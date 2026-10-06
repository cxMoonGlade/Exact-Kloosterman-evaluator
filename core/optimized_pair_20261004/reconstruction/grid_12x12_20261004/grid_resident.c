/* Grid-capacity extension of D15: owned public FLINT objects for the D11 ordered trace.
   No Python layout, borrowed coefficient pointer, or cross-query math cache. */
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"
#include "fmpz_mat.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

_Static_assert(CHAR_BIT == 8 && sizeof(long) == 8 && sizeof(slong) == 8,
               "D15 requires the registered LP64 platform");
_Static_assert(sizeof(ulong) == 8 && sizeof(size_t) == 8 && sizeof(void *) == 8,
               "D15 requires 64-bit words, sizes and pointers");
_Static_assert(sizeof(fmpz) == 8 && sizeof(mp_limb_t) == 8 && sizeof(int) == 4,
               "D15 requires the registered FLINT/GMP layout");
_Static_assert(FLINT_BITS == 64 && GMP_NAIL_BITS == 0,
               "D15 requires 64-bit GMP limbs without nails");

enum { OK, SHAPE, CAPACITY, CANONICAL, STATE, ALLOC, ALGEBRA, ABI,
       NOT_AVAILABLE, COUNTER };
enum { READY=1, RUNNING, DONE, FAILED };
enum { TRACE, INITIAL_A, LEFT, RIGHT, PRODUCT, JUMP_BEFORE, JUMP_AFTER,
       FINAL_MATRIX };
#define COUNTS(X) \
 X(creates) X(runs) X(exports) X(steps) X(full_products) X(terminal_products) \
 X(zero_products) X(real_products) X(one_imag_products) X(gaussian_products) \
 X(ks_calls) X(mod_g_calls) X(compose_schedules) X(vec_calls) \
 X(compose_components) X(compose_arguments) X(constant_components) \
 X(inflate_components) X(inflate_arguments) X(argument_squares) \
 X(hinv_builds) X(hinv_checks) X(modular_copy_coefficients) \
 X(modular_to_integer_coefficients) X(wire_input_fields) X(wire_output_fields) \
 X(wire_input_bytes) X(wire_output_bytes) X(debug_bytes) \
 X(poly_inits) X(poly_clears) X(poly_mat_inits) X(poly_mat_clears) \
 X(ctx_inits) X(ctx_clears) X(bridge_allocs) X(bridge_frees) \
 X(max_ks_bits) X(max_ks_span_bits) X(max_ks_limbs) X(max_vec_slots)
#define COUNT_ENUM(s) C_##s,
enum { COUNTS(COUNT_ENUM) C_COUNT };
#undef COUNT_ENUM
#define COUNT_NAME(s) #s,
static const char *const count_names[] = { COUNTS(COUNT_NAME) };
#undef COUNT_NAME
enum { T_CONTROL, T_PARSE, T_CONTEXT, T_COPY, T_HINV, T_COMPOSE_PREPARE,
       T_COMPOSE_KERNEL, T_INFLATE, T_JUMP, T_KS_PREPARE, T_KS_KERNEL,
       T_COMBINE_REDUCE, T_DEBUG, T_EXPORT, T_CLEANUP, T_COUNT };
static const char *const time_names[] = {
    "control", "parse", "context_setup", "owner_copy", "hinv_setup",
    "compose_prepare", "compose_kernel", "inflate", "jump", "ks_prepare",
    "ks_kernel", "combine_reduce", "debug", "export", "cleanup"
};
typedef struct {
    long abi, state, status, N, L, d, P, e, audit, expected_products;
    long steps_done, hinv_created;
} metadata;
typedef struct { fmpz_mod_poly_t p; int live; } mpoly;
typedef struct { fmpz_poly_mat_t p; int live; } imatrix;
typedef struct {
    long rows, cols, initialized;
    fmpz_mod_poly_struct *parts;
} matrix;
typedef struct {
    long rows, cols, parts, fields;
    size_t length;
    unsigned char *bytes;
} snapshot;
typedef struct {
    long kind, phase, digit, terminal, rows, inner, cols, branch;
    long native_products, mod_g_calls, ks_max_bits, ks_max_span_bits;
    snapshot left, right, product, before, after;
} step_record;
typedef struct {
    metadata meta;
    long counts[C_COUNT];
    double seconds[T_COUNT], total;
    long kind, index, rows, cols, parts, degree, bits, fields;
    size_t length;
    unsigned char *bytes;
} result;
typedef struct {
    metadata meta;
    long counts[C_COUNT], H, limbs;
    double seconds[T_COUNT], total, phase_started, api_started;
    int active, context_live, scalar_live;
    fmpz_mod_ctx_t ctx;
    fmpz_t scalar;
    mpoly g, x, hinv, trace[2];
    mpoly *jump, *pending;
    matrix *A, *current, *shifted, *temporary, *base;
    ulong *limb_buffer;
    step_record *step;
    snapshot initial;
    result *receipt;
    size_t debug_limit;
} owner;

static double now_seconds(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec * 1e-9;
}
static void phase(owner *o, int next) {
    double t=now_seconds();
    if (o->active>=0) o->seconds[o->active]+=t-o->phase_started;
    o->active=next; o->phase_started=t;
}
static void api_begin(owner *o) {
    o->api_started=now_seconds(); o->phase_started=o->api_started;
    o->active=T_CONTROL;
}
static void api_end(owner *o) {
    phase(o,-1); o->total+=o->phase_started-o->api_started;
}
static int add_count(owner *o, int key, long amount) {
    if (amount<0 || o->counts[key]>LONG_MAX-amount) {
        if (!o->meta.status) o->meta.status=COUNTER;
        return COUNTER;
    }
    o->counts[key]+=amount; return OK;
}
#define COUNT(o,key,n) do { int cs_=add_count((o),C_##key,(n)); \
    if (cs_) return cs_; } while (0)
static void max_count(owner *o, int key, long value) {
    if (o->counts[key]<value) o->counts[key]=value;
}
static int size_add(size_t a, size_t b, size_t *out) {
    if (a>(size_t)PTRDIFF_MAX || b>(size_t)PTRDIFF_MAX-a) return 0;
    *out=a+b; return 1;
}
static int size_mul(size_t a, size_t b, size_t *out) {
    if (b && a>(size_t)PTRDIFF_MAX/b) return 0;
    *out=a*b; return 1;
}
static int long_add(long a, long b, long *out) {
    if (a<0 || b<0 || a>LONG_MAX-b) return 0;
    *out=a+b; return 1;
}
static int long_mul(long a, long b, long *out) {
    if (a<0 || b<0 || (b && a>LONG_MAX/b)) return 0;
    *out=a*b; return 1;
}
static long bit_length(ulong v) {
    long n=0; while (v) { n++; v>>=1; } return n;
}
static long products(long v) {
    long bits=bit_length((ulong)v), ones=0;
    for (ulong t=(ulong)v;t;t>>=1) ones+=(long)(t&1);
    return bits+ones-2;
}
static int cells(long rows, long cols, size_t item) {
    long n; size_t bytes;
    return long_mul(rows,cols,&n) && size_mul((size_t)n,item,&bytes);
}
static int coeff_capacity(long len) {
    long reserve;
    return long_mul(len,4,&reserve) && cells(reserve,1,sizeof(fmpz));
}
static int limb_capacity(long bits, long *count) {
    long n, extra, span; size_t bytes;
    if (bits<0) return 0;
    n=bits ? (bits-1)/64+1 : 0;
    if (!long_add(n,2,&extra) || extra>INT_MAX
        || !long_mul(extra,64,&span)
        || !size_mul((size_t)extra,sizeof(mp_limb_t),&bytes)) return 0;
    if (count) *count=n;
    return 1;
}
static int matching_abi(void) {
    return !strcmp(FLINT_VERSION,"3.6.0") && !strcmp(flint_version,"3.6.0");
}
int qcb_resident_abi(void) { return 1; }
const char *qcb_resident_header_version(void) { return FLINT_VERSION; }
const char *qcb_resident_runtime_version(void) { return flint_version; }
size_t qcb_resident_slong_size(void) { return sizeof(slong); }

/* All public-domain/size decisions precede owner allocation and byte reads. */
static int dimensions(long N, long L, long d, long P, int audit,
                      size_t gb, size_t ab, size_t *debug_limit) {
    long square, twice, fields, S, H, c, dc, baby, blocks, bc, sum, t;
    size_t gbytes, abytes, bytes, total=0;
    if (N<1 || N>128 || L<2 || L>128 || P<1 || P>8201 || d<1 || d>N
        || N%d || (audit!=0 && audit!=1)) return SHAPE;
    H=P/8+(P%8!=0); S=products(d)+products(N/d);
    if (!long_mul(L,L,&square) || !long_mul(2,square,&twice)
        || !long_mul(twice,d,&fields)
        || !size_mul((size_t)(d+1),(size_t)H,&gbytes)
        || !size_mul((size_t)fields,(size_t)H,&abytes)) return CAPACITY;
    if (gb!=gbytes || ab!=abytes) return SHAPE;
    if (!size_add(gbytes,abytes,&bytes) || !coeff_capacity(2*(d+1))
        || !cells(twice,1,sizeof(fmpz_mod_poly_struct))
        || !cells(fields,1,sizeof(fmpz))
        || !cells(S,1,sizeof(step_record))) return CAPACITY;
    /* Maximum actual vector batch, with both parts and one jump. */
    c=twice+1;
    if (!long_mul(d,c,&dc)) return CAPACITY;
    baby=(long)n_sqrt((ulong)dc)+1; blocks=(d+1)/baby+1;
    if (!long_mul(blocks,c,&bc) || !cells(baby,d,sizeof(fmpz))
        || !cells(bc,baby,sizeof(fmpz)) || !cells(bc,d,sizeof(fmpz))
        || !cells(baby,1,sizeof(fmpz *))
        || !cells(2*c+1,1,sizeof(fmpz_mod_poly_struct))
        || !cells((2*c+1)*d,1,sizeof(fmpz))) return CAPACITY;
    if (!long_mul(6*S,square,&sum) || !long_mul(sum,d,&sum)
        || !long_mul(2*S,d,&t) || !long_add(sum,t,&sum)
        || !long_add(sum,fields,&sum)
        || !size_mul((size_t)sum,(size_t)H,&total)) return CAPACITY;
    *debug_limit=total;
    return OK;
}

/* Receipt/export allocations are deliberately outside these owner counters. */
static void *bridge_alloc(owner *o, size_t count, size_t item) {
    size_t bytes; void *p;
    if (!size_mul(count,item,&bytes)) { o->meta.status=CAPACITY; return NULL; }
    if (!bytes) return NULL;
    p=calloc(count,item);
    if (!p) { o->meta.status=ALLOC; return NULL; }
    if (add_count(o,C_bridge_allocs,1)) { free(p); return NULL; }
    return p;
}
static void bridge_free(owner *o, void *p) {
    if (p) { free(p); add_count(o,C_bridge_frees,1); }
}
static void mp_init(owner *o, mpoly *p) {
    fmpz_mod_poly_init(p->p,o->ctx); p->live=1;
    add_count(o,C_poly_inits,1);
}
static void mp_clear(owner *o, mpoly *p) {
    if (p->live) {
        fmpz_mod_poly_clear(p->p,o->ctx); p->live=0;
        add_count(o,C_poly_clears,1);
    }
}
static int mp_copy(owner *o, fmpz_mod_poly_t out, const fmpz_mod_poly_t in) {
    fmpz_mod_poly_set(out,in,o->ctx);
    COUNT(o,modular_copy_coefficients,fmpz_mod_poly_length(in,o->ctx));
    return o->meta.status;
}
static matrix *matrix_new(owner *o, long rows, long cols) {
    matrix *m; long n;
    if (!long_mul(rows,cols,&n) || !long_mul(n,2,&n)
        || !cells(n,1,sizeof(fmpz_mod_poly_struct))) {
        o->meta.status=CAPACITY; return NULL;
    }
    m=bridge_alloc(o,1,sizeof(*m)); if (!m) return NULL;
    m->rows=rows; m->cols=cols;
    m->parts=bridge_alloc(o,(size_t)n,sizeof(fmpz_mod_poly_struct));
    if (!m->parts) { bridge_free(o,m); return NULL; }
    for (long j=0;j<n;j++) {
        fmpz_mod_poly_init(m->parts+j,o->ctx); m->initialized++;
        add_count(o,C_poly_inits,1);
    }
    return m;
}
static void matrix_clear(owner *o, matrix **slot) {
    matrix *m=*slot;
    if (m) {
        for (long j=0;j<m->initialized;j++) {
            fmpz_mod_poly_clear(m->parts+j,o->ctx);
            add_count(o,C_poly_clears,1);
        }
        bridge_free(o,m->parts); bridge_free(o,m); *slot=NULL;
    }
}
static matrix *matrix_copy(owner *o, const matrix *in) {
    matrix *out=matrix_new(o,in->rows,in->cols);
    if (!out) return NULL;
    for (long j=0;j<in->initialized;j++) {
        if (mp_copy(o,out->parts+j,in->parts+j)) {
            matrix_clear(o,&out); return NULL;
        }
    }
    return out;
}
static int parse_field(owner *o, const unsigned char *bytes, fmpz_t value) {
    memset(o->limb_buffer,0,(size_t)o->limbs*sizeof(ulong));
    for (long j=0;j<o->H;j++) {
        unsigned char b=bytes[j];
        if (o->meta.P%8 && j==o->H-1 && (b>>(o->meta.P%8))) return CANONICAL;
        o->limb_buffer[j/8]|=(ulong)b<<(8*(j%8));
    }
    fmpz_set_ui_array(value,o->limb_buffer,o->limbs);
    COUNT(o,wire_input_fields,1); COUNT(o,wire_input_bytes,o->H);
    return OK;
}
static int write_field(owner *o, unsigned char *bytes, const fmpz_t value) {
    if (fmpz_sgn(value)<0 || fmpz_bits(value)>(ulong)o->meta.P) return ALGEBRA;
    fmpz_get_ui_array(o->limb_buffer,o->limbs,value);
    for (long j=0;j<o->H;j++) bytes[j]=(unsigned char)(o->limb_buffer[j/8]>>(8*(j%8)));
    return OK;
}
static int snapshot_make(owner *o, snapshot *s, long rows, long cols,
                         long parts, const fmpz_mod_poly_struct *polys) {
    long fields, count; size_t length, end;
    int status;
    if (s->bytes || !long_mul(rows,cols,&count) || !long_mul(count,parts,&count)
        || !long_mul(count,o->meta.d,&fields)
        || !size_mul((size_t)fields,(size_t)o->H,&length)
        || !size_add((size_t)o->counts[C_debug_bytes],length,&end)
        || end>o->debug_limit) return CAPACITY;
    s->rows=rows; s->cols=cols; s->parts=parts; s->fields=fields; s->length=length;
    s->bytes=bridge_alloc(o,length,1); if (!s->bytes) return o->meta.status;
    size_t offset=0;
    for (long j=0;j<count;j++) for (long k=0;k<o->meta.d;k++) {
        fmpz_mod_poly_get_coeff_fmpz(o->scalar,polys+j,k,o->ctx);
        status=write_field(o,s->bytes+offset,o->scalar); if (status) return status;
        offset+=(size_t)o->H;
    }
    COUNT(o,debug_bytes,(long)length);
    return OK;
}
static int snapshot_matrix(owner *o, snapshot *s, const matrix *m) {
    return snapshot_make(o,s,m->rows,m->cols,2,m->parts);
}
static void snapshot_clear(owner *o, snapshot *s) {
    bridge_free(o,s->bytes); memset(s,0,sizeof(*s));
}

static int im_init(owner *o, imatrix *m, long rows, long cols) {
    if (!cells(rows,cols,sizeof(fmpz_poly_struct))
        || !cells(rows,cols,sizeof(fmpz))) return CAPACITY;
    fmpz_poly_mat_init(m->p,rows,cols); m->live=1;
    COUNT(o,poly_mat_inits,1); return OK;
}
static void im_clear(owner *o, imatrix *m) {
    if (m->live) {
        fmpz_poly_mat_clear(m->p); m->live=0;
        add_count(o,C_poly_mat_clears,1);
    }
}
/* All fit requests here have a bounded polynomial role, not unbounded reuse. */
static int fit_mod(owner *o, fmpz_mod_poly_t p, long target) {
    long max=target>p->alloc ? target : p->alloc, reserve;
    if (!long_mul(max,2,&reserve) || !coeff_capacity(reserve)) return CAPACITY;
    fmpz_mod_poly_fit_length(p,target,o->ctx); return OK;
}
static int mod_rem(owner *o, fmpz_mod_poly_t out,
                   const fmpz_mod_poly_t input, int matrix_output) {
    fmpz_t factor; int status=OK;
    fmpz_init(factor);
    fmpz_mod_poly_rem_f(factor,out,input,o->g.p,o->ctx);
    if (matrix_output) add_count(o,C_mod_g_calls,1);
    if (!fmpz_is_one(factor) || fmpz_mod_poly_degree(out,o->ctx)>=o->meta.d)
        status=ALGEBRA;
    fmpz_clear(factor); return status ? status : (int)o->meta.status;
}
static int branch_of(owner *o, const matrix *a, const matrix *b,
                     int *ai, int *bi) {
    int az=1,bz=1;
    *ai=*bi=0;
    for (long j=0;j<a->initialized;j++) if (!fmpz_mod_poly_is_zero(a->parts+j,o->ctx)) {
        az=0; if (j%2) *ai=1;
    }
    for (long j=0;j<b->initialized;j++) if (!fmpz_mod_poly_is_zero(b->parts+j,o->ctx)) {
        bz=0; if (j%2) *bi=1;
    }
    return az || bz ? 0 : *ai && *bi ? 3 : *ai || *bi ? 2 : 1;
}
static int integer_matrix(owner *o, imatrix *out, const matrix *in,
                           int part, int terminal, int right) {
    long L=o->meta.L;
    long rows=terminal ? (right ? L*L : 1) : L;
    long cols=terminal ? (right ? 1 : L*L) : L;
    int status=im_init(o,out,rows,cols); if (status) return status;
    for (long i=0;i<L;i++) for (long j=0;j<L;j++) {
        long ii=terminal ? (right ? i*L+j : 0) : i;
        long jj=terminal ? (right ? 0 : i*L+j) : j;
        long source=terminal && right ? j*L+i : i*L+j;
        const fmpz_mod_poly_struct *p=in->parts+2*source+part;
        fmpz_mod_poly_get_fmpz_poly(fmpz_poly_mat_entry(out->p,ii,jj),p,o->ctx);
        COUNT(o,modular_to_integer_coefficients,fmpz_mod_poly_length(p,o->ctx));
    }
    return OK;
}
/* Capacity mirrors actual KS max lengths/width and guards all unpack rounding.
   Its span is an upper bound, not measured GMP allocation or output bit count. */
static int ks_capacity(owner *o, const imatrix *a, const imatrix *b,
                       long *width, long *span, long *limbs) {
    long la=fmpz_poly_mat_max_length(a->p), lb=fmpz_poly_mat_max_length(b->p);
    long ba=fmpz_poly_mat_max_bits(a->p), bb=fmpz_poly_mat_max_bits(b->p);
    long inner=fmpz_poly_mat_nrows(b->p), d=o->meta.d;
    long v,V,T,qa,packed,rounded,maxlimbs=0,z;
    if (la<0 || lb<0 || la>d || lb>d || ba<0 || bb<0
        || ba>o->meta.P+1 || bb>o->meta.P+1 || inner<1) return ALGEBRA;
    if (!long_add(ba,bb,&v)
        || !long_add(v,bit_length((ulong)(la<lb?la:lb)),&v)
        || !long_add(v,bit_length((ulong)inner),&v)
        || !long_mul(2,o->meta.P+1,&V) || !long_add(V,1,&V)
        || !long_add(V,bit_length((ulong)d),&V)
        || !long_add(V,bit_length((ulong)inner),&V)
        || v<1 || v>V || !long_mul(2*d,V,&T)
        || !long_add(T,V-1,&rounded)) return CAPACITY;
    qa=(T-1)/v+1;
    if (!long_add(qa,1,&z) || !coeff_capacity(z)
        || !long_mul(qa,v,&packed) || !limb_capacity(packed,&maxlimbs)
        || !long_add(T,v-1,&rounded) || !limb_capacity(T,&z)) return CAPACITY;
    if (z>maxlimbs) maxlimbs=z;
    if (!long_mul(la,v,&packed) || !limb_capacity(packed,&z)) return CAPACITY;
    if (z>maxlimbs) maxlimbs=z;
    if (!long_mul(lb,v,&packed) || !limb_capacity(packed,&z)) return CAPACITY;
    if (z>maxlimbs) maxlimbs=z;
    if (!limb_capacity(o->meta.P+1,&z)) return CAPACITY;
    if (z>maxlimbs) maxlimbs=z;
    if (!coeff_capacity(2*d)) return CAPACITY;
    *width=v; *span=T; *limbs=maxlimbs; return OK;
}
static int ks_product(owner *o, imatrix *out, const imatrix *a,
                       const imatrix *b, step_record *s) {
    long width,span,limbs; int status;
    phase(o,T_KS_PREPARE);
    if (fmpz_poly_mat_ncols(a->p)!=fmpz_poly_mat_nrows(b->p)) return ALGEBRA;
    status=ks_capacity(o,a,b,&width,&span,&limbs); if (status) return status;
    status=im_init(o,out,fmpz_poly_mat_nrows(a->p),fmpz_poly_mat_ncols(b->p));
    if (status) return status;
    phase(o,T_KS_KERNEL);
    fmpz_poly_mat_mul_KS(out->p,a->p,b->p);
    COUNT(o,ks_calls,1); s->native_products++;
    max_count(o,C_max_ks_bits,width); max_count(o,C_max_ks_span_bits,span);
    max_count(o,C_max_ks_limbs,limbs);
    if (s->ks_max_bits<width) s->ks_max_bits=width;
    if (s->ks_max_span_bits<span) s->ks_max_span_bits=span;
    phase(o,T_KS_PREPARE);
    if (fmpz_poly_mat_max_length(out->p)>2*o->meta.d-1) return ALGEBRA;
    return OK;
}

/* The branch is fixed by inputs. A zero output imaginary part still gets rem. */
static int matrix_product(owner *o, const matrix *a, const matrix *b,
                           step_record *s, matrix **output) {
    enum { AR, AI, BR, BI, SA, SB, U, V, W, E, I, IM_COUNT };
    imatrix m[IM_COUNT]={0};
    mpoly converted={0};
    int status=OK, ai=0, bi=0, branch;
    long outdim=s->terminal ? 1 : o->meta.L;
    matrix *out=NULL;
    *output=NULL;
    phase(o,T_KS_PREPARE);
    if (!a || !b || a->rows!=o->meta.L || a->cols!=o->meta.L
        || b->rows!=o->meta.L || b->cols!=o->meta.L) return ALGEBRA;
    branch=branch_of(o,a,b,&ai,&bi); s->branch=branch;
    out=matrix_new(o,outdim,outdim); if (!out) return (int)o->meta.status;
    if (!branch) goto success;
#define MP_TRY(expr) do { status=(expr); if (status) goto finish; } while (0)
    MP_TRY(integer_matrix(o,&m[AR],a,0,(int)s->terminal,0));
    MP_TRY(integer_matrix(o,&m[BR],b,0,(int)s->terminal,1));
    if (ai) MP_TRY(integer_matrix(o,&m[AI],a,1,(int)s->terminal,0));
    if (bi) MP_TRY(integer_matrix(o,&m[BI],b,1,(int)s->terminal,1));
    MP_TRY(ks_product(o,&m[U],&m[AR],&m[BR],s));
    if (ai && bi) {
        MP_TRY(ks_product(o,&m[V],&m[AI],&m[BI],s));
        phase(o,T_KS_PREPARE);
        MP_TRY(im_init(o,&m[SA],fmpz_poly_mat_nrows(m[AR].p),fmpz_poly_mat_ncols(m[AR].p)));
        MP_TRY(im_init(o,&m[SB],fmpz_poly_mat_nrows(m[BR].p),fmpz_poly_mat_ncols(m[BR].p)));
        fmpz_poly_mat_add(m[SA].p,m[AR].p,m[AI].p);
        fmpz_poly_mat_add(m[SB].p,m[BR].p,m[BI].p);
        MP_TRY(ks_product(o,&m[W],&m[SA],&m[SB],s));
        phase(o,T_COMBINE_REDUCE);
        MP_TRY(im_init(o,&m[E],outdim,outdim));
        MP_TRY(im_init(o,&m[I],outdim,outdim));
        fmpz_poly_mat_sub(m[E].p,m[U].p,m[V].p);
        fmpz_poly_mat_sub(m[I].p,m[W].p,m[U].p);
        /* Public matrix sub is entrywise and permits its output as an input. */
        fmpz_poly_mat_sub(m[I].p,m[I].p,m[V].p);
    } else if (ai || bi) {
        MP_TRY(ks_product(o,&m[I],ai ? &m[AI] : &m[AR],bi ? &m[BI] : &m[BR],s));
    }
    phase(o,T_COMBINE_REDUCE);
    mp_init(o,&converted);
    MP_TRY(fit_mod(o,converted.p,2*o->meta.d));
    for (long row=0;row<outdim;row++) for (long col=0;col<outdim;col++) {
        for (int part=0;part<(ai||bi ? 2 : 1);part++) {
            imatrix *src=part ? &m[I] : ai&&bi ? &m[E] : &m[U];
            fmpz_mod_poly_set_fmpz_poly(converted.p,fmpz_poly_mat_entry(src->p,row,col),o->ctx);
            MP_TRY(mod_rem(o,out->parts+2*(row*outdim+col)+part,converted.p,1));
            s->mod_g_calls++;
        }
    }
success:
    status=(int)o->meta.status;
    if (!status) { *output=out; out=NULL; }
finish:
    phase(o,T_CLEANUP);
    mp_clear(o,&converted);
    for (int j=IM_COUNT-1;j>=0;j--) im_clear(o,&m[j]);
    matrix_clear(o,&out);
    phase(o,T_CONTROL);
#undef MP_TRY
    return status ? status : (int)o->meta.status;
}

static int hinv_setup(owner *o) {
    mpoly reversed={0},check={0};
    fmpz_t factor;
    int status=OK;
    if (o->hinv.live) return OK;
    phase(o,T_HINV);
    fmpz_init(factor); mp_init(o,&reversed); mp_init(o,&check); mp_init(o,&o->hinv);
    status=fit_mod(o,reversed.p,o->meta.d+1); if (status) goto finish;
    status=fit_mod(o,check.p,2*(o->meta.d+1)); if (status) goto finish;
    status=fit_mod(o,o->hinv.p,o->meta.d+1); if (status) goto finish;
    fmpz_mod_poly_reverse(reversed.p,o->g.p,o->meta.d+1,o->ctx);
    fmpz_mod_poly_inv_series_f(factor,o->hinv.p,reversed.p,o->meta.d+1,o->ctx);
    add_count(o,C_hinv_builds,1);
    if (!fmpz_is_one(factor)) { status=ALGEBRA; goto finish; }
    fmpz_mod_poly_mullow(check.p,reversed.p,o->hinv.p,o->meta.d+1,o->ctx);
    add_count(o,C_hinv_checks,1);
    if (!fmpz_mod_poly_is_one(check.p,o->ctx)) { status=ALGEBRA; goto finish; }
    o->meta.hinv_created=1;
finish:
    phase(o,T_CLEANUP);
    mp_clear(o,&check); mp_clear(o,&reversed); fmpz_clear(factor);
    phase(o,T_CONTROL);
    return status ? status : (int)o->meta.status;
}
static int vector_capacity(owner *o, long c, long *slots) {
    long d=o->meta.d, dc, baby, blocks, bc, s0,s1,s2,total,bits,t;
    size_t bytes, extra;
    if (c<1 || c>2*o->meta.L*o->meta.L+1 || !long_mul(d,c,&dc)) return CAPACITY;
    baby=(long)n_sqrt((ulong)dc)+1; blocks=(d+1)/baby+1;
    if (!long_mul(blocks,c,&bc) || !long_mul(baby,d,&s0)
        || !long_mul(bc,baby,&s1) || !long_mul(bc,d,&s2)
        || !long_add(s0,s1,&total) || !long_add(total,s2,&total)
        || !cells(baby,d,sizeof(fmpz)) || !cells(bc,baby,sizeof(fmpz))
        || !cells(bc,d,sizeof(fmpz)) || !cells(total,1,sizeof(fmpz))
        || !cells(baby,1,sizeof(fmpz *)) || !cells(2*d,1,sizeof(fmpz))
        || !cells(2*c+1,1,sizeof(fmpz_mod_poly_struct))
        || !cells((2*c+1)*d,1,sizeof(fmpz)) || !cells(c,1,sizeof(long))
        || !cells(3*d-2,1,sizeof(fmpz))) return CAPACITY;
    if (!size_mul((size_t)total,sizeof(fmpz),&bytes)
        || !size_mul((size_t)baby,sizeof(fmpz *),&extra)
        || !size_add(bytes,extra,&bytes)
        || !size_mul((size_t)(2*d),sizeof(fmpz),&extra)
        || !size_add(bytes,extra,&bytes)) return CAPACITY;
    if (!long_mul(2,o->meta.P,&bits)
        || !long_add(bits,bit_length((ulong)(baby-1)),&bits)
        || !long_add(bits,1,&bits) || !limb_capacity(bits,&t)) return CAPACITY;
    *slots=total;
    return OK;
}
static int inflate_one(owner *o, fmpz_mod_poly_t out,
                        const fmpz_mod_poly_t input, mpoly *scratch) {
    long len=fmpz_mod_poly_length(input,o->ctx);
    long target=len<=1 ? len : 2*(len-1)+1;
    int status=fit_mod(o,scratch->p,target); if (status) return status;
    fmpz_mod_poly_inflate(scratch->p,input,2,o->ctx);
    return mod_rem(o,out,scratch->p,0);
}
static int compose_matrix(owner *o, const matrix *source,
                           int include_jump, int first) {
    fmpz_mod_poly_struct *inputs=NULL,*results=NULL;
    long *positions=NULL, input_init=0,result_init=0,c=0,u=0,vec_slots=0;
    mpoly inflated={0};
    int status=OK;
    if (!source || o->shifted || !o->jump || !o->jump->live || o->meta.d<2)
        return STATE;
    COUNT(o,compose_schedules,1);
    phase(o,first ? T_INFLATE : T_COMPOSE_PREPARE);
    o->shifted=matrix_new(o,o->meta.L,o->meta.L);
    if (!o->shifted) return (int)o->meta.status;
    if (include_jump && !o->pending->live) mp_init(o,o->pending);
    if (first) {
        mp_init(o,&inflated);
        for (long j=0;j<source->initialized;j++) {
            if (!fmpz_mod_poly_is_zero(source->parts+j,o->ctx)) {
                status=inflate_one(o,o->shifted->parts+j,source->parts+j,&inflated);
                if (status) goto finish;
                add_count(o,C_inflate_components,1);
            }
        }
        if (include_jump) {
            status=inflate_one(o,o->pending->p,o->jump->p,&inflated);
            if (status) goto finish;
            add_count(o,C_inflate_arguments,1);
        }
        goto finish;
    }
    /* Keep the original first-general-schedule setup, including empty batches. */
    status=hinv_setup(o); if (status) goto finish;
    phase(o,T_COMPOSE_PREPARE);
    for (long j=0;j<source->initialized;j++)
        if (!fmpz_mod_poly_is_zero(source->parts+j,o->ctx)) u++;
    c=u+include_jump;
    add_count(o,C_compose_components,u);
    add_count(o,C_compose_arguments,include_jump);
    if (!c) goto finish;
    status=vector_capacity(o,c,&vec_slots); if (status) goto finish;
    inputs=bridge_alloc(o,(size_t)(c+1),sizeof(*inputs));
    results=bridge_alloc(o,(size_t)c,sizeof(*results));
    positions=bridge_alloc(o,(size_t)c,sizeof(*positions));
    if (!inputs || !results || !positions) { status=(int)o->meta.status; goto finish; }
    for (long j=0;j<=c;j++) {
        fmpz_mod_poly_init(inputs+j,o->ctx); input_init++; add_count(o,C_poly_inits,1);
    }
    for (long j=0;j<c;j++) {
        fmpz_mod_poly_init(results+j,o->ctx); result_init++; add_count(o,C_poly_inits,1);
    }
    long at=0;
    for (long j=0;j<source->initialized;j++) {
        if (!fmpz_mod_poly_is_zero(source->parts+j,o->ctx)) {
            positions[at]=j;
            status=mp_copy(o,inputs+at,source->parts+j); if (status) goto finish;
            at++;
        }
    }
    if (include_jump) {
        positions[at]=-1;
        status=mp_copy(o,inputs+at,o->jump->p); if (status) goto finish;
        at++;
    }
    if (at!=c || !fmpz_mod_poly_is_zero(inputs+c,o->ctx)) { status=ALGEBRA; goto finish; }
    phase(o,T_COMPOSE_KERNEL);
    if (fmpz_mod_poly_is_zero(o->jump->p,o->ctx)) {
        for (long j=0;j<c;j++) {
            fmpz_mod_poly_get_coeff_fmpz(o->scalar,inputs+j,0,o->ctx);
            fmpz_mod_poly_set_coeff_fmpz(results+j,0,o->scalar,o->ctx);
            add_count(o,C_constant_components,1);
        }
    } else {
        fmpz_mod_poly_compose_mod_brent_kung_vec_preinv(results,inputs,c+1,c,
            o->jump->p,o->g.p,o->hinv.p,o->ctx);
        add_count(o,C_vec_calls,1);
        max_count(o,C_max_vec_slots,vec_slots);
    }
    phase(o,T_COMPOSE_PREPARE);
    for (long j=0;j<c;j++) {
        if (fmpz_mod_poly_degree(results+j,o->ctx)>=o->meta.d) { status=ALGEBRA; goto finish; }
        if (positions[j]<0) status=mp_copy(o,o->pending->p,results+j);
        else status=mp_copy(o,o->shifted->parts+positions[j],results+j);
        if (status) goto finish;
    }
finish:
    phase(o,T_CLEANUP);
    for (long j=0;j<result_init;j++) {
        fmpz_mod_poly_clear(results+j,o->ctx); add_count(o,C_poly_clears,1);
    }
    for (long j=0;j<input_init;j++) {
        fmpz_mod_poly_clear(inputs+j,o->ctx); add_count(o,C_poly_clears,1);
    }
    bridge_free(o,positions); bridge_free(o,results); bridge_free(o,inputs);
    mp_clear(o,&inflated);
    phase(o,T_CONTROL);
    return status ? status : (int)o->meta.status;
}

static int initial_jump(owner *o) {
    mpoly raw={0}; int status=OK;
    phase(o,T_JUMP);
    o->jump=bridge_alloc(o,1,sizeof(*o->jump));
    o->pending=bridge_alloc(o,1,sizeof(*o->pending));
    if (!o->jump || !o->pending) return (int)o->meta.status;
    mp_init(o,&o->x); mp_init(o,o->jump); mp_init(o,&raw);
    fmpz_mod_poly_set_coeff_ui(raw.p,1,1,o->ctx);
    status=mod_rem(o,o->x.p,raw.p,0); if (status) goto finish;
    fmpz_mod_poly_mulmod(o->jump->p,o->x.p,o->x.p,o->g.p,o->ctx);
    add_count(o,C_argument_squares,1);
finish:
    phase(o,T_CLEANUP); mp_clear(o,&raw); phase(o,T_CONTROL);
    return status ? status : (int)o->meta.status;
}
static int square_jump(owner *o) {
    mpoly next={0}; int status=OK;
    phase(o,T_JUMP); mp_init(o,&next);
    fmpz_mod_poly_mulmod(next.p,o->jump->p,o->jump->p,o->g.p,o->ctx);
    add_count(o,C_argument_squares,1);
    status=mp_copy(o,o->jump->p,next.p);
    phase(o,T_CLEANUP); mp_clear(o,&next); phase(o,T_CONTROL);
    return status ? status : (int)o->meta.status;
}

/* Audit operands are logical L by L; only the multiplication kernel flattens. */
static int do_step(owner *o, const matrix *left, const matrix *right,
                    int kind, int digit) {
    long index=o->meta.steps_done;
    int status=OK;
    step_record *s;
    if (index<0 || index>=o->meta.expected_products || o->temporary) return ALGEBRA;
    s=o->step+index;
    s->kind=kind; s->phase=kind<2 ? 0 : 1; s->digit=digit;
    s->terminal=index+1==o->meta.expected_products;
    s->rows=s->cols=s->terminal ? 1 : o->meta.L;
    s->inner=s->terminal ? o->meta.L*o->meta.L : o->meta.L;
    if (o->meta.audit) {
        phase(o,T_DEBUG);
        status=snapshot_matrix(o,&s->left,left); if (status) return status;
        status=snapshot_matrix(o,&s->right,right); if (status) return status;
        if (kind<2) {
            status=snapshot_make(o,&s->before,1,1,1,o->jump->p);
            if (status) return status;
        }
    }
    status=matrix_product(o,left,right,s,&o->temporary); if (status) return status;
    if (o->meta.audit) {
        phase(o,T_DEBUG);
        status=snapshot_matrix(o,&s->product,o->temporary); if (status) return status;
    }
    phase(o,T_COMBINE_REDUCE);
    if (s->terminal) {
        for (int j=0;j<2;j++) {
            status=mp_copy(o,o->trace[j].p,o->temporary->parts+j);
            if (status) return status;
        }
    }
    /* The old current has exactly one owning role, even for its square. */
    phase(o,T_CLEANUP); matrix_clear(o,&o->current);
    if (s->terminal) matrix_clear(o,&o->temporary);
    else { o->current=o->temporary; o->temporary=NULL; }
    matrix_clear(o,&o->shifted);
    phase(o,T_COMBINE_REDUCE);
    if (kind==0) {
        mpoly *old=o->jump;
        if (!o->pending || !o->pending->live) return ALGEBRA;
        o->jump=o->pending; o->pending=old;
        phase(o,T_CLEANUP); mp_clear(o,o->pending);
    } else if (kind==1) {
        status=square_jump(o); if (status) return status;
    }
    if (o->meta.audit && kind<2) {
        phase(o,T_DEBUG);
        status=snapshot_make(o,&s->after,1,1,1,o->jump->p); if (status) return status;
    }
    phase(o,T_CONTROL);
    COUNT(o,steps,1); o->meta.steps_done++;
    if (s->terminal) { COUNT(o,terminal_products,1); }
    else { COUNT(o,full_products,1); }
    static const int keys[]={C_zero_products,C_real_products,C_one_imag_products,C_gaussian_products};
    return add_count(o,keys[s->branch],1);
}

int qcb_resident_run(void *opaque) {
    owner *o=opaque;
    int status=OK;
    if (!o) return STATE;
    api_begin(o);
    if (o->meta.state!=READY) { api_end(o); return STATE; }
    o->meta.state=RUNNING; add_count(o,C_runs,1);
    if (o->meta.d>1) {
        status=initial_jump(o); if (status) goto finish;
        int first=1;
        for (long bit=bit_length((ulong)o->meta.d)-2;bit>=0;bit--) {
            int digit=(int)(((ulong)o->meta.d>>bit)&1);
            status=compose_matrix(o,o->current,1,first); if (status) goto finish;
            status=do_step(o,o->shifted,o->current,0,digit); if (status) goto finish;
            first=0;
            if (digit) {
                status=compose_matrix(o,o->A,0,0); if (status) goto finish;
                status=do_step(o,o->shifted,o->current,1,1); if (status) goto finish;
            }
        }
        phase(o,T_JUMP);
        if (!fmpz_mod_poly_equal(o->jump->p,o->x.p,o->ctx)) { status=ALGEBRA; goto finish; }
        phase(o,T_CLEANUP);
        mp_clear(o,&o->hinv); mp_clear(o,&o->x);
        mp_clear(o,o->jump); mp_clear(o,o->pending);
        bridge_free(o,o->jump); bridge_free(o,o->pending);
        o->jump=o->pending=NULL;
    }
    if (o->meta.e>1) {
        phase(o,T_COPY);
        if (!o->current) { status=ALGEBRA; goto finish; }
        o->base=matrix_copy(o,o->current);
        if (!o->base) { status=(int)o->meta.status; goto finish; }
        for (long bit=bit_length((ulong)o->meta.e)-2;bit>=0;bit--) {
            int digit=(int)(((ulong)o->meta.e>>bit)&1);
            status=do_step(o,o->current,o->current,2,digit); if (status) goto finish;
            if (digit) {
                status=do_step(o,o->base,o->current,3,1); if (status) goto finish;
            }
        }
    }
    if (o->meta.steps_done!=o->meta.expected_products) { status=ALGEBRA; goto finish; }
    if (!o->meta.expected_products) {
        phase(o,T_COMBINE_REDUCE);
        for (long j=0;j<o->meta.L;j++) for (int part=0;part<2;part++)
            fmpz_mod_poly_add(o->trace[part].p,o->trace[part].p,
                o->A->parts+2*(j*o->meta.L+j)+part,o->ctx);
    }
finish:
    if (!status) status=(int)o->meta.status;
    if (status) {
        phase(o,T_CLEANUP);
        matrix_clear(o,&o->temporary); matrix_clear(o,&o->shifted);
        if (o->pending) mp_clear(o,o->pending);
    }
    o->meta.status=status; o->meta.state=status ? FAILED : DONE;
    api_end(o);
    return status;
}

static void owner_cleanup(owner *o) {
    phase(o,T_CLEANUP);
    snapshot_clear(o,&o->initial);
    if (o->step) {
        for (long j=0;j<o->meta.expected_products;j++) {
            snapshot_clear(o,&o->step[j].left); snapshot_clear(o,&o->step[j].right);
            snapshot_clear(o,&o->step[j].product); snapshot_clear(o,&o->step[j].before);
            snapshot_clear(o,&o->step[j].after);
        }
        bridge_free(o,o->step); o->step=NULL;
    }
    matrix_clear(o,&o->temporary); matrix_clear(o,&o->shifted);
    matrix_clear(o,&o->current); matrix_clear(o,&o->base); matrix_clear(o,&o->A);
    if (o->jump) { mp_clear(o,o->jump); bridge_free(o,o->jump); o->jump=NULL; }
    if (o->pending) { mp_clear(o,o->pending); bridge_free(o,o->pending); o->pending=NULL; }
    mp_clear(o,&o->hinv); mp_clear(o,&o->x);
    mp_clear(o,&o->trace[0]); mp_clear(o,&o->trace[1]); mp_clear(o,&o->g);
    if (o->scalar_live) { fmpz_clear(o->scalar); o->scalar_live=0; }
    if (o->context_live) {
        fmpz_mod_ctx_clear(o->ctx); o->context_live=0; add_count(o,C_ctx_clears,1);
    }
    bridge_free(o,o->limb_buffer); o->limb_buffer=NULL;
}

int qcb_resident_create(long N, long L, long d, long P, int audit,
    const unsigned char *g_wire, size_t g_bytes,
    const unsigned char *A_wire, size_t A_bytes, void **owner_out) {
    owner *o=NULL;
    size_t debug_limit=0,offset=0;
    double entered=now_seconds();
    int status;
    if (!owner_out) return SHAPE;
    *owner_out=NULL;
    status=dimensions(N,L,d,P,audit,g_bytes,A_bytes,&debug_limit);
    if (status) return status;
    if (!g_wire || !A_wire) return SHAPE;
    if (!matching_abi()) return ABI;
    o=calloc(1,sizeof(*o)); if (!o) return ALLOC;
    o->active=T_CONTROL; o->api_started=o->phase_started=entered;
    o->counts[C_bridge_allocs]=1;
    o->meta=(metadata){1,READY,0,N,L,d,P,N/d,audit,products(d)+products(N/d),0,0};
    o->H=P/8+(P%8!=0); o->limbs=P/64+(P%64!=0); o->debug_limit=debug_limit;
    o->receipt=calloc(1,sizeof(*o->receipt));
    if (!o->receipt) { status=ALLOC; goto finish; }
    o->limb_buffer=bridge_alloc(o,(size_t)o->limbs,sizeof(ulong));
    if (!o->limb_buffer) { status=(int)o->meta.status; goto finish; }
    if (o->meta.expected_products) {
        o->step=bridge_alloc(o,(size_t)o->meta.expected_products,sizeof(step_record));
        if (!o->step) { status=(int)o->meta.status; goto finish; }
    }
    phase(o,T_CONTEXT);
    fmpz_init(o->scalar); o->scalar_live=1;
    fmpz_one(o->scalar); fmpz_mul_2exp(o->scalar,o->scalar,(ulong)P);
    fmpz_mod_ctx_init(o->ctx,o->scalar); o->context_live=1; add_count(o,C_ctx_inits,1);
    mp_init(o,&o->g); mp_init(o,&o->trace[0]); mp_init(o,&o->trace[1]);
    o->A=matrix_new(o,L,L); if (!o->A) { status=(int)o->meta.status; goto finish; }
    phase(o,T_PARSE);
    status=fit_mod(o,o->g.p,d+1); if (status) goto finish;
    for (long k=0;k<=d;k++) {
        status=parse_field(o,g_wire+offset,o->scalar); if (status) goto finish;
        if (k==d && !fmpz_is_one(o->scalar)) { status=CANONICAL; goto finish; }
        fmpz_mod_poly_set_coeff_fmpz(o->g.p,k,o->scalar,o->ctx); offset+=(size_t)o->H;
    }
    if (offset!=g_bytes || fmpz_mod_poly_degree(o->g.p,o->ctx)!=d
        || !fmpz_mod_poly_is_monic(o->g.p,o->ctx)) { status=CANONICAL; goto finish; }
    offset=0;
    for (long j=0;j<2*L*L;j++) {
        status=fit_mod(o,o->A->parts+j,d); if (status) goto finish;
        for (long k=0;k<d;k++) {
            status=parse_field(o,A_wire+offset,o->scalar); if (status) goto finish;
            fmpz_mod_poly_set_coeff_fmpz(o->A->parts+j,k,o->scalar,o->ctx);
            offset+=(size_t)o->H;
        }
    }
    if (offset!=A_bytes) { status=SHAPE; goto finish; }
    phase(o,T_COPY);
    o->current=matrix_copy(o,o->A);
    if (!o->current) { status=(int)o->meta.status; goto finish; }
    if (audit) {
        phase(o,T_DEBUG);
        status=snapshot_matrix(o,&o->initial,o->A); if (status) goto finish;
    }
    status=add_count(o,C_creates,1);
finish:
    if (!status) status=(int)o->meta.status;
    if (status) {
        o->meta.status=status;
        owner_cleanup(o); free(o->receipt);
        add_count(o,C_bridge_frees,1); api_end(o); free(o);
    } else {
        api_end(o); *owner_out=o;
    }
    return status;
}

static long meta_value(const metadata *m, const char *key) {
    if (!m || !key) return -1;
#define META(name) if (!strcmp(key,#name)) return m->name
    META(abi); META(state); META(status); META(N); META(L); META(d); META(P);
    META(e); META(audit); META(expected_products); META(steps_done); META(hinv_created);
#undef META
    return -1;
}
long qcb_resident_meta(const void *opaque, const char *key) {
    const owner *o=opaque; return o ? meta_value(&o->meta,key) : -1;
}
static long count_value(const long *counts, const char *key) {
    if (!counts || !key) return -1;
    for (int j=0;j<C_COUNT;j++) if (!strcmp(key,count_names[j])) return counts[j];
    return -1;
}
static double time_value(const double *seconds, double total, const char *key) {
    if (!seconds || !key) return -1;
    if (!strcmp(key,"total")) return total;
    for (int j=0;j<T_COUNT;j++) if (!strcmp(key,time_names[j])) return seconds[j];
    return -1;
}
long qcb_resident_count(const void *opaque, const char *key) {
    const owner *o=opaque; return o ? count_value(o->counts,key) : -1;
}
double qcb_resident_seconds(const void *opaque, const char *key) {
    const owner *o=opaque; return o ? time_value(o->seconds,o->total,key) : -1;
}
long qcb_resident_step_meta(const void *opaque, long index, const char *key) {
    const owner *o=opaque;
    const step_record *s;
    /* FAILED retains only committed metadata, never a half-step or export. */
    if (!o || !key || index<0 || index>=o->meta.steps_done)
        return -1;
    s=o->step+index;
#define STEP(name) if (!strcmp(key,#name)) return s->name
    STEP(kind); STEP(phase); STEP(digit); STEP(terminal); STEP(rows); STEP(inner);
    STEP(cols); STEP(branch); STEP(native_products); STEP(mod_g_calls);
    STEP(ks_max_bits); STEP(ks_max_span_bits);
#undef STEP
    return -1;
}
static void result_metadata(owner *o, result *r) {
    r->meta=o->meta;
    memcpy(r->counts,o->counts,sizeof(r->counts));
    memcpy(r->seconds,o->seconds,sizeof(r->seconds)); r->total=o->total;
}
void qcb_resident_result_free(void *opaque) {
    result *r=opaque;
    if (r) { free(r->bytes); free(r); }
}
const unsigned char *qcb_resident_result_data(const void *opaque) {
    const result *r=opaque; return r ? r->bytes : NULL;
}
size_t qcb_resident_result_size(const void *opaque) {
    const result *r=opaque; return r ? r->length : 0;
}
long qcb_resident_result_meta(const void *opaque, const char *key) {
    const result *r=opaque;
    if (!r || !key) return -1;
#define RESULT(name) if (!strcmp(key,#name)) return r->name
    RESULT(kind); RESULT(index); RESULT(rows); RESULT(cols); RESULT(parts);
    RESULT(degree); RESULT(bits); RESULT(fields);
#undef RESULT
    if (!strcmp(key,"bytes")) return (long)r->length;
    return meta_value(&r->meta,key);
}
long qcb_resident_result_count(const void *opaque, const char *key) {
    const result *r=opaque; return r ? count_value(r->counts,key) : -1;
}
double qcb_resident_result_seconds(const void *opaque, const char *key) {
    const result *r=opaque; return r ? time_value(r->seconds,r->total,key) : -1;
}

int qcb_resident_export(void *opaque, int kind, long index, void **result_out) {
    owner *o=opaque;
    result *r=NULL;
    const snapshot *s=NULL;
    int status=OK;
    if (!result_out) return SHAPE;
    *result_out=NULL;
    if (!o) return STATE;
    api_begin(o);
    if (o->meta.state!=DONE) { status=STATE; goto finish; }
    if (kind<TRACE || kind>FINAL_MATRIX || index<0) { status=SHAPE; goto finish; }
    if (kind>=LEFT && kind<=JUMP_AFTER) {
        if (index>=o->meta.steps_done) { status=SHAPE; goto finish; }
    } else if (index) { status=SHAPE; goto finish; }
    if (kind!=TRACE && !o->meta.audit) { status=NOT_AVAILABLE; goto finish; }
    if (kind==INITIAL_A) s=&o->initial;
    else if (kind==FINAL_MATRIX) {
        if (o->meta.expected_products) { status=NOT_AVAILABLE; goto finish; }
        s=&o->initial;
    } else if (kind>=LEFT && kind<=JUMP_AFTER) {
        const step_record *p=o->step+index;
        if (kind==LEFT) s=&p->left;
        else if (kind==RIGHT) s=&p->right;
        else if (kind==PRODUCT) s=&p->product;
        else if (kind==JUMP_BEFORE) s=&p->before;
        else s=&p->after;
        if (!s->bytes) { status=NOT_AVAILABLE; goto finish; }
    }
    phase(o,T_EXPORT);
    r=calloc(1,sizeof(*r)); if (!r) { status=ALLOC; goto finish; }
    r->kind=kind; r->index=index; r->degree=o->meta.d; r->bits=o->meta.P;
    if (s) {
        r->rows=s->rows; r->cols=s->cols; r->parts=s->parts;
        r->fields=s->fields; r->length=s->length;
    } else {
        r->rows=r->cols=1; r->parts=2; r->fields=2*o->meta.d;
        if (!size_mul((size_t)r->fields,(size_t)o->H,&r->length)) {
            status=CAPACITY; goto finish;
        }
    }
    r->bytes=malloc(r->length); if (!r->bytes) { status=ALLOC; goto finish; }
    if (s) memcpy(r->bytes,s->bytes,s->length);
    else {
        size_t offset=0;
        for (int part=0;part<2;part++) for (long k=0;k<o->meta.d;k++) {
            fmpz_mod_poly_get_coeff_fmpz(o->scalar,o->trace[part].p,k,o->ctx);
            status=write_field(o,r->bytes+offset,o->scalar); if (status) goto finish;
            offset+=(size_t)o->H;
        }
    }
    status=add_count(o,C_exports,1); if (status) goto finish;
    status=add_count(o,C_wire_output_fields,r->fields); if (status) goto finish;
    status=add_count(o,C_wire_output_bytes,(long)r->length);
finish:
    if (status) {
        phase(o,T_CLEANUP); qcb_resident_result_free(r); r=NULL;
    }
    api_end(o);
    if (r) { result_metadata(o,r); *result_out=r; }
    return status;
}

int qcb_resident_destroy(void **opaque, void **receipt_out) {
    owner *o;
    result *r;
    if (receipt_out) *receipt_out=NULL;
    if (!opaque || !*opaque) return OK;
    o=*opaque; api_begin(o);
    r=o->receipt; o->receipt=NULL;
    owner_cleanup(o);
    /* This last owner allocation is about to be released; the receipt is not
       in bridge_allocs and has an independent public free operation. */
    add_count(o,C_bridge_frees,1);
    api_end(o);
    r->kind=-1; r->index=-1;
    result_metadata(o,r);
    double cleanup_tail=o->phase_started;
    free(o); *opaque=NULL;
    cleanup_tail=now_seconds()-cleanup_tail;
    r->seconds[T_CLEANUP]+=cleanup_tail; r->total+=cleanup_tail;
    if (receipt_out) *receipt_out=r; else qcb_resident_result_free(r);
    return OK;
}
