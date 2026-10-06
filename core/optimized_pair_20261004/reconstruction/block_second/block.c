/* D13: owned finite block elimination; no Python object is borrowed. */
#include <limits.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_mat.h"

enum { OK, ARGUMENT, CAPACITY, WIRE_LENGTH, NONCANONICAL, ALLOCATION,
       GAMMA, FILTER, HAT_OR_REMAINDER, NONINTEGRAL, NONUNIT, INTERNAL, ABI };
#define COUNTERS(X) \
 X(hat_count) X(rhs_count) X(hat_slots) X(unit_slots) \
 X(normalized_rhs_slots) X(prefix_slots) X(output_slots) \
 X(fall_divisions) X(fall_products) X(unit_inversions) X(ratio_products) \
 X(contexts) X(rhs_contexts) \
 X(block_groups) X(constant_groups) X(short_groups) X(successor_groups) X(terminal_groups) \
 X(block_coefficient_slots) X(block_successor_slots) X(taylor_products) \
 X(short_unit_inverses) X(short_power_products) X(triangular_dot_terms) \
 X(triangular_diagonal_products) X(diagonal_power_products) X(long_basis_products) \
 X(correction_products) X(long_state_products) X(block_divisions) \
 X(b_integrality_slots) X(remainder_checked_slots) X(power2_checked_slots) X(block_shifted_slots) \
 X(matrix_inits) X(matrix_clears) X(matrix_slots_initialized) X(matrix_slots_cleared) \
 X(matrix_coefficient_gets) X(matrix_coefficient_sets) X(matrix_products) \
 X(matrix_input_entries) X(matrix_output_entries) X(audit_next_slots) \
 X(H_filter_slots) X(gamma_valuation_checks) X(h0_rhs_products) X(h1_rhs_products) \
 X(second_input_matrix_slots) X(native_factorizations) \
 X(input_fields_parsed) X(input_bytes_parsed) X(output_fields_written) \
 X(output_bytes_written) X(array_import_calls) X(array_export_calls) \
 X(scalar_coefficient_gets) X(scalar_coefficient_sets) \
 X(integer_poly_get_calls) X(integer_poly_get_slots) \
 X(integer_poly_set_calls) X(integer_poly_set_slots) \
 X(bulk_right_shift_calls) X(bulk_right_shift_input_slots) \
 X(canonical_checked_fields) X(filtered_checked_slots) X(audit_fields_written) \
 X(fmpz_inits) X(fmpz_clears) X(poly_inits) X(poly_clears) \
 X(context_inits) X(context_clears) X(workspace_alloc_calls) \
 X(workspace_free_calls) X(workspace_allocated_bytes) X(workspace_freed_bytes)
#define ENUM_COUNT(s) CT_##s,
enum { COUNTERS(ENUM_COUNT) CT_COUNT };
#undef ENUM_COUNT
#define NAME_COUNT(s) #s,
static const char *const count_names[] = { COUNTERS(NAME_COUNT) };
#undef NAME_COUNT
enum { ST_PARSE, ST_CONTEXTS, ST_RHS, ST_UNITS, ST_RATIOS, ST_PREPARE,
       ST_TAYLOR, ST_SHORT, ST_LONG, ST_RESIDUAL, ST_ASSEMBLY,
       ST_SERIALIZE, ST_CLEANUP, ST_COUNT };
static const char *const stage_names[] = {
    "parse_check", "contexts", "RHS", "units", "ratios", "block_prepare",
    "taylor_matmul", "short_solve", "long_reconstruct", "residual_check_shift",
    "assembly", "serialize", "workspace_cleanup"
};

typedef struct {
    long status, mode, Q, L, n, node, k, next_bits, residual_bits, audit, field_bytes;
    long counts[CT_COUNT];
    double seconds[ST_COUNT], total;
    uint8_t *bytes;
    size_t length;
} block_result;
typedef struct { fmpz_poly_t p; int live; } ipoly;
typedef struct {
    fmpz_mod_poly_t p;
    const fmpz_mod_ctx_struct *ctx;
    slong bits;
    int live;
} mpoly;
typedef struct {
    fmpz_mod_mat_t p;
    const fmpz_mod_ctx_struct *ctx;
    long rows, cols, bits;
    int live;
} matrix;
enum { MAT_D, MAT_C, MAT_J, MAT_E, MAT_b, MAT_B, MAT_LOW, MAT_CORR, MAT_COUNT };
enum { BP_R, BP_INV, BP_Y, BP_POWER, BP_WORK, BP_U, BP_UTEMP,
       BP_BJ, BP_NUM, BP_QUOT, BP_REM, BP_COUNT };
typedef struct { void *pointer; size_t bytes; } allocation;
typedef struct {
    long q, l, r, k[2], jslots, uslots, wslots, aslots, nslots, B, limbs;
    size_t input_fields, output_fields, input_bytes, output_bytes;
} dimensions;
enum { MP_FALL, MP_HAT, MP_QUOTIENT, MP_REMAINDER, MP_LINEAR,
       MP_NUMERATOR, MP_DENOMINATOR, MP_INVERSE, MP_BLOCK, MP_RATIO,
       MP_DIVIDED, MP_NEXT, MP_HM, MP_HMPOWER, MP_S0, MP_S1, MP_COUNT };
typedef struct {
    dimensions d;
    block_result *result;
    int status, mode, audit, fixed_init;
    fmpz_t gamma, temp[7];
    ipoly raw, hats_dummy;
    ipoly *hats, *H, *units, *norms, *snapshots;
    mpoly scratch[MP_COUNT], *states, *p0, *p1, *p0low;
    fmpz_mod_ctx_struct *ctx;
    long ctx_init, sparse[3], sparse_count;
    fmpz_mod_ctx_t rhs_ctx[2];
    int rhs_init;
    fmpz *prefixes, *output, *weights;
    long prefixes_init, output_init, weights_init;
    matrix matrices[MAT_COUNT];
    mpoly block_poly[BP_COUNT], *rho_powers;
    fmpz *next_snapshots;
    long next_init, next_used, state_count;
    int active_stage;
    double active_started;
    long *row_starts;
    uint8_t *written;
    ulong *limb_buffer;
    allocation blocks[24];
    int block_count;
    const uint8_t *input;
    size_t input_length, input_offset, output_offset;
} workspace;

static double now_seconds(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + (double)t.tv_nsec * 1e-9;
}
static int size_add(size_t a, size_t b, size_t *out) {
    if (a > SIZE_MAX-b) return 0;
    *out = a+b; return 1;
}
static int size_mul(size_t a, size_t b, size_t *out) {
    if (b && a > SIZE_MAX/b) return 0;
    *out = a*b; return 1;
}
static int long_product(long a, long b, long *out) {
    if (a < 0 || b < 0 || (b && a > LONG_MAX/b)) return 0;
    *out = a*b; return 1;
}
static int plan_array(size_t count, size_t item, size_t *total) {
    size_t bytes;
    return size_mul(count,item,&bytes) && bytes<=LONG_MAX
        && size_add(*total,bytes,total) && *total<=LONG_MAX;
}
static int add_count(workspace *w, int key, long amount) {
    long *value = &w->result->counts[key];
    if (amount < 0 || *value > LONG_MAX-amount) {
        w->status = CAPACITY; return CAPACITY;
    }
    *value += amount; return OK;
}
#define COUNT(w, key, amount) do { \
    int count_status_ = add_count((w), CT_##key, (amount)); \
    if (count_status_) return count_status_; \
} while (0)

int qcb_block_abi(void) { return 2; }
const char *qcb_block_header_version(void) { return FLINT_VERSION; }
const char *qcb_block_runtime_version(void) { return flint_version; }
long qcb_block_word_size(const char *key) {
    if (!key) return -1;
    if (!strcmp(key,"char_bits")) return CHAR_BIT;
    if (!strcmp(key,"flint_bits")) return FLINT_BITS;
    if (!strcmp(key,"long_bytes")) return sizeof(long);
    if (!strcmp(key,"slong_bytes")) return sizeof(slong);
    if (!strcmp(key,"ulong_bytes")) return sizeof(ulong);
    if (!strcmp(key,"size_t_bytes")) return sizeof(size_t);
    return -1;
}
long qcb_block_result_meta(const void *opaque, const char *key) {
    const block_result *r = opaque;
    if (!r || !key) return -1;
#define META(s) if (!strcmp(key,#s)) return r->s
    META(status); META(mode); META(Q); META(L); META(n); META(node);
    META(k); META(next_bits); META(residual_bits); META(audit); META(field_bytes);
#undef META
    return -1;
}
long qcb_block_result_count(const void *opaque, const char *key) {
    const block_result *r = opaque;
    int i;
    if (!r || !key) return -1;
    for (i=0;i<CT_COUNT;i++) if (!strcmp(key,count_names[i])) return r->counts[i];
    return -1;
}
double qcb_block_result_seconds(const void *opaque, const char *key) {
    const block_result *r = opaque;
    int i;
    if (!r || !key) return -1;
    if (!strcmp(key,"total")) return r->total;
    for (i=0;i<ST_COUNT;i++) if (!strcmp(key,stage_names[i])) return r->seconds[i];
    return -1;
}
size_t qcb_block_result_output_length(const void *opaque) {
    const block_result *r = opaque; return r ? r->length : 0;
}
const uint8_t *qcb_block_result_bytes(const void *opaque) {
    const block_result *r = opaque; return r ? r->bytes : NULL;
}
void qcb_block_result_free(void *opaque) {
    block_result *r = opaque;
    if (r) { free(r->bytes); free(r); }
}
static int matching_abi(void) {
    return CHAR_BIT==8 && FLINT_BITS==64 && sizeof(long)==8 && sizeof(slong)==8
        && sizeof(ulong)==8 && sizeof(size_t)==8
        && !strcmp(FLINT_VERSION,"3.6.0") && !strcmp(flint_version,"3.6.0");
}

/* All sizes used by either layout are checked before reading the first byte. */
static int batch_dimensions(dimensions *d, long q, long l, int audit, int mode) {
    long v, z, half, s, max_count;
    size_t fields, extra, allocation_bound=0;
    if (l<2 || q<1 || l>q/2 || (audit!=0 && audit!=1)) return ARGUMENT;
    if (q>LONG_MAX/16) return CAPACITY;
    memset(d,0,sizeof(*d));
    d->q=q; d->l=l; d->r=(q-2)/(2*l)+2;
    d->B=q/8+(q%8!=0); d->limbs=q/64+(q%64!=0);
    for (int h=0;h<2;h++) {
        long t=q-l*h, a=t/(2*l), b=t%(2*l);
        d->k[h]=l*a+(b<l?b:l);
    }
    if (d->k[0]+d->k[1]!=q) return INTERNAL;
    if (!long_product(l,q,&d->aslots) || !long_product(l,2*q-l,&d->jslots)
        || !long_product(d->r,q,&v)) return CAPACITY;
    half=d->r%2 ? (d->r-1)/2 : d->r/2;
    if (!long_product(half,d->r%2 ? d->r : d->r-1,&z)) return CAPACITY;
    d->uslots=v-z;
    if (!long_product(d->r-1,q-1,&v)
        || !long_product(d->r-1,d->r-2,&z)
        || !long_product(z,l,&z)) return CAPACITY;
    d->wslots=v-z;
    half=q%2 ? (q-1)/2 : q/2;
    if (!long_product(half,q%2 ? q : q-1,&v)) return CAPACITY;
    s=1+(q-1)%l;
    if (v<s || !long_product(l,v-s,&max_count)) return CAPACITY;
    /* Aggregate initialization/conversion counters also have ample finite bounds. */
    if (!long_product(l,q,&v) || !long_product(v,q,&v)
        || !long_product(v,512,&max_count)) return CAPACITY;
    d->nslots=0;
    for (int h=0;h<2;h++) {
        long t=q-l*h, D=(t-1)/(2*l), a, b;
        if (!long_product(D,t,&a) || !long_product(D,D+1,&b)
            || !long_product(b,l,&b) || a<b
            || !long_product(a-b,l,&a) || d->nslots>LONG_MAX-a) return CAPACITY;
        d->nslots+=a;
    }
    if (!size_mul((size_t)d->r,(size_t)q,&fields)
        || !size_add(fields,mode ? (size_t)d->jslots : (size_t)(3*q),&fields)
        || !size_add(fields,1,&fields)) return CAPACITY;
    d->input_fields=fields;
    fields=(size_t)d->aslots;
    if (audit) {
        if (!size_add((size_t)d->jslots,(size_t)d->uslots,&extra)
            || !size_add(extra,(size_t)d->wslots,&extra)
            || !size_add(extra,(size_t)d->aslots,&extra)
            || !size_add(extra,(size_t)d->nslots,&extra)
            || !size_add(fields,extra,&fields)) return CAPACITY;
    }
    d->output_fields=fields;
    if (!size_mul(d->input_fields,(size_t)d->B,&d->input_bytes)
        || !size_mul(fields,(size_t)d->B,&d->output_bytes)
        || d->input_bytes>LONG_MAX || d->output_bytes>LONG_MAX) return CAPACITY;
    /* Bound every explicit workspace array, including its aggregate byte counter,
       before parsing any field. Weights use at most r entries. */
    if (!plan_array((size_t)d->limbs,sizeof(ulong),&allocation_bound)
        || !plan_array((size_t)(2*d->r),sizeof(ipoly),&allocation_bound)
        || !plan_array((size_t)(2*l),sizeof(ipoly),&allocation_bound)
        || !plan_array((size_t)(2*l),sizeof(mpoly),&allocation_bound)
        || !plan_array((size_t)(q+1),sizeof(fmpz_mod_ctx_struct),&allocation_bound)
        || !plan_array((size_t)d->aslots,sizeof(fmpz),&allocation_bound)
        || !plan_array((size_t)d->aslots,sizeof(fmpz),&allocation_bound)
        || !plan_array((size_t)d->r,sizeof(fmpz),&allocation_bound)
        || !plan_array((size_t)d->aslots,sizeof(uint8_t),&allocation_bound)
        || !plan_array((size_t)l,sizeof(long),&allocation_bound)) return CAPACITY;
    if (!mode && (!plan_array(3,sizeof(ipoly),&allocation_bound)
        || !plan_array((size_t)(l+1),sizeof(mpoly),&allocation_bound)
        || !plan_array((size_t)(l+1),sizeof(mpoly),&allocation_bound)
        || !plan_array((size_t)l,sizeof(mpoly),&allocation_bound))) return CAPACITY;
    if (audit && !plan_array((size_t)(d->r-1),sizeof(ipoly),&allocation_bound)) return CAPACITY;
    if (!plan_array((size_t)(l+1),sizeof(mpoly),&allocation_bound)
        || !plan_array((size_t)d->aslots,16*sizeof(fmpz),&allocation_bound)
        || !plan_array((size_t)(q+l),16*sizeof(fmpz *),&allocation_bound)
        || (audit && !plan_array((size_t)d->nslots,sizeof(fmpz),&allocation_bound))) return CAPACITY;
    return OK;
}
static int block_dimensions(dimensions *d, long n, long l, long node) {
    long k, next, v, z;
    size_t bound=0;
    if (n<1 || l<2 || node<0) return ARGUMENT;
    if (n>LONG_MAX/16) return CAPACITY;
    memset(d,0,sizeof(*d));
    k=l<n?l:n; next=l<=(n-1)/2?n-2*l:0;
    d->q=n; d->l=l; d->k[0]=k;
    d->B=n/8+(n%8!=0); d->limbs=n/64+(n%64!=0);
    if (!long_product(l,n,&v) || !long_product(l,k,&d->aslots)
        || !long_product(l,next,&d->nslots)
        || !long_product(v,n,&z) || !long_product(z,512,&z)
        || l==LONG_MAX || !long_product(n,l+1,&z)
        || !long_product(k,k,&z) || !long_product(k,l+1,&z)) return CAPACITY;
    if (!size_add((size_t)v,k>1?(size_t)(n-1):0,&d->input_fields)
        || !size_add((size_t)d->aslots,(size_t)d->nslots,&d->output_fields)
        || !size_mul(d->input_fields,(size_t)d->B,&d->input_bytes)
        || !size_mul(d->output_fields,(size_t)d->B,&d->output_bytes)
        || d->input_bytes>LONG_MAX || d->output_bytes>LONG_MAX) return CAPACITY;
    if (!plan_array((size_t)d->limbs,sizeof(ulong),&bound)
        || !plan_array((size_t)(n+1),sizeof(fmpz_mod_ctx_struct),&bound)
        || !plan_array((size_t)l,sizeof(ipoly),&bound)
        || !plan_array((size_t)l,sizeof(mpoly),&bound)
        || !plan_array((size_t)d->aslots,sizeof(fmpz),&bound)) return CAPACITY;
    if (next && !plan_array((size_t)(l+1),sizeof(mpoly),&bound)) return CAPACITY;
    /* Conservative simultaneous matrix storage, including row pointers. */
    if (!plan_array((size_t)v,16*sizeof(fmpz),&bound)
        || !plan_array((size_t)(n+l),16*sizeof(fmpz *),&bound)) return CAPACITY;
    return OK;
}
static void *workspace_alloc(workspace *w, size_t count, size_t item) {
    size_t bytes;
    void *p;
    if (!size_mul(count,item,&bytes) || bytes>LONG_MAX || w->block_count>=24
        || w->result->counts[CT_workspace_allocated_bytes]>LONG_MAX-(long)bytes) {
        w->status=CAPACITY; return NULL;
    }
    if (add_count(w,CT_workspace_alloc_calls,1)) return NULL;
    p=calloc(count,item);
    if (!p) {
        /* The failed allocation acquired zero bytes; free(NULL) is its matched
           no-op release, so actual allocation/free attempts remain auditable. */
        free(p); add_count(w,CT_workspace_free_calls,1);
        w->status=ALLOCATION; return NULL;
    }
    w->blocks[w->block_count++]=(allocation){p,bytes};
    add_count(w,CT_workspace_allocated_bytes,(long)bytes);
    return p;
}
static void ip_init(workspace *w, ipoly *p) {
    fmpz_poly_init(p->p); p->live=1; add_count(w,CT_poly_inits,1);
}
static void ip_clear(workspace *w, ipoly *p) {
    if (p && p->live) {
        fmpz_poly_clear(p->p); p->live=0; add_count(w,CT_poly_clears,1);
    }
}
static void mp_clear(workspace *w, mpoly *p) {
    if (p && p->live) {
        fmpz_mod_poly_clear(p->p,p->ctx); p->live=0;
        add_count(w,CT_poly_clears,1);
    }
}
static void mp_init_at(workspace *w, mpoly *p, const fmpz_mod_ctx_struct *ctx, long bits) {
    mp_clear(w,p);
    fmpz_mod_poly_init(p->p,ctx); p->ctx=ctx; p->bits=bits; p->live=1;
    add_count(w,CT_poly_inits,1);
}
static void mp_init(workspace *w, mpoly *p, long bits) {
    mp_init_at(w,p,w->ctx+bits,bits);
}
static void iget(workspace *w, fmpz_t value, const ipoly *p, long index) {
    fmpz_poly_get_coeff_fmpz(value,p->p,index);
    add_count(w,CT_scalar_coefficient_gets,1);
}
static void iset(workspace *w, ipoly *p, long index, const fmpz_t value) {
    fmpz_poly_set_coeff_fmpz(p->p,index,value);
    add_count(w,CT_scalar_coefficient_sets,1);
}
static void mget(workspace *w, fmpz_t value, const mpoly *p, long index) {
    fmpz_mod_poly_get_coeff_fmpz(value,p->p,index,p->ctx);
    add_count(w,CT_scalar_coefficient_gets,1);
}
static void mset(workspace *w, mpoly *p, long index, const fmpz_t value) {
    fmpz_mod_poly_set_coeff_fmpz(p->p,index,value,p->ctx);
    add_count(w,CT_scalar_coefficient_sets,1);
}
static void integer_get(workspace *w, ipoly *out, const mpoly *in) {
    fmpz_mod_poly_get_fmpz_poly(out->p,in->p,in->ctx);
    add_count(w,CT_integer_poly_get_calls,1);
    add_count(w,CT_integer_poly_get_slots,fmpz_mod_poly_length(in->p,in->ctx));
}
static void integer_set(workspace *w, mpoly *out, const ipoly *in) {
    fmpz_mod_poly_set_fmpz_poly(out->p,in->p,out->ctx);
    add_count(w,CT_integer_poly_set_calls,1);
    add_count(w,CT_integer_poly_set_slots,fmpz_poly_length(in->p));
}
static void bulk_shift(workspace *w, ipoly *p, long shift) {
    add_count(w,CT_bulk_right_shift_calls,1);
    add_count(w,CT_bulk_right_shift_input_slots,fmpz_poly_length(p->p));
    fmpz_poly_scalar_fdiv_2exp(p->p,p->p,(ulong)shift);
}
static int canonical(workspace *w, const fmpz_t x, long bits) {
    COUNT(w,canonical_checked_fields,1);
    return fmpz_sgn(x)<0 || fmpz_bits(x)>(ulong)bits ? NONCANONICAL : OK;
}
static int filtered(workspace *w, const fmpz_t x, long bits, long degree) {
    COUNT(w,filtered_checked_slots,1);
    fmpz_fdiv_r_2exp(w->temp[6],x,(ulong)(degree<bits?degree:bits));
    return fmpz_is_zero(w->temp[6]) ? OK : FILTER;
}
static int parse_field(workspace *w, fmpz_t out, long bits) {
    long b, B=w->d.B;
    if (w->input_offset>w->input_length || (size_t)B>w->input_length-w->input_offset)
        return WIRE_LENGTH;
    memset(w->limb_buffer,0,(size_t)w->d.limbs*sizeof(ulong));
    for (b=0;b<B;b++) {
        uint8_t byte=w->input[w->input_offset+(size_t)b];
        if (b>=(bits+7)/8 && byte) return NONCANONICAL;
        if (bits%8 && b==bits/8 && (byte>>(bits%8))) return NONCANONICAL;
        w->limb_buffer[b/8] |= ((ulong)byte) << (8*(b%8));
    }
    fmpz_set_ui_array(out,w->limb_buffer,w->d.limbs);
    COUNT(w,array_import_calls,1);
    COUNT(w,canonical_checked_fields,1);
    w->input_offset+=(size_t)B;
    COUNT(w,input_fields_parsed,1); COUNT(w,input_bytes_parsed,B);
    return OK;
}
static int parse_poly(workspace *w, ipoly *p, long length, long bits, int slot_key) {
    int status;
    ip_init(w,p);
    for (long j=0;j<length;j++) {
        status=parse_field(w,w->temp[0],bits); if (status) return status;
        iset(w,p,j,w->temp[0]);
        if (slot_key>=0) { status=add_count(w,slot_key,1); if (status) return status; }
    }
    return w->status;
}
static int write_field(workspace *w, const fmpz_t value, long bits, int audit) {
    int status=canonical(w,value,bits);
    long b, B=w->d.B;
    if (status) return status;
    if (w->output_offset>w->d.output_bytes || (size_t)B>w->d.output_bytes-w->output_offset)
        return INTERNAL;
    fmpz_get_ui_array(w->limb_buffer,w->d.limbs,value);
    COUNT(w,array_export_calls,1);
    for (b=0;b<B;b++)
        w->result->bytes[w->output_offset+(size_t)b]=(uint8_t)(w->limb_buffer[b/8]>>(8*(b%8)));
    w->output_offset+=(size_t)B;
    COUNT(w,output_fields_written,1); COUNT(w,output_bytes_written,B);
    if (audit) { COUNT(w,audit_fields_written,1); }
    return w->status;
}
static int write_poly(workspace *w, const ipoly *p, long length, long bits) {
    for (long j=0;j<length;j++) {
        int status;
        iget(w,w->temp[0],p,j);
        status=write_field(w,w->temp[0],bits,1); if (status) return status;
    }
    return w->status;
}

static int fixed_init(workspace *w) {
    fmpz_init(w->gamma);
    for (int j=0;j<7;j++) fmpz_init(w->temp[j]);
    w->fixed_init=1; COUNT(w,fmpz_inits,8);
    ip_init(w,&w->raw);
    w->limb_buffer=workspace_alloc(w,(size_t)w->d.limbs,sizeof(ulong));
    return w->limb_buffer ? OK : w->status;
}
static int init_contexts(workspace *w) {
    long q=w->d.q;
    w->ctx=workspace_alloc(w,(size_t)(q+1),sizeof(fmpz_mod_ctx_struct));
    if (!w->ctx) return w->status;
    for (long b=1;b<=q;b++) {
        fmpz_one(w->temp[0]); fmpz_mul_2exp(w->temp[0],w->temp[0],(ulong)b);
        fmpz_mod_ctx_init(w->ctx+b,w->temp[0]); w->ctx_init=b;
        COUNT(w,context_inits,1); COUNT(w,contexts,1);
    }
    return OK;
}
static int vector_init(workspace *w, fmpz **out, long length, long *initialized) {
    *out=workspace_alloc(w,(size_t)length,sizeof(fmpz));
    if (!*out) return w->status;
    for (long j=0;j<length;j++) {
        fmpz_init(*out+j); (*initialized)++; COUNT(w,fmpz_inits,1);
    }
    return OK;
}
static int batch_storage(workspace *w) {
    long r=w->d.r, l=w->d.l;
    w->state_count=2*l;
    w->hats=workspace_alloc(w,(size_t)r,sizeof(ipoly));
    w->units=workspace_alloc(w,(size_t)r,sizeof(ipoly));
    w->norms=workspace_alloc(w,(size_t)(2*l),sizeof(ipoly));
    w->states=workspace_alloc(w,(size_t)(2*l),sizeof(mpoly));
    if (!w->hats || !w->units || !w->norms || !w->states) return w->status;
    if (w->mode==0) {
        w->H=workspace_alloc(w,3,sizeof(ipoly));
        w->p0=workspace_alloc(w,(size_t)(l+1),sizeof(mpoly));
        w->p1=workspace_alloc(w,(size_t)(l+1),sizeof(mpoly));
        w->p0low=workspace_alloc(w,(size_t)l,sizeof(mpoly));
        if (!w->H || !w->p0 || !w->p1 || !w->p0low) return w->status;
    }
    if (w->audit) {
        w->snapshots=workspace_alloc(w,(size_t)(r-1),sizeof(ipoly));
        if (!w->snapshots) return w->status;
    }
    return OK;
}
static int check_gamma(workspace *w) {
    COUNT(w,gamma_valuation_checks,1);
    fmpz_fdiv_r_2exp(w->temp[0],w->gamma,2);
    return fmpz_equal_ui(w->temp[0],2) ? OK : GAMMA;
}
static int parse_batch(workspace *w) {
    int status;
    long q=w->d.q, r=w->d.r, l=w->d.l;
    status=parse_field(w,w->gamma,q); if (status) return status;
    status=check_gamma(w); if (status) return status;
    for (long d=0;d<r;d++) {
        status=parse_poly(w,w->hats+d,q,q,CT_hat_slots); if (status) return status;
        COUNT(w,hat_count,1);
    }
    if (w->mode==0) {
        for (int h=0;h<3;h++) {
            status=parse_poly(w,w->H+h,q,q,-1); if (status) return status;
        }
    } else {
        for (long i=0;i<2*l;i++) {
            long bits=q-l*(i%2);
            status=parse_poly(w,w->norms+i,bits,bits,-1); if (status) return status;
        }
    }
    return w->input_offset==w->input_length ? w->status : WIRE_LENGTH;
}
static void make_linear(workspace *w, mpoly *out, long bits, long node) {
    mp_init(w,out,bits);
    fmpz_set_si(w->temp[0],-node); mset(w,out,0,w->temp[0]);
    fmpz_one(w->temp[0]); mset(w,out,1,w->temp[0]);
}
static void copy_low(workspace *w, mpoly *out, const ipoly *input, long bits) {
    fmpz_poly_set(w->raw.p,input->p);
    fmpz_poly_truncate(w->raw.p,bits);
    integer_set(w,out,&w->raw);
}

/* Identical powers and products to the frozen _normalized_rhs. */
static int construct_rhs(workspace *w) {
    long q=w->d.q, l=w->d.l, low=q-l;
    int status;
    mpoly *hm=&w->scratch[MP_HM], *hm_power=&w->scratch[MP_HMPOWER];
    mpoly *s0=&w->scratch[MP_S0], *s1=&w->scratch[MP_S1];
    for (int h=0;h<3;h++) for (long k=0;k<q;k++) {
        iget(w,w->temp[0],w->H+h,k);
        COUNT(w,H_filter_slots,1);
        status=filtered(w,w->temp[0],q,k+(h==2)); if (status) return status;
    }
    for (int h=0;h<2;h++) {
        long bits=q-l*h;
        fmpz_one(w->temp[0]); fmpz_mul_2exp(w->temp[0],w->temp[0],(ulong)bits);
        fmpz_mod_ctx_init(w->rhs_ctx[h],w->temp[0]); w->rhs_init++;
        COUNT(w,context_inits,1); COUNT(w,rhs_contexts,1);
    }
    for (long i=0;i<=l;i++) {
        mp_init_at(w,w->p0+i,w->rhs_ctx[0],q);
        mp_init_at(w,w->p1+i,w->rhs_ctx[1],low);
    }
    fmpz_one(w->temp[0]); mset(w,w->p0,0,w->temp[0]); mset(w,w->p1,0,w->temp[0]);
    integer_set(w,w->p0+1,w->H+1);
    fmpz_poly_zero(w->raw.p);
    for (long k=0;k<low;k++) {
        iget(w,w->temp[0],w->H+2,k);
        fmpz_fdiv_q_2exp(w->temp[0],w->temp[0],1);
        fmpz_fdiv_r_2exp(w->temp[0],w->temp[0],(ulong)low);
        iset(w,&w->raw,k,w->temp[0]);
    }
    integer_set(w,w->p1+1,&w->raw);
    for (long i=2;i<=l;i++) {
        fmpz_mod_poly_mullow(w->p0[i].p,w->p0[i-1].p,w->p0[1].p,q,w->rhs_ctx[0]);
        COUNT(w,h0_rhs_products,1);
        fmpz_mod_poly_mullow(w->p1[i].p,w->p1[i-1].p,w->p1[1].p,low,w->rhs_ctx[1]);
        COUNT(w,h1_rhs_products,1);
    }
    for (long i=0;i<l;i++) {
        mp_init_at(w,w->p0low+i,w->rhs_ctx[1],low);
        integer_get(w,&w->raw,w->p0+i); fmpz_poly_truncate(w->raw.p,low);
        integer_set(w,w->p0low+i,&w->raw);
    }
    mp_init_at(w,hm,w->rhs_ctx[0],q); integer_set(w,hm,w->H);
    mp_init_at(w,hm_power,w->rhs_ctx[0],q);
    fmpz_one(w->temp[0]); mset(w,hm_power,0,w->temp[0]);
    mp_init_at(w,s0,w->rhs_ctx[0],q); mp_init_at(w,s1,w->rhs_ctx[1],low);
    /* temp[4]=gamma/2 mod 2^low, temp[5] its current power. */
    fmpz_fdiv_q_2exp(w->temp[4],w->gamma,1);
    fmpz_fdiv_r_2exp(w->temp[4],w->temp[4],(ulong)low); fmpz_one(w->temp[5]);
    for (long i=0;i<l;i++) {
        if (i==1) fmpz_mod_poly_set(hm_power->p,hm->p,w->rhs_ctx[0]);
        else if (i>1) {
            fmpz_mod_poly_mullow(s0->p,hm_power->p,hm->p,q,w->rhs_ctx[0]);
            COUNT(w,h0_rhs_products,1);
            fmpz_mod_poly_set(hm_power->p,s0->p,w->rhs_ctx[0]);
        }
        if (!i) {
            fmpz_mod_poly_set(s0->p,w->p0[l].p,w->rhs_ctx[0]);
            fmpz_mod_poly_set(s1->p,w->p1[l].p,w->rhs_ctx[1]);
        } else {
            fmpz_mod_poly_mullow(s0->p,w->p0[l-i].p,hm_power->p,q,w->rhs_ctx[0]);
            COUNT(w,h0_rhs_products,1);
            fmpz_mod_poly_mullow(s1->p,w->p1[l-i].p,w->p0low[i].p,low,w->rhs_ctx[1]);
            COUNT(w,h1_rhs_products,1);
        }
        ip_init(w,w->norms+2*i); integer_get(w,w->norms+2*i,s0);
        ip_init(w,w->norms+2*i+1);
        for (long k=0;k<low;k++) {
            mget(w,w->temp[0],s1,k);
            fmpz_mul(w->temp[0],w->temp[0],w->temp[5]);
            fmpz_fdiv_r_2exp(w->temp[0],w->temp[0],(ulong)low);
            iset(w,w->norms+2*i+1,k,w->temp[0]);
        }
        if (i+1<l) {
            fmpz_mul(w->temp[5],w->temp[5],w->temp[4]);
            fmpz_fdiv_r_2exp(w->temp[5],w->temp[5],(ulong)low);
        }
    }
    return w->status;
}

static int start_states(workspace *w) {
    long q=w->d.q, l=w->d.l;
    for (long c=0;c<2*l;c++) {
        int status;
        long bits=q-l*(c%2);
        if (fmpz_poly_length(w->norms[c].p)>bits) return INTERNAL;
        for (long k=0;k<bits;k++) {
            iget(w,w->temp[0],w->norms+c,k);
            status=canonical(w,w->temp[0],bits); if (status) return status;
            status=filtered(w,w->temp[0],bits,k); if (status) return status;
            COUNT(w,normalized_rhs_slots,1);
        }
        COUNT(w,rhs_count,1);
        mp_init(w,w->states+c,bits); integer_set(w,w->states+c,w->norms+c);
        if (!w->audit) ip_clear(w,w->norms+c);
    }
    return w->status;
}
static int extract_units(workspace *w) {
    long q=w->d.q, r=w->d.r;
    mpoly *fall=&w->scratch[MP_FALL], *hat=&w->scratch[MP_HAT];
    mpoly *quot=&w->scratch[MP_QUOTIENT], *rem=&w->scratch[MP_REMAINDER];
    mpoly *linear=&w->scratch[MP_LINEAR];
    int status;
    if (fmpz_poly_length(w->hats[0].p)!=1) return HAT_OR_REMAINDER;
    iget(w,w->temp[0],w->hats,0);
    if (!fmpz_is_one(w->temp[0])) return HAT_OR_REMAINDER;
    mp_init(w,fall,q); fmpz_one(w->temp[0]); mset(w,fall,0,w->temp[0]);
    mp_init(w,hat,q); mp_init(w,quot,q); mp_init(w,rem,q);
    for (long d=0;d<r;d++) {
        integer_set(w,hat,w->hats+d);
        COUNT(w,fall_divisions,1);
        fmpz_mod_poly_divrem_f(w->temp[1],quot->p,rem->p,hat->p,fall->p,w->ctx+q);
        if (!fmpz_is_one(w->temp[1])) return NONUNIT;
        if (fmpz_mod_poly_length(rem->p,rem->ctx)
            || fmpz_mod_poly_length(quot->p,quot->ctx)>q-d) return HAT_OR_REMAINDER;
        ip_init(w,w->units+d); integer_get(w,w->units+d,quot);
        for (long k=0;k<q-d;k++) {
            iget(w,w->temp[0],w->units+d,k);
            fmpz_fdiv_r_2exp(w->temp[1],w->temp[0],(ulong)d);
            if (!fmpz_is_zero(w->temp[1])) return NONINTEGRAL;
        }
        bulk_shift(w,w->units+d,d);
        for (long k=0;k<q-d;k++) {
            iget(w,w->temp[0],w->units+d,k);
            if (k==0 && !fmpz_is_odd(w->temp[0])) return NONUNIT;
            status=filtered(w,w->temp[0],q-d,k); if (status) return status;
            COUNT(w,unit_slots,1);
        }
        if (d+1<r) {
            make_linear(w,linear,q,d);
            fmpz_mod_poly_mul(fall->p,fall->p,linear->p,w->ctx+q);
            COUNT(w,fall_products,1);
        }
    }
    return w->status;
}

/* One active bucket accounts for every early return from the block code. */
static void stage(workspace *w, int key) {
    double t=now_seconds();
    if (w->active_stage>=0) w->result->seconds[w->active_stage]+=t-w->active_started;
    w->active_stage=key; w->active_started=t;
}
static void mat_clear(workspace *w, matrix *p) {
    if (p->live) {
        fmpz_mod_mat_clear(p->p,p->ctx); p->live=0;
        add_count(w,CT_matrix_clears,1);
        add_count(w,CT_matrix_slots_cleared,p->rows*p->cols);
    }
}
static int mat_init(workspace *w, matrix *p, long rows, long cols, long bits) {
    long slots;
    size_t bytes;
    if (rows<1 || cols<1 || bits<1 || bits>w->d.q
        || !long_product(rows,cols,&slots)
        || !size_mul((size_t)slots,sizeof(fmpz),&bytes) || bytes>LONG_MAX
        || !size_mul((size_t)rows,sizeof(fmpz *),&bytes) || bytes>LONG_MAX) return CAPACITY;
    mat_clear(w,p);
    p->rows=rows; p->cols=cols; p->bits=bits; p->ctx=w->ctx+bits;
    fmpz_mod_mat_init(p->p,rows,cols,p->ctx); p->live=1;
    COUNT(w,matrix_inits,1); COUNT(w,matrix_slots_initialized,slots);
    return w->status;
}
static void mat_get(workspace *w, fmpz_t value, const matrix *p, long row, long col) {
    if (!p->live || row<0 || row>=p->rows || col<0 || col>=p->cols) {
        w->status=INTERNAL; fmpz_zero(value); return;
    }
    fmpz_mod_mat_get_entry(value,p->p,row,col,p->ctx);
    add_count(w,CT_matrix_coefficient_gets,1);
}
static int mat_set(workspace *w, matrix *p, long row, long col, const fmpz_t value) {
    int status;
    if (!p->live || row<0 || row>=p->rows || col<0 || col>=p->cols) return INTERNAL;
    status=canonical(w,value,p->bits); if (status) return status;
    fmpz_mod_mat_set_entry(p->p,row,col,value,p->ctx);
    COUNT(w,matrix_coefficient_sets,1); return w->status;
}
static int mat_product(workspace *w, matrix *out, const matrix *a,
                       const matrix *b, int bucket) {
    if (!out->live || !a->live || !b->live || out==a || out==b
        || out->ctx!=a->ctx || a->ctx!=b->ctx || a->cols!=b->rows
        || out->rows!=a->rows || out->cols!=b->cols) return INTERNAL;
    COUNT(w,matrix_products,1);
    COUNT(w,matrix_input_entries,a->rows*a->cols);
    COUNT(w,matrix_input_entries,b->rows*b->cols);
    COUNT(w,matrix_output_entries,out->rows*out->cols);
    stage(w,bucket);
    fmpz_mod_mat_mul(out->p,a->p,b->p,out->ctx);
    return w->status;
}
static void clear_block_objects(workspace *w) {
    for (int j=MAT_COUNT-1;j>=0;j--) mat_clear(w,&w->matrices[j]);
    for (int j=BP_COUNT-1;j>=0;j--) mp_clear(w,&w->block_poly[j]);
    if (w->rho_powers) for (long j=w->d.l;j>=0;j--) mp_clear(w,w->rho_powers+j);
}
static int block_storage(workspace *w, int successor_possible) {
    if (successor_possible) {
        w->rho_powers=workspace_alloc(w,(size_t)(w->d.l+1),sizeof(mpoly));
        if (!w->rho_powers) return w->status;
    }
    if (w->audit && w->d.nslots) {
        int status=vector_init(w,&w->next_snapshots,w->d.nslots,&w->next_init);
        if (status) return status;
    }
    return w->status;
}

/* Actual output extraction: production keeps b intact; the bare entry returns a. */
static int emit_block_coefficients(workspace *w, long n, long k, long node,
                                    int h, long stride, long first) {
    long l=w->d.l, q=w->d.q;
    matrix *b=&w->matrices[MAT_b];
    stage(w,w->mode==2?ST_SHORT:ST_ASSEMBLY);
    for (long i=0;i<l;i++) for (long j=0;j<k;j++) {
        int status;
        if (k==1) {
            mget(w,w->temp[0],w->states+stride*i+first,0);
            COUNT(w,b_integrality_slots,1); /* j=0 is an exact, vacuous check. */
        } else mat_get(w,w->temp[0],b,j,i);
        if (w->mode==2) {
            fmpz_fdiv_q_2exp(w->output+i*k+j,w->temp[0],(ulong)j);
            status=canonical(w,w->output+i*k+j,n-j);
        } else {
            long K=l*node+j, offset=i*q+(h?w->d.k[0]:0)+K;
            if (K>=w->d.k[h] || offset<0 || offset>=w->d.aslots) return INTERNAL;
            fmpz_mul_2exp(w->prefixes+offset,w->temp[0],(ulong)(l*node+l*h));
            status=canonical(w,w->prefixes+offset,q-l*node);
            if (!status) status=filtered(w,w->prefixes+offset,q-l*node,K+l*h);
            if (!status) status=add_count(w,CT_prefix_slots,1);
        }
        if (status) return status;
        COUNT(w,block_coefficient_slots,1);
    }
    return w->status;
}

/* The only block mathematics implementation, shared by all three entrypoints. */
static int block_core(workspace *w, long n, long node, int h,
                      long stride, long first, const mpoly *ratio) {
    long l=w->d.l, k=l<n?l:n, next=l<=(n-1)/2?n-2*l:0;
    long m=next?n-l:0;
    int status;
    matrix *D=&w->matrices[MAT_D], *C=&w->matrices[MAT_C], *J=&w->matrices[MAT_J];
    matrix *E=&w->matrices[MAT_E], *b=&w->matrices[MAT_b];
    matrix *B=&w->matrices[MAT_B], *low=&w->matrices[MAT_LOW], *corr=&w->matrices[MAT_CORR];
    mpoly *R=&w->block_poly[BP_R], *inverse=&w->block_poly[BP_INV];
    mpoly *Y=&w->block_poly[BP_Y], *power=&w->block_poly[BP_POWER];
    mpoly *work=&w->block_poly[BP_WORK], *U=&w->block_poly[BP_U];
    mpoly *utemp=&w->block_poly[BP_UTEMP], *bj=&w->block_poly[BP_BJ];
    mpoly *num=&w->block_poly[BP_NUM], *quot=&w->block_poly[BP_QUOT], *rem=&w->block_poly[BP_REM];
    stage(w,ST_PREPARE);
    COUNT(w,block_groups,1);
    if (next) { COUNT(w,successor_groups,1); } else { COUNT(w,terminal_groups,1); }
    for (long i=0;i<l;i++) {
        const mpoly *state=w->states+stride*i+first;
        if (!state->live || state->bits!=n || state->ctx!=w->ctx+n
            || fmpz_mod_poly_length(state->p,state->ctx)>n) return INTERNAL;
    }
    if (k==1) {
        COUNT(w,constant_groups,1);
        status=emit_block_coefficients(w,n,k,node,h,stride,first);
        stage(w,-1); return status;
    }
    COUNT(w,short_groups,1);
    if (!ratio || !ratio->live || ratio->bits<n-1 || ratio->ctx!=w->ctx+ratio->bits)
        return INTERNAL;
    status=mat_init(w,D,k,n,n); if (status) return status;
    status=mat_init(w,C,n,l+1,n); if (status) return status;
    status=mat_init(w,J,k,l+1,n); if (status) return status;
    fmpz_one(w->temp[0]); status=mat_set(w,D,0,0,w->temp[0]); if (status) return status;
    for (long j=0;j+1<n;j++) for (long s=0;s<k;s++) {
        mat_get(w,w->temp[0],D,s,j);
        fmpz_mul_ui(w->temp[0],w->temp[0],(ulong)node);
        if (s) { mat_get(w,w->temp[1],D,s-1,j); fmpz_add(w->temp[0],w->temp[0],w->temp[1]); }
        fmpz_fdiv_r_2exp(w->temp[0],w->temp[0],(ulong)n);
        status=mat_set(w,D,s,j+1,w->temp[0]); if (status) return status;
    }
    for (long s=0;s<n;s++) {
        if (s<n-1) mget(w,w->temp[0],ratio,s); else fmpz_zero(w->temp[0]);
        fmpz_fdiv_r_2exp(w->temp[0],w->temp[0],(ulong)(n-1));
        status=mat_set(w,C,s,0,w->temp[0]); if (status) return status;
        for (long i=0;i<l;i++) {
            mget(w,w->temp[0],w->states+stride*i+first,s);
            status=mat_set(w,C,s,i+1,w->temp[0]); if (status) return status;
        }
    }
    COUNT(w,taylor_products,1);
    status=mat_product(w,J,D,C,ST_TAYLOR); if (status) return status;
    stage(w,ST_SHORT);
    status=mat_init(w,E,k,k,n); if (status) return status;
    status=mat_init(w,b,k,l,n); if (status) return status;
    mp_init(w,R,n); mp_init(w,inverse,n); mp_init(w,Y,n);
    mp_init(w,power,n); mp_init(w,work,n);
    for (long s=0;s<k;s++) {
        mat_get(w,w->temp[0],J,s,0); mset(w,R,s,w->temp[0]);
    }
    mget(w,w->temp[4],R,0);
    if (!fmpz_is_odd(w->temp[4])) return NONUNIT;
    COUNT(w,short_unit_inverses,1);
    fmpz_mod_poly_inv_series_f(w->temp[1],inverse->p,R->p,k,R->ctx);
    if (!fmpz_is_one(w->temp[1])) return NONUNIT;
    for (long s=1;s<k;s++) { mget(w,w->temp[0],inverse,s-1); mset(w,Y,s,w->temp[0]); }
    fmpz_one(w->temp[0]); status=mat_set(w,E,0,0,w->temp[0]); if (status) return status;
    fmpz_mod_poly_set(power->p,Y->p,Y->ctx);
    for (long j=1;j<k;j++) {
        if (j>1) {
            fmpz_mod_poly_mullow(work->p,power->p,Y->p,k,Y->ctx);
            COUNT(w,short_power_products,1);
            fmpz_mod_poly_set(power->p,work->p,work->ctx);
        }
        for (long s=0;s<k;s++) {
            mget(w,w->temp[0],power,s);
            status=mat_set(w,E,s,j,w->temp[0]); if (status) return status;
        }
    }
    fmpz_one(w->temp[5]); /* R(0)^s, the inverse of E[s,s]. */
    for (long s=0;s<k;s++) {
        mat_get(w,w->temp[0],E,s,s);
        if (!fmpz_is_odd(w->temp[0])) return NONUNIT;
        if (s) {
            fmpz_mul(w->temp[5],w->temp[5],w->temp[4]);
            fmpz_fdiv_r_2exp(w->temp[5],w->temp[5],(ulong)n);
            COUNT(w,diagonal_power_products,1);
        }
        for (long i=0;i<l;i++) {
            mat_get(w,w->temp[2],J,s,i+1);
            for (long j=0;j<s;j++) {
                mat_get(w,w->temp[0],E,s,j); mat_get(w,w->temp[1],b,j,i);
                fmpz_mul(w->temp[0],w->temp[0],w->temp[1]);
                fmpz_sub(w->temp[2],w->temp[2],w->temp[0]);
                COUNT(w,triangular_dot_terms,1);
            }
            fmpz_mul(w->temp[2],w->temp[2],w->temp[5]);
            fmpz_fdiv_r_2exp(w->temp[2],w->temp[2],(ulong)n);
            COUNT(w,triangular_diagonal_products,1);
            COUNT(w,b_integrality_slots,1);
            fmpz_fdiv_r_2exp(w->temp[0],w->temp[2],(ulong)s);
            if (!fmpz_is_zero(w->temp[0])) return NONINTEGRAL;
            status=mat_set(w,b,s,i,w->temp[2]); if (status) return status;
        }
    }
    status=emit_block_coefficients(w,n,k,node,h,stride,first); if (status) return status;
    stage(w,ST_PREPARE);
    mat_clear(w,D); mat_clear(w,C); mat_clear(w,J); mat_clear(w,E);
    mp_clear(w,R); mp_clear(w,inverse); mp_clear(w,Y); mp_clear(w,power); mp_clear(w,work);
    if (next) {
        if (!w->rho_powers || k!=l) return INTERNAL;
        status=mat_init(w,B,m,l,m); if (status) return status;
        status=mat_init(w,low,l,l,m); if (status) return status;
        status=mat_init(w,corr,m,l,m); if (status) return status;
        for (long j=1;j<=l;j++) mp_init(w,w->rho_powers+j,m);
        integer_get(w,&w->raw,ratio); fmpz_poly_truncate(w->raw.p,m);
        integer_set(w,w->rho_powers+1,&w->raw);
        for (long j=2;j<=l;j++) {
            fmpz_mod_poly_mullow(w->rho_powers[j].p,w->rho_powers[j-1].p,w->rho_powers[1].p,m,w->ctx+m);
            COUNT(w,long_basis_products,1);
        }
        mp_init(w,U,m); mp_init(w,utemp,m); mp_init(w,bj,m); mp_init(w,work,m);
        fmpz_set_si(w->temp[0],-node); mset(w,U,0,w->temp[0]);
        fmpz_one(w->temp[0]); mset(w,U,1,w->temp[0]);
        fmpz_mod_poly_set(work->p,U->p,U->ctx); /* fixed X-node */
        for (long s=0;s<m;s++) {
            mget(w,w->temp[0],w->rho_powers+l,s);
            status=mat_set(w,B,s,0,w->temp[0]); if (status) return status;
        }
        for (long j=1;j<=l;j++) {
            if (j>1) {
                fmpz_mod_poly_mullow(utemp->p,U->p,work->p,m,U->ctx);
                COUNT(w,long_basis_products,1);
                fmpz_mod_poly_set(U->p,utemp->p,utemp->ctx);
            }
            if (j<l) {
                fmpz_mod_poly_mullow(bj->p,U->p,w->rho_powers[l-j].p,m,U->ctx);
                COUNT(w,long_basis_products,1);
                for (long s=0;s<m;s++) {
                    mget(w,w->temp[0],bj,s);
                    status=mat_set(w,B,s,j,w->temp[0]); if (status) return status;
                }
            }
        }
        for (long j=0;j<l;j++) for (long i=0;i<l;i++) {
            mat_get(w,w->temp[0],b,j,i);
            fmpz_fdiv_r_2exp(w->temp[0],w->temp[0],(ulong)m);
            status=mat_set(w,low,j,i,w->temp[0]); if (status) return status;
        }
        COUNT(w,correction_products,1);
        status=mat_product(w,corr,B,low,ST_LONG); if (status) return status;
        for (long i=0;i<l;i++) {
            mpoly *state=w->states+stride*i+first;
            stage(w,ST_PREPARE);
            mp_init(w,num,m); mp_init(w,quot,m); mp_init(w,rem,m);
            integer_get(w,&w->raw,state); fmpz_poly_truncate(w->raw.p,m);
            integer_set(w,bj,&w->raw); /* disjoint input for rho^L multiplication */
            stage(w,ST_LONG);
            fmpz_mod_poly_mullow(num->p,w->rho_powers[l].p,bj->p,m,num->ctx);
            COUNT(w,long_state_products,1);
            for (long s=0;s<m;s++) {
                mget(w,w->temp[0],num,s); mat_get(w,w->temp[1],corr,s,i);
                fmpz_sub(w->temp[0],w->temp[0],w->temp[1]);
                fmpz_fdiv_r_2exp(w->temp[0],w->temp[0],(ulong)m);
                mset(w,num,s,w->temp[0]);
            }
            stage(w,ST_RESIDUAL);
            COUNT(w,block_divisions,1);
            fmpz_mod_poly_divrem_f(w->temp[1],quot->p,rem->p,num->p,U->p,num->ctx);
            if (!fmpz_is_one(w->temp[1])) return NONUNIT;
            if (fmpz_mod_poly_length(rem->p,rem->ctx)>l
                || fmpz_mod_poly_length(quot->p,quot->ctx)>next) return INTERNAL;
            for (long s=0;s<l;s++) {
                mget(w,w->temp[0],rem,s); COUNT(w,remainder_checked_slots,1);
                if (!fmpz_is_zero(w->temp[0])) return HAT_OR_REMAINDER;
            }
            integer_get(w,&w->raw,quot);
            for (long s=0;s<next;s++) {
                iget(w,w->temp[0],&w->raw,s); COUNT(w,power2_checked_slots,1);
                fmpz_fdiv_r_2exp(w->temp[1],w->temp[0],(ulong)l);
                if (!fmpz_is_zero(w->temp[1])) return NONINTEGRAL;
            }
            bulk_shift(w,&w->raw,l); COUNT(w,block_shifted_slots,next);
            fmpz_poly_truncate(w->raw.p,next);
            for (long s=0;s<next;s++) {
                iget(w,w->temp[0],&w->raw,s);
                status=canonical(w,w->temp[0],next); if (status) return status;
                status=filtered(w,w->temp[0],next,s); if (status) return status;
            }
            mp_init(w,state,next); integer_set(w,state,&w->raw);
            for (long s=0;s<next;s++) {
                COUNT(w,block_successor_slots,1);
                if (w->audit) {
                    if (w->next_used>=w->d.nslots) return INTERNAL;
                    mget(w,w->next_snapshots+w->next_used,state,s); w->next_used++;
                }
            }
        }
    }
    stage(w,ST_PREPARE); clear_block_objects(w); stage(w,-1);
    return w->status;
}

static int solve_blocks(workspace *w) {
    long q=w->d.q, l=w->d.l;
    mpoly *num=&w->scratch[MP_NUMERATOR], *den=&w->scratch[MP_DENOMINATOR];
    mpoly *inverse=&w->scratch[MP_INVERSE], *ratio=&w->scratch[MP_BLOCK];
    int status;
    stage(w,ST_PREPARE);
    status=vector_init(w,&w->prefixes,w->d.aslots,&w->prefixes_init); if (status) return status;
    status=block_storage(w,q>2*l); if (status) return status;
    for (long d=0;d<=(w->d.k[0]-1)/l;d++) {
        long n0=q-2*l*d;
        if (n0>1) {
            long bits=n0-1;
            stage(w,ST_RATIOS);
            if (d+1>=w->d.r || bits>q-d-1) return INTERNAL;
            mp_init(w,num,bits); mp_init(w,den,bits);
            copy_low(w,num,w->units+d,bits); copy_low(w,den,w->units+d+1,bits);
            mget(w,w->temp[0],den,0); if (!fmpz_is_odd(w->temp[0])) return NONUNIT;
            mp_init(w,inverse,bits); mp_init(w,ratio,bits);
            COUNT(w,unit_inversions,1);
            fmpz_mod_poly_inv_series_f(w->temp[1],inverse->p,den->p,bits,den->ctx);
            if (!fmpz_is_one(w->temp[1])) return NONUNIT;
            COUNT(w,ratio_products,1);
            fmpz_mod_poly_mullow(ratio->p,num->p,inverse->p,bits,ratio->ctx);
            if (w->audit) { ip_init(w,w->snapshots+d); integer_get(w,w->snapshots+d,ratio); }
        }
        for (int h=0;h<2;h++) {
            long n=n0-l*h;
            if (n<=0) continue;
            status=block_core(w,n,d,h,2,h,n>1?ratio:NULL); if (status) return status;
        }
    }
    stage(w,-1);
    if (w->audit && w->next_used!=w->d.nslots) return INTERNAL;
    return w->status;
}

static int assemble(workspace *w) {
    long q=w->d.q, l=w->d.l, maxd=(w->d.k[0]-1)/l;
    int status;
    status=vector_init(w,&w->output,w->d.aslots,&w->output_init); if (status) return status;
    status=vector_init(w,&w->weights,maxd+1,&w->weights_init); if (status) return status;
    w->written=workspace_alloc(w,(size_t)w->d.aslots,sizeof(uint8_t));
    if (!w->written) return w->status;
    w->row_starts=workspace_alloc(w,(size_t)l,sizeof(long));
    if (!w->row_starts) return w->status;
    {
        long accumulated=0;
        for (long j=0;j<l;j++) {
            w->row_starts[j]=l*accumulated;
            accumulated+=(q-1-j)/l+1;
        }
        if (accumulated!=q) return INTERNAL;
    }
    fmpz_one(w->temp[3]);
    for (long i=0;i<l;i++) {
        fmpz_mul(w->temp[3],w->temp[3],w->gamma);
        fmpz_fdiv_r_2exp(w->temp[3],w->temp[3],(ulong)q);
    }
    fmpz_one(w->weights);
    for (long d=1;d<=maxd;d++) {
        fmpz_mul(w->weights+d,w->weights+d-1,w->temp[3]);
        fmpz_fdiv_r_2exp(w->weights+d,w->weights+d,(ulong)q);
    }
    for (long i=0;i<l;i++) for (int h=0;h<2;h++) for (long K=0;K<w->d.k[h];K++) {
        long d=K/l, j=K%l, degree=2*d+h;
        long count=(q-1-j)/l+1, offset, prefix=i*q+(h?w->d.k[0]:0)+K;
        if (l*degree+j>=q || degree>=count) return INTERNAL;
        offset=w->row_starts[j]+i*count+degree;
        if (offset<0 || offset>=w->d.aslots || w->written[offset]) return INTERNAL;
        fmpz_mul(w->output+offset,w->weights+d,w->prefixes+prefix);
        fmpz_fdiv_r_2exp(w->output+offset,w->output+offset,(ulong)q);
        status=filtered(w,w->output+offset,q,l*degree+j); if (status) return status;
        w->written[offset]=1; COUNT(w,output_slots,1);
    }
    for (long j=0;j<w->d.aslots;j++) if (!w->written[j]) return INTERNAL;
    return w->status;
}
static int serialize_batch(workspace *w) {
    long q=w->d.q, l=w->d.l;
    int status;
    w->result->bytes=malloc(w->d.output_bytes);
    if (!w->result->bytes) return ALLOCATION;
    for (long j=0;j<w->d.aslots;j++) {
        status=write_field(w,w->output+j,q,0); if (status) return status;
    }
    if (w->audit) {
        for (long c=0;c<2*l;c++) {
            long bits=q-l*(c%2);
            status=write_poly(w,w->norms+c,bits,bits); if (status) return status;
        }
        for (long d=0;d<w->d.r;d++) {
            status=write_poly(w,w->units+d,q-d,q-d); if (status) return status;
        }
        for (long d=0;d<w->d.r-1;d++) {
            long bits=q-2*l*d-1;
            if (!w->snapshots[d].live) return INTERNAL;
            status=write_poly(w,w->snapshots+d,bits,bits); if (status) return status;
        }
        for (long i=0;i<l;i++) for (int h=0;h<2;h++) for (long K=0;K<w->d.k[h];K++) {
            long offset=i*q+(h?w->d.k[0]:0)+K;
            status=write_field(w,w->prefixes+offset,q-l*(K/l),1); if (status) return status;
        }
        long position=0;
        for (long d=0;d<=(w->d.k[0]-1)/l;d++) for (int h=0;h<2;h++) {
            long n=q-l*h-2*l*d;
            if (n<=2*l) continue;
            long next=n-2*l;
            for (long i=0;i<l;i++) for (long s=0;s<next;s++) {
                if (position>=w->next_used) return INTERNAL;
                status=write_field(w,w->next_snapshots+position,next,1); if (status) return status;
                position++; COUNT(w,audit_next_slots,1);
            }
        }
        if (position!=w->d.nslots || position!=w->next_used) return INTERNAL;
    }
    if (w->output_offset!=w->d.output_bytes) return INTERNAL;
    w->result->length=w->output_offset;
    return w->status;
}
static void workspace_clear(workspace *w) {
    long l=w->d.l, r=w->d.r;
    clear_block_objects(w);
    for (long j=0;j<w->next_init;j++) { fmpz_clear(w->next_snapshots+j); add_count(w,CT_fmpz_clears,1); }
    for (int j=MP_COUNT-1;j>=0;j--) mp_clear(w,&w->scratch[j]);
    for (long j=w->state_count;j>0;j--) {
        if (w->states) mp_clear(w,w->states+j-1);
        if (w->norms) ip_clear(w,w->norms+j-1);
    }
    for (long j=l+1;j>0;j--) {
        if (w->p0) mp_clear(w,w->p0+j-1);
        if (w->p1) mp_clear(w,w->p1+j-1);
        if (w->p0low && j<=l) mp_clear(w,w->p0low+j-1);
    }
    for (long j=r;j>0;j--) {
        if (w->hats) ip_clear(w,w->hats+j-1);
        if (w->units) ip_clear(w,w->units+j-1);
        if (w->snapshots && j<r) ip_clear(w,w->snapshots+j-1);
    }
    if (w->H) for (int j=2;j>=0;j--) ip_clear(w,w->H+j);

    ip_clear(w,&w->raw); ip_clear(w,&w->hats_dummy);
    for (long j=0;j<w->prefixes_init;j++) { fmpz_clear(w->prefixes+j); add_count(w,CT_fmpz_clears,1); }
    for (long j=0;j<w->output_init;j++) { fmpz_clear(w->output+j); add_count(w,CT_fmpz_clears,1); }
    for (long j=0;j<w->weights_init;j++) { fmpz_clear(w->weights+j); add_count(w,CT_fmpz_clears,1); }
    if (w->fixed_init) {
        fmpz_clear(w->gamma); for (int j=0;j<7;j++) fmpz_clear(w->temp[j]);
        add_count(w,CT_fmpz_clears,8);
    }
    for (int h=w->rhs_init-1;h>=0;h--) {
        fmpz_mod_ctx_clear(w->rhs_ctx[h]); add_count(w,CT_context_clears,1);
    }
    for (long b=w->ctx_init;b>=1;b--) {
        fmpz_mod_ctx_clear(w->ctx+b); add_count(w,CT_context_clears,1);
    }
    for (long j=w->sparse_count;j>0;j--) {
        fmpz_mod_ctx_clear(w->ctx+w->sparse[j-1]); add_count(w,CT_context_clears,1);
    }
    for (int j=w->block_count-1;j>=0;j--) {
        free(w->blocks[j].pointer);
        add_count(w,CT_workspace_free_calls,1);
        add_count(w,CT_workspace_freed_bytes,(long)w->blocks[j].bytes);
    }
}
static int finish(workspace *w, void **out, int status, double started) {
    stage(w,-1);
    double clean=now_seconds();
    workspace_clear(w);
    w->result->seconds[ST_CLEANUP]+=now_seconds()-clean;
    if (!status) status=w->status;
    if (status) { free(w->result->bytes); w->result->bytes=NULL; w->result->length=0; }
    w->result->status=status;
    w->result->total=now_seconds()-started;
    *out=w->result;
    return status;
}
static int batch_entry(long Q, long L, int audit, const uint8_t *data,
                       size_t length, void **out, int mode) {
    double began=now_seconds(), t;
    workspace w;
    int status;
    if (!out) return ARGUMENT;
    *out=NULL;
    if (!matching_abi()) return ABI;
    memset(&w,0,sizeof(w)); w.active_stage=-1;
    status=batch_dimensions(&w.d,Q,L,audit,mode); if (status) return status;
    if (length!=w.d.input_bytes) return WIRE_LENGTH;
    if (!data) return ARGUMENT;
    w.result=calloc(1,sizeof(block_result)); if (!w.result) return ALLOCATION;
    w.result->mode=mode; w.result->Q=Q; w.result->L=L;
    w.result->audit=audit; w.result->field_bytes=w.d.B;
    w.mode=mode; w.audit=audit; w.input=data; w.input_length=length;
    t=now_seconds();
    status=fixed_init(&w);
    if (!status) status=batch_storage(&w);
    if (!status) status=parse_batch(&w);
    w.result->seconds[ST_PARSE]+=now_seconds()-t;
    if (status) return finish(&w,out,status,began);
    t=now_seconds(); status=init_contexts(&w);
    w.result->seconds[ST_CONTEXTS]+=now_seconds()-t;
    if (status) return finish(&w,out,status,began);
    t=now_seconds();
    if (!mode) status=construct_rhs(&w);
    if (!status) status=start_states(&w);
    w.result->seconds[ST_RHS]+=now_seconds()-t;
    if (status) return finish(&w,out,status,began);
    t=now_seconds(); status=extract_units(&w);
    w.result->seconds[ST_UNITS]+=now_seconds()-t;
    if (status) return finish(&w,out,status,began);
    status=solve_blocks(&w);
    if (status) return finish(&w,out,status,began);
    t=now_seconds(); status=assemble(&w);
    w.result->seconds[ST_ASSEMBLY]+=now_seconds()-t;
    if (status) return finish(&w,out,status,began);
    t=now_seconds(); status=serialize_batch(&w);
    w.result->seconds[ST_SERIALIZE]+=now_seconds()-t;
    return finish(&w,out,status,began);
}
int qcb_block_solve(long Q, long L, int audit, const uint8_t *data,
                     size_t length, void **out) {
    return batch_entry(Q,L,audit,data,length,out,0);
}
int qcb_block_solve_rhs(long Q, long L, int audit, const uint8_t *data,
                         size_t length, void **out) {
    return batch_entry(Q,L,audit,data,length,out,1);
}

int qcb_block_apply(long n, long L, long node, const uint8_t *data,
                    size_t length, void **out) {
    workspace w;
    double began=now_seconds(), t;
    long k, next;
    int status;
    mpoly *ratio;
    if (!out) return ARGUMENT;
    *out=NULL;
    if (!matching_abi()) return ABI;
    memset(&w,0,sizeof(w)); w.active_stage=-1;
    status=block_dimensions(&w.d,n,L,node); if (status) return status;
    if (length!=w.d.input_bytes) return WIRE_LENGTH;
    if (!data) return ARGUMENT;
    k=L<n?L:n; next=L<=(n-1)/2?n-2*L:0;
    w.result=calloc(1,sizeof(block_result)); if (!w.result) return ALLOCATION;
    w.result->mode=2; w.result->L=L; w.result->n=n; w.result->node=node;
    w.result->k=k; w.result->next_bits=next; w.result->residual_bits=next?n-L:0;
    w.result->field_bytes=w.d.B;
    w.mode=2; w.input=data; w.input_length=length; w.state_count=L;
    t=now_seconds(); status=fixed_init(&w);
    if (!status) {
        w.norms=workspace_alloc(&w,(size_t)L,sizeof(ipoly));
        w.states=workspace_alloc(&w,(size_t)L,sizeof(mpoly));
        if (!w.norms || !w.states) status=w.status;
    }
    for (long i=0;!status && i<L;i++) status=parse_poly(&w,w.norms+i,n,n,-1);
    if (!status && k>1) status=parse_poly(&w,&w.hats_dummy,n-1,n-1,-1);
    if (!status && w.input_offset!=length) status=WIRE_LENGTH;
    for (long i=0;!status && i<L;i++) for (long s=0;!status && s<n;s++) {
        iget(&w,w.temp[0],w.norms+i,s);
        status=filtered(&w,w.temp[0],n,s);
    }
    if (!status && k>1) {
        for (long s=0;!status && s<n-1;s++) {
            iget(&w,w.temp[0],&w.hats_dummy,s);
            if (!s && !fmpz_is_odd(w.temp[0])) status=NONUNIT;
            if (!status) status=filtered(&w,w.temp[0],n-1,s);
        }
    }
    w.result->seconds[ST_PARSE]+=now_seconds()-t;
    if (status) return finish(&w,out,status,began);
    t=now_seconds();
    w.ctx=workspace_alloc(&w,(size_t)(n+1),sizeof(fmpz_mod_ctx_struct));
    if (!w.ctx) status=w.status;
    else {
        long precisions[3]={n,next?n-L:0,next};
        for (int j=0;j<3 && precisions[j];j++) {
            long bits=precisions[j];
            fmpz_one(w.temp[0]); fmpz_mul_2exp(w.temp[0],w.temp[0],(ulong)bits);
            fmpz_mod_ctx_init(w.ctx+bits,w.temp[0]); w.sparse[w.sparse_count++]=bits;
            add_count(&w,CT_context_inits,1); add_count(&w,CT_contexts,1);
        }
    }
    w.result->seconds[ST_CONTEXTS]+=now_seconds()-t;
    if (status) return finish(&w,out,status,began);
    stage(&w,ST_PREPARE);
    status=vector_init(&w,&w.output,w.d.aslots,&w.output_init);
    if (!status) status=block_storage(&w,next>0);
    if (!status) {
        for (long i=0;i<L;i++) {
            mp_init(&w,w.states+i,n); integer_set(&w,w.states+i,w.norms+i);
            ip_clear(&w,w.norms+i);
        }
        ratio=&w.scratch[MP_RATIO];
        if (k>1) {
            /* Rho is lifted as an integer polynomial into the n-bit context. */
            mp_init(&w,ratio,n); integer_set(&w,ratio,&w.hats_dummy);
        }
        status=block_core(&w,n,node,0,1,0,k>1?ratio:NULL);
    }
    stage(&w,-1);
    if (status) return finish(&w,out,status,began);
    t=now_seconds();
    w.result->bytes=malloc(w.d.output_bytes);
    if (!w.result->bytes) status=ALLOCATION;
    for (long i=0;!status && i<L;i++) for (long j=0;!status && j<k;j++)
        status=write_field(&w,w.output+i*k+j,n-j,0);
    for (long i=0;!status && i<L;i++) for (long s=0;!status && s<next;s++) {
        mget(&w,w.temp[0],w.states+i,s);
        status=write_field(&w,w.temp[0],next,0);
    }
    if (!status && w.output_offset!=w.d.output_bytes) status=INTERNAL;
    if (!status) w.result->length=w.output_offset;
    w.result->seconds[ST_SERIALIZE]+=now_seconds()-t;
    return finish(&w,out,status,began);
}

