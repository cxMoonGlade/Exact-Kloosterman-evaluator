/* D5 finite first transform (4Q/3Q/2Q). All mathematical objects are owned here. */
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

enum { ST_SETUP, ST_PARSE, ST_KERNEL, ST_TREE, ST_WEIGHTS, ST_SOLVE,
       ST_SERIALIZE, ST_CLEANUP, ST_COUNT };
enum { CT_Q, CT_R, CT_RHS, CT_NODES, CT_KERNEL_STEPS, CT_KERNEL_INVERSES,
       CT_UNIT_INVERSES, CT_WEIGHT_INVERSES, CT_TREE_PRODUCTS, CT_CONVOLUTIONS,
       CT_PRODUCTS, CT_ADDITIONS, CT_NODE_VALUES, CT_OUTPUT_SLOTS,
       CT_W_SLOTS, CT_C_SLOTS, CT_MAX_BITS, CT_KERNEL_PRODUCTS, CT_COUNT };

typedef struct {
    char **strings;
    long counts[CT_COUNT];
    double seconds[ST_COUNT];
} first_result;

typedef struct {
    fmpz_mod_poly_t W, C;
    slong left, right, span;
} tree_node;

typedef struct {
    slong q, r, nodes, values_init, nodes_init, term_count;
    int fixed_init;
    fmpz_t mod_q, mod_2, mod_3, mod_now, tmp, tmp2, tmp3, numerator, inverse, coefficient;
    fmpz_poly_t phi, H[3];
    fmpz_mod_ctx_t ctx_q, ctx_2, ctx_3;
    fmpz_mod_poly_t A3, H3, T3, scratch2, output;
    fmpz *A, *unit_inverse, *weight, *odd_factorial;
    slong *valuation, *leaves, *terms;
    tree_node *tree;
} first_workspace;

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

int qcb_first_abi(void) { return 2; }
const char *qcb_first_header_version(void) { return FLINT_VERSION; }
const char *qcb_first_runtime_version(void) { return flint_version; }
size_t qcb_first_slong_size(void) { return sizeof(slong); }

static int long_product(long a, long b, long *out) {
    if (a < 0 || b < 0 || (b && a > LONG_MAX / b)) return 0;
    *out = a * b;
    return 1;
}

static int array_size(long length, size_t item) {
    return length >= 0 && (size_t)length <= SIZE_MAX / item;
}

static slong valuation_ui(ulong x) {
    slong count = 0;
    while (!(x & 1)) { count++; x >>= 1; }
    return count; /* only called for strictly positive x */
}

/* Validate the complete wire before fmpz_poly_set_str, whose own parser is lax.
   Coefficients may carry a minus sign here; canonical checking rejects them.
   A canonical Q-bit coefficient has at most Q decimal digits. */
static int valid_text(const char *text, slong q) {
    const unsigned char *p;
    char *end;
    long length, j, digits;
    if (!text || !*text) return 0;
    if (!strcmp(text, "0")) return 1;
    if (!isdigit((unsigned char)*text)) return 0;
    errno = 0;
    length = strtol(text, &end, 10);
    if (errno || length < 1 || length > q || end[0] != ' ' || end[1] != ' ') return 0;
    if (!array_size(length, sizeof(fmpz))) return 0;
    p = (const unsigned char *)(end + 2);
    for (j = 0; j < length; j++) {
        if (*p == '-') p++;
        if (!isdigit(*p)) return 0;
        digits = 0;
        do {
            if (++digits > q) return 0;
            p++;
        } while (isdigit(*p));
        if (j + 1 == length) return *p == '\0';
        if (*p != ' ') return 0;
        p++;
    }
    return 0;
}

static int parse_canonical(fmpz_poly_t polynomial, const char *text,
                           first_workspace *w) {
    slong j;
    if (!valid_text(text, w->q) || fmpz_poly_set_str(polynomial, text)) return 3;
    if (fmpz_poly_degree(polynomial) >= w->q) return 1;
    for (j = 0; j < fmpz_poly_length(polynomial); j++) {
        fmpz_poly_get_coeff_fmpz(w->tmp, polynomial, j);
        if (fmpz_sgn(w->tmp) < 0 || fmpz_cmp(w->tmp, w->mod_q) >= 0) return 1;
    }
    return 0;
}

static void workspace_clear(first_workspace *w) {
    slong j;
    for (j = 0; j < w->nodes_init; j++) {
        fmpz_mod_poly_clear(w->tree[j].W, w->ctx_2);
        fmpz_mod_poly_clear(w->tree[j].C, w->ctx_2);
    }
    for (j = 0; j < w->values_init; j++) {
        fmpz_clear(w->A + j);
        fmpz_clear(w->unit_inverse + j);
        fmpz_clear(w->weight + j);
        fmpz_clear(w->odd_factorial + j);
    }
    free(w->tree); free(w->A); free(w->unit_inverse); free(w->weight);
    free(w->odd_factorial); free(w->valuation); free(w->leaves); free(w->terms);
    if (w->fixed_init) {
        fmpz_mod_poly_clear(w->A3, w->ctx_3);
        fmpz_mod_poly_clear(w->H3, w->ctx_3);
        fmpz_mod_poly_clear(w->T3, w->ctx_3);
        fmpz_mod_poly_clear(w->scratch2, w->ctx_2);
        fmpz_mod_poly_clear(w->output, w->ctx_q);
        fmpz_mod_ctx_clear(w->ctx_q);
        fmpz_mod_ctx_clear(w->ctx_2);
        fmpz_mod_ctx_clear(w->ctx_3);
        fmpz_poly_clear(w->phi);
        for (j = 0; j < 3; j++) fmpz_poly_clear(w->H[j]);
        fmpz_clear(w->mod_q); fmpz_clear(w->mod_2); fmpz_clear(w->mod_3);
        fmpz_clear(w->mod_now); fmpz_clear(w->tmp); fmpz_clear(w->tmp2);
        fmpz_clear(w->tmp3); fmpz_clear(w->numerator); fmpz_clear(w->inverse);
        fmpz_clear(w->coefficient);
    }
}

static int workspace_init(first_workspace *w, slong q, slong r, slong nodes) {
    slong j;
    w->q = q; w->r = r; w->nodes = nodes;
    fmpz_init(w->mod_q); fmpz_init(w->mod_2); fmpz_init(w->mod_3);
    fmpz_init(w->mod_now); fmpz_init(w->tmp); fmpz_init(w->tmp2);
    fmpz_init(w->tmp3); fmpz_init(w->numerator); fmpz_init(w->inverse);
    fmpz_init(w->coefficient);
    fmpz_one(w->mod_q); fmpz_mul_2exp(w->mod_q, w->mod_q, (ulong)q);
    fmpz_one(w->mod_2); fmpz_mul_2exp(w->mod_2, w->mod_2, (ulong)(2*q));
    fmpz_one(w->mod_3); fmpz_mul_2exp(w->mod_3, w->mod_3, (ulong)(3*q));
    fmpz_mod_ctx_init(w->ctx_q, w->mod_q);
    fmpz_mod_ctx_init(w->ctx_2, w->mod_2);
    fmpz_mod_ctx_init(w->ctx_3, w->mod_3);
    fmpz_poly_init(w->phi);
    for (j = 0; j < 3; j++) fmpz_poly_init(w->H[j]);
    fmpz_mod_poly_init(w->A3, w->ctx_3);
    fmpz_mod_poly_init(w->H3, w->ctx_3);
    fmpz_mod_poly_init(w->T3, w->ctx_3);
    fmpz_mod_poly_init(w->scratch2, w->ctx_2);
    fmpz_mod_poly_init(w->output, w->ctx_q);
    w->fixed_init = 1;
    w->A = calloc((size_t)q, sizeof(fmpz));
    w->unit_inverse = calloc((size_t)q, sizeof(fmpz));
    w->weight = calloc((size_t)q, sizeof(fmpz));
    w->odd_factorial = calloc((size_t)q, sizeof(fmpz));
    w->valuation = calloc((size_t)q, sizeof(slong));
    w->leaves = calloc((size_t)q, sizeof(slong));
    w->terms = calloc((size_t)q, sizeof(slong));
    w->tree = calloc((size_t)nodes, sizeof(tree_node));
    if (!w->A || !w->unit_inverse || !w->weight || !w->odd_factorial ||
        !w->valuation || !w->leaves || !w->terms || !w->tree) return 5;
    for (j = 0; j < q; j++) {
        fmpz_init(w->A+j); fmpz_init(w->unit_inverse+j);
        fmpz_init(w->weight+j); fmpz_init(w->odd_factorial+j);
        w->values_init++;
    }
    for (j = 0; j < nodes; j++) {
        fmpz_mod_poly_init(w->tree[j].W, w->ctx_2);
        fmpz_mod_poly_init(w->tree[j].C, w->ctx_2);
        w->nodes_init++;
    }
    return 0;
}

/* Return 4 for a failed odd inverse, 6 for a nonexact binary division. */
static int checked_shift(fmpz_t out, const fmpz_t value, slong shift, fmpz_t scratch) {
    if (shift < 0) return 6;
    fmpz_fdiv_r_2exp(scratch, value, (ulong)shift);
    if (!fmpz_is_zero(scratch)) return 6;
    fmpz_fdiv_q_2exp(out, value, (ulong)shift);
    return 0;
}

static int checked_inverse(fmpz_t out, const fmpz_t value, const fmpz_t modulus) {
    if (!fmpz_is_odd(value) || !fmpz_invmod(out, value, modulus)) return 4;
    return 0;
}

static int construct_kernel(first_workspace *w, first_result *result) {
    slong m, j, shift, bits = 4*w->q;
    ulong odd;
    int status;
    fmpz_set(w->A, w->mod_q);
    for (m = 1; m < w->q; m++) {
        fmpz_set(w->numerator, w->A + m - 1);
        for (j = 0; j < w->term_count && w->terms[j] <= m; j++) {
            slong d = w->terms[j];
            fmpz_poly_get_coeff_fmpz(w->coefficient, w->phi, d);
            fmpz_addmul(w->numerator, w->coefficient, w->A + m - d);
            result->counts[CT_KERNEL_PRODUCTS]++;
        }
        /* Older A coefficients are known to at least this current precision. */
        fmpz_fdiv_r_2exp(w->numerator, w->numerator, (ulong)bits);
        shift = valuation_ui((ulong)m);
        status = checked_shift(w->numerator, w->numerator, shift, w->tmp);
        if (status) return status;
        bits -= shift;
        w->valuation[m] = w->valuation[m-1] + shift;
        odd = (ulong)m >> shift;
        if (odd != 1) {
            fmpz_set_ui(w->tmp, odd);
            fmpz_one(w->mod_now);
            fmpz_mul_2exp(w->mod_now, w->mod_now, (ulong)bits);
            status = checked_inverse(w->inverse, w->tmp, w->mod_now);
            if (status) return status;
            fmpz_mul(w->numerator, w->numerator, w->inverse);
            result->counts[CT_KERNEL_INVERSES]++;
        }
        fmpz_fdiv_r_2exp(w->A + m, w->numerator, (ulong)bits);
        result->counts[CT_KERNEL_STEPS]++;
    }
    for (m = 0; m < w->q; m++) {
        shift = w->q - w->valuation[m];
        status = checked_shift(w->tmp2, w->A + m, shift, w->tmp);
        if (status) return status;
        status = checked_inverse(w->unit_inverse + m, w->tmp2, w->mod_2);
        if (status) return status;
        result->counts[CT_UNIT_INVERSES]++;
        /* Explicit precision conversion, after extracting the odd unit. */
        fmpz_fdiv_r_2exp(w->A + m, w->A + m, (ulong)(3*w->q));
        fmpz_mod_poly_set_coeff_fmpz(w->A3, m, w->A + m, w->ctx_3);
    }
    return 0;
}

static slong construct_tree(first_workspace *w, slong lo, slong hi,
                            slong *next, first_result *result) {
    slong index = (*next)++, middle;
    tree_node *node = w->tree + index;
    node->span = hi - lo;
    result->counts[CT_W_SLOTS] += node->span + 1;
    result->counts[CT_C_SLOTS] += node->span;
    if (hi - lo == 1) {
        w->leaves[lo] = index;
        node->left = node->right = -1;
        fmpz_mod_poly_set_coeff_si(node->W, 0, -lo, w->ctx_2);
        fmpz_mod_poly_set_coeff_ui(node->W, 1, 1, w->ctx_2);
    } else {
        middle = lo + (hi - lo)/2;
        node->left = construct_tree(w, lo, middle, next, result);
        node->right = construct_tree(w, middle, hi, next, result);
        fmpz_mod_poly_mul(node->W, w->tree[node->left].W,
                          w->tree[node->right].W, w->ctx_2);
        result->counts[CT_TREE_PRODUCTS]++;
    }
    return index;
}

static int construct_weights(first_workspace *w, first_result *result) {
    slong m, shift;
    int status;
    fmpz_one(w->odd_factorial);
    for (m = 1; m < w->q; m++) {
        ulong odd = (ulong)m >> valuation_ui((ulong)m);
        fmpz_mul_ui(w->odd_factorial + m, w->odd_factorial + m - 1, odd);
        fmpz_fdiv_r_2exp(w->odd_factorial + m, w->odd_factorial + m, (ulong)(2*w->q));
    }
    for (m = 0; m < w->q; m++) {
        shift = w->q - w->valuation[m] - w->valuation[w->q-1-m];
        if (shift < 1) return 6;
        fmpz_mul(w->tmp, w->odd_factorial + m, w->odd_factorial + w->q-1-m);
        if ((w->q-1-m) & 1) fmpz_neg(w->tmp, w->tmp);
        fmpz_fdiv_r_2exp(w->tmp, w->tmp, (ulong)(2*w->q));
        status = checked_inverse(w->weight + m, w->tmp, w->mod_2);
        if (status) return status;
        fmpz_mul_2exp(w->weight + m, w->weight + m, (ulong)shift);
        fmpz_fdiv_r_2exp(w->weight + m, w->weight + m, (ulong)(2*w->q));
        result->counts[CT_WEIGHT_INVERSES]++;
    }
    return 0;
}

static void merge_tree(first_workspace *w, slong index, first_result *result) {
    tree_node *node = w->tree + index;
    if (node->left < 0) return;
    merge_tree(w, node->left, result);
    merge_tree(w, node->right, result);
    fmpz_mod_poly_mul(node->C, w->tree[node->left].C,
                      w->tree[node->right].W, w->ctx_2);
    fmpz_mod_poly_mul(w->scratch2, w->tree[node->right].C,
                      w->tree[node->left].W, w->ctx_2);
    fmpz_mod_poly_add(node->C, node->C, w->scratch2, w->ctx_2);
    result->counts[CT_PRODUCTS] += 2;
    result->counts[CT_ADDITIONS]++;
}

static int solve_column(first_workspace *w, slong column, first_result *result) {
    slong m;
    int status;
    if (column >= w->r) {
        fmpz_mod_poly_set_fmpz_poly(w->H3, w->H[column-w->r], w->ctx_3);
        fmpz_mod_poly_mullow(w->T3, w->A3, w->H3, w->q, w->ctx_3);
        result->counts[CT_CONVOLUTIONS]++;
    }
    for (m = 0; m < w->q; m++) {
        if (column < w->r) {
            if (m < column) fmpz_zero(w->tmp2);
            else {
                fmpz_mul_2exp(w->tmp2, w->A + m-column, (ulong)column);
                fmpz_fdiv_r_2exp(w->tmp2, w->tmp2, (ulong)(3*w->q));
            }
        } else {
            fmpz_mod_poly_get_coeff_fmpz(w->tmp2, w->T3, m, w->ctx_3);
        }
        status = checked_shift(w->tmp2, w->tmp2, w->q-w->valuation[m], w->tmp);
        if (status) return status;
        fmpz_mul(w->tmp2, w->tmp2, w->unit_inverse + m);
        fmpz_fdiv_r_2exp(w->tmp2, w->tmp2, (ulong)(2*w->q));
        fmpz_mul(w->tmp2, w->tmp2, w->weight + m);
        fmpz_fdiv_r_2exp(w->tmp2, w->tmp2, (ulong)(2*w->q));
        fmpz_mod_poly_set_fmpz(w->tree[w->leaves[m]].C, w->tmp2, w->ctx_2);
        result->counts[CT_NODE_VALUES]++;
    }
    merge_tree(w, 0, result);
    if (fmpz_mod_poly_degree(w->tree[0].C, w->ctx_2) >= w->q) return 1;
    fmpz_mod_poly_zero(w->output, w->ctx_q);
    for (m = 0; m < w->q; m++) {
        fmpz_mod_poly_get_coeff_fmpz(w->tmp2, w->tree[0].C, m, w->ctx_2);
        status = checked_shift(w->tmp2, w->tmp2, w->q, w->tmp);
        if (status) return status;
        fmpz_fdiv_r_2exp(w->tmp2, w->tmp2, (ulong)w->q);
        fmpz_mod_poly_set_coeff_fmpz(w->output, m, w->tmp2, w->ctx_q);
    }
    return 0;
}

void qcb_first_result_free(void *opaque) {
    first_result *result = opaque;
    long j;
    if (!result) return;
    if (result->strings) {
        for (j = 0; j < result->counts[CT_RHS]; j++)
            if (result->strings[j]) flint_free(result->strings[j]);
        free(result->strings);
    }
    free(result);
}

/* Status: 1 contract, 2 size overflow, 3 wire parse, 4 odd inverse,
   5 adapter allocation, 6 nonexact binary division. FLINT OOM can abort. */
int qcb_first_solve(long Q, long r, const char *phi,
                    const char *const H[3], void **output) {
    first_workspace w;
    first_result *result = NULL;
    long nodes, count, output_slots, product_count, tree_bound, levels = 1;
    ulong shifted;
    slong j, next = 0;
    int status = 0;
    double started = now(), cleanup_started;
    if (!output) return 1;
    *output = NULL;
    if (Q < 2 || r < 1 || r > Q || !phi || !H || !H[0] || !H[1] || !H[2]) return 1;
    if (Q > LONG_MAX/16) return 2;
    nodes = 2*Q-1; count = r+3;
    shifted = (ulong)(Q-1);
    while (shifted) { levels++; shifted >>= 1; }
    if (!long_product(Q, count, &output_slots) ||
        !long_product(2*(Q-1), count, &product_count) ||
        !long_product(2*Q, levels, &tree_bound) ||
        !array_size(Q, sizeof(fmpz)) || !array_size(Q, sizeof(slong)) ||
        !array_size(nodes, sizeof(tree_node)) || !array_size(count, sizeof(char *)) ||
        !array_size(output_slots, sizeof(fmpz)) || !array_size(tree_bound, sizeof(fmpz))) return 2;
    memset(&w, 0, sizeof(w));
    result = calloc(1, sizeof(first_result));
    if (!result) return 5;
    result->counts[CT_Q] = Q; result->counts[CT_R] = r;
    result->counts[CT_RHS] = count; result->counts[CT_NODES] = nodes;
    result->counts[CT_OUTPUT_SLOTS] = output_slots;
    /* Explicit kernel numerator*inverse can use 8Q bits; this also bounds
       3Q polynomial-product accumulators for Q>=2. Not allocator telemetry. */
    result->counts[CT_MAX_BITS] = 8*Q;
    result->strings = calloc((size_t)count, sizeof(char *));
    if (!result->strings) { status = 5; goto finish; }
    status = workspace_init(&w, Q, r, nodes);
    if (status) goto finish;
    result->seconds[ST_SETUP] = now()-started;
    started = now();
    status = parse_canonical(w.phi, phi, &w);
    if (status) goto finish;
    for (j = 0; j < 3; j++) {
        status = parse_canonical(w.H[j], H[j], &w);
        if (status) goto finish;
    }
    fmpz_poly_get_coeff_fmpz(w.tmp, w.phi, 0);
    fmpz_poly_get_coeff_fmpz(w.tmp2, w.phi, 1);
    if (!fmpz_is_zero(w.tmp) || !fmpz_is_one(w.tmp2)) { status = 1; goto finish; }
    for (j = 2; j < Q; j++) {
        fmpz_poly_get_coeff_fmpz(w.tmp, w.phi, j);
        if (!fmpz_is_zero(w.tmp)) {
            if ((j & (j-1)) || fmpz_is_odd(w.tmp)) { status = 1; goto finish; }
            w.terms[w.term_count++] = j;
        }
    }
    result->seconds[ST_PARSE] = now()-started;
    started = now();
    status = construct_kernel(&w, result);
    if (status) goto finish;
    result->seconds[ST_KERNEL] = now()-started;
    started = now();
    construct_tree(&w, 0, Q, &next, result);
    if (next != nodes) { status = 1; goto finish; }
    result->seconds[ST_TREE] = now()-started;
    started = now();
    status = construct_weights(&w, result);
    if (status) goto finish;
    result->seconds[ST_WEIGHTS] = now()-started;
    for (j = 0; j < count; j++) {
        started = now();
        status = solve_column(&w, j, result);
        if (status) goto finish;
        result->seconds[ST_SOLVE] += now()-started;
        started = now();
        result->strings[j] = fmpz_mod_poly_get_str(w.output, w.ctx_q);
        if (!result->strings[j]) { status = 5; goto finish; }
        result->seconds[ST_SERIALIZE] += now()-started;
    }
finish:
    cleanup_started = now();
    workspace_clear(&w);
    result->seconds[ST_CLEANUP] = now()-cleanup_started;
    if (status) qcb_first_result_free(result);
    else *output = result;
    return status;
}

const void *qcb_first_result_strings(const void *opaque) {
    return opaque ? ((const first_result *)opaque)->strings : NULL;
}
const long *qcb_first_result_counts(const void *opaque) {
    return opaque ? ((const first_result *)opaque)->counts : NULL;
}
const double *qcb_first_result_seconds(const void *opaque) {
    return opaque ? ((const first_result *)opaque)->seconds : NULL;
}
