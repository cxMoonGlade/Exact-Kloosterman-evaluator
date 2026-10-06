/* Query-owned public FLINT modular-composition objects; no Python ABI access. */
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
#include "ulong_extras.h"

typedef struct {
    fmpz_mod_ctx_t ctx;
    fmpz_mod_poly_t modulus, inverse;
    slong degree, bits;
    double setup[3]; /* context/parse, inverse acquisition, inverse check */
} batch_context;

typedef struct {
    char **strings;
    long dimensions[5]; /* count, baby width, block count, work slots, bit bound */
    double seconds[4]; /* parse/allocation, compose, serialization, cleanup */
} batch_result;

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

int qcb_batch_abi(void) { return 1; }
const char *qcb_batch_header_version(void) { return FLINT_VERSION; }
const char *qcb_batch_runtime_version(void) { return flint_version; }
size_t qcb_batch_slong_size(void) { return sizeof(slong); }

static int valid_text(const char *text) {
    char *end;
    long count, j;
    const unsigned char *p;
    if (!text || !*text) return 0;
    if (!strcmp(text, "0")) return 1;
    if (!isdigit((unsigned char)*text)) return 0;
    errno = 0;
    count = strtol(text, &end, 10);
    if (errno || count < 1 || end[0] != ' ' || end[1] != ' ') return 0;
    if ((size_t)count > SIZE_MAX / sizeof(fmpz)) return 0;
    p = (const unsigned char *)(end + 2);
    for (j = 0; j < count; j++) {
        if (*p == '-') p++;
        if (!isdigit(*p)) return 0;
        do { p++; } while (isdigit(*p));
        if (j + 1 == count) return *p == '\0';
        if (*p != ' ') return 0;
        p++;
    }
    return 0;
}

static int parse_poly(fmpz_mod_poly_t out, const char *text,
                      const fmpz_mod_ctx_t ctx) {
    fmpz_poly_t integer;
    int ok;
    if (!valid_text(text)) return 0;
    fmpz_poly_init(integer);
    ok = !fmpz_poly_set_str(integer, text);
    if (ok) fmpz_mod_poly_set_fmpz_poly(out, integer, ctx);
    fmpz_poly_clear(integer);
    return ok;
}

void qcb_batch_destroy(void *opaque) {
    batch_context *c = opaque;
    if (c) {
        fmpz_mod_poly_clear(c->modulus, c->ctx);
        fmpz_mod_poly_clear(c->inverse, c->ctx);
        fmpz_mod_ctx_clear(c->ctx);
        free(c);
    }
}

/* Status 1: invalid input; 2: size overflow; 3: parse; 4: inverse;
   5: adapter allocator. FLINT allocation failure may terminate the process. */
int qcb_batch_create(long bits, const char *modulus,
                      const char *inverse_override, void **output) {
    batch_context *c;
    fmpz_t number;
    fmpz_mod_poly_t reversed, check;
    double started = now();
    int status = 0;
    if (!output) return 1;
    *output = NULL;
    if (bits < 1 || bits > LONG_MAX / 4 || !modulus) return 1;
    c = calloc(1, sizeof(batch_context));
    if (!c) return 5;
    fmpz_init(number);
    fmpz_one(number);
    fmpz_mul_2exp(number, number, (ulong)bits);
    fmpz_mod_ctx_init(c->ctx, number);
    fmpz_clear(number);
    fmpz_mod_poly_init(c->modulus, c->ctx);
    fmpz_mod_poly_init(c->inverse, c->ctx);
    fmpz_mod_poly_init(reversed, c->ctx);
    fmpz_mod_poly_init(check, c->ctx);
    c->bits = bits;
    if (!parse_poly(c->modulus, modulus, c->ctx)) { status = 3; goto finish; }
    c->degree = fmpz_mod_poly_degree(c->modulus, c->ctx);
    if (c->degree < 2 || c->degree == LONG_MAX ||
        !fmpz_mod_poly_is_monic(c->modulus, c->ctx)) { status = 1; goto finish; }
    c->setup[0] = now() - started;
    started = now();
    fmpz_mod_poly_reverse(reversed, c->modulus, c->degree + 1, c->ctx);
    if (inverse_override) {
        if (!parse_poly(c->inverse, inverse_override, c->ctx)) { status = 3; goto finish; }
        if (fmpz_mod_poly_degree(c->inverse, c->ctx) > c->degree) { status = 1; goto finish; }
    } else {
        fmpz_mod_poly_inv_series(c->inverse, reversed, c->degree + 1, c->ctx);
    }
    c->setup[1] = now() - started;
    started = now();
    fmpz_mod_poly_mullow(check, reversed, c->inverse, c->degree + 1, c->ctx);
    if (!fmpz_mod_poly_is_one(check, c->ctx)) { status = 4; goto finish; }
    c->setup[2] = now() - started;
finish:
    fmpz_mod_poly_clear(reversed, c->ctx);
    fmpz_mod_poly_clear(check, c->ctx);
    if (status) qcb_batch_destroy(c); else *output = c;
    return status;
}

const double *qcb_batch_setup_seconds(const void *opaque) {
    return ((const batch_context *)opaque)->setup;
}

void qcb_batch_result_free(void *opaque) {
    batch_result *r = opaque;
    if (r) {
        if (r->strings) {
            long j;
            for (j = 0; j < r->dimensions[0]; j++)
                if (r->strings[j]) flint_free(r->strings[j]);
            free(r->strings);
        }
        free(r);
    }
}

static int add_slots(size_t *sum, slong rows, slong cols) {
    size_t product;
    if (rows < 1 || cols < 1 || rows > LONG_MAX / cols ||
        (size_t)rows > SIZE_MAX / (size_t)cols) return 0;
    product = (size_t)rows * (size_t)cols;
    if (product > SIZE_MAX / sizeof(fmpz) ||
        *sum > SIZE_MAX / sizeof(fmpz) - product) return 0;
    *sum += product;
    return *sum <= (size_t)LONG_MAX;
}

int qcb_batch_compose(void *opaque, long count, const char * const *polys,
                       const char *argument, void **output) {
    batch_context *c = opaque;
    batch_result *r = NULL;
    fmpz_mod_poly_struct *inputs = NULL, *results = NULL;
    fmpz_mod_poly_t inner;
    fmpz_t constant;
    slong s, k, j, input_init = 0, result_init = 0, log_s = 0;
    size_t slots = 0;
    ulong log_value;
    double started = now(), cleanup_started;
    int status = 0;
    if (!output) return 1;
    *output = NULL;
    if (!c || count < 1 || !polys || !argument) return 1;
    if (count == LONG_MAX || c->degree > LONG_MAX / count ||
        (size_t)(count + 1) > SIZE_MAX / sizeof(fmpz_mod_poly_struct) ||
        (size_t)count > SIZE_MAX / sizeof(char *)) return 2;
    s = (slong)n_sqrt((ulong)(c->degree * count)) + 1;
    k = (c->degree + 1) / s + 1;
    if (k > LONG_MAX / count ||
        !add_slots(&slots, s, c->degree) ||
        !add_slots(&slots, k * count, s) ||
        !add_slots(&slots, k * count, c->degree)) return 2;
    log_value = (ulong)(s - 1);
    while (log_value) { log_s++; log_value >>= 1; }
    if (c->bits > (LONG_MAX - log_s - 1) / 2) return 2;
    r = calloc(1, sizeof(batch_result));
    if (!r) return 5;
    r->dimensions[0] = count; r->dimensions[1] = s;
    r->dimensions[2] = k; r->dimensions[3] = (long)slots;
    r->dimensions[4] = 2 * c->bits + log_s + 1;
    r->strings = calloc((size_t)count, sizeof(char *));
    inputs = calloc((size_t)count + 1, sizeof(fmpz_mod_poly_struct));
    results = calloc((size_t)count, sizeof(fmpz_mod_poly_struct));
    fmpz_mod_poly_init(inner, c->ctx);
    fmpz_init(constant);
    if (!r->strings || !inputs || !results) { status = 5; goto finish; }
    for (j = 0; j <= count; j++) {
        fmpz_mod_poly_init(inputs + j, c->ctx); input_init++;
    }
    for (j = 0; j < count; j++) {
        fmpz_mod_poly_init(results + j, c->ctx); result_init++;
    }
    if (!parse_poly(inner, argument, c->ctx)) { status = 3; goto finish; }
    if (fmpz_mod_poly_degree(inner, c->ctx) >= c->degree) { status = 1; goto finish; }
    for (j = 0; j < count; j++) {
        if (!parse_poly(inputs + j, polys[j], c->ctx)) { status = 3; goto finish; }
        if (fmpz_mod_poly_degree(inputs + j, c->ctx) >= c->degree) { status = 1; goto finish; }
    }
    r->seconds[0] = now() - started;
    started = now();
    if (fmpz_mod_poly_is_zero(inner, c->ctx)) {
        for (j = 0; j < count; j++) {
            fmpz_mod_poly_get_coeff_fmpz(constant, inputs + j, 0, c->ctx);
            fmpz_mod_poly_set_coeff_fmpz(results + j, 0, constant, c->ctx);
        }
        r->dimensions[1] = r->dimensions[2] = r->dimensions[3] = 0;
    } else {
        /* Separate arrays, and a final zero input, meet the stricter docs. */
        fmpz_mod_poly_compose_mod_brent_kung_vec_preinv(results, inputs,
            count + 1, count, inner, c->modulus, c->inverse, c->ctx);
    }
    r->seconds[1] = now() - started;
    started = now();
    for (j = 0; j < count; j++) {
        r->strings[j] = fmpz_mod_poly_get_str(results + j, c->ctx);
        if (!r->strings[j]) { status = 5; goto finish; }
    }
    r->seconds[2] = now() - started;
finish:
    cleanup_started = now();
    for (j = 0; j < input_init; j++) fmpz_mod_poly_clear(inputs + j, c->ctx);
    for (j = 0; j < result_init; j++) fmpz_mod_poly_clear(results + j, c->ctx);
    free(inputs); free(results);
    fmpz_mod_poly_clear(inner, c->ctx);
    fmpz_clear(constant);
    if (status) qcb_batch_result_free(r);
    else { r->seconds[3] = now() - cleanup_started; *output = r; }
    return status;
}

void * const *qcb_batch_result_strings(const void *opaque) {
    return (void * const *)((const batch_result *)opaque)->strings;
}
const long *qcb_batch_result_dimensions(const void *opaque) {
    return ((const batch_result *)opaque)->dimensions;
}
const double *qcb_batch_result_seconds(const void *opaque) {
    return ((const batch_result *)opaque)->seconds;
}
