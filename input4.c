#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t res, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gens, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gens, const fmpz_mpoly_ctx_t ctx);

int main(void) {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_vec_t test_case1, test_case2;
    fmpz_mpoly_vec_t basis1, basis2;
    const char *var[] = {"x", "y", "z"};

    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);

    // Test case 1
    fmpz_mpoly_vec_init(test_case1, 0, ctx);
    fmpz_mpoly_t poly1, poly2;
    fmpz_mpoly_init(poly1, ctx);
    fmpz_mpoly_init(poly2, ctx);
    fmpz_mpoly_set_str_pretty(poly1, "2*x+3*y+4*z-5", var, ctx);
    fmpz_mpoly_set_str_pretty(poly2, "3*x+4*y+5*z-2", var, ctx);
    fmpz_mpoly_vec_append(test_case1, poly1, ctx);
    fmpz_mpoly_vec_append(test_case1, poly2, ctx);

    fmpz_mpoly_vec_init(basis1, 0, ctx);
    buchberger_naive(basis1, test_case1, ctx);

    for (slong i = 0; i < basis1->length; i++) {
        char *poly_str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(basis1, i), var, ctx);
        fprintf(stderr, "%s\n", poly_str);
        flint_free(poly_str);
    }

    // Test case 2
    fmpz_mpoly_vec_init(test_case2, 0, ctx);
    fmpz_mpoly_vec_append(test_case2, poly1, ctx);
    fmpz_mpoly_vec_append(test_case2, poly2, ctx);

    fmpz_mpoly_vec_init(basis2, 0, ctx);
    buchberger_reduced(basis2, test_case2, ctx);

    for (slong i = 0; i < basis2->length; i++) {
        char *poly_str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(basis2, i), var, ctx);
        fprintf(stderr, "%s\n", poly_str);
        flint_free(poly_str);
    }

    fmpz_mpoly_clear(poly1, ctx);
    fmpz_mpoly_clear(poly2, ctx);
    fmpz_mpoly_vec_clear(test_case1, ctx);
    fmpz_mpoly_vec_clear(test_case2, ctx);
    fmpz_mpoly_vec_clear(basis1, ctx);
    fmpz_mpoly_vec_clear(basis2, ctx);
    fmpz_mpoly_ctx_clear(ctx);

    return 0;
}
