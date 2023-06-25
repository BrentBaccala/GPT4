#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t p1, const fmpz_mpoly_t p2, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t p, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);

void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);

void buchberger_reduced(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);

void reduce_by_matching_term(fmpz_mpoly_t poly_to_reduce, const fmpz_mpoly_t matching_term, const fmpz_mpoly_t poly_to_reduce_with, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t gcd, temp1, temp2;
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp1, ctx);
    fmpz_mpoly_init(temp2, ctx);

    const char *var[] = {"x", "y", "z"};
    char *matching_term_str = fmpz_mpoly_get_str_pretty(matching_term, var, ctx);
    char *poly_to_reduce_str = fmpz_mpoly_get_str_pretty(poly_to_reduce, var, ctx);
    char *poly_to_reduce_with_str = fmpz_mpoly_get_str_pretty(poly_to_reduce_with, var, ctx);

    flint_fprintf(stderr, "Matching term: %s\n", matching_term_str);
    flint_fprintf(stderr, "Polynomial to reduce: %s\n", poly_to_reduce_str);
    flint_fprintf(stderr, "Polynomial to reduce with: %s\n", poly_to_reduce_with_str);

    fmpz_mpoly_leadterm(temp1, poly_to_reduce_with, ctx);
    fmpz_mpoly_gcd(gcd, matching_term, temp1, ctx);
    fmpz_mpoly_mul(temp1, poly_to_reduce_with, matching_term, ctx);
    fmpz_mpoly_div(temp2, temp1, gcd, ctx);
    fmpz_mpoly_sub(poly_to_reduce, poly_to_reduce, temp2, ctx);

    flint_free(matching_term_str);
    flint_free(poly_to_reduce_str);
    flint_free(poly_to_reduce_with_str);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}
