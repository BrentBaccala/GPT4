#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_redux, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t out_vec, const fmpz_mpoly_vec_t gen_vec, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t out_vec, const fmpz_mpoly_vec_t gen_vec, const fmpz_mpoly_ctx_t ctx);

void reduce_by_matching_term(fmpz_mpoly_t poly_to_reduce, const fmpz_mpoly_t matching_term, const fmpz_mpoly_t reducer_poly, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t gcd_term, temp_poly;
    char *matching_term_str, *poly_to_reduce_str, *reducer_poly_str;
    const char *varnames[] = {"x", "y", "z"};

    fmpz_mpoly_init(gcd_term, ctx);
    fmpz_mpoly_init(temp_poly, ctx);

    matching_term_str = fmpz_mpoly_get_str_pretty(matching_term, varnames, ctx);
    poly_to_reduce_str = fmpz_mpoly_get_str_pretty(poly_to_reduce, varnames, ctx);
    reducer_poly_str = fmpz_mpoly_get_str_pretty(reducer_poly, varnames, ctx);

    flint_fprintf(stderr, "Matching term: %s\n", matching_term_str);
    flint_fprintf(stderr, "Poly to reduce: %s\n", poly_to_reduce_str);
    flint_fprintf(stderr, "Reducer poly: %s\n", reducer_poly_str);

    fmpz_mpoly_leadterm(gcd_term, reducer_poly, ctx);
    fmpz_mpoly_gcd(gcd_term, gcd_term, matching_term, ctx);
    fmpz_mpoly_mul(temp_poly, reducer_poly, gcd_term, ctx);
    fmpz_mpoly_sub(poly_to_reduce, poly_to_reduce, temp_poly, ctx);

    fmpz_mpoly_clear(gcd_term, ctx);
    fmpz_mpoly_clear(temp_poly, ctx);
    flint_free(matching_term_str);
    flint_free(poly_to_reduce_str);
    flint_free(reducer_poly_str);
}
