#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t leadterm, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const int flag, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx);

void reduce_by_matching_term(fmpz_mpoly_t poly, const fmpz_mpoly_t matching_term, const fmpz_mpoly_t reducer, const fmpz_mpoly_ctx_t ctx) {
    const char* vars[] = {"x", "y", "z"};

    char* poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    char* matching_term_str = fmpz_mpoly_get_str_pretty(matching_term, vars, ctx);
    char* reducer_str = fmpz_mpoly_get_str_pretty(reducer, vars, ctx);
    fprintf(stderr, "Reducing polynomial: %s\nMatching term: %s\nReducer: %s\n", poly_str, matching_term_str, reducer_str);
    
    free(poly_str);
    free(matching_term_str);
    free(reducer_str);

    fmpz_mpoly_t lead_reducer, gcd, temp;
    fmpz_mpoly_init(lead_reducer, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);
    
    fmpz_mpoly_leadterm(lead_reducer, reducer, ctx);
    fmpz_mpoly_gcd(gcd, matching_term, lead_reducer, ctx);
    fmpz_mpoly_mul(temp, reducer, matching_term, ctx);
    fmpz_mpoly_divides(temp, temp, gcd, ctx);
    fmpz_mpoly_sub(poly, poly, temp, ctx);

    fmpz_mpoly_clear(lead_reducer, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}

