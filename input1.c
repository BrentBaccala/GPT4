#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t output, const fmpz_mpoly_t input, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t output, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduce, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t output, fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t output, fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);

void reduce_by_matching_term(fmpz_mpoly_t poly, const fmpz_mpoly_t match_term, const fmpz_mpoly_t reduce_with, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t gcd_term, temp_poly;
    fmpz_mpoly_init(gcd_term, ctx);
    fmpz_mpoly_init(temp_poly, ctx);

    // Log the matching term
    const char *varnames[] = {"x", "y", "z"};
    char *match_term_str = fmpz_mpoly_get_str_pretty(match_term, varnames, ctx);
    flint_fprintf(stderr, "Matching term: %s\n", match_term_str);
    flint_free(match_term_str);

    // Compute GCD of the matching term and the leading term of the polynomial to reduce with
    fmpz_mpoly_t reduce_with_lead;
    fmpz_mpoly_init(reduce_with_lead, ctx);
    fmpz_mpoly_leadterm(reduce_with_lead, reduce_with, ctx);
    fmpz_mpoly_gcd(gcd_term, match_term, reduce_with_lead, ctx);

    // Multiply the polynomial to reduce with by the GCD, store result in temp_poly
    fmpz_mpoly_mul(temp_poly, reduce_with, gcd_term, ctx);

    // Subtract temp_poly from the polynomial to be reduced
    fmpz_mpoly_sub(poly, poly, temp_poly, ctx);

    fmpz_mpoly_clear(gcd_term, ctx);
    fmpz_mpoly_clear(temp_poly, ctx);
    fmpz_mpoly_clear(reduce_with_lead, ctx);
}
