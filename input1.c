#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t out, const fmpz_mpoly_t in, const fmpz_mpoly_ctx_t ctx);

void construct_s_pair(fmpz_mpoly_t out, const fmpz_mpoly_t p1, const fmpz_mpoly_t p2, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t inout, const fmpz_mpoly_vec_t vec, int lead_reduce, const fmpz_mpoly_ctx_t ctx);

void buchberger_naive(fmpz_mpoly_vec_t out, const fmpz_mpoly_vec_t in, const fmpz_mpoly_ctx_t ctx);

void buchberger_reduced(fmpz_mpoly_vec_t out, const fmpz_mpoly_vec_t in, const fmpz_mpoly_ctx_t ctx);

void reduce_by_matching_term(fmpz_mpoly_t inout, const fmpz_mpoly_t match_term, const fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx);

void reduce_by_matching_term(fmpz_mpoly_t inout, const fmpz_mpoly_t match_term, const fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t gcd, lcm, reduced_p;
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(lcm, ctx);
    fmpz_mpoly_init(reduced_p, ctx);

    // Log the matching term, the polynomial to be reduced, and the polynomial to reduce with
    const char* varnames[] = {"x", "y", "z"};
    char* match_term_str = fmpz_mpoly_get_str_pretty(match_term, varnames, ctx);
    char* inout_str = fmpz_mpoly_get_str_pretty(inout, varnames, ctx);
    char* p_str = fmpz_mpoly_get_str_pretty(p, varnames, ctx);
    flint_fprintf(stderr, "Matching term: %s\n", match_term_str);
    flint_fprintf(stderr, "Polynomial to be reduced: %s\n", inout_str);
    flint_fprintf(stderr, "Polynomial to reduce with: %s\n", p_str);
    flint_free(match_term_str);
    flint_free(inout_str);
    flint_free(p_str);

    // Compute the GCD of the matching term and the leading term of the polynomial to reduce with
    fmpz_mpoly_t lead_term_p;
    fmpz_mpoly_init(lead_term_p, ctx);
    fmpz_mpoly_leadterm(lead_term_p, p, ctx);
    fmpz_mpoly_gcd(gcd, match_term, lead_term_p, ctx);

    // Multiply the polynomial to reduce with by the GCD
    fmpz_mpoly_lcm(lcm, match_term, lead_term_p, ctx);
    fmpz_mpoly_mul(reduced_p, p, gcd, ctx);

    // Subtract the result from the polynomial to be reduced
    fmpz_mpoly_sub(inout, inout, reduced_p, ctx);

    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(lcm, ctx);
    fmpz_mpoly_clear(reduced_p, ctx);
    fmpz_mpoly_clear(lead_term_p, ctx);
}
