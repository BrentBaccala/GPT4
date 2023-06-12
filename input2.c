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

void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx) {
    // Log the start of the construction
    fprintf(stderr, "Constructing s-pair for:\n");
    char **vars = {"x", "y", "z"};  // adjust according to the number of variables in your polynomial ring
    char *poly1_str = fmpz_mpoly_get_str_pretty(poly1, vars, ctx);
    char *poly2_str = fmpz_mpoly_get_str_pretty(poly2, vars, ctx);
    fprintf(stderr, "Poly1: %s\n", poly1_str);
    fprintf(stderr, "Poly2: %s\n", poly2_str);
    flint_free(poly1_str);
    flint_free(poly2_str);

    // Construct the leading terms of the input polynomials
    fmpz_mpoly_t lt1, lt2;
    fmpz_mpoly_init(lt1, ctx);
    fmpz_mpoly_init(lt2, ctx);
    fmpz_mpoly_leadterm(lt1, poly1, ctx);
    fmpz_mpoly_leadterm(lt2, poly2, ctx);

    // Construct the GCD of the leading terms
    fmpz_mpoly_t gcd;
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_gcd(gcd, lt1, lt2, ctx);

    // Multiply a copy of each input polynomial by the leading term of the other input polynomial
    fmpz_mpoly_t poly1_copy, poly2_copy;
    fmpz_mpoly_init(poly1_copy, ctx);
    fmpz_mpoly_init(poly2_copy, ctx);
    fmpz_mpoly_mul(poly1_copy, poly1, lt2, ctx);
    fmpz_mpoly_mul(poly2_copy, poly2, lt1, ctx);

    // Divide each copy by the GCD
    fmpz_mpoly_div(poly1_copy, poly1_copy, gcd, ctx);
    fmpz_mpoly_div(poly2_copy, poly2_copy, gcd, ctx);

    // Subtract the two copies to get the result
    fmpz_mpoly_sub(s_pair, poly1_copy, poly2_copy, ctx);

    // Log the constructed s-pair
    char *s_pair_str = fmpz_mpoly_get_str_pretty(s_pair, vars, ctx);
    fprintf(stderr, "Constructed s-pair: %s\n", s_pair_str);
    flint_free(s_pair_str);

    // Clean up
    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(poly1_copy, ctx);
    fmpz_mpoly_clear(poly2_copy, ctx);
}
