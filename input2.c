#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t f, const fmpz_mpoly_t g, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_get_coeff_vars_ui(lt, f, 0, NULL, 0, ctx);
}

void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t f, const fmpz_mpoly_t g, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lc_f, lc_g, gcd_lcm, tmp_f, tmp_g;
    const char *vars[] = {"x", "y", "z"};

    // Log input polynomials
    fprintf(stderr, "Beginning construction of the s-pair\n");
    fprintf(stderr, "f: %s\n", fmpz_mpoly_get_str_pretty(f, vars, ctx));
    fprintf(stderr, "g: %s\n", fmpz_mpoly_get_str_pretty(g, vars, ctx));

    // Initialize and construct leading terms
    fmpz_mpoly_init(lc_f, ctx);
    fmpz_mpoly_init(lc_g, ctx);
    fmpz_mpoly_init(gcd_lcm, ctx);
    fmpz_mpoly_init(tmp_f, ctx);
    fmpz_mpoly_init(tmp_g, ctx);
    
    fmpz_mpoly_leadterm(lc_f, f, ctx);
    fmpz_mpoly_leadterm(lc_g, g, ctx);

    // Construct GCD of leading terms
    fmpz_mpoly_gcd(gcd_lcm, lc_f, lc_g, ctx);

    // Multiply input polynomials by leading terms of the other
    fmpz_mpoly_mul(tmp_f, f, lc_g, ctx);
    fmpz_mpoly_mul(tmp_g, g, lc_f, ctx);

    // Divide by GCD
    fmpz_mpoly_divexact(tmp_f, tmp_f, gcd_lcm, ctx);
    fmpz_mpoly_divexact(tmp_g, tmp_g, gcd_lcm, ctx);

    // Subtract to get the s-pair
    fmpz_mpoly_sub(s_pair, tmp_f, tmp_g, ctx);

    // Log the result
    fprintf(stderr, "Constructed an s-pair for f and g:\n");
    fprintf(stderr, "s-pair: %s\n", fmpz_mpoly_get_str_pretty(s_pair, vars, ctx));

    // Clear memory
    fmpz_mpoly_clear(lc_f, ctx);
    fmpz_mpoly_clear(lc_g, ctx);
    fmpz_mpoly_clear(gcd_lcm, ctx);
    fmpz_mpoly_clear(tmp_f, ctx);
    fmpz_mpoly_clear(tmp_g, ctx);
}
