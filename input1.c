#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t res, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t inout, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t inout, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    slong i, j;
    fmpz_mpoly_t leadterm_inout, leadterm_vec, temp, gcd;
    const char *vars[] = {"x", "y", "z"};

    fmpz_mpoly_init(leadterm_inout, ctx);
    fmpz_mpoly_init(leadterm_vec, ctx);
    fmpz_mpoly_init(temp, ctx);
    fmpz_mpoly_init(gcd, ctx);

    char *inout_str = fmpz_mpoly_get_str_pretty(inout, vars, ctx);
    flint_fprintf(stderr, "Starting reduction: %s\n", inout_str);
    flint_free(inout_str);

    fmpz_mpoly_leadterm(leadterm_inout, inout, ctx);

    for (i = 0; i < vec->length && !fmpz_mpoly_is_zero(inout, ctx); i++) {
        fmpz_mpoly_set(temp, fmpz_mpoly_vec_entry(vec, i), ctx);

        fmpz_mpoly_leadterm(leadterm_vec, temp, ctx);

        if (lead_reduction) {
            while (fmpz_mpoly_divides(gcd, leadterm_inout, leadterm_vec, ctx) && !fmpz_mpoly_is_zero(inout, ctx)) {
                fmpz_mpoly_divexact(temp, inout, gcd, ctx);
                fmpz_mpoly_mul(temp, temp, fmpz_mpoly_vec_entry(vec, i), ctx);
                fmpz_mpoly_sub(inout, inout, temp, ctx);
                fmpz_mpoly_leadterm(leadterm_inout, inout, ctx);
            }
        } else {
            for (j = 0; j < fmpz_mpoly_length(inout, ctx); j++) {
                fmpz_mpoly_get_term_coeff_fmpz(gcd, inout, j, ctx);
                if (fmpz_mpoly_divides(gcd, gcd, leadterm_vec, ctx)) {
                    char *term_str = fmpz_mpoly_get_term_str_pretty(inout, j, vars, ctx);
                    flint_fprintf(stderr, "Matching term: %s\n", term_str);
                    flint_free(term_str);

                    fmpz_mpoly_get_term_monomial(temp, inout, j, ctx);
                    fmpz_mpoly_mul(temp, temp, fmpz_mpoly_vec_entry(vec, i), ctx);
                    fmpz_mpoly_sub(inout, inout, temp, ctx);
                }
            }
        }
    }

    inout_str = fmpz_mpoly_get_str_pretty(inout, vars, ctx);
    flint_fprintf(stderr, "Reduced polynomial: %s\n", inout_str);
    flint_free(inout_str);

    fmpz_mpoly_clear(leadterm_inout, ctx);
    fmpz_mpoly_clear(leadterm_vec, ctx);
    fmpz_mpoly_clear(temp, ctx);
    fmpz_mpoly_clear(gcd, ctx);
}
