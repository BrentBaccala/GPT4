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

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t leadterm, gcd_leadterm, temp_poly;
    int found_divisor;
    slong i, nvars;
    char **varnames;

    nvars = ctx->minfo->nvars;
    varnames = (char **) flint_malloc(nvars * sizeof(char *));
    for (i = 0; i < nvars; i++) {
        varnames[i] = (char *) flint_malloc(10 * sizeof(char));
        sprintf(varnames[i], "x%ld", i + 1);
    }

    fprintf(stderr, "Input: %s\n", fmpz_mpoly_get_str_pretty(poly, (const char **) varnames, ctx));

    fmpz_mpoly_init(leadterm, ctx);
    fmpz_mpoly_init(gcd_leadterm, ctx);
    fmpz_mpoly_init(temp_poly, ctx);

    do {
        found_divisor = 0;
        fmpz_mpoly_leadterm(leadterm, poly, ctx);

        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_t vec_poly;
            fmpz_mpoly_init(vec_poly, ctx);
            fmpz_mpoly_set(vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);

            if (fmpz_mpoly_divides(gcd_leadterm, leadterm, vec_poly, ctx)) {
                found_divisor = 1;
                fprintf(stderr, "Divisor found: %s\n", fmpz_mpoly_get_str_pretty(vec_poly, (const char **) varnames, ctx));

                fmpz_mpoly_mul(temp_poly, vec_poly, gcd_leadterm, ctx);
                fmpz_mpoly_sub(poly, poly, temp_poly, ctx);
                fmpz_mpoly_leadterm(leadterm, poly, ctx);
                break;
            }

            fmpz_mpoly_clear(vec_poly, ctx);
        }
    } while (found_divisor && !fmpz_mpoly_is_zero(poly, ctx));

    fprintf(stderr, "Output: %s\n", fmpz_mpoly_get_str_pretty(poly, (const char **) varnames, ctx));

    for (i = 0; i < nvars; i++) {
        flint_free(varnames[i]);
    }
    flint_free(varnames);

    fmpz_mpoly_clear(leadterm, ctx);
    fmpz_mpoly_clear(gcd_leadterm, ctx);
    fmpz_mpoly_clear(temp_poly, ctx);
}
