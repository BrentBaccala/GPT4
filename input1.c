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
    slong i;
    int changed;
    fmpz_mpoly_t leadterm, temp_poly, gcd_term;

    char **varnames = (char **) flint_malloc(ctx->minfo->nvars * sizeof(char *));
    for (i = 0; i < ctx->minfo->nvars; i++) {
        varnames[i] = (char *) flint_malloc(16 * sizeof(char));
        flint_sprintf(varnames[i], "x%wd", i + 1);
    }

    char *poly_str = fmpz_mpoly_get_str_pretty(poly, (const char **) varnames, ctx);
    flint_fprintf(stderr, "Input polynomial: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_init(leadterm, ctx);
    fmpz_mpoly_init(temp_poly, ctx);
    fmpz_mpoly_init(gcd_term, ctx);

    do {
        changed = 0;
        fmpz_mpoly_leadterm(leadterm, poly, ctx);

        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_t vec_poly;
            fmpz_mpoly_init(vec_poly, ctx);
            fmpz_mpoly_set(vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);

            if (fmpz_mpoly_divides(temp_poly, leadterm, fmpz_mpoly_vec_entry(vec, i), ctx)) {
                fmpz_mpoly_gcd(gcd_term, leadterm, fmpz_mpoly_vec_entry(vec, i), ctx);
                fmpz_mpoly_mul(temp_poly, vec_poly, gcd_term, ctx);
                fmpz_mpoly_sub(poly, poly, temp_poly, ctx);

                poly_str = fmpz_mpoly_get_str_pretty(poly, (const char **) varnames, ctx);
                flint_fprintf(stderr, "Reduced by %wd-th polynomial: %s\n", i, poly_str);
                flint_free(poly_str);

                fmpz_mpoly_leadterm(leadterm, poly, ctx);
                changed = 1;
                break;
            }

            fmpz_mpoly_clear(vec_poly);
        }
    } while (changed && !fmpz_mpoly_is_zero(poly, ctx));

    poly_str = fmpz_mpoly_get_str_pretty(poly, (const char **) varnames, ctx);
    flint_fprintf(stderr, "Resulting polynomial: %s\n", poly_str);
    flint_free(poly_str);

    for (i = 0; i < ctx->minfo->nvars; i++) {
        flint_free(varnames[i]);
    }
    flint_free(varnames);

    fmpz_mpoly_clear(leadterm);
    fmpz_mpoly_clear(temp_poly);
    fmpz_mpoly_clear(gcd_term);
}
