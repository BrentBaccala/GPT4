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
    slong i, nvars;
    char **var_names;
    fmpz_mpoly_t lt_poly, temp_poly, gcd_term;

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(temp_poly, ctx);
    fmpz_mpoly_init(gcd_term, ctx);

    nvars = ctx->minfo->nvars;
    var_names = (char **) malloc(nvars * sizeof(char *));
    for (i = 0; i < nvars; i++) {
        var_names[i] = (char *) malloc(16 * sizeof(char));
        snprintf(var_names[i], 16, "x%ld", i + 1);
    }

    char *poly_str = fmpz_mpoly_get_str_pretty(poly, (const char **) var_names, ctx);
    fprintf(stderr, "Input polynomial: %s\n", poly_str);
    flint_free(poly_str);

    int progress;
    do {
        progress = 0;
        for (i = 0; i < vec->length && !fmpz_mpoly_is_zero(poly, ctx); i++) {
            fmpz_mpoly_t vec_poly;
            fmpz_mpoly_init_set(vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);

            fmpz_mpoly_leadterm(lt_poly, vec_poly, ctx);
            if (fmpz_mpoly_divides(gcd_term, poly, lt_poly, ctx)) {
                poly_str = fmpz_mpoly_get_str_pretty(vec_poly, (const char **) var_names, ctx);
                fprintf(stderr, "Divisor found: %s\n", poly_str);
                flint_free(poly_str);

                fmpz_mpoly_mul(temp_poly, vec_poly, gcd_term, ctx);
                fmpz_mpoly_sub(poly, poly, temp_poly, ctx);
                progress = 1;

                poly_str = fmpz_mpoly_get_str_pretty(poly, (const char **) var_names, ctx);
                fprintf(stderr, "Current polynomial: %s\n", poly_str);
                flint_free(poly_str);
            }

            fmpz_mpoly_clear(vec_poly, ctx);
        }
    } while (progress && !fmpz_mpoly_is_zero(poly, ctx));

    poly_str = fmpz_mpoly_get_str_pretty(poly, (const char **) var_names, ctx);
    fprintf(stderr, "Reduced polynomial: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(temp_poly, ctx);
    fmpz_mpoly_clear(gcd_term, ctx);

    for (i = 0; i < nvars; i++) {
        free(var_names[i]);
    }
    free(var_names);
}
