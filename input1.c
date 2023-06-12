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
    int reduced;
    fmpz_mpoly_t temp1, temp2, lead_term;

    fmpz_mpoly_init(temp1, ctx);
    fmpz_mpoly_init(temp2, ctx);
    fmpz_mpoly_init(lead_term, ctx);

    char **var_names = (char **)malloc(ctx->minfo->nvars * sizeof(char *));
    for (slong i = 0; i < ctx->minfo->nvars; i++) {
        var_names[i] = (char *)malloc(4 * sizeof(char));
        snprintf(var_names[i], 4, "x%ld", i + 1);
    }

    char *poly_str = fmpz_mpoly_get_str_pretty(poly, var_names, ctx);
    fprintf(stderr, "Input polynomial: %s\n", poly_str);
    flint_free(poly_str);

    do {
        reduced = 0;
        for (slong i = 0; i < vec->length; i++) {
            fmpz_mpoly_t vec_poly;
            fmpz_mpoly_init(vec_poly, ctx);
            fmpz_mpoly_set(vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);

            fmpz_mpoly_leadterm(lead_term, vec_poly, ctx);

            if (fmpz_mpoly_divides(temp1, poly, lead_term, ctx)) {
                fmpz_mpoly_gcd(temp2, lead_term, poly, ctx);
                fmpz_mpoly_divexact(temp1, poly, temp2, ctx);
                fmpz_mpoly_mul(temp1, temp1, vec_poly, ctx);
                fmpz_mpoly_sub(poly, poly, temp1, ctx);

                char *vec_poly_str = fmpz_mpoly_get_str_pretty(vec_poly, var_names, ctx);
                fprintf(stderr, "Reduced by: %s\n", vec_poly_str);
                flint_free(vec_poly_str);

                reduced = 1;
                break;
            }

            fmpz_mpoly_clear(vec_poly, ctx);
        }
    } while (reduced && !fmpz_mpoly_is_zero(poly, ctx));

    poly_str = fmpz_mpoly_get_str_pretty(poly, var_names, ctx);
    fprintf(stderr, "Reduced polynomial: %s\n", poly_str);
    flint_free(poly_str);

    for (slong i = 0; i < ctx->minfo->nvars; i++) {
        free(var_names[i]);
    }
    free(var_names);

    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
    fmpz_mpoly_clear(lead_term, ctx);
}
