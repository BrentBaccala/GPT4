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
    int found;
    slong i;
    fmpz_mpoly_t lead_poly, g, temp_poly;

    fmpz_mpoly_init(lead_poly, ctx);
    fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_init(temp_poly, ctx);

    do {
        found = 0;
        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_t current_poly;
            fmpz_mpoly_init(current_poly, ctx);
            fmpz_mpoly_set(current_poly, fmpz_mpoly_vec_entry(vec, i), ctx);

            fmpz_mpoly_leadterm(lead_poly, current_poly, ctx);
            if (fmpz_mpoly_divides(g, poly, lead_poly, ctx)) {
                found = 1;
                fmpz_mpoly_mul(temp_poly, current_poly, g, ctx);
                fmpz_mpoly_sub(poly, poly, temp_poly, ctx);

                char *poly_str = fmpz_mpoly_get_str_pretty(poly, NULL, ctx);
                fprintf(stderr, "Reduced poly: %s\n", poly_str);
                flint_free(poly_str);

                if (fmpz_mpoly_is_zero(poly, ctx)) {
                    break;
                }
            }
            fmpz_mpoly_clear(current_poly, ctx);
        }
    } while (found && !fmpz_mpoly_is_zero(poly, ctx));

    fmpz_mpoly_clear(lead_poly, ctx);
    fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_clear(temp_poly, ctx);
}

