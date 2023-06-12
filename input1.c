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
    fmpz_mpoly_t lead_poly, divisor, temp;
    int reduction_occurred;

    fmpz_mpoly_init(lead_poly, ctx);
    fmpz_mpoly_init(divisor, ctx);
    fmpz_mpoly_init(temp, ctx);

    fprintf(stderr, "Input polynomial: %s\n", fmpz_mpoly_get_str_pretty(poly, NULL, ctx));

    do {
        reduction_occurred = 0;
        for (slong i = 0; i < vec->length && !fmpz_mpoly_is_zero(poly, ctx); i++) {
            fmpz_mpoly_t current_poly = fmpz_mpoly_vec_entry(vec, i);
            fmpz_mpoly_leadterm(lead_poly, current_poly, ctx);

            if (fmpz_mpoly_divides(temp, poly, lead_poly, ctx)) {
                fprintf(stderr, "Reducing by %s\n", fmpz_mpoly_get_str_pretty(current_poly, NULL, ctx));
                fmpz_mpoly_mul(temp, temp, current_poly, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                reduction_occurred = 1;
                break;
            }
        }
    } while (reduction_occurred);

    fprintf(stderr, "Reduced polynomial: %s\n", fmpz_mpoly_get_str_pretty(poly, NULL, ctx));

    fmpz_mpoly_clear(lead_poly, ctx);
    fmpz_mpoly_clear(divisor, ctx);
    fmpz_mpoly_clear(temp, ctx);
}
