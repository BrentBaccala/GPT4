#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void buchberger_reduced(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gens, const fmpz_mpoly_ctx_t ctx) {
    slong i, j;
    int changed, is_zero;
    fmpz_mpoly_t poly, reduced_poly;
    fmpz_mpoly_vec_t temp_vec;

    buchberger_naive(res, gens, ctx);

    fmpz_mpoly_init(poly, ctx);
    fmpz_mpoly_init(reduced_poly, ctx);
    fmpz_mpoly_vec_init(temp_vec, 0, ctx);

    do {
        changed = 0;
        for (i = 0; i < res->length; i++) {
            fmpz_mpoly_set(poly, fmpz_mpoly_vec_entry(res, i), ctx);
            fmpz_mpoly_vec_set(temp_vec, res, ctx);
            fmpz_mpoly_vec_remove_index(temp_vec, i, ctx);

            reduce_by_vector(reduced_poly, poly, temp_vec, ctx);
            is_zero = fmpz_mpoly_is_zero(reduced_poly, ctx);

            if (!is_zero && !fmpz_mpoly_equal(poly, reduced_poly, ctx)) {
                fmpz_mpoly_swap(reduced_poly, fmpz_mpoly_vec_entry(res, i), ctx);
                changed = 1;
            } else if (is_zero) {
                fmpz_mpoly_vec_remove_index(res, i, ctx);
                i--;
                changed = 1;
            }
        }
    } while (changed);

    fmpz_mpoly_clear(poly, ctx);
    fmpz_mpoly_clear(reduced_poly, ctx);
    fmpz_mpoly_vec_clear(temp_vec, ctx);
}
