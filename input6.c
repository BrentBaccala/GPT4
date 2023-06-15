#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void buchberger_reduced(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gens, const fmpz_mpoly_ctx_t ctx) {
    slong i, j;
    fmpz_mpoly_vec_t gens_copy;
    fmpz_mpoly_t tmp1, tmp2, red;

    fmpz_mpoly_vec_init(gens_copy, 0, ctx);
    for (i = 0; i < gens->length; i++) {
        fmpz_mpoly_vec_append(gens_copy, fmpz_mpoly_vec_entry(gens, i), ctx);
    }

    buchberger_naive(res, gens_copy, ctx);
    fmpz_mpoly_vec_clear(gens_copy, ctx);

    fmpz_mpoly_init(tmp1, ctx);
    fmpz_mpoly_init(tmp2, ctx);
    fmpz_mpoly_init(red, ctx);

    i = 0;
    while (i < res->length) {
        fmpz_mpoly_set(tmp1, fmpz_mpoly_vec_entry(res, i), ctx);
        for (j = 0; j < res->length; j++) {
            if (j != i) {
                fmpz_mpoly_set(tmp2, fmpz_mpoly_vec_entry(res, j), ctx);
                reduce_by_vector(red, tmp1, tmp2, ctx);
                fmpz_mpoly_swap(tmp1, red, ctx);
            }
        }

        if (fmpz_mpoly_is_zero(tmp1, ctx)) {
            fmpz_mpoly_vec_remove_entry(res, i, ctx);
        } else {
            fmpz_mpoly_swap(fmpz_mpoly_vec_entry(res, i), tmp1, ctx);
            i++;
        }
    }

    fmpz_mpoly_clear(tmp1, ctx);
    fmpz_mpoly_clear(tmp2, ctx);
    fmpz_mpoly_clear(red, ctx);
}
