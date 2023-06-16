#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_red_flag, const fmpz_mpoly_ctx_t ctx);

void buchberger_naive(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gen, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_vec_t basis, spairs;
    fmpz_mpoly_vec_init(basis, 0, ctx);
    fmpz_mpoly_vec_init(spairs, 0, ctx);

    for (slong i = 0; i < gen->length; i++) {
        fmpz_mpoly_vec_append(basis, fmpz_mpoly_vec_entry(gen, i), ctx);
    }

    for (slong i = 0; i < basis->length; i++) {
        for (slong j = i + 1; j < basis->length; j++) {
            fmpz_mpoly_t s_pair;
            fmpz_mpoly_init(s_pair, ctx);
            construct_s_pair(s_pair, fmpz_mpoly_vec_entry(basis, i), fmpz_mpoly_vec_entry(basis, j), ctx);
            fmpz_mpoly_vec_append(spairs, s_pair, ctx);
            fmpz_mpoly_clear(s_pair, ctx);
        }
    }

    slong idx = 0;
    while (idx < spairs->length) {
        fmpz_mpoly_t red_poly;
        fmpz_mpoly_init(red_poly, ctx);
        fmpz_mpoly_set(red_poly, fmpz_mpoly_vec_entry(spairs, idx), ctx);
        reduce_by_vector(red_poly, basis, 1, ctx);

        if (!fmpz_mpoly_is_zero(red_poly, ctx)) {
            fmpz_mpoly_vec_append(basis, red_poly, ctx);
            for (slong i = 0; i < basis->length - 1; i++) {
                fmpz_mpoly_t new_s_pair;
                fmpz_mpoly_init(new_s_pair, ctx);
                construct_s_pair(new_s_pair, fmpz_mpoly_vec_entry(basis, i), fmpz_mpoly_vec_entry(basis, basis->length - 1), ctx);
                fmpz_mpoly_vec_append(spairs, new_s_pair, ctx);
                fmpz_mpoly_clear(new_s_pair, ctx);
            }
        }
        fmpz_mpoly_clear(red_poly, ctx);
        idx++;
    }

    for (slong i = 0; i < basis->length; i++) {
        fmpz_mpoly_vec_append(res, fmpz_mpoly_vec_entry(basis, i), ctx);
    }

    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(spairs, ctx);
}
