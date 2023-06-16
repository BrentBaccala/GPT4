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

void buchberger_reduced(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gens, const fmpz_mpoly_ctx_t ctx) {
    buchberger_naive(res, gens, ctx);

    slong i, j;
    fmpz_mpoly_t reduced_poly;
    fmpz_mpoly_vec_t temp_vec;

    fmpz_mpoly_init(reduced_poly, ctx);
    fmpz_mpoly_vec_init(temp_vec, 0, ctx);

    for (i = 0; i < res->length; i++) {
        fmpz_mpoly_set(reduced_poly, fmpz_mpoly_vec_entry(res, i), ctx);
        reduce_by_vector(reduced_poly, res, ctx);

        if (!fmpz_mpoly_is_zero(reduced_poly, ctx)) {
            fmpz_mpoly_swap(fmpz_mpoly_vec_entry(res, i), reduced_poly, ctx);
        } else {
            for (j = 0; j < res->length; j++) {
                if (j != i) {
                    fmpz_mpoly_vec_append(temp_vec, fmpz_mpoly_vec_entry(res, j), ctx);
                }
            }
            fmpz_mpoly_vec_swap(res, temp_vec, ctx);
            fmpz_mpoly_vec_clear(temp_vec, ctx);
            fmpz_mpoly_vec_init(temp_vec, 0, ctx);
            i--;
        }
    }

    fmpz_mpoly_clear(reduced_poly, ctx);
    fmpz_mpoly_vec_clear(temp_vec, ctx);
}
