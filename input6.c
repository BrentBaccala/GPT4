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
void buchberger_reduced(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gens, const fmpz_mpoly_ctx_t ctx);

void buchberger_reduced(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gens, const fmpz_mpoly_ctx_t ctx) {
    buchberger_naive(res, gens, ctx);

    slong i, j, n = res->length;
    fmpz_mpoly_vec_t temp_vec;
    fmpz_mpoly_vec_init(temp_vec, 0, ctx);

    for (i = 0; i < n; i++) {
        fmpz_mpoly_t current_poly;
        fmpz_mpoly_init(current_poly, ctx);
        fmpz_mpoly_set(current_poly, fmpz_mpoly_vec_entry(res, i), ctx);

        fmpz_mpoly_vec_t reduced_basis;
        fmpz_mpoly_vec_init(reduced_basis, 0, ctx);

        for (j = 0; j < n; j++) {
            if (j != i) {
                fmpz_mpoly_vec_append(reduced_basis, fmpz_mpoly_vec_entry(res, j), ctx);
            }
        }

        reduce_by_vector(current_poly, reduced_basis, ctx);

        if (fmpz_mpoly_is_zero(current_poly, ctx)) {
            fmpz_mpoly_vec_swap_entry(res, i, n - 1, ctx);
            res->length--;
            n--;
            i--;
        } else {
            fmpz_mpoly_swap(fmpz_mpoly_vec_entry(res, i), current_poly, ctx);
        }

        fmpz_mpoly_clear(current_poly, ctx);
        fmpz_mpoly_vec_clear(reduced_basis, ctx);
    }

    fmpz_mpoly_vec_clear(temp_vec, ctx);
}
