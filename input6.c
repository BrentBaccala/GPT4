#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

void construct_s_pair(fmpz_mpoly_t spair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);

void buchberger_naive(fmpz_mpoly_vec_t basis, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx);

void buchberger_reduced(fmpz_mpoly_vec_t reduced_basis, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx) {
    slong i, j, len;
    fmpz_mpoly_vec_t basis, temp_vec;
    fmpz_mpoly_t temp_poly;

    buchberger_naive(basis, generators, ctx);
    len = basis->length;

    fmpz_mpoly_vec_init(reduced_basis, 0, ctx);
    fmpz_mpoly_vec_init(temp_vec, 0, ctx);
    fmpz_mpoly_init(temp_poly, ctx);

    for (i = 0; i < len; i++) {
        fmpz_mpoly_set(temp_poly, fmpz_mpoly_vec_entry(basis, i), ctx);

        fmpz_mpoly_vec_clear(temp_vec, ctx);
        fmpz_mpoly_vec_init(temp_vec, 0, ctx);

        for (j = 0; j < len; j++) {
            if (j != i) {
                fmpz_mpoly_vec_append(temp_vec, fmpz_mpoly_vec_entry(basis, j), ctx);
            }
        }

        reduce_by_vector(temp_poly, temp_vec, 1, ctx);

        if (!fmpz_mpoly_is_zero(temp_poly, ctx)) {
            fmpz_mpoly_vec_append(reduced_basis, temp_poly, ctx);
        }
    }

    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(temp_vec, ctx);
    fmpz_mpoly_clear(temp_poly, ctx);
}
