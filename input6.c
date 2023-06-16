#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t output, const fmpz_mpoly_t input, const fmpz_mpoly_ctx_t ctx);

void construct_s_pair(fmpz_mpoly_t output, const fmpz_mpoly_t f1, const fmpz_mpoly_t f2, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const int lead_reduction_flag, const fmpz_mpoly_ctx_t ctx);

void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);

void buchberger_reduced(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx) {
    slong i, j, n;
    fmpz_mpoly_vec_t basis, reduced_basis;
    fmpz_mpoly_t temp_poly;

    buchberger_naive(output, input, ctx);
    fmpz_mpoly_vec_init(reduced_basis, 0, ctx);
    fmpz_mpoly_init(temp_poly, ctx);

    for (i = 0; i < output->length; i++) {
        fmpz_mpoly_set(temp_poly, fmpz_mpoly_vec_entry(output, i), ctx);
        fmpz_mpoly_vec_init(basis, output->length - 1, ctx);
        n = 0;

        for (j = 0; j < output->length; j++) {
            if (j != i) {
                fmpz_mpoly_set(fmpz_mpoly_vec_entry(basis, n), fmpz_mpoly_vec_entry(output, j), ctx);
                n++;
            }
        }

        reduce_by_vector(temp_poly, basis, 1, ctx);

        if (!fmpz_mpoly_is_zero(temp_poly, ctx)) {
            fmpz_mpoly_vec_append(reduced_basis, temp_poly, ctx);
        }

        fmpz_mpoly_vec_clear(basis, ctx);
    }

    fmpz_mpoly_vec_swap(output, reduced_basis, ctx);
    fmpz_mpoly_vec_clear(reduced_basis, ctx);
    fmpz_mpoly_clear(temp_poly, ctx);
}
