#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t output, const fmpz_mpoly_t input, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t output, const fmpz_mpoly_t input1, const fmpz_mpoly_t input2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);

void buchberger_reduced(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_vec_t basis;
    fmpz_mpoly_vec_init(basis, 0, ctx);

    buchberger_naive(basis, input, ctx);

    for (slong i = 0; i < basis->length; ++i) {
        fmpz_mpoly_vec_t other_polynomials;
        fmpz_mpoly_vec_init(other_polynomials, 0, ctx);

        for (slong j = 0; j < basis->length; ++j) {
            if (j != i) {
                fmpz_mpoly_vec_append(other_polynomials, fmpz_mpoly_vec_entry(basis, j), ctx);
            }
        }

        fmpz_mpoly_t reduced_poly;
        fmpz_mpoly_init(reduced_poly, ctx);
        fmpz_mpoly_set(reduced_poly, fmpz_mpoly_vec_entry(basis, i), ctx);

        reduce_by_vector(reduced_poly, other_polynomials, 0, ctx);
        fmpz_mpoly_swap(fmpz_mpoly_vec_entry(basis, i), reduced_poly, ctx);

        fmpz_mpoly_clear(reduced_poly, ctx);
        fmpz_mpoly_vec_clear(other_polynomials, ctx);
    }

    fmpz_mpoly_vec_init(output, 0, ctx);
    for (slong i = 0; i < basis->length; ++i) {
        if (!fmpz_mpoly_is_zero(fmpz_mpoly_vec_entry(basis, i), ctx)) {
            fmpz_mpoly_vec_append(output, fmpz_mpoly_vec_entry(basis, i), ctx);
        }
    }

    fmpz_mpoly_vec_clear(basis, ctx);
}
