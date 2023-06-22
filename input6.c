#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t out, const fmpz_mpoly_t in, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t out, const fmpz_mpoly_t a, const fmpz_mpoly_t b, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t out, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t out, const fmpz_mpoly_vec_t in, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t out, const fmpz_mpoly_vec_t in, const fmpz_mpoly_ctx_t ctx);

void buchberger_reduced(fmpz_mpoly_vec_t out, const fmpz_mpoly_vec_t in, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_vec_t basis, reduced_basis;
    fmpz_mpoly_vec_init(basis, 0, ctx);
    fmpz_mpoly_vec_init(reduced_basis, 0, ctx);

    buchberger_naive(basis, in, ctx);

    for (slong i = 0; i < basis->length; i++) {
        fmpz_mpoly_t reduced_poly;
        fmpz_mpoly_vec_t others;

        fmpz_mpoly_init(reduced_poly, ctx);
        fmpz_mpoly_vec_init(others, 0, ctx);

        for (slong j = 0; j < basis->length; j++) {
            if (j != i) {
                fmpz_mpoly_vec_append(others, fmpz_mpoly_vec_entry(basis, j), ctx);
            }
        }

        fmpz_mpoly_set(reduced_poly, fmpz_mpoly_vec_entry(basis, i), ctx);
        reduce_by_vector(reduced_poly, others, 0, ctx);

        if (fmpz_mpoly_is_zero(reduced_poly, ctx)) {
            fmpz_mpoly_clear(reduced_poly, ctx);
        } else {
            fmpz_mpoly_swap(reduced_poly, fmpz_mpoly_vec_entry(basis, i), ctx);
            fmpz_mpoly_vec_append(reduced_basis, reduced_poly, ctx);
            fmpz_mpoly_clear(reduced_poly, ctx);
        }

        fmpz_mpoly_vec_clear(others, ctx);
    }

    fmpz_mpoly_vec_swap(out, reduced_basis, ctx);
    fmpz_mpoly_vec_clear(reduced_basis, ctx);
    fmpz_mpoly_vec_clear(basis, ctx);
}
