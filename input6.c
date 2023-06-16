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
    fmpz_mpoly_vec_t basis;
    fmpz_mpoly_vec_init(basis, 0, ctx);
    buchberger_naive(basis, gens, ctx);

    for (slong i = 0; i < basis->length; ++i) {
        fmpz_mpoly_t poly, reduced_poly;
        fmpz_mpoly_init(poly, ctx);
        fmpz_mpoly_init(reduced_poly, ctx);
        fmpz_mpoly_set(poly, fmpz_mpoly_vec_entry(basis, i), ctx);

        fmpz_mpoly_vec_t others;
        fmpz_mpoly_vec_init(others, 0, ctx);
        for (slong j = 0; j < basis->length; ++j) {
            if (j != i) {
                fmpz_mpoly_vec_append(others, fmpz_mpoly_vec_entry(basis, j), ctx);
            }
        }

        reduce_by_vector(poly, others, ctx);

        if (fmpz_mpoly_is_zero(poly, ctx)) {
            fmpz_mpoly_vec_t new_basis;
            fmpz_mpoly_vec_init(new_basis, 0, ctx);
            for (slong k = 0; k < basis->length; ++k) {
                if (k != i) {
                    fmpz_mpoly_vec_append(new_basis, fmpz_mpoly_vec_entry(basis, k), ctx);
                }
            }
            fmpz_mpoly_vec_swap(basis, new_basis, ctx);
            fmpz_mpoly_vec_clear(new_basis, ctx);
            --i;
        } else if (!fmpz_mpoly_equal(poly, fmpz_mpoly_vec_entry(basis, i), ctx)) {
            fmpz_mpoly_swap(poly, fmpz_mpoly_vec_entry(basis, i), ctx);
        }

        fmpz_mpoly_vec_clear(others, ctx);
        fmpz_mpoly_clear(poly, ctx);
        fmpz_mpoly_clear(reduced_poly, ctx);
    }

    fmpz_mpoly_vec_swap(res, basis, ctx);
    fmpz_mpoly_vec_clear(basis, ctx);
}
