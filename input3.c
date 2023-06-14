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

void buchberger_naive(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gens, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_vec_t basis, s_pairs;
    fmpz_mpoly_vec_init(basis, 0, ctx);
    fmpz_mpoly_vec_init(s_pairs, 0, ctx);
    
    for (slong i = 0; i < gens->length; i++) {
        fmpz_mpoly_vec_append(basis, fmpz_mpoly_vec_entry(gens, i), ctx);
    }

    for (slong i = 0; i < gens->length; i++) {
        for (slong j = i + 1; j < gens->length; j++) {
            fmpz_mpoly_t s_pair;
            fmpz_mpoly_init(s_pair, ctx);
            construct_s_pair(s_pair, fmpz_mpoly_vec_entry(gens, i), fmpz_mpoly_vec_entry(gens, j), ctx);
            fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
            fmpz_mpoly_clear(s_pair, ctx);
        }
    }

    slong idx = 0;
    while (idx < s_pairs->length) {
        fmpz_mpoly_t reduced;
        fmpz_mpoly_init(reduced, ctx);
        fmpz_mpoly_set(reduced, fmpz_mpoly_vec_entry(s_pairs, idx), ctx);
        reduce_by_vector(reduced, basis, ctx);

        if (!fmpz_mpoly_is_zero(reduced, ctx)) {
            fmpz_mpoly_vec_append(basis, reduced, ctx);
            slong basis_length = basis->length;
            for (slong i = 0; i < basis_length - 1; i++) {
                fmpz_mpoly_t s_pair;
                fmpz_mpoly_init(s_pair, ctx);
                construct_s_pair(s_pair, fmpz_mpoly_vec_entry(basis, i), fmpz_mpoly_vec_entry(basis, basis_length - 1), ctx);
                fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
                fmpz_mpoly_clear(s_pair, ctx);
            }
        }

        fmpz_mpoly_clear(reduced, ctx);
        idx++;
    }

    for (slong i = 0; i < basis->length; i++) {
        fmpz_mpoly_vec_append(res, fmpz_mpoly_vec_entry(basis, i), ctx);
    }

    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(s_pairs, ctx);
}
