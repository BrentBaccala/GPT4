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

    for (slong i = 0; i < gens->length; i++)
    {
        fmpz_mpoly_vec_append(basis, fmpz_mpoly_vec_entry(gens, i), ctx);
        for (slong j = i + 1; j < gens->length; j++)
        {
            fmpz_mpoly_t sp;
            fmpz_mpoly_init(sp, ctx);
            construct_s_pair(sp, fmpz_mpoly_vec_entry(gens, i), fmpz_mpoly_vec_entry(gens, j), ctx);
            fmpz_mpoly_vec_append(s_pairs, sp, ctx);
            fmpz_mpoly_clear(sp, ctx);
        }
    }

    for (slong i = 0; i < s_pairs->length; i++)
    {
        fmpz_mpoly_t r;
        fmpz_mpoly_init(r, ctx);
        fmpz_mpoly_set(r, fmpz_mpoly_vec_entry(s_pairs, i), ctx);
        reduce_by_vector(r, basis, ctx);

        if (!fmpz_mpoly_is_zero(r, ctx))
        {
            fmpz_mpoly_vec_append(basis, r, ctx);
            for (slong j = 0; j < basis->length - 1; j++)
            {
                fmpz_mpoly_t sp;
                fmpz_mpoly_init(sp, ctx);
                construct_s_pair(sp, fmpz_mpoly_vec_entry(basis, j), fmpz_mpoly_vec_entry(basis, basis->length - 1), ctx);
                fmpz_mpoly_vec_append(s_pairs, sp, ctx);
                fmpz_mpoly_clear(sp, ctx);
            }
        }
        fmpz_mpoly_clear(r, ctx);
    }

    for (slong i = 0; i < basis->length; i++)
    {
        fmpz_mpoly_vec_append(res, fmpz_mpoly_vec_entry(basis, i), ctx);
    }

    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(s_pairs, ctx);
}
