#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t);
void construct_s_pair(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t);
void reduce_by_vector(fmpz_mpoly_t, const fmpz_mpoly_vec_t, int, const fmpz_mpoly_ctx_t);

void buchberger_naive(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_vec_t B, S;
    fmpz_mpoly_t s, g;
    slong i, j;

    fmpz_mpoly_vec_init(B, 0, ctx);
    fmpz_mpoly_vec_init(S, 0, ctx);
    fmpz_mpoly_init(s, ctx);
    fmpz_mpoly_init(g, ctx);

    for (i = 0; i < F->length; i++)
    {
        fmpz_mpoly_set(g, fmpz_mpoly_vec_entry(F, i), ctx);
        fmpz_mpoly_vec_append(B, g, ctx);
    }

    for (i = 0; i < B->length; i++)
    {
        for (j = i + 1; j < B->length; j++)
        {
            construct_s_pair(s, fmpz_mpoly_vec_entry(B, i), fmpz_mpoly_vec_entry(B, j), ctx);
            fmpz_mpoly_vec_append(S, s, ctx);
        }
    }

    for (i = 0; i < S->length; i++)
    {
        fmpz_mpoly_set(s, fmpz_mpoly_vec_entry(S, i), ctx);
        reduce_by_vector(s, B, 1, ctx);

        if (!fmpz_mpoly_is_zero(s, ctx))
        {
            fmpz_mpoly_vec_append(B, s, ctx);

            for (j = 0; j < B->length - 1; j++)
            {
                construct_s_pair(g, fmpz_mpoly_vec_entry(B, j), s, ctx);
                fmpz_mpoly_vec_append(S, g, ctx);
            }
        }
    }

    for (i = 0; i < B->length; i++)
    {
        fmpz_mpoly_set(g, fmpz_mpoly_vec_entry(B, i), ctx);
        fmpz_mpoly_vec_append(G, g, ctx);
    }

    fmpz_mpoly_clear(s, ctx);
    fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_vec_clear(B, ctx);
    fmpz_mpoly_vec_clear(S, ctx);
}
