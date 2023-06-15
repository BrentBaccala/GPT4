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
    fmpz_mpoly_t s_pair, reduced, poly1, poly2;
    slong i, j, k;

    fmpz_mpoly_vec_init(basis, 0, ctx);
    fmpz_mpoly_vec_init(s_pairs, 0, ctx);

    for (i = 0; i < gens->length; i++)
    {
        fmpz_mpoly_vec_append(basis, fmpz_mpoly_vec_entry(gens, i), ctx);
    }

    for (i = 0; i < gens->length - 1; i++)
    {
        for (j = i + 1; j < gens->length; j++)
        {
            fmpz_mpoly_init(poly1, ctx);
            fmpz_mpoly_init(poly2, ctx);
            fmpz_mpoly_set(poly1, fmpz_mpoly_vec_entry(gens, i), ctx);
            fmpz_mpoly_set(poly2, fmpz_mpoly_vec_entry(gens, j), ctx);
            fmpz_mpoly_init(s_pair, ctx);
            construct_s_pair(s_pair, poly1, poly2, ctx);
            fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
            fmpz_mpoly_clear(poly1, ctx);
            fmpz_mpoly_clear(poly2, ctx);
            fmpz_mpoly_clear(s_pair, ctx);
        }
    }

    for (i = 0; i < s_pairs->length; i++)
    {
        fmpz_mpoly_init(reduced, ctx);
        fmpz_mpoly_set(reduced, fmpz_mpoly_vec_entry(s_pairs, i), ctx);
        reduce_by_vector(reduced, basis, ctx);

        if (!fmpz_mpoly_is_zero(reduced, ctx))
        {
            fmpz_mpoly_vec_append(basis, reduced, ctx);

            for (j = 0; j < basis->length - 1; j++)
            {
                fmpz_mpoly_init(poly1, ctx);
                fmpz_mpoly_init(poly2, ctx);
                fmpz_mpoly_set(poly1, fmpz_mpoly_vec_entry(basis, j), ctx);
                fmpz_mpoly_set(poly2, fmpz_mpoly_vec_entry(basis, basis->length - 1), ctx);
                fmpz_mpoly_init(s_pair, ctx);
                construct_s_pair(s_pair, poly1, poly2, ctx);

                int is_unique = 1;
                for (k = 0; k < s_pairs->length; k++)
                {
                    if (fmpz_mpoly_equal(s_pair, fmpz_mpoly_vec_entry(s_pairs, k), ctx))
                    {
                        is_unique = 0;
                        break;
                    }
                }

                if (is_unique)
                {
                    fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
                }

                fmpz_mpoly_clear(poly1, ctx);
                fmpz_mpoly_clear(poly2, ctx);
                fmpz_mpoly_clear(s_pair, ctx);
            }
        }

        fmpz_mpoly_clear(reduced, ctx);
    }

    for (i = 0; i < basis->length; i++)
    {
        fmpz_mpoly_vec_append(res, fmpz_mpoly_vec_entry(basis, i), ctx);
    }

    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(s_pairs, ctx);
}
