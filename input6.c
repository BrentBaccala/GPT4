#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void buchberger_reduced(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gens, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_vec_t naive_basis;
    fmpz_mpoly_vec_init(naive_basis, 0, ctx);
    buchberger_naive(naive_basis, gens, ctx);

    slong i, j;
    slong len = naive_basis->length;
    fmpz_mpoly_t temp_poly, reduced_poly;
    fmpz_mpoly_vec_t temp_vec;
    fmpz_mpoly_vec_init(temp_vec, 0, ctx);

    for (i = 0; i < len; i++)
    {
        fmpz_mpoly_init(temp_poly, ctx);
        fmpz_mpoly_set(temp_poly, fmpz_mpoly_vec_entry(naive_basis, i), ctx);
        fmpz_mpoly_vec_append(temp_vec, temp_poly, ctx);
    }

    i = 0;
    while (i < temp_vec->length)
    {
        fmpz_mpoly_init(reduced_poly, ctx);
        fmpz_mpoly_set(reduced_poly, fmpz_mpoly_vec_entry(temp_vec, i), ctx);
        reduce_by_vector(reduced_poly, temp_vec, ctx);

        int changed = !fmpz_mpoly_equal(reduced_poly, fmpz_mpoly_vec_entry(temp_vec, i), ctx);

        if (fmpz_mpoly_is_zero(reduced_poly, ctx))
        {
            fmpz_mpoly_vec_swap_entry(temp_vec, i, temp_vec->length - 1, ctx);
            fmpz_mpoly_vec_set_length(temp_vec, temp_vec->length - 1, ctx);
        }
        else
        {
            if (changed)
            {
                fmpz_mpoly_set(fmpz_mpoly_vec_entry(temp_vec, i), reduced_poly, ctx);
            }
            i++;
        }

        fmpz_mpoly_clear(reduced_poly, ctx);
    }

    fmpz_mpoly_vec_swap(res, temp_vec, ctx);
    fmpz_mpoly_vec_clear(temp_vec, ctx);
    fmpz_mpoly_vec_clear(naive_basis, ctx);
}
