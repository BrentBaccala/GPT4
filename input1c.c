#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t);
void construct_s_pair(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t);
void reduce_by_vector(fmpz_mpoly_t, const fmpz_mpoly_vec_t, const int, const fmpz_mpoly_ctx_t);
void buchberger_naive(fmpz_mpoly_vec_t, const fmpz_mpoly_vec_t, const fmpz_mpoly_ctx_t);
void buchberger_reduced(fmpz_mpoly_vec_t, const fmpz_mpoly_vec_t, const fmpz_mpoly_ctx_t);

void fmpz_mpoly_get_term(fmpz_mpoly_t output, const fmpz_mpoly_t input, slong term_num, const fmpz_mpoly_ctx_t ctx)
{
    if (term_num == 0)
    {
        fmpz_mpoly_leadterm(output, input, ctx);
    }
    else
    {
        fmpz_mpoly_t temp;
        fmpz_mpoly_init(temp, ctx);
        fmpz_mpoly_set(temp, input, ctx);
        fmpz_mpoly_truncate(temp, term_num, ctx);
        fmpz_mpoly_sub(temp, input, temp, ctx);
        fmpz_mpoly_leadterm(output, temp, ctx);
        fmpz_mpoly_clear(temp, ctx);
    }
}
