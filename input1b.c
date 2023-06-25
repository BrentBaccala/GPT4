#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);

int fmpz_mpoly_is_divisible(const fmpz_mpoly_t dividend, const fmpz_mpoly_t divisor, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t dummy;
    int divisible;

    fmpz_mpoly_init(dummy, ctx);
    divisible = fmpz_mpoly_divides(dummy, dividend, divisor, ctx);
    fmpz_mpoly_clear(dummy, ctx);

    return divisible;
}
