#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t res, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

void construct_s_pair(fmpz_mpoly_t spair, const fmpz_mpoly_t p1, const fmpz_mpoly_t p2, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t f, const fmpz_mpoly_vec_t F, int lead_reduc, const fmpz_mpoly_ctx_t ctx);

void buchberger_naive(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);

void buchberger_reduced(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_get_term(fmpz_mpoly_t term, const fmpz_mpoly_t poly, slong n, const fmpz_mpoly_ctx_t ctx) {
    if (n == 0) {
        fmpz_mpoly_leadterm(term, poly, ctx);
    } else {
        fmpz_mpoly_t truncated_poly;
        fmpz_mpoly_init(truncated_poly, ctx);
        fmpz_mpoly_set(truncated_poly, poly, ctx);
        fmpz_mpoly_truncate(truncated_poly, n, ctx);
        fmpz_mpoly_sub(truncated_poly, poly, truncated_poly, ctx);
        fmpz_mpoly_leadterm(term, truncated_poly, ctx);
        fmpz_mpoly_clear(truncated_poly, ctx);
    }
}
