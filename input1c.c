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

void fmpz_mpoly_get_term(fmpz_mpoly_t rop, const fmpz_mpoly_t poly, slong term_num, const fmpz_mpoly_ctx_t ctx) {
    if (term_num == 0) {
        fmpz_mpoly_leadterm(rop, poly, ctx);
    } else {
        fmpz_mpoly_t truncated_poly, diff_poly;
        fmpz_mpoly_init(truncated_poly, ctx);
        fmpz_mpoly_init(diff_poly, ctx);
        fmpz_mpoly_truncate(truncated_poly, poly, term_num, ctx);
        fmpz_mpoly_sub(diff_poly, poly, truncated_poly, ctx);
        fmpz_mpoly_leadterm(rop, diff_poly, ctx);
        fmpz_mpoly_clear(truncated_poly, ctx);
        fmpz_mpoly_clear(diff_poly, ctx);
    }
}
