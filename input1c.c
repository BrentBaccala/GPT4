#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t output, const fmpz_mpoly_t input, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t output, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t inout, const fmpz_mpoly_vec_t polys, const int lead_reducing, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);

void fmpz_mpoly_get_term(fmpz_mpoly_t output, const fmpz_mpoly_t input, ulong term_num, const fmpz_mpoly_ctx_t ctx) {
    if (term_num == 0) {
        fmpz_mpoly_leadterm(output, input, ctx);
    } else {
        fmpz_mpoly_t truncated_input, temp_diff;
        fmpz_mpoly_init(truncated_input, ctx);
        fmpz_mpoly_init(temp_diff, ctx);

        fmpz_mpoly_truncate(truncated_input, input, term_num, ctx);
        fmpz_mpoly_sub(temp_diff, input, truncated_input, ctx);
        fmpz_mpoly_leadterm(output, temp_diff, ctx);

        fmpz_mpoly_clear(truncated_input, ctx);
        fmpz_mpoly_clear(temp_diff, ctx);
    }
}
