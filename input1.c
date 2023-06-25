#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t output, const fmpz_mpoly_t input, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t output, const fmpz_mpoly_t input1, const fmpz_mpoly_t input2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);
void reduce_by_matching_term(fmpz_mpoly_t poly, const fmpz_mpoly_t matching_term, const fmpz_mpoly_t reducer, const fmpz_mpoly_ctx_t ctx);

void reduce_by_matching_term(fmpz_mpoly_t poly, const fmpz_mpoly_t matching_term, const fmpz_mpoly_t reducer, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t gcd, temp;
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);
    fmpz_mpoly_gcd(gcd, matching_term, fmpz_mpoly_leadterm(temp, reducer, ctx), ctx);
    fmpz_mpoly_mul(temp, reducer, gcd, ctx);
    fmpz_mpoly_sub(poly, poly, temp, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lead_term, vec_lead_term, temp;
    fmpz_mpoly_init(lead_term, ctx);
    fmpz_mpoly_init(vec_lead_term, ctx);
    fmpz_mpoly_init(temp, ctx);

    fmpz_mpoly_leadterm(lead_term, poly, ctx);

    for (slong i = 0; i < vec->length; i++) {
        fmpz_mpoly_set(temp, fmpz_mpoly_vec_entry(vec, i), ctx);
        fmpz_mpoly_leadterm(vec_lead_term, temp, ctx);

        if (lead_reduction) {
            if (fmpz_mpoly_divides(lead_term, lead_term, vec_lead_term, ctx)) {
                reduce_by_matching_term(poly, lead_term, temp, ctx);
                fmpz_mpoly_leadterm(lead_term, poly, ctx);
                if (fmpz_mpoly_is_zero(poly, ctx)) {
                    break;
                }
            }
        } else {
            for (slong j = 0; j < poly->length; j++) {
                fmpz_mpoly_t term;
                fmpz_mpoly_init(term, ctx);
                fmpz_mpoly_get_term_coeff_fmpz(term, poly, j, ctx);

                if (fmpz_mpoly_divides(term, term, vec_lead_term, ctx)) {
                    reduce_by_matching_term(poly, term, temp, ctx);
                }

                fmpz_mpoly_clear(term, ctx);
            }
        }
    }

    fmpz_mpoly_clear(lead_term, ctx);
    fmpz_mpoly_clear(vec_lead_term, ctx);
    fmpz_mpoly_clear(temp, ctx);
}
