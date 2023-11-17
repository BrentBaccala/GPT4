#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t);
void construct_s_pair(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t);
void reduce_by_vector(fmpz_mpoly_t, const fmpz_mpoly_vec_t, int, const fmpz_mpoly_ctx_t);
void buchberger_naive(fmpz_mpoly_vec_t, const fmpz_mpoly_vec_t, const fmpz_mpoly_ctx_t);
void buchberger_reduced(fmpz_mpoly_vec_t, const fmpz_mpoly_vec_t, const fmpz_mpoly_ctx_t);
void fmpz_mpoly_get_term(fmpz_mpoly_t, const fmpz_mpoly_t, slong, const fmpz_mpoly_ctx_t);
void reduce_by_matching_term(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t leadterm, vec_leadterm;
    slong i, j;
    const char * vars[] = {"x", "y", "z"};

    flint_printf("Input to reduce_by_vector: %s\n", fmpz_mpoly_get_str_pretty(poly, vars, ctx));

    fmpz_mpoly_init(leadterm, ctx);
    fmpz_mpoly_init(vec_leadterm, ctx);

    fmpz_mpoly_leadterm(leadterm, poly, ctx);

    for (i = 0; i < vec->length; i++) {
        fmpz_mpoly_t vec_poly;
        fmpz_mpoly_init(vec_poly, ctx);
        fmpz_mpoly_set(vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);

        fmpz_mpoly_leadterm(vec_leadterm, vec_poly, ctx);

        if (lead_reduction) {
            if (fmpz_mpoly_is_divisible(leadterm, vec_leadterm, ctx)) {
                reduce_by_matching_term(poly, leadterm, vec_poly, ctx);
                if (fmpz_mpoly_is_zero(poly, ctx)) break;
                fmpz_mpoly_leadterm(leadterm, poly, ctx);
                i = -1; // restart the loop
            }
        } else {
            for (j = 0; j < fmpz_mpoly_length(poly, ctx); j++) {
                fmpz_mpoly_t term;
                fmpz_mpoly_init(term, ctx);
                fmpz_mpoly_get_term(term, poly, j, ctx);
                if (fmpz_mpoly_is_divisible(term, vec_leadterm, ctx)) {
                    reduce_by_matching_term(poly, term, vec_poly, ctx);
                    if (fmpz_mpoly_is_zero(poly, ctx)) break;
                }
                fmpz_mpoly_clear(term, ctx);
            }
        }

        fmpz_mpoly_clear(vec_poly, ctx);
    }

    fmpz_mpoly_clear(leadterm, ctx);
    fmpz_mpoly_clear(vec_leadterm, ctx);

    flint_printf("Output of reduce_by_vector: %s\n", fmpz_mpoly_get_str_pretty(poly, vars, ctx));
}
