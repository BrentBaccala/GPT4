#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t output, const fmpz_mpoly_t input, const fmpz_mpoly_ctx_t ctx);

void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);

void reduce_by_matching_term(fmpz_mpoly_t poly, const fmpz_mpoly_t matching_term, const fmpz_mpoly_t reducer, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t reducer_vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t reducer_vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lead_term, lead_term_reducer, temp_poly;
    slong i, j;
    const char *varnames[] = {"x", "y", "z"};

    // Log input/output polynomial
    char *poly_str = fmpz_mpoly_get_str_pretty(poly, varnames, ctx);
    flint_fprintf(stderr, "Input polynomial: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_init(lead_term, ctx);
    fmpz_mpoly_init(lead_term_reducer, ctx);
    fmpz_mpoly_init(temp_poly, ctx);

    while (!fmpz_mpoly_is_zero(poly, ctx)) {
        fmpz_mpoly_leadterm(lead_term, poly, ctx);

        int reduced = 0;
        for (i = 0; i < reducer_vec->length; i++) {
            fmpz_mpoly_t reducer = fmpz_mpoly_vec_entry(reducer_vec, i);
            fmpz_mpoly_leadterm(lead_term_reducer, reducer, ctx);

            if (lead_reduction) {
                if (fmpz_mpoly_is_divisible(lead_term, lead_term_reducer, ctx)) {
                    reduce_by_matching_term(poly, lead_term, reducer, ctx);
                    reduced = 1;
                    break;
                }
            } else {
                for (j = 0; j < fmpz_mpoly_length(poly, ctx); j++) {
                    fmpz_mpoly_t term_j;
                    fmpz_mpoly_init(term_j, ctx);
                    fmpz_mpoly_get_term_monomial(term_j, poly, j, ctx);

                    if (fmpz_mpoly_is_divisible(term_j, lead_term_reducer, ctx)) {
                        reduce_by_matching_term(poly, term_j, reducer, ctx);
                        reduced = 1;
                    }
                    fmpz_mpoly_clear(term_j, ctx);
                }
            }
        }

        if (!reduced) {
            break;
        }
    }

    // Log the result of the reduction
    poly_str = fmpz_mpoly_get_str_pretty(poly, varnames, ctx);
    flint_fprintf(stderr, "Reduced polynomial: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_clear(lead_term, ctx);
    fmpz_mpoly_clear(lead_term_reducer, ctx);
    fmpz_mpoly_clear(temp_poly, ctx);
}
