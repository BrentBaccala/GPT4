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
int fmpz_mpoly_is_divisible(const fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lead_term, lead_term2;
    const char *vars[] = {"x", "y", "z"};

    char *poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    flint_fprintf(stderr, "Reducing: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_init(lead_term, ctx);
    fmpz_mpoly_init(lead_term2, ctx);
    fmpz_mpoly_leadterm(lead_term, poly, ctx);

    for (slong i = 0; i < vec->length; i++) {
        fmpz_mpoly_t poly_i;
        fmpz_mpoly_init_set(poly_i, fmpz_mpoly_vec_entry(vec, i), ctx);
        fmpz_mpoly_leadterm(lead_term2, poly_i, ctx);

        if (lead_reduction) {
            while (fmpz_mpoly_is_divisible(lead_term, lead_term2, ctx)) {
                reduce_by_matching_term(poly, lead_term, poly_i, ctx);
                if (fmpz_mpoly_is_zero(poly, ctx)) {
                    break;
                }
                fmpz_mpoly_leadterm(lead_term, poly, ctx);
            }
        } else {
            slong n_terms = fmpz_mpoly_length(poly, ctx);
            for (slong j = 0; j < n_terms; j++) {
                fmpz_mpoly_t term_j;
                fmpz_mpoly_init(term_j, ctx);
                fmpz_mpoly_get_term(term_j, poly, j, ctx);
                if (fmpz_mpoly_is_divisible(term_j, lead_term2, ctx)) {
                    reduce_by_matching_term(poly, term_j, poly_i, ctx);
                }
                fmpz_mpoly_clear(term_j, ctx);
            }
        }
        fmpz_mpoly_clear(poly_i, ctx);
    }

    fmpz_mpoly_clear(lead_term, ctx);
    fmpz_mpoly_clear(lead_term2, ctx);

    poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    flint_fprintf(stderr, "Reduced: %s\n", poly_str);
    flint_free(poly_str);
}
