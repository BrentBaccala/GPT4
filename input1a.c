#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t p1, const fmpz_mpoly_t p2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction,
                      const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);

void reduce_by_matching_term(fmpz_mpoly_t poly, const fmpz_mpoly_t term, const fmpz_mpoly_t to_reduce_with, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    const char *var_names[] = {"x", "y", "z"};
    char *poly_str = fmpz_mpoly_get_str_pretty(poly, var_names, ctx);
    flint_fprintf(stderr, "Reduce by vector: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_t lt_poly;
    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_leadterm(lt_poly, poly, ctx);

    slong i, j;
    for (i = 0; i < vec->length && !fmpz_mpoly_is_zero(poly, ctx); i++) {
        fmpz_mpoly_t lt_vec;
        fmpz_mpoly_init(lt_vec, ctx);
        fmpz_mpoly_t vec_poly;
        fmpz_mpoly_init(vec_poly, ctx);
        fmpz_mpoly_set(vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);
        fmpz_mpoly_leadterm(lt_vec, vec_poly, ctx);

        if (lead_reduction) {
            if (fmpz_mpoly_divides(lt_poly, lt_poly, lt_vec, ctx)) {
                reduce_by_matching_term(poly, lt_poly, vec_poly, ctx);
                fmpz_mpoly_leadterm(lt_poly, poly, ctx);
            }
        } else {
            for (j = 0; j < fmpz_mpoly_length(poly, ctx); j++) {
                fmpz_mpoly_t term;
                fmpz_mpoly_init(term, ctx);
                fmpz_mpoly_get_term_monomial(term, poly, j, ctx);

                if (fmpz_mpoly_divides(term, term, lt_vec, ctx)) {
                    reduce_by_matching_term(poly, term, vec_poly, ctx);
                }
                fmpz_mpoly_clear(term, ctx);
            }
        }

        fmpz_mpoly_clear(lt_vec, ctx);
        fmpz_mpoly_clear(vec_poly, ctx);
    }

    poly_str = fmpz_mpoly_get_str_pretty(poly, var_names, ctx);
    flint_fprintf(stderr, "Reduced by vector: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_clear(lt_poly, ctx);
}
