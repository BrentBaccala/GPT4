#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t res, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t res, const fmpz_mpoly_t a, const fmpz_mpoly_t b, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx);
void fmpz_mpoly_get_term(fmpz_mpoly_t res, const fmpz_mpoly_t poly, slong n, const fmpz_mpoly_ctx_t ctx);
void reduce_by_matching_term(fmpz_mpoly_t poly, const fmpz_mpoly_t term, const fmpz_mpoly_t reducer, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    slong i, j;
    fmpz_mpoly_t lt_poly, lt_vec_poly;
    const char *vars[] = {"x", "y", "z"};
    char *poly_str;

    // Log input/output polynomial
    poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    flint_fprintf(stderr, "Input/output polynomial: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_vec_poly, ctx);

    while (!fmpz_mpoly_is_zero(poly, ctx)) {
        fmpz_mpoly_leadterm(lt_poly, poly, ctx);

        int reduced = 0;
        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_t vec_poly;
            fmpz_mpoly_init(vec_poly, ctx);
            fmpz_mpoly_set(vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);

            fmpz_mpoly_leadterm(lt_vec_poly, vec_poly, ctx);

            if (lead_reduction) {
                if (fmpz_mpoly_is_divisible(lt_poly, lt_vec_poly, ctx)) {
                    reduce_by_matching_term(poly, lt_poly, vec_poly, ctx);
                    reduced = 1;
                    break;
                }
            } else {
                for (j = 0; j < fmpz_mpoly_length(poly, ctx); j++) {
                    fmpz_mpoly_t term;
                    fmpz_mpoly_init(term, ctx);
                    fmpz_mpoly_get_term(term, poly, j, ctx);

                    if (fmpz_mpoly_is_divisible(term, lt_vec_poly, ctx)) {
                        reduce_by_matching_term(poly, term, vec_poly, ctx);
                    }
                    
                    fmpz_mpoly_clear(term, ctx);
                }
            }

            fmpz_mpoly_clear(vec_poly, ctx);
        }

        if (!reduced) {
            break;
        }
    }

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(lt_vec_poly, ctx);

    // Log reduced polynomial
    poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    flint_fprintf(stderr, "Reduced polynomial: %s\n", poly_str);
    flint_free(poly_str);
}
