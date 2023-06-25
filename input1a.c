#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t poly_vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);
void reduce_by_matching_term(fmpz_mpoly_t poly, const fmpz_mpoly_t matching_term, const fmpz_mpoly_t reducer, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t poly_vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lt_poly, lt_vec_entry;
    const char *var_names[] = {"x", "y", "z"};

    char *poly_str = fmpz_mpoly_get_str_pretty(poly, var_names, ctx);
    flint_fprintf(stderr, "Input poly: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_vec_entry, ctx);

    while (!fmpz_mpoly_is_zero(poly, ctx)) {
        fmpz_mpoly_leadterm(lt_poly, poly, ctx);

        int found_matching_term= 0;
        for (slong i = 0; i < poly_vec->length; ++i) {
            fmpz_mpoly_t vec_entry;
            fmpz_mpoly_init(vec_entry, ctx);
            fmpz_mpoly_set(vec_entry, fmpz_mpoly_vec_entry(poly_vec, i), ctx);

            fmpz_mpoly_leadterm(lt_vec_entry, vec_entry, ctx);

            if (lead_reduction) {
                if (fmpz_mpoly_is_divisible(lt_poly, lt_vec_entry, ctx)) {
                    reduce_by_matching_term(poly, lt_poly, vec_entry, ctx);
                    found_matching_term = 1;
                    break;
                }
            } else {
                fmpz_mpoly_t temp_poly;
                fmpz_mpoly_init(temp_poly, ctx);
                fmpz_mpoly_set(temp_poly, poly, ctx);

                for (slong term_idx = 0; term_idx < fmpz_mpoly_length(temp_poly, ctx); ++term_idx) {
                    fmpz_mpoly_t term;
                    fmpz_mpoly_init(term, ctx);
                    fmpz_mpoly_set(term, temp_poly, ctx);
                    fmpz_mpoly_truncate(term, term_idx + 1, ctx);
                    fmpz_mpoly_sub(term, temp_poly, term, ctx);
                    fmpz_mpoly_leadterm(term, term, ctx);

                    if (fmpz_mpoly_is_divisible(term, lt_vec_entry, ctx)) {
                        reduce_by_matching_term(poly, term, vec_entry, ctx);
                        found_matching_term = 1;
                    }

                    fmpz_mpoly_clear(term, ctx);
                }

                fmpz_mpoly_clear(temp_poly, ctx);
            }

            fmpz_mpoly_clear(vec_entry, ctx);
            if (found_matching_term) {
                break;
            }
        }

        if (!found_matching_term) {
            break;
        }
    }

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(lt_vec_entry, ctx);

    poly_str = fmpz_mpoly_get_str_pretty(poly, var_names, ctx);
    flint_fprintf(stderr, "Reduced poly: %s\n", poly_str);
    flint_free(poly_str);
}
