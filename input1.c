#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t output, const fmpz_mpoly_t input, const fmpz_mpoly_ctx_t ctx);

void construct_s_pair(fmpz_mpoly_t output, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t io_poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);

void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);

void buchberger_reduced(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t io_poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    const char *var_names[3] = {"x", "y", "z"};
    char *poly_str = fmpz_mpoly_get_str_pretty(io_poly, var_names, ctx);
    flint_fprintf(stderr, "Input/output polynomial: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_t lt_io_poly;
    fmpz_mpoly_init(lt_io_poly, ctx);
    fmpz_mpoly_leadterm(lt_io_poly, io_poly, ctx);

    slong i, j;
    for (i = 0; i < vec->length; i++) {
        fmpz_mpoly_t poly_vec, lt_poly_vec;
        fmpz_mpoly_init(poly_vec, ctx);
        fmpz_mpoly_set(poly_vec, fmpz_mpoly_vec_entry(vec, i), ctx);
        fmpz_mpoly_init(lt_poly_vec, ctx);
        fmpz_mpoly_leadterm(lt_poly_vec, poly_vec, ctx);

        if (lead_reduction) {
            if (fmpz_mpoly_divides(lt_io_poly, lt_io_poly, lt_poly_vec, ctx)) {
                fmpz_mpoly_t temp_poly;
                fmpz_mpoly_init(temp_poly, ctx);
                fmpz_mpoly_mul(temp_poly, poly_vec, lt_io_poly, ctx);
                fmpz_mpoly_sub(io_poly, io_poly, temp_poly, ctx);
                fmpz_mpoly_clear(temp_poly, ctx);
                fmpz_mpoly_leadterm(lt_io_poly, io_poly, ctx);
            }
        } else {
            for (j = 0; j < fmpz_mpoly_length(io_poly, ctx); j++) {
                fmpz_mpoly_t term_j, lt_term_j;
                fmpz_mpoly_init(term_j, ctx);
                fmpz_mpoly_set(term_j, io_poly, ctx);
                fmpz_mpoly_truncate(term_j, j + 1, ctx);
                fmpz_mpoly_init(lt_term_j, ctx);
                fmpz_mpoly_leadterm(lt_term_j, term_j, ctx);

                if (fmpz_mpoly_divides(lt_term_j, lt_term_j, lt_poly_vec, ctx)) {
                    poly_str = fmpz_mpoly_get_str_pretty(lt_term_j, var_names, ctx);
                    flint_fprintf(stderr, "Matching term: %s\n", poly_str);
                    flint_free(poly_str);
                    
                    fmpz_mpoly_t temp_poly;
                    fmpz_mpoly_init(temp_poly, ctx);
                    fmpz_mpoly_mul(temp_poly, poly_vec, lt_term_j, ctx);
                    fmpz_mpoly_sub(io_poly, io_poly, temp_poly, ctx);
                    fmpz_mpoly_clear(temp_poly, ctx);
                }

                fmpz_mpoly_clear(term_j, ctx);
                fmpz_mpoly_clear(lt_term_j, ctx);
            }
        }

        fmpz_mpoly_clear(poly_vec, ctx);
        fmpz_mpoly_clear(lt_poly_vec, ctx);

        if (fmpz_mpoly_is_zero(io_poly, ctx)) {
            break;
        }
    }

    poly_str = fmpz_mpoly_get_str_pretty(io_poly, var_names, ctx);
    flint_fprintf(stderr, "Result: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_clear(lt_io_poly, ctx);
}
