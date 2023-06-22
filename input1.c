#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lt_poly, lt_veci, gcd, temp;
    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_veci, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    const char *vars[] = {"x", "y", "z"};
    char *poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    flint_fprintf(stderr, "Input/output polynomial: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_leadterm(lt_poly, poly, ctx);

    for (slong i = 0; i < vec->length; i++) {
        fmpz_mpoly_t veci;
        fmpz_mpoly_init(veci, ctx);
        fmpz_mpoly_set(veci, fmpz_mpoly_vec_entry(vec, i), ctx);
        fmpz_mpoly_leadterm(lt_veci, veci, ctx);

        int done = 0;
        while (!done) {
            if (lead_reduction) {
                if (fmpz_mpoly_divides(gcd, lt_poly, lt_veci, ctx)) {
                    char *lt_str = fmpz_mpoly_get_str_pretty(lt_poly, vars, ctx);
                    flint_fprintf(stderr, "Matching term: %s\n", lt_str);
                    flint_free(lt_str);

                    fmpz_mpoly_mul(temp, veci, gcd, ctx);
                    fmpz_mpoly_sub(poly, poly, temp, ctx);
                    fmpz_mpoly_leadterm(lt_poly, poly, ctx);
                } else {
                    done = 1;
                }
            } else {
                done = 1;
                for (slong j = 0; j < fmpz_mpoly_length(poly, ctx); j++) {
                    fmpz_mpoly_t term_j;
                    fmpz_mpoly_init(term_j, ctx);
                    fmpz_mpoly_set(term_j, poly, ctx);
                    fmpz_mpoly_truncate(term_j, j + 1, ctx);

                    fmpz_mpoly_leadterm(term_j, term_j, ctx);
                    if (fmpz_mpoly_divides(gcd, term_j, lt_veci, ctx)) {
                        char *lt_str = fmpz_mpoly_get_str_pretty(term_j, vars, ctx);
                        flint_fprintf(stderr, "Matching term: %s\n", lt_str);
                        flint_free(lt_str);

                        fmpz_mpoly_mul(temp, veci, gcd, ctx);
                        fmpz_mpoly_sub(poly, poly, temp, ctx);
                    }
                    fmpz_mpoly_clear(term_j, ctx);
                }
            }

            if (fmpz_mpoly_is_zero(poly, ctx)) {
                break;
            }
        }
        fmpz_mpoly_clear(veci, ctx);
    }

    poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    flint_fprintf(stderr, "Reduced polynomial: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(lt_veci, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}
