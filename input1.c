#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

void construct_s_pair(fmpz_mpoly_t spair, const fmpz_mpoly_t p1, const fmpz_mpoly_t p2, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);

void buchberger_naive(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);

void buchberger_reduced(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    const char *vars[] = {"x", "y", "z"};
    char *s;
    slong i, j;
    fmpz_mpoly_t lt, lt_vec, tmp;
    int found_matching_term;

    s = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    flint_fprintf(stderr, "Input/output polynomial: %s\n", s);
    flint_free(s);

    fmpz_mpoly_init(lt, ctx);
    fmpz_mpoly_init(lt_vec, ctx);
    fmpz_mpoly_init(tmp, ctx);

    fmpz_mpoly_leadterm(lt, poly, ctx);

    while (!fmpz_mpoly_is_zero(poly, ctx)) {
        found_matching_term = 0;
        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_set(tmp, fmpz_mpoly_vec_entry(vec, i), ctx);
            fmpz_mpoly_leadterm(lt_vec, tmp, ctx);

            if (lead_reduction) {
                if (fmpz_mpoly_divides(tmp, lt, lt_vec, ctx)) {
                    found_matching_term = 1;
                    s = fmpz_mpoly_get_str_pretty(lt_vec, vars, ctx);
                    flint_fprintf(stderr, "Matching term: %s\n", s);
                    flint_free(s);
                }
            } else {
                for (j = 0; j < fmpz_mpoly_length(poly, ctx); j++) {
                    fmpz_mpoly_get_term_monomial(lt, poly, j, ctx);
                    if (fmpz_mpoly_divides(tmp, lt, lt_vec, ctx)) {
                        found_matching_term = 1;
                        s = fmpz_mpoly_get_str_pretty(lt, vars, ctx);
                        flint_fprintf(stderr, "Matching term: %s\n", s);
                        flint_free(s);
                    }
                }
            }

            if (found_matching_term) {
                fmpz_mpoly_gcd(lt, lt, lt_vec, ctx);
                fmpz_mpoly_divexact(tmp, lt, lt_vec, ctx);
                fmpz_mpoly_mul(tmp, tmp, fmpz_mpoly_vec_entry(vec, i), ctx);
                fmpz_mpoly_sub(poly, poly, tmp, ctx);

                if (lead_reduction) {
                    fmpz_mpoly_leadterm(lt, poly, ctx);
                } else {
                    break;
                }
            }
        }

        if (!found_matching_term) {
            break;
        }
    }

    s = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    flint_fprintf(stderr, "Reduced polynomial: %s\n", s);
    flint_free(s);

    fmpz_mpoly_clear(lt, ctx);
    fmpz_mpoly_clear(lt_vec, ctx);
    fmpz_mpoly_clear(tmp, ctx);
}
