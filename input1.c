#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t, const fmpz_mpoly_t, fmpz_mpoly_ctx_t);
void construct_s_pair(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t, fmpz_mpoly_ctx_t);
void reduce_by_vector(fmpz_mpoly_t, const fmpz_mpoly_vec_t, int, fmpz_mpoly_ctx_t);
void buchberger_naive(fmpz_mpoly_vec_t, const fmpz_mpoly_vec_t, fmpz_mpoly_ctx_t);
void buchberger_reduced(fmpz_mpoly_vec_t, const fmpz_mpoly_vec_t, fmpz_mpoly_ctx_t);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, fmpz_mpoly_ctx_t ctx) {
    const char *vars[] = {"x", "y", "z"};
    char *poly_str;
    slong i, j;
    fmpz_mpoly_t leadterm_poly, vec_leadterm, term, gcd, temp;

    fmpz_mpoly_init(leadterm_poly, ctx);
    fmpz_mpoly_init(vec_leadterm, ctx);
    fmpz_mpoly_init(term, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    flint_fprintf(stderr, "Input polynomial: %s\n", poly_str);
    flint_free(poly_str);

    while (!fmpz_mpoly_is_zero(poly, ctx)) {
        fmpz_mpoly_leadterm(leadterm_poly, poly, ctx);

        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_t vec_poly;

            fmpz_mpoly_init(vec_poly, ctx);
            fmpz_mpoly_set(vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);

            fmpz_mpoly_leadterm(vec_leadterm, vec_poly, ctx);

            if (lead_reduction) {
                if (fmpz_mpoly_divides(leadterm_poly, leadterm_poly, vec_leadterm, ctx)) {
                    fmpz_mpoly_gcd(gcd, leadterm_poly, vec_leadterm, ctx);
                    fmpz_mpoly_mul(temp, vec_poly, gcd, ctx);
                    fmpz_mpoly_sub(poly, poly, temp, ctx);

                    if (fmpz_mpoly_is_zero(poly, ctx)) {
                        break;
                    }

                    fmpz_mpoly_leadterm(leadterm_poly, poly, ctx);
                }
            } else {
                for (j = 0; j < fmpz_mpoly_length(poly, ctx); j++) {
                    fmpz_mpoly_get_term(term, poly, j, ctx);

                    if (fmpz_mpoly_divides(term, term, vec_leadterm, ctx)) {
                        fmpz_mpoly_gcd(gcd, term, vec_leadterm, ctx);
                        fmpz_mpoly_mul(temp, vec_poly, gcd, ctx);
                        fmpz_mpoly_sub(poly, poly, temp, ctx);
                    }
                }
            }

            fmpz_mpoly_clear(vec_poly, ctx);
        }

        if (!lead_reduction) {
            break;
        }
    }

    poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    flint_fprintf(stderr, "Reduced polynomial: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_clear(leadterm_poly, ctx);
    fmpz_mpoly_clear(vec_leadterm, ctx);
    fmpz_mpoly_clear(term, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}
