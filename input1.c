#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t res, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gens, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gens, const fmpz_mpoly_ctx_t ctx);


void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    slong i;
    int is_reducible;
    fmpz_mpoly_t lt_poly, lt_vec_poly, gcd_lt, temp;

    // Log input/output polynomial
    char *poly_str = fmpz_mpoly_get_str_pretty(poly, (const char *[]){"x", "y", "z"}, ctx);
    fprintf(stderr, "Input/output polynomial: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_vec_poly, ctx);
    fmpz_mpoly_init(gcd_lt, ctx);
    fmpz_mpoly_init(temp, ctx);

    do {
        is_reducible = 0;
        fmpz_mpoly_leadterm(lt_poly, poly, ctx);

        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_t vec_poly;
            fmpz_mpoly_init(vec_poly, ctx);
            fmpz_mpoly_set(vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);

            fmpz_mpoly_leadterm(lt_vec_poly, vec_poly, ctx);

            if (fmpz_mpoly_divides(gcd_lt, lt_poly, lt_vec_poly, ctx)) {
                is_reducible = 1;

                // Log the division
                char *vec_poly_str = fmpz_mpoly_get_str_pretty(vec_poly, (const char *[]){"x", "y", "z"}, ctx);
                fprintf(stderr, "Divides: %s\n", vec_poly_str);
                flint_free(vec_poly_str);

                fmpz_mpoly_mul(temp, vec_poly, gcd_lt, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);

                // Recompute the leading term
                fmpz_mpoly_leadterm(lt_poly, poly, ctx);
                break;
            }

            fmpz_mpoly_clear(vec_poly, ctx);
        }
    } while (is_reducible && !fmpz_mpoly_is_zero(poly, ctx));

    // Log the result of the reduction
    char *result_str = fmpz_mpoly_get_str_pretty(poly, (const char *[]){"x", "y", "z"}, ctx);
    fprintf(stderr, "Reduced polynomial: %s\n", result_str);
    flint_free(result_str);

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(lt_vec_poly, ctx);
    fmpz_mpoly_clear(gcd_lt, ctx);
    fmpz_mpoly_clear(temp, ctx);
}
