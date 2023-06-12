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

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t leadterm_poly, temp_poly, gcd_poly, vec_poly;

    fmpz_mpoly_init(leadterm_poly, ctx);
    fmpz_mpoly_init(temp_poly, ctx);
    fmpz_mpoly_init(gcd_poly, ctx);
    fmpz_mpoly_init(vec_poly, ctx);

    char **var_strs = (char *[]){"x", "y", "z"}; // You can modify this to match the number of variables in your polynomial ring
    char *poly_str = fmpz_mpoly_get_str_pretty(poly, var_strs, ctx);
    fprintf(stderr, "Input polynomial: %s\n", poly_str);
    flint_free(poly_str);

    int reduced;
    do {
        reduced = 0;
        for (slong i = 0; i < vec->length; i++) {
            fmpz_mpoly_set(vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);
            fmpz_mpoly_leadterm(leadterm_poly, vec_poly, ctx);
            
            if (fmpz_mpoly_divides(temp_poly, poly, leadterm_poly, ctx)) {
                fmpz_mpoly_gcd(gcd_poly, leadterm_poly, temp_poly, ctx);
                fmpz_mpoly_mul(temp_poly, vec_poly, gcd_poly, ctx);
                fmpz_mpoly_sub(poly, poly, temp_poly, ctx);

                poly_str = fmpz_mpoly_get_str_pretty(poly, var_strs, ctx);
                fprintf(stderr, "Reduced by vector polynomial %ld: %s\n", i, poly_str);
                flint_free(poly_str);

                reduced = 1;
                break;
            }
        }
    } while (reduced && !fmpz_mpoly_is_zero(poly, ctx));

    poly_str = fmpz_mpoly_get_str_pretty(poly, var_strs, ctx);
    fprintf(stderr, "Result after reduction: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_clear(leadterm_poly, ctx);
    fmpz_mpoly_clear(temp_poly, ctx);
    fmpz_mpoly_clear(gcd_poly, ctx);
    fmpz_mpoly_clear(vec_poly, ctx);
}
