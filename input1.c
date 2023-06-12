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
    char *poly_str;
    const char *var_names[] = {"x", "y", "z"};
    fmpz_mpoly_t lt_poly, lt_vec_poly, temp_poly;
    slong i;

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_vec_poly, ctx);
    fmpz_mpoly_init(temp_poly, ctx);

    poly_str = fmpz_mpoly_get_str_pretty(poly, var_names, ctx);
    flint_printf("Starting reduction with polynomial %s\n", poly_str);
    flint_free(poly_str);

    int reduced;
    do {
        reduced = 0;
        fmpz_mpoly_leadterm(lt_poly, poly, ctx);

        for (i = 0; i < vec->length && !fmpz_mpoly_is_zero(poly, ctx); i++) {
            fmpz_mpoly_set(lt_vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);
            fmpz_mpoly_leadterm(lt_vec_poly, lt_vec_poly, ctx);

            if (fmpz_mpoly_divides(temp_poly, lt_poly, lt_vec_poly, ctx)) {
                reduced = 1;
                poly_str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(vec, i), var_names, ctx);
                flint_printf("Divisible by %s\n", poly_str);
                flint_free(poly_str);

                fmpz_mpoly_mul(temp_poly, temp_poly, fmpz_mpoly_vec_entry(vec, i), ctx);
                fmpz_mpoly_sub(poly, poly, temp_poly, ctx);
                fmpz_mpoly_leadterm(lt_poly, poly, ctx);
            }
        }
    } while (reduced && !fmpz_mpoly_is_zero(poly, ctx));

    poly_str = fmpz_mpoly_get_str_pretty(poly, var_names, ctx);
    flint_printf("Reduction result: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(lt_vec_poly, ctx);
    fmpz_mpoly_clear(temp_poly, ctx);
}
