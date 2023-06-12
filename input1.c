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
    int reduced;
    fmpz_mpoly_t lt_poly, lt_vec_poly, gcd, temp;

    const char *vars[] = {"x", "y", "z"};
    char *poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    fprintf(stderr, "Starting reduction: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_vec_poly, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    do {
        reduced = 0;
        fmpz_mpoly_leadterm(lt_poly, poly, ctx);

        for (slong i = 0; i < vec->length; i++) {
            fmpz_mpoly_t vec_poly;
            fmpz_mpoly_init_set(vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);
            fmpz_mpoly_leadterm(lt_vec_poly, vec_poly, ctx);

            if (fmpz_mpoly_divides(gcd, lt_poly, lt_vec_poly, ctx)) {
                fprintf(stderr, "Reducing by: %s\n", fmpz_mpoly_get_str_pretty(vec_poly, vars, ctx));
                fmpz_mpoly_mul(temp, vec_poly, gcd, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                fmpz_mpoly_leadterm(lt_poly, poly, ctx);
                reduced = 1;
                break;
            }

            fmpz_mpoly_clear(vec_poly, ctx);
        }
    } while (reduced && !fmpz_mpoly_is_zero(poly, ctx));

    poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    fprintf(stderr, "Reduced: %s\n", poly_str);
    flint_free(poly_str);

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(lt_vec_poly, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}
