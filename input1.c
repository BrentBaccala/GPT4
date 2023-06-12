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
    const char *vars[] = {"x", "y", "z"};
    fmpz_mpoly_t leadterm, gcd, temp;

    fmpz_mpoly_init(leadterm, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    fprintf(stderr, "Input: %s\n", fmpz_mpoly_get_str_pretty(poly, vars, ctx));

    int reduced;
    do {
        reduced = 0;
        fmpz_mpoly_leadterm(leadterm, poly, ctx);

        for (slong i = 0; i < vec->length && !reduced; i++) {
            fmpz_mpoly_t vec_poly;

            fmpz_mpoly_init(vec_poly, ctx);
            fmpz_mpoly_set(vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);

            if (fmpz_mpoly_divides(gcd, leadterm, vec_poly, ctx)) {
                fprintf(stderr, "Reducing by: %s\n", fmpz_mpoly_get_str_pretty(vec_poly, vars, ctx));

                fmpz_mpoly_mul(temp, vec_poly, gcd, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                fmpz_mpoly_leadterm(leadterm, poly, ctx);
                reduced = 1;
            }

            fmpz_mpoly_clear(vec_poly, ctx);
        }
    } while (reduced && !fmpz_mpoly_is_zero(poly, ctx));

    fprintf(stderr, "Reduced output: %s\n", fmpz_mpoly_get_str_pretty(poly, vars, ctx));

    fmpz_mpoly_clear(leadterm, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}
