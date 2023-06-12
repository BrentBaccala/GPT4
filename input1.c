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
    slong i;
    int is_reducible;
    fmpz_mpoly_t leadterm, divisor, gcd, temp;
    
    fmpz_mpoly_init(leadterm, ctx);
    fmpz_mpoly_init(divisor, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    do {
        is_reducible = 0;
        fmpz_mpoly_leadterm(leadterm, poly, ctx);

        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_set(divisor, fmpz_mpoly_vec_entry(vec, i), ctx);
            fmpz_mpoly_leadterm(gcd, divisor, ctx);

            if (fmpz_mpoly_divides(gcd, leadterm, gcd, ctx)) {
                is_reducible = 1;
                fmpz_mpoly_divexact(temp, leadterm, gcd, ctx);
                fmpz_mpoly_mul(temp, temp, divisor, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                break;
            }
        }
    } while (is_reducible && !fmpz_mpoly_is_zero(poly, ctx));

    fmpz_mpoly_clear(leadterm, ctx);
    fmpz_mpoly_clear(divisor, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}
