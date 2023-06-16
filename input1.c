#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);

void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx);

void buchberger_reduced(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t inout, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    const char *vars[] = {"x", "y", "z"};
    char *inout_str = fmpz_mpoly_get_str_pretty(inout, vars, ctx);
    flint_fprintf(stderr, "Starting reduction with polynomial: %s\n", inout_str);
    free(inout_str);

    fmpz_mpoly_t lt_inout, lt_vec, gcd, temp;
    fmpz_mpoly_init(lt_inout, ctx);
    fmpz_mpoly_init(lt_vec, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    int changed;
    do {
        changed = 0;
        fmpz_mpoly_leadterm(lt_inout, inout, ctx);
        for (slong i = 0; i < vec->length; ++i) {
            fmpz_mpoly_t vec_poly;
            fmpz_mpoly_init(vec_poly, ctx);
            fmpz_mpoly_set(vec_poly, fmpz_mpoly_vec_entry(vec, i), ctx);

            fmpz_mpoly_leadterm(lt_vec, vec_poly, ctx);
            if (fmpz_mpoly_divides(gcd, lt_inout, lt_vec, ctx)) {
                flint_fprintf(stderr, "Divisor found\n");
                fmpz_mpoly_mul(temp, vec_poly, gcd, ctx);
                fmpz_mpoly_sub(inout, inout, temp, ctx);
                fmpz_mpoly_leadterm(lt_inout, inout, ctx);
                changed = 1;
                break;
            }
            fmpz_mpoly_clear(vec_poly, ctx);
        }
    } while (changed && !fmpz_mpoly_is_zero(inout, ctx));

    inout_str = fmpz_mpoly_get_str_pretty(inout, vars, ctx);
    flint_fprintf(stderr, "Reduction result: %s\n", inout_str);
    free(inout_str);

    fmpz_mpoly_clear(lt_inout, ctx);
    fmpz_mpoly_clear(lt_vec, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}
