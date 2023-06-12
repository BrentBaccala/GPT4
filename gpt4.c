#include <stdio.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_set(lt, poly, ctx);
    fmpz_mpoly_truncate(lt, 1, ctx);
}

void construct_s_pair(fmpz_mpoly_t spair, const fmpz_mpoly_t poly1,
                      const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lt1, lt2, gcd_lt, temp1, temp2;

    fmpz_mpoly_init(lt1, ctx);
    fmpz_mpoly_init(lt2, ctx);
    fmpz_mpoly_init(gcd_lt, ctx);
    fmpz_mpoly_init(temp1, ctx);
    fmpz_mpoly_init(temp2, ctx);

    fmpz_mpoly_leadterm(lt1, poly1, ctx);
    fmpz_mpoly_leadterm(lt2, poly2, ctx);
    fmpz_mpoly_gcd(gcd_lt, lt1, lt2, ctx);

    fmpz_mpoly_mul(temp1, poly1, lt2, ctx);
    fmpz_mpoly_mul(temp2, poly2, lt1, ctx);
    fmpz_mpoly_div(temp1, temp1, gcd_lt, ctx);
    fmpz_mpoly_div(temp2, temp2, gcd_lt, ctx);

    fmpz_mpoly_sub(spair, temp1, temp2, ctx);

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd_lt, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    int done;
    slong i;
    fmpz_mpoly_t lt_poly, lt_vec, gcd_lt, temp;

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_vec, ctx);
    fmpz_mpoly_init(gcd_lt, ctx);
    fmpz_mpoly_init(temp, ctx);

    do {
        done = 1;
        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_leadterm(lt_poly, poly, ctx);
            fmpz_mpoly_leadterm(lt_vec, fmpz_mpoly_vec_entry(vec, i), ctx);

            if (fmpz_mpoly_divides(gcd_lt, lt_poly, lt_vec, ctx)) {
                fmpz_mpoly_mul(temp, fmpz_mpoly_vec_entry(vec, i), gcd_lt, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                done = 0;
            }
        }
    } while (!done);

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(lt_vec, ctx);
    fmpz_mpoly_clear(gcd_lt, ctx);
    fmpz_mpoly_clear(temp, ctx);
}

void buchberger_naive(fmpz_mpoly_vec_t basis, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx) {
    slong i, j;
    fmpz_mpoly_vec_t spairs;
    fmpz_mpoly_t spair;

    fmpz_mpoly_vec_init(spairs, 0, ctx);
    fmpz_mpoly_vec_set(basis, generators, ctx);
    fmpz_mpoly_init(spair, ctx);

    for (i = 0; i < basis->length; i++) {
        for (j = i + 1; j < basis->length; j++) {
            construct_s_pair(spair, fmpz_mpoly_vec_entry(basis, i), fmpz_mpoly_vec_entry(basis, j), ctx);
            if (!fmpz_mpoly_is_zero(spair, ctx)) {
                fmpz_mpoly_vec_append(spairs, spair, ctx);
            }
        }
    }

    for (i = 0; i < spairs->length; i++) {
        reduce_by_vector(fmpz_mpoly_vec_entry(spairs, i), basis, ctx);
        if (!fmpz_mpoly_is_zero(fmpz_mpoly_vec_entry(spairs, i), ctx)) {
            fmpz_mpoly_vec_append(basis, fmpz_mpoly_vec_entry(spairs, i), ctx);

            for (j = 0; j < basis->length - 1; j++) {
                construct_s_pair(spair, fmpz_mpoly_vec_entry(basis, j), fmpz_mpoly_vec_entry(basis, basis->length - 1), ctx);
                if (!fmpz_mpoly_is_zero(spair, ctx)) {
                    fmpz_mpoly_vec_append(spairs, spair, ctx);
                }
            }
        }
    }

    fmpz_mpoly_clear(spair, ctx);
    fmpz_mpoly_vec_clear(spairs, ctx);
}

void test_case_1() {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_vec_t generators, basis;
    const char *vars[] = {"x", "y", "z"};

    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);
    fmpz_mpoly_vec_init(generators, 2, ctx);
    fmpz_mpoly_vec_init(basis, 0, ctx);

    fmpz_mpoly_set_str_pretty(fmpz_mpoly_vec_entry(generators, 0), "2*x+3*y+4*z-5", vars, ctx);
    fmpz_mpoly_set_str_pretty(fmpz_mpoly_vec_entry(generators, 1), "3*x+4*y+5*z-2", vars, ctx);

    buchberger_naive(basis, generators, ctx);

    for (slong i = 0; i < basis->length; i++) {
        char *poly_str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(basis, i), vars, ctx);
        printf("%s\n", poly_str);
        flint_free(poly_str);
    }

    fmpz_mpoly_vec_clear(generators, ctx);
    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

int main() {
    test_case_1();
    return 0;
}
