#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_set(lt, poly, ctx);
    fmpz_mpoly_truncate(lt, 1, ctx);
}

void construct_s_pair(fmpz_mpoly_t spair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lt1, lt2, gcd, temp1, temp2;

    fmpz_mpoly_init(lt1, ctx);
    fmpz_mpoly_init(lt2, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp1, ctx);
    fmpz_mpoly_init(temp2, ctx);

    fmpz_mpoly_leadterm(lt1, poly1, ctx);
    fmpz_mpoly_leadterm(lt2, poly2, ctx);
    fmpz_mpoly_gcd(gcd, lt1, lt2, ctx);

    fmpz_mpoly_mul(temp1, poly1, lt2, ctx);
    fmpz_mpoly_mul(temp2, poly2, lt1, ctx);
    fmpz_mpoly_divexact(temp1, temp1, gcd, ctx);
    fmpz_mpoly_divexact(temp2, temp2, gcd, ctx);

    fmpz_mpoly_sub(spair, temp1, temp2, ctx);

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    int changed;
    fmpz_mpoly_t lt_poly, lt_vec_poly, gcd, temp;

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_vec_poly, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    do {
        changed = 0;
        for (slong i = 0; i < vec->length; i++) {
            fmpz_mpoly_t vec_poly = fmpz_mpoly_vec_entry(vec, i);

            fmpz_mpoly_leadterm(lt_poly, poly, ctx);
            fmpz_mpoly_leadterm(lt_vec_poly, vec_poly, ctx);

            if (fmpz_mpoly_divides(gcd, lt_poly, lt_vec_poly, ctx)) {
                fmpz_mpoly_mul(temp, vec_poly, gcd, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                changed = 1;
            }
        }
    } while (changed);

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(lt_vec_poly, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}

void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_vec_t basis, spairs;
    fmpz_mpoly_t spair, reduced_spair;

    fmpz_mpoly_vec_init(basis, 0, ctx);
    fmpz_mpoly_vec_init(spairs, 0, ctx);
    fmpz_mpoly_init(spair, ctx);
    fmpz_mpoly_init(reduced_spair, ctx);

    fmpz_mpoly_vec_set(output, generators, ctx);

    for (slong i = 0; i < generators->length - 1; i++) {
        for (slong j = i + 1; j < generators->length; j++) {
            construct_s_pair(spair, fmpz_mpoly_vec_entry(generators, i), fmpz_mpoly_vec_entry(generators, j), ctx);
            fmpz_mpoly_vec_append(spairs, spair, ctx);
        }
    }

    for (slong i = 0; i < spairs->length; i++) {
        fmpz_mpoly_set(reduced_spair, fmpz_mpoly_vec_entry(spairs, i), ctx);
        reduce_by_vector(reduced_spair, basis, ctx);

        if (!fmpz_mpoly_is_zero(reduced_spair, ctx)) {
            fmpz_mpoly_vec_append(basis, reduced_spair, ctx);

            for(slong j = 0; j < basis->length - 1; j++) {
                construct_s_pair(spair, fmpz_mpoly_vec_entry(basis, j), reduced_spair, ctx);
                fmpz_mpoly_vec_append(spairs, spair, ctx);
            }
        }
    }

    fmpz_mpoly_vec_set(output, basis, ctx);

    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(spairs, ctx);
    fmpz_mpoly_clear(spair, ctx);
    fmpz_mpoly_clear(reduced_spair, ctx);
}

void test_case_1() {
    const char* vars[] = {"x", "y", "z"};
    fmpz_mpoly_t poly1, poly2;
    fmpz_mpoly_vec_t generators, basis;
    fmpz_mpoly_ctx_t ctx;

    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);

    fmpz_mpoly_init(poly1, ctx);
    fmpz_mpoly_init(poly2, ctx);

    fmpz_mpoly_set_str_pretty(poly1, "2*x + 3*y + 4*z - 5", vars, ctx);
    fmpz_mpoly_set_str_pretty(poly2, "3*x + 4*y + 5*z - 2", vars, ctx);

    fmpz_mpoly_vec_init(generators, 2, ctx);
    fmpz_mpoly_vec_init(basis, 0, ctx);

    fmpz_mpoly_vec_append(generators, poly1, ctx);
    fmpz_mpoly_vec_append(generators, poly2, ctx);

    buchberger_naive(basis, generators, ctx);

    for (slong i = 0; i < basis->length; i++) {
        char *poly_str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(basis, i), vars, ctx);
        printf("%s\n", poly_str);
        flint_free(poly_str);
    }

    fmpz_mpoly_clear(poly1, ctx);
    fmpz_mpoly_clear(poly2, ctx);
    fmpz_mpoly_vec_clear(generators, ctx);
    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

int main() {
    test_case_1();
    return 0;
}
