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

void construct_s_pair(const fmpz_mpoly_t a, const fmpz_mpoly_t b,
                      fmpz_mpoly_t s_pair, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t a_lt, b_lt, gcd;

    fmpz_mpoly_init(a_lt, ctx);
    fmpz_mpoly_init(b_lt, ctx);
    fmpz_mpoly_init(gcd, ctx);

    fmpz_mpoly_leadterm(a_lt, a, ctx);
    fmpz_mpoly_leadterm(b_lt, b, ctx);
    fmpz_mpoly_gcd(gcd, a_lt, b_lt, ctx);

    fmpz_mpoly_mul(a_lt, a, b_lt, ctx);
    fmpz_mpoly_mul(b_lt, b, a_lt, ctx);

    fmpz_mpoly_div(a_lt, a_lt, gcd, ctx);
    fmpz_mpoly_div(b_lt, b_lt, gcd, ctx);

    fmpz_mpoly_sub(s_pair, a_lt, b_lt, ctx);

    fmpz_mpoly_clear(a_lt, ctx);
    fmpz_mpoly_clear(b_lt, ctx);
    fmpz_mpoly_clear(gcd, ctx);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    int found;
    slong i;
    fmpz_mpoly_t lt, gcd, temp;

    fmpz_mpoly_init(lt, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    do {
        found = 0;
        for (i = 0; i < vec->length && !fmpz_mpoly_is_zero(poly, ctx); i++) {
            fmpz_mpoly_leadterm(lt, fmpz_mpoly_vec_entry(vec, i), ctx);
            fmpz_mpoly_gcd(gcd, lt, poly, ctx);

            if (fmpz_mpoly_equal(gcd, lt, ctx)) {
                found = 1;
                fmpz_mpoly_div(temp, poly, gcd, ctx);
                fmpz_mpoly_mul(temp, temp, fmpz_mpoly_vec_entry(vec, i), ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
            }
        }
    } while (found && !fmpz_mpoly_is_zero(poly, ctx));

    fmpz_mpoly_clear(lt, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}

void buchberger_naive(fmpz_mpoly_vec_t basis, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx) {
    slong i, j;
    fmpz_mpoly_t s_pair, temp;
    fmpz_mpoly_vec_t s_pairs;

    fmpz_mpoly_init(s_pair, ctx);
    fmpz_mpoly_init(temp, ctx);
    fmpz_mpoly_vec_init(s_pairs, 0, ctx);

    fmpz_mpoly_vec_set(basis, generators, ctx);

    for (i = 0; i < basis->length; i++) {
        for (j = i + 1; j < basis->length; j++) {
            construct_s_pair(fmpz_mpoly_vec_entry(basis, i), fmpz_mpoly_vec_entry(basis, j), s_pair, ctx);
            if (!fmpz_mpoly_is_zero(s_pair, ctx)) {
                fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
            }
        }
    }

    for (i = 0; i < s_pairs->length; i++) {
        reduce_by_vector(s_pair, basis, ctx);
        if (!fmpz_mpoly_is_zero(s_pair, ctx)) {
            fmpz_mpoly_vec_append(basis, s_pair, ctx);
            for (j = 0; j < basis->length - 1; j++) {
                construct_s_pair(fmpz_mpoly_vec_entry(basis, j), s_pair, temp, ctx);
                if (!fmpz_mpoly_is_zero(temp, ctx)) {
                    fmpz_mpoly_vec_append(s_pairs, temp, ctx);
                }
            }
        }
    }

    fmpz_mpoly_clear(s_pair, ctx);
    fmpz_mpoly_clear(temp, ctx);
    fmpz_mpoly_vec_clear(s_pairs, ctx);
}

void test_case_1() {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_vec_t generators, basis;
    fmpz_mpoly_t poly1, poly2;
    const char* vars[] = {"x", "y", "z"};

    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);

    fmpz_mpoly_vec_init(generators, 0, ctx);
    fmpz_mpoly_vec_init(basis, 0, ctx);

    fmpz_mpoly_init(poly1, ctx);
    fmpz_mpoly_init(poly2, ctx);

    fmpz_mpoly_set_str_pretty(poly1, "2*x + 3*y + 4*z - 5", vars, ctx);
    fmpz_mpoly_set_str_pretty(poly2, "3*x + 4*y + 5*z - 2", vars, ctx);

    fmpz_mpoly_vec_append(generators, poly1, ctx);
    fmpz_mpoly_vec_append(generators, poly2, ctx);

    buchberger_naive(basis, generators, ctx);

    printf("Basis:\n");
    for (slong i = 0; i < basis->length; i++) {
        char* basis_str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(basis, i), vars, ctx);
        printf("%ld: %s\n", i, basis_str);
        flint_free(basis_str);
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
