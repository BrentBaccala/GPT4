#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t res, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_set(res, poly, ctx);
    fmpz_mpoly_truncate(res, 1, ctx);
}

void construct_s_pair(fmpz_mpoly_t spair, const fmpz_mpoly_t a, const fmpz_mpoly_t b, const fmpz_mpoly_ctx_t ctx) {
    fprintf(stderr, "Constructing s-pair\n");

    fmpz_mpoly_t a_lt, b_lt, gcd, a_copy, b_copy;
    fmpz_mpoly_init(a_lt, ctx);
    fmpz_mpoly_init(b_lt, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(a_copy, ctx);
    fmpz_mpoly_init(b_copy, ctx);

    fmpz_mpoly_leadterm(a_lt, a, ctx);
    fmpz_mpoly_leadterm(b_lt, b, ctx);
    fmpz_mpoly_gcd(gcd, a_lt, b_lt, ctx);

    fmpz_mpoly_mul(a_copy, a, b_lt, ctx);
    fmpz_mpoly_divexact(a_copy, a_copy, gcd, ctx);

    fmpz_mpoly_mul(b_copy, b, a_lt, ctx);
    fmpz_mpoly_divexact(b_copy, b_copy, gcd, ctx);

    fmpz_mpoly_sub(spair, a_copy, b_copy, ctx);

    fprintf(stderr, "Constructed s-pair\n");

    fmpz_mpoly_clear(a_lt, ctx);
    fmpz_mpoly_clear(b_lt, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(a_copy, ctx);
    fmpz_mpoly_clear(b_copy, ctx);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t poly_vec, const fmpz_mpoly_ctx_t ctx) {
    int changed;
    fmpz_mpoly_t tmp, gcd, tmp_poly;
    fmpz_mpoly_init(tmp, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(tmp_poly, ctx);
    do {
        changed = 0;
        for (slong i = 0; i < poly_vec->length; i++) {
            const fmpz_mpoly_t vec_poly = fmpz_mpoly_vec_entry(poly_vec, i);
            fmpz_mpoly_leadterm(tmp, vec_poly, ctx);
            if (fmpz_mpoly_divides(gcd, fmpz_mpoly_leadterm(poly, ctx), tmp, ctx)) {
                fmpz_mpoly_mul(tmp_poly, vec_poly, gcd, ctx);
                fmpz_mpoly_sub(poly, poly, tmp_poly, ctx);
                changed = 1;
            }
        }
    } while (changed);
    fmpz_mpoly_clear(tmp, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(tmp_poly, ctx);
}

void buchberger_naive(fmpz_mpoly_vec_t output_basis, const fmpz_mpoly_vec_t input_generators, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_vec_t basis, s_pairs;
    fmpz_mpoly_vec_init(basis, input_generators->length, ctx);
    fmpz_mpoly_vec_init(s_pairs, 0, ctx);

    fmpz_mpoly_vec_set(basis, input_generators, ctx);

    for (slong i = 0; i < basis->length; i++) {
        for (slong j = i + 1; j < basis->length; j++) {
            fmpz_mpoly_t s_pair;
            fmpz_mpoly_init(s_pair, ctx);

            construct_s_pair(s_pair, fmpz_mpoly_vec_entry(basis, i), fmpz_mpoly_vec_entry(basis, j), ctx);

            if (!fmpz_mpoly_is_zero(s_pair, ctx)) {
                fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
            }
            fmpz_mpoly_clear(s_pair, ctx);
        }
    }

    slong index = 0;
    while (index < s_pairs->length) {
        fmpz_mpoly_t reduced_s_pair;
        fmpz_mpoly_init(reduced_s_pair, ctx);
        fmpz_mpoly_set(reduced_s_pair, fmpz_mpoly_vec_entry(s_pairs, index), ctx);

        reduce_by_vector(reduced_s_pair, basis, ctx);

        if (!fmpz_mpoly_is_zero(reduced_s_pair, ctx)) {
            fmpz_mpoly_vec_append(basis, reduced_s_pair, ctx);

            for (slong i = 0; i < basis->length - 1; i++) {
                fmpz_mpoly_t new_s_pair;
                fmpz_mpoly_init(new_s_pair, ctx);
                construct_s_pair(new_s_pair, fmpz_mpoly_vec_entry(basis, i), fmpz_mpoly_vec_entry(basis, basis->length - 1), ctx);

                if (!fmpz_mpoly_is_zero(new_s_pair, ctx)) {
                    fmpz_mpoly_vec_append(s_pairs, new_s_pair, ctx);
                }
                fmpz_mpoly_clear(new_s_pair, ctx);
            }
        }
        fmpz_mpoly_clear(reduced_s_pair, ctx);
        index++;
    }

    fmpz_mpoly_vec_swap(output_basis, basis, ctx);
    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(s_pairs, ctx);
}

void test_case_1() {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);

    char *vars[] = {"x", "y", "z"};

    fmpz_mpoly_t f1, f2;
    fmpz_mpoly_init(f1, ctx);
    fmpz_mpoly_init(f2, ctx);
    fmpz_mpoly_set_str_pretty(f1, "2*x + 3*y + 4*z - 5", vars, ctx);
    fmpz_mpoly_set_str_pretty(f2, "3*x + 4*y + 5*z - 2", vars, ctx);

    fmpz_mpoly_vec_t generators, basis;
    fmpz_mpoly_vec_init(generators, 2, ctx);
    fmpz_mpoly_vec_init(basis, 0, ctx);
    fmpz_mpoly_vec_set_coeff(generators, 0, f1, ctx);
    fmpz_mpoly_vec_set_coeff(generators, 1, f2, ctx);

    buchberger_naive(basis, generators, ctx);

    for (slong i = 0; i < basis->length; i++) {
        char *str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(basis, i), vars, ctx);
        printf("%s\n", str);
        flint_free(str);
    }

    fmpz_mpoly_clear(f1, ctx);
    fmpz_mpoly_clear(f2, ctx);
    fmpz_mpoly_vec_clear(generators, ctx);
    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

int main() {
    test_case_1();
    return 0;
}
