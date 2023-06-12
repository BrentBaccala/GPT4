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

void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t p1, const fmpz_mpoly_t p2, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lt1, lt2, gcd, temp1, temp2;

    fmpz_mpoly_init(lt1, ctx);
    fmpz_mpoly_init(lt2, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp1, ctx);
    fmpz_mpoly_init(temp2, ctx);

    fmpz_mpoly_leadterm(lt1, p1, ctx);
    fmpz_mpoly_leadterm(lt2, p2, ctx);
    fmpz_mpoly_gcd(gcd, lt1, lt2, ctx);

    fmpz_mpoly_mul(temp1, p1, lt2, ctx);
    fmpz_mpoly_divexact(temp1, temp1, gcd, ctx);
    fmpz_mpoly_mul(temp2, p2, lt1, ctx);
    fmpz_mpoly_divexact(temp2, temp2, gcd, ctx);

    fmpz_mpoly_sub(s_pair, temp1, temp2, ctx);

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    slong i;
    int changed;
    fmpz_mpoly_t lt_poly, lt_vec, gcd, temp;

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_vec, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    do {
        changed = 0;
        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_leadterm(lt_poly, poly, ctx);
            fmpz_mpoly_leadterm(lt_vec, fmpz_mpoly_vec_entry(vec, i), ctx);

            if (fmpz_mpoly_divides(gcd, lt_poly, lt_vec, ctx)) {
                fmpz_mpoly_mul(temp, fmpz_mpoly_vec_entry(vec, i), gcd, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                changed = 1;
            }
        }
    } while (changed);

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(lt_vec, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}

void buchberger_naive(fmpz_mpoly_vec_t basis, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_vec_t s_pairs;
    fmpz_mpoly_t s_pair;
    slong i, j;

    fmpz_mpoly_vec_init(basis, 0, ctx);
    fmpz_mpoly_vec_init(s_pairs, 0, ctx);
    fmpz_mpoly_init(s_pair, ctx);

    for (i = 0; i < generators->length; i++) {
        fmpz_mpoly_vec_append(basis, fmpz_mpoly_vec_entry(generators, i), ctx);
    }

    for (i = 0; i < generators->length; i++) {
        for (j = i + 1; j < generators->length; j++) {
            construct_s_pair(s_pair, fmpz_mpoly_vec_entry(generators, i), fmpz_mpoly_vec_entry(generators, j), ctx);
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
                construct_s_pair(s_pair, fmpz_mpoly_vec_entry(basis, j), fmpz_mpoly_vec_entry(basis, basis->length - 1), ctx);
                if (!fmpz_mpoly_is_zero(s_pair, ctx)) {
                    fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
                }
            }
        }
    }

    fmpz_mpoly_clear(s_pair, ctx);
    fmpz_mpoly_vec_clear(s_pairs, ctx);
}

void test_case_1() {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_vec_t generators, basis;
    fmpz_mpoly_t p1, p2;
    const char* vars[] = {"x", "y", "z"};

    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);
    fmpz_mpoly_vec_init(generators, 0, ctx);
    fmpz_mpoly_vec_init(basis, 0, ctx);

    fmpz_mpoly_init(p1, ctx);
    fmpz_mpoly_init(p2, ctx);

    fmpz_mpoly_set_str_pretty(p1, "2*x + 3*y + 4*z - 5", vars, ctx);
    fmpz_mpoly_set_str_pretty(p2, "3*x + 4*y + 5*z - 2", vars, ctx);

    fmpz_mpoly_vec_append(generators, p1, ctx);
    fmpz_mpoly_vec_append(generators, p2, ctx);

    buchberger_naive(basis, generators, ctx);

    for (slong i = 0; i < basis->length; i++) {
        char *str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(basis, i), vars, ctx);
        printf("%s\n", str);
        flint_free(str);
    }

    fmpz_mpoly_clear(p1, ctx);
    fmpz_mpoly_clear(p2, ctx);
    fmpz_mpoly_vec_clear(generators, ctx);
    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

int main() {
    test_case_1();
    return 0;
}
