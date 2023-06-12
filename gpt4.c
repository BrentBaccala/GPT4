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

void construct_s_pair(fmpz_mpoly_t res, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lt1, lt2, lt_gcd, temp1, temp2;

    fmpz_mpoly_init(lt1, ctx);
    fmpz_mpoly_init(lt2, ctx);
    fmpz_mpoly_init(lt_gcd, ctx);
    fmpz_mpoly_init(temp1, ctx);
    fmpz_mpoly_init(temp2, ctx);

    fmpz_mpoly_leadterm(lt1, poly1, ctx);
    fmpz_mpoly_leadterm(lt2, poly2, ctx);
    fmpz_mpoly_gcd(lt_gcd, lt1, lt2, ctx);

    fmpz_mpoly_mul(temp1, poly1, lt2, ctx);
    fmpz_mpoly_mul(temp2, poly2, lt1, ctx);
    fmpz_mpoly_divexact(temp1, temp1, lt_gcd, ctx);
    fmpz_mpoly_divexact(temp2, temp2, lt_gcd, ctx);

    fmpz_mpoly_sub(res, temp1, temp2, ctx);

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(lt_gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    int changed;
    fmpz_mpoly_t lt_poly, lt_vec, temp, gcd;

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_vec, ctx);
    fmpz_mpoly_init(temp, ctx);
    fmpz_mpoly_init(gcd, ctx);

    do {
        changed = 0;
        fmpz_mpoly_leadterm(lt_poly, poly, ctx);

        for (slong i = 0; i < vec->length; i++) {
            fmpz_mpoly_t vec_entry;
            fmpz_mpoly_init(vec_entry, ctx);
            fmpz_mpoly_set(vec_entry, fmpz_mpoly_vec_entry(vec, i), ctx);
            fmpz_mpoly_leadterm(lt_vec, vec_entry, ctx);

            if (fmpz_mpoly_divides(gcd, lt_poly, lt_vec, ctx)) {
                fmpz_mpoly_mul(temp, vec_entry, gcd, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                fmpz_mpoly_leadterm(lt_poly, poly, ctx);
                changed = 1;
            }

            fmpz_mpoly_clear(vec_entry, ctx);
        }
    } while (changed);

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(lt_vec, ctx);
    fmpz_mpoly_clear(temp, ctx);
    fmpz_mpoly_clear(gcd, ctx);
}

void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_vec_t basis, s_pairs;
    fmpz_mpoly_vec_init(basis, generators->length, ctx);
    fmpz_mpoly_vec_init(s_pairs, 0, ctx);

    fmpz_mpoly_vec_set(basis, generators, ctx);

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

    for (slong i = 0; i < s_pairs->length; i++) {
        fmpz_mpoly_t reduced_s_pair;
        fmpz_mpoly_init(reduced_s_pair, ctx);
        fmpz_mpoly_set(reduced_s_pair, fmpz_mpoly_vec_entry(s_pairs, i), ctx);
        reduce_by_vector(reduced_s_pair, basis, ctx);

        if (!fmpz_mpoly_is_zero(reduced_s_pair, ctx)) {
            fmpz_mpoly_vec_append(basis, reduced_s_pair, ctx);

            for (slong j = 0; j < basis->length - 1; j++) {
                fmpz_mpoly_t s_pair;
                fmpz_mpoly_init(s_pair, ctx);
                construct_s_pair(s_pair, fmpz_mpoly_vec_entry(basis, j), fmpz_mpoly_vec_entry(basis, basis->length - 1), ctx);

                if (!fmpz_mpoly_is_zero(s_pair, ctx)) {
                    fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
                }

                fmpz_mpoly_clear(s_pair, ctx);
            }
        }

        fmpz_mpoly_clear(reduced_s_pair, ctx);
    }

    fmpz_mpoly_vec_set(output, basis, ctx);

    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(s_pairs, ctx);
}

void test_case_1(void) {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_vec_t generators, output;
    const char *varnames[] = {"x", "y", "z"};

    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);
    fmpz_mpoly_vec_init(generators, 2, ctx);
    fmpz_mpoly_vec_init(output, 0, ctx);

    fmpz_mpoly_set_str_pretty(fmpz_mpoly_vec_entry(generators, 0), "2*x+3*y+4*z-5", varnames, ctx);
    fmpz_mpoly_set_str_pretty(fmpz_mpoly_vec_entry(generators, 1), "3*x+4*y+5*z-2", varnames, ctx);

    buchberger_naive(output, generators, ctx);

    for (slong i = 0; i < output->length; i++) {
        char *str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(output, i), varnames, ctx);
        printf("%s\n", str);
        flint_free(str);
    }

    fmpz_mpoly_vec_clear(generators, ctx);
    fmpz_mpoly_vec_clear(output, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

int main(void) {
    test_case_1();
    return 0;
}
