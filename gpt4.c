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

void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lt1, lt2, gcd, temp1, temp2;
    fmpz_mpoly_init(lt1, ctx);
    fmpz_mpoly_init(lt2, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp1, ctx);
    fmpz_mpoly_init(temp2, ctx);

    fprintf(stderr, "Constructing s-pair...\n");
    fmpz_mpoly_leadterm(lt1, poly1, ctx);
    fmpz_mpoly_leadterm(lt2, poly2, ctx);
    fmpz_mpoly_gcd(gcd, lt1, lt2, ctx);

    fmpz_mpoly_mul(temp1, poly1, lt2, ctx);
    fmpz_mpoly_mul(temp2, poly2, lt1, ctx);
    fmpz_mpoly_divexact(temp1, temp1, gcd, ctx);
    fmpz_mpoly_divexact(temp2, temp2, gcd, ctx);

    fmpz_mpoly_sub(s_pair, temp1, temp2, ctx);
    fprintf(stderr, "Constructed s-pair for input polynomials.\n");

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    slong i;
    int changed = 1;
    fmpz_mpoly_t lt_poly, lt_vec, gcd, temp;

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_vec, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    while (changed) {
        changed = 0;
        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_leadterm(lt_poly, poly, ctx);
            fmpz_mpoly_leadterm(lt_vec, fmpz_mpoly_vec_entry(vec, i), ctx);

            if (fmpz_mpoly_divides(gcd, lt_poly, lt_vec, ctx)) {
                fprintf(stderr, "Reducing...\n");
                fmpz_mpoly_mul(temp, fmpz_mpoly_vec_entry(vec, i), gcd, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                changed = 1;
            }
        }
    }

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(lt_vec, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}

void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_vec_t basis, s_pairs;
    slong i, j;
    fmpz_mpoly_t s_pair, reduced_s_pair;

    fmpz_mpoly_vec_init(basis, 0, ctx);
    fmpz_mpoly_vec_init(s_pairs, 0, ctx);
    fmpz_mpoly_init(s_pair, ctx);
    fmpz_mpoly_init(reduced_s_pair, ctx);

    for (i = 0; i < generators->length; i++) {
        fmpz_mpoly_vec_append(basis, fmpz_mpoly_vec_entry(generators, i), ctx);
    }

    for (i = 0; i < basis->length - 1; i++) {
        for (j = i + 1; j < basis->length; j++) {
            construct_s_pair(s_pair, fmpz_mpoly_vec_entry(basis, i), fmpz_mpoly_vec_entry(basis, j), ctx);
            if (!fmpz_mpoly_is_zero(s_pair, ctx)) {
                fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
            }
        }
    }

    for (i = 0; i <s_pairs->length; i++) {
        fmpz_mpoly_set(reduced_s_pair, fmpz_mpoly_vec_entry(s_pairs, i), ctx);
        reduce_by_vector(reduced_s_pair, basis, ctx);
        if (!fmpz_mpoly_is_zero(reduced_s_pair, ctx)) {
            fmpz_mpoly_vec_append(basis, reduced_s_pair, ctx);
            for (j = 0; j < basis->length - 1; j++) {
                construct_s_pair(s_pair, fmpz_mpoly_vec_entry(basis, j), reduced_s_pair, ctx);
                if (!fmpz_mpoly_is_zero(s_pair, ctx)) {
                    fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
                }
            }
        }
    }

    for (i = 0; i < basis->length; i++) {
        fmpz_mpoly_vec_append(output, fmpz_mpoly_vec_entry(basis, i), ctx);
    }

    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(s_pairs, ctx);
    fmpz_mpoly_clear(s_pair, ctx);
    fmpz_mpoly_clear(reduced_s_pair, ctx);
}

void test_case_1() {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_vec_t generators, output;
    fmpz_mpoly_t poly1, poly2;
    const char *varnames[] = {"x", "y", "z"};

    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);
    fmpz_mpoly_vec_init(generators, 0, ctx);
    fmpz_mpoly_vec_init(output, 0, ctx);

    fmpz_mpoly_init(poly1, ctx);
    fmpz_mpoly_init(poly2, ctx);

    fmpz_mpoly_set_str_pretty(poly1, "2*x+3*y+4*z-5", varnames, ctx);
    fmpz_mpoly_set_str_pretty(poly2, "3*x+4*y+5*z-2", varnames, ctx);

    fmpz_mpoly_vec_append(generators, poly1, ctx);
    fmpz_mpoly_vec_append(generators, poly2, ctx);

    buchberger_naive(output, generators, ctx);

    for (slong i = 0; i < output->length; i++) {
        char *str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(output, i), varnames, ctx);
        printf("%s\n", str);
        flint_free(str);
    }

    fmpz_mpoly_vec_clear(generators, ctx);
    fmpz_mpoly_vec_clear(output, ctx);
    fmpz_mpoly_clear(poly1, ctx);
    fmpz_mpoly_clear(poly2, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

int main() {
    test_case_1();
    return 0;
}
