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

void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx) {
    fprintf(stderr, "Beginning construction of s-pair.\n");
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
    fmpz_mpoly_div(temp1, temp1, gcd, ctx);
    fmpz_mpoly_div(temp2, temp2, gcd, ctx);
    fmpz_mpoly_sub(s_pair, temp1, temp2, ctx);

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);

    fprintf(stderr, "Constructed s-pair for input polynomials.\n");
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    int changed = 1;
    fmpz_mpoly_t lt1, lt2, gcd, temp;
    fmpz_mpoly_init(lt1, ctx);
    fmpz_mpoly_init(lt2, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    while (changed && !fmpz_mpoly_is_zero(poly, ctx)) {
        changed = 0;
        for (slong i = 0; i < vec->length; ++i) {
            fmpz_mpoly_t curr_poly = *fmpz_mpoly_vec_entry(vec, i);
            fmpz_mpoly_leadterm(lt1, poly, ctx);
            fmpz_mpoly_leadterm(lt2, curr_poly, ctx);

            if (fmpz_mpoly_divides(gcd, lt1, lt2, ctx)) {
                fprintf(stderr, "Reducing by vector entry %ld.\n", i);
                fmpz_mpoly_mul(temp, curr_poly, gcd, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                changed = 1;
                break;
            }
        }
    }

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}

void buchberger_naive(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gens, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_vec_t basis, s_pairs;
    fmpz_mpoly_vec_init(basis, 0, ctx);
    fmpz_mpoly_vec_init(s_pairs, 0, ctx);

    for (slong i = 0; i < gens->length; ++i) {
        fmpz_mpoly_vec_append(basis, fmpz_mpoly_vec_entry(gens, i), ctx);
    }

    for (slong i = 0; i < gens->length - 1; ++i) {
        for (slong j = i + 1; j < gens->length; ++j) {
            fmpz_mpoly_t s_pair;
            fmpz_mpoly_init(s_pair, ctx);
            construct_s_pair(s_pair, fmpz_mpoly_vec_entry(gens, i), fmpz_mpoly_vec_entry(gens, j), ctx);
            fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
            fmpz_mpoly_clear(s_pair, ctx);
        }
    }

    for (slong i = 0;i < s_pairs->length; ++i) {
        fmpz_mpoly_t reduced_s_pair;
        fmpz_mpoly_init(reduced_s_pair, ctx);
        fmpz_mpoly_set(reduced_s_pair, fmpz_mpoly_vec_entry(s_pairs, i), ctx);
        reduce_by_vector(reduced_s_pair, basis, ctx);

        if (!fmpz_mpoly_is_zero(reduced_s_pair, ctx)) {
            fmpz_mpoly_vec_append(basis, reduced_s_pair, ctx);
            for (slong j = 0; j < basis->length - 1; ++j) {
                fmpz_mpoly_t new_s_pair;
                fmpz_mpoly_init(new_s_pair, ctx);
                construct_s_pair(new_s_pair, fmpz_mpoly_vec_entry(basis, j), reduced_s_pair, ctx);
                fmpz_mpoly_vec_append(s_pairs, new_s_pair, ctx);
                fmpz_mpoly_clear(new_s_pair, ctx);
            }
        }
        fmpz_mpoly_clear(reduced_s_pair, ctx);
    }

    for (slong i = 0; i < basis->length; ++i) {
        fmpz_mpoly_vec_append(res, fmpz_mpoly_vec_entry(basis, i), ctx);
    }

    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(s_pairs, ctx);
}

void test_case_1() {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);
    const char* vars[] = {"x", "y", "z"};

    fmpz_mpoly_vec_t gens, basis;
    fmpz_mpoly_vec_init(gens, 0, ctx);
    fmpz_mpoly_vec_init(basis, 0, ctx);

    fmpz_mpoly_t poly1, poly2;
    fmpz_mpoly_init(poly1, ctx);
    fmpz_mpoly_init(poly2, ctx);
    fmpz_mpoly_set_str_pretty(poly1, "2*x+3*y+4*z-5", vars, ctx);
    fmpz_mpoly_set_str_pretty(poly2, "3*x+4*y+5*z-2", vars, ctx);

    fmpz_mpoly_vec_append(gens, poly1, ctx);
    fmpz_mpoly_vec_append(gens, poly2, ctx);

    buchberger_naive(basis, gens, ctx);

    for (slong i = 0; i < basis->length; ++i) {
        char* str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(basis, i), vars, ctx);
        printf("%s\n", str);
        flint_free(str);
    }

    fmpz_mpoly_clear(poly1, ctx);
    fmpz_mpoly_clear(poly2, ctx);
    fmpz_mpoly_vec_clear(gens, ctx);
    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

int main() {
    test_case_1();
    return 0;
}
