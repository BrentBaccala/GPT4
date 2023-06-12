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

void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t p1, const fmpz_mpoly_t p2, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lt_p1, lt_p2, gcd, temp1, temp2;
    fmpz_mpoly_init(lt_p1, ctx);
    fmpz_mpoly_init(lt_p2, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp1, ctx);
    fmpz_mpoly_init(temp2, ctx);

    fprintf(stderr, "Constructing S-pair\n");

    fmpz_mpoly_leadterm(lt_p1, p1, ctx);
    fmpz_mpoly_leadterm(lt_p2, p2, ctx);
    fmpz_mpoly_gcd(gcd, lt_p1, lt_p2, ctx);

    fmpz_mpoly_mul(temp1, p1, lt_p2, ctx);
    fmpz_mpoly_mul(temp2, p2, lt_p1, ctx);

    fmpz_mpoly_divexact(temp1, temp1, gcd, ctx);
    fmpz_mpoly_divexact(temp2, temp2, gcd, ctx);

    fmpz_mpoly_sub(s_pair, temp1, temp2, ctx);

    fprintf(stderr, "Constructed S-pair for input polynomials\n");

    fmpz_mpoly_clear(lt_p1, ctx);
    fmpz_mpoly_clear(lt_p2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    int changed;
    fmpz_mpoly_t lt_poly, lt_vec_i, gcd, temp;

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_vec_i, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    do {
        changed = 0;
        for (slong i = 0; i < fmpz_mpoly_vec_length(vec, ctx); i++) {
            fmpz_mpoly_leadterm(lt_poly, poly, ctx);
            fmpz_mpoly_leadterm(lt_vec_i, fmpz_mpoly_vec_entry(vec, i), ctx);

            if (fmpz_mpoly_divides(gcd, lt_poly, lt_vec_i, ctx)) {
                fprintf(stderr, "Reducing by vector\n");
                fmpz_mpoly_mul(temp, fmpz_mpoly_vec_entry(vec, i), gcd, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                changed = 1;
                break;
            }
        }
    } while (changed && !fmpz_mpoly_is_zero(poly, ctx));

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(lt_vec_i, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}

void buchberger_naive(fmpz_mpoly_vec_t basis, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_vec_t s_pairs;
    fmpz_mpoly_t s_pair, reduced_s_pair;

    fmpz_mpoly_vec_init(basis, 0, ctx);
    fmpz_mpoly_vec_init(s_pairs, 0, ctx);
    fmpz_mpoly_init(s_pair, ctx);
    fmpz_mpoly_init(reduced_s_pair, ctx);

    fmpz_mpoly_vec_set(basis, generators, ctx);

    for (slong i = 0; i < fmpz_mpoly_vec_length(generators, ctx); i++) {
        for (slong j = i + 1; j < fmpz_mpoly_vec_length(generators, ctx); j++) {
            construct_s_pair(s_pair, fmpz_mpoly_vec_entry(generators, i), fmpz_mpoly_vec_entry(generators, j), ctx);
            if (!fmpz_mpoly_is_zero(s_pair, ctx)) {
                fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
            }
        }
    }

    for (slong i = 0; i < fmpz_mpoly_vec_length(s_pairs, ctx); i++) {
        fmpz_mpoly_set(reduced_s_pair, fmpz_mpoly_vec_entry(s_pairs, i), ctx);
        reduce_by_vector(reduced_s_pair, basis, ctx);

        if (!fmpz_mpoly_is_zero(reduced_s_pair, ctx)) {
            fmpz_mpoly_vec_append(basis, reduced_s_pair, ctx);

            for (slong j = 0; j < fmpz_mpoly_vec_length(basis, ctx) - 1; j++) {
                construct_s_pair(s_pair, reduced_s_pair, fmpz_mpoly_vec_entry(basis, j), ctx);
                if (!fmpz_mpoly_is_zero(s_pair, ctx)) {
                    fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
                }
            }
        }
    }

    fmpz_mpoly_vec_clear(s_pairs, ctx);
    fmpz_mpoly_clear(s_pair, ctx);
    fmpz_mpoly_clear(reduced_s_pair, ctx);
}

void test_case_1() {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t p1, p2;
    fmpz_mpoly_vec_t generators, basis;
    const char* vars[] = {"x", "y", "z"};

    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);
    fmpz_mpoly_init(p1, ctx);
    fmpz_mpoly_init(p2, ctx);
    fmpz_mpoly_vec_init(generators, 0, ctx);
    fmpz_mpoly_vec_init(basis, 0, ctx);

    fmpz_mpoly_set_str_pretty(p1, "2*x+3*y+4*z-5", vars, ctx);
    fmpz_mpoly_set_str_pretty(p2, "3*x+4*y+5*z-2", vars, ctx);

    fmpz_mpoly_vec_append(generators, p1, ctx);
    fmpz_mpoly_vec_append(generators, p2, ctx);

    buchberger_naive(basis, generators, ctx);

    for (slong i = 0; i < fmpz_mpoly_vec_length(basis, ctx); i++) {
        char* str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(basis, i), vars, ctx);
        printf("%s\n", str);
        flint_free(str);
    }

    fmpz_mpoly_vec_clear(generators, ctx);
    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_clear(p1, ctx);
    fmpz_mpoly_clear(p2, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

int main() {
    test_case_1();
    return 0;
}
