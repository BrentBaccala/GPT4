#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(const fmpz_mpoly_t poly, fmpz_mpoly_t leadterm,
                         const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_set(leadterm, poly, ctx);
    fmpz_mpoly_truncate(leadterm, 1, ctx);
}

void construct_s_pair(const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                      fmpz_mpoly_t s_pair, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lt1, lt2, gcd, temp1, temp2;

    fmpz_mpoly_init(lt1, ctx);
    fmpz_mpoly_init(lt2, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp1, ctx);
    fmpz_mpoly_init(temp2, ctx);

    fmpz_mpoly_leadterm(poly1, lt1, ctx);
    fmpz_mpoly_leadterm(poly2, lt2, ctx);
    fmpz_mpoly_gcd(gcd, lt1, lt2, ctx);

    fmpz_mpoly_mul(temp1, poly1, lt2, ctx);
    fmpz_mpoly_divexact(temp1, temp1, gcd, ctx);

    fmpz_mpoly_mul(temp2, poly2, lt1, ctx);
    fmpz_mpoly_divexact(temp2, temp2, gcd, ctx);

    fmpz_mpoly_sub(s_pair, temp1, temp2, ctx);

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);

    fprintf(stderr, "Constructed s-pair for %s and %s\n", poly1, poly2);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec,
                      const fmpz_mpoly_ctx_t ctx) {
    int changed;
    slong i;
    fmpz_mpoly_t lt_poly, lt_vec, gcd, temp;

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_vec, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    do {
        changed = 0;
        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_leadterm(poly, lt_poly, ctx);
            fmpz_mpoly_leadterm(fmpz_mpoly_vec_entry(vec, i), lt_vec, ctx);

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

fmpz_mpoly_vec_t buchberger_naive(const fmpz_mpoly_vec_t generators,
                                  const fmpz_mpoly_ctx_t ctx) {
    slong i, j;
    fmpz_mpoly_vec_t basis, s_pairs;
    fmpz_mpoly_t s_pair;

    fmpz_mpoly_vec_init(basis, generators->length, ctx);
    fmpz_mpoly_vec_set(basis, generators, ctx);

    fmpz_mpoly_vec_init(s_pairs, 0, ctx);
    fmpz_mpoly_init(s_pair, ctx);

    for (i = 0; i < generators->length; i++) {
        for (j = i + 1; j < generators->length; j++) {
            construct_s_pair(fmpz_mpoly_vec_entry(generators, i),
                             fmpz_mpoly_vec_entry(generators, j), s_pair, ctx);
            fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
        }
    }

    for (i = 0; i < s_pairs->length; i++) {
        reduce_by_vector(fmpz_mpoly_vec_entry(s_pairs, i), basis, ctx);

        if (!fmpz_mpoly_is_zero(fmpz_mpoly_vec_entry(s_pairs, i), ctx)) {
            fmpz_mpoly_vec_append(basis, fmpz_mpoly_vec_entry(s_pairs, i), ctx);

            for (j = 0; j < basis->length - 1; j++) {
                construct_s_pair(fmpz_mpoly_vec_entry(basis, j),
                                 fmpz_mpoly_vec_entry(basis, basis->length - 1),
                                 s_pair, ctx);
                fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
            }
        }
    }

    fmpz_mpoly_clear(s_pair, ctx);
    fmpz_mpoly_vec_clear(s_pairs, ctx);

    return basis;
}

void test_case_1() {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_vec_t generators, basis;
    fmpz_mpoly_t poly1, poly2;
    const char *vars[] = {"x", "y", "z"};

    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);

    fmpz_mpoly_init(poly1, ctx);
    fmpz_mpoly_init(poly2, ctx);
    fmpz_mpoly_set_str_pretty(poly1, "2*x + 3*y + 4*z - 5", vars, ctx);
    fmpz_mpoly_set_str_pretty(poly2, "3*x + 4*y + 5*z - 2", vars, ctx);

    fmpz_mpoly_vec_init(generators, 2, ctx);
    fmpz_mpoly_vec_set_entry(generators, 0, poly1, ctx);
    fmpz_mpoly_vec_set_entry(generators, 1, poly2, ctx);

    basis = buchberger_naive(generators, ctx);

    for (slong i = 0; i < basis->length; i++) {
        fmpz_mpoly_print_pretty(fmpz_mpoly_vec_entry(basis, i), vars, ctx);
        printf("\n");
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
