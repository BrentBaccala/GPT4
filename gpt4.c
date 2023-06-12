#include <stdio.h>
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_vec.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_set(lt, poly, ctx);
    fmpz_mpoly_truncate(lt, 1, ctx);
}

void construct_s_pair(fmpz_mpoly_t spair, const fmpz_mpoly_t p1, const fmpz_mpoly_t p2, const fmpz_mpoly_ctx_t ctx) {
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
    fmpz_mpoly_mul(temp2, p2, lt1, ctx);

    fmpz_mpoly_divexact(temp1, temp1, gcd, ctx);
    fmpz_mpoly_divexact(temp2, temp2, gcd, ctx);

    fmpz_mpoly_sub(spair, temp1, temp2, ctx);

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}

void reduce_by_vector(fmpz_mpoly_t p, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    int changed;
    slong i;
    fmpz_mpoly_t lt_p, lt_vi, temp, gcd;

    fmpz_mpoly_init(lt_p, ctx);
    fmpz_mpoly_init(lt_vi, ctx);
    fmpz_mpoly_init(temp, ctx);
    fmpz_mpoly_init(gcd, ctx);

    do {
        changed = 0;
        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_t vi = fmpz_mpoly_vec_entry(vec, i);
            fmpz_mpoly_leadterm(lt_p, p, ctx);
            fmpz_mpoly_leadterm(lt_vi, vi, ctx);
            if (fmpz_mpoly_divides(gcd, lt_p, lt_vi, ctx)) {
                fmpz_mpoly_mul(temp, vi, gcd, ctx);
                fmpz_mpoly_sub(p, p, temp, ctx);
                changed = 1;
                break;
            }
        }
    } while (changed);

    fmpz_mpoly_clear(lt_p, ctx);
    fmpz_mpoly_clear(lt_vi, ctx);
    fmpz_mpoly_clear(temp, ctx);
    fmpz_mpoly_clear(gcd, ctx);
}

void buchberger_naive(fmpz_mpoly_vec_t basis, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_vec_t spairs;
    fmpz_mpoly_t spair;
    slong i, j;

    fmpz_mpoly_vec_init(spairs, 0, ctx);
    fmpz_mpoly_init(spair, ctx);

    fmpz_mpoly_vec_set(basis, generators, ctx);

    for (i = 0; i < generators->length - 1; i++) {
        for (j = i + 1; j < generators->length; j++) {
            construct_s_pair(spair, fmpz_mpoly_vec_entry(generators, i), fmpz_mpoly_vec_entry(generators, j), ctx);
            if (!fmpz_mpoly_is_zero(spair, ctx)) {
                fmpz_mpoly_vec_append(spairs, spair, ctx);
            }
        }
    }

    for (i = 0; i < spairs->length; i++) {
        fmpz_mpoly_t ri = fmpz_mpoly_vec_entry(spairs, i);
        reduce_by_vector(ri, basis, ctx);
        if (!fmpz_mpoly_is_zero(ri, ctx)) {
            fmpz_mpoly_vec_append(basis, ri, ctx);
            for (j = 0; j < basis->length - 1; j++) {
                construct_s_pair(spair, fmpz_mpoly_vec_entry(basis, j), ri, ctx);
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

    fmpz_mpoly_set_str_pretty(fmpz_mpoly_vec_entry(generators, 0), "2*x + 3*y + 4*z - 5", vars, ctx);
    fmpz_mpoly_set_str_pretty(fmpz_mpoly_vec_entry(generators, 1), "3*x + 4*y + 5*z - 2", vars, ctx);

    buchberger_naive(basis, generators, ctx);

    for (slong i = 0; i < basis->length; i++) {
        char *str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(basis, i), vars, ctx);
        printf("Basis element %ld: %s\n", i, str);
        flint_free(str);
    }

    fmpz_mpoly_vec_clear(generators, ctx);
    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

int main() {
    test_case_1();
    return 0;
}
