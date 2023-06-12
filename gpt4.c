#include <stdio.h>
#include "flint/fmpz_mpoly.h"
#include "utils_flint.h"

void fmpz_mpoly_leadterm(const fmpz_mpoly_t poly, fmpz_mpoly_t leadterm, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_set(leadterm, poly, ctx);
    fmpz_mpoly_truncate(leadterm, 1, ctx);
}

void construct_s_pair(const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, fmpz_mpoly_t s_pair, const fmpz_mpoly_ctx_t ctx) {
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
    fmpz_mpoly_mul(temp2, poly2, lt1, ctx);
    fmpz_mpoly_div(temp1, temp1, gcd, ctx);
    fmpz_mpoly_div(temp2, temp2, gcd, ctx);
    fmpz_mpoly_sub(s_pair, temp1, temp2, ctx);

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    int reduced;
    slong i, len;
    fmpz_mpoly_t temp, lt1, lt2, gcd;

    fmpz_mpoly_init(temp, ctx);
    fmpz_mpoly_init(lt1, ctx);
    fmpz_mpoly_init(lt2, ctx);
    fmpz_mpoly_init(gcd, ctx);

    do {
        reduced = 0;
        len = fmpz_mpoly_vec_length(vec, ctx);

        for (i = 0; i < len; i++) {
            fmpz_mpoly_t divisor = fmpz_mpoly_vec_entry(vec, i, ctx);

            fmpz_mpoly_leadterm(poly, lt1, ctx);
            fmpz_mpoly_leadterm(divisor, lt2, ctx);
            fmpz_mpoly_gcd(gcd, lt1, lt2, ctx);

            if (fmpz_mpoly_equal(gcd, lt2, ctx)) {
                fmpz_mpoly_mul(temp, poly, lt2, ctx);
                fmpz_mpoly_div(temp, temp, gcd, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                reduced = 1;
                break;
            }
        }
    } while (reduced);

    fmpz_mpoly_clear(temp, ctx);
    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
}

void buchberger_naive(const fmpz_mpoly_vec_t generators, fmpz_mpoly_vec_t basis, const fmpz_mpoly_ctx_t ctx) {
    slong i, j, len;
    fmpz_mpoly_t s_pair;

    fmpz_mpoly_init(s_pair, ctx);
    fmpz_mpoly_vec_set(basis, generators, ctx);
    len = fmpz_mpoly_vec_length(generators, ctx);

    for (i = 0; i < len; i++) {
        for (j = i + 1; j < len; j++) {
            fmpz_mpoly_t poly1 = fmpz_mpoly_vec_entry(generators, i, ctx);
            fmpz_mpoly_t poly2 = fmpz_mpoly_vec_entry(generators, j, ctx);
            construct_s_pair(poly1, poly2, s_pair, ctx);
            reduce_by_vector(s_pair, basis, ctx);

            if (!fmpz_mpoly_is_zero(s_pair, ctx)) {
                fmpz_mpoly_vec_append(basis, s_pair, ctx);
                len = fmpz_mpoly_vec_length(generators, ctx);
            }
        }
    }

    fmpz_mpoly_clear(s_pair, ctx);
}

void test_case_1() {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_vec_t generators, basis;
    fmpz_mpoly_t poly1, poly2;
    const char *var_names[] = {"x", "y", "z"};

    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);
    fmpz_mpoly_vec_init(generators, 2, ctx);
    fmpz_mpoly_vec_init(basis, 0, ctx);
    fmpz_mpoly_init(poly1, ctx);
    fmpz_mpoly_init(poly2, ctx);

    fmpz_mpoly_set_str_pretty(poly1, "2*x+3*y+4*z-5", var_names, ctx);
    fmpz_mpoly_set_str_pretty(poly2, "3*x+4*y+5*z-2", var_names, ctx);
    fmpz_mpoly_vec_set_coeff(generators, 0, poly1, ctx);
    fmpz_mpoly_vec_set_coeff(generators, 1, poly2, ctx);

    buchberger_naive(generators, basis, ctx);

    fmpz_mpoly_vec_print_pretty(basis, var_names, ctx);

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
