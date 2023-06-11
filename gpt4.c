#include <stdio.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_factor.h>
#include <calcium/utils_flint.h>

void fmpz_mpoly_leadterm(const fmpz_mpoly_t poly, fmpz_mpoly_t lt, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_set(lt, poly, ctx);
    fmpz_mpoly_truncate(lt, 1, ctx);
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
    fmpz_mpoly_div(temp1, temp1, gcd, ctx);

    fmpz_mpoly_mul(temp2, poly2, lt1, ctx);
    fmpz_mpoly_div(temp2, temp2, gcd, ctx);

    fmpz_mpoly_sub(s_pair, temp1, temp2, ctx);

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    int changed;
    fmpz_mpoly_t lt_poly, lt_vec, gcd, temp;

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(lt_vec, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    do {
        changed = 0;
        fmpz_mpoly_leadterm(poly, lt_poly, ctx);

        for (slong i = 0; i < vec->length; i++) {
            fmpz_mpoly_t vec_poly = *fmpz_mpoly_vec_entry(vec, i);
            fmpz_mpoly_leadterm(vec_poly, lt_vec, ctx);

            if (fmpz_mpoly_divides(gcd, lt_poly, lt_vec, ctx)) {
                fmpz_mpoly_mul(temp, vec_poly, gcd, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                changed = 1;
                break;
            }
        }
    } while (changed);

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(lt_vec, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}

void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_vec_t basis, s_pairs;
    fmpz_mpoly_vec_init(basis, ctx);
    fmpz_mpoly_vec_init(s_pairs, ctx);

    fmpz_mpoly_vec_set(basis, generators, ctx);

    for (slong i = 0; i < basis->length; i++) {
        for (slong j = i + 1; j < basis->length; j++) {
            fmpz_mpoly_t s_pair;
            fmpz_mpoly_init(s_pair, ctx);
            construct_s_pair(*fmpz_mpoly_vec_entry(basis, i), *fmpz_mpoly_vec_entry(basis, j), s_pair, ctx);

            if (!fmpz_mpoly_is_zero(s_pair, ctx)) {
                fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
            }

            fmpz_mpoly_clear(s_pair, ctx);
        }
    }

    for (slong i = 0; i < s_pairs->length; i++) {
        fmpz_mpoly_t reduced_s_pair = *fmpz_mpoly_vec_entry(s_pairs, i);
        reduce_by_vector(reduced_s_pair, basis, ctx);

        if (!fmpz_mpoly_is_zero(reduced_s_pair, ctx)) {
            fmpz_mpoly_vec_append(basis, reduced_s_pair, ctx);

            for (slong j = 0; j < basis->length - 1; j++)                {
                fmpz_mpoly_t new_s_pair;
                fmpz_mpoly_init(new_s_pair, ctx);
                construct_s_pair(*fmpz_mpoly_vec_entry(basis, j), reduced_s_pair, new_s_pair, ctx);

                if (!fmpz_mpoly_is_zero(new_s_pair, ctx)) {
                    fmpz_mpoly_vec_append(s_pairs, new_s_pair, ctx);
                }

                fmpz_mpoly_clear(new_s_pair, ctx);
            }
        }
    }

    fmpz_mpoly_vec_set(output, basis, ctx);

    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(s_pairs, ctx);
}
