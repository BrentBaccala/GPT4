#include <stdio.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_vec_init(fmpz_mpoly_vec_t vec, slong len, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_vec_init(vec, len, ctx);
}

void fmpz_mpoly_leadterm(fmpz_mpoly_t res, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set(res, poly, ctx);
    if (fmpz_mpoly_length(poly, ctx) > 1)
        fmpz_mpoly_truncate(res, 1, ctx);
}

void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
{
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
    fmpz_mpoly_divexact(temp1, temp1, gcd, ctx);
    fmpz_mpoly_divexact(temp2, temp2, gcd, ctx);
    fmpz_mpoly_sub(s_pair, temp1, temp2, ctx);

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_t lt1, lt2, gcd, temp;

    fmpz_mpoly_init(lt1, ctx);
    fmpz_mpoly_init(lt2, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    int found = 1;
    while (found)
    {
        found = 0;
        for (i = 0; i < vec->length; i++)
        {
            fmpz_mpoly_t vec_poly = fmpz_mpoly_vec_entry(vec, i);
            fmpz_mpoly_leadterm(lt1, poly, ctx);
            fmpz_mpoly_leadterm(lt2, vec_poly, ctx);
            if (fmpz_mpoly_divides(gcd, lt1, lt2, ctx))
            {
                fmpz_mpoly_mul(temp, vec_poly, gcd, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                found = 1;
                break;
            }
        }
    }

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}

void buchberger_naive(fmpz_mpoly_vec_t basis, const fmpz_mpoly_vec_t generators, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_vec_init(basis, 0, ctx);
    fmpz_mpoly_vec_set(basis, generators, ctx);

    fmpz_mpoly_vec_t s_pairs;
    fmpz_mpoly_vec_init(s_pairs, 0, ctx);

    for (slong i = 0; i < generators->length - 1; i++)
    {
        for (slong j = i + 1; j < generators->length; j++)
        {
            fmpz_mpoly_t s_pair;
            fmpz_mpoly_init(s_pair, ctx);
            construct_s_pair(s_pair, fmpz_mpoly_vec_entry(generators, i), fmpz_mpoly_vec_entry(generators, j), ctx);
            fmpz_mpoly_vec_append(s_pairs, s_pair, ctx);
            fmpz_mpoly_clear(s_pair, ctx);
        }
    }

    for (slong i = 0; i < s_pairs->length; i++)
    {
        fmpz_mpoly_t reduced_s_pair;
        fmpz_mpoly_init(reduced_s_pair, ctx);
        fmpz_mpoly_set(reduced_s_pair, fmpz_mpoly_vec_entry(s_pairs, i), ctx);
        reduce_by_vector(reduced_s_pair, basis, ctx);

        if (!fmpz_mpoly_is_zero(reduced_s_pair, ctx))
        {
            fmpz_mpoly_vec_append(basis, reduced_s_pair, ctx);
            for (slong j = 0; j < basis->length - 1; j++)
            {
                fmpz_mpoly_t new_s_pair;
                fmpz_mpoly_init(new_s_pair, ctx);
                construct_s_pair(new_s_pair, fmpz_mpoly_vec_entry(basis, j), reduced_s_pair, ctx);
                fmpz_mpoly_vec_append(s_pairs, new_s_pair, ctx);
                fmpz_mpoly_clear(new_s_pair, ctx);
            }
        }
        fmpz_mpoly_clear(reduced_s_pair, ctx);
    }

    fmpz_mpoly_vec_clear(s_pairs, ctx);
}

void test_case_1()
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX, fmpz_QQ_ctx_mp());
    fmpz_mpoly_vec_t generators;
    fmpz_mpoly_vec_init(generators, 2, ctx);

    fmpz_mpoly_t poly1, poly2;
    fmpz_mpoly_init(poly1, ctx);
    fmpz_mpoly_init(poly2, ctx);

    fmpz_mpoly_set_str_pretty(poly1, "2*x+3*y+4*z-5", "x y z", ctx);
    fmpz_mpoly_set_str_pretty(poly2, "3*x+4*y+5*z-2", "x y z", ctx);

    fmpz_mpoly_vec_set_coeff(generators, 0, poly1, ctx);
    fmpz_mpoly_vec_set_coeff(generators, 1, poly2, ctx);

    fmpz_mpoly_vec_t basis;
    buchberger_naive(basis, generators, ctx);

    for (slong i = 0; i < basis->length; i++)
    {
        fmpz_mpoly_t basis_poly = fmpz_mpoly_vec_entry(basis, i);
        fmpz_mpoly_print_pretty(basis_poly, "x y z", ctx);
        printf("\n");
    }

    fmpz_mpoly_clear(poly1, ctx);
    fmpz_mpoly_clear(poly2, ctx);
    fmpz_mpoly_vec_clear(generators, ctx);
    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

int main()
{
    test_case_1();
    return 0;
}
