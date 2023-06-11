#include "flint.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_vec.h"
#include "utils_flint.h"

void fmpz_mpoly_leadterm(const fmpz_mpoly_t poly, fmpz_mpoly_t leadterm, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set(leadterm, poly, ctx);
    fmpz_mpoly_truncate(leadterm, 1, ctx);
}

void construct_s_pair(const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, fmpz_mpoly_t spair, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t leadterm1, leadterm2, gcd, mul1, mul2;

    fmpz_mpoly_init(leadterm1, ctx);
    fmpz_mpoly_init(leadterm2, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(mul1, ctx);
    fmpz_mpoly_init(mul2, ctx);

    fmpz_mpoly_leadterm(poly1, leadterm1, ctx);
    fmpz_mpoly_leadterm(poly2, leadterm2, ctx);

    fmpz_mpoly_gcd(gcd, leadterm1, leadterm2, ctx);

    fmpz_mpoly_mul(mul1, poly1, leadterm2, ctx);
    fmpz_mpoly_mul(mul2, poly2, leadterm1, ctx);

    fmpz_mpoly_div(mul1, mul1, gcd, ctx);
    fmpz_mpoly_div(mul2, mul2, gcd, ctx);

    fmpz_mpoly_sub(spair, mul1, mul2, ctx);

    fmpz_mpoly_clear(leadterm1, ctx);
    fmpz_mpoly_clear(leadterm2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(mul1, ctx);
    fmpz_mpoly_clear(mul2, ctx);
}

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx)
{
    int changed;
    fmpz_mpoly_t leadterm_poly, leadterm_vec, gcd, mul;

    fmpz_mpoly_init(leadterm_poly, ctx);
    fmpz_mpoly_init(leadterm_vec, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(mul, ctx);

    do
    {
        changed = 0;
        for (slong i = 0; i < vec->length; i++)
        {
            fmpz_mpoly_leadterm(poly, leadterm_poly, ctx);
            fmpz_mpoly_leadterm(vec->entries[i], leadterm_vec, ctx);

            if (fmpz_mpoly_divides(gcd, leadterm_poly, leadterm_vec, ctx))
            {
                fmpz_mpoly_mul(mul, vec->entries[i], gcd, ctx);
                fmpz_mpoly_sub(poly, poly, mul, ctx);
                changed = 1;
            }
        }
    } while (changed);

    fmpz_mpoly_clear(leadterm_poly, ctx);
    fmpz_mpoly_clear(leadterm_vec, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(mul, ctx);
}

void buchberger_naive(const fmpz_mpoly_vec_t generators, fmpz_mpoly_vec_t basis, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_vec_t spairs;
    fmpz_mpoly_t spair;

    fmpz_mpoly_vec_init(basis, generators->length, ctx);
    fmpz_mpoly_vec_set(basis, generators, ctx);

    fmpz_mpoly_vec_init(spairs, 0, ctx);
    fmpz_mpoly_init(spair, ctx);

    for (slong i = 0; i < generators->length; i++)
    {
        for (slong j = i + 1; j < generators->length; j++)
        {
            construct_s_pair(generators->entries[i], generators->entries[j], spair, ctx);
            if (!fmpz_mpoly_is_zero(spair, ctx))
            {
                fmpz_mpoly_vec_append(spairs, spair, ctx);
            }
        }
    }

    for (slong i = 0; i < spairs->length; i++)
    {
        reduce_by_vector(spairs->entries[i], basis, ctx);
        if (!fmpz_mpoly_is_zero(spairs->entries[i], ctx))
        {
            fmpz_mpoly_vec_append(basis, spairs->entries[i], ctx);
            for (slong j = 0; j < basis->length - 1; j++)
            {
                construct_s_pair(basis->entries[j], basis->entries[basis->length - 1], spair, ctx);
                if (!fmpz_mpoly_is_zero(spair, ctx))
                {
                    fmpz_mpoly_vec_append(spairs, spair, ctx);
                }
            }
        }
    }

    fmpz_mpoly_clear(spair, ctx);
    fmpz_mpoly_vec_clear(spairs, ctx);
}
