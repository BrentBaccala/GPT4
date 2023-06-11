#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_vec.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_set(lt, poly, ctx);
    fmpz_mpoly_truncate(lt, 1, ctx);
}

void construct_s_pair(const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, fmpz_mpoly_t spair, const fmpz_mpoly_ctx_t ctx) {
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

    fmpz_mpoly_sub(spair, temp1, temp2, ctx);

    flint_printf("Constructed s-pair for two polynomials\n");

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}

void reduce_by_set(fmpz_mpoly_t inout_poly, const fmpz_mpoly_vec_t in_vec, const fmpz_mpoly_ctx_t ctx) {
    int changed;
    slong i;
    fmpz_mpoly_t lt_inout, lt_vec, gcd, temp;

    fmpz_mpoly_init(lt_inout, ctx);
    fmpz_mpoly_init(lt_vec, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    do {
        changed = 0;

        for (i = 0; i < in_vec->length; i++) {
            fmpz_mpoly_leadterm(lt_inout, inout_poly, ctx);
            fmpz_mpoly_leadterm(lt_vec, in_vec->coeffs + i, ctx);

            if (fmpz_mpoly_divides(gcd, lt_inout, lt_vec, ctx)) {
                fmpz_mpoly_mul(temp, in_vec->coeffs + i, gcd, ctx);
                fmpz_mpoly_sub(inout_poly, inout_poly, temp, ctx);
                changed = 1;
            }
        }
    } while (changed);

    fmpz_mpoly_clear(lt_inout, ctx);
    fmpz_mpoly_clear(lt_vec, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}
