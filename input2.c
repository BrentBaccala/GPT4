#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);

void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t lt1, lt2, gcd, temp1, temp2;
    const char *vars[] = {"x", "y", "z"};

    fmpz_mpoly_init(lt1, ctx);
    fmpz_mpoly_init(lt2, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp1, ctx);
    fmpz_mpoly_init(temp2, ctx);

    fprintf(stderr, "Constructing s-pair...\n");
    char *str_poly1 = fmpz_mpoly_get_str_pretty(poly1, vars, ctx);
    char *str_poly2 = fmpz_mpoly_get_str_pretty(poly2, vars, ctx);
    fprintf(stderr, "Input polynomials: \n%s\n%s\n", str_poly1, str_poly2);
    flint_free(str_poly1);
    flint_free(str_poly2);

    fmpz_mpoly_leadterm(lt1, poly1, ctx);
    fmpz_mpoly_leadterm(lt2, poly2, ctx);

    fmpz_mpoly_gcd(gcd, lt1, lt2, ctx);

    fmpz_mpoly_mul(temp1, poly1, lt2, ctx);
    fmpz_mpoly_mul(temp2, poly2, lt1, ctx);

    fmpz_mpoly_div(temp1, temp1, gcd, ctx);
    fmpz_mpoly_div(temp2, temp2, gcd, ctx);

    fmpz_mpoly_sub(s_pair, temp1, temp2, ctx);

    char *str_s_pair = fmpz_mpoly_get_str_pretty(s_pair, vars, ctx);
    fprintf(stderr, "Constructed s-pair for input polynomials: \n%s\n", str_s_pair);
    flint_free(str_s_pair);

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}
