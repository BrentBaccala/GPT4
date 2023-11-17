#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t);
void construct_s_pair(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t);

void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t leadterm1, leadterm2, gcd, temp1, temp2;
    const char * vars[] = {"x", "y", "z"};

    fmpz_mpoly_init(leadterm1, ctx);
    fmpz_mpoly_init(leadterm2, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp1, ctx);
    fmpz_mpoly_init(temp2, ctx);

    fmpz_mpoly_leadterm(leadterm1, poly1, ctx);
    fmpz_mpoly_leadterm(leadterm2, poly2, ctx);

    fmpz_mpoly_gcd(gcd, leadterm1, leadterm2, ctx);

    fmpz_mpoly_mul(temp1, poly1, leadterm2, ctx);
    fmpz_mpoly_mul(temp2, poly2, leadterm1, ctx);

    fmpz_mpoly_div(temp1, temp1, gcd, ctx);
    fmpz_mpoly_div(temp2, temp2, gcd, ctx);

    fmpz_mpoly_sub(s_pair, temp1, temp2, ctx);

    fprintf(stderr, "Beginning construction of the s-pair for:\n%s\n%s\n", fmpz_mpoly_get_str_pretty(poly1, vars, ctx), fmpz_mpoly_get_str_pretty(poly2, vars, ctx));
    fprintf(stderr, "Constructed s-pair:\n%s\n", fmpz_mpoly_get_str_pretty(s_pair, vars, ctx));

    fmpz_mpoly_clear(leadterm1, ctx);
    fmpz_mpoly_clear(leadterm2, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp1, ctx);
    fmpz_mpoly_clear(temp2, ctx);
}
