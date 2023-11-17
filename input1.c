#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduc, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t out, const fmpz_mpoly_vec_t in, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t out, const fmpz_mpoly_vec_t in, const fmpz_mpoly_ctx_t ctx);

void reduce_by_matching_term(fmpz_mpoly_t poly, const fmpz_mpoly_t term, const fmpz_mpoly_t reducer, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lt, gcd, temp;
    char *poly_str, *term_str, *reducer_str;
    const char *vars[] = {"x", "y", "z"};

    fmpz_mpoly_init(lt, ctx);
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(temp, ctx);

    fmpz_mpoly_leadterm(lt, reducer, ctx);
    fmpz_mpoly_gcd(gcd, term, lt, ctx);

    poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    term_str = fmpz_mpoly_get_str_pretty(term, vars, ctx);
    reducer_str = fmpz_mpoly_get_str_pretty(reducer, vars, ctx);

    fprintf(stderr, "Reducing polynomial %s by term %s with reducer %s\n", poly_str, term_str, reducer_str);

    fmpz_mpoly_mul(temp, reducer, term, ctx);
    fmpz_mpoly_divexact(temp, temp, gcd, ctx);
    fmpz_mpoly_sub(poly, poly, temp, ctx);

    flint_free(poly_str);
    flint_free(term_str);
    flint_free(reducer_str);

    fmpz_mpoly_clear(lt, ctx);
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(temp, ctx);
}
