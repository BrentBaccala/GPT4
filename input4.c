#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t res, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s_pair, const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t res, const fmpz_mpoly_vec_t gens, const fmpz_mpoly_ctx_t ctx);

void test_case_1() {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);

    fmpz_mpoly_t poly1, poly2;
    fmpz_mpoly_init(poly1, ctx);
    fmpz_mpoly_init(poly2, ctx);

    const char *var_names[] = {"x", "y", "z"};
    fmpz_mpoly_set_str_pretty(poly1, "2*x+3*y+4*z-5", var_names, ctx);
    fmpz_mpoly_set_str_pretty(poly2, "3*x+4*y+5*z-2", var_names, ctx);

    fmpz_mpoly_vec_t vec;
    fmpz_mpoly_vec_init(vec, 0, ctx);
    fmpz_mpoly_vec_append(vec, poly1, ctx);
    fmpz_mpoly_vec_append(vec, poly2, ctx);

    fmpz_mpoly_vec_t res;
    fmpz_mpoly_vec_init(res, 0, ctx);
    buchberger_naive(res, vec, ctx);

    for (slong i = 0; i < res->length; i++) {
        char *poly_str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(res, i), var_names, ctx);
        printf("%s\n", poly_str);
        free(poly_str);
    }

    fmpz_mpoly_vec_clear(vec, ctx);
    fmpz_mpoly_vec_clear(res, ctx);
    fmpz_mpoly_clear(poly1, ctx);
    fmpz_mpoly_clear(poly2, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

int main() {
    test_case_1();
    return 0;
}
