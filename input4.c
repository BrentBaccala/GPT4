#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t out, const fmpz_mpoly_t in, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t out, const fmpz_mpoly_t in1, const fmpz_mpoly_t in2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t out, const fmpz_mpoly_vec_t vec, int flag, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t out, const fmpz_mpoly_vec_t in, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t out, const fmpz_mpoly_vec_t in, const fmpz_mpoly_ctx_t ctx);

int main() {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX);

    const char* vars[] = {"x", "y", "z"};

    fmpz_mpoly_vec_t test_case_1;
    fmpz_mpoly_vec_init(test_case_1, 0, ctx);
    fmpz_mpoly_t poly1, poly2;
    fmpz_mpoly_init(poly1, ctx);
    fmpz_mpoly_init(poly2, ctx);
    fmpz_mpoly_set_str_pretty(poly1, "2*x+3*y+4*z-5", vars, ctx);
    fmpz_mpoly_set_str_pretty(poly2, "3*x+4*y+5*z-2", vars, ctx);
    fmpz_mpoly_vec_append(test_case_1, poly1, ctx);
    fmpz_mpoly_vec_append(test_case_1, poly2, ctx);

    fmpz_mpoly_vec_t basis;
    fmpz_mpoly_vec_init(basis, 0, ctx);
    buchberger_naive(basis, test_case_1, ctx);

    for (int i = 0; i < basis->length; i++) {
        char* str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(basis, i), vars, ctx);
        fprintf(stderr, "%s\n", str);
        flint_free(str);
    }

    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(test_case_1, ctx);
    fmpz_mpoly_clear(poly1, ctx);
    fmpz_mpoly_clear(poly2, ctx);

    fmpz_mpoly_vec_init(test_case_1, 0, ctx);
    fmpz_mpoly_init(poly1, ctx);
    fmpz_mpoly_init(poly2, ctx);
    fmpz_mpoly_set_str_pretty(poly1, "2*x+3*y+4*z-5", vars, ctx);
    fmpz_mpoly_set_str_pretty(poly2, "3*x+4*y+5*z-2", vars, ctx);
    fmpz_mpoly_vec_append(test_case_1, poly1, ctx);
    fmpz_mpoly_vec_append(test_case_1, poly2, ctx);

    fmpz_mpoly_vec_init(basis, 0, ctx);
    buchberger_reduced(basis, test_case_1, ctx);

    for (int i = 0; i < basis->length; i++) {
        char* str = fmpz_mpoly_get_str_pretty(fmpz_mpoly_vec_entry(basis, i), vars, ctx);
        fprintf(stderr, "%s\n", str);
        flint_free(str);
    }

    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(test_case_1, ctx);
    fmpz_mpoly_clear(poly1, ctx);
    fmpz_mpoly_clear(poly2, ctx);

    fmpz_mpoly_ctx_clear(ctx);

    return 0;
}
