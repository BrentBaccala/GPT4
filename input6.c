#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

/* Assume these function declarations */
void fmpz_mpoly_leadterm(fmpz_mpoly_t output, const fmpz_mpoly_t input, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t output, const fmpz_mpoly_t input1, const fmpz_mpoly_t input2, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduce, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx);

void buchberger_reduced(fmpz_mpoly_vec_t output, const fmpz_mpoly_vec_t input, const fmpz_mpoly_ctx_t ctx) {
    slong i, j;
    fmpz_mpoly_t poly, tmp;
    fmpz_mpoly_vec_t basis, new_basis;

    /* Initialize polynomial vector and polynomial */
    fmpz_mpoly_vec_init(basis, 0, ctx);
    fmpz_mpoly_vec_init(new_basis, 0, ctx);
    fmpz_mpoly_init(poly, ctx);
    fmpz_mpoly_init(tmp, ctx);

    /* Construct a naive Groebner basis */
    buchberger_naive(basis, input, ctx);

    /* Loop over all of the polynomials in the basis */
    for (i = 0; i < basis->length; i++) {
        /* Construct the polynomial vector of all the other polynomials in the basis */
        for (j = 0; j < basis->length; j++) {
            if (i != j) {
                fmpz_mpoly_set(tmp, fmpz_mpoly_vec_entry(basis, j), ctx);
                fmpz_mpoly_vec_append(new_basis, tmp, ctx);
            }
        }

        /* Reduce the selected polynomial by the polynomial vector composed of all the other polynomials */
        fmpz_mpoly_set(poly, fmpz_mpoly_vec_entry(basis, i), ctx);
        reduce_by_vector(poly, new_basis, 0, ctx);

        /* Replace the polynomial in the basis vector by its reduction */
        if (!fmpz_mpoly_is_zero(poly, ctx)) {
            fmpz_mpoly_vec_append(output, poly, ctx);
        }

        /* Clear and reset new_basis for next loop */
        fmpz_mpoly_vec_clear(new_basis, ctx);
        fmpz_mpoly_vec_init(new_basis, 0, ctx);
    }

    /* Clean up */
    fmpz_mpoly_clear(poly, ctx);
    fmpz_mpoly_clear(tmp, ctx);
    fmpz_mpoly_vec_clear(basis, ctx);
    fmpz_mpoly_vec_clear(new_basis, ctx);
}
