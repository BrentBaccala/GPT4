#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t spair, const fmpz_mpoly_t f, const fmpz_mpoly_t g, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t r, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);

static void reduce_by_matching_term(fmpz_mpoly_t r, const fmpz_mpoly_t matching_term, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t gcd, tmp;
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(tmp, ctx);

    fmpz_mpoly_gcd(gcd, matching_term, fmpz_mpoly_leadterm(tmp, poly, ctx), ctx);
    fmpz_mpoly_mul(tmp, poly, gcd, ctx);
    fmpz_mpoly_sub(r, r, tmp, ctx);

    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(tmp, ctx);
}

void reduce_by_vector(fmpz_mpoly_t r, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lt, lt_vec_i;
    fmpz_mpoly_init(lt, ctx);
    fmpz_mpoly_init(lt_vec_i, ctx);

    fmpz_mpoly_leadterm(lt, r, ctx);

    for (slong i = 0; i < vec->length; i++) {
        fmpz_mpoly_t poly_i;
        fmpz_mpoly_init(poly_i, ctx);
        fmpz_mpoly_set(poly_i, fmpz_mpoly_vec_entry(vec, i), ctx);
        fmpz_mpoly_leadterm(lt_vec_i, poly_i, ctx);

        if (lead_reduction) {
            if (fmpz_mpoly_divides(lt, lt, lt_vec_i, ctx)) {
                reduce_by_matching_term(r, lt, poly_i, ctx);
            }
        } else {
            for (slong j = 0; j < r->length; j++) {
                fmpz_mpoly_t term_j;
                fmpz_mpoly_init(term_j, ctx);
                fmpz_mpoly_truncate(term_j, r, j, ctx);
                fmpz_mpoly_sub(term_j, r, term_j, ctx);
                fmpz_mpoly_leadterm(term_j, term_j, ctx);

                if (fmpz_mpoly_divides(lt_vec_i, term_j, lt_vec_i, ctx)) {
                    reduce_by_matching_term(r, term_j, poly_i, ctx);
                }
                fmpz_mpoly_clear(term_j, ctx);
            }
        }

        if (fmpz_mpoly_is_zero(r, ctx)) {
            break;
        }
        fmpz_mpoly_leadterm(lt, r, ctx);

        fmpz_mpoly_clear(poly_i, ctx);
    }

    fmpz_mpoly_clear(lt, ctx);
    fmpz_mpoly_clear(lt_vec_i, ctx);
}
