#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t);
void construct_s_pair(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t);
void reduce_by_vector(fmpz_mpoly_t, const fmpz_mpoly_vec_t, int, const fmpz_mpoly_ctx_t);
void buchberger_naive(fmpz_mpoly_vec_t, const fmpz_mpoly_vec_t, const fmpz_mpoly_ctx_t);
void buchberger_reduced(fmpz_mpoly_vec_t, const fmpz_mpoly_vec_t, const fmpz_mpoly_ctx_t);

void reduce_by_vector(fmpz_mpoly_t f, const fmpz_mpoly_vec_t V, int lead_reduc, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t f_lead, g_lead, g, term;
    slong i, j, n = V->length;
    char *f_str;
    const char *vars[] = {"x", "y", "z"};
    
    f_str = fmpz_mpoly_get_str_pretty(f, vars, ctx);
    flint_fprintf(stderr, "Input: %s\n", f_str);
    flint_free(f_str);

    fmpz_mpoly_init(f_lead, ctx);
    fmpz_mpoly_init(g_lead, ctx);
    fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_init(term, ctx);

    fmpz_mpoly_leadterm(f_lead, f, ctx);

    for (i = 0; i < n; i++) {
        fmpz_mpoly_set(g, fmpz_mpoly_vec_entry(V, i), ctx);
        fmpz_mpoly_leadterm(g_lead, g, ctx);

        if (lead_reduc) {
            if (fmpz_mpoly_divides(term, f_lead, g_lead, ctx)) {
                fmpz_mpoly_mul(term, g, term, ctx);
                fmpz_mpoly_sub(f, f, term, ctx);
            }
        } else {
            for (j = 0; j < fmpz_mpoly_length(f, ctx); j++) {
                fmpz_mpoly_get_term_coeff_fmpz(term, f, j, ctx);
                if (fmpz_mpoly_divides(term, term, g_lead, ctx)) {
                    fmpz_mpoly_mul(term, g, term, ctx);
                    fmpz_mpoly_sub(f, f, term, ctx);
                }
            }
        }

        if (fmpz_mpoly_is_zero(f, ctx))
            break;

        fmpz_mpoly_leadterm(f_lead, f, ctx);
    }

    f_str = fmpz_mpoly_get_str_pretty(f, vars, ctx);
    flint_fprintf(stderr, "Output: %s\n", f_str);
    flint_free(f_str);

    fmpz_mpoly_clear(f_lead, ctx);
    fmpz_mpoly_clear(g_lead, ctx);
    fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_clear(term, ctx);
}
