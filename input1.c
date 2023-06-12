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

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx) {
    slong i;
    int changed;
    fmpz_mpoly_t lt_poly, tmp_poly;
    char **var_names;
    slong nvars = ctx->minfo->nvars;

    var_names = (char **) malloc(nvars * sizeof(char *));
    for (i = 0; i < nvars; i++) {
        var_names[i] = (char *) malloc(16 * sizeof(char));
        snprintf(var_names[i], 16, "x%ld", i + 1);
    }

    fprintf(stderr, "Reducing: %s\n", fmpz_mpoly_get_str_pretty(poly, var_names, ctx));

    fmpz_mpoly_init(lt_poly, ctx);
    fmpz_mpoly_init(tmp_poly, ctx);

    do {
        changed = 0;
        fmpz_mpoly_leadterm(lt_poly, poly, ctx);

        for (i = 0; i < vec->length; i++) {
            fmpz_mpoly_t vec_entry;
            fmpz_mpoly_init_set(vec_entry, fmpz_mpoly_vec_entry(vec, i), ctx);

            if (fmpz_mpoly_divides(tmp_poly, lt_poly, vec_entry, ctx)) {
                fprintf(stderr, "Divisible by: %s\n", fmpz_mpoly_get_str_pretty(vec_entry, var_names, ctx));

                fmpz_mpoly_mul(tmp_poly, tmp_poly, vec_entry, ctx);
                fmpz_mpoly_sub(poly, poly, tmp_poly, ctx);
                fmpz_mpoly_leadterm(lt_poly, poly, ctx);
                changed = 1;
                break;
            }

            fmpz_mpoly_clear(vec_entry, ctx);
        }
    } while (changed && !fmpz_mpoly_is_zero(poly, ctx));

    fprintf(stderr, "Reduced to: %s\n", fmpz_mpoly_get_str_pretty(poly, var_names, ctx));

    fmpz_mpoly_clear(lt_poly, ctx);
    fmpz_mpoly_clear(tmp_poly, ctx);

    for (i = 0; i < nvars; i++) {
        free(var_names[i]);
    }
    free(var_names);
}
