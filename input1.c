#include <stdio.h>
#include <stdlib.h>
#include "flint/flint.h"
#include "flint/fmpz_mpoly.h"
#include "flint/fmpz_mpoly_factor.h"
#include "calcium/utils_flint.h"

void fmpz_mpoly_leadterm(fmpz_mpoly_t lt, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);
void construct_s_pair(fmpz_mpoly_t s, const fmpz_mpoly_t a, const fmpz_mpoly_t b, const fmpz_mpoly_ctx_t ctx);
void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx);
void buchberger_naive(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);
void buchberger_reduced(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx);

void reduce_by_vector(fmpz_mpoly_t poly, const fmpz_mpoly_vec_t vec, int lead_reduction, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_t lt1, lt2, temp;
    const char *vars[] = {"x", "y", "z"};

    flint_printf("Input/output polynomial: ");
    char *poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    flint_fprintf(stderr, "%s\n", poly_str);
    free(poly_str);

    fmpz_mpoly_init(lt1, ctx);
    fmpz_mpoly_init(lt2, ctx);
    fmpz_mpoly_init(temp, ctx);

    fmpz_mpoly_leadterm(lt1, poly, ctx);

    for (slong i = 0; i < vec->length; i++) {
        fmpz_mpoly_set(temp, fmpz_mpoly_vec_entry(vec, i), ctx);
        fmpz_mpoly_leadterm(lt2, temp, ctx);

        if (lead_reduction) {
            if (fmpz_mpoly_divides(temp, lt1, lt2, ctx)) {
                fmpz_mpoly_gcd(lt1, lt1, lt2, ctx);
                fmpz_mpoly_mul(temp, temp, lt1, ctx);
                fmpz_mpoly_sub(poly, poly, temp, ctx);
                if (fmpz_mpoly_is_zero(poly, ctx)) {
                    break;
                }
                fmpz_mpoly_leadterm(lt1, poly, ctx);
            }
        } else {
            for (slong j = fmpz_mpoly_length(poly, ctx) - 1; j >= 0; j--) {
                fmpz_mpoly_t term;
                fmpz_mpoly_init(term, ctx);
                fmpz_mpoly_truncate(poly, j + 1, ctx);
                fmpz_mpoly_sub(term, poly, term, ctx);
                fmpz_mpoly_leadterm(term, term, ctx);
                if (fmpz_mpoly_divides(lt2, term, lt2, ctx)) {
                    fmpz_mpoly_gcd(lt2, term, lt2, ctx);
                    fmpz_mpoly_mul(temp, temp, lt2, ctx);
                    fmpz_mpoly_sub(poly, poly, temp, ctx);
                }
                fmpz_mpoly_clear(term, ctx);
            }
        }
    }

    fmpz_mpoly_clear(lt1, ctx);
    fmpz_mpoly_clear(lt2, ctx);
    fmpz_mpoly_clear(temp, ctx);

    flint_printf("Reduced polynomial: ");
    poly_str = fmpz_mpoly_get_str_pretty(poly, vars, ctx);
    flint_fprintf(stderr, "%s\n", poly_str);
    free(poly_str);
}
