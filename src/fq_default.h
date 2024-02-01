/*
    Copyright (C) 2021 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_DEFAULT_H
#define FQ_DEFAULT_H

#ifdef FQ_DEFAULT_INLINES_C
#define FQ_DEFAULT_INLINE
#else
#define FQ_DEFAULT_INLINE static inline
#endif

#include "nmod.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fq.h"
#include "fq_nmod.h"
#include "fq_zech.h"
#include "gr.h"

#define FQ_DEFAULT_FQ_ZECH  1
#define FQ_DEFAULT_FQ_NMOD  2
#define FQ_DEFAULT_FQ       3
#define FQ_DEFAULT_NMOD     4
#define FQ_DEFAULT_FMPZ_MOD 5

/* Data types and context ****************************************************/

#ifdef __cplusplus
extern "C" {
#endif

typedef union fq_default_struct
{
    fq_t fq;
    fq_nmod_t fq_nmod;
    fq_zech_t fq_zech;
    ulong nmod;
    fmpz_t fmpz_mod;
} fq_default_struct;

typedef fq_default_struct fq_default_t[1];

typedef struct
{
    gr_ctx_t gr_ctx;
    int _gr_inited;
    int type;
    union ctx
    {
        fq_ctx_t fq_UNUSED;
        fq_nmod_ctx_t fq_nmod_UNUSED;
        fq_zech_ctx_t fq_zech_UNUSED;
        struct {
            nmod_t mod;
            mp_limb_t a;    /* minpoly is x - a */
        } nmod;
        struct {
            fmpz_mod_ctx_t mod;
            fmpz_t a;       /* minpoly is x - a */
        } fmpz_mod;
    } ctx;
} fq_default_ctx_struct;

typedef fq_default_ctx_struct fq_default_ctx_t[1];

#define FQ_DEFAULT_TYPE(ctx) ((ctx)->type)
#define FQ_DEFAULT_GR_CTX(ctx) ((gr_ctx_struct *) ((ctx)->gr_ctx))

#define FQ_DEFAULT_CTX_FQ(ctx) ((fq_ctx_struct *) GR_CTX_DATA_AS_PTR(FQ_DEFAULT_GR_CTX(ctx)))
#define FQ_DEFAULT_CTX_FQ_ZECH(ctx) ((fq_zech_ctx_struct *) GR_CTX_DATA_AS_PTR(FQ_DEFAULT_GR_CTX(ctx)))
#define FQ_DEFAULT_CTX_FQ_NMOD(ctx) ((fq_nmod_ctx_struct *) GR_CTX_DATA_AS_PTR(FQ_DEFAULT_GR_CTX(ctx)))

#define FQ_DEFAULT_CTX_NMOD_OLD(ctx) ((ctx)->ctx.nmod)
#define FQ_DEFAULT_CTX_NMOD(ctx) ((nmod_t *) FQ_DEFAULT_GR_CTX(ctx))[0]

#define FQ_DEFAULT_CTX_FMPZ_MOD_OLD(ctx) ((ctx)->ctx.fmpz_mod)
#define FQ_DEFAULT_CTX_FMPZ_MOD(ctx) ((fmpz_mod_ctx_struct *) GR_CTX_DATA_AS_PTR(FQ_DEFAULT_GR_CTX(ctx)))

typedef struct
{
    fmpz_mod_ctx_struct * ctx;
    truth_t is_prime;
    fmpz a;    /* when used as finite field with defining polynomial x - a */
}
_gr_fmpz_mod_ctx_struct;

typedef struct
{
    nmod_t nmod;
    ulong a;   /* when used as finite field with defining polynomial x - a */
}
_gr_nmod_ctx_struct;

#define NMOD_CTX_A(ring_ctx) (&((((_gr_nmod_ctx_struct *)(ring_ctx))->a)))
#define FMPZ_MOD_CTX_A(ring_ctx) (&((((_gr_fmpz_mod_ctx_struct *)(ring_ctx))->a)))

#define FQ_DEFAULT_CTX_NMOD_A(ctx) NMOD_CTX_A(ctx)
#define FQ_DEFAULT_CTX_FMPZ_MOD_A(ctx) FMPZ_MOD_CTX_A(ctx)


void fq_default_ctx_init_type(fq_default_ctx_t ctx,
                            const fmpz_t p, slong d, const char *var, int type);

FQ_DEFAULT_INLINE void fq_default_ctx_init(fq_default_ctx_t ctx,
                                      const fmpz_t p, slong d, const char *var)
{
    fq_default_ctx_init_type(ctx, p, d, var, 0);
}

void fq_default_ctx_init_modulus_type(fq_default_ctx_t ctx,
                const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx,
                                                   const char * var, int type);

void fq_default_ctx_init_modulus(fq_default_ctx_t ctx,
      const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var);


void fq_default_ctx_init_modulus_nmod_type(fq_default_ctx_t ctx,
                        const nmod_poly_t modulus, const char * var, int type);


void fq_default_ctx_init_modulus_nmod(fq_default_ctx_t ctx,
                                  const nmod_poly_t modulus, const char * var);

FQ_DEFAULT_INLINE void fq_default_ctx_clear(fq_default_ctx_t ctx)
{
    if (!ctx->_gr_inited)
        flint_abort();

    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_ctx_clear(FQ_DEFAULT_CTX_FMPZ_MOD_OLD(ctx).mod);
        fmpz_clear(FQ_DEFAULT_CTX_FMPZ_MOD_OLD(ctx).a);
    }

    if (ctx->_gr_inited)
        gr_ctx_clear(FQ_DEFAULT_GR_CTX(ctx));
}

FQ_DEFAULT_INLINE int fq_default_ctx_type(const fq_default_ctx_t ctx)
{
    return FQ_DEFAULT_TYPE(ctx);
}

FQ_DEFAULT_INLINE slong fq_default_ctx_degree(const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_ctx_degree(FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_ctx_degree(FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        return 1;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        return 1;
    }
    else
    {
        return fq_ctx_degree(FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_INLINE void fq_default_ctx_prime(fmpz_t prime,
                                                    const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fmpz_set_ui(prime, fq_zech_ctx_prime(FQ_DEFAULT_CTX_FQ_ZECH(ctx)));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fmpz_set_ui(prime, fq_nmod_ctx_prime(FQ_DEFAULT_CTX_FQ_NMOD(ctx)));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        fmpz_set_ui(prime, FQ_DEFAULT_CTX_NMOD(ctx).n);
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_set(prime, fmpz_mod_ctx_modulus(FQ_DEFAULT_CTX_FMPZ_MOD(ctx)));
    }
    else
    {
        fmpz_set(prime, fq_ctx_prime(FQ_DEFAULT_CTX_FQ(ctx)));
    }
}

void fq_default_ctx_modulus(fmpz_mod_poly_t p,
		                                   const fq_default_ctx_t ctx);

FQ_DEFAULT_INLINE void fq_default_ctx_order(fmpz_t f,
		                                    const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fmpz_set_ui(f, fq_zech_ctx_order_ui(FQ_DEFAULT_CTX_FQ_ZECH(ctx)));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_ctx_order(f, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        fmpz_set_ui(f, FQ_DEFAULT_CTX_NMOD(ctx).n);
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_set(f, fmpz_mod_ctx_modulus(FQ_DEFAULT_CTX_FMPZ_MOD(ctx)));
    }
    else
    {
        fq_ctx_order(f, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

#ifdef FLINT_HAVE_FILE
int fq_default_ctx_fprint(FILE * file, const fq_default_ctx_t ctx);
#endif

void fq_default_ctx_print(const fq_default_ctx_t ctx);

/* Memory management  *********************************************************/

FQ_DEFAULT_INLINE void fq_default_init(fq_default_t rop,
		                                    const fq_default_ctx_t ctx)
{
    gr_init(rop, FQ_DEFAULT_GR_CTX(ctx));
}

FQ_DEFAULT_INLINE void fq_default_init2(fq_default_t rop,
		                                    const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_init2(rop->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_init2(rop->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        rop->nmod = 0;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_init(rop->fmpz_mod);
    }
    else
    {
        fq_init2(rop->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_INLINE void fq_default_clear(fq_default_t rop,
		                                    const fq_default_ctx_t ctx)
{
    gr_clear(rop, FQ_DEFAULT_GR_CTX(ctx));
}

/* Predicates ****************************************************************/

FQ_DEFAULT_INLINE int fq_default_is_invertible(const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_is_invertible(op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_is_invertible(op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        return op->nmod != 0;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        return !fmpz_is_zero(op->fmpz_mod);
    }
    else
    {
        return fq_is_invertible(op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/* Basic arithmetic **********************************************************/

/* GR_IGNORE(...) because we assume that basic arithmetic is well-implemented */

FQ_DEFAULT_INLINE void fq_default_add(fq_default_t rop, const fq_default_t op1,
		            const fq_default_t op2, const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_add(rop, op1, op2, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_sub(fq_default_t rop, const fq_default_t op1,
		            const fq_default_t op2, const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_sub(rop, op1, op2, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_sub_one(fq_default_t rop,
		            const fq_default_t op1, const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_sub_ui(rop, op1, 1, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_neg(fq_default_t rop,
		            const fq_default_t op1, const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_neg(rop, op1, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_mul(fq_default_t rop, const fq_default_t op1,
                            const fq_default_t op2, const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_mul(rop, op1, op2, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_mul_fmpz(fq_default_t rop,
             const fq_default_t op, const fmpz_t x, const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_mul_fmpz(rop, op, x, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_mul_si(fq_default_t rop,
		    const fq_default_t op, slong x, const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_mul_si(rop, op, x, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_mul_ui(fq_default_t rop,
		    const fq_default_t op, ulong x, const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_mul_ui(rop, op, x, FQ_DEFAULT_GR_CTX(ctx)));
}

/* todo: verify all implementations define mul_fmpz, mul_si, mul_ui, sqr, pow_fmpz, pow_ui */

FQ_DEFAULT_INLINE void fq_default_sqr(fq_default_t rop,
                             const fq_default_t op, const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_sqr(rop, op, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_inv(fq_default_t rop,
		            const fq_default_t op1, const fq_default_ctx_t ctx)
{
    GR_MUST_SUCCEED(gr_inv(rop, op1, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_div(fq_default_t rop, fq_default_t op1,
                                  fq_default_t op2, const fq_default_ctx_t ctx)
{
    GR_MUST_SUCCEED(gr_div(rop, op1, op2, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_pow(fq_default_t rop,
	    const fq_default_t op1, const fmpz_t e, const fq_default_ctx_t ctx)
{
    GR_MUST_SUCCEED(gr_pow_fmpz(rop, op1, e, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_pow_ui(fq_default_t rop,
              const fq_default_t op, const ulong e, const fq_default_ctx_t ctx)
{
    GR_MUST_SUCCEED(gr_pow_ui(rop, op, e, FQ_DEFAULT_GR_CTX(ctx)));
}

/* Roots *********************************************************************/

FQ_DEFAULT_INLINE int fq_default_sqrt(fq_default_t rop,
		             const fq_default_t op, const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_sqrt(rop->fq_zech, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_sqrt(rop->fq_nmod, op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        if (op->nmod == 0)
        {
            rop->nmod = 0;
            return 1;
        }
        else
        {
            rop->nmod = n_sqrtmod(op->nmod, FQ_DEFAULT_CTX_NMOD(ctx).n);
            return rop->nmod != 0;
        }
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_sqrtmod(rop->fmpz_mod, op->fmpz_mod,
                                  fmpz_mod_ctx_modulus(FQ_DEFAULT_CTX_FMPZ_MOD(ctx)));
    }
    else
    {
        return fq_sqrt(rop->fq, op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_INLINE void fq_default_pth_root(fq_default_t rop,
		            const fq_default_t op1, const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_pth_root(rop->fq_zech, op1->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_pth_root(rop->fq_nmod, op1->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        rop->nmod = op1->nmod;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_set(rop->fmpz_mod, op1->fmpz_mod);
    }
    else
    {
        fq_pth_root(rop->fq, op1->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_INLINE int fq_default_is_square(const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_is_square(op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_is_square(op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        return op->nmod == 0 || n_sqrtmod(op->nmod, FQ_DEFAULT_CTX_NMOD(ctx).n) != 0;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        int res;
        fmpz_t t;
        fmpz_init(t);
        res = fmpz_sqrtmod(t, op->fmpz_mod, fmpz_mod_ctx_modulus(FQ_DEFAULT_CTX_FMPZ_MOD(ctx)));
        fmpz_clear(t);
        return res;
    }
    else
    {
        return fq_is_square(op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/* Randomisation *************************************************************/

FQ_DEFAULT_INLINE void fq_default_randtest(fq_default_t rop,
		                flint_rand_t state, const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_randtest(rop->fq_zech, state, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_randtest(rop->fq_nmod, state, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        rop->nmod = n_randint(state, FQ_DEFAULT_CTX_NMOD(ctx).n);
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_rand(rop->fmpz_mod, state, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_randtest(rop->fq, state, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_INLINE void fq_default_randtest_not_zero(fq_default_t rop,
		                flint_rand_t state, const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_randtest_not_zero(rop->fq_zech, state, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_randtest_not_zero(rop->fq_nmod, state, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        rop->nmod = n_randint(state, FQ_DEFAULT_CTX_NMOD(ctx).n - 1) + 1;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_rand_not_zero(rop->fmpz_mod, state, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_randtest_not_zero(rop->fq, state, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_INLINE void fq_default_rand(fq_default_t rop,
		                flint_rand_t state, const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_rand(rop->fq_zech, state, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_rand(rop->fq_nmod, state, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        rop->nmod = n_randint(state, FQ_DEFAULT_CTX_NMOD(ctx).n);
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_rand(rop->fmpz_mod, state, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_rand(rop->fq, state, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_INLINE void fq_default_rand_not_zero(fq_default_t rop,
		                flint_rand_t state, const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_rand_not_zero(rop->fq_zech, state, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_rand_not_zero(rop->fq_nmod, state, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        rop->nmod = n_randint(state, FQ_DEFAULT_CTX_NMOD(ctx).n - 1) + 1;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_mod_rand_not_zero(rop->fmpz_mod, state, FQ_DEFAULT_CTX_FMPZ_MOD(ctx));
    }
    else
    {
        fq_rand_not_zero(rop->fq, state, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/* Comparison ****************************************************************/

FQ_DEFAULT_INLINE int fq_default_equal(const fq_default_t op1,
		            const fq_default_t op2, const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_equal(op1->fq_zech, op2->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_equal(op1->fq_nmod, op2->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        return op1->nmod == op2->nmod;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_equal(op1->fmpz_mod, op2->fmpz_mod);
    }
    else
    {
        return fq_equal(op1->fq, op2->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_INLINE int fq_default_is_zero(const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_is_zero(op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_is_zero(op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        return op->nmod == 0;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_is_zero(op->fmpz_mod);
    }
    else
    {
        return fq_is_zero(op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_INLINE int fq_default_is_one(const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_is_one(op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_is_one(op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        return op->nmod == 1;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_is_one(op->fmpz_mod);
    }
    else
    {
        return fq_is_one(op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/* Assignments and conversions ***********************************************/

FQ_DEFAULT_INLINE void fq_default_set(fq_default_t rop,
		             const fq_default_t op, const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_set(rop, op, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_set_fmpz(fq_default_t rop,
		                    const fmpz_t x, const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_set_fmpz(rop, x, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_set_ui(fq_default_t rop,
		                     const ulong x, const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_set_ui(rop, x, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_set_si(fq_default_t rop,
		                     const slong x, const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_set_si(rop, x, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_zero(fq_default_t rop,
		                                    const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_zero(rop, FQ_DEFAULT_GR_CTX(ctx)));
}

FQ_DEFAULT_INLINE void fq_default_one(fq_default_t rop,
		                                    const fq_default_ctx_t ctx)
{
    GR_IGNORE(gr_one(rop, FQ_DEFAULT_GR_CTX(ctx)));
}

/* todo: should this swap the fq_default_structs? */
FQ_DEFAULT_INLINE void fq_default_swap(fq_default_t op1,
		                  fq_default_t op2, const fq_default_ctx_t ctx)
{
    gr_swap(op1, op2, FQ_DEFAULT_GR_CTX(ctx));
}

FQ_DEFAULT_INLINE void fq_default_gen(fq_default_t rop,
		                                    const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_gen(rop->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_gen(rop->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        rop->nmod = *FQ_DEFAULT_CTX_NMOD_A(ctx);
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_set(rop->fmpz_mod, FQ_DEFAULT_CTX_FMPZ_MOD_A(ctx));
    }
    else
    {
        fq_gen(rop->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_INLINE void fq_default_get_nmod_poly(nmod_poly_t poly,
                             const fq_default_t op, const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_get_nmod_poly(poly, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_get_nmod_poly(poly, op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        nmod_poly_fit_length(poly, 1);
        poly->length = (op->nmod != 0);
        poly->coeffs[0] = op->nmod;
    }
    else
    {
        flint_throw(FLINT_ERROR, "Impossible conversion\n");
    }
}

FQ_DEFAULT_INLINE void fq_default_set_nmod_poly(fq_default_t op,
                            const nmod_poly_t poly, const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_set_nmod_poly(op->fq_zech, poly, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_set_nmod_poly(op->fq_nmod, poly, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        op->nmod = nmod_poly_evaluate_nmod(poly, *FQ_DEFAULT_CTX_NMOD_A(ctx));
    }
    else
    {
        flint_throw(FLINT_ERROR, "Impossible conversion\n");
    }
}

FQ_DEFAULT_INLINE int fq_default_get_fmpz(fmpz_t z, const fq_default_t op,
                                                    const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_get_fmpz(z, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_get_fmpz(z, op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        fmpz_set_ui(z, op->nmod);
        return 1;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_set(z, op->fmpz_mod);
        return 1;
    }
    else
    {
        return fq_get_fmpz(z, op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

void fq_default_get_fmpz_mod_poly(fmpz_mod_poly_t poly,
                            const fq_default_t op, const fq_default_ctx_t ctx);

void fq_default_set_fmpz_mod_poly(fq_default_t op,
                       const fmpz_mod_poly_t poly, const fq_default_ctx_t ctx);

void fq_default_get_fmpz_poly(fmpz_poly_t poly,
                            const fq_default_t op, const fq_default_ctx_t ctx);

void fq_default_set_fmpz_poly(fq_default_t op,
                           const fmpz_poly_t poly, const fq_default_ctx_t ctx);

FQ_DEFAULT_INLINE void fq_default_get_coeff_fmpz(fmpz_t c,
                          fq_default_t op, slong n, const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        nmod_poly_t p;
        ulong c0;
        nmod_poly_init(p, fq_zech_ctx_prime(FQ_DEFAULT_CTX_FQ_ZECH(ctx)));
        fq_zech_get_nmod_poly(p, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
        c0 = nmod_poly_get_coeff_ui(p, n);
        fmpz_set_ui(c, c0);
        nmod_poly_clear(p);
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        ulong c0 = nmod_poly_get_coeff_ui((nmod_poly_struct *) op->fq_nmod, n);
        fmpz_set_ui(c, c0);
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        if (n != 0)
            fmpz_zero(c);
        else
            fmpz_set_ui(c, op->nmod);
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        if (n != 0)
            fmpz_zero(c);
        else
            fmpz_set(c, op->fmpz_mod);
    }
    else
    {
        fmpz_mod_ctx_t mod_ctx;
        fmpz_mod_ctx_init(mod_ctx, fq_ctx_prime(FQ_DEFAULT_CTX_FQ(ctx)));
        fmpz_mod_poly_get_coeff_fmpz(c,
		                  (fmpz_mod_poly_struct *) op->fq, n, mod_ctx);
        fmpz_mod_ctx_clear(mod_ctx);
    }
}

/* Output ********************************************************************/

#ifdef FLINT_HAVE_FILE
int fq_default_fprint(FILE * file, const fq_default_t op, const fq_default_ctx_t ctx);
int fq_default_fprint_pretty(FILE * file, const fq_default_t op, const fq_default_ctx_t ctx);
#endif

void fq_default_print(const fq_default_t op, const fq_default_ctx_t ctx);
void fq_default_print_pretty(const fq_default_t op, const fq_default_ctx_t ctx);

FQ_DEFAULT_INLINE char * fq_default_get_str(const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_get_str(op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_get_str(op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        fmpz_t f;
        char* s;
        fmpz_init_set_ui(f, op->nmod);
        s = fmpz_get_str(NULL, 10, f);
        fmpz_clear(f);
        return s;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_get_str(NULL, 10, op->fmpz_mod);
    }
    else
    {
        return fq_get_str(op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_INLINE char * fq_default_get_str_pretty(const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        return fq_zech_get_str_pretty(op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        return fq_nmod_get_str_pretty(op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        fmpz_t f;
        char* s;
        fmpz_init_set_ui(f, op->nmod);
        s = fmpz_get_str(NULL, 10, f);
        fmpz_clear(f);
        return s;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        return fmpz_get_str(NULL, 10, op->fmpz_mod);
    }
    else
    {
        return fq_get_str_pretty(op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

/* Special functions *********************************************************/

FQ_DEFAULT_INLINE void fq_default_trace(fmpz_t rop,
		             const fq_default_t op, const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_trace(rop, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_trace(rop, op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        fmpz_set_ui(rop, op->nmod);
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_set(rop, op->fmpz_mod);
    }
    else
    {
        fq_trace(rop, op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_INLINE void fq_default_frobenius(fq_default_t rop,
		    const fq_default_t op, slong e, const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_frobenius(rop->fq_zech, op->fq_zech, e, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_frobenius(rop->fq_nmod, op->fq_nmod, e, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        rop->nmod = op->nmod;
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_set(rop->fmpz_mod, op->fmpz_mod);
    }
    else
    {
        fq_frobenius(rop->fq, op->fq, e, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

FQ_DEFAULT_INLINE void fq_default_norm(fmpz_t rop,
		             const fq_default_t op, const fq_default_ctx_t ctx)
{
    if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_ZECH)
    {
        fq_zech_norm(rop, op->fq_zech, FQ_DEFAULT_CTX_FQ_ZECH(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FQ_NMOD)
    {
        fq_nmod_norm(rop, op->fq_nmod, FQ_DEFAULT_CTX_FQ_NMOD(ctx));
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_NMOD)
    {
        fmpz_set_ui(rop, op->nmod);
    }
    else if (FQ_DEFAULT_TYPE(ctx) == FQ_DEFAULT_FMPZ_MOD)
    {
        fmpz_set(rop, op->fmpz_mod);
    }
    else
    {
        fq_norm(rop, op->fq, FQ_DEFAULT_CTX_FQ(ctx));
    }
}

#ifdef __cplusplus
}
#endif

#endif
