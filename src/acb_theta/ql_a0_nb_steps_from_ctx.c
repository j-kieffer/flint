/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong
acb_theta_ql_a0_nb_steps_from_ctx(slong * s, const acb_theta_ql_ctx_t ctx, slong prec)
{
    slong g = acb_theta_ql_ctx_g(ctx);
    slong lp = ACB_THETA_LOW_PREC;
    arb_t x, t;
    slong res;

    arb_init(x);
    arb_init(t);

	/* Compute s */
    for (*s = g - 1; *s >= 1; *s--)
    {
        arb_mul_2exp_si(x, arb_mat_entry(acb_theta_ql_ctx_y(ctx), *s - 1, *s - 1),
            FLINT_MAX(1, 3 + *s - g));
        if (arb_lt(cmp, arb_mat_entry(acb_theta_ql_ctx_y(ctx), *s, *s)))
        {
            break;
        }
    }

    FLINT_ASSERT(*s >= 0 && *s < g);

	/* Decide number of steps heuristically */
    arb_sqr(x, arb_mat_entry(C, *s, *s), lp);
    arb_const_log2(t, lp);
    arb_div(x, x, t, lp);
    arb_div_si(x, x, prec, lp);
    arb_log(x, x, lp);
    arb_div(x, x, t, lp);

    if (!arb_is_finite(x) || arf_cmpabs_2exp_si(arb_midref(x), FLINT_BITS - 4) > 0)
    {
        arb_clear(x);
        arb_clear(t);
        return 0;
    }

    res =  -arf_get_si(arb_midref(x), ARF_RND_NEAR);
    if (*s == 0)
    {
        if (g == 1)
        {
            res -= 7;
        }
        else if (g == 2)
        {
            res -= 3;
        }
        else if (g <= 5)
        {
            res -= 1;
        }
    }
    else
    {
        res += 1;
    }
    res = FLINT_MAX(0, res);

    arb_clear(x);
    arb_clear(t);
    return res;
}
