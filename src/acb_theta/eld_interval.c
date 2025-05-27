/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb_theta.h"

#define ACB_THETA_ELD_MAX_ERR 100

static int
acb_theta_arf_get_si_safe(slong * m, const arf_t x, arf_rnd_t rnd)
{
    if (!arf_is_finite(x))
    {
        return 0;
    }
    else if (arf_cmpabs_2exp_si(x, FLINT_BITS - 4) > 0)
    {
        return 0;
    }
    else
    {
        *m = arf_get_si(x, rnd);
        return 1;
    }
}

int
acb_theta_eld_interval(slong * min, slong * mid, slong * max,
    const arb_t c, const arf_t R2, const arb_t v, slong prec)
{
    slong e;
    arb_t ctr, y;
    arf_t rad, b;
    int res;

    arb_init(ctr);
    arb_init(y);
    arf_init(rad);
    arf_init(b);

    /* Compute center and mid */
    arb_neg(ctr, v);
    arb_div(ctr, ctr, c, prec);

    arf_set_mag(b, arb_radref(ctr));
    res = acb_theta_arf_get_si_safe(&e, b, ARF_RND_NEAR);
    if (res)
    {
        res = (e <= ACB_THETA_ELD_MAX_ERR);
    }

    res = res && acb_theta_arf_get_si_safe(mid, arb_midref(ctr), ARF_RND_NEAR);

    if (arf_cmp_si(R2, 0) < 0)
    {
        *min = *mid + 1;
        *max = *mid;
    }
    else
    {
        arb_set_arf(y, R2);
        arb_sqrt(y, y, prec);
        arb_div(y, y, c, prec);
        arb_get_ubound_arf(rad, y, prec);

        arb_set_arf(y, rad);
        arb_add(y, ctr, y, prec);
        arb_get_ubound_arf(b, y, prec);
        res = res && acb_theta_arf_get_si_safe(max, b, ARF_RND_FLOOR);

        arb_set_arf(y, rad);
        arb_sub(y, ctr, y, prec);
        arb_get_lbound_arf(b, y, prec);
        res = res && acb_theta_arf_get_si_safe(min, b, ARF_RND_CEIL);
    }

    arb_clear(ctr);
    arb_clear(y);
    arf_clear(rad);
    arf_clear(b);
    return res;
}
