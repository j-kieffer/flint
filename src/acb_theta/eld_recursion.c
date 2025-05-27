/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_mat.h"
#include "acb_theta.h"

void
acb_theta_eld_recursion(arb_t next_R2, arb_ptr next_v, const arb_mat_t cho,
    const arb_t R2, arb_srcptr v, slong n, slong g, slong prec)
{
    arb_t x;
    slong k;

    arb_init(x);

    /* Set next_R2 to R2 - (v + gamma * n)^2 */
    arb_mul_si(x, arb_mat_entry(cho, g - 1, g - 1), n, prec);
    arb_add(x, x, &v[g - 1], prec);
    arb_sqr(x, x, prec);
    arb_sub(next_R2, R2, x, prec);
    arb_nonnegative_part(next_R2, next_R2);

    /* Set next_v */
    for (k = 0; k < g - 1; k++)
    {
        arb_set(&next_v[k], &v[k]);
        arb_addmul_si(&next_v[k], arb_mat_entry(cho, k, g - 1), n, prec);
    }

    arb_clear(x);
}
