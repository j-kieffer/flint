/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_eld_interval, state)
{
    slong iter;

    /* Test: check that min, mid, max satisfy the required inequalities */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong prec = ACB_THETA_LOW_PREC;
        slong bits = 1 + n_randint(state, 5);
        acb_mat_t tau;
        arb_mat_t cho, yinv;
        acb_ptr z;
        arb_ptr v, next_v;
        arf_t R2;
        arb_t x, t;
        slong min, mid, max;
        int res;

        acb_mat_init(tau, g, g);
        arb_mat_init(cho, g, g);
        arb_mat_init(yinv, g, g);
        z = _acb_vec_init(g);
        v = _arb_vec_init(g);
        next_v = _arb_vec_init(g - 1);
        arf_init(R2);
        arb_init(x);
        arb_init(t);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_siegel_cho_yinv(cho, yinv, tau, prec);
        acb_siegel_randtest_vec(z, state, g, prec);
        _acb_vec_get_real(v, z, g);
        arf_randtest_special(R2, state, bits, 2);

        res = acb_theta_eld_interval(&min, &mid, &max,
            arb_mat_entry(cho, g - 1, g - 1), R2, &v[g - 1], prec);

        if (res)
        {
            /* min - 1 and max + 1 do not satisfy the inequality */
            arb_set_arf(x, R2);
            acb_theta_eld_recursion(t, next_v, cho, x, v, min - 1, g, prec);

            if (arb_is_positive(t))
            {
                flint_printf("FAIL (min - 1)\n");
                flint_abort();
            }

            acb_theta_eld_recursion(t, next_v, cho, x, v, max + 1, g, prec);

            if (arb_is_positive(t))
            {
                flint_printf("FAIL (max + 1)\n");
                flint_abort();
            }

            if (min <= max)
            {
                /* mid is between them */
                if (min > mid || mid > max)
                {
                    flint_printf("FAIL (inequalities)\n");
                    flint_abort();
                }
            }
        }

        acb_mat_clear(tau);
        arb_mat_clear(cho);
        arb_mat_clear(yinv);
        _acb_vec_clear(z, g);
        _arb_vec_clear(v, g);
        _arb_vec_clear(next_v, g - 1);
        arf_clear(R2);
        arb_clear(x);
        arb_clear(t);
    }

    TEST_FUNCTION_END(state);
}

