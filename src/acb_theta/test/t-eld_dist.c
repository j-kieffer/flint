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
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_eld_dist, state)
{
    slong iter;

    /* Test: ellipsoid must have points, and the points are not too close */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = ACB_THETA_LOW_PREC;
        slong hprec = 200;
        slong bits = n_randint(state, 5);
        acb_mat_t tau;
        arb_mat_t cho, yinv;
        acb_ptr z;
        arb_ptr v;
        arb_t d;
        acb_theta_eld_t E;
        arf_t R2;
        int b;
        int omit_zero = iter % 2;

        acb_mat_init(tau, g, g);
        arb_mat_init(cho, g, g);
        arb_mat_init(yinv, g, g);
        z = _acb_vec_init(g);
        v = _arb_vec_init(g);
        acb_theta_eld_init(E, g, g);
        arf_init(R2);

        acb_siegel_randtest_reduced(tau, state, hprec, bits);
        acb_siegel_cho_yinv(cho, yinv, tau, hprec);
        acb_siegel_randtest_vec(z, state, g, hprec);
        _acb_vec_get_real(v, z, g);

        acb_theta_eld_dist(d, v, cho, omit_zero, prec);

        arb_get_ubound_arf(R2, d, prec);
        b = acb_theta_eld_set(E, cho, R2, v);

        if (b && acb_theta_eld_nb_pts(E) == 0)
        {
            flint_printf("FAIL (no points)\n");
            flint_printf("g = %wd, omit_zero = %wd, d = ", g, omit_zero);
            arb_printd(d, 5);
            flint_printf("\n");
            arb_mat_printd(cho, 5);
            _arb_vec_printd(v, g, 5);
            flint_abort();
        }

        if (b)
        {
            slong * pts;
            slong k;
            arb_t t;
            arb_t dmin;

            pts = flint_malloc(g * acb_theta_eld_nb_pts(E) * sizeof(slong));
            arb_init(t);
            arb_init(dmin);

            arb_pos_inf(dmin);
            acb_theta_eld_points(pts, E);

            for (k = 0; k < acb_theta_eld_nb_pts(E); k++)
            {
                if (omit_zero && acb_theta_eld_pt_is_zero(pts + k * g, g))
                {
                    continue;
                }
                acb_theta_eld_dist_pt(t, v, cho, pts + k * g, prec);
                arb_min(dmin, dmin, t, prec);
                if (arb_lt(t, d))
                {
                    flint_printf("FAIL (point too close)\n");
                    flint_printf("g = %wd, omit_zero = %wd, d = ", g, omit_zero);
                    arb_printd(d, 5);
                    flint_printf("\nt = ");
                    arb_printd(t, 5);
                    flint_printf("\n");
                    arb_mat_printd(cho, 5);
                    _arb_vec_printd(v, g, 5);
                    for (slong j = 0; j < g; j++)
                    {
                        flint_printf(" %wd", pts[k * g + j]);
                    }
                    flint_printf("\n");
                    flint_abort();
                }
            }

            if (arb_gt(dmin, d))
            {
                flint_printf("FAIL (points too far)\n");
            }

            flint_free(pts);
            arb_clear(t);
            arb_clear(dmin);
        }

        acb_mat_clear(tau);
        arb_mat_clear(cho);
        arb_mat_clear(yinv);
        _acb_vec_clear(z, g);
        _arb_vec_clear(v, g);
        acb_theta_eld_clear(E);
        arf_clear(R2);
    }

    TEST_FUNCTION_END(state);
}
