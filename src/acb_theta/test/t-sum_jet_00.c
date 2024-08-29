/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_sum_jet_00, state)
{
    slong iter;

    /* Test: matches jet_naive_00 */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong prec = 100 + n_randint(state, 100);
        slong ord = n_randint(state, 4);
        slong mag_bits = n_randint(state, 4);
        slong nbz = 1 + n_randint(state, 4);
        slong nbth = acb_theta_jet_nb(ord, g);
        acb_mat_t tau;
        acb_ptr zs;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_struct * vec;
        acb_ptr th1, th2;
        slong j;

        acb_mat_init(tau, g, g);
        zs = _acb_vec_init(nbz * g);
        acb_theta_ctx_tau_init(ctx_tau, g);
        vec = acb_theta_ctx_z_vec_init(nbz, g);
        th1 = _acb_vec_init(nbz * nbth);
        th2 = _acb_vec_init(nbz * nbth);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_siegel_randtest_vec(zs, state, nbz * g, prec);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        for (j = 0; j < nbz; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
        }

        acb_theta_sum_jet_00(th1, vec, nbz, ctx_tau, ord, prec);
        for (j = 0; j < nbz; j++)
        {
            acb_theta_jet_naive_00(th2 + j * nbth, zs + j * g, tau, ord, prec);
        }

        if (!_acb_vec_overlaps(th1, th2, nbz * nbth))
        {
            flint_printf("FAIL\n");
            flint_printf("\n\ng=%wd, ord=%wd\n", g, ord);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(zs, nbz * g, 5);
            flint_printf("th1 : ");
            _acb_vec_printd(th1, nbz * nbth, 5);
            flint_printf("th2 : ");
            _acb_vec_printd(th2, nbz * nbth, 5);
            flint_printf("Difference: ");
            _acb_vec_sub(th1, th1, th2, nbz * nbth, prec);
            _acb_vec_printd(th1, nbz * nbth, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(zs, nbz * g);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(vec, nbz);
        _acb_vec_clear(th1, nbz * nbth);
        _acb_vec_clear(th2, nbz * nbth);
    }

    TEST_FUNCTION_END(state);
}