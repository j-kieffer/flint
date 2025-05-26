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

static int
acb_theta_slong_vec_is_zero(const slong * n, slong g)
{
    slong k;

    for (k = 0; k < g; k++)
    {
        if (n[k] != 0)
        {
            return 0;
        }
    }
    return 1;
}

static void
acb_theta_eld_dist_unif(arb_t d, const arb_mat_t cho, int omit_zero, slong prec)
{
    slong g = arb_mat_nrows(cho);
    arb_ptr v;
    slong k;

    v = _arb_vec_init(g);

    for (k = 0; k < g; k++)
    {
        arb_zero_pm_one(&v[k]);
        if (!omit_zero)
        {
            arb_mul_2exp_si(&v[k], &v[k], -1);
        }
    }
    arb_mat_vector_mul_col(v, cho, v, prec);
    arb_dot(d, NULL, 0, v, 1, v, 1, g, prec);

    _arb_vec_clear(v, g);
}

static void
acb_theta_eld_dist_pt(arb_t d, arb_srcptr v, const arb_mat_t cho,
    const slong * n, slong prec)
{
    slong g = arb_mat_nrows(cho);
    arb_ptr w;
    slong k;

    w = _arb_vec_init(g);

    for (k = 0; k < g; k++)
    {
        arb_set_si(&w[k], n[k]);
    }
    arb_mat_vector_mul_col(w, cho, w, prec);
    _arb_vec_add(w, w, v, g, prec);
    arb_dot(d, NULL, 0, w, 1, w, 1, g, prec);

    _arb_vec_clear(w, g);
}

static void
acb_theta_eld_dist_ubound(arf_t u, arb_srcptr v, const arb_mat_t cho,
    int omit_zero, slong prec)
{
    slong g = arb_mat_nrows(cho);
    slong nb = 1 << g;
    arb_mat_t choinv;
    arb_ptr x;
    slong * approx_up;
    slong * pt;
    arb_t d;
    arf_t b;
    slong j, k;
    int r = 1;

    arb_mat_init(choinv, g, g);
    x = _arb_vec_init(g);
    approx_up = flint_malloc(g * sizeof(slong));
    pt = flint_malloc(g * sizeof(slong));
    arb_init(d);
    arf_init(b);

    arb_mat_one(choinv);
    arb_mat_solve_triu(choinv, cho, choinv, 0, prec);
    arb_mat_vector_mul_col(x, choinv, v, prec);
    r = _arb_vec_is_finite(x, g);

    for (k = 0; (k < g) && r; k++)
    {
        r = (arf_cmpabs_2exp_si(arb_midref(&x[k]), 30) <= 0);
        if (r)
        {
            approx_up[k] = - arf_get_si(arb_midref(&x[k]), ARF_RND_FLOOR);
        }
    }

    arf_pos_inf(u);
    if (r)
    {
        for (k = 0; k < nb; k++)
        {
            for (j = 0; j < g; j++)
            {
                pt[j] = approx_up[j];
                if (k & (1 << j))
                {
                    pt[j] -= 1;
                }
            }
            if (!omit_zero || !acb_theta_slong_vec_is_zero(pt, g))
            {
                acb_theta_eld_dist_pt(d, v, cho, pt, prec);
                arb_get_ubound_arf(b, d, prec);
                arf_min(u, u, b);
            }
        }
    }

    arb_mat_clear(choinv);
    _arb_vec_clear(x, g);
    flint_free(approx_up);
    flint_free(pt);
    arb_clear(d);
    arf_clear(b);
}

void
acb_theta_eld_dist(arb_t d, arb_srcptr v, const arb_mat_t cho, int omit_zero, slong prec)
{
    slong g = arb_mat_nrows(cho);
    acb_theta_eld_t E;
    slong nb;
    slong * pts;
    arf_t u;
    arb_t x;
    slong k;
    int b;

    acb_theta_eld_init(E, g, g);
    arf_init(u);
    arb_init(x);

    acb_theta_eld_dist_ubound(u, v, cho, omit_zero, prec);
    b = acb_theta_eld_set(E, cho, u, v);

    if (b)
    {
        nb = acb_theta_eld_nb_pts(E);
        pts = flint_malloc(nb * g * sizeof(slong));
        acb_theta_eld_points(pts, E);

        arb_pos_inf(d);
        for (k = 0; k < nb; k++)
        {
            if (!omit_zero || !acb_theta_slong_vec_is_zero(pts + k * g, g))
            {
                acb_theta_eld_dist_pt(x, v, cho, pts + k * g, prec);
                arb_min(d, d, x, prec);
            }
        }

        flint_free(pts);
    }
    else
    {
        acb_theta_eld_dist_unif(d, cho, omit_zero, prec);
    }
    arb_nonnegative_part(d, d);

    acb_theta_eld_clear(E);
    arf_clear(u);
    arb_clear(x);
}
