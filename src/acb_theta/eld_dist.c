/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_mat.h"
#include "acb_theta.h"

#define ACB_THETA_ELD_MAX_ERR 100

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
acb_theta_eld_next_dmax(arb_t next_dmax, const arb_t dmax, const arb_t gamma,
    const arb_t v, slong k, slong prec)
{
    arb_t x;

    arb_init(x);

    /* Set next_dmax to dmax - (v + gamma*k)^2 */
    arb_mul_si(x, gamma, k, prec);
    arb_add(x, x, v, prec);
    arb_sqr(x, x, prec);
    arb_sub(next_dmax, dmax, x, prec);
    arb_nonnegative_part(next_dmax, next_dmax);

    arb_clear(x);
}

static void
acb_theta_eld_dist_rec(arb_t d, const arb_mat_t cho, arb_srcptr v,
    const arb_t dmax, slong g, int omit_zero, slong prec)
{
    arb_t x, ctr;
    arf_t rad;
    slong min, mid, max;
    int res = 1;

    arb_init(x);
    arb_init(ctr);
    arf_init(rad);

    /* Get ctr and rad for possible values of n_{g-1} */
    arb_sqrt(x, dmax, prec);
    arb_div(x, x, arb_mat_entry(cho, g - 1, g - 1), prec);
    arb_get_ubound_arf(rad, x, prec);
    arb_div(ctr, &v[g - 1], arb_mat_entry(cho, g - 1, g - 1), prec);
    arb_neg(ctr, ctr);
    res = acb_theta_eld_interval(&min, &mid, &max, ctr, rad);

    if (res && min > max)
    {
        /* Distance is bigger than dmax */
        arb_set(d, dmax);
    }
    else if (res && g == 1)
    {
        /* Can compute distance directly */
        if (!omit_zero || mid != 0)
        {
            arb_sub_si(x, ctr, mid, prec);
            arb_sqr(d, x, prec);
        }
        else
        {
            arb_sub_si(x, ctr, mid + 1, prec);
            arb_sqr(d, x, prec);
        }

        if (!omit_zero || mid + 1 != 0)
        {
            arb_sub_si(x, ctr, mid + 1, prec);
            arb_sqr(x, x, prec);
            arb_min(d, d, x, prec);
        }
        if (!omit_zero || mid - 1 != 0)
        {
            arb_sub_si(x, ctr, mid - 1, prec);
            arb_sqr(x, x, prec);
            arb_min(d, d, x, prec);
        }
        arb_nonnegative_part(d, d);
        arb_sqr(x, arb_mat_entry(cho, 0, 0), prec);
        arb_mul(d, d, x, prec);
    }
    else if (res)
    {
        /* Induction */
        arb_ptr v_mid, v_diff, next_v;
        arb_t next_dmax;
        slong k;

        v_mid = _arb_vec_init(g - 1);
        v_diff = _arb_vec_init(g - 1);
        next_v = _arb_vec_init(g - 1);
        arb_init(next_dmax);

        arb_set(d, dmax);

        /* Set v_mid, v_diff */
        for (k = 0; k < g - 1; k++)
        {
            arb_set(&v_diff[k], arb_mat_entry(cho, k, g - 1));
            arb_mul_si(&v_mid[k], &v_diff[k], mid, prec);
        }
        _arb_vec_add(v_mid, v_mid, v, g - 1, prec);

        /* Right loop */
        _arb_vec_set(next_v, v_mid, g - 1);
        k = 0;
        while (res && mid + k < max + 1)
        {
            acb_theta_eld_next_dmax(next_dmax, d,
                arb_mat_entry(cho, g - 1, g - 1), &v[g - 1], mid + k, prec);
            acb_theta_eld_dist_rec(x, cho, next_v, next_dmax, g - 1,
                omit_zero && (mid + k == 0), prec);
            arb_add(x, x, d, prec);
            arb_sub(x, x, next_dmax, prec);
            arb_min(d, d, x, prec);
            arb_nonnegative_part(d, d);

            /* Recompute min, max */
            arb_sqrt(x, d, prec);
            arb_div(x, x, arb_mat_entry(cho, g - 1, g - 1), prec);
            arb_get_ubound_arf(rad, x, prec);
            res = acb_theta_eld_interval(&min, &mid, &max, ctr, rad);

            k++;
            _arb_vec_add(next_v, next_v, v_diff, g - 1, prec);
        }

        /* Left loop */
        _arb_vec_set(next_v, v_mid, g - 1);
        k = 0;
        while (res && mid - (k + 1) > min - 1)
        {
            _arb_vec_sub(next_v, next_v, v_diff, g - 1, prec);

            acb_theta_eld_next_dmax(next_dmax, d,
                arb_mat_entry(cho, g - 1, g - 1), &v[g - 1], mid - (k + 1), prec);
            acb_theta_eld_dist_rec(x, cho, next_v, next_dmax, g - 1,
                omit_zero && (mid - (k + 1)) == 0, prec);
            arb_add(x, x, d, prec);
            arb_sub(x, x, next_dmax, prec);
            arb_min(d, d, x, prec);
            arb_nonnegative_part(d, d);

            /* Recompute min, max */
            arb_sqrt(x, d, prec);
            arb_div(x, x, arb_mat_entry(cho, g - 1, g - 1), prec);
            arb_get_ubound_arf(rad, x, prec);
            res = acb_theta_eld_interval(&min, &mid, &max, ctr, rad);
            k++;
        }

        _arb_vec_clear(v_mid, g - 1);
        _arb_vec_clear(v_diff, g - 1);
        _arb_vec_clear(next_v, g - 1);
        arb_clear(next_dmax);
    }

    if (!res)
    {
        /* Could not compute values of n_{g-1} or induction failed, return uniform bound */
        acb_theta_eld_dist_unif(d, cho, omit_zero, prec);
    }

    arb_clear(x);
    arb_clear(ctr);
    arf_clear(rad);
}

static void
acb_theta_eld_dist_ubound(arf_t u, arb_srcptr v, const arb_mat_t cho,
    int omit_zero, slong prec)
{
    slong g = arb_mat_nrows(cho);
    arb_mat_t choinv;
    arb_ptr x;
    slong * pt;
    arb_t d;
    slong k;
    int r = 1;

    arb_mat_init(choinv, g, g);
    x = _arb_vec_init(g);
    pt = flint_malloc(g * sizeof(slong));
    arb_init(d);

    arb_mat_one(choinv);
    arb_mat_solve_triu(choinv, cho, choinv, 0, prec);
    arb_mat_vector_mul_col(x, choinv, v, prec);
    r = _arb_vec_is_finite(x, g);

    for (k = 0; (k < g) && r; k++)
    {
        r = (arf_cmpabs_2exp_si(arb_midref(&x[k]), 30) <= 0);
        if (r)
        {
            pt[k] = - arf_get_si(arb_midref(&x[k]), ARF_RND_NEAR);
        }
    }

    if (r)
    {
        if (omit_zero && acb_theta_eld_pt_is_zero(pt, g))
        {
            pt[0] = 1;
        }
        acb_theta_eld_dist_pt(d, v, cho, pt, prec);
    }
    else
    {
        acb_theta_eld_dist_unif(d, cho, omit_zero, prec);
    }
    arb_get_ubound_arf(u, d, prec);

    arb_mat_clear(choinv);
    _arb_vec_clear(x, g);
    flint_free(pt);
    arb_clear(d);
}

void
acb_theta_eld_dist(arb_t d, arb_srcptr v, const arb_mat_t cho, int omit_zero, slong prec)
{
    slong g = arb_mat_nrows(cho);
    arf_t u;
    arb_t dmax;

    arf_init(u);
    arb_init(dmax);

    acb_theta_eld_dist_ubound(u, v, cho, omit_zero, prec);
    arb_set_arf(dmax, u);
    acb_theta_eld_dist_rec(d, cho, v, dmax, g, omit_zero, prec);
    arb_nonnegative_part(d, d);

    arf_clear(u);
    arb_clear(dmax);
}
