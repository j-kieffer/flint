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

/* Returns an interval d that contains Dist(v, cho Z^g)^2 for every vector v */
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

/* Returns an upper bound u on Dist(v, cho Z^g)^2 */
static void
acb_theta_eld_dist_ubound(arf_t dmax, arb_srcptr v, const arb_mat_t cho,
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
    arb_get_ubound_arf(dmax, d, prec);

    arb_mat_clear(choinv);
    _arb_vec_clear(x, g);
    flint_free(pt);
    arb_clear(d);
}

/* Main recursive function */
static void
acb_theta_eld_dist_rec(arb_t d, const arb_mat_t cho, arb_srcptr v,
    const arf_t dmax, slong g, int omit_zero, slong prec);

/* Helper function for eld_dist */
static int
acb_theta_eld_dist_update(slong * min, slong * mid, slong * max, arb_t d,
    arf_t d_up, arf_t d_low, const arb_mat_t cho, arb_srcptr v, slong n,
    slong g, int omit_zero, slong prec)
{
    arb_ptr next_v;
    arb_t next_dmax, x;
    arf_t t;
    int update;
    int res = 1;

    next_v = _arb_vec_init(g - 1);
    arb_init(next_dmax);
    arb_init(x);
    arf_init(t);

    acb_theta_eld_recursion(next_dmax, next_v, cho, d, v, n, g, prec);
    arb_get_ubound_arf(t, next_dmax, prec);
    acb_theta_eld_dist_rec(x, cho, next_v, t, g - 1, omit_zero && (n == 0), prec);
    arb_sub(x, x, next_dmax, prec);
    arb_add(x, x, d, prec);

    arb_get_lbound_arf(t, x, prec);
    if (arf_cmp(t, d_low) < 0)
    {
        update = 1;
        if (arf_cmp_si(t, 0) <= 0)
        {
            arf_zero(d_low);
        }
        else
        {
            arf_set(d_low, t);
        }
    }

    arb_get_ubound_arf(t, x, prec);
    if (arf_cmp(t, d_up) < 0)
    {
        update = 1;
        arf_set(d_up, t);

        /* Recompute min, max */
        res = acb_theta_eld_interval(min, mid, max,
            arb_mat_entry(cho, g - 1, g - 1), d_up, &v[g - 1], prec);
    }

    if (update)
    {
        arb_set_arf(d, d_low);
        arb_set_arf(x, d_up);
        arb_union(d, d, x, prec);
    }

    _arb_vec_clear(next_v, g - 1);
    arb_clear(next_dmax);
    arb_clear(x);
    arf_clear(t);
    return res;
}

static void
acb_theta_eld_dist_rec(arb_t d, const arb_mat_t cho, arb_srcptr v,
    const arf_t dmax, slong g, int omit_zero, slong prec)
{
    slong min, mid, max;
    int res = 1;

    /* Get possible interval for n_{g-1} */
    res = acb_theta_eld_interval(&min, &mid, &max,
        arb_mat_entry(cho, g - 1, g - 1), dmax, &v[g - 1], prec);

    if (res && min > max)
    {
        /* Distance is bigger than dmax */
        arb_set_arf(d, dmax);
        return;
    }
    else if (res && g == 1)
    {
        /* Can compute distance directly */
        arb_t x;

        arb_init(x);

        arb_set_si(x, (omit_zero && mid == 0 ? mid + 1 : mid));
        arb_mul(x, x, arb_mat_entry(cho, g - 1, g - 1), prec);
        arb_add(x, x, &v[g - 1], prec);
        arb_sqr(d, x, prec);

        if (!omit_zero || mid + 1 != 0)
        {
            arb_mul_si(x, arb_mat_entry(cho, g - 1, g - 1), mid + 1, prec);
            arb_add(x, x, &v[g - 1], prec);
            arb_sqr(x, x, prec);
            arb_min(d, d, x, prec);
        }
        if (!omit_zero || mid - 1 != 0)
        {
            arb_mul_si(x, arb_mat_entry(cho, g - 1, g - 1), mid - 1, prec);
            arb_add(x, x, &v[g - 1], prec);
            arb_sqr(x, x, prec);
            arb_min(d, d, x, prec);
        }

        arb_clear(x);
    }
    else if (res)
    {
        /* Induction */
        arf_t d_low, d_up; // manipulating arf's reduces precision losses
        slong k;

        arf_init(d_low);
        arf_init(d_up);

        arb_set_arf(d, dmax);
        arf_set(d_low, dmax);
        arf_set(d_up, dmax);

        /* Right loop */
        k = 0;
        while (res && mid + k < max + 1)
        {
            acb_theta_eld_dist_update(&min, &mid, &max, d, d_up, d_low,
                cho, v, mid + k, g, omit_zero, prec);
            k++;
        }

        /* Left loop */
        k = 0;
        while (res && mid - (k + 1) > min - 1)
        {
            acb_theta_eld_dist_update(&min, &mid, &max, d, d_up, d_low,
                cho, v, mid - (k + 1), g, omit_zero, prec);
            k++;
        }

        arf_clear(d_low);
        arf_clear(d_up);
    }

    if (!res)
    {
        /* Could not compute values of n_{g-1} or induction failed, return uniform bound */
        acb_theta_eld_dist_unif(d, cho, omit_zero, prec);
    }
}

void
acb_theta_eld_dist(arb_t d, arb_srcptr v, const arb_mat_t cho, int omit_zero, slong prec)
{
    slong g = arb_mat_nrows(cho);
    arf_t dmax;

    arf_init(dmax);

    acb_theta_eld_dist_ubound(dmax, v, cho, omit_zero, prec);
    acb_theta_eld_dist_rec(d, cho, v, dmax, g, omit_zero, prec);

    arf_clear(dmax);
}
