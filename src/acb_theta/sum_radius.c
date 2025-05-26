/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "arb.h"
#include "arb_mat.h"
#include "arb_hypgeom.h"
#include "acb_theta.h"

/* Assuming a >= 0, return R2 such that x - (a/2)*log(x)\geq b for all
   x\geq R2, and R2 is close to the smallest possible */

static void
acb_theta_invert_lin_plus_log(arf_t R2, slong a, const arb_t b, slong prec)
{
    arb_t x, y, t;
    arf_t z;
    slong k;

    arb_init(x);
    arb_init(y);
    arb_init(t);
    arf_init(z);

    if (a == 0)
    {
        arb_get_ubound_arf(R2, b, prec);
        goto exit;
    }

    /* minimum is at x=a/2 */
    arb_set_si(x, a);
    arb_div_si(x, x, 2, prec);
    arb_log(y, x, prec);
    arb_mul(y, y, x, prec);
    arb_sub(y, x, y, prec);

    /* x = max(a, 2*(b - min) + a log 2) is always large enough; then iterate
       function a few times */
    arb_sub(y, b, y, prec);
    arb_const_log2(t, prec);
    arb_mul_2exp_si(t, t, -1);
    arb_mul_si(t, t, a, prec);
    arb_add(y, y, t, prec);
    arb_max(y, y, x, prec);
    arb_mul_si(x, y, 2, prec);
    arb_get_ubound_arf(z, x, prec);
    arb_set_arf(x, z);

    for (k = 0; k < 4; k++)
    {
        arb_log(y, x, prec);
        arb_mul_si(y, y, a, prec);
        arb_div_si(y, y, 2, prec);
        arb_add(x, b, y, prec);
        arb_get_ubound_arf(z, x, prec);
        arb_set_arf(x, z);
    }

    arb_get_ubound_arf(R2, x, prec);
    goto exit;

  exit:
    {
        arb_clear(x);
        arb_clear(y);
        arb_clear(t);
        arf_clear(z);
    }
}

static void
acb_theta_sum_radius_v1(arf_t R2, const arb_mat_t cho, slong ord, slong prec)
{
    slong g = arb_mat_nrows(cho);
    slong lp = ACB_THETA_LOW_PREC;
    arb_t b, temp, sqrt2pi;
    arf_t cmp;
    slong k;

    arb_init(b);
    arb_init(temp);
    arb_init(sqrt2pi);
    arf_init(cmp);

    arb_const_pi(sqrt2pi, lp);
    arb_mul_2exp_si(sqrt2pi, sqrt2pi, 1);
    arb_sqrt(sqrt2pi, sqrt2pi, lp);

    /* Set b such that
       (1 + 8/sqrt(pi)) * prod_j (1 + sqrt(2pi)/c_j) * b \leq 2^(-prec) */
    arb_const_sqrt_pi(temp, lp);
    arb_inv(temp, temp, lp);
    arb_mul_2exp_si(temp, temp, 3);
    arb_add_si(b, temp, 1, lp);
    for (k = 0; k < g; k++)
    {
        arb_div(temp, sqrt2pi, arb_mat_entry(cho, k, k), lp);
        arb_add_si(temp, temp, 1, lp);
        arb_mul(b, b, temp, lp);
    }
    arb_inv(b, b, lp);
    arb_mul_2exp_si(b, b, -prec);

    /* Solve R2^((g-1)/2+ord) exp(-R2) \leq b */
    arb_log(b, b, lp);
    arb_neg(b, b);
    acb_theta_invert_lin_plus_log(R2, g - 1 + 2 * ord, b, lp);

    /* Max with 4, 2*ord for formula to be valid */
    arf_set_si(cmp, FLINT_MAX(4, 2 * ord));
    arf_max(R2, R2, cmp);

    arb_clear(b);
    arb_clear(temp);
    arf_clear(cmp);
}

static void
acb_theta_eld_shortest(arb_t rho, const arb_mat_t cho, slong prec)
{
    slong g = arb_mat_nrows(cho);
    arb_ptr zero;

    zero = _arb_vec_init(g);

    acb_theta_eld_dist(rho, zero, cho, 1, prec);

    _arb_vec_clear(zero, g);
}

static void
acb_theta_radius_rhs_v2(arb_t res, const arb_t rho, const arf_t R2, slong ord,
    slong g, slong prec)
{
    arb_t x, t;

    arb_init(x);
    arb_init(t);

    /* Set x = (R - rho/2)^2 */
    arb_set_arf(x, R2);
    arb_sqrt(x, x, prec);
    arb_mul_2exp_si(t, rho, -1);
    arb_sub(x, x, t, prec);
    arb_sqr(x, x, prec);

    /* Set x to incomplete Gamma */
    arb_set_si(t, g + ord);
    arb_mul_2exp_si(t, t, -1);
    arb_hypgeom_gamma_upper(x, t, x, 0, prec);

    /* Multiply by the right factors */
    arb_set_si(t, 2);
    arb_div(t, t, rho, prec);
    arb_pow_ui(t, t, g, prec);
    arb_mul_si(t, t, g, prec);
    arb_mul_2exp_si(t, t, -1);
    arb_mul(x, x, t, prec);

    arb_set(res, x);
    arb_clear(x);
    arb_clear(t);
}

static void
acb_theta_sum_radius_v2(arf_t R2, const arf_t R2max, const arb_mat_t cho,
    slong ord, slong prec)
{
    slong g = arb_mat_nrows(cho);
    slong lp = ACB_THETA_LOW_PREC;
    arb_t rho;
    arb_t u;
    arf_t newR2, R2min;

    arb_init(rho);
    arb_init(u);
    arf_init(newR2);
    arf_init(R2min);

    acb_theta_eld_shortest(rho, cho, lp);
    arf_set(R2, R2max);
    arf_set(newR2, R2);

    /* Set R2min */
    arb_set_si(u, g * g + 8 * ord);
    arb_sqrt(u, u, lp);
    arb_add_si(u, u, g + 2 * ord, lp);
    arb_sqrt(u, u, lp);
    arb_add(u, u, rho, lp);
    arb_mul_2exp_si(u, u, -1);
    arb_sqr(u, u, lp);
    arb_get_ubound_arf(R2min, u, lp);

    while (1)
    {
        /* Let newR2 = R2 * (1 - 1/10g) */
        arb_set_si(u, - 10 * g);
        arb_inv(u, u, lp);
        arb_add_si(u, u, 1, lp);
        arb_mul_arf(u, u, R2, lp);
        arb_get_ubound_arf(newR2, u, lp);
        if (arf_cmp(newR2, R2min) < 0)
        {
            arf_set(newR2, R2min);
        }
        if (arf_cmp(newR2, R2) >= 0)
        {
            break;
        }

        /* Check if newR2 still works, otherwise stop */
        acb_theta_radius_rhs_v2(u, rho, newR2, ord, g, lp);
        arb_mul_2exp_si(u, u, prec);
        arb_sub_si(u, u, 1, lp);
        if (arb_is_positive(u))
        {
            break;
        }
        else
        {
            arf_set(R2, newR2);
        }
    }

    arb_clear(rho);
    arb_clear(u);
    arf_clear(newR2);
    arf_clear(R2min);
}

void
acb_theta_sum_radius(arf_t R2, arf_t eps, const arb_mat_t cho, slong ord, slong prec)
{
    arf_t u;
    arf_init(u);

    acb_theta_sum_radius_v1(R2, cho, ord, prec);
    acb_theta_sum_radius_v2(u, R2, cho, ord, prec);

    arf_min(R2, R2, u);

    /* Set error 2^(-prec) */
    arf_one(eps);
    arf_mul_2exp_si(eps, eps, -prec);

    arf_clear(u);
}
