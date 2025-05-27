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

static void
acb_theta_slong_vec_max(slong * r, const slong * v1, const slong * v2, slong d)
{
    slong k;
    for (k = 0; k < d; k++)
    {
        r[k] = FLINT_MAX(v1[k], v2[k]);
    }
}

/* Main recursive function in dimension d > 1 */
static int
acb_theta_eld_set_rec(acb_theta_eld_t E, const arb_mat_t cho,
    const arf_t R2, arb_srcptr v, slong * coords, slong prec);

/* Helper function for eld_set_rec */
static int
acb_theta_eld_set_child(acb_theta_eld_t E, acb_theta_eld_t child,
    const arb_mat_t cho, const arf_t R2, arb_srcptr v,
    slong n, slong * coords, slong d, slong prec)
{
    arb_ptr next_v;
    arb_t x;
    arf_t next_R2;
    int res;

    next_v = _arb_vec_init(d - 1);
    arb_init(x);
    arf_init(next_R2);

    arb_set_arf(x, R2);
    acb_theta_eld_recursion(x, next_v, cho, x, v, n, d, prec);
    arb_get_ubound_arf(next_R2, x, prec);

    coords[d - 1] = n;
    res = acb_theta_eld_set_rec(child, cho, next_R2, next_v, coords, prec);
    if (res)
    {
        E->nb_pts += child->nb_pts;
        E->nb_border += child->nb_border;
        acb_theta_slong_vec_max(E->box, E->box, child->box, d - 1);
    }

    _arb_vec_clear(next_v, d - 1);
    arb_clear(x);
    arf_clear(next_R2);
    return res;
}

static int
acb_theta_eld_set_rec(acb_theta_eld_t E, const arb_mat_t cho,
    const arf_t R2, arb_srcptr v, slong * coords, slong prec)
{
    slong d = E->dim;
    slong g = E->ambient_dim;
    slong min, mid, max, nr, nl;
    slong k;
    int res;

    res = acb_theta_eld_interval(&min, &mid, &max,
        arb_mat_entry(cho, d - 1, d - 1), R2, &v[d - 1], prec);
    if (!res)
    {
        return 0;
    }

    for (k = 0; k < g - d; k++)
    {
        E->last_coords[k] = coords[d + k];
    }
    E->min = min;
    E->mid = mid;
    E->max = max;
    E->nb_pts = 0;
    E->nb_border = (d == 1 ? 2 : 0);
    E->box[d - 1] = FLINT_MAX(max, -min);
    for (k = 0; k < d - 1; k++)
    {
        E->box[k] = 0;
    }

    /* Induction only if d > 1 and min <= max */
    if (min > max)
    {
        E->box[d - 1] = 0;
        return 1;
    }
    else if (d == 1)
    {
        E->nb_pts = max - min + 1;
        return 1;
    }

    nr = max - mid + 1;
    nl = mid - min;

    /* Initialize children */
    if (nr > 0) /* should always be the case */
    {
        E->rchildren = flint_malloc(nr * sizeof(struct acb_theta_eld_struct));
        E->nr = nr;
        for (k = 0; k < nr; k++)
        {
            acb_theta_eld_init(&E->rchildren[k], d - 1, g);
        }
    }
    if (nl > 0)
    {
        E->lchildren = flint_malloc(nl * sizeof(struct acb_theta_eld_struct));
        E->nl = nl;
        for (k = 0; k < nl; k++)
        {
            acb_theta_eld_init(&E->lchildren[k], d - 1, g);
        }
    }

    /* Loop over children */
    for (k = 0; (k < nr) && res; k++)
    {
        res = acb_theta_eld_set_child(E, &E->rchildren[k], cho, R2, v,
            mid + k, coords, d, prec);
    }
    for (k = 0; (k < nl) && res; k++)
    {
        res = acb_theta_eld_set_child(E, &E->lchildren[k], cho, R2, v,
            mid - (k + 1), coords, d, prec);
    }

    return res;
}

int
acb_theta_eld_set(acb_theta_eld_t E, const arb_mat_t cho, const arf_t R2,
    arb_srcptr v, slong prec)
{
    slong d = E->dim;
    slong g = E->ambient_dim;
    slong * coords;
    int res;

    coords = flint_malloc(g * sizeof(slong));
    acb_theta_eld_clear(E);
    acb_theta_eld_init(E, d, g);

    res = acb_theta_eld_set_rec(E, cho, R2, v, coords, prec);

    flint_free(coords);
    return res;
}
