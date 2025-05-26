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
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_eld_shortest(arb_t rho, const arb_mat_t cho, slong prec)
{
    slong g = arb_mat_nrows(cho);
    arb_ptr zero;

    zero = _arb_vec_init(g);

    acb_theta_eld_dist(rho, zero, cho, 1, prec);
    arb_sqrt(rho, rho, prec);

    _arb_vec_clear(zero, g);
}

