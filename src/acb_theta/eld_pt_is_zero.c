/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
acb_theta_eld_pt_is_zero(const slong * pt, slong g)
{
    slong k;

    for (k = 0; k < g; k++)
    {
        if (pt[k] != 0)
        {
            return 0;
        }
    }
    return 1;
}
