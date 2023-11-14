/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_exp(acb_t r, const acb_t z, slong prec)
{
#define a acb_realref(z)
#define b acb_imagref(z)
    if (arb_is_zero(b))
    {
        arb_exp(acb_realref(r), a, prec);
        arb_zero(acb_imagref(r));
    }
    else if (arb_is_zero(a))
    {
        arb_sin_cos(acb_imagref(r), acb_realref(r), b, prec);
    }
    else
    {
        arb_t t, u, v;

        arb_init(t);
        arb_init(u);
        arb_init(v);

        arb_exp(t, a, prec);
        arb_sin_cos(u, v, b, prec);

        arb_mul(acb_realref(r), t, v, prec);
        arb_mul(acb_imagref(r), t, u, prec);

        arb_clear(t);
        arb_clear(u);
        arb_clear(v);
    }
#undef a
#undef b
}

