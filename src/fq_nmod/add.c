/*
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "fq_nmod.h"

void fq_nmod_add(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
{
    slong max = FLINT_MAX(op1->length, op2->length);

    nmod_poly_fit_length(rop, max);

    _nmod_poly_add(rop->coeffs,
                   op1->coeffs, op1->length, op2->coeffs, op2->length,
                   rop->mod);

    _nmod_poly_set_length(rop, max);
    _nmod_poly_normalise(rop);
}
