/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq.h"
#include "fmpq_poly.h"

void fmpq_poly_set_coeff_fmpq(fmpq_poly_t poly, slong n, const fmpq_t x)
{
    slong len = poly->length;
    const int replace = (n < len && !fmpz_is_zero(poly->coeffs + n));

    if (!replace && fmpq_is_zero(x))
        return;

    if (n + 1 > len)
    {
        fmpq_poly_fit_length(poly, n + 1);
        _fmpq_poly_set_length(poly, n + 1);
        flint_mpn_zero((mp_ptr) poly->coeffs + len, (n + 1) - len);
        len = n + 1;
    }

    if (replace)
    {
        fmpz_t c;

        fmpz_init(c);

        fmpz_zero(poly->coeffs + n);
        _fmpz_poly_content(c, poly->coeffs, len);
        _fmpz_vec_scalar_mul_fmpz(poly->coeffs, poly->coeffs, len, fmpq_denref(x));
        fmpz_mul(c, c, fmpq_denref(x));

        fmpz_mul(poly->coeffs + n, fmpq_numref(x), poly->den);
        fmpz_gcd(c, c, poly->coeffs + n);
        fmpz_mul(poly->den, poly->den, fmpq_denref(x));

        if (!fmpz_is_one(c))
            fmpz_gcd(c, c, poly->den);
        if (!fmpz_is_one(c))
        {
            _fmpz_vec_scalar_divexact_fmpz(poly->coeffs, poly->coeffs, len, c);
            fmpz_divexact(poly->den, poly->den, c);
        }

        _fmpq_poly_normalise(poly);
        fmpz_clear(c);
    }
    else
    {
        fmpz_t d, t;

        fmpz_init(d);
        fmpz_init(t);

        fmpz_gcd(d, poly->den, fmpq_denref(x));
        fmpz_divexact(t, fmpq_denref(x), d);

        _fmpz_vec_scalar_mul_fmpz(poly->coeffs, poly->coeffs, len, t);
        fmpz_set(poly->coeffs + n, fmpq_numref(x));
        fmpz_mul(poly->coeffs + n, poly->coeffs + n, poly->den);
        fmpz_divexact(poly->coeffs + n, poly->coeffs + n, d);

        fmpz_mul(poly->den, poly->den, t);

        fmpz_clear(d);
        fmpz_clear(t);
    }
}

void fmpq_poly_set_coeff_fmpz(fmpq_poly_t poly, slong n, const fmpz_t x)
{
    slong len = poly->length;
    const int replace = (n < len && !fmpz_is_zero(poly->coeffs + n));

    if (!replace && fmpz_is_zero(x))
        return;

    if (n + 1 > len)
    {
        fmpq_poly_fit_length(poly, n + 1);
        _fmpq_poly_set_length(poly, n + 1);
        flint_mpn_zero((mp_ptr) poly->coeffs + len, (n + 1) - len);
    }

    if (*poly->den == WORD(1))
    {
        fmpz_set(poly->coeffs + n, x);
        if (replace)
            _fmpq_poly_normalise(poly);
    }
    else
    {
        fmpz_mul(poly->coeffs + n, poly->den, x);
        if (replace)
            fmpq_poly_canonicalise(poly);
    }
}

void fmpq_poly_set_coeff_si(fmpq_poly_t poly, slong n, slong x)
{
    slong len = poly->length;
    const int replace = (n < len && !fmpz_is_zero(poly->coeffs + n));

    if (!replace && (x == WORD(0)))
        return;

    if (n + 1 > len)
    {
        fmpq_poly_fit_length(poly, n + 1);
        _fmpq_poly_set_length(poly, n + 1);
        flint_mpn_zero((mp_ptr) poly->coeffs + len, (n + 1) - len);
    }

    if (*poly->den == WORD(1))
    {
        fmpz_set_si(poly->coeffs + n, x);
        if (replace)
            _fmpq_poly_normalise(poly);
    }
    else
    {
        fmpz_mul_si(poly->coeffs + n, poly->den, x);
        if (replace)
            fmpq_poly_canonicalise(poly);
    }
}

void fmpq_poly_set_coeff_ui(fmpq_poly_t poly, slong n, ulong x)
{
    slong len = poly->length;
    const int replace = (n < len && !fmpz_is_zero(poly->coeffs + n));

    if (!replace && (x == UWORD(0)))
        return;

    if (n + 1 > len)
    {
        fmpq_poly_fit_length(poly, n + 1);
        _fmpq_poly_set_length(poly, n + 1);
        flint_mpn_zero((mp_ptr) poly->coeffs + len, (n + 1) - len);
    }

    if (*poly->den == WORD(1))
    {
        fmpz_set_ui(poly->coeffs + n, x);
        if (replace)
            _fmpq_poly_normalise(poly);
    }
    else
    {
        fmpz_mul_ui(poly->coeffs + n, poly->den, x);
        if (replace)
            fmpq_poly_canonicalise(poly);
    }
}
