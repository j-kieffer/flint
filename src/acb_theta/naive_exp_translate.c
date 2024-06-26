/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_naive_exp_translate(acb_t c, acb_ptr exp_z, const acb_mat_t exp_tau_div_4,
	ulong a, slong prec)
{
	slong g = acb_mat_nrows(exp_tau_div_4);
	acb_t s;
	slong j, k;

	acb_init(s);

	acb_one(c);
	for (j = 0; j < g; j++)
	{
		/* do nothing if aj = 0 */
		if (!((a >> j) & 1))
		{
			continue;
		}

		/* cofactor */
		for (k = j; k < g; k++)
		{
			if ((a >> k) & 1)
			{
				acb_mul(c, c, acb_mat_entry(exp_tau_div_4, j, k), prec);
			}
		}

		/* translation of z */
		acb_one(&exp_z[j]);
		for (k = 0; k < g; k++)
		{
			if (!((a >> k) & 1))
			{
				continue;
			}

			if (k < j)
			{
				acb_set(s, acb_mat_entry(exp_tau_div_4, k, j));
			}
			else if (k == j)
			{
				acb_sqr(s, acb_mat_entry(exp_tau_div_4, k, k), prec);
			}
			else
			{
				acb_set(s, acb_mat_entry(exp_tau_div_4, j, k));
			}
			acb_mul(&exp_z[j], &exp_z[j], s, prec);
		}
	}
	acb_inv(c, c, prec);

	acb_clear(s);
}
