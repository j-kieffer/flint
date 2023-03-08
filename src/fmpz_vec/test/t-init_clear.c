/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    
    flint_printf("init/clear....");
    fflush(stdout);
    
    

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a;
        slong j, len = n_randint(state, 100) + 1;

        a = _fmpz_vec_init(len);
        for (j = 0; j < len; j++)
            fmpz_zero(a + j);

        _fmpz_vec_clear(a, len);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
