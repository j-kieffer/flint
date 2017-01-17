/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qsieve.h"

void qsieve_linalg_init(qs_t qs_inf)
{
    slong i, num_primes;

    qs_inf->extra_rels = 64; /* number of opportunities to factor n */
    qs_inf->max_factors = 60; /* maximum number of factors a (merged) relation can have */

    /* allow as many dups as relations */
    num_primes = qs_inf->num_primes;
    qs_inf->num_primes += qs_inf->ks_primes;

    qs_inf->buffer_size = 2*(qs_inf->num_primes + qs_inf->extra_rels);
    qs_inf->small = flint_malloc(qs_inf->small_primes*sizeof(mp_limb_t));
    qs_inf->factor = flint_malloc(qs_inf->max_factors*sizeof(fac_t));
    qs_inf->matrix = flint_malloc((qs_inf->buffer_size)*sizeof(la_col_t));
    qs_inf->unmerged = qs_inf->matrix + qs_inf->buffer_size;
    qs_inf->Y_arr = flint_malloc(qs_inf->buffer_size*sizeof(fmpz));
    qs_inf->curr_rel = qs_inf->relation
                     = flint_malloc(2*qs_inf->buffer_size*qs_inf->max_factors*sizeof(slong));

    for (i = 0; i < qs_inf->buffer_size; i++)
    {
        fmpz_init(qs_inf->Y_arr + i);
        qs_inf->matrix[i].weight = 0;
        qs_inf->matrix[i].data = NULL;
    }

    qs_inf->prime_count = flint_malloc(qs_inf->num_primes*sizeof(slong));

    qs_inf->num_primes = num_primes;
    qs_inf->num_unmerged = 0;
    qs_inf->columns = 0;
    qs_inf->num_relations = 0;

    /* parameter related to partials */
    qs_inf->full_relation = 0;
    qs_inf->edges = 0;
    qs_inf->vertices = 0;
    qs_inf->components = 1;
    qs_inf->num_cycles = 0;

    qs_inf->table_size = 10000;
    qs_inf->hash_table = flint_calloc(1 << 25, sizeof(mp_limb_t));
    qs_inf->table = flint_malloc(qs_inf->table_size * sizeof(hash_t));
}

/* re-initialize all the linear algebra parameter */

void qsieve_linalg_re_init(qs_t qs_inf)
{
    slong i;
    qs_inf->curr_rel = qs_inf->relation;

    for (i = 0; i < qs_inf->buffer_size; i++)
    {
        fmpz_init(qs_inf->Y_arr + i);
        qs_inf->matrix[i].weight = 0;
        qs_inf->matrix[i].data = NULL;
    }

    qs_inf->num_unmerged = 0;
    qs_inf->columns = 0;
    qs_inf->num_relations = 0;

    /* parameters related to partials */
    qs_inf->full_relation = 0;
    qs_inf->edges = 0;
    qs_inf->vertices = 0;
    qs_inf->components = 1;
    qs_inf->num_cycles = 0;

    memset(qs_inf->hash_table, 0, (1 << 25) * sizeof(mp_limb_t));
}

/* increase size of different array after factor base increment*/

void qsieve_linalg_re_alloc(qs_t qs_inf)
{
    slong num_primes;

    num_primes = qs_inf->num_primes;
    qs_inf->num_primes += qs_inf->ks_primes;
    qs_inf->buffer_size = 2*(qs_inf->num_primes + qs_inf->extra_rels);
    qs_inf->matrix = flint_realloc(qs_inf->matrix, qs_inf->buffer_size*sizeof(la_col_t));
    qs_inf->unmerged = qs_inf->matrix + qs_inf->buffer_size;
    qs_inf->Y_arr = flint_realloc(qs_inf->Y_arr, qs_inf->buffer_size*sizeof(fmpz));
    qs_inf->curr_rel = qs_inf->relation
                     = flint_realloc(qs_inf->relation, 2*qs_inf->buffer_size*qs_inf->max_factors*sizeof(slong));

    qs_inf->prime_count = flint_realloc(qs_inf->prime_count, qs_inf->num_primes*sizeof(slong));
    qs_inf->num_primes = num_primes;
    qsieve_linalg_init(qs_inf);
}