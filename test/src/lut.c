/* Hello emacs, this is -*- c -*- */
/* $Id: lut.c,v 1.1 2001/07/13 15:48:02 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

#include "lut.h"
#include "real.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static real_t* lut_table; /* The lookup table itself */
static real_t  lut_limit; /* The upper value of the table */
static real_t  lut_precision; /* The step precision */
static real_t* lut_greatest; /* The greatest value on the table. Remember the
				lowest is always lut_table!! */

void config_lut(const real_t prec, const real_t lim)
{
  real_t* tmp;
  real_t it; /* iterator */

  int size = (int)ceilf(lim/prec) + 1;

  lut_precision = prec;
  lut_limit = lim;
  lut_table = (real_t*) calloc (size, sizeof(real_t));
  if ( lut_table == NULL ) {
    perror("[lut] LUT space");
    exit(EXIT_FAILURE);
  }

  tmp = lut_table; /* just for manipulation, look: */
  /* This will actually build the lookup table */
  for (it=0; it<lim; it+=prec) *(tmp++) = FUNC_TANH(it);

  /* Now we check, because, sometimes the precision is _not_ enough short to
     fill all SIZE positions */
  if ( lut_table[size-1] < lut_table[size-2] ) {
    lut_greatest = &lut_table[size-2];
    lut_limit -= lut_precision;
  }
  else {
    lut_greatest = &lut_table[size-1];
  }

  fprintf(stderr, "[lut] Single-Sided Lookup Table Configured\n");

  return;
}

void clear_lut(void) 
{
  free(lut_table);
  lut_precision = 0;
}

real_t lookup(const real_t val)
{
  int pos; /* controls the table position */
  real_t absolute; /* a temporary real_t holder */

  /* Most of the time we get extremes, than I have to test for them first. */
  if ( (absolute=FUNC_ABS(val)) >= lut_limit ) {
    if (val > 0) return +1;
    else return -1;
  }

  /* If we're here, I have to calculate the position anyhow */
  pos = (int)FUNC_RINT(absolute/lut_precision);
  if (val < 0) return -1*lut_table[pos];
  
  return lut_table[pos];
}



