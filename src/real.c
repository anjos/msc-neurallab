/* Hello emacs, this is -*- c -*- */
/* $Id: real.c,v 1.1 2001/07/13 15:48:02 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

#include "real.h"

/* These implementations are usually slower than the ones from LIBC, but you
   can use them just to profile your applications. */

/* Returns the modules of the real value given by the first argument */
real_t modulus (const real_t arg) 
{
  if ( arg < 0 ) return -arg;
  return arg;
}

/* Returns the closest integer to the given argument */
int realint (const real_t arg)
{
  int closest = (int)arg;
  if ( arg-closest > (closest+1)-arg ) return closest+1;
  return closest;
}
