/* Hello emacs, this is -*- c -*- */
/* $Id: lut.h,v 1.1 2001/07/13 15:48:02 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

/* This module controls some aspects of lookup table processing. It's mainly
   use nowadays is for neural networking. */

#ifndef LUT_H_
#define LUT_H_

#include "real.h"

/* This function will configure a single side lookup table function. If it was
   already configured in a previous step, use clear_lut() before re-configuring
   the table. This will free any required resources.  The first argument is the
   precision of the table and the second represents the upper limit. So,
   spacing occurs from 0 to limit.  During LUT configuration, the size
   accouting goes like this: if the limit/precision is an integer, the table
   is exactly calculated, but if limit/precision is _not_ an integer, than,
   the table size is rounded upwards to the next integer using ceilf. The
   number accounting, in that case may not be the desired result. */
void config_lut(const real_t, const real_t);

/* This function clears the lookup table, freeing all required resources */
void clear_lut(void);

/* This function will perform the lookup of a real_t value */
real_t lookup(const real_t);

#endif /* LUT_H_ */
