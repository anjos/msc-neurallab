/* Hello emacs, this is -*- c -*- */
/* $Id: real.h,v 1.1 2001/07/13 15:48:02 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

/* This module defines types and functions common to all other modules. */

#ifndef REAL_H_
#define REAL_H_

/* Choose 64 or 32-bit precision. All related function calls will auto-adjust
   for that choice. */
#define _PRECISION_32BIT

/* Do _not_ fiddle from here on */
#ifdef _PRECISION_32BIT
typedef float real_t; /* This will define what's a real to us */
#define FUNC_TANH tanhf /* This will define what is a tanh for us */
#define FUNC_ABS fabsf /* This will define what is a abs() for us */
#define FUNC_POW powf /* how to apply a power */
#define FUNC_RINT rintf /* This is how numbers are rounded */
#define FUNC_FLOOR floorf
#define FUNC_CEIL ceilf
#define FSCANF_TYPE "%e" /* How to read reals */
#define PREC_POPT_REAL POPT_ARG_STRING /* The input precision description */
#endif /* _PRECISION_32BIT */

#ifdef _PRECISION_64BIT
typedef double real_t; /* This will define what's a real to us */
#define FUNC_TANH tanh /* This will define what is a tanh for us */
#define FUNC_ABS fabs /* This will define what is a abs() for us */
#define FUNC_POW pow /* how to apply a power */
#define FUNC_RINT rint /* This is how numbers are rounded */
#define FUNC_FLOOR floor
#define FUNC_CEIL ceil
#define FSCANF_TYPE "%le" /* How to read reals */
#define PREC_POPT_REAL POPT_ARG_STRING /* The input precision description */
#endif /* _PRECISION_64BIT */

/* What is a boolean for us */
typedef enum bool_t {FALSE = 0, TRUE = 1} bool_t;

/* This functions define some basic stuff we have to use */

/* Returns the modules of the real value given by the first argument */
real_t modulus (const real_t);

/* Returns the closest integer to the given argument. Watch-out! This function
   does _not_ have provisions to deal with int overflow! */
int realint (const real_t);

#endif /* REAL_H_ */
