/* Hello emacs, this is -*- c -*- */
/* $Id: analysis.c,v 1.1 2001/07/13 15:48:02 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "real.h"
#include "analysis.h"

/* This function takes the outputs from the neural network at the first
   argument and a threshold at the second argument and finds the efficiencies
   for class1xclass2 discrimination. From those it calculates the SP value and
   returns it. We expect that class 1 has a target _greater_ than class2, WATCH
   OUT FOR THAT!.
   This function will update the last 2 variables with the values of
   efficiency. */
real_t get_sp (const outputdata_t*, const real_t, real_t*, real_t*);

/* Calculate the efficiency for both class 1 and 2 given a threshold */
void effic (const outputdata_t*, const real_t, real_t*, real_t*);

/* This function takes the outputs from the neural network at the first
   argument computes the integral of the curve "eff1 (y) versus
   `3rd.arg'*(1-eff2) (x)". This integral is subtracted from the integral of
   y=1 between 0 and `3rd.arg*(1-eff2[last])' (the last efficiency value for
   eff2). This is what's output. The integration method uses the Simpson's 1/3
   rule (Newton-Cotes, with 2nd. order approximation).  This criterium should
   be a better one than the sp value, since it computes the integral over the
   whole discrimination area instead of a single point. Of course it can't be
   negative.  */
void eval_area (const outputdata_t*, analysis_t*);

/* This function will check on an integer. If it's even, the output is true,
   casen not, the output is false (0). */
int iseven (const int);

/* Now we compute the area under that curve using the Simpon's 1/3 Rule
   Repeated (of course) over all intervals. The formula is:
   I = (h/3) * { f(x0) + f(xm) + 4*[f(x1) + f(x3) + ... + f(xm-1)] + 2*[f(x2) +
   + f(x4) + ... + f(xm-2)] }.
   I've got to find out the even and odd points. The value of m is
   ANALYSIS_NO_THRESHOLDS, if ANALYSIS_NO_THRESHOLDS is even or
   ANALYSIS_NO_THRESHOLDS-1 in the case it's odd, since this integration method
   requires m to be even. The parameters are: y values of the curve, the number
   of y values (checked for evenness), the step h used to extract those
   values. */ 
real_t int_simpson (const real_t*, const int, const real_t);

void analyze (analysis_t* ap, const outputdata_t* out,
	      const analysis_parameter_t* params)
{
  register int i; /* iterator */
  real_t step;
  real_t csp; /* the current SP value just calculated */
  real_t sse1, sse2; /* holds the SSE of each class */

  /* WE EXPECT THAT CLASS 1 HAS A TARGET VALUE GREATER THAN CLASS 2 */
  step = (out->target1 - out->target2) / ANALYSIS_NO_THRESHOLDS; 

  if (step <= 0) {
    fprintf(stderr, "[analyze] target 1 is less than target 2\n");
    fprintf(stderr, "[analyze] target 1 is ");
    fprintf(stderr, "%e", out->target1);
    fprintf(stderr, "\n[analyze] target 2 is ");
    fprintf(stderr, "%e", out->target2);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
  }

  /* Execute all routines */
  for (i=0; i<=ANALYSIS_NO_THRESHOLDS; ++i) {
    ap->threshold[i] = out->target2 + ( i * step );
    csp = get_sp(out, ap->threshold[i], &ap->eff1[i], &ap->eff2[i]);
    if (csp > ap->bestsp.value ) { /* got better */
      ap->bestsp.value = csp;
      ap->bestsp.index = i;
    }
  }

  /* Now evaluate the area under 1 and the efficiency curve */
  ap->area = 0;
  if (params->eval_area) eval_area (out, ap);

  /* Now evaluate the SSE and MSE */
  sse1 = 0;
  sse2 = 0;
  for (i=0; i<out->n1; ++i) sse1 += FUNC_POW(out->target1 - out->class1[i], 2);
  for (i=0; i<out->n2; ++i) sse2 += FUNC_POW(out->target2 - out->class2[i], 2);
  ap->sse = sse1 + sse2;
  ap->mse = ((sse1 / out->n1) + (sse2 / out->n2)) / 2;

  return;
}

/* This function takes the outputs from the neural network at the first
   argument and a threshold at the second argument and finds the efficiencies
   for class1xclass2 discrimination. From those it calculates the SP value and
   returns it. We expect that class 1 has a target _greater_ than class2, WATCH
   OUT FOR THAT!*/
real_t get_sp (const outputdata_t* out, const real_t t, real_t* eff1,
	       real_t* eff2)
{
  /* calculates the efficiency for both classes */
  effic(out, t, eff1, eff2); 

  /* return the SP product */
  return ((*eff1) + (*eff2)) * ((*eff1)*(*eff2));
}

/* Calculate the efficiency for both class 1 and 2 given a threshold */
void effic (const outputdata_t* out, const real_t t, real_t* eff1, real_t* eff2)
{
  /* an iterator */
  register int i;

  /* some counters */
  int right1 = 0;
  int right2 = 0;
  
  /* Count the right ones we got */
  for (i=0; i<out->n1; ++i) if (out->class1[i] >= t) ++right1;
  for (i=0; i<out->n2; ++i) if (out->class2[i] <= t) ++right2;

  (*eff1) = ((real_t)right1/(real_t)out->n1);
  (*eff2) = ((real_t)right2/(real_t)out->n2);
  
  return;
}

/* Copies src int dest */
void copy_analysis (const analysis_t* src, analysis_t* dest)
{
  register int i; /* iterator */

  for (i=0; i<=ANALYSIS_NO_THRESHOLDS; ++i) {
    dest->eff1[i] = src->eff1[i];
    dest->eff2[i] = src->eff2[i];
    dest->threshold[i] = src->threshold[i];
  }

  dest->area = src->area;
  dest->sse = src->sse;
  dest->mse = src->mse;
  dest->bestsp.value = src->bestsp.value;
  dest->bestsp.index = src->bestsp.index;

  return;
}

/* This function takes the outputs from the neural network at the first
   argument computes the integral of the curve "eff1 (y) versus
   (1-eff2) (x)". This integral is subtracted from the integral of
   y=1 between 0 and `(1-eff2[last])' (the last efficiency value for
   eff2). This is what's output. The integration method uses the Simpson's 1/3
   rule (Newton-Cotes, with 2nd. order approximation).  This criterium should
   be a better one than the sp value, since it computes the integral over the
   whole discrimination area instead of a single point. Of course it can't be
   negative.  */
void eval_area (const outputdata_t* out, analysis_t* ap)
{
  real_t under; /* the area under the curve */
  int ninterval; /* the number of intervals used */

  /* Decides the number of intervals we may use in order to apply the Simpson
     1/3 Integration Rule. ninterval must be even */
  ninterval = ANALYSIS_NO_THRESHOLDS;
  if ( ! iseven(ninterval) ) --ninterval;
  
  /* Calculate the area under the curve. The interval (or step) is taken from
     the first one since the curve is calculated backwards. */
  under = int_simpson (ap->eff1, ninterval, 
		       (1-ap->eff2[0])/(real_t)ANALYSIS_NO_THRESHOLDS);

  /* Now we calculate the difference from y=1 and return */
  ap->area = 1 - ap->eff2[0] - under;

  return;
}

/* This function will check on an integer. If it's even, the output is true,
   casen not, the output is false (0). */
int iseven (const int n)
{
  div_t d;
  d = div(n,2);
  if ( d.rem == 0 ) return 1;
  else return 0;
}
 
/* Now we compute the area under that curve using the Simpon's 1/3 Rule
   Repeated (of course) over all intervals. The formula is:
   I = (h/3) * { f(x0) + f(xm) + 4*[f(x1) + f(x3) + ... + f(xm-1)] + 2*[f(x2) +
   + f(x4) + ... + f(xm-2)] }.
   I've got to find out the even and odd points. The value of m is
   ANALYSIS_NO_THRESHOLDS, if ANALYSIS_NO_THRESHOLDS is even or
   ANALYSIS_NO_THRESHOLDS-1 in the case it's odd, since this integration method
   requires m to be even. The parameters are: y values of the curve, the number
   of y values (checked for evenness), the step h used to extract those
   values. */ 
real_t int_simpson (const real_t* y, const int ny, const real_t h)
{
  register int i; /* iterator */
  real_t oddacc = 0.; /* accumulator for odd positions */
  real_t evenacc = 0.; /* accumulator for even positions */

  for (i=1; i<ny; ++i) { /* for the remaining values */
    if ( iseven(i) ) evenacc += y[i];
    else oddacc += y[i];
  }

  return ( (h/3)*(y[0] + y[ny] + 4*oddacc + 2*evenacc) );
}

/* Print efficiencies on file indicated by the first argument */
void print_efficiencies (const char* fname, const analysis_t* ap)
{
  register int i; /* iterator */
  FILE* fp;

  /* open file */
  if ( NULL == (fp=fopen(fname,"w") ) ) {
    perror("[main] switching to stdout");
    fp = stdout;
  }
  
  /* write efficiencies */
  for (i=0; i<=ANALYSIS_NO_THRESHOLDS; ++i)
    fprintf(fp, "%e %e\n", ap->eff1[i], ap->eff2[i]);

  if (fp != stdout) fclose(fp);

  return;
}










