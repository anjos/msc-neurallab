/* Hello emacs, this is -*- c -*- */
/* $Id: analysis.h,v 1.1 2001/07/13 15:48:02 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

/* This module contains functions for output data analysis. */

#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "real.h"
#include "event.h"

/* this will keep the number of thresholds to analyze */
#define ANALYSIS_NO_THRESHOLDS 200

typedef struct sp_t {
  real_t value; /* the value of SP found at this run */
  int index; /* the index from analysis_t that guarantees it */
} sp_t;

typedef struct analysis_t {
  sp_t bestsp;
  /* the efficiency for class 1 of the current run */
  real_t eff1[ANALYSIS_NO_THRESHOLDS+1];
  /* the efficiency for class 1 of the current run */
  real_t eff2[ANALYSIS_NO_THRESHOLDS+1];
  /* the thresholds that guarantee those */
  real_t threshold[ANALYSIS_NO_THRESHOLDS+1];
  real_t area; /* the area under y=1 and the efficiency curve by Simpon's
		  integration */
  real_t sse; /* the sum of square errors */
  real_t mse; /* the mean square error */
} analysis_t;

typedef struct analysis_parameter_t {
  bool_t eval_sp;
  bool_t eval_area;
} analysis_parameter_t;

/* This function will analyze the output data at the second argument, placing
   the results at the first one.
   The analysis consists of the following steps: 
   1. For the ANALYSIS_NO_THRESHOLDS, the procedure will scan the interval
   between the targets, computing the efficiency between class 1 and class 2
   discrimination;
   2. During this procedure, if it will save the best SP products it finds in
   bestsp, keeping the indexes for that value and the value itself if the last
   argument is true on "eval_sp";
   3. Will evaluate the area under the ROC if the last argument is true on
   "eval_area";
   4. Will calculate the MSE and SSE of input. The MSE is evaluated using the
   SSE of both classes, dividing each by the number of entries in each class
   and evaluating the mean between both sets.
   *** WE EXPECT THAT CLASS 1 HAS TARGETS GREATER THAN CLASS 2 */
void analyze (analysis_t*, const outputdata_t*, const analysis_parameter_t*);

/* Print efficiencies on file indicated by the first argument */
void print_efficiencies (const char*, const analysis_t*);

/* Copies a analysis_t into another */
void copy_analysis (const analysis_t*, analysis_t*);

#endif /* ANALYSIS_H_ */


