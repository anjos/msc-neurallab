/* Hello emacs, this is -*- c -*- */

/* $Id: parameter.h,v 1.1 2001/07/13 15:48:02 andre Exp $ */

/* This modules controls program initialization. It uses the popt library to
   setup argument options and so. Be aware! */

/* André Rabello <Andre.Rabello@ufrj.br> */

#ifndef PARAMETER_H_ /* our doorkeeper */
#define PARAMETER_H_

#include <stdlib.h>
#include "neural.h"
#include "real.h"

#define MAX_FILENAME_SIZE 80

typedef struct parameter_t {
  char ifilename[MAX_FILENAME_SIZE]; /* The input events filename */
  char cfilename[MAX_FILENAME_SIZE]; /* The network configuration filename */
  FILE* cfile; /* A pointer (if non-null) to the configuration file */
  char ofilename[MAX_FILENAME_SIZE]; /* Where I should store the data */
  char hint[MAX_FILENAME_SIZE]; /* The hint for extra output (appending str) */
  bool_t relevance; /* Do I have to calculate the input relevance to net
		       performance?. Default is NO. */
  bool_t importance; /* Do I have to calculate the input importance to net
			performance?. Default is NO. */
  int act_method; /* This holds the current activation method */
  real_t lut_prec; /* This will configure the lookup table precision */
  real_t lut_limit; /* This will configure the lookup table limits */
} parameter_t;

/* This function will take argc and argv from main() and will parse them
   composing the parameters of the current run. The arguments should specify a
   configuration file with the number of input and hidden units of the neural
   network and all weights. Also, the arguments should point to a input file,
   for network testing and in which way the activation of each neuron is
   performed. */
parameter_t* getparam (int argc, const char** argv);

/* The next function will simply de-allocate any pre-allocated space, close
   opened files, etc. */
void termparam (parameter_t*);

#endif /* PARAMETER_H_ */

