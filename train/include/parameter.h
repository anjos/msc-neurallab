/* Hello emacs, this is -*- c -*- */

/* $Id: parameter.h,v 1.1 2001/07/13 15:48:02 andre Exp $ */

/* This modules controls program initialization. It uses the popt library to
   setup argument options and so. Be aware! */

/* André Rabello <Andre.Rabello@ufrj.br> */

#ifndef PARAMETER_H_ /* our doorkeeper */
#define PARAMETER_H_

#include <stdlib.h>
#include "analysis.h"
#include "real.h"
#include "event.h"

#define MAX_FILENAME_SIZE 80

typedef struct parameter_t {
  char train_filename[MAX_FILENAME_SIZE]; /* The training data set */
  char test_filename[MAX_FILENAME_SIZE]; /* The test data set */
  char eff_filename[MAX_FILENAME_SIZE]; /* The file where the efficiencies are
					   going to be stored */
  char config_filename[MAX_FILENAME_SIZE]; /* Where to put config. data */
  FILE* cfile; /* A pointer (if non-null) to the configuration file */

  char run_filename[MAX_FILENAME_SIZE]; /* Where I should store the run data */
  FILE* runfile; /* A pointer (if non-null) to the output file */
  
  real_t lr; /* The learning rate (LR) */
  real_t momentum; /* The momentum parameter */
  real_t lr_decay; /* LR decay rate */

  long int batch; /* The size of the batch (patterns per update) */
  long int epoch; /* The size of the epoch (updates per parameter change) */
  long int maxsteps; /* The number of steps to train the network */
  int ninputs; /* The number of input nodes */
  int nhidden; /* The number of hidden layers */

  /* Parameters for network analysis */
  analysis_parameter_t anpar;

  /* Variables for data manipulation only */
  inputdata_t* train;
  inputdata_t* test;

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









