/* Hello emacs, this is -*- c -*- */
/* $Id: parameter.c,v 1.2 2001/08/02 03:46:30 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

#include <popt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "parameter.h"
#include "event.h"

/* Some "local global" stuff: The next variable declares a parameter the valid
   option flags. They have to be global or local to main. */
static struct poptOption opttable[] = {
  /* longName, shortName, argInfo, arg, val, descript, argDescript */
  {"maxsteps", 'a', POPT_ARG_LONG, NULL, 'a',
   "defines the maximum number of steps to train [default: 500]", 
   "long int (>0)"}, 
  {"config-file", 'c', POPT_ARG_STRING, NULL, 'c',
   "output filename for configuration storage [default: config.dat]",
   "filename"},  
  {"lr-decay", 'd', PREC_POPT_REAL, NULL, 'd',
   "defines the Learning Rate decay [default: 1.0]", "float (>0)"},
  {"eff-file", 'e', POPT_ARG_STRING, NULL, 'e',
   "output filename for training efficinency storage [default: effic.dat]",
   "filename"},
  {"no-inputs", 'i', POPT_ARG_LONG, NULL, 'i',
   "the number of input nodes [default: 10]", "int (>0)"},
  {"no-hidden", 'j', POPT_ARG_LONG, NULL, 'j',
   "the number of hidden nodes [default: 5]", "int (>0)"},
  {"learn-rate", 'l', PREC_POPT_REAL, NULL, 'l',
   "defines the Learning Rate [default: 0.1]", "float (>0)"},
  {"momentum", 'm', PREC_POPT_REAL, NULL, 'm',
   "defines the training momentum [default: 0]", "float (>0)"},
  {"run-file", 'n', POPT_ARG_STRING, NULL, 'n',
   "output filename for running data storage [default: stdout]",
   "filename"},
  {"batch", 'p', POPT_ARG_LONG, NULL, 'p',
   "defines the batch size (patterns per update) [default: 10]", "long int (>0)"},
  {"epoch", 'q', POPT_ARG_LONG, NULL, 'q',
   "defines the epoch size (updates per parameter update) [default: 100]", "long int (>0)"},
  {"train-file", 'r', POPT_ARG_STRING, NULL, 'r', 
   "input data filename for training data [default: NO DEFAULT]", "filename"},
  {"use-sp", 's', POPT_ARG_NONE, NULL, 's',
   "Use the SP product value to find best network [DEFAULT] (faster)", NULL},
  {"test-file", 't', POPT_ARG_STRING, NULL, 't',
   "input data filename for test data [default: NO DEFAULT]", "filename"},
  {"use-area", 'x', POPT_ARG_NONE, NULL, 'x',
   "Use the ROC area value to find best network (better results)", NULL},
  POPT_AUTOHELP
  {NULL, 0, 0, NULL, 0} };

/* This function sets the default values for the parameters at startup. */
void initparam (parameter_t*);

/* The next function just copies the contentes one parameter variable into the
   other. */
parameter_t* cpparam (parameter_t*, const parameter_t*);

/* This function tests some aspectes of the parameters you selected, not
   allowing strange things to happen. If, for instance, you forget the training
   and test sets file, it warns you... */
void test_parameters (parameter_t*);

/* This will dump to the file named p->config_filename, the configuration of
   this run. */
void dump_config (parameter_t*);

/* **********************
   FUNCTION IMPLEMENATION
   ********************** */

void initparam (parameter_t* p)
{
  strncpy(p->train_filename, "", MAX_FILENAME_SIZE);
  strncpy(p->test_filename, "", MAX_FILENAME_SIZE);
  strncpy(p->eff_filename, "effic.dat", MAX_FILENAME_SIZE);
  strncpy(p->config_filename, "config.dat", MAX_FILENAME_SIZE);
  strncpy(p->run_filename, "", MAX_FILENAME_SIZE);

  /* Dump configuration to default file and run data to screen */
  p->runfile = stdout;
  p->cfile = NULL;

  p->anpar.eval_sp = FALSE;
  p->anpar.eval_area = FALSE;

  /* Network training parameters */
  p->lr = 0.1;
  p->momentum = 0;
  p->lr_decay = 1;
  p->batch = 10;
  p->epoch = 100;
  p->maxsteps = 500;
  
  /* Network architecture is (10-5-1) by default */
  p->ninputs = 10;
  p->nhidden = 5;

}

parameter_t* getparam (int argc, const char** argv)
{
  parameter_t p; /* This is what will be used inside the function */
  parameter_t* retval; /* This will temporarily hold the returned values */
  int c = 0; /* for returned values by poptGetNextOpt() */
  poptContext context; /* context for parsing command-line options */

  /* Initialize the parameters */
  initparam(&p);

  /* Create our current context */
  context = poptGetContext(NULL, argc, argv, opttable, 0);
  
  /* If the user signaled not to understand the app, give usage to him/her */
  if (argc < 2) {
    poptPrintUsage(context, stderr, 0);
    exit(EXIT_SUCCESS);
  }

  /* Now do options processing */
  while ( (c=poptGetNextOpt(context)) >= 0 ) {
    switch (c) {
    case 'a':
      p.maxsteps = strtol(poptGetOptArg(context),NULL,10);
      if ( p.maxsteps <= 0 ) {
	fprintf(stderr, "[parameter] Invalid value for STEPS => %ld\n",
		p.maxsteps);
	exit(EXIT_FAILURE);
      }
      break;

    case 'c':
      strncpy(p.config_filename, poptGetOptArg(context), MAX_FILENAME_SIZE);
      break;

    case 'd': 
      p.lr_decay = (real_t)strtod(poptGetOptArg(context),NULL);
      if ( p.lr_decay < 0 ) {
	fprintf(stderr, "[parameter] Invalid value for LR Decay => %e\n", 
		p.lr_decay);
	exit(EXIT_FAILURE);
      }
      break;	
      
    case 'e':
      strncpy(p.eff_filename, poptGetOptArg(context), MAX_FILENAME_SIZE);
      break;
      
    case 'i':
      p.ninputs = (int)strtol(poptGetOptArg(context),NULL,10);
      if ( p.ninputs <= 0 ) {
	fprintf(stderr, "[parameter] Invalid value for # inputs => %d\n", 
		p.ninputs);
	exit(EXIT_FAILURE);
      }
      break;	

    case 'j':
      p.nhidden = (int)strtol(poptGetOptArg(context),NULL,10);
      if ( p.nhidden <= 0 ) {
	fprintf(stderr, "[parameter] Invalid value for # hidden => %d\n",
		p.nhidden);
	exit(EXIT_FAILURE);
      }
      break;	

    case 'l':
      p.lr = (real_t)strtod(poptGetOptArg(context),NULL);
      if ( p.lr <= 0 ) {
	fprintf(stderr, "[parameter] Invalid value for LR => %e\n",
		p.lr);
	exit(EXIT_FAILURE);
      }
      break;

    case 'm':
      p.momentum = (real_t)strtod(poptGetOptArg(context),NULL);
      if ( p.momentum < 0 ) {
	fprintf(stderr, "[parameter] Invalid value for momentum => %e\n",
		p.momentum);
	exit(EXIT_FAILURE);
      }
      break;

    case 'n':
      strncpy(p.run_filename, poptGetOptArg(context), MAX_FILENAME_SIZE);
      if ( (p.runfile=fopen(p.run_filename, "w")) == NULL ) {
	perror("[parameter] runfile");
	exit(EXIT_FAILURE);
      }
      break;

    case 'p':
      p.batch = (int)strtol(poptGetOptArg(context),NULL,10);
      if ( p.batch < 1 ) {
	fprintf(stderr, "[parameter] Invalid value for batch => %ld\n",
		p.batch);
	exit(EXIT_FAILURE);
      }
      break;	

    case 'q':
      p.epoch = (int)strtol(poptGetOptArg(context),NULL,10);
      if ( p.epoch < 1 ) {
	fprintf(stderr, "[parameter] Invalid value for epoch => %ld\n",
		p.epoch);
	exit(EXIT_FAILURE);
      }
      break;

    case 'r':
      strncpy(p.train_filename, poptGetOptArg(context), MAX_FILENAME_SIZE);
      break;

    case 's':
      p.anpar.eval_sp = TRUE;
      break;

    case 't':
      strncpy(p.test_filename, poptGetOptArg(context), MAX_FILENAME_SIZE);
      break;

    case 'x':
      p.anpar.eval_area = TRUE;
      break;

    default:
      fprintf(stderr, "[parameter] Return value _not_ expected => %d\n", c);
      exit(EXIT_FAILURE);
    }
  }

  /* check for errors during option processing */
  if ( c < -1 ) {
    /* an error occurred during option processing */
    fprintf(stderr, "[parameter] %s: %s\n", 
	    poptBadOption(context, POPT_BADOPTION_NOALIAS), poptStrerror(c));
    exit(EXIT_FAILURE);
  }

  poptFreeContext(context);

  /* Readin data files */
  p.train = fread_dat(p.train_filename, p.ninputs);
  p.test = fread_dat(p.test_filename, p.ninputs);

  /* Now we must copy the contents from the static memory into the dynamic */
  retval = (parameter_t*) malloc (sizeof(parameter_t));
  cpparam(retval, &p);

  /* Test some aspects of parameters */
  test_parameters(retval);

  /* Dump parameters into file */
  dump_config(retval);

  /* If we got here, we can safely return */
  return (retval);
} /* getparam() */

parameter_t* cpparam (parameter_t* dest, const parameter_t* src)
{
  dest->lr = src->lr;
  dest->momentum = src->momentum;
  dest->lr_decay = src->lr_decay;
  dest->epoch = src->epoch;
  dest->batch = src->batch;
  dest->maxsteps = src->maxsteps;
  dest->ninputs = src->ninputs;
  dest->nhidden = src->nhidden;
  dest->runfile = src->runfile;
  dest->train = src->train;
  dest->test = src->test;
  dest->anpar.eval_sp = src->anpar.eval_sp;
  dest->anpar.eval_area = src->anpar.eval_area;
  strncpy(dest->train_filename, src->train_filename, MAX_FILENAME_SIZE);
  strncpy(dest->test_filename, src->test_filename, MAX_FILENAME_SIZE);
  strncpy(dest->eff_filename, src->eff_filename, MAX_FILENAME_SIZE);
  strncpy(dest->config_filename, src->config_filename, MAX_FILENAME_SIZE);
  strncpy(dest->run_filename, src->run_filename, MAX_FILENAME_SIZE);
  return dest;
}

void termparam (parameter_t* p)
{
  if ( p->runfile != stdout || p->runfile != NULL ) fclose(p->runfile);
  if ( p->cfile != stdout || p->cfile != NULL ) fclose(p->cfile);

  /* The input data has to be freed as well */
  free_inputdata_t(p->train);
  free_inputdata_t(p->test);

  free(p);
  return;
}

/* This function tests some aspectes of the parameters you selected, not
   allowing strange things to happen. If, for instance, you forget the training
   and test sets file, it warns you... */
void test_parameters (parameter_t* p)
{
  if (p->lr > 2) {
    fprintf(stderr,"[parameter] Your LR is greater than 2! (= %e)\n", p->lr); 
    fprintf(stderr,"            Strange things might happen...\n");
  }

  if (p->momentum > 1) {
    fprintf(stderr,"[parameter] Your momentum is greater than 1! (= %e)\n",
	    p->momentum);
    fprintf(stderr,"            Strange things might happen...\n");
  }

  if (p->anpar.eval_sp && p->anpar.eval_area) {
    fprintf(stderr, "[parameter] Can't locate network by SP an Area under ROC");
    fprintf(stderr, "\n[parameter] I Will use \"SP Product\"");
    p->anpar.eval_area = FALSE;
  }

  if (!p->anpar.eval_sp && !p->anpar.eval_area) {
    fprintf(stderr, "[parameter] Next time you have to select a method for");
    fprintf(stderr, " network saving...\n"); 
    fprintf(stderr, "[parameter] I Will use \"SP Product\" this time!");
    p->anpar.eval_sp = TRUE;
  }

  return;
}

/* This will dump to the file named p->config_filename, the configuration of
   this run. */
void dump_config (parameter_t* p)
{
  time_t current_time = time(NULL); /* just a dummy variable for timing */ 

  if ( (p->cfile=fopen(p->config_filename,"w")) == NULL ) {
    perror("[parameter] configuration file");
    exit(EXIT_FAILURE);
  }

  /* SDB version */
  fprintf(p->cfile, "[sdb]\nfloat version = 0.1\n");
  
  fprintf(p->cfile, "\n[files]\n");
  fprintf(p->cfile, "string date = %s", ctime(&current_time));

  fprintf(p->cfile, "string train_file = %s\n", p->train_filename);
  fprintf(p->cfile, "string test_file = %s\n", p->test_filename);
  fprintf(p->cfile, "string effic_file = %s\n", p->eff_filename);
  fprintf(p->cfile, "string runfile = %s\n",
	  (p->runfile==stdout)?"stdout":p->run_filename);

  fprintf(p->cfile, "\n[parameters]\n");
  fprintf(p->cfile, "float learning_rate = %.4e\n", p->lr); 
  fprintf(p->cfile, "float momentum = %.4e\n", p->momentum);
  fprintf(p->cfile, "float lr_decay = %.4e\n", p->lr_decay);
  fprintf(p->cfile, "int steps = %ld\n", p->maxsteps);
  fprintf(p->cfile, "int batch = %ld\n", p->batch);
  fprintf(p->cfile, "int epoch = %ld\n", p->epoch);
  fprintf(p->cfile, "int input_dimension = %d\n", p->ninputs);
  fprintf(p->cfile, "int hidden_dimension = %d\n", p->nhidden);
  fprintf(p->cfile, "string network_saving = ");
  if (p->anpar.eval_sp) fprintf(p->cfile, "SP\n");
  if (p->anpar.eval_area) fprintf(p->cfile, "Area\n");

  /* Do _not_ close the configuration file, somethings will be added later */
  /* fclose(p->cfile); */

  return;
}
