/* Hello emacs, this is -*- c -*- */
/* $Id: parameter.c,v 1.1 2001/07/13 15:48:02 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

#include <popt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "parameter.h"

/* Some "local global" stuff: The next variable declares a parameter the valid
   option flags. They have to be global or local to main. */
static struct poptOption opttable[] = {
  /* longName, shortName, argInfo, arg, val, descript, argDescript */
  {"input-file", 'i', POPT_ARG_STRING, NULL, 'i', 
   "input data vectors filename [default: stdin]", "filename"},
  {"output-file", 'o', POPT_ARG_STRING, NULL, 'o',
   "output data filename [default: stdout]", "filename"},
  {"hint", 'k', POPT_ARG_STRING, NULL, 'k',
   "hint output filename [default: '']", "string"},
  {"net-file", 'n', POPT_ARG_STRING, NULL, 'n',
   "neural network configuration filename [default: config.net]", "filename"},
  {"act-builtin", 'f', POPT_ARG_NONE, NULL, 'f', 
   "selects builtin tanh() function as the activation method", NULL},
  {"act-lookup", 'l', POPT_ARG_NONE, NULL, 'l',
   "use a lookup table to implement tanh()", NULL},
  {"lut-precision", 'p', PREC_POPT_REAL, NULL, 'p',
   "defines the LUT precision [default: 0.01]", "float (>0)"},
  {"lut-end", 'e', PREC_POPT_REAL, NULL, 'e',
   "defines the ending point for the LUT [default: 6.0]", "float (>0)"},
  {"relevance", 'r', POPT_ARG_NONE, NULL, 'r',
   "evaluate the input relevance to net performance (MSE)", NULL},
  {"importance", 't', POPT_ARG_NONE, NULL, 't',
   "evaluate the input importance to net performance (EFFICIENCY)", NULL},
  POPT_AUTOHELP
  {NULL, 0, 0, NULL, 0} };

/* This function sets the default values for the parameters at startup. */
void initparam (parameter_t*);

/* The next function just copies the contentes one parameter variable into the
   other. */
parameter_t* cpparam (parameter_t*, const parameter_t*);

/* **********************
   FUNCTION IMPLEMENATION
   ********************** */

void initparam (parameter_t* p)
{
  strncpy(p->ifilename, "", MAX_FILENAME_SIZE);
  strncpy(p->cfilename, "config.net", MAX_FILENAME_SIZE);
  /* give hint for other output files => 'hint'-type.dat */
  strncpy(p->hint, "", MAX_FILENAME_SIZE);
  p->cfile = NULL;
  strncpy(p->ofilename, "", MAX_FILENAME_SIZE);

  p->act_method = 'f'; /* this means use 'function' */
  p->lut_prec = 0.01; /* The default precision for the lookup table */
  p->lut_limit = 6.0; /* The limits are set from -6 till +6, after that, the
			 results are truncated to -1 or +1 respectively */

  p->relevance = FALSE; /* do not calculate the input relevance by default */
  p->importance = FALSE; /* do not calculate the input importance by default */
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
    case 'l': /* use lookup */
    case 'f': /* use builtin function */
      p.act_method = c;
      break;

    case 'p': /* The precision can't be less or equal zero! */
      p.lut_prec = (real_t)strtod(poptGetOptArg(context),NULL);
      if ( p.lut_prec <= 0 ) {
	fprintf(stderr, "[main] Invalid table precision => %e\n",
		p.lut_prec);
	exit(EXIT_FAILURE);
      }
      break;

    case 'e': /* The table boundary */
      p.lut_limit = (real_t)strtod(poptGetOptArg(context),NULL);
      if ( p.lut_limit <= 0 ) {
	fprintf(stderr, "[main] Invalid table limit => %e\n",
		p.lut_limit);
	exit(EXIT_FAILURE);
      }
      break;	

    case 'i': /* The input file has been changed */
      strncpy(p.ifilename, poptGetOptArg(context), MAX_FILENAME_SIZE);
      break;

    case 'k': /* The hint for extra output files has been changed */
      strncpy(p.hint, poptGetOptArg(context), MAX_FILENAME_SIZE);
      break;

    case 'o': /* The output file has been changed */
      strncpy(p.ofilename, poptGetOptArg(context), MAX_FILENAME_SIZE);
      break;

    case 'n': /* The neural network config file has been changed */
      strncpy(p.cfilename, poptGetOptArg(context), MAX_FILENAME_SIZE);
      if ( (p.cfile=fopen(p.cfilename, "r")) == NULL ) {
	perror("[main] configfile");
	exit(EXIT_FAILURE);
      }
      break;

    case 'r': /* calculate the input relevance to net performance */
      p.relevance = TRUE;
      break;

    case 't': /* calculate the input relevance to net performance */
      p.importance = TRUE;
      break;

    default:
      fprintf(stderr, "[main] Return value _not_ expected => %d\n", c);
      exit(EXIT_FAILURE);
    }
  }

  /* check for errors during option processing */
  if ( c < -1 ) {
    /* an error occurred during option processing */
    fprintf(stderr, "[main] %s: %s\n",
	    poptBadOption(context, POPT_BADOPTION_NOALIAS), poptStrerror(c));
    exit(EXIT_FAILURE);
  }

  poptFreeContext(context);

  /* Now we must copy the contents from the static memory into the dynamic */
  retval = (parameter_t*) malloc (sizeof(parameter_t));
  cpparam(retval, &p);

  /* If we got here, we can safely return */
  return (retval);
} /* getparam() */

parameter_t* cpparam (parameter_t* dest, const parameter_t* src)
{
  strncpy(dest->ifilename, src->ifilename, MAX_FILENAME_SIZE);
  strncpy(dest->ofilename, src->ofilename, MAX_FILENAME_SIZE);
  strncpy(dest->cfilename, src->cfilename, MAX_FILENAME_SIZE);
  strncpy(dest->hint, src->hint, MAX_FILENAME_SIZE);
  dest->cfile = src ->cfile;
  dest->act_method = src->act_method;
  dest->lut_prec = src->lut_prec;
  dest->lut_limit = src->lut_limit;
  dest->relevance = src->relevance;
  dest->importance = src->importance;
  return dest;
}

void termparam (parameter_t* p)
{
  if ( p->cfile != NULL ) fclose(p->cfile);
  free(p);
  return;
}
