/* Hello emacs, this is -*- c -*- */
/* $Id: main.c,v 1.1 2001/07/13 15:48:02 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

/* This program will run an Artificial Neural Network with one output and a
   configurable number of input and hidden nodes. The activation funtion will
   always be a tgh() and the output function, the identity line. */

#include <stdlib.h>
#include <stdio.h>
#include <sys/times.h>
#include <time.h>
#include <string.h>

#include "parameter.h"
#include "neural.h"
#include "lut.h"
#include "real.h"
#include "event.h"

int main (int argc, char** argv)
{
  register int i; /* dump iterator */
  parameter_t* parameters;
  neuralnet_t* net;

  inputdata_t* input;
  outputdata_t* output;

  clock_t acc, start; /* for time measurement */

  /* This will initialize and load the command line parameters */
  parameters = getparam(argc, (const char**)argv);

  /* This will print some informative information */
#ifdef _PRECISION_32BIT
  fprintf(stderr, "[main] Using 32-bit logic\n");
#endif
#ifdef _PRECISION_64BIT
  fprintf(stderr, "[main] Using 64-bit logic\n");
#endif

  /* Configure neural net and lookup table, if needed */
  net = loadnet(parameters->cfile);
  change_act(parameters->act_method);
  if (parameters->act_method == 'l') {
    config_lut(parameters->lut_prec, parameters->lut_limit);
    fprintf(stderr, "[main] LUT precision: %e, limit: %e\n",
	    parameters->lut_prec, parameters->lut_limit);
  }

  /* Readin the input events */
  input = fread_dat(parameters->ifilename, net->input_dim);

  /* Now we start processing extracting the network output */
  start = clock();
  for (i=0; i<100; ++i) output = batch_propagate(input, net);
  acc = clock() - start;

  /* Save the output to file */
  fwrite_output(parameters->ofilename, output);

  /* Calculate the relevance of input over net performance */
  if (parameters->relevance) {
    register int j; /* iterator */
    real_t* r; /* placeholder */
    FILE* rfile;
    char* filename; /* for the output filename */
    asprintf(&filename, "%s-relevance.dat", parameters->hint);

    if ( (rfile=fopen(filename, "w")) == NULL ) {
      perror ("[main] relevance file");
      exit (EXIT_FAILURE);
    }
    r = relevance(input, output, net);
    for (j=0; j<net->input_dim; ++j) fprintf(rfile, "%e\n", r[j]);
    free(r);
    free(filename);
    fclose(rfile);
    fprintf(stderr, "[main] Finished relevance evaluation\n");
  }

  /* Calculate the importance of input over net performance */
  if (parameters->importance) {
    register int j; /* iterator */
    real_t* imp; /* placeholder */
    FILE* impfile;
    char* filename; /* for the output filename */
    asprintf(&filename, "%s-importance.dat", parameters->hint);
    
    if ( (impfile=fopen(filename, "w")) == NULL ) {
      perror ("[main] importance file");
      exit (EXIT_FAILURE);
    }
    imp = importance(input, output, net);
    for (j=0; j<net->input_dim; ++j) fprintf(impfile, "%e\n", imp[j]);
    free(imp);
    free(filename);
    fclose(impfile);
    fprintf(stderr, "[main] Finished importance evaluation\n");
  }

  /* This will output some usefull comments to the user */
  fprintf(stderr, "[main] The processing for %s is now finished.\n",
	  (strcmp(parameters->ifilename,"")==0)?"stdin":parameters->ifilename);
  fprintf(stderr, "[main] The net outputs were%s saved%s%s.\n",
	  (strcmp(parameters->ofilename,"")==0)?" not":"", 
	  (strcmp(parameters->ofilename,"")==0)?"":" at ", 
	  (strcmp(parameters->ofilename,"")==0)?"":parameters->ofilename);
  fprintf(stderr, "[main] The total processing time was %e seconds.\n",
	  ((double)acc)/CLOCKS_PER_SEC);
  fprintf(stderr, "[main] The mean of the processing time is %e us.\n",
	  (10000*(double)acc)/((input->n1+input->n2)*CLOCKS_PER_SEC));
  fprintf(stderr, "[main] %d Events were analysed in total.\n",
	  (input->n1+input->n2)*100);
  fprintf(stderr, "Bye...\n");

  /* This will terminate all internal run parameters */
  termparam(parameters);
  clear_lut();
  free_inputdata_t(input); 
  free_outputdata_t(output);

  exit(EXIT_SUCCESS);

} /* main() */
  
