#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <f2c.h>

#include "jetnet.h"
#include "real.h"
#include "event.h"
#include "analysis.h"
#include "parameter.h"

/* Trains a neural network */

/* This function trains the network with the input data set */
void train_net (const inputdata_t*);

/* This function will place a random input event in OIN and its output vector
   (the target) in OUT. */
void place_random (const inputdata_t*);

/* This will place a generic const event into the right positions for training
   or testing. */
void place_event (const event_t*, const int);

/* Initializes JETNET parameters and newly created network */
void init_net (const parameter_t*);

/* This will test the network and will dump to the screen the results of such
   testing. The results will be thrown over the second argument. */
void test_net (const inputdata_t*, analysis_t*, const analysis_parameter_t*);

/* This will save the network output into a file named after the first
   input. The data that will be used should be given as a second parameter. */
void save_output(const char*, const inputdata_t*);

/* Saves the network on to a file */
void save_net(void);

/* This funtion just returns a valid long based on a string. It's an
   implemenatation of atol(), but with the verification of strtol(). */
double to_valid_double(const char*);

/* Converts a string into a valid double or issue an error message */
long to_valid_long(const char*);

extern struct jndat1_1_ jndat1_;

int main (int argc, char** argv)
{
  analysis_t train_an; /* where the train analysis is going to be stored */
  analysis_t train_best; /* where the train analysis is going to be stored */
  analysis_t an; /* where the analysis is going to be stored */
  analysis_t best; /* where I keep the smaller area so far */

  /* The parameters I shall use */
  parameter_t* params;

  /* These are training iterators */
  int step;
  int lastsaved = 0; /* the last saved step */
  div_t runtest; /* shall I run the tests? */

  /* Get parameters from input flags */
  params = getparam(argc, (const char**)argv);

  /* Initializing the random number genearator, for this program */
  srand(time(0));

  best.bestsp.value = 0; /* just to initialize it... */
  best.area = 1;

  init_net(params); /* Initializes the network and its configuration */

  fprintf(stderr, " => Wait while I'm processing...\n");
  fprintf(stderr, "%7d/%7ld", 0, params->maxsteps);

  /* Ok, agora vamos ao treinamento */
  for (step=1; step <= params->maxsteps; ++step) {

    train_net(params->train);
    runtest = div(step, jndat1_1.mstjn[1]);
    if (!runtest.rem) { /* run test analysis, will dump to screen! */
      fprintf(stderr, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
      fprintf(stderr, "%7d/%7ld", step, params->maxsteps); 

      /* tests net with train file and save the resulting data */
      test_net(params->train, &train_an, &(params->anpar));
      fprintf(params->runfile, "%d \n", runtest.quot);
      fprintf(params->runfile, " %e %e %e %e %e %e %d\n",
	      train_an.bestsp.value,
	      train_an.threshold[train_an.bestsp.index],
	      train_an.mse,
	      train_an.eff1[train_an.bestsp.index],
	      train_an.eff2[train_an.bestsp.index],
	      train_an.area, -1);

      /* tests net with test file and save the resulting data */
      test_net(params->test, &an, &(params->anpar));

      /* save unconditionally the results to file */
      fprintf(params->runfile, " %e %e %e %e %e %e ", an.bestsp.value,
	      an.threshold[an.bestsp.index], an.mse, an.eff1[an.bestsp.index],
	      an.eff2[an.bestsp.index], an.area);

      /* If the net is good, save it */
      if (( params->anpar.eval_area && (an.area < best.area) ) ||
	  (params->anpar.eval_sp && (an.bestsp.value > best.bestsp.value) )) {
	/* Save all */
	copy_analysis(&an, &best);
	copy_analysis(&train_an, &train_best);
	print_efficiencies(params->eff_filename, &best);
	save_net();
	lastsaved = step;
	fprintf(params->runfile, "%d\n", 1); /* network saved flag */
      }

      else {
	fprintf(params->runfile, "%d\n", 0); /* net NOT saved flag */
      }
    }
  }

  fprintf(stderr, "\n");

  /* Now I will save the network outputs into some files*/
  /*    save_output("jetnet-out-train.dat", params->train); */
  /*    save_output("jetnet-out-test.dat", params->test); */

  /* This will print some relevant information into file */
  fprintf(params->cfile, " +- The LAST network was SAVED at step:: %d\n",
	  lastsaved);

  fprintf(params->cfile, " +-+- <<TRAIN RESULTS>>\n");
  fprintf(params->cfile, "   +- SP:: %.2f\n", train_best.bestsp.value);
  fprintf(params->cfile, "   +- THRESH:: %.3f\n",
	  train_best.threshold[train_best.bestsp.index]);
  fprintf(params->cfile, "   +- MSE:: %.2e\n", train_best.mse);
  fprintf(params->cfile, "   +- PE:: %.2e\n", train_best.area);
  fprintf(params->cfile, "   +- Eff1 (%%):: %.2f\n",
	  100*train_best.eff1[train_best.bestsp.index]);
  fprintf(params->cfile, "   +- Eff2 (%%):: %.2f\n",
	  100*train_best.eff2[train_best.bestsp.index]);
  fprintf(params->cfile, "   +- BackRate2 (kHz):: %.2f\n", 
	  25*(1-train_best.eff2[train_best.bestsp.index]));

  fprintf(params->cfile, " +-+- <<TEST RESULTS>>\n");
  fprintf(params->cfile, "   +- SP:: %.2f\n", best.bestsp.value);
  fprintf(params->cfile, "   +- THRESH:: %.3f\n", 
	  best.threshold[best.bestsp.index]);
  fprintf(params->cfile, "   +- MSE:: %.2e\n", best.mse);
  fprintf(params->cfile, "   +- PE:: %.2e\n", best.area);
  fprintf(params->cfile, "   +- Eff1 (%%):: %.2f\n",
	  100*best.eff1[best.bestsp.index]);
  fprintf(params->cfile, "   +- Eff2 (%%):: %.2f\n",
	  100*best.eff2[best.bestsp.index]);
  fprintf(params->cfile, "   +- BackRate2 (kHz):: %.2f\n", 
	  25*(1-best.eff2[best.bestsp.index]));
  
  /* Free the input data */
  termparam(params);

  exit(EXIT_SUCCESS);
}

/* Initializes JETNET parameters and newly created network */
void init_net (const parameter_t* params)
{ 
  /* Network configuration per se */
  jndat1_1.mstjn[1] = params->batch; /* number of inputs before updating the
					net */ 
  jndat1_1.mstjn[8] = params->epoch; /* number of weight updates before
					updating the training parameters */
  jndat1_1.mstjn[2] = 2; /* overall activation function = tanh() */
  jndat1_1.mstjn[9] = params->ninputs; /* number of nodes in input layer */
  jndat1_1.mstjn[10] = params->nhidden; /* number of nodes in hidden layer */
  jndat1_1.mstjn[11] = 1; /* number of nodes at output layer */

  /* Training parametrisation */
  jndat1_1.parjn[0] = params->lr; /* learning rate */
  jndat1_1.parjn[1] = params->momentum; /* momentum */
  jndat1_1.parjn[3] = 1.; /* window for weight initialization */
  jndat1_1.parjn[10] = params->lr_decay; /* Learning Rate decay */

  /* Initialize network */
  jninit_();
  fprintf(stderr, "Network Initialization done!\n");
  return;
}

/* This function trains the network with the input data set */
void train_net (const inputdata_t* d)
{
  place_random(d);
  jntral_();
}

void test_net (const inputdata_t* d, analysis_t* an, 
	       const analysis_parameter_t* anpar)
{
  register int i; /* iterator */
  outputdata_t out; /* the outputs from the network */

  out.n1 = d->n1;
  out.target1 = d->class1[0].target; /* ignores target variance:read event.h */
  if ( NULL == (out.class1 = (real_t*) alloca (d->n1 * sizeof(real_t))) ) {
    perror("[test_net] output allocation for class 1");
    exit(EXIT_FAILURE);
  }

  out.n2 = d->n2;
  out.target2 = d->class2[0].target; /* ignores target variance:read event.h */
  if ( NULL == (out.class2 = (real_t*) alloca (d->n2 * sizeof(real_t))) ) {
    perror("[test_net] output allocation for class 2");
    exit(EXIT_FAILURE);
  }

  for(i=0; i<d->n1; ++i) { /* test events of class 1 */
    place_event(&d->class1[i], d->event_size);
    jntest_(); /* call training routine */
    out.class1[i] = jndat1_1.out[0];
  }

  for(i=0; i<d->n2; ++i) { /* test events of class 2 */
    place_event(&d->class2[i], d->event_size);
    jntest_(); /* call training routine */
    out.class2[i] = jndat1_1.out[0];
  }

  /* Now I got the outputs from the network, I just have to analyze them */
  analyze(an, &out, anpar);

  return;
}

/* This will save the network output into a file named after the first
   input. The data that will be used should be given as a second parameter. */
void save_output(const char* fname, const inputdata_t* d) 
{
  register int i; /* iterator */
  outputdata_t out; /* the outputs from the network */
  FILE* fp;
  /* long int save = 10; This is the current file number */

  out.n1 = d->n1;
  out.target1 = d->class1[0].target; /* ignores target variance:read event.h */
  if ( NULL == (out.class1 = (real_t*) alloca (d->n1 * sizeof(real_t))) ) {
    perror("[test_net] output allocation for class 1");
    exit(EXIT_FAILURE);
  }

  out.n2 = d->n2;
  out.target2 = d->class2[0].target; /* ignores target variance:read event.h */
  if ( NULL == (out.class2 = (real_t*) alloca (d->n2 * sizeof(real_t))) ) {
    perror("[test_net] output allocation for class 2");
    exit(EXIT_FAILURE);
  }

  for(i=0; i<d->n1; ++i) { /* test events of class 1 */
    place_event(&d->class1[i], d->event_size);
    jntest_(); /* call training routine */
    out.class1[i] = jndat1_1.out[0];
  }

  for(i=0; i<d->n2; ++i) { /* test events of class 2 */
    place_event(&d->class2[i], d->event_size);
    jntest_(); /* call training routine */
    out.class2[i] = jndat1_1.out[0];
  }

  /* Now I should dump the data, but only if the output filename is not "" */
  if ( strcmp(fname, "") == 0 )  return;

  if ( (fp = fopen(fname,"w") ) == NULL ) {
    perror ("[main] Net output");
    fprintf(stderr, "[main] Can't output %s for writing...\n", fname);
    exit(EXIT_FAILURE);
  }

  fprintf(stderr, "[main] Saving net output at %s... ", fname);
  for(i=0; i<d->n1; ++i) fprintf(fp, "%e\n", out.class1[i]);
  for(i=0; i<d->n2; ++i) fprintf(fp, "%e\n", out.class2[i]);
  fprintf(stderr, "done!");
  fclose(fp);

  return;
}
  

/* Saves the network on to a file */
void save_net(void)
{ /* put network on file */
  long int save = 10; /* This is the current file number */
  jndump_(&save);
}

/* This function will place a random input event in OIN and its output vector
   (the target) in OUT. */
void place_random (const inputdata_t* d)
{
  place_event(get_random(d), d->event_size);
}

/* This will place a generic const event into the right positions for training
   or testing. */
void place_event (const event_t* evp, const int size)
{
  int i; /* iterator */
  for(i=0; i<size; ++i) jndat1_1.oin[i] = evp->data[i];
  jndat1_1.out[0] = evp->target;
}

/* Converts a string int a valid double or issue an error message */
double to_valid_double(const char* str)
{
  char** invalid_number = (char**)malloc(sizeof(char*));
  double number;

  *invalid_number = NULL;

  number = strtod(str,invalid_number);
  if (**invalid_number != '\0') {
    fprintf(stderr,"(util) -%s- is not a valid float.\n", str);
    exit(EXIT_FAILURE);
  }

  return(number);
}

/* Converts a string into a valid double or issue an error message */
long to_valid_long(const char* str)
{
  char** invalid_number = (char**)malloc(sizeof(char*));
  long number;

  *invalid_number = NULL;

  number = strtol(str,invalid_number,10);
  if (**invalid_number != '\0') {
    fprintf(stderr,"(parameter.c) -%s- is not a valid integer.\n", str);
    exit(EXIT_FAILURE);
  }

  return(number);
}
