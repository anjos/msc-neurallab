/* Hello emacs, this is -*- c -*- */
/* $Id: neural.c,v 1.1 2001/07/13 15:48:02 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "neural.h"
#include "lut.h"
#include "real.h"
#include "event.h"
#include "analysis.h"

/* The next declaration will mark the current activation type. In order to
   change it, use change_act() */
static real_t (*global_act) (real_t) = &FUNC_TANH;

/* PROTOTYPES */

/* This will calculate the means of each variable over the patterns hold in
   `input' (1st. argument). The results will be stored over the _pre-allocated_
   space pointed by the second argument. This is a fair process and doesn't
   take in account the number of patterns in each class.
   It's accomplished by extracting the means of each class separately and then
   summing and dividing both by 2. */
real_t* mean (const inputdata_t*, real_t*);

/* IMPLEMENTATION */

void change_act(const char c)
{
  switch (c) {
  case 'f':
    global_act = &FUNC_TANH;
    fprintf(stderr, "[neural] Activation Function changed to BUILTIN\n");
    break;
  case 'l':
    global_act = &lookup;
    fprintf(stderr, "[neural] Activation Function changed to LOOKUP\n");
    break;
  default:
    global_act = &FUNC_TANH;
    fprintf(stderr, "[neural] Activation Function changed to BUILTIN\n");
  }
  return;
}

neuralnet_t* loadnet (FILE* fp)
{
  int i,j; /* iterator */
  neuralnet_t* net;

  net = (neuralnet_t*) calloc (1, sizeof(neuralnet_t));
  if (net == NULL) {
    perror("[neural] net");
    exit(EXIT_FAILURE);
  }

  /* first, read the dimensionality of the input patterns */
  fscanf(fp, "%d", &net->input_dim);

  /* Read the number of hidden neurons */
  fscanf(fp, "%d", &net->nhidden);
  
  /* Now we configure the output neuron */
  net->out.weights = (real_t*) calloc (net->nhidden, sizeof(real_t));
  if (net->out.weights == NULL) {
    perror("[neural] net->out.weights");
    exit(EXIT_FAILURE);
  }

  /* And configure also the hidden layers */
  net->hidden = (neuron_t*) calloc (net->nhidden, sizeof(neuron_t));
  for (i=0; i < net->nhidden; ++i) {
    net->hidden[i].weights = (real_t*) calloc (net->input_dim, sizeof(real_t));
    if (net->out.weights == NULL) {
      perror("[neural] net->hidden.weights");
      exit(EXIT_FAILURE);
    }
  }

  /* Now we have to read the weight values for each synapse on the way. This
     means all the weights and biases for the net->nhidden neurons and also for
     the output layer. It's expected that the file contain a matrix of floats
     arranged rowwise with the weights and bias (respectively) of all such
     neurons. The last weights (and a bias) are taken to be the ones that
     connect the hidden layer to the output layer. */
  for (i=0; i < net->nhidden; ++i) {
    for (j=0; j< net->input_dim; ++j)
      fscanf(fp, FSCANF_TYPE, &net->hidden[i].weights[j]);
    /* read the bias for this neuron */
    fscanf(fp, FSCANF_TYPE, &net->hidden[i].bias); 
  }

  for (i=0; i < net->nhidden; ++i)
    fscanf(fp, FSCANF_TYPE, &net->out.weights[i]);
  /* read the bias for this neuron */
  fscanf(fp, FSCANF_TYPE, &net->out.bias);

  fprintf(stderr, "[neural] Network %d-%d-1 loaded from file\n",
	  net->input_dim, net->nhidden);

  /* Now we can safely finish */
  return net;
}

void freenet (neuralnet_t* net)
{
  int i;
  for(i=0; i<net->nhidden; ++i) free(net->hidden[i].weights);
  free(net->hidden);
  free(net->out.weights);
  free(net);
}

void activate (neuron_t* n)
{
  n->out = global_act(n->u);
}

real_t propagate (const real_t* pat, neuralnet_t* net)
{
  register int i; /* iterators */
  register int j;

  /* We start with the hidden layers */
  net->out.u = net->out.bias;
  for (i=0; i<net->nhidden; ++i) {
    net->hidden[i].u = net->hidden[i].bias;
    for (j=0; j<net->input_dim; ++j) {
      net->hidden[i].u += net->hidden[i].weights[j] * pat[j];
    }
    activate(&net->hidden[i]);
    net->out.u += net->hidden[i].out * net->out.weights[i];
  }

  activate(&net->out);

  return net->out.out;
}

/* This function takes as argument all the inputs of all classes and process
    them in batch. The first argument is a pointer to the events, the second,
    the network. This function will return all the outputs, placed in a newly
    allocated space. */
outputdata_t* batch_propagate (const inputdata_t* input, neuralnet_t* net)
{
  register int i; /* iterator */
  outputdata_t* retval;
  
  /* Initialize the output data */
  retval = (outputdata_t*) malloc (sizeof(outputdata_t));
  retval->n1 = input->n1;
  retval->n2 = input->n2;
  retval->target1 = input->class1[0].target;
  retval->target2 = input->class2[0].target;
  
  /* Get the output */
  retval->class1 = (real_t*) malloc (input->n1 * sizeof(real_t));
  for (i=0; i<input->n1; ++i) 
    retval->class1[i] = propagate(input->class1[i].data, net);

  retval->class2 = (real_t*) malloc (input->n2 * sizeof(real_t));
  for (i=0; i<input->n2; ++i) 
    retval->class2[i] = propagate(input->class2[i].data, net);

  /* Return the newly allocated data */
  return retval;
}

/* Evaluate the relevance of each input to network performance. This is done by
   substituting each input by its mean over the whole set. The entry with the
   variable substituted with the mean is subtracted by the real entry and this
   is done for all entries. Such subtractions are raised to the power of 2, in
   order for the signal to disappear and are summed over all entries. The
   result is divided by the number of input patterns and that's the relevance
   of that input. The arguments of this function are: the input patterns, the
   outupt patterns previously evaluated and the neural network itself. */
real_t* relevance (const inputdata_t* input, const outputdata_t* output, 
		   neuralnet_t* net)
{
  register int i,j; /* iterators */
  real_t changed; /* the output with the value substituted by the mean */

  real_t* m; /* These are for the means */
  real_t* data; /* and these, a placeholder for the current event */
  real_t* retval; /* this is what we will return to the caller */

  m = (real_t*) alloca (net->input_dim*sizeof(real_t));
  data = (real_t*) alloca (net->input_dim*sizeof(real_t));

  retval = (real_t*) calloc (net->input_dim, sizeof(real_t));

  /* First calculate the means. This is a fair process and doesn't depend on
     the number of patterns in each class. */
  m = mean(input,m);

  /* now we calculate the relevance for all classes. */
  for (i=0; i<net->input_dim; ++i) {
    for(j=0; j<input->n1; ++j) {
      memcpy(data, input->class1[j].data, net->input_dim*sizeof(real_t));
      data[i] = m[i]; /* substitute the value by its mean */
      changed = propagate(data, net);
      retval[i] += FUNC_POW(output->class1[j]-changed, 2);
    }
    for(j=0; j<input->n2; ++j) {
      memcpy(data, input->class2[j].data, net->input_dim*sizeof(real_t));
      data[i] = m[i]; /* substitute the value by its mean */
      changed = propagate(data, net);
      retval[i] += FUNC_POW(output->class2[j]-changed, 2);
    }
    retval[i] /= (input->n1 + input->n2);
  }

  return retval;
}

/* This function will evaluate the _importance_. This will measure how much a
   variable is important to discrimination, more or less in the same sense of
   _relevance_ to MSE. The comparison point is the area between 1 and the
   'class1 efficiency versus class2 error'. Such area will increase more in the
   case a variable is very important to discrimination. Else, the area may
   decrease, in the case the variable is not helping with the discrimination at
   all. The area is calculated using the Simpsons's Rule (1/3). It's a
   quadratic approximation of the curve. */
real_t* importance (const inputdata_t* input, const outputdata_t* output, 
		    neuralnet_t* net)
{ 
  register int i,j; /* iterators */
  analysis_t an, changed; /* the analysis structures */
  real_t* m; /* These are for the means */
  real_t* retval; /* this is what we will return to the caller */
  analysis_parameter_t anpar; /* Do not eval SP, but do AREA */

  anpar.eval_sp = FALSE;
  anpar.eval_area = TRUE;

  m = (real_t*) alloca (input->event_size*sizeof(real_t));
  retval = (real_t*) calloc (input->event_size, sizeof(real_t));

  /* First calculate the means. This is a fair process and doesn't depend on
     the number of patterns in each class. */
  m = mean(input,m);

  /* We firstly calculate the area between 1 and the efficiency curve for our
     current discriminator. */
  analyze(&an, output, &anpar);

  /* now we calculate the importance. This is done by replicating the
     inputdata_t variable and substituting each value by its mean. */
  for (i=0; i<input->event_size; ++i) { /* for all input positions */
      inputdata_t* newinput = copy_inputdata_t (input);
      outputdata_t* newoutput;
      for (j=0; j<input->n1; ++j)
	newinput->class1[j].data[i] = m[i]; /* substitute the mean */
      for (j=0; j<input->n2; ++j)
	newinput->class2[j].data[i] = m[i]; /* substitute the mean */

      /* re-analyze the data */
      newoutput = batch_propagate(newinput, net);
      analyze(&changed, newoutput, &anpar);

      /* Evaluates the difference. Note that if the difference area increases,
	 it means that the variable was important to discrimination, else, the
	 variable is indifferent (0) or is confusing the discrimination process
	 (<0). */
      retval[i] = changed.area - an.area;

      /* Free used data */
      free_inputdata_t(newinput);
      free_outputdata_t(newoutput);
  }
      
  return retval;  
} 

/* This will calculate the means of each variable over the patterns hold in
   `input' (1st. argument). The results will be stored over the _pre-allocated_
   space pointed by the second argument. This is a fair process and doesn't
   take in account the number of patterns in each class.
   It's accomplished by extracting the means of each class separately and then
   summing and dividing both by 2. */
real_t* mean (const inputdata_t* input, real_t* d)
{
  register int i, j; /* iterators */
  real_t cm1; /* The current mean for class 1 */
  real_t cm2; /* The current mean for class 2 */

  for (i=0; i<input->event_size; ++i) {
    cm1 = 0;
    cm2 = 0;
    for (j=0; j<input->n1; ++j) cm1 += input->class1[j].data[i];
    cm1 /= input->n1;
    for (j=0; j<input->n2; ++j) cm2 += input->class2[j].data[i];
    cm2 /= input->n2;
    d[i] = (cm1 + cm2)/2;
  }

  return d;
}
