/* Hello emacs, this is -*- c -*- */
/* $Id: neural.h,v 1.1 2001/07/13 15:48:02 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

/* This module controls Neural Network prototyping and associated 
   functionality. */

#ifndef NEURAL_H_ /* The door keeper */
#define NEURAL_H_

#include "real.h"
#include "event.h"
#include <stdio.h>

/* This is what is a neural network for us: */
/* 1: The neuron */
typedef struct neuron_t {
  real_t* weights; /* These are the sinapses */
  real_t bias; /* The bias of the current neuron */
  real_t u; /* The non-activated sum of all weighted inputs */
  real_t out; /* The output, i.e., u after activation */
} neuron_t;

/* 2: The network */
typedef struct neuralnet_t {
  int input_dim; /* The number of dimensions per input: fixed to all hidden
		    neurons. */
  int nhidden; /* The number of neurons at the hidden layer */
  neuron_t* hidden; /* The hidden layer */
  neuron_t out; /* The _only_ output layer */
} neuralnet_t;

/* The next function will create a neural network based on the contents of the
   given file pointer. */
neuralnet_t* loadnet (FILE*);

/* The next function will free all the resources allocated for the current
   neural network */
void freenet (neuralnet_t*);

/* This function will change the activation function into one of the known
   activation types:
   'f' -> use builtin tanhf();
   'l' -> use a configurable lookup table().
   The lookup table configuration should be performed manually using
   config_lut().
   By default, the builtin tanhf() is used. */
void change_act(const char);

/* The following function will activate the given neuron. The activation itself
   is masked from the user. It can be implemented with a lookup structure or
   applying a common hyperbolic tangent function implementation. */
void activate (neuron_t*);

/* The next function will apply a set of values to the neural network and will
   collect its output. The return value is the output value from the net. The
   arguments are a vector containing the exact number of input patterns
   required and the neural network itself. The input vector, at this process,
   is changed by the input bias of the neural network, be warned! */
real_t propagate (const real_t*, neuralnet_t*);

/* This function takes as argument all the inputs of all classes and process
    them in batch. The first argument is a pointer to the events, the second,
    the network. This function will return all the outputs, placed in a newly
    allocated space. */
outputdata_t* batch_propagate (const inputdata_t*, neuralnet_t*);

/* Evaluate the relevance of each input to network performance. This is done by
   substituting each input by its mean over the whole set. The entry with the
   variable substituted with the mean is subtracted by the real entry and this
   is done for all entries. Such subtractions are raised to the power of 2, in
   order for the signal to disappear and are summed over all entries. The
   result is divided by the number of input patterns and that's the relevance
   of that input. The arguments of this function are: the input patterns, the
   outupt patterns previously evaluated and the neural network itself. */
real_t* relevance (const inputdata_t*, const outputdata_t*, neuralnet_t*);

/* This function will evaluate the _importance_. This will measure how much a
   variable is important to discrimination, more or less in the same sense of
   _relevance_ to MSE. The comparison point is the area between 1 and the
   'class1 efficiency versus class2 error'. Such area will increase more in the
   case a variable is very important to discrimination. Else, the area may
   decrease, in the case the variable is not helping with the discrimination at
   all. The area is calculated using the Simpsons's Rule (1/3). It's a
   quadratic approximation of the curve. */
real_t* importance (const inputdata_t*, const outputdata_t*, neuralnet_t*);
#endif /* NEURAL_H_ */

