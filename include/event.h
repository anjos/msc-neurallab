/* Hello emacs, this is -*- c -*- */
/* $Id: event.h,v 1.1 2001/07/13 15:48:02 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

/* This module defines what is an event for our neural nets. The associated
   functions are also here. */

#ifndef EVENT_H_
#define EVENT_H_

#include "real.h"

typedef struct event_t {
  real_t* data; /* the data */
  real_t target; /* the target for this data */
} event_t;

typedef struct inputdata_t {
  int event_size; /* event size */
  int n1; /* number of events of class 1 */
  event_t* class1;
  int n2; /* number of events on class 2 */
  event_t* class2;
} inputdata_t;

typedef struct outputdata_t {
  int n1; /* number of events of class 1 */
  real_t* class1; /* the class 1 output from the net */
  int n2; /* number of events on class 2 */
  real_t* class2; /* the class 2 output from the net */
  /* This program will not work with localized targets, though the inputdata_t
     has entries for that */
  real_t target1; /* the global target of class 1 */
  real_t target2; /* the global target of class 2 */
} outputdata_t;

/* This function will readin the input data indicated by the filename of
   parameter 1 and will place the results on a newly allocated buffer. This
   procedure will automatically divide the classes of events by the target
   value which is expected to by > 0 for class 1 and < 0 for class 2. The
   second parameter indicates the event size. */
inputdata_t* fread_dat (const char*, const int);

/* This will write the output data into file. The scheme is first the real
   network output, then the target. The parameters are the output filename and
   the output data itself (it won't be modified). If the filename
   (1st. parameter) is "", then the output is redirected to stdout. */
void fwrite_output (char*, const outputdata_t*);

/* This function returns a random event of alternate classes from within the
   parameter 1. It uses rand for that job. An internal (static) variable keeps
   the class number of the last function call, guaranteeing that every other
   class of events will be returned all the time. If you would like this
   fuction to _really_ genearete random sequences, you should initialize its
   seed by issuing at the main routine: srand(time(0)); The function does the
   sorting as explained in Numerical Recipes in C. (the _right_ way) */
const event_t* get_random (const inputdata_t*);

/* This will create a copy the a inputdata_t* and will return it */
inputdata_t* copy_inputdata_t (const inputdata_t*);

/* This liberates the memory used by the input data. */
void free_inputdata_t (inputdata_t*);

/* This liberates the memory used by the output data. */
void free_outputdata_t (outputdata_t*);

#endif /* EVENT_H_ */


