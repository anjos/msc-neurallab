/* Hello emacs, this is -*- c -*- */
/* $Id: event.c,v 1.1 2001/07/13 15:48:02 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

#include "event.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* This function will readin the input data indicated by the filename of
   parameter 1 and will place the results on a newly allocated buffer. This
   procedure will automatically divide the classes of events by the target
   value which is expected to by > 0 for class 1 and < 0 for class 2. The
   second parameter indicates the event size. */
inputdata_t* fread_dat (const char* fname, const int evsize)
{
  int i; /* iterator */
  int end; /* end of file indicator */
  real_t* temp; /* a temporary holder for data */
  real_t ctarget; /* the target of the current event being read-in */
  FILE* fp; /* a file pointer to the newly opened file */
  inputdata_t* d; /* what will be returned */

  d = (inputdata_t*) malloc (sizeof(inputdata_t));

  /* Initialize the structure inputdata_t. Say "bye" if you had anything in
     there... */
  d->n1 = 0;
  d->n2 = 0;
  d->class1 = NULL;
  d->class2 = NULL;
  d->event_size = evsize;

  if ( (fp=fopen(fname,"r")) == NULL ) {
    perror("[fread_dat] System said");
    fprintf(stderr, "[fread_dat] Can't open file %s for reading\n", fname);
    exit(EXIT_FAILURE);
  }

  temp = (real_t*) malloc (d->event_size*sizeof(real_t));

  for(;;) { /* forever */
    for (i=0; i<d->event_size; ++i) {
      end = fscanf(fp, FSCANF_TYPE, &temp[i]);
    }
    end = fscanf(fp, FSCANF_TYPE, &ctarget);
    if (end==EOF) break;

    if (ctarget > 0) {
      d->class1 = (event_t*) realloc (d->class1, (d->n1+1)*sizeof(event_t));
      d->class1[d->n1].data = (real_t*) malloc (d->event_size*sizeof(real_t));
      memcpy(d->class1[d->n1].data, temp, d->event_size*sizeof(real_t));
      d->class1[d->n1].target = ctarget;
      ++d->n1;
    }
    else {
      d->class2 = (event_t*) realloc (d->class2, (d->n2+1)*sizeof(event_t));
      d->class2[d->n2].data = (real_t*) malloc (d->event_size*sizeof(real_t));
      memcpy(d->class2[d->n2].data, temp, d->event_size*sizeof(real_t));
      d->class2[d->n2].target = ctarget;
      ++d->n2;
    }
  }
  
  fclose(fp); /* the file is _not_ needed anymore... */

  fprintf(stderr,
	  "[fread_dat] Read %d events for class 1 and %d events for class 2\n",
	  d->n1, d->n2);

  return d;
}

/* This function returns a random event of alternate classes from within the
   parameter 1. It uses rand for that job. An internal (static) variable keeps
   the class number of the last function call. */
const event_t* get_random (const inputdata_t* d) 
{
  static int last = 2; /* This is our marker of last class called */
  div_t v; /* the event number of this time */

  if (last == 1) { /* get an event of class 2 */
    last = 2;
    v=div(rand(),d->n2);
    return &d->class2[v.rem];
  }
  else { /* get an event of class 1 */
    last = 1;
    v=div(rand(),d->n1);
    return &d->class1[v.rem];
  }

}

/* This liberates the memory used by the input data. */
void free_inputdata_t (inputdata_t* d)
{
  register int i; /* iterator */
  
  for (i=0; i<d->n1; ++i) 
    free(d->class1[i].data);
  for (i=0; i<d->n2; ++i) 
    free(d->class2[i].data);
  free(d->class1);
  free(d->class2);

  free(d); /* finally, we free the allocated resources */

  return;
}

/* This will create a copy the a inputdata_t* and will return it */
inputdata_t* copy_inputdata_t (const inputdata_t* d)
{
  register int i; /* iterator */
  inputdata_t* retval;

  retval = (inputdata_t*) malloc (sizeof(inputdata_t));
  retval->event_size = d->event_size;
  retval->n1 = d->n1;
  retval->n2 = d->n2;
  
  retval->class1 = (event_t*) malloc (d->n1*sizeof(event_t));
  for (i=0; i<d->n1; ++i) {
    retval->class1[i].data = (real_t*) malloc (d->event_size*sizeof(real_t));
    memcpy(retval->class1[i].data, d->class1[i].data,
	   d->event_size*sizeof(real_t));
    retval->class1[i].target = d->class1[i].target;
  }
  retval->class2 = (event_t*) malloc (d->n2*sizeof(event_t));
  for (i=0; i<d->n2; ++i) {
    retval->class2[i].data = (real_t*) malloc (d->event_size*sizeof(real_t));
    memcpy(retval->class2[i].data, d->class2[i].data, 
	   d->event_size*sizeof(real_t));
    retval->class2[i].target = d->class2[i].target;
  }

  return retval;
}
  

/* This liberates the memory used by the output data. */
void free_outputdata_t (outputdata_t* d)
{
  free(d->class1);
  free(d->class2);
  free(d);
  
  return;
}

/* This will write the output data into file. The scheme is first the real
   network output, then the target. The parameters are the output filename and
   the output data itself (it won't be modified). If the filename
   (1st. parameter) is "", then the output is redirected to stdout. */
void fwrite_output (char* fname, const outputdata_t* d)
{
  register int i; /* iterator */
  FILE* fp; /* the input file pointer */

  /* open the input file, if needed */
  if ( strcmp(fname, "") == 0 ) fp = stdout;
  else {
    if ( (fp=fopen(fname,"w")) == NULL ) {
      perror("[fwrite_output] System said");
      fprintf(stderr, "[fwrite_output] Can't open file %s for writing\n", fname);
      exit(EXIT_FAILURE);
    }
  }
  
  for (i=0; i<d->n1; ++i) 
    fprintf(fp, "%e %.1e\n", d->class1[i], d->target1);
  for (i=0; i<d->n2; ++i) 
    fprintf(fp, "%e %.1e\n", d->class2[i], d->target2);
  
  /* close the output file */
  if (fp != stdout) fclose(fp);

  return;
}
