/* Copyright 2019,2021,2024 IPB, Universite de Bordeaux, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
**
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
**
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
**
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : test_scotch_context.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Clement BARTHELEMY                      **/
/**                                                        **/
/**   FUNCTION   : This module tests the thread import     **/
/**                feature of the library Context object.  **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 25 aug 2019     **/
/**                                 to   : 30 jul 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE               600
#endif /* _XOPEN_SOURCE */
#ifndef __USE_XOPEN2K
#define __USE_XOPEN2K                             /* For POSIX pthread_barrier_t */
#endif /* __USE_XOPEN2K */

#include <stdio.h>
#if (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H))
#include <stdint.h>
#endif /* (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H)) */
#include <stdlib.h>
#include <string.h>

#include "../libscotch/module.h"
#include "../libscotch/common.h"
#include "../libscotch/context.h"
#include "../libscotch/graph.h"
#include "scotch.h"

/*
**  The type and structure definitions.
*/

/*+ The thread-specific data block that represents what a thread has access to. +*/

typedef struct TestThreadData_ {
  SCOTCH_Context *          contptr;              /* Pointer to the context */
  int                       thrdnum;
} TestThreadData;


/**********************************/
/*                                */
/* The fake threaded API routine. */
/*                                */
/**********************************/

/* This routine checks that a context is proprely
** bound to a graph and tests it.
** It returns:
** - 0   : if the binding succeeded.
** - !0  : on error.
*/

static
void
scotchGraphDoUsefulStuff (
ThreadDescriptor * restrict const descptr,
const Graph * restrict const      grafptr)
{
  const int                 thrdnbr = threadNbr (descptr);
  const int                 thrdnum = threadNum (descptr);

  printf ("(%d/%d) : " GNUMSTRING "\n", thrdnum, thrdnbr, grafptr->vertnbr);
}

/* Split work routine. */

void
scotchSplit (
Context * const             contptr,              /* Sub-context                  */
const int                   spltnum,              /* Rank of sub-context (0 or 1) */
Graph * restrict const      grafptr)
{
  const int                 thrdnbr = contextThreadNbr (contptr);

  printf ("Sub-context %d (%d) : " GNUMSTRING ", " GNUMSTRING "\n", spltnum, thrdnbr, grafptr->vertnbr, contextIntRandVal (contptr, grafptr->vertnbr));

  printf ("Launching multi-threaded work in sub-context\n");

  contextThreadLaunch (contptr,
                       (ThreadFunc) scotchGraphDoUsefulStuff,
                       grafptr);
}

/* Context-dependent main task routine. */

int
SCOTCH_graphDoUsefulStuff (
const SCOTCH_Graph * restrict const libgrafptr)
{
  CONTEXTDECL        (libgrafptr);

  if (CONTEXTINIT (libgrafptr)) {
    errorPrint (STRINGIFY (SCOTCH_contextGraphTest) ": cannot initialize context");
    return     (1);
  }

  printf ("Launching multi-threaded work in provided context\n");

  contextThreadLaunch (CONTEXTGETDATA (libgrafptr), /* Fake libScotch routine uses worker threads to perform some task in parallel */
                       (ThreadFunc) scotchGraphDoUsefulStuff,
                       (void *) CONTEXTGETOBJECT (libgrafptr));

  printf ("Splitting provided context\n");

  contextThreadLaunchSplit (CONTEXTGETDATA (libgrafptr),
                            (ContextSplitFunc) scotchSplit,
                            (void *) CONTEXTGETOBJECT (libgrafptr));

  printf ("Back to provided context\n");

  CONTEXTEXIT (libgrafptr);

  return (0);
}

/********************************/
/*                              */
/* The threaded import routine. */
/*                              */
/********************************/

#ifdef COMMON_PTHREAD

static
void *
testContextImport (
TestThreadData *            thrdptr)
{
  pthread_t           thidval;

  thidval = pthread_self ();
  pthread_detach (thidval);                       /* Pretend not to care about the worker threads */

  SCOTCH_contextThreadImport2 (thrdptr->contptr, thrdptr->thrdnum); /* All threads should call contextThreadImport2(); worker threads are captured */

  return (NULL);                                  /* Worker threads live their own life after the context is freed */
}

#endif /* COMMON_PTHREAD */

/*********************/
/*                   */
/* The main routine. */
/*                   */
/*********************/

int
main (
int                 argc,
char *              argv[])
{
  FILE *              fileptr;
  SCOTCH_Context      contdat;
  TestThreadData *    thrdtab;                    /* Array of thread local data           */
  int                 thrdnbr;
  int                 thrdnum;
  SCOTCH_Graph        graftab[2];                 /* Original graph and context container */
  int *               coretab;

  SCOTCH_errorProg (argv[0]);

  if (argc != 2) {
    SCOTCH_errorPrint ("usage: %s graph_file", argv[0]);
    exit (EXIT_FAILURE);
  }

#ifdef SCOTCH_PTHREAD_NUMBER
  thrdnbr = SCOTCH_PTHREAD_NUMBER;                /* If prescribed number defined at compile time, use it as default */
#else /* SCOTCH_PTHREAD_NUMBER */
  thrdnbr = 3;
#endif /* SCOTCH_PTHREAD_NUMBER */

  SCOTCH_graphInit (&graftab[0]);

  if ((fileptr = fopen (argv[1], "r")) == NULL) {
    SCOTCH_errorPrint ("main: cannot open file");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphLoad (&graftab[0], fileptr, -1, 0) != 0) { /* Read source graph */
    SCOTCH_errorPrint ("main: cannot load graph");
    exit (EXIT_FAILURE);
  }

  fclose (fileptr);

  printf ("Using a default context\n");

  if (SCOTCH_contextInit (&contdat) != 0) {
    SCOTCH_errorPrint ("main: cannot initialize context");
    exit (EXIT_FAILURE);
  }

  SCOTCH_graphInit (&graftab[1]);
  if (SCOTCH_contextBindGraph (&contdat, &graftab[0], &graftab[1]) != 0) { /* graftab[1] is the context graph */
    SCOTCH_errorPrint ("main: cannot bind context (1)");
    exit (EXIT_FAILURE);
  }

  SCOTCH_graphDoUsefulStuff (&graftab[1]);        /* Do some useful work using the context graph */

  SCOTCH_graphExit   (&graftab[1]);               /* Free the context graph before its bound context */
  SCOTCH_contextExit (&contdat);

  printf ("Using a tailored context (1)\n");

  SCOTCH_contextInit (&contdat);
  SCOTCH_contextThreadSpawn (&contdat, thrdnbr, NULL);

  SCOTCH_graphInit (&graftab[1]);
  if (SCOTCH_contextBindGraph (&contdat, &graftab[0], &graftab[1]) != 0) { /* graftab[1] is the context graph */
    SCOTCH_errorPrint ("main: cannot bind context (2)");
    exit (EXIT_FAILURE);
  }

  SCOTCH_graphDoUsefulStuff (&graftab[1]);        /* Do some useful work using the context graph */

  SCOTCH_graphExit   (&graftab[1]);
  SCOTCH_contextExit (&contdat);

  printf ("Using a tailored context (2)\n");

  if ((coretab = (int *) malloc (thrdnbr * sizeof (int))) == NULL) { /* Create core placement array */
    SCOTCH_errorPrint ("main: out of memory (1)");
    exit (EXIT_FAILURE);
  }
  coretab[0] = 0;
  for (thrdnum = 1; thrdnum < thrdnbr; thrdnum ++)
    coretab[thrdnum] = thrdnbr - thrdnum;

  SCOTCH_contextInit (&contdat);
  SCOTCH_contextThreadSpawn (&contdat, thrdnbr, coretab);
  if (SCOTCH_contextRandomClone (&contdat) != 0) {
    SCOTCH_errorPrint ("main: cannot clone random context");
    exit (EXIT_FAILURE);
  }

  free (coretab);

  SCOTCH_graphInit (&graftab[1]);
  if (SCOTCH_contextBindGraph (&contdat, &graftab[0], &graftab[1]) != 0) { /* graftab[1] is the context graph */
    SCOTCH_errorPrint ("main: cannot bind context (3)");
    exit (EXIT_FAILURE);
  }

  SCOTCH_graphDoUsefulStuff (&graftab[1]);        /* Do some useful work using the context graph */

  SCOTCH_graphExit   (&graftab[1]);
  SCOTCH_contextExit (&contdat);

  printf ("Using an imported context\n");

  if ((thrdtab = (TestThreadData *) malloc (thrdnbr * sizeof (TestThreadData))) == NULL) { /* Create a fake environment for threads */
    SCOTCH_errorPrint ("main: out of memory (2)");
    exit (EXIT_FAILURE);
  }
  for (thrdnum = 0; thrdnum < thrdnbr; thrdnum ++) {
    thrdtab[thrdnum].contptr = &contdat;          /* All threads know where the context is              */
    thrdtab[thrdnum].thrdnum = thrdnum;           /* All threads know their index in their thread group */
  }

  SCOTCH_contextInit (&contdat);                  /* Leader thread creates the context                */
  SCOTCH_contextThreadImport1 (&contdat, thrdnbr); /* Leader thread calls "Import1" before any others */

#ifdef COMMON_PTHREAD
  for (thrdnum = 1; thrdnum < thrdnbr; thrdnum ++) { /* Then all other threads run their part */
    pthread_t           thidval;

    if (pthread_create (&thidval, NULL, (void * (*) (void *)) testContextImport, (void *) &thrdtab[thrdnum]) != 0) {
      SCOTCH_errorPrint ("main: cannot launch thread (%d)", thrdnum);
      exit (EXIT_FAILURE);
    }
  }
#endif /* COMMON_PTHREAD */

  SCOTCH_contextThreadImport2 (&contdat, 0);      /* Leader thread calls "Import2" along with all other threads */

  SCOTCH_graphInit (&graftab[1]);
  if (SCOTCH_contextBindGraph (&contdat, &graftab[0], &graftab[1]) != 0) { /* graftab[1] is the context graph */
    SCOTCH_errorPrint ("main: cannot bind context (4)");
    exit (EXIT_FAILURE);
  }

  SCOTCH_graphDoUsefulStuff (&graftab[1]);        /* Do some useful work using the context graph */

  SCOTCH_graphExit   (&graftab[1]);
  SCOTCH_contextExit (&contdat);                  /* Worker threads are released when the context is destroyed */

  free (thrdtab);

  SCOTCH_graphExit (&graftab[0]);

  exit (EXIT_SUCCESS);
}
