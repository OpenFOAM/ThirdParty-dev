/* Copyright 2012,2014,2015,2018,2019,2021 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_common_thread.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the thread            **/
/**                management module.                      **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 04 nov 2012     **/
/**                                 to   : 10 jul 2018     **/
/**                # Version 7.0  : from : 21 aug 2019     **/
/**                                 to   : 31 aug 2021     **/
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

#include "../libscotch/module.h"
#include "../libscotch/common.h"
#include "../libscotch/common_thread.h"
#include "../libscotch/common_thread_system.h"

#define COMPVAL(n)                  (((n) * ((n) + 1)) / 2)

/*
**  The static and global variables.
*/

static int                  C_erroval = 0;        /* Global error value */

/*
**  The type and structure definitions.
*/

/*+ The thread-specific data block. +*/

typedef struct TestData_ {
  int                       reduval;              /*+ Value to reduce +*/
  int                       scanval[2];           /*+ Values for scan +*/
} TestData;

/*+ The block data structure +*/

typedef struct TestGroup_ {
  TestData *                datatab;              /*+ Thread test data              +*/
  int                       redusum;              /*+ Value to compare reduction to +*/
} TestGroup;

/*************************/
/*                       */
/* The threaded routine. */
/*                       */
/*************************/

static
void
testReduce (
TestData * restrict const   vlocptr,              /* Pointer to local value  */
TestData * restrict const   vremptr,              /* Pointer to remote value */
void *                      dataptr)              /* Pointer to shared data  */
{
  if (dataptr != &C_erroval)
    errorPrint ("testReduce: invalid data pointer");

  vlocptr->reduval += vremptr->reduval;
}

static
void
testScan (
TestData * restrict const   vlocptr,              /* Pointer to local value  */
TestData * restrict const   vremptr,              /* Pointer to remote value */
const int                   srcpval,              /* Source phase value      */
const int                   dstpval,              /* Destination phase value */
void *                      dataptr)              /* Pointer to shared data  */
{
  if (dataptr != &C_erroval)
    errorPrint ("testScan: invalid data pointer");

  vlocptr->scanval[dstpval] = vlocptr->scanval[srcpval] + ((vremptr == NULL) ? 0 : vremptr->scanval[srcpval]);
}

static
void
testThreads (
ThreadDescriptor * restrict const descptr,
volatile TestGroup * restrict     grouptr)
{
  const int           thrdnum = threadNum (descptr);
  TestData * const    dataptr = &grouptr->datatab[thrdnum];

  printf ("%d: running\n", thrdnum);

  threadBarrier (descptr);
  fflush (stdout);

  if (thrdnum == 0)
    printf ("Performing reduction\n");

  dataptr->reduval = thrdnum + 1;
  threadReduce (descptr, (void *) dataptr, sizeof (TestData), (ThreadReduceFunc) testReduce, 0, &C_erroval);

  if ((thrdnum == 0) &&                           /* Test reduction result on thread 0 */
      (dataptr->reduval != grouptr->redusum)) {
    errorPrint ("testThreads: invalid reduction operator (0)");
    C_erroval = 1;
  }

  threadBarrier (descptr);

  if (thrdnum == 0)
    printf ("Performing scan\n");

  dataptr->scanval[0] = 1 + thrdnum;
  threadScan (descptr, (void *) dataptr, sizeof (TestData), (ThreadScanFunc) testScan, &C_erroval);

  if (dataptr->scanval[0] != COMPVAL (thrdnum + 1)) {
    errorPrint ("testThreads: invalid scan operator (%d)", thrdnum);
    C_erroval = 1;
  }

  threadBarrier (descptr);                        /* Final barrier before freeing work array */
}

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
  ThreadContext       contdat;
  TestGroup           groudat;
  int                 thrdnbr;

  errorProg (argv[0]);

#ifdef SCOTCH_PTHREAD_NUMBER
  thrdnbr = SCOTCH_PTHREAD_NUMBER;                /* If prescribed number defined at compile time, use it as default */
#else /* SCOTCH_PTHREAD_NUMBER */
  thrdnbr = -1;                                   /* Else take the number of cores at run time */
#endif /* SCOTCH_PTHREAD_NUMBER */
  thrdnbr = envGetInt ("SCOTCH_PTHREAD_NUMBER", thrdnbr);
  if (thrdnbr < 1)
    thrdnbr = threadSystemCoreNbr ();

  if (threadContextInit (&contdat, thrdnbr, NULL) != 0) {
    errorPrint ("main: cannot initialize thread context");
    exit       (EXIT_FAILURE);
  }

  thrdnbr = threadContextNbr (&contdat);
  printf ("%d threads in context\n", thrdnbr);

  if ((groudat.datatab = malloc (thrdnbr * sizeof (TestData))) == NULL) {
    errorPrint ("main: out of memory");
    exit       (EXIT_FAILURE);
  }
  groudat.redusum = COMPVAL (thrdnbr);

  threadLaunch (&contdat, (ThreadFunc) testThreads, (void *) &groudat);

  free (groudat.datatab);

  threadContextExit (&contdat);

  exit ((C_erroval == 0) ? EXIT_SUCCESS : EXIT_FAILURE);
}
