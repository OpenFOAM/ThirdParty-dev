/* Copyright 2016,2018 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_fibo.c                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the fibonacci heap    **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 23 aug 2016     **/
/**                                 to     22 may 2018     **/
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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../libscotch/module.h"
#include "../libscotch/common.h"
#include "../libscotch/fibo.h"

/*
**  The type and structure definitions.
*/

/* The fibo cell structure. */

typedef struct TestFibo_ {
  FiboNode                  nodedat;              /* TRICK: FIRST */
  int                       randval;
} TestFibo;

/***************************/
/*                         */
/* The comparison routine. */
/*                         */
/***************************/

static
int
testFiboCmpFunc (
const FiboNode *            nod0ptr,
const FiboNode *            nod1ptr)
{
  int                 ran0val;
  int                 ran1val;

  ran0val = ((TestFibo *) nod0ptr)->randval;
  ran1val = ((TestFibo *) nod1ptr)->randval;
  if (ran0val < ran1val)
    return (-1);
  return (1);
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
  FiboHeap            fibodat;
  TestFibo *          nodetab;
  TestFibo *          nodeptr;
  int                 nodesiz;
  int                 nodemax;
  int                 nodenbr;
  int                 nodenum;
  int                 nodetmp;
  int                 randval;
  int                 randtmp;
  int                 passnbr;
  int                 passnum;

  SCOTCH_errorProg (argv[0]);

  if (fiboHeapInit (&fibodat, testFiboCmpFunc) != 0) {
    SCOTCH_errorPrint ("main: cannot initialize Fibonacci heap");
    exit (EXIT_FAILURE);
  }

  intRandInit ();                                 /* Initialize random generator */

  nodesiz = 100;
  passnbr = -1;
  switch (argc) {
    case 4 :
      intRandSeed (MAX (0, atoi (argv[3])));
    case 3 :
      passnbr = MAX (1, atoi (argv[2]));
    case 2 :
      nodesiz = MAX (1, atoi (argv[1]));
    case 1 :
      break;
    default :
      SCOTCH_errorPrint ("usage: %s [nodenbr [passnbr [seed]]]", argv[0]);
      exit (EXIT_FAILURE);
  }
  if (passnbr < 0)
    passnbr = 10 * nodesiz;

  if ((nodetab = malloc (nodesiz * sizeof (TestFibo))) == NULL) {
    SCOTCH_errorPrint ("main: out of memory");
    exit (EXIT_FAILURE);
  }
  for (nodenum = 0; nodenum < nodesiz; nodenum ++) /* Initialize node array */
    nodetab[nodenum].randval = -1;

  nodemax = nodesiz - 1;                          /* Set maximum index */
  nodenbr = 0;                                    /* Array is empty    */
  for (passnum = 0; passnum < passnbr; passnum ++) {
    switch (intRandVal (6)) {
      case 0 :                                    /* Add node */
      case 1 :
      case 2 :                                    /* More additions than deletions on average */
        if (nodenbr >= nodemax)
          break;
        for (nodenum = 0; nodenum <= nodemax; nodenum ++) { /* Search for a free slot */
          if (nodetab[nodenum].randval < 0)
            break;
        }
        if (nodenum > nodemax) {
          SCOTCH_errorPrint ("main: invalid node array (1)");
          exit (EXIT_FAILURE);
        }
        nodetab[nodenum].randval = abs (intRandVal (INTVALMAX));
        fiboHeapAdd (&fibodat, (FiboNode *) &nodetab[nodenum]);
        nodenbr ++;
        break;
      case 3 :                                    /* Remove arbitrary node */
        if (nodenbr <= 0)
          break;
        nodetmp = intRandVal (nodenbr);
        for (nodenum = 0; ; nodenum ++) {         /* Search for non-empty slot */
          if (nodenum > nodemax) {
            SCOTCH_errorPrint ("main: invalid node array (2)");
            exit (EXIT_FAILURE);
          }
          if (nodetab[nodenum].randval >= 0) {
            if (-- nodetmp <= 0)
              break;
          }
        }
        fiboHeapDel (&fibodat, (FiboNode *) &nodetab[nodenum]);
        nodetab[nodenum].randval = -1;
        nodenbr --;
        break;
      case 4 :                                    /* Remove minimum node */
        if (nodenbr <= 0)
          break;
        nodeptr = (TestFibo *) fiboHeapMin (&fibodat);
        randval = nodeptr->randval;               /* Keep node key value   */
        fiboHeapDel (&fibodat, (FiboNode *) nodeptr); /* Remove node first */
        nodeptr->randval = -1;
        nodenbr --;
        for (nodenum = 0; nodenum <= nodemax; nodenum ++) { /* Check if smaller node exists */
          if ((nodetab[nodenum].randval >= 0) &&
              (nodetab[nodenum].randval <  randval)) {
            SCOTCH_errorPrint ("main: node is not of minimum key");
          }
        }
        break;
      case 5 :                                    /* Decrease value of arbitrary node */
        if (nodenbr <= 0)
          break;
        nodetmp = intRandVal (nodenbr);
        for (nodenum = 0; ; nodenum ++) {         /* Search for non-empty slot */
          if (nodenum > nodemax) {
            SCOTCH_errorPrint ("main: invalid node array (3)");
            exit (EXIT_FAILURE);
          }
          if (nodetab[nodenum].randval >= 0) {
            if (-- nodetmp <= 0)
              break;
          }
        }
        if (nodetab[nodenum].randval <= 0)        /* Cannot decrease smallest value */
          break;
        randtmp = intRandVal (nodetab[nodenum].randval + 1);
        if (randtmp > nodetab[nodenum].randval)
          break;
        nodetab[nodenum].randval -= randtmp;
        fiboHeapDecrease (&fibodat, (FiboNode *) &nodetab[nodenum]);
        break;
    }
  }

  fiboHeapExit (&fibodat);
  free         (nodetab);

  exit (EXIT_SUCCESS);
}
