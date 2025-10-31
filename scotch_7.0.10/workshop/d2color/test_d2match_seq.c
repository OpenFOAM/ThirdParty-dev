/* Copyright 2021 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_d2match_seq.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests distance-2 graph      **/
/**                matching algorithms.                    **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 05 jun 2021     **/
/**                                 to   : 05 jun 2021     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "graph.h"

/*
**  The type and structure definitions.
*/

/*+ The work queue. +*/

typedef struct Queue_ {
  Gnum *                    queutab;              /*+ Pointer to queue +*/
  Gnum                      queusiz;              /*+ Size of queue    +*/
  Gnum                      qhedidx;              /*+ Queue head index +*/
  Gnum                      qtalidx;              /*+ Queue tail index +*/
} Queue;

/*
**  The function prototypes.
*/

/**********************/
/*                    */
/* The test routines. */
/*                    */
/**********************/

static
Gnum
graphD2match2 (
const Graph * restrict const  grafptr,
const Gnum * restrict         randtax,
const Gnum                    vertnum)
{
  Gnum                randmax;
  Gnum                veraend;
  Gnum                edgenum;

  const Gnum * restrict const verttax = grafptr->verttax;
  const Gnum * restrict const vendtax = grafptr->vendtax;
  const Gnum * restrict const edgetax = grafptr->edgetax;

  veraend = vertnum;
  randmax = randtax[vertnum];
  for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
    Gnum                vertend;
    Gnum                raveval;

    vertend = edgetax[edgenum];
    raveval = randtax[vertend];
    if ((raveval > randmax) ||
        ((raveval == randmax) && (vertend > veraend))) {
      randmax = raveval;
      veraend = vertend;
    }
  }
  return (veraend);
}

int
graphD2match (
Graph * restrict const      grafptr,
Gnum * restrict const       matetab)
{
  IntRandContext      randdat;
  Gnum * restrict     randtax;
  Gnum * restrict     vmaxtax;
  Gnum * restrict     matetax;
  Queue               queudat;
  Gnum                stepnbr;
  Gnum                vrmnnbr;
  Gnum                vertnum;
  Clock               timedat;                    /* Timing variable */

  const Gnum                  baseval = grafptr->baseval;
  const Gnum                  vertnbr = grafptr->vertnbr;
  const Gnum                  vertnnd = grafptr->vertnnd;
  const Gnum * restrict const verttax = grafptr->verttax;
  const Gnum * restrict const vendtax = grafptr->vendtax;
  const Gnum * restrict const edgetax = grafptr->edgetax;

  if (vertnbr <= 0)                               /* Do not consider empty graphs */
    return (0);

  if (memAllocGroup ((void **) (void *)
                     &queudat.queutab, (size_t) (vertnbr * sizeof (Gnum)),
                     &randtax,         (size_t) (vertnbr * sizeof (Gnum)),
                     &vmaxtax,         (size_t) (vertnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("graphD2match: out of memory");
    return (-1);
  }
  matetax  = matetab - baseval;
  randtax -= baseval;
  vmaxtax -= baseval;
  queudat.queusiz = vertnbr;                      /* Initialize work queue */
  queudat.qhedidx =                               /* Work queue is empty   */
  queudat.qtalidx = 0;

  intRandProc (&randdat, 1);                      /* Use a fixed random seed for the sake of reproducibility */
  intRandSeed (&randdat, 1);

  clockInit  (&timedat);
  clockStart (&timedat);

  for (vertnum = baseval; vertnum < vertnnd; vertnum ++)
    randtax[vertnum] = intRandVal2 (&randdat) & (((Gunum) -1) >> 1); /* Always get positive numbers */

  for (vertnum = baseval; vertnum < vertnnd; vertnum ++) {
    vmaxtax[vertnum] = graphD2match2 (grafptr, randtax, vertnum);

    if (vmaxtax[vertnum] == vertnum) {            /* If vertex is locally dominant     */
      queudat.queutab[queudat.qtalidx] = vertnum; /* Enqueue it for further processing */
      queudat.qtalidx ++;                         /* No overflow possible yet          */
    }
  }

  for (vrmnnbr = vertnbr, stepnbr = 0; ; ) {      /* Outer loop on steps */
    Gnum                qtalidx;                  /* Saved tail index    */
    Gnum                vemanbr;

    vemanbr = 0;
    qtalidx = queudat.qtalidx;                    /* Keep current tail index                */
    while (queudat.qhedidx != qtalidx) {          /* As long as there are enqueued vertices */
      Gnum                vertnum;
      Gnum                edgenum;
      Gnum                vertend;

      vertnum = queudat.queutab[queudat.qhedidx]; /* Dequeue vertex to consider */
      queudat.qhedidx = (queudat.qhedidx + 1) % queudat.queusiz;

      if (randtax[vertnum] < 0)                   /* If vertex already matched, skip it */
        continue;
      if (vmaxtax[vertnum] != vertnum)            /* If vertex not locally dominant, skip it */
        continue;

      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                vertend;

        vertend = edgetax[edgenum];
        if (vmaxtax[vertend] != vertnum)
          break;
      }
      if (edgenum < vendtax[vertnum])             /* If vertex is not distance-2 locally dominant, skip it */
        continue;

      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) { /* Select first unmated neighbor */
        vertend = edgetax[edgenum];
        if (randtax[vertend] >= 0)
          break;
      }

      if (edgenum < vendtax[vertnum]) {           /* If suitable mate found                       */
        matetax[vertend] = vertnum;               /* Mate neighbor mate with vertex               */
        randtax[vertend] = -1;                    /* Do no longer consider neighbor for dominance */
        vemanbr ++;
        queudat.queutab[queudat.qtalidx] = vertend; /* Enqueue mate for further processing */
        queudat.qtalidx = (queudat.qtalidx + 1) % queudat.queusiz;
      }
      else
        vertend = vertnum;
      matetax[vertnum] = vertend;                 /* Mate vertex with its neighbor              */
      randtax[vertnum] = -1;                      /* Do no longer consider vertex for dominance */
      vemanbr ++;
      queudat.queutab[queudat.qtalidx] = vertnum; /* Enqueue vertex for further processing */
      queudat.qtalidx = (queudat.qtalidx + 1) % queudat.queusiz;
    }
    printf ("Step " GNUMSTRING " : " GNUMSTRING "\n", stepnbr, vemanbr);

    stepnbr ++;                                   /* One more step done */
    vrmnnbr -= vemanbr;
    if (vrmnnbr <= 0)                             /* Exit outer loop if nothing more to do */
      break;

    qtalidx = queudat.qtalidx;                    /* Keep current tail index                */
    while (queudat.qhedidx != qtalidx) {          /* As long as there are enqueued vertices */
      Gnum                vertnum;
      Gnum                edgenum;

      vertnum = queudat.queutab[queudat.qhedidx]; /* Dequeue vertex whose neighbors are to recompute */
      queudat.qhedidx = (queudat.qhedidx + 1) % queudat.queusiz;
      vmaxtax[vertnum] = graphD2match2 (grafptr, randtax, vertnum);

      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                vertend;

        vertend = edgetax[edgenum];
        vmaxtax[vertend] = graphD2match2 (grafptr, randtax, vertend);
        queudat.queutab[queudat.qtalidx] = vmaxtax[vertend]; /* Enqueue it as vertex of interest */
        queudat.qtalidx = (queudat.qtalidx + 1) % queudat.queusiz;
      }
    }
  }

  clockStop (&timedat);
  printf    ("Time: %g\n", (double) clockVal (&timedat));

  memFree (queudat.queutab);                      /* Free group leader */

  return (stepnbr);
}

/*********************/
/*                   */
/* The main routine. */
/*                   */
/*********************/

int
main (
int                         argc,
char *                      argv[])
{
  Graph               grafdat;
  Gnum * restrict     matetab;
  FILE *              fileptr;
  Gnum                vertnum;

  if (argc != 3) {
    errorPrint ("usage: %s graph matching", argv[0]);
    exit       (EXIT_FAILURE);
  }

  graphInit (&grafdat);
  if ((fileptr = fopen (argv[1], "r")) == NULL) {
    errorPrint ("main: cannot open graph file");
    exit       (EXIT_FAILURE);
  }
  if (graphLoad (&grafdat, fileptr, -1, 0) != 0) {
    errorPrint ("main: cannot read graph");
    exit       (EXIT_FAILURE);
  }
  fclose (fileptr);

  if ((matetab = memAlloc (grafdat.vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("main: out of memory (1)");
    exit       (EXIT_FAILURE);
  }

  memSet (matetab, ~0, grafdat.vertnbr * sizeof (Gnum)); /* Test for internal errors */

  graphD2match (&grafdat, matetab);

  if ((fileptr = fopen (argv[2], "w")) == NULL) {
    errorPrint ("main: cannot open matching file");
    exit       (EXIT_FAILURE);
  }
  fprintf (fileptr, GNUMSTRING "\n", grafdat.vertnbr);
  for (vertnum = 0; vertnum < grafdat.vertnbr; vertnum ++)
    fprintf (fileptr, GNUMSTRING "\n" GNUMSTRING "\n", vertnum + grafdat.baseval, matetab[vertnum]);
  fclose (fileptr);

  for (vertnum = 0; vertnum < grafdat.vertnbr; vertnum ++) { /* Un-based loop on matetab */
    if ((matetab[vertnum] < 0) || (matetab[vertnum] >= grafdat.vertnbr)) {
      errorPrint ("main: invalid matching");
      exit       (EXIT_FAILURE);
    }
  }

  memFree   (matetab);
  graphExit (&grafdat);

  exit (EXIT_SUCCESS);
}
