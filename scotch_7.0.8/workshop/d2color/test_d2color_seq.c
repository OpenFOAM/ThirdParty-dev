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
/**   NAME       : test_d2color_seq.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests distance-2 graph      **/
/**                coloring algorithms.                    **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 01 jun 2021     **/
/**                                 to   : 03 jun 2021     **/
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
graphD2color2 (
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
graphD2color (
Graph * restrict const      grafptr,
Gnum * restrict const       colotab)
{
  IntRandContext      randdat;
  Gnum * restrict     randtax;
  Gnum * restrict     vmaxtax;
  Gnum * restrict     colotax;
  Queue               queudat;
  Gnum                colonbr;
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
    errorPrint ("graphD2color: out of memory");
    return (-1);
  }
  colotax  = colotab - baseval;
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
    vmaxtax[vertnum] = graphD2color2 (grafptr, randtax, vertnum);

    if (vmaxtax[vertnum] == vertnum) {            /* If vertex is locally dominant     */
      queudat.queutab[queudat.qtalidx] = vertnum; /* Enqueue it for further processing */
      queudat.qtalidx ++;                         /* No overflow possible yet          */
    }
  }

  for (vrmnnbr = vertnbr, colonbr = 0; ; ) {      /* Outer loop on colors */
    Gnum                qtalidx;                  /* Saved tail index     */
    Gnum                veconbr;

    veconbr = 0;
    qtalidx = queudat.qtalidx;                    /* Keep current tail index                */
    while (queudat.qhedidx != qtalidx) {          /* As long as there are enqueued vertices */
      Gnum                vertnum;
      Gnum                edgenum;

      vertnum = queudat.queutab[queudat.qhedidx]; /* Dequeue vertex to consider */
      queudat.qhedidx = (queudat.qhedidx + 1) % queudat.queusiz;

      if (randtax[vertnum] < 0)                   /* If vertex already colored, skip it */
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

      colotax[vertnum] = colonbr;                 /* Set vertex color                           */
      randtax[vertnum] = -1;                      /* Do no longer consider vertex for dominance */
      veconbr ++;
      queudat.queutab[queudat.qtalidx] = vertnum; /* Enqueue it for further processing */
      queudat.qtalidx = (queudat.qtalidx + 1) % queudat.queusiz;
    }
    printf ("Color " GNUMSTRING " : " GNUMSTRING "\n", colonbr, veconbr);

    colonbr ++;                                   /* One more color used */
    vrmnnbr -= veconbr;
    if (vrmnnbr <= 0)                             /* Exit outer loop if nothing more to do */
      break;

    qtalidx = queudat.qtalidx;                    /* Keep current tail index                */
    while (queudat.qhedidx != qtalidx) {          /* As long as there are enqueued vertices */
      Gnum                vertnum;
      Gnum                edgenum;

      vertnum = queudat.queutab[queudat.qhedidx]; /* Dequeue vertex whose neighbors are to recompute */
      queudat.qhedidx = (queudat.qhedidx + 1) % queudat.queusiz;
      vmaxtax[vertnum] = graphD2color2 (grafptr, randtax, vertnum);

      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                vertend;

        vertend = edgetax[edgenum];
        vmaxtax[vertend] = graphD2color2 (grafptr, randtax, vertend);
        queudat.queutab[queudat.qtalidx] = vmaxtax[vertend]; /* Enqueue it as vertex of interest */
        queudat.qtalidx = (queudat.qtalidx + 1) % queudat.queusiz;
      }
    }
  }

  clockStop (&timedat);
  printf    ("Time: %g\n", (double) clockVal (&timedat));

  memFree (queudat.queutab);                      /* Free group leader */

  return (colonbr);
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
  Gnum * restrict     colotab;
  Gnum                colonbr;
  Gnum                colonum;
  Gnum * restrict     chektax;
  FILE *              fileptr;
  Gnum                vertnum;

  if (argc != 3) {
    errorPrint ("usage: %s graph mapping", argv[0]);
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

  if ((colotab = memAlloc (grafdat.vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("main: out of memory (1)");
    exit       (EXIT_FAILURE);
  }

  memSet (colotab, ~0, grafdat.vertnbr * sizeof (Gnum)); /* Test for internal errors */

  colonbr = graphD2color (&grafdat, colotab);

  if ((fileptr = fopen (argv[2], "w")) == NULL) {
    errorPrint ("main: cannot open mapping file");
    exit       (EXIT_FAILURE);
  }
  fprintf (fileptr, GNUMSTRING "\n", grafdat.vertnbr);
  for (vertnum = 0; vertnum < grafdat.vertnbr; vertnum ++)
    fprintf (fileptr, GNUMSTRING "\n" GNUMSTRING "\n", vertnum + grafdat.baseval, colotab[vertnum]);
  fclose (fileptr);

  for (vertnum = 0; vertnum < grafdat.vertnbr; vertnum ++) { /* Un-based loop on colotab */
    if ((colotab[vertnum] < 0) || (colotab[vertnum] >= colonbr)) {
      errorPrint ("main: invalid coloring (1)");
      exit       (EXIT_FAILURE);
    }
  }

  if ((chektax = memAlloc (grafdat.vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("main: out of memory (2)");
    exit       (EXIT_FAILURE);
  }
  memSet (chektax, ~0, grafdat.vertnbr * sizeof (Gnum));
  chektax -= grafdat.baseval;

  for (colonum = 0; colonum < colonbr; colonum ++) {
    Gnum                vertnum;

    for (vertnum = grafdat.baseval; vertnum < grafdat.vertnnd; vertnum ++) {
      Gnum                edgenum;

      if (colotab[vertnum - grafdat.baseval] != colonum)
        continue;

      if (chektax[vertnum] == colonum) {
        errorPrint ("main: invalid coloring (2)");
        exit       (EXIT_FAILURE);
      }
      chektax[vertnum] = colonum;

      for (edgenum = grafdat.verttax[vertnum]; edgenum < grafdat.vendtax[vertnum]; edgenum ++) {
        Gnum                vertend;

        vertend = grafdat.edgetax[edgenum];
        if (chektax[vertend] == colonum) {
          errorPrint ("main: invalid coloring (3)");
          exit       (EXIT_FAILURE);
        }
        chektax[vertend] = colonum;
      }
    }
  }

  memFree   (chektax + grafdat.baseval);
  memFree   (colotab);
  graphExit (&grafdat);

  exit (EXIT_SUCCESS);
}
