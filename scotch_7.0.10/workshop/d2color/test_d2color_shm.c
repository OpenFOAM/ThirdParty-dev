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
/**   NAME       : test_d2color_shm.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests shared-memory         **/
/**                distance-2 graph coloring algorithms.   **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 02 jun 2021     **/
/**                                 to   : 02 jun 2021     **/
/**                                                        **/
/************************************************************/

#define CENTRALIZED_RANDOM

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "common_thread.h"
#include "common_thread_system.h"
#include "graph.h"

/*
**  The type and structure definitions.
*/

/*+ The thread-specific data block. +*/

typedef struct ColoData_ {
  int                       seedval;              /*+ Random seed for all threads +*/
  Gnum                      veconbr;              /*+ Number of colored vertices  +*/
} ColoData;

/*+ The block data structure +*/

typedef struct ColoGroup_ {
  ColoData *                datatab;              /*+ Thread test data              +*/
  Graph *                   grafptr;              /*+ Pointer to graph              +*/
  Gnum *                    randtax;              /*+ Array of vertex random values +*/
  Gnum *                    vmaxtax;              /*+ Array of maximum vertices     +*/
  Gnum *                    colotax;              /*+ Vertex coloring array         +*/
  Gnum                      colonbr;              /*+ Number of colors              +*/
  int                       revaval;              /*+ Return value for threads      +*/
} ColoGroup;

/*************************/
/*                       */
/* The threaded routine. */
/*                       */
/*************************/

static
void
graphD2colorReduce (
ColoData * restrict const   vlocptr,              /* Pointer to local value  */
ColoData * restrict const   vremptr)              /* Pointer to remote value */
{
  vlocptr->veconbr += vremptr->veconbr;
}

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

void
graphD2colorShm (
ThreadDescriptor * restrict const descptr,
volatile ColoGroup * restrict     grouptr)
{
  IntRandContext      randdat;
  Gnum * restrict     vnumtab;
  Gnum                vnumnbr;
  Gnum * restrict     queutab;
  Gnum                queunbr;
  Gnum                colonbr;
  Gnum                vrmnnbr;                    /*+ Global number of vertices to process +*/
  Gnum                vertbas;
  Gnum                vertnnd;
  Gnum                vertnum;

  const Graph * restrict const  grafptr = grouptr->grafptr;
  const Gnum * restrict const   verttax = grafptr->verttax;
  const Gnum * restrict const   vendtax = grafptr->vendtax;
  const Gnum * restrict const   edgetax = grafptr->edgetax;
  const int                     thrdnbr = threadNbr (descptr);
  const int                     thrdnum = threadNum (descptr);
  ColoData * const              dataptr = &grouptr->datatab[thrdnum];
  Gnum *                        randtax = grouptr->randtax;
  Gnum *                        vmaxtax = grouptr->vmaxtax;
  Gnum *                        colotax = grouptr->colotax;

  if (grafptr->vertnbr <= 0)                      /* Do not consider empty graphs */
    return;

  vertbas = grafptr->baseval + DATASCAN (grafptr->vertnbr, thrdnbr, thrdnum);
  vertnnd = grafptr->baseval + DATASCAN (grafptr->vertnbr, thrdnbr, thrdnum + 1);
  vnumnbr = vertnnd - vertbas;

  if (memAllocGroup ((void **) (void *)
                     &vnumtab, (size_t) (vnumnbr * sizeof (Gnum)),
                     &queutab, (size_t) (vnumnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("graphD2colorShm: out of memory");
    grouptr->revaval = 1;                         /* No need of lock since always same value written */
  }

  intRandProc (&randdat, 1);                      /* Use a fixed random seed for the sake of reproducibility */
  intRandSeed (&randdat, dataptr->seedval);       /* Each thread has its own seed                            */

#ifndef CENTRALIZED_RANDOM
  for (vertnum = vertbas; vertnum < vertnnd; vertnum ++)
    randtax[vertnum] = intRandVal2 (&randdat) & (((Gunum) -1) >> 1); /* Always get positive numbers */
#endif /* CENTRALIZED_RANDOM */

  threadBarrier (descptr);                        /* Wait until all random numbers assigned */
  if (grouptr->revaval != 0)                      /* In case of error, all threads abort    */
    return;

  for (vertnum = vertbas; vertnum < vertnnd; vertnum ++) { /* Compute maximum vertex of local vertices */
    vmaxtax[vertnum] = graphD2color2 (grafptr, randtax, vertnum);
    vnumtab[vertnum - vertbas] = vertnum;         /* Fill list of remaining vertices */
  }

  threadBarrier (descptr);                        /* Wait until all maximum vertices computed */

  for (colonbr = 0, vrmnnbr = grafptr->vertnbr; ; ) { /* Outer loop on colors */
    Gnum                vnumidx;
    Gnum                vnumtmp;

    vnumtmp = 0;
    queunbr = 0;
    for (vnumidx = 0; vnumidx < vnumnbr; vnumidx ++) {
      Gnum                vertnum;
      Gnum                edgenum;

      vertnum = vnumtab[vnumidx];                 /* Get vertex number from index array */
      if (vmaxtax[vertnum] != vertnum) {          /* If vertex is not locally dominant  */
        vnumtab[vnumtmp ++] = vertnum;            /* Re-enqueue it                      */
        continue;
      }

      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) {
        Gnum                vertend;

        vertend = edgetax[edgenum];
	if (vmaxtax[vertend] != vertnum)
          break;
      }
      if (edgenum < vendtax[vertnum]) {           /* If vertex is not distance-2 dominant */
        vnumtab[vnumtmp ++] = vertnum;            /* Re-enqueue it                        */
        continue;
      }

      colotax[vertnum] = colonbr;
      randtax[vertnum] = -1;
      queutab[queunbr ++] = vertnum;              /* Enqueue it for further processing */
    }
    dataptr->veconbr = vnumnbr - vnumtmp;         /* Number of locally colored vertices is difference in index array size */
    vnumnbr = vnumtmp;                            /* Set new size of vertex number array                                  */

    colonbr ++;                                   /* One more color used */
    threadReduce (descptr, (void *) dataptr, sizeof (ColoData), (ThreadReduceFunc) graphD2colorReduce, 0);
    if (thrdnum == 0)
      printf ("Color " GNUMSTRING " : " GNUMSTRING "\n", colonbr, dataptr->veconbr);
    vrmnnbr -= grouptr->datatab[0].veconbr;       /* Update global number of remaining vertices */
    if (vrmnnbr <= 0)                             /* Exit outer loop if nothing more to do      */
      break;

    while (queunbr > 0) {
      Gnum                vertnum;
      Gnum                edgenum;

      vertnum = queutab[-- queunbr];

      vmaxtax[vertnum] = graphD2color2 (grafptr, randtax, vertnum); /* Recompute maximum vertex of vertex itself */

      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) { /* For all distance-1 neighbors of colored vertex */
        Gnum                vertend;

        vertend = edgetax[edgenum];
        vmaxtax[vertend] = graphD2color2 (grafptr, randtax, vertend); /* Recompute maximum vertex of distance-1 neighbor */
      }
    }

    threadBarrier (descptr);                      /* Wait until all maximum vertex updates computed */
  }

  memFree (vnumtab);                              /* Free group leader */

  if (thrdnum == 0)                               /* Set number of colors */
    grouptr->colonbr = colonbr;
}

int
graphD2color (
Graph * restrict const      grafptr,
Gnum * restrict const       colotab,
int                         thrdnbr)              /*+ Desired nunber of threads, or -1 +*/
{
  ThreadContext       contdat;
  IntRandContext      randdat;
  ColoGroup           groudat;
  Clock               timedat;                    /* Timing variable */
  Gnum                thrdnum;
#ifdef CENTRALIZED_RANDOM
  Gnum                vertnum;
#endif /* CENTRALIZED_RANDOM */

  const Gnum                baseval = grafptr->baseval;
  const Gnum                vertnbr = grafptr->vertnbr;

  if (thrdnbr < 1)
    thrdnbr = threadSystemCoreNbr ();

  if (threadContextInit (&contdat, thrdnbr, NULL) != 0) {
    SCOTCH_errorPrint ("graphD2color: cannot initialize thread context");
    return (1);
  }

  thrdnbr = threadContextNbr (&contdat);
  printf ("%d threads in context\n", thrdnbr);

  if (memAllocGroup ((void **) (void *)
                     &groudat.datatab, (size_t) (thrdnbr * sizeof (ColoData)),
                     &groudat.randtax, (size_t) (vertnbr * sizeof (Gnum)),
                     &groudat.vmaxtax, (size_t) (vertnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("graphD2color: out of memory");
    return (1);
  }
  groudat.grafptr  = grafptr;
  groudat.colotax  = colotab - baseval;
  groudat.colonbr  = 0;                           /* In case of empty graph */
  groudat.randtax -= baseval;
  groudat.vmaxtax -= baseval;
  groudat.revaval  = 0;                           /* Assume everything will go well */

  intRandProc (&randdat, 1);                      /* Use a fixed random seed for the sake of reproducibility */
  intRandSeed (&randdat, 1);

#ifdef CENTRALIZED_RANDOM
  for (vertnum = grafptr->baseval; vertnum < grafptr->vertnnd; vertnum ++)
    groudat.randtax[vertnum] = intRandVal2 (&randdat) & (((Gunum) -1) >> 1); /* Always get positive numbers */
#endif /* CENTRALIZED_RANDOM */

  clockInit  (&timedat);
  clockStart (&timedat);

  for (thrdnum = 0; thrdnum < thrdnbr; thrdnum ++)
    groudat.datatab[thrdnum].seedval = intRandVal2 (&randdat) & (((Gunum) -1) >> 1); /* Always get positive numbers */

  threadLaunch (&contdat, (ThreadFunc) graphD2colorShm, (void *) &groudat);

  memFree (groudat.datatab);

  clockStop (&timedat);
  printf    ("time: %g\n", (double) clockVal (&timedat));

  threadContextExit (&contdat);

  return ((groudat.revaval != 0) ? -1 : groudat.colonbr);
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
  int                 thrdnbr;

  if ((argc < 3) || (argc > 4)) {
    errorPrint ("usage: %s graph mapping thrdnbr", argv[0]);
    exit       (EXIT_FAILURE);
  }

  thrdnbr = -1;
  if ((argc == 4) && ((thrdnbr = atoi (argv[3])) <= 0)) {
    errorPrint ("main: invalid number of threads");
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

  colonbr = graphD2color (&grafdat, colotab, thrdnbr);

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
