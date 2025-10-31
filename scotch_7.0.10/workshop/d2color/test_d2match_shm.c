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
/**   NAME       : test_d2match_shm.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests shared-memory         **/
/**                distance-2 graph matching algorithms.   **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 05 jun 2021     **/
/**                                 to   : 05 jun 2021     **/
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

typedef struct MateData_ {
  int                       seedval;              /*+ Random seed for all threads +*/
  Gnum                      vemanbr;              /*+ Number of matched vertices  +*/
} MateData;

/*+ The block data structure +*/

typedef struct MateGroup_ {
  MateData *                datatab;              /*+ Thread test data               +*/
  Graph *                   grafptr;              /*+ Pointer to graph               +*/
  Gnum *                    randtax;              /*+ Array of vertex random values  +*/
  Gnum *                    vmaxtax;              /*+ Array of maximum vertex values +*/
  Gnum *                    matetax;              /*+ Vertex matching array          +*/
  Gnum                      stepnbr;              /*+ Number of steps                +*/
  int                       revaval;              /*+ Return value for threads       +*/
} MateGroup;

/*************************/
/*                       */
/* The threaded routine. */
/*                       */
/*************************/

static
void
graphD2matchReduce (
MateData * restrict const   vlocptr,              /* Pointer to local value  */
MateData * restrict const   vremptr)              /* Pointer to remote value */
{
  vlocptr->vemanbr += vremptr->vemanbr;
}

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

void
graphD2matchShm (
ThreadDescriptor * restrict const descptr,
volatile MateGroup * restrict     grouptr)
{
  IntRandContext      randdat;
  Gnum * restrict     vnumtab;
  Gnum                vnumnbr;
  Gnum * restrict     queutab;
  Gnum                queunbr;
  Gnum                stepnbr;
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
  MateData * const              dataptr = &grouptr->datatab[thrdnum];
  Gnum *                        randtax = grouptr->randtax;
  Gnum *                        vmaxtax = grouptr->vmaxtax;
  Gnum *                        matetax = grouptr->matetax;

  if (grafptr->vertnbr <= 0)                      /* Do not consider empty graphs */
    return;

  vertbas = grafptr->baseval + DATASCAN (grafptr->vertnbr, thrdnbr, thrdnum);
  vertnnd = grafptr->baseval + DATASCAN (grafptr->vertnbr, thrdnbr, thrdnum + 1);
  vnumnbr = vertnnd - vertbas;

  if (memAllocGroup ((void **) (void *)
                     &vnumtab, (size_t) (vnumnbr * sizeof (Gnum)),
                     &queutab, (size_t) (vnumnbr * sizeof (Gnum) * 2), NULL) == NULL) { /* A local vertex and its non-local mate can be enqueued */
    errorPrint ("graphD2matchShm: out of memory");
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
    vmaxtax[vertnum] = graphD2match2 (grafptr, randtax, vertnum);
    vnumtab[vertnum - vertbas] = vertnum;         /* Fill list of remaining vertices */
  }

  threadBarrier (descptr);                        /* Wait until all maximum vertices computed */

  for (stepnbr = 0, vrmnnbr = grafptr->vertnbr; ; ) { /* Outer loop on steps */
    Gnum                vnumidx;
    Gnum                vnumtmp;
    Gnum                vemanbr;

    vemanbr = 0;
    queunbr = 0;
    for (vnumidx = vnumtmp = 0; vnumidx < vnumnbr; vnumidx ++) {
      Gnum                vertnum;
      Gnum                edgenum;
      Gnum                vertend;

      vertnum = vnumtab[vnumidx];                 /* Get vertex number from index array */

      if (randtax[vertnum] < 0)                   /* If local vertex has been mated by a distant neighbor, skip it */
	continue;

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

      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) { /* Select first unmated neighbor */
        vertend = edgetax[edgenum];
        if (randtax[vertend] >= 0)
          break;
      }

      if (edgenum < vendtax[vertnum]) {           /* If suitable mate found                       */
        matetax[vertend] = vertnum;               /* Mate neighbor mate with vertex               */
        randtax[vertend] = -1;                    /* Do no longer consider neighbor for dominance */
        vemanbr ++;
        queutab[queunbr ++] = vertend;            /* Enqueue mate for further processing */
      }
      else
        vertend = vertnum;
      matetax[vertnum] = vertend;                 /* Mate vertex with its neighbor              */
      randtax[vertnum] = -1;                      /* Do no longer consider vertex for dominance */
      vemanbr ++;
      queutab[queunbr ++] = vertnum;              /* Enqueue vertex for further processing */
    }
    vnumnbr = vnumtmp;                            /* Set new size of vertex number array */

    stepnbr ++;                                   /* One more step done    */
    dataptr->vemanbr = vemanbr;                   /* Prepare for reduction */
    threadReduce (descptr, (void *) dataptr, sizeof (MateData), (ThreadReduceFunc) graphD2matchReduce, 0);
    if (thrdnum == 0)
      printf ("Mater " GNUMSTRING " : " GNUMSTRING "\n", stepnbr, dataptr->vemanbr);
    vrmnnbr -= grouptr->datatab[0].vemanbr;       /* Update global number of remaining vertices */
    if (vrmnnbr <= 0)                             /* Exit outer loop if nothing more to do      */
      break;

    while (queunbr > 0) {
      Gnum                vertnum;
      Gnum                edgenum;

      vertnum = queutab[-- queunbr];

      vmaxtax[vertnum] = graphD2match2 (grafptr, randtax, vertnum); /* Recompute maximum vertex of vertex itself */

      for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++) { /* For all distance-1 neighbors of colored vertex */
        Gnum                vertend;

        vertend = edgetax[edgenum];
        vmaxtax[vertend] = graphD2match2 (grafptr, randtax, vertend); /* Recompute maximum vertex of distance-1 neighbor */
      }
    }

    threadBarrier (descptr);                      /* Wait until all maximum vertex updates computed */
  }

  memFree (vnumtab);                              /* Free group leader */

  if (thrdnum == 0)                               /* Set number of colors */
    grouptr->stepnbr = stepnbr;
}

int
graphD2match (
Graph * restrict const      grafptr,
Gnum * restrict const       matetab,
int                         thrdnbr)              /*+ Desired nunber of threads, or -1 +*/
{
  ThreadContext       contdat;
  IntRandContext      randdat;
  MateGroup           groudat;
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
    SCOTCH_errorPrint ("graphD2match: cannot initialize thread context");
    return (1);
  }

  thrdnbr = threadContextNbr (&contdat);
  printf ("%d threads in context\n", thrdnbr);

  if (memAllocGroup ((void **) (void *)
                     &groudat.datatab, (size_t) (thrdnbr * sizeof (MateData)),
                     &groudat.randtax, (size_t) (vertnbr * sizeof (Gnum)),
                     &groudat.vmaxtax, (size_t) (vertnbr * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("graphD2match: out of memory");
    return (1);
  }
  groudat.grafptr  = grafptr;
  groudat.matetax  = matetab - baseval;
  groudat.randtax -= baseval;
  groudat.vmaxtax -= baseval;
  groudat.stepnbr  = 0;                           /* In case of empty graph         */
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

  threadLaunch (&contdat, (ThreadFunc) graphD2matchShm, (void *) &groudat);

  memFree (groudat.datatab);

  clockStop (&timedat);
  printf    ("time: %g\n", (double) clockVal (&timedat));

  threadContextExit (&contdat);

  return ((groudat.revaval != 0) ? -1 : groudat.stepnbr);
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
  int                 thrdnbr;

  if ((argc < 3) || (argc > 4)) {
    errorPrint ("usage: %s graph matching thrdnbr", argv[0]);
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

  if ((matetab = memAlloc (grafdat.vertnbr * sizeof (Gnum))) == NULL) {
    errorPrint ("main: out of memory (1)");
    exit       (EXIT_FAILURE);
  }

  memSet (matetab, ~0, grafdat.vertnbr * sizeof (Gnum)); /* Test for internal errors */

  graphD2match (&grafdat, matetab, thrdnbr);

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
