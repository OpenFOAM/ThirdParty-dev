/* Copyright 2014,2015,2018 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_scotch_graph_coarsen.c             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the sequential graph coarsening         **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 15 oct 2014     **/
/**                                 to     22 may 2018     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include <stdio.h>
#if (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H))
#include <stdint.h>
#endif /* (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H)) */
#include <stdlib.h>
#include <string.h>

#include "scotch.h"

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
  SCOTCH_Num              baseval;                /* Base value                */
  SCOTCH_Graph            finegrafdat;            /* Fine graph                */
  SCOTCH_Num              finevertnbr;            /* Number of fine vertices   */
  SCOTCH_Num              finevertnnd;            /* Based end of vertex array */
  SCOTCH_Num              finevertnum;
  SCOTCH_Num *            fineverttax;
  SCOTCH_Num *            finevendtax;
  SCOTCH_Num *            fineedgetax;
  SCOTCH_Num *            finematetax;
  SCOTCH_Graph            coargrafdat;            /* Coarse graph              */
  SCOTCH_Num *            coarmulttab;            /* Multinode array           */
  SCOTCH_Num              coarvertnbr;            /* Number of coarse vertices */
  SCOTCH_Num              coarvertnum;
  SCOTCH_Num              coarverttmp;
  SCOTCH_Num              coaredgenbr;
  FILE *                  fileptr;

  SCOTCH_errorProg (argv[0]);

  if (argc != 2) {
    SCOTCH_errorPrint ("usage: %s graph_file", argv[0]);
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphInit (&finegrafdat) != 0) {     /* Initialize fine source graph */
    SCOTCH_errorPrint ("main: cannot initialize graph");
    exit (EXIT_FAILURE);
  }

  if ((fileptr = fopen (argv[1], "r")) == NULL) { /* Open fine graph file */
    SCOTCH_errorPrint ("main: cannot open file");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphLoad (&finegrafdat, fileptr, -1, 0) != 0) { /* Read fine source graph */
    SCOTCH_errorPrint ("main: cannot load graph");
    exit (EXIT_FAILURE);
  }

  fclose (fileptr);

  SCOTCH_graphData (&finegrafdat, &baseval,
                    &finevertnbr, &fineverttax, &finevendtax, NULL, NULL,
                    NULL, &fineedgetax, NULL);
  fineverttax -= baseval;
  finevendtax -= baseval;
  fineedgetax -= baseval;

  if ((finematetax = malloc (finevertnbr * sizeof (SCOTCH_Num))) == NULL) { /* Array is un-based */
    SCOTCH_errorPrint ("main: out of memory (1)");
    exit (EXIT_FAILURE);
  }

  coarvertnbr = 0;                                /* Set minimum number of vertices in graph */
  if (SCOTCH_graphCoarsenMatch (&finegrafdat, &coarvertnbr, 1.0, SCOTCH_COARSENNOMERGE, finematetax) != 0) {
    SCOTCH_errorPrint ("main: cannot compute matching");
    exit (EXIT_FAILURE);
  }

  finematetax -= baseval;                         /* Array is based from now on */
  for (finevertnum = baseval, finevertnnd = finevertnbr + baseval, coarverttmp = finevertnbr;
       finevertnum < finevertnnd; finevertnum ++) {
    SCOTCH_Num          finevertend;

    finevertend = finematetax[finevertnum];
    if ((finevertend <  baseval) ||
        (finevertend >= finevertnnd)) {
      SCOTCH_errorPrint ("main: invalid mate array (1)");
      exit (EXIT_FAILURE);
    }
    if (finevertend > finevertnum) {              /* If mates make a multinode (counted once) */
      coarverttmp --;                             /* One less coarse vertex                   */
      if (finematetax[finevertend] != finevertnum) {
        SCOTCH_errorPrint ("main: invalid mate array (2)");
        exit (EXIT_FAILURE);
      }
    }
  }
  if (coarverttmp != coarvertnbr) {
    SCOTCH_errorPrint ("main: invalid mate array (3)");
    exit (EXIT_FAILURE);
  }

  printf ("Mate array has " SCOTCH_NUMSTRING " vertices\n",
          coarvertnbr);
  printf ("Graph mated with a ratio of %lg\n", (double) coarvertnbr / (double) finevertnbr);

  if ((coarmulttab = malloc (coarvertnbr * 2 * sizeof (SCOTCH_Num))) == NULL) {
    SCOTCH_errorPrint ("main: out of memory (2)");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphCoarsenBuild (&finegrafdat, coarvertnbr, finematetax + baseval, &coargrafdat, coarmulttab) != 0) {
    SCOTCH_errorPrint ("main: cannot compute coarse graph (1)");
    exit (EXIT_FAILURE);
  }

  SCOTCH_graphSize (&coargrafdat, &coarvertnbr, &coaredgenbr);
  printf ("Coarse graph has " SCOTCH_NUMSTRING " vertices and " SCOTCH_NUMSTRING " edges\n",
          coarvertnbr,
          coaredgenbr);
  printf ("Graph coarsened with a ratio of %lg\n", (double) coarvertnbr / (double) finevertnbr);

  for (coarvertnum = 0, coarverttmp = finevertnbr, finevertnnd = finevertnbr + baseval;
       coarvertnum < coarvertnbr; coarvertnum ++) {
    SCOTCH_Num          finevertnum0;
    SCOTCH_Num          finevertnum1;

    finevertnum0 = coarmulttab[2 * coarvertnum];
    finevertnum1 = coarmulttab[2 * coarvertnum + 1];
    if ((finevertnum0 <  baseval)     ||
        (finevertnum0 >= finevertnnd) ||
        (finevertnum1 <  baseval)     ||
        (finevertnum1 >= finevertnnd)) {
      SCOTCH_errorPrint ("main: invalid multinode array (1)");
      exit (EXIT_FAILURE);
    }
    if (finevertnum0 != finevertnum1)
      coarverttmp --;
  }
  if (coarverttmp != coarvertnbr) {
    SCOTCH_errorPrint ("main: invalid multinode array (2)");
    exit (EXIT_FAILURE);
  }

  SCOTCH_graphExit (&coargrafdat);
  free             (coarmulttab);                 /* Free first multinode array */

  memset (finematetax + baseval, ~0, finevertnbr * sizeof (SCOTCH_Num)); /* Set based array to "-1"'s */

  for (finevertnum = baseval, finevertnnd = finevertnbr + baseval, coarvertnbr = 0; /* Compute simple matching */
       finevertnum < finevertnnd; finevertnum ++) {
    SCOTCH_Num          finematenum;

    if (finematetax[finevertnum] >= 0)            /* If vertex already matched */
      continue;

    coarvertnbr ++;                               /* One more coarse vertex will be created */

    finematenum = finevertnum;                    /* Assume we mate with ourselves     */
    if (fineverttax[finevertnum] == finevendtax[finevertnum]) { /* If isolated vertex  */
      while (1) {                                 /* Use first free vertex as a mate   */
        if (++ finevertnum >= finevertnnd) {      /* If no other free vertex available */
          finematetax[finematenum] = finematenum; /* Mate isolated vertex with itself  */
          goto endloop;                           /* Exit nested loops                 */
        }
        if (finematetax[finevertnum] < 0)         /* If vertex is free, keep it as mate */
          break;                                  /* Increment performed in outer loop  */
      }
    }
    else {
      SCOTCH_Num          fineedgenum;

      for (fineedgenum = fineverttax[finevertnum]; /* Find a suitable mate */
           fineedgenum < finevendtax[finevertnum]; fineedgenum ++) {
        SCOTCH_Num        finevertend;

        finevertend = fineedgetax[fineedgenum];   /* FInd end vertex     */
        if (finematetax[finevertend] < 0) {       /* If vertex is free   */
          finematenum = finevertend;              /* Keep vertex as mate */
          break;
        }
      }
    }

    finematetax[finevertnum] = finematenum;
    finematetax[finematenum] = finevertnum;
  }
endloop:
  printf ("Mate array has " SCOTCH_NUMSTRING " vertices\n",
          coarvertnbr);
  printf ("Graph mated with a ratio of %lg\n", (double) coarvertnbr / (double) finevertnbr);

  if ((coarmulttab = malloc (coarvertnbr * 2 * sizeof (SCOTCH_Num))) == NULL) {
    SCOTCH_errorPrint ("main: out of memory (3)");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphCoarsenBuild (&finegrafdat, coarvertnbr, finematetax + baseval, &coargrafdat, coarmulttab) != 0) {
    SCOTCH_errorPrint ("main: cannot compute coarse graph (2)");
    exit (EXIT_FAILURE);
  }

  SCOTCH_graphSize (&coargrafdat, &coarvertnbr, &coaredgenbr);
  printf ("Coarse graph has " SCOTCH_NUMSTRING " vertices and " SCOTCH_NUMSTRING " edges\n",
          coarvertnbr,
          coaredgenbr);
  printf ("Graph coarsened with a ratio of %lg\n", (double) coarvertnbr / (double) finevertnbr);

  SCOTCH_graphExit (&coargrafdat);
  free (coarmulttab);                             /* Free second multinode array */
  free (finematetax + baseval);

  if ((coarmulttab = malloc (finevertnbr * 2 * sizeof (SCOTCH_Num))) == NULL) { /* Number of coarse vertices and ratio not bounded */
    SCOTCH_errorPrint ("main: out of memory (3)");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphCoarsen (&finegrafdat, 1, 1.0, SCOTCH_COARSENNONE, &coargrafdat, coarmulttab) != 0) {
    SCOTCH_errorPrint ("main: cannot coarsen graph");
    exit (EXIT_FAILURE);
  }

  SCOTCH_graphSize (&coargrafdat, &coarvertnbr, NULL);
  printf ("Coarse graph has " SCOTCH_NUMSTRING " vertices and " SCOTCH_NUMSTRING " edges\n",
          coarvertnbr,
          coaredgenbr);
  printf ("Graph coarsened with a ratio of %lg\n", (double) coarvertnbr / (double) finevertnbr);

  for (coarvertnum = 0, coarverttmp = finevertnbr, finevertnnd = finevertnbr + baseval;
       coarvertnum < coarvertnbr; coarvertnum ++) {
    SCOTCH_Num          finevertnum0;
    SCOTCH_Num          finevertnum1;

    finevertnum0 = coarmulttab[2 * coarvertnum];
    finevertnum1 = coarmulttab[2 * coarvertnum + 1];
    if ((finevertnum0 <  baseval)     ||
        (finevertnum0 >= finevertnnd) ||
        (finevertnum1 <  baseval)     ||
        (finevertnum1 >= finevertnnd)) {
      SCOTCH_errorPrint ("main: invalid multinode array (3)");
      exit (EXIT_FAILURE);
    }
    if (finevertnum0 != finevertnum1)
      coarverttmp --;
  }
  if (coarverttmp != coarvertnbr) {
    SCOTCH_errorPrint ("main: invalid multinode array (4)");
    exit (EXIT_FAILURE);
  }

  SCOTCH_graphExit (&coargrafdat);
  free             (coarmulttab);                 /* Free third multinode array */
  SCOTCH_graphExit (&finegrafdat);

  exit (EXIT_SUCCESS);
}
