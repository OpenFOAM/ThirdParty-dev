/* Copyright 2019 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_scotch_libmetis.c                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the libscotchmetis routines.            **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 15 may 2019     **/
/**                                 to     19 may 2019     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include <math.h>
#include <stdio.h>
#if (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H))
#include <stdint.h>
#endif /* (((defined __STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || (defined HAVE_STDINT_H)) */
#include <stdlib.h>
#include <string.h>

#include "scotch.h"
#include "metis.h"                                /* Our "metis.h" file */

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
  FILE *                  fileptr;
  SCOTCH_Graph            grafdat;
  SCOTCH_Num              baseval;
  SCOTCH_Num              vertnbr;
  SCOTCH_Num *            verttab;
  SCOTCH_Num *            velotab;
  SCOTCH_Num *            edgetab;
  SCOTCH_Num *            edlotab;
  SCOTCH_Num              edgecut;
  SCOTCH_Num *            parttab;
  SCOTCH_Num *            peritab;
#if (SCOTCH_METIS_VERSION == 3)
  SCOTCH_Num              fwgtval;
#endif /* (SCOTCH_METIS_VERSION == 3) */

  const SCOTCH_Num          foptval = 0;
  const SCOTCH_Num          partnbr = 9;
#if (SCOTCH_METIS_VERSION == 5)
  const SCOTCH_Num          awgttab[9] = { 2, 2, 1, 2, 2, 3, 1, 1, 1 };
  const SCOTCH_Num          nconval = 1;
  const double              kbaltab[1] = { 0.05 };
#endif /* (SCOTCH_METIS_VERSION == 5) */

  SCOTCH_errorProg (argv[0]);

  if (argc != 2) {
    SCOTCH_errorPrint ("usage: %s graph_file", argv[0]);
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphInit (&grafdat) != 0) {         /* Initialize source graph */
    SCOTCH_errorPrint ("main: cannot initialize graph");
    exit (EXIT_FAILURE);
  }

  if ((fileptr = fopen (argv[1], "r")) == NULL) { /* Read a square 2D grid graph */
    SCOTCH_errorPrint ("main: cannot open file (1)");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphLoad (&grafdat, fileptr, -1, 0) != 0) { /* Read source graph */
    SCOTCH_errorPrint ("main: cannot load graph");
    exit (EXIT_FAILURE);
  }

  fclose (fileptr);

  SCOTCH_graphData (&grafdat, &baseval, &vertnbr, &verttab, NULL, &velotab, NULL, NULL, &edgetab, &edlotab);

  if (((parttab = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL) ||
      ((peritab = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL)) {
    SCOTCH_errorPrint ("main: out of memory");
    exit (EXIT_FAILURE);
  }

#if (SCOTCH_METIS_VERSION == 3)
  fwgtval = ((velotab != NULL) ? 2 : 0) +
            ((edlotab != NULL) ? 1 : 0);

  if (METIS_PartGraphKway (&vertnbr, verttab, edgetab, velotab, edlotab,
                           &fwgtval, &baseval, &partnbr, &foptval, &edgecut, parttab) != METIS_OK) {
    SCOTCH_errorPrint ("main: error in METIS_V3_PartGraphKway");
    exit (EXIT_FAILURE);
  }

  if (METIS_PartGraphRecursive (&vertnbr, verttab, edgetab, velotab, edlotab,
                                &fwgtval, &baseval, &partnbr, &foptval, &edgecut, parttab) != METIS_OK) {
    SCOTCH_errorPrint ("main: error in METIS_V3_PartGraphRecursive");
    exit (EXIT_FAILURE);
  }

  fwgtval &= ~2;                                  /* Take vertex load array as communication volume array */
  if (METIS_PartGraphVKway (&vertnbr, verttab, edgetab, NULL, velotab,
                            &fwgtval, &baseval, &partnbr, &foptval, &edgecut, parttab) != METIS_OK) {
    SCOTCH_errorPrint ("main: error in METIS_V3_PartGraphKway");
    exit (EXIT_FAILURE);
  }

  if (METIS_EdgeND (&vertnbr, verttab, edgetab, &baseval, &foptval, peritab, parttab) != METIS_OK) {
    SCOTCH_errorPrint ("main: error in METIS_V3_EdgeND");
    exit (EXIT_FAILURE);
  }

  if (METIS_NodeND (&vertnbr, verttab, edgetab, &baseval, &foptval, peritab, parttab) != METIS_OK) {
    SCOTCH_errorPrint ("main: error in METIS_V3_NodeND");
    exit (EXIT_FAILURE);
  }

  if (METIS_NodeWND (&vertnbr, verttab, edgetab, velotab, &baseval, &foptval, peritab, parttab) != METIS_OK) {
    SCOTCH_errorPrint ("main: error in METIS_V3_NodeWND");
    exit (EXIT_FAILURE);
  }
#endif /* (SCOTCH_METIS_VERSION == 3) */

#if (SCOTCH_METIS_VERSION == 5)
  SCOTCH_graphBase (&grafdat, 0);                 /* MeTiS v5 no longer handles graph base */

  if (METIS_PartGraphKway (&vertnbr, &nconval, verttab, edgetab, velotab, NULL, edlotab,
                           &partnbr, awgttab, kbaltab, &foptval, &edgecut, parttab) != METIS_OK) {
    SCOTCH_errorPrint ("main: error in METIS_V5_PartGraphKway");
    exit (EXIT_FAILURE);
  }

  if (METIS_PartGraphRecursive (&vertnbr, &nconval, verttab, edgetab, velotab, NULL, edlotab,
                                &partnbr, awgttab, kbaltab, &foptval, &edgecut, parttab) != METIS_OK) {
    SCOTCH_errorPrint ("main: error in METIS_V5_PartGraphRecursive");
    exit (EXIT_FAILURE);
  }

  if (METIS_NodeND (&vertnbr, verttab, edgetab, velotab, &foptval, peritab, parttab) != METIS_OK) {
    SCOTCH_errorPrint ("main: error in METIS_V5_NodeND");
    exit (EXIT_FAILURE);
  }
#endif /* (SCOTCH_METIS_VERSION == 5) */

  free (peritab);
  free (parttab);
  SCOTCH_graphExit (&grafdat);

  exit (EXIT_SUCCESS);
}
