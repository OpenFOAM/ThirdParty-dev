/* Copyright 2020 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_scotch_libesmumps.c                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the libesmumps routines.                **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 19 jan 2020     **/
/**                                 to   : 19 jan 2020     **/
/**                # Version 6.1  : from : 22 feb 2020     **/
/**                                 to   : 05 sep 2020     **/
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
#include "esmumps.h"

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
  SCOTCH_Num              vertnbr;
  SCOTCH_Num              vertnum;
  SCOTCH_Num *            verttab;
  SCOTCH_Num *            velotab;
  SCOTCH_Num              edgenbr;
  SCOTCH_Num *            edgetab;
  SCOTCH_Num *            elentab;
  SCOTCH_Num *            lasttab;
  SCOTCH_Num *            lentab;
  SCOTCH_Num *            nvtab;
  SCOTCH_Num *            petab;

  SCOTCH_errorProg (argv[0]);

  if (argc != 2) {
    SCOTCH_errorPrint ("usage: %s graph_file", argv[0]);
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphInit (&grafdat) != 0) {         /* Initialize source graph */
    SCOTCH_errorPrint ("main: cannot initialize graph");
    exit (EXIT_FAILURE);
  }

  if ((fileptr = fopen (argv[1], "r")) == NULL) {
    SCOTCH_errorPrint ("main: cannot open file (1)");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphLoad (&grafdat, fileptr, 1, 0) != 0) { /* Read source graph with base 1 for esmumps */
    SCOTCH_errorPrint ("main: cannot load graph");
    exit (EXIT_FAILURE);
  }

  fclose (fileptr);

  SCOTCH_graphData (&grafdat, NULL, &vertnbr, &verttab, NULL, &velotab, NULL, &edgenbr, &edgetab, NULL);

  if (((elentab = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL) ||
      ((lasttab = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL) ||
      ((lentab  = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL) ||
      ((nvtab   = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL) ||
      ((petab   = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL)) {
    SCOTCH_errorPrint ("main: out of memory");
    exit (EXIT_FAILURE);
  }

  memcpy (petab, verttab, vertnbr * sizeof (SCOTCH_Num)); /* Prepare graph topology arrays */
  for (vertnum = 0; vertnum < vertnbr; vertnum ++)
    lentab[vertnum] = verttab[vertnum + 1] - verttab[vertnum];

  if (esmumps (vertnbr, 0, petab, edgenbr + 1, lentab, edgetab, nvtab, elentab, lasttab) != 0) {
    SCOTCH_errorPrint ("main: cannot run esmumps");
    exit (EXIT_FAILURE);
  }

#ifdef ESMUMPS_HAS_ESMUMPSV
  if (velotab != NULL) {
    memcpy (petab, verttab, vertnbr * sizeof (SCOTCH_Num)); /* Prepare graph topology arrays */
    for (vertnum = 0; vertnum < vertnbr; vertnum ++)
      lentab[vertnum] = verttab[vertnum + 1] - verttab[vertnum];

    memcpy (nvtab, velotab, vertnbr * sizeof (SCOTCH_Num)); /* Prepare node multiplicity array */

    if (esmumpsv (vertnbr, 0, petab, edgenbr + 1, lentab, edgetab, nvtab, elentab, lasttab) != 0) {
      SCOTCH_errorPrint ("main: cannot run esmumpsv");
      exit (EXIT_FAILURE);
    }
  }
#endif /* ESMUMPS_HAS_ESMUMPSV */

  free (petab);
  free (nvtab);
  free (lentab);
  free (lasttab);
  free (elentab);
  SCOTCH_graphExit (&grafdat);

  exit (EXIT_SUCCESS);
}
