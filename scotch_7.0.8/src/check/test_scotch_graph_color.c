/* Copyright 2012,2014,2018,2025 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_scotch_graph_color.c               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the SCOTCH_graphColor() routine.        **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 06 jan 2012     **/
/**                                 to   : 22 may 2018     **/
/**                # Version 7.0  : from : 04 jul 2025     **/
/**                                 to   : 04 jul 2025     **/
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

#include "../libscotch/module.h"
#include "../libscotch/common.h"

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
  FILE *              fileptr;
  SCOTCH_Graph        grafdat;
  SCOTCH_Num          vertnbr;
  SCOTCH_Num          vertnum;
  SCOTCH_Num          colonbr;
  SCOTCH_Num          colonum;
  SCOTCH_Num *        colotab;
  SCOTCH_Num *        cnbrtab;

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
    SCOTCH_errorPrint ("main: cannot open file");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphLoad (&grafdat, fileptr, -1, 0) != 0) { /* Read source graph */
    SCOTCH_errorPrint ("main: cannot load graph");
    exit (EXIT_FAILURE);
  }

  fclose (fileptr);

  SCOTCH_graphSize (&grafdat, &vertnbr, NULL);

  if ((colotab = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL) {
    SCOTCH_errorPrint ("main: out of memory (1)");
    exit (EXIT_FAILURE);
  }

  if ((cnbrtab = malloc (vertnbr * sizeof (SCOTCH_Num))) == NULL) {
    SCOTCH_errorPrint ("main: out of memory (1)");
    exit (EXIT_FAILURE);
  }
  memset (cnbrtab, 0, vertnbr * sizeof (SCOTCH_Num));

  if (SCOTCH_graphColor (&grafdat, colotab, &colonbr, 0) != 0) {
    SCOTCH_errorPrint ("main: cannot color graph");
    exit (EXIT_FAILURE);
  }

  printf ("Number of colors: %ld\n", (long) colonbr);

  for (vertnum = 0; vertnum < vertnbr; vertnum ++) /* Sum-up color histogram */
    cnbrtab[colotab[vertnum]] ++;

  for (colonum = 0; colonum < colonbr; colonum ++)
    printf ("Color %5ld: %ld\n",
            (long) colonum,
            (long) cnbrtab[colonum]);

  free (cnbrtab);
  free (colotab);
  SCOTCH_graphExit (&grafdat);

  exit (EXIT_SUCCESS);
}
