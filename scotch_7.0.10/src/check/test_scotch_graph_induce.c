/* Copyright 2019,2025 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_graph_induce.c                     **/
/**                                                        **/
/**   AUTHOR     : Amaury JACQUES (v6.0)                   **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operations of     **/
/**                the SCOTCH_graphInduceList() and        **/
/**                SCOTCH_graphInducePart() routines.      **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 15 apr 2019     **/
/**                                 to   : 16 apr 2019     **/
/**                # Version 7.0  : from : 13 sep 2019     **/
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
  SCOTCH_Num          baseval;
  SCOTCH_Num          orgvertnbr;
  SCOTCH_Num          orgvertnum;
  SCOTCH_GraphPart2 * orgparttab;
  SCOTCH_Graph        orggrafdat;
  SCOTCH_Graph        indgrafdat;
  SCOTCH_Num *        indlisttab;
  SCOTCH_Num          indvertnbr;
  SCOTCH_Num          indvertnum;

  if (argc != 2) {
    SCOTCH_errorPrint ("usage: %s graph_file", argv[0]);
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphInit (&orggrafdat) != 0) {
    SCOTCH_errorPrint ("main: cannot initialize graph (1)");
    exit (EXIT_FAILURE);
  }

  if ((fileptr = fopen (argv[1], "r")) == NULL) {
    SCOTCH_errorPrint ("main: cannot open file");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphLoad (&orggrafdat, fileptr, -1, 0) != 0) { /* Read source graph */
    SCOTCH_errorPrint ("main: cannot load graph");
    exit (EXIT_FAILURE);
  }

  fclose (fileptr);

  SCOTCH_graphData (&orggrafdat, &baseval, &orgvertnbr, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  if ((orgparttab = malloc (orgvertnbr * sizeof (SCOTCH_GraphPart2))) == NULL) {
    SCOTCH_errorPrint ("main: out of memory (1)");
    exit (EXIT_FAILURE);
  }
  if ((indlisttab = malloc (orgvertnbr * sizeof (SCOTCH_Num))) == NULL) { /* TRICK: size is orgvertnbr */
    SCOTCH_errorPrint ("main: out of memory (2)");
    exit (EXIT_FAILURE);
  }

  SCOTCH_randomReset ();                          /* Initialize global random generator */

  for (orgvertnum = 0; orgvertnum < orgvertnbr; orgvertnum ++) /* Random permutation of all original graph vertices */
    indlisttab[orgvertnum] = orgvertnum + baseval;
  for (orgvertnum = 0; orgvertnum < (orgvertnbr - 1); orgvertnum ++) {
    SCOTCH_Num          randval;
    SCOTCH_Num          verttmp;

    randval = SCOTCH_randomVal (orgvertnbr - orgvertnum);
    verttmp = indlisttab[orgvertnum + randval];
    indlisttab[orgvertnum + randval] = indlisttab[orgvertnum];
    indlisttab[orgvertnum] = verttmp;
  }

  indvertnbr = (orgvertnbr + 1) / 2;              /* Keep only half of the original vertices */

  memset (orgparttab, 0, orgvertnbr * sizeof (SCOTCH_GraphPart2));
  for (indvertnum = 0; indvertnum < indvertnbr; indvertnum ++) /* Flag kept vertices as belonging to part 1 */
    orgparttab[indlisttab[indvertnum] - baseval] = 1;

  if (SCOTCH_graphInit (&indgrafdat) != 0) {
    SCOTCH_errorPrint ("main: cannot initialize graph (2)");
    exit (EXIT_FAILURE);
  }
  if (SCOTCH_graphInduceList (&orggrafdat, indvertnbr, indlisttab, &indgrafdat) != 0) {
    SCOTCH_errorPrint ("main: cannot induce graph (1)");
    exit (EXIT_FAILURE);
  }
  if (SCOTCH_graphCheck (&indgrafdat) != 0) {
    SCOTCH_errorPrint ("main: invalid induced graph (1)");
    exit (EXIT_FAILURE);
  }
  SCOTCH_graphExit (&indgrafdat);

  if (SCOTCH_graphInit (&indgrafdat) != 0) {
    SCOTCH_errorPrint ("main: cannot initialize graph (3)");
    exit (EXIT_FAILURE);
  }
  if (SCOTCH_graphInducePart (&orggrafdat, indvertnbr, orgparttab, 1, &indgrafdat) != 0) {
    SCOTCH_errorPrint ("main: cannot induce graph (2)");
    exit (EXIT_FAILURE);
  }
  if (SCOTCH_graphCheck (&indgrafdat) != 0) {
    SCOTCH_errorPrint ("main: invalid induced graph (2)");
    exit (EXIT_FAILURE);
  }
  SCOTCH_graphExit (&indgrafdat);

  free (indlisttab);
  free (orgparttab);

  SCOTCH_graphExit (&orggrafdat);

  exit (EXIT_SUCCESS);
}
