/* Copyright 2015,2016,2018,2019 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_scotch_arch_deco.c                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the SCOTCH_arch*() routines.            **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 28 dec 2015     **/
/**                                 to     25 sep 2019     **/
/**                                                        **/
/**   DATES      : # This code derives from that of        **/
/**                  test_scotch_arch.c                    **/
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

#define ARCHNBR                     6

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
  SCOTCH_Num          vertnbr;
  SCOTCH_Graph        grafdat;
  SCOTCH_Strat        stradat;
  SCOTCH_Arch         archtab[ARCHNBR];
  int                 archnbr = 0;
  int                 i;
  SCOTCH_Num          listnbr = 5;
  SCOTCH_Num          listtab[] = { 0, 4, 1, 5, 7 };

  SCOTCH_errorProg (argv[0]);

  if (argc != 3) {
    SCOTCH_errorPrint ("usage: %s graph_file architecture_file", argv[0]);
    exit (EXIT_FAILURE);
  }

  SCOTCH_randomReset ();

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
  if (vertnbr < 8) {
    SCOTCH_errorPrint ("main: graph is too small");
    exit (EXIT_FAILURE);
  }

  if ((fileptr = fopen (argv[2], "w")) == NULL) {
    SCOTCH_errorPrint ("main: cannot open file");
    exit (EXIT_FAILURE);
  }

  SCOTCH_stratInit (&stradat);

  for (i = 1, archnbr = 0; i <= 2; i ++) {
    int                 j;
    int                 o;

    if (i == 1)
      o = SCOTCH_archBuild0 (&archtab[archnbr], &grafdat, vertnbr, NULL, &stradat);
    else
      o = SCOTCH_archBuild2 (&archtab[archnbr], &grafdat, vertnbr, NULL);
    if (o != 0) {
      SCOTCH_errorPrint ("main: cannot create decomposition-described architecture (%d)", 2 * i - 1);
      exit (EXIT_FAILURE);
    }

    if (i == 1)
      o = SCOTCH_archBuild0 (&archtab[archnbr + 1], &grafdat, listnbr, listtab, &stradat);
    else
      o = SCOTCH_archBuild2 (&archtab[archnbr + 1], &grafdat, listnbr, listtab);
    if (o != 0) {
      SCOTCH_errorPrint ("main: cannot create decomposition-described architecture (%d)", 2 * i);
      exit (EXIT_FAILURE);
    }

    if (i == 2) {
      if (SCOTCH_archSub (&archtab[archnbr + 2], &archtab[archnbr], listnbr, listtab) != 0) {
        SCOTCH_errorPrint ("main: cannot create sub-architecture (%d)", i);
        exit (EXIT_FAILURE);
      }
    }

    for (j = 0; j < (1 + i); j ++) {
      if (SCOTCH_archSave (&archtab[archnbr + j], fileptr) != 0) {
        SCOTCH_errorPrint ("main: cannot save architecture (%d)", archnbr + 1 + j);
        exit (EXIT_FAILURE);
      }
    }

    archnbr += (1 + i);
  }

  fclose (fileptr);

  SCOTCH_stratExit (&stradat);

  for (i = archnbr - 1; i >= 0; i --)             /* Destroy architectures in reverse order to destroy sub-architectures first */
    SCOTCH_archExit (&archtab[i]);

  SCOTCH_graphExit (&grafdat);

  if ((fileptr = fopen (argv[2], "r")) == NULL) { /* Read all architectures from file where they were written to */
    SCOTCH_errorPrint ("main: cannot open file (2)");
    exit (EXIT_FAILURE);
  }

  for (i = 0; i < archnbr; i ++) {
    if ((SCOTCH_archInit (&archtab[i])          != 0) ||
        (SCOTCH_archLoad (&archtab[i], fileptr) != 0)) {
      SCOTCH_errorPrint ("main: cannot load architecture (%d)", 1 + i);
      exit (EXIT_FAILURE);
    }
  }

  fclose (fileptr);

  for (i = 0; i < archnbr; i ++)                  /* Destroy architectures in any order, as they are now all autonomous from each other */
    SCOTCH_archExit (&archtab[i]);

  exit (EXIT_SUCCESS);
}
