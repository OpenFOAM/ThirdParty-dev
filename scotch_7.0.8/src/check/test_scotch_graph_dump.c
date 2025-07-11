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
/**   NAME       : test_scotch_graph_dump.c                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the SCOTCH_graphDump() routine.         **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 08 jan 2019     **/
/**                                 to   : 08 jan 2019     **/
/**                # Version 6.1  : from : 24 jun 2021     **/
/**                                 to   : 24 jun 2021     **/
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

int                         testGraphBuild (SCOTCH_Graph *);

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

  SCOTCH_errorProg (argv[0]);

  if (argc != 2) {
    SCOTCH_errorPrint ("usage: %s graph_file", argv[0]);
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphInit (&grafdat) != 0) {         /* Initialize source graph */
    SCOTCH_errorPrint ("main: cannot initialize graph");
    exit (EXIT_FAILURE);
  }

  if (testGraphBuild (&grafdat) != 0) {           /* Build source graph */
    SCOTCH_errorPrint ("main: cannot build graph");
    exit (EXIT_FAILURE);
  }

  if ((fileptr = fopen (argv[1], "w")) == NULL) {
    SCOTCH_errorPrint ("main: cannot open file");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphSave (&grafdat, fileptr) != 0) { /* Save source graph */
    SCOTCH_errorPrint ("main: cannot save graph");
    exit (EXIT_FAILURE);
  }

  fclose (fileptr);

  SCOTCH_graphExit (&grafdat);

  exit (EXIT_SUCCESS);
}

/* The testGraphBuild() code will be pasted after this line */

