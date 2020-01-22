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
/**   NAME       : test_scotch_arch.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the SCOTCH_arch*() routines.            **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 25 jun 2014     **/
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

#define ARCHNBR                     20

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
  SCOTCH_Arch         archtab[ARCHNBR];
  SCOTCH_Num          vnumnbr = 5;
  SCOTCH_Num          vnumtab[8] = { 0, 4, 1, 5, 7, 11, 13, 15 };
  SCOTCH_Num          dimnnbr = 5;
  SCOTCH_Num          dimntab[5] = { 3, 3, 5, 4, 8 };
  SCOTCH_Num          levlnbr = 3;
  SCOTCH_Num          sizetab[3] = { 6, 3, 4 };
  SCOTCH_Num          linktab[3] = { 20, 5, 1 };
  int                 archnbr = 0;
  int                 i;

  SCOTCH_errorProg (argv[0]);

  SCOTCH_randomReset ();

  if (SCOTCH_archHcub (&archtab[archnbr ++], 4) != 0) { /* TRICK: must be in rank 0 to create sub-architecture from it */
    SCOTCH_errorPrint ("main: cannot create hcub architecture");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_archCmplt (&archtab[archnbr ++], 8) != 0) {
    SCOTCH_errorPrint ("main: cannot create cmplt architecture");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_archMesh2 (&archtab[archnbr ++], 2, 4) != 0) {
    SCOTCH_errorPrint ("main: cannot create mesh2D architecture");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_archMesh3 (&archtab[archnbr ++], 3, 4, 5) != 0) {
    SCOTCH_errorPrint ("main: cannot create mesh3D architecture");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_archMeshX (&archtab[archnbr ++], dimnnbr, dimntab) != 0) {
    SCOTCH_errorPrint ("main: cannot create meshXD architecture");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_archTleaf (&archtab[archnbr ++], levlnbr, sizetab, linktab) != 0) {
    SCOTCH_errorPrint ("main: cannot create tleaf architecture");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_archTorus2 (&archtab[archnbr ++], 2, 4) != 0) {
    SCOTCH_errorPrint ("main: cannot create torus2D architecture");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_archTorus3 (&archtab[archnbr ++], 3, 4, 5) != 0) {
    SCOTCH_errorPrint ("main: cannot create torus3D architecture");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_archTorusX (&archtab[archnbr ++], dimnnbr, dimntab) != 0) {
    SCOTCH_errorPrint ("main: cannot create torusXD architecture");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_archSub (&archtab[archnbr ++], &archtab[0], vnumnbr + 3, vnumtab) != 0) { /* TRICK: create sub-architecture of hypercube at rank 0 */
    SCOTCH_errorPrint ("main: cannot create sub-architecture (1)");
    exit (EXIT_FAILURE);
  }

  for (i = 0; i < archnbr; i ++) {
    if (SCOTCH_archSub (&archtab[i + archnbr], &archtab[i], vnumnbr, vnumtab) != 0) { /* Create sub-architectures (of sub-architecture, too) */
      SCOTCH_errorPrint ("main: cannot create sub-architecture (%d)", 2 + i);
      exit (EXIT_FAILURE);
    }
  }
  archnbr *= 2;                                   /* Number of architectures has doubled */

  if ((fileptr = fopen (argv[1], "w+")) == NULL) { /* Write all architectures to file */
    SCOTCH_errorPrint ("main: cannot open file (1)");
    exit (EXIT_FAILURE);
  }

  for (i = 0; i < archnbr; i ++) {                /* Save all architectures to same file */
    if (SCOTCH_archSave (&archtab[i], fileptr) != 0) {
      SCOTCH_errorPrint ("main: cannot save architecture (%d)", 1 + i);
      exit (EXIT_FAILURE);
    }
  }

  fclose (fileptr);

  for (i = archnbr - 1; i >= 0; i --)             /* Destroy architectures in reverse order to destroy sub-architectures first */
    SCOTCH_archExit (&archtab[i]);

  if ((fileptr = fopen (argv[1], "r")) == NULL) { /* Read all architectures from file where they were written to */
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
