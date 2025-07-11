/* Copyright 2014,2015,2018,2021,2023-2025 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_scotch_dgraph_check.c              **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Clement BARTHELEMY                      **/
/**                                                        **/
/**   FUNCTION   : This module tests the operation of      **/
/**                the SCOTCH_dgraphCheck() routine.       **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 28 sep 2014     **/
/**                                 to   : 22 may 2018     **/
/**                # Version 6.1  : from : 28 dec 2021     **/
/**                                 to   : 28 dec 2021     **/
/**                # Version 7.0  : from : 03 jul 2023     **/
/**                                 to   : 04 jul 2025     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include <mpi.h>

#include "../libscotch/module.h"
#include "../libscotch/common.h"

#include "ptscotch.h"

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
  MPI_Comm            proccomm;
  int                 procglbnbr;                 /* Number of processes sharing graph data */
  int                 proclocnum;                 /* Number of this process                 */
  SCOTCH_Dgraph       grafdat;
  FILE *              file;
#ifdef SCOTCH_PTHREAD
  int                 thrdreqlvl;
  int                 thrdprolvl;
#endif /* SCOTCH_PTHREAD */

  SCOTCH_errorProg (argv[0]);

#ifdef SCOTCH_PTHREAD
  thrdreqlvl = MPI_THREAD_MULTIPLE;
  if (MPI_Init_thread (&argc, &argv, thrdreqlvl, &thrdprolvl) != MPI_SUCCESS)
    SCOTCH_errorPrint ("main: Cannot initialize (1)");
#else /* SCOTCH_PTHREAD */
  if (MPI_Init (&argc, &argv) != MPI_SUCCESS)
    SCOTCH_errorPrint ("main: Cannot initialize (2)");
#endif /* SCOTCH_PTHREAD */

  if (argc != 2) {
    SCOTCH_errorPrint ("usage: %s graph_file", argv[0]);
    exit (EXIT_FAILURE);
  }

  proccomm = MPI_COMM_WORLD;
  MPI_Comm_size (proccomm, &procglbnbr);          /* Get communicator data */
  MPI_Comm_rank (proccomm, &proclocnum);

#ifdef SCOTCH_CHECK_NOAUTO
  fprintf (stderr, "Proc %2d of %2d, pid %d\n", proclocnum, procglbnbr, getpid ());

  if (proclocnum == 0) {                          /* Synchronize on keybord input */
    char           c;

    printf ("Waiting for key press...\n");
    scanf ("%c", &c);
  }
#endif /* SCOTCH_CHECK_NOAUTO */

  if (MPI_Barrier (proccomm) != MPI_SUCCESS) {    /* Synchronize for debug */
    SCOTCH_errorPrint ("main: cannot communicate (1)");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_dgraphInit (&grafdat, proccomm) != 0) { /* Initialize source graph */
    SCOTCH_errorPrint ("main: cannot initialize graph");
    exit (EXIT_FAILURE);
  }

  file = NULL;
  if ((proclocnum == 0) &&
      ((file = fopen (argv[1], "r")) == NULL)) {
    SCOTCH_errorPrint ("main: cannot open graph file");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_dgraphLoad (&grafdat, file, -1, 0) != 0) {
    SCOTCH_errorPrint ("main: cannot load graph");
    exit (EXIT_FAILURE);
  }

  if (file != NULL)
    fclose (file);

  if (SCOTCH_dgraphCheck (&grafdat) != 0) {
    SCOTCH_errorPrint ("main: invalid graph");
    exit (EXIT_FAILURE);
  }

  if (MPI_Barrier (proccomm) != MPI_SUCCESS) {    /* Synchronize for debug */
    SCOTCH_errorPrint ("main: cannot communicate (2)");
    exit (EXIT_FAILURE);
  }

  SCOTCH_dgraphExit (&grafdat);

  MPI_Finalize ();
  exit (EXIT_SUCCESS);
}
