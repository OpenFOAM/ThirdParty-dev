/* Copyright 2024-2025 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : adm2dgr.c                               **/
/**                                                        **/
/**   AUTHOR     : Marc FUENTES                            **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Converts an auxiliary distributed mesh  **/
/**                file into a Scotch distributed graph    **/
/**                file that represents the dual element   **/
/**                graph of the mesh.                      **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 11 jan 2024     **/
/**                                 to   : 16 aug 2025     **/
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

#include "../libscotch/module.h"
#include "../libscotch/common.h"

#include "ptscotch.h"
#include "parmetis.h"                             /* Our "parmetis.h" file */
#include "adm2dgr.h"

/*
**  The static variables.
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* File array              */
                              { FILEMODER },
                              { FILEMODEW } };

static const char *         C_usageList[] = {
  "adm2dgr <input msh file> <out dgr file>",
  "  -c  : Set the minimum number of shared nodes to define adjacency",
  "  -h  : Display this help",
  "  -V  : Print program version and copyright",
  NULL };

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
  MPI_Comm            proccomm;
  int                 proclocnum;
  int                 procglbnbr;
  int                 protglbnum;
  SCOTCH_Dmesh        meshdat;
  SCOTCH_Dgraph       grafdat;
  SCOTCH_Num          baseval;
  SCOTCH_Num          noconbr;                    /* Number of common nodes                         */
  SCOTCH_Num          velmlocnbr;                 /* Local number of element vertices               */
  SCOTCH_Num *        velmloctab;                 /* Array of element vertex indices [velmloctab+1] */
  SCOTCH_Num *        eelmloctab;                 /* Array of element edge neighbors                */
  SCOTCH_Num *        prelvrttab;                 /* Array of element-to-process distribution       */
  SCOTCH_Num          preldspadj;                 /* Adjustment value for element-to-process array  */
  SCOTCH_Num          vertlocnbr;                 /* Number of local vertices in dual graph         */
  SCOTCH_Num *        vertloctab;                 /* Dual graph local vertex array                  */
  SCOTCH_Num          edgelocnbr;                 /* Number of local edges in dual graph            */
  SCOTCH_Num *        edgeloctab;                 /* Dual graph local edge array                    */
#ifdef SCOTCH_PTHREAD_MPI
  int                 thrdreqlvl;
  int                 thrdprolvl;
#endif /* SCOTCH_PTHREAD */
  int                 procnum;
  int                 i;

  errorProg ("adm2dgr");

#ifdef SCOTCH_PTHREAD_MPI
  thrdreqlvl = MPI_THREAD_MULTIPLE;
  if (MPI_Init_thread (&argc, &argv, thrdreqlvl, &thrdprolvl) != MPI_SUCCESS)
    SCOTCH_errorPrint ("main: Cannot initialize (1)");
  if (thrdreqlvl > thrdprolvl)
    SCOTCH_errorPrint ("main: MPI implementation is not thread-safe: recompile without SCOTCH_PTHREAD");
#else /* SCOTCH_PTHREAD_MPI */
  if (MPI_Init (&argc, &argv) != MPI_SUCCESS)
    SCOTCH_errorPrint ("main: Cannot initialize (2)");
#endif /* SCOTCH_PTHREAD_MPI */

  proccomm = MPI_COMM_WORLD;
  MPI_Comm_rank (proccomm, &proclocnum);
  MPI_Comm_size (proccomm, &procglbnbr);
  protglbnum = 0;
  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (EXIT_SUCCESS);
  }

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  noconbr = 1;
  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes */
    if ((argv[i][0] != '+') &&                    /* If found a file name      */
        ((argv[i][0] != '-') || (argv[i][1] == '\0'))) {
      if (C_fileNum < C_FILEARGNBR)               /* A file name has been given */
        fileBlockName (C_fileTab, C_fileNum ++) = argv[i];
      else
        SCOTCH_errorPrint ("main: too many file names given");
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (EXIT_SUCCESS);
        case 'C' :                                /* Number of common nodes */
        case 'c' :
          noconbr = atoi (&argv[i][2]);
          if ((noconbr < 1) || (noconbr > 3))
            SCOTCH_errorPrint ("main: invalid number of common nodes");
          break;
        case 'V' :
        case 'v' :
          fprintf (stderr, "adm2dgr, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, SCOTCH_COPYRIGHT_STRING "\n");
          fprintf (stderr, SCOTCH_LICENSE_STRING "\n");
          return  (EXIT_SUCCESS);
        default :
          SCOTCH_errorPrint ("main: unprocessed option '%s'", argv[i]);
      }
    }
  }

  fileBlockOpenDist (C_fileTab, C_FILENBR, procglbnbr, proclocnum, protglbnum); /* Open all files */

  if (SCOTCH_dmeshInit (&meshdat, proccomm) != 0)
    SCOTCH_errorPrint ("main: cannot initialize mesh");

  if (SCOTCH_dmeshLoad (&meshdat, C_filepntradminp, -1, 0) != 0)
    SCOTCH_errorPrint ("main: cannot read input mesh file");

  SCOTCH_dmeshData (&meshdat, &baseval, NULL, &velmlocnbr, &velmloctab, NULL, NULL, &eelmloctab, NULL, NULL);

  if ((prelvrttab = memAlloc ((procglbnbr + 1) * sizeof (SCOTCH_Num))) == NULL)
    errorPrint ("main: out of memory");
  if (MPI_Allgather (&velmlocnbr, 1, SCOTCH_NUM_MPI,
                     prelvrttab,  1, SCOTCH_NUM_MPI, proccomm) != MPI_SUCCESS) {
    errorPrint ("dmeshBuildAdm: communication error");
    return (1);
  }
  preldspadj = baseval;                           /* Use the new base to build displacement array */
  for (procnum = 0; procnum < procglbnbr; procnum ++) { /* Build element-to-process array         */
    SCOTCH_Num          preldspval;

    preldspval = prelvrttab[procnum];
    prelvrttab[procnum] = preldspadj;
    preldspadj += preldspval;
  }
  prelvrttab[procnum] = preldspadj;               /* Set end of displacement array */

  if (SCOTCH_ParMETIS_V3_Mesh2Dual (prelvrttab, velmloctab, eelmloctab, &baseval, &noconbr,
                                    &vertloctab, &edgeloctab, &proccomm) != METIS_OK)
    SCOTCH_errorPrint ("main: cannot compute dual graph");

  vertlocnbr = velmlocnbr;
  edgelocnbr = vertloctab[vertlocnbr] - baseval;

  SCOTCH_dgraphInit (&grafdat, MPI_COMM_WORLD);

  if (SCOTCH_dgraphBuild (&grafdat, baseval,
                          vertlocnbr, vertlocnbr, vertloctab, vertloctab + 1, NULL, NULL,
                          edgelocnbr, edgelocnbr, edgeloctab, NULL, NULL) != 0) {
    SCOTCH_errorPrint ("main: cannot build dual graph");
  }
  SCOTCH_dgraphCheck (&grafdat);
  SCOTCH_dgraphSave  (&grafdat, C_filepntrdgrout);
  SCOTCH_dgraphExit  (&grafdat);
  SCOTCH_dmeshExit   (&meshdat);

  free (prelvrttab);

  fileBlockClose (C_fileTab, C_FILENBR);

  MPI_Finalize ();

  return (EXIT_SUCCESS);
}
