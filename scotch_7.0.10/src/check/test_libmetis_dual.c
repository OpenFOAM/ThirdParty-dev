/* Copyright 2021,2024,2025 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_libmetis_dual.c                    **/
/**                                                        **/
/**   AUTHOR     : Marc FUENTES                            **/
/**                Francois PELLEGRINI                     **/
/**                Clement BARTHELEMY                      **/
/**                                                        **/
/**   FUNCTION   : This module tests the operations of     **/
/**                libscotchmetis dual graph routines.     **/
/**                                                        **/
/**   DATES      : # Version 6.1  : from : 10 feb 2021     **/
/**                                 to   : 17 jul 2021     **/
/**                # Version 7.0  : from : 08 aug 2024     **/
/**                                 to   : 04 jul 2025     **/
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
  SCOTCH_Num          vertnum;
  SCOTCH_Num          edgenum;
  SCOTCH_Num *        xadj;
  SCOTCH_Num *        adjncy;
  SCOTCH_Num          ncommon = 1;
  SCOTCH_Num          baseval = 0;
  SCOTCH_Num          xadnbr;
  SCOTCH_Num          adjnbr;
  SCOTCH_Num *        epart;
  SCOTCH_Num *        npart;
  SCOTCH_Num          nparts;
  SCOTCH_Num          objval;
  SCOTCH_Num          options[METIS_NOPTIONS];
  SCOTCH_Num          ne         = 6;
  SCOTCH_Num          nn         = 7;
  double              tpwgt[3]   = { 0.75, 0.125, 0.125 };
  SCOTCH_Num          vgwt[6]    = { 1, 2, 1, 1, 1, 1 };
  SCOTCH_Num          eptr[]     = { 0, 3, 6, 9, 12, 15, 18 };
  SCOTCH_Num          eind[]     = { 0, 1, 2, 0, 1, 5, 1, 5, 4, 1, 4, 6, 1, 6, 3, 1, 3, 2 };
  SCOTCH_Num          xadj_c[]   = { 0, 5, 10, 15, 20, 25, 30 };
  SCOTCH_Num          adjncy_c[] = { 1, 2, 3, 4, 5, 0, 2, 3, 4, 5, 0, 1, 3, 4, 5, 0, 1, 2, 4, 5, 0, 1, 2, 3, 5, 0, 1, 2, 3, 4 };

  SCOTCH_errorProg (argv[0]);

  xadnbr = sizeof (xadj_c) / sizeof (SCOTCH_Num);
  adjnbr = sizeof (adjncy_c) / sizeof (SCOTCH_Num);

  if (argc >= 2) {
    SCOTCH_Num          edgenbr;
    SCOTCH_Num          velmnum;

    baseval = atoi (argv[1]);
    edgenbr = eptr[ne];                           /* Preserve un-based number of edges */

    for (velmnum = 0; velmnum <= ne; velmnum ++)
      eptr[velmnum] += baseval;
    for (edgenum = 0; edgenum < edgenbr; edgenum ++)
      eind[edgenum] += baseval;
    for (vertnum = 0; vertnum < xadnbr; vertnum ++)
      xadj_c[vertnum] += baseval;
    for (edgenum = 0; edgenum < adjnbr; edgenum ++)
      adjncy_c[edgenum] += baseval;
  }

  SCOTCHMETISNAMEC (METIS_MeshToDual) (&ne, &nn, eptr, eind, &ncommon, &baseval, &xadj, &adjncy);
  if (xadj == NULL) {
    SCOTCH_errorPrint ("main: error in METIS_MeshToDual");
    exit (EXIT_FAILURE);
  }

  for (vertnum = 0; vertnum < xadnbr; vertnum ++) {
    if (xadj[vertnum] != xadj_c[vertnum]) {
      SCOTCH_errorPrint ("main: invalid vertex array");
      exit (EXIT_FAILURE);
    }
  }
  for (edgenum = 0; edgenum < adjnbr; edgenum ++) {
    if (adjncy[edgenum] != adjncy_c[edgenum]) {
      SCOTCH_errorPrint ("main: invalid edge array");
      exit (EXIT_FAILURE);
    }
  }

  free (xadj);                                    /* Free memory arrays allocated by METIS_MeshToDual() */
  free (adjncy);

  if (((epart = (SCOTCH_Num *) malloc (ne * sizeof (SCOTCH_Num))) == NULL) ||
      ((npart = (SCOTCH_Num *) malloc (nn * sizeof (SCOTCH_Num))) == NULL)) {
    SCOTCH_errorPrint ("main: out of memory");
    exit (EXIT_FAILURE);
  }

  nparts = 3;
  ncommon ++;
  SCOTCHMETISNAMEC (METIS_SetDefaultOptions) (options);
  options[METIS_OPTION_NUMBERING] = baseval;

  SCOTCHMETISNAMEC (METIS_PartMeshDual) (&ne, &nn, eptr, eind, vgwt, NULL, &ncommon, &nparts, tpwgt, options, &objval, epart, npart);
  if (objval < 0) {
    SCOTCH_errorPrint ("main: error in METIS_PartMeshDual");
    exit (EXIT_FAILURE);
  }

  free (npart);
  free (epart);

  exit (EXIT_SUCCESS);
}
