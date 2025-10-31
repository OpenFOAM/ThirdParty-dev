/* Copyright 2004,2007,2009,2012,2015,2018,2020,2022,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : main_mumps.c                            **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This is the test module for the MUMPS   **/
/**                interface routine.                      **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 17 may 2001     **/
/**                                 to   : 17 may 2001     **/
/**                # Version 1.0  : from : 17 jun 2005     **/
/**                                 to   : 17 jun 2005     **/
/**                # Version 5.1  : from : 22 jan 2009     **/
/**                                 to   : 22 jan 2009     **/
/**                # Version 6.0  : from : 01 dec 2012     **/
/**                                 to   : 22 jan 2020     **/
/**                # Version 7.0  : from : 10 dec 2022     **/
/**                                 to   : 10 sep 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "graph.h"
#include "esmumps.h"

void                        ESMUMPSF            (const INT * const, const INT * const, INT * const, const INT * const, INT * const, INT * const, INT * const, INT * const, INT * const, INT * const);

/******************************/
/*                            */
/* This is the main function. */
/*                            */
/******************************/

int
main (
int                         argc,
char *                      argv[])
{
  Graph               grafdat;                    /* Graph to load */
  INT                 vertnbr;
  INT *               verttab;
  INT                 edgenbr;
  INT *               edgetab;
  INT *               lentab;
  INT *               nvtab;
  INT *               elentab;
  INT *               lasttab;
  INT                 pfree;
  INT                 ncmpa;
  INT                 vertnum;
  FILE *              stream;

  if (argc != 2) {
    errorPrint ("test_esmumps: usage: test_esmumps graph_file");
    exit       (EXIT_FAILURE);
  }

  graphInit (&grafdat);
  if ((stream = fopen (argv[1], "r")) == NULL) {
    errorPrint ("test_esmumps: cannot open graph file");
    graphExit  (&grafdat);
    exit       (EXIT_FAILURE);
  }
  if (graphLoad (&grafdat, stream, 1, 3) != 0) {  /* Base graph with base value 1, no loads */
    errorPrint ("test_esmumps: cannot load graph file");
    graphExit  (&grafdat);
    exit       (EXIT_FAILURE);
  }
  fclose (stream);

  graphData (&grafdat, NULL, &vertnbr, &verttab, NULL, NULL, NULL, &edgenbr, &edgetab, NULL);

  nvtab   =                                       /* Assume an error */
  elentab =
  lasttab = NULL;
  if (((lentab  = malloc (vertnbr * sizeof (INT))) == NULL) ||
      ((nvtab   = malloc (vertnbr * sizeof (INT))) == NULL) ||
      ((elentab = malloc (vertnbr * sizeof (INT))) == NULL) ||
      ((lasttab = malloc (vertnbr * sizeof (INT))) == NULL)) {
    errorPrint ("test_esmumps: out of memory");
    free       (lentab);
    free       (nvtab);
    free       (elentab);
    free       (lasttab);
    graphExit  (&grafdat);
    exit       (EXIT_FAILURE);
  }

  for (vertnum = 0; vertnum < vertnbr; vertnum ++) {
    if (verttab[vertnum] == verttab[vertnum + 1]) {
      lentab[vertnum] = 0;
      verttab[vertnum] = 0;                       /* Graph structure no longer valid in Emilio */
    }
    else
      lentab[vertnum] = verttab[vertnum + 1] - verttab[vertnum];
  }

  pfree = edgenbr + 1;
  ESMUMPSF (&vertnbr, &edgenbr, verttab, &pfree,
            lentab, edgetab, nvtab, elentab, lasttab, &ncmpa);

  free      (lentab);
  free      (nvtab);
  free      (elentab);
  free      (lasttab);
  graphExit (&grafdat);

  if (ncmpa < 0) {
    errorPrint ("test_esmumps: error in ESMUMPSF (%d)", ncmpa);
    exit       (EXIT_FAILURE);
  }

  exit (EXIT_SUCCESS);
}
