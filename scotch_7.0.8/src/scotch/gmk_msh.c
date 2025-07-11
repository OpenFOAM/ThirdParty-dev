/* Copyright 2004,2007,2008,2010-2012,2014,2018,2019,2021,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : gmk_msh.c                               **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Marc FUENTES (v6.1)                     **/
/**                                                        **/
/**   FUNCTION   : Part of a mesh-to-graph converter.      **/
/**                This module contains the main function. **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 21 jan 2004     **/
/**                                 to   : 21 jan 2004     **/
/**                # Version 5.0  : from : 23 dec 2007     **/
/**                                 to   : 16 mar 2008     **/
/**                # Version 5.1  : from : 01 jul 2010     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 24 sep 2019     **/
/**                # Version 6.1  : from : 17 feb 2021     **/
/**                                 to   : 28 feb 2021     **/
/**                # Version 7.0  : from : 21 jan 2023     **/
/**                                 to   : 21 jan 2023     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "gmk_msh.h"

/*
**  The static and global variables.
*/

static int                  C_fileNum    = 0;     /* Number of file in arg list */
static File                 C_fileTab[2] = {      /* File array                 */
                              { FILEMODER },
                              { FILEMODEW } };

static const char *         C_usageList[] = {
  "gmk_msh [<input mesh file> [<output graph file>]] <options>",
  "  -d<noconbr>  : Build a dual graph instead of a primal graph, based on a minimum number of common neighbors (default: 1)",
  "  -h           : Display this help",
  "  -V           : Print program version and copyright",
  NULL };

/*****************************/
/*                           */
/* This is the main function */
/*                           */
/*****************************/

int
main (
int                         argc,
char *                      argv[])
{
  SCOTCH_Mesh         meshdat;
  SCOTCH_Graph        grafdat;
  int                 noconbr;                    /* Flag greater than zero if dual graph wanted */
  int                 i;

  errorProg ("gmk_msh");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (EXIT_SUCCESS);
  }

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  noconbr = 0;                                    /* No dual graph by default                         */
  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_fileNum < C_FILEARGNBR)               /* File name has been given                         */
        fileBlockName (C_fileTab, C_fileNum ++) = argv[i];
      else
        errorPrint ("main: too many file names given");
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'd' :                                /* Compute dual graph */
        case 'D' :
          if (argv[i][2] == '\0')                 /* If no argument, set value to 1 */
            noconbr = 1;
          else if ((noconbr = atoi (argv[i] + 2)) < 1) /* Get number of common points */
            errorPrint ("main: invalid number of common points '%s'", argv[i] + 2);
          break;
        case 'H' :                               /* Give help */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (EXIT_SUCCESS);
        case 'V' :
          fprintf (stderr, "gmk_msh, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, SCOTCH_COPYRIGHT_STRING "\n");
          fprintf (stderr, SCOTCH_LICENSE_STRING "\n");
          return  (EXIT_SUCCESS);
        default :
          errorPrint ("main: unprocessed option '%s'", argv[i]);
      }
    }
  }

  fileBlockOpen (C_fileTab, C_FILENBR);           /* Open all files */

  SCOTCH_meshInit  (&meshdat);
  SCOTCH_graphInit (&grafdat);

  SCOTCH_meshLoad  (&meshdat, C_filepntrmshinp, -1);
  SCOTCH_meshCheck (&meshdat);

  if (noconbr > 0)
    SCOTCH_meshGraphDual (&meshdat, &grafdat, noconbr);
  else
    SCOTCH_meshGraph (&meshdat, &grafdat);

#ifdef SCOTCH_DEBUG_ALL
  if (SCOTCH_graphCheck (&grafdat) != 0)
    errorPrint ("main: bad graph structure");
#endif /* SCOTCH_DEBUG_ALL */
  SCOTCH_graphSave (&grafdat, C_filepntrgrfout);

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  SCOTCH_graphExit (&grafdat);
  SCOTCH_meshExit  (&meshdat);

  return (EXIT_SUCCESS);
}
