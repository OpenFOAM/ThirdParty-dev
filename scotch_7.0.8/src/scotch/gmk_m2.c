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
/**   NAME       : gmk_m2.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Creates the source graph for 2D mesh    **/
/**                graphs.                                 **/
/**                                                        **/
/**   DATES      : # Version 2.0  : from : 30 oct 1994     **/
/**                                 to   : 08 nov 1994     **/
/**                # Version 3.0  : from : 11 jul 1995     **/
/**                                 to   : 02 oct 1995     **/
/**                # Version 3.2  : from : 03 jun 1997     **/
/**                                 to   : 03 jun 1997     **/
/**                # Version 3.3  : from : 06 oct 1998     **/
/**                                 to   : 06 oct 1998     **/
/**                # Version 3.4  : from : 03 feb 2000     **/
/**                                 to   : 18 may 2004     **/
/**                # Version 5.0  : from : 13 dec 2007     **/
/**                                 to   : 16 mar 2008     **/
/**                # Version 5.1  : from : 01 jul 2010     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 24 sep 2019     **/
/**                # Version 6.1  : from : 02 apr 2021     **/
/**                                 to   : 02 apr 2021     **/
/**                # Version 7.0  : from : 02 apr 2021     **/
/**                                 to   : 21 jan 2023     **/
/**                                                        **/
/**   NOTES      : # The vertices of the (dX,dY) mesh are  **/
/**                  numbered as terminals so that         **/
/**                  t(0,0) = 0, t(1,0) = 1,               **/
/**                  t(dX - 1, 0) = dX - 1, t(0,1) = dX,   **/
/**                  and t(x,y) = (y * dX) + x.            **/
/**                                                        **/
/**                # The bit mask array represents the     **/
/**                  possible neighbors of a vertex, in    **/
/**                  the following order:                  **/
/**                         1     2     4                  **/
/**                         8   (16)   32                  **/
/**                        64   128   256                  **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "gmk_m2.h"

/*
**  The static definitions.
*/

static int                  C_paraNum = 0;        /* Number of parameters       */
static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* The file array          */
                              { FILEMODEW },
                              { FILEMODEW } };

static const char *         C_usageList[] = {
  "gmk_m2 <dimX> [<dimY> [<output source file>]] <options>",
  "  -b<val>   : Set base value for output (0 or 1)",
  "  -e        : Build a 8-neighbor grid rather than a 4-neighbor one",
  "  -g<file>  : Output the geometry to <file>",
  "  -h        : Display this help",
  "  -t        : Build a torus rather than a grid",
  "  -V        : Print program version and copyright",
  "  -y        : Invert y coordinate in geometry",
  NULL };

/****************************************/
/*                                      */
/* The main routine, which computes the */
/* source graph description.            */
/*                                      */
/****************************************/

int
main (
int                         argc,
char *                      argv[])
{
  int                 flagval;                    /* Process flags          */
  SCOTCH_Num          baseval;                    /* Base value             */
  SCOTCH_Num          edgenbr;                    /* Number of edges (arcs) */
  SCOTCH_Num          d[2] = { 1, 1 };            /* Mesh dimensions        */
  SCOTCH_Num          c[2];                       /* Vertex coordinates     */
  int                 cbitmsk;                    /* Coordinate flag mask   */
  int                 i;

  errorProg ("gmk_m2");

  flagval = C_FLAGDEFAULT;                        /* Set default flags */
  baseval = 0;
  cbitmsk = 2 + 8 + 32 + 128;                     /* 4-point mask */

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (EXIT_SUCCESS);
  }

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_paraNum < 2) {                        /* If number of parameters not reached              */
        if ((d[C_paraNum ++] = atoi (argv[i])) < 1) /* Get the dimension                              */
          errorPrint ("main: invalid dimension '%s'", argv[i]);
        continue;                                 /* Process the other parameters */
      }
      if (C_fileNum < C_FILEARGNBR)               /* A file name has been given */
        fileBlockName (C_fileTab, C_fileNum ++) = argv[i];
      else
        errorPrint ("main: too many file names given");
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'B' :                                /* Set base value */
        case 'b' :
          baseval = (SCOTCH_Num) atol (&argv[i][2]);
          if ((baseval < 0) || (baseval > 1)) {
            errorPrint ("main: invalid base value '" SCOTCH_NUMSTRING "'", (SCOTCH_Num) baseval);
          }
          break;
        case 'E' :                                /* Build a finite-element grid */
        case 'e' :
          flagval |= C_FLAGELEM;
          cbitmsk  = 1 + 2 + 4 + 8 + 32 + 64 + 128 + 256; /* 8-point mask */
          break;
        case 'G' :                                /* Output the geometry */
        case 'g' :
          flagval |= C_FLAGGEOOUT;
          if (argv[i][2] != '\0')
            C_filenamegeoout = &argv[i][2];
          break;
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (EXIT_SUCCESS);
        case 'T' :                                /* Build a torus */
        case 't' :
          flagval |= C_FLAGTORUS;
          break;
        case 'V' :
          fprintf (stderr, "gmk_m2, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, SCOTCH_COPYRIGHT_STRING "\n");
          fprintf (stderr, SCOTCH_LICENSE_STRING "\n");
          return  (EXIT_SUCCESS);
        case 'Y' :                                /* Invert the y coordinate in geometry file */
        case 'y' :
          flagval |= C_FLAGGEOINVY;
          break;
        default :
          errorPrint ("main: unprocessed option '%s'", argv[i]);
      }
    }
  }

  fileBlockOpen (C_fileTab, C_FILENBR);           /* Open all files */

  if ((flagval & C_FLAGTORUS) != 0) {             /* Build a torus           */
    edgenbr = ((flagval & C_FLAGELEM) != 0)       /* Compute number of edges */
              ? (8 * (d[0] * d[1]) - ((d[0] < 3) ? (6 * d[1]) : 0) - ((d[1] < 3) ? (6 * d[0]) : 0))
              : (4 * (d[0] * d[1]) - ((d[0] < 3) ? (2 * d[1]) : 0) - ((d[1] < 3) ? (2 * d[0]) : 0));

    fprintf (C_filepntrsrcout, "0\n" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\n" SCOTCH_NUMSTRING "\t000\n",
             (SCOTCH_Num) (d[0] * d[1]),          /* Print number of vertices */
             (SCOTCH_Num) edgenbr,
             (SCOTCH_Num) baseval);

    for (c[1] = 0; c[1] < d[1]; c[1] ++) {        /* Output neighbor list */
      for (c[0] = 0; c[0] < d[0]; c[0] ++) {
        int                 cbitval;

        cbitval = cbitmsk;
	if (d[0] <= 2) {                          /* Remove ascending arcs if duplicates */
          cbitval &= ~(4 + 32 + 256);
          if (d[0] <= 1)                          /* Remove descending arcs if no width */
            cbitval &= ~(1 + 8 + 64);
        }
        if (d[1] <= 2) {
          cbitval &= ~(64 + 128 + 256);
          if (d[1] <= 1)
            cbitval &= ~(1 + 2 + 4);
        }

        C_edgeOutput (d, c, baseval, cbitval, C_filepntrsrcout);
      }
    }
  }
  else {                                          /* Build a grid            */
    edgenbr = ((flagval & C_FLAGELEM) != 0)       /* Compute number of edges */
              ? (8 * (d[0] * d[1]) - 6 * (d[0] + d[1]) + 4)
              : (4 * (d[0] * d[1]) - 2 * (d[0] + d[1]));

    fprintf (C_filepntrsrcout, "0\n" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\n" SCOTCH_NUMSTRING "\t000\n",
             (SCOTCH_Num) (d[0] * d[1]),          /* Print number of vertices */
             (SCOTCH_Num) edgenbr,
             (SCOTCH_Num) baseval);

    for (c[1] = 0; c[1] < d[1]; c[1] ++) {        /* Output neighbor list */
      for (c[0] = 0; c[0] < d[0]; c[0] ++) {
        int                 cbitval;

        cbitval = cbitmsk;
        if (c[0] <= 0)
          cbitval &= ~(1 + 8 + 64);
        if (c[0] >= (d[0] - 1))
          cbitval &= ~(4 + 32 + 256);
        if (c[1] <= 0)
          cbitval &= ~(1 + 2 + 4);
        if (c[1] >= (d[1] - 1))
          cbitval &= ~(64 + 128 + 256);

        C_edgeOutput (d, c, baseval, cbitval, C_filepntrsrcout);
      }
    }
  }

  if (flagval & C_FLAGGEOOUT) {                   /* If geometry is wanted                */
   fprintf (C_filepntrgeoout, "2\n" SCOTCH_NUMSTRING "\n", /* Output geometry file header */
            (SCOTCH_Num) (d[0] * d[1]));

    for (c[1] = 0; c[1] < d[1]; c[1] ++) {        /* Output mesh coordinates */
      for (c[0] = 0; c[0] < d[0]; c[0] ++)
        fprintf (C_filepntrgeoout, SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\n",
                 (SCOTCH_Num) (c[1] * d[0] + c[0] + baseval),
                 (SCOTCH_Num) c[0],
                 (SCOTCH_Num) (((flagval & C_FLAGGEOINVY) != 0) ? (d[1] - 1 - c[1]) : c[1]));
    }
  }

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  return (EXIT_SUCCESS);
}

/* This routine writes the neighbors of a 2D grid vertex.
** It returns:
** - void  : in all cases.
*/

void
C_edgeOutput (
const SCOTCH_Num * const    d,
const SCOTCH_Num * const    c,
const SCOTCH_Num            baseval,
const int                   cbitval,
FILE * const                fileptr)
{
  SCOTCH_Num          n[2];                       /* Neighbor loop indices */
  int                 degrval;
  int                 cbitnum;

  for (cbitnum = 256, degrval = 0; cbitnum != 0; cbitnum >>= 1) { /* Count effective number of neighbors */
    if ((cbitval & cbitnum) != 0)
      degrval ++;
  }

  fprintf (fileptr, "%d", degrval);               /* Write number of neighbors */

  for (n[1] = c[1] - 1, cbitnum = 1; n[1] <= c[1] + 1; n[1] ++) { /* Output the effective neighbors */
    for (n[0] = c[0] - 1; n[0] <= c[0] + 1; n[0] ++, cbitnum <<= 1) {
      if ((cbitval & cbitnum) != 0)
        fprintf (fileptr, "\t" SCOTCH_NUMSTRING,
                 (SCOTCH_Num) (((n[1] + d[1]) % d[1]) * d[0] + ((n[0] + d[0]) % d[0]) + baseval));
    }
  }

  fprintf (fileptr, "\n");
}
