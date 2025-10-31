/* Copyright 2004,2007,2008,2010-2012,2014,2018-2020,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : gotst.c                                 **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Bruno MARCUSSEAU (v3.1)                 **/
/**                                                        **/
/**   FUNCTION   : Graph symbolic factorizer.              **/
/**                This module contains the main function. **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 27 jan 2004     **/
/**                                 to   : 28 nov 2005     **/
/**                # Version 5.0  : from : 25 jun 2007     **/
/**                                 to   : 16 mar 2008     **/
/**                # Version 5.1  : from : 20 apr 2010     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 24 sep 2019     **/
/**                # Version 6.1  : from : 01 jun 2020     **/
/**                                 to   : 14 jun 2020     **/
/**                # Version 7.0  : from : 21 jan 2023     **/
/**                                 to   : 21 jan 2023     **/
/**                                                        **/
/**   NOTES      : # The cost analysis routine leaves the  **/
/**                  memory management to malloc and free  **/
/**                  because it is assumed to be the most  **/
/**                  efficient to merge free blocks and    **/
/**                  reallocate them.                      **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "gotst.h"

/*
**  The static and global variables.
*/

static int                  C_fileNum    = 0;     /* Number of file in arg list  */
static File                 C_fileTab[3] = {      /* File array                  */
                              { FILEMODER },
                              { FILEMODER },
                              { FILEMODEW } };

static const char *         C_usageList[] = {
  "gotst [<input graph file> [<input ordering file> [<output data file>]]] <options>",
  "  -h  : Display this help",
  "  -v  : Do not account for vertex weights",
  "  -V  : Print program version and copyright",
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
  SCOTCH_Graph        grafdat;
  SCOTCH_Num          vertnbr;
  SCOTCH_Num *        verttab;
  SCOTCH_Num *        vendtab;
  SCOTCH_Num *        velotab;
  SCOTCH_Num          edgenbr;
  SCOTCH_Num *        edgetab;
  SCOTCH_Num          baseval;
  SCOTCH_Ordering     ordedat;
  SCOTCH_Num *        permtab;
  SCOTCH_Num *        peritab;
  int                 flagval;
  int                 i;

  flagval = C_FLAGNONE;                           /* Default behavior */

  errorProg ("gotst");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (EXIT_SUCCESS);
  }

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_fileNum < C_FILEARGNBR)               /* File name has been given                         */
        fileBlockName (C_fileTab, C_fileNum ++) = argv[i];
      else
        errorPrint ("main: too many file names given");
    }
    else {                                       /* If found an option name */
      switch (argv[i][1]) {
        case 'H' :                               /* Give help */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (EXIT_SUCCESS);
        case 'v' :
          flagval |= C_FLAGNOVERT;
          break;
        case 'V' :
          fprintf (stderr, "gotst, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, SCOTCH_COPYRIGHT_STRING "\n");
          fprintf (stderr, SCOTCH_LICENSE_STRING "\n");
          return  (EXIT_SUCCESS);
        default :
          errorPrint ("main: unprocessed option '%s'", argv[i]);
      }
    }
  }

  fileBlockOpen (C_fileTab, C_FILENBR);           /* Open all files */

  SCOTCH_graphInit (&grafdat);
  SCOTCH_graphLoad (&grafdat, C_filepntrgrfinp, -1, ((flagval & C_FLAGNOVERT) != 0) ? 3 : 2); /* Read source graph with or without vertex weights */
  SCOTCH_graphData (&grafdat, &baseval, &vertnbr, &verttab, &vendtab, &velotab, NULL, &edgenbr, &edgetab, NULL);
#ifdef SCOTCH_DEBUG_ALL
  if (vendtab != (verttab + 1))
    errorPrint ("main: graph should be compact");
#endif /* SCOTCH_DEBUG_ALL */
  if (memAllocGroup ((void **) (void *)
                     &peritab, (size_t) (vertnbr * sizeof (SCOTCH_Num)),
                     &permtab, (size_t) (vertnbr * sizeof (SCOTCH_Num)), NULL) == NULL)
    errorPrint ("main: out of memory");

  SCOTCH_graphOrderInit (&grafdat, &ordedat, permtab, peritab, NULL, NULL, NULL);
  SCOTCH_graphOrderLoad (&grafdat, &ordedat, C_filepntrordinp);
  if (SCOTCH_graphOrderCheck (&grafdat, &ordedat) != 0)
    errorPrint ("main: invalid ordering");

  factorView (baseval, vertnbr, verttab, velotab, edgenbr, edgetab, permtab, peritab, C_filepntrdatout);

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  memFree               (peritab);
  SCOTCH_graphOrderExit (&grafdat, &ordedat);
  SCOTCH_graphExit      (&grafdat);

  return (EXIT_SUCCESS);
}

/*************************************/
/*                                   */
/* These routines compute statistics */
/* on orderings.                     */
/*                                   */
/*************************************/

static
int
factorView (
const SCOTCH_Num              baseval,
const SCOTCH_Num              vertnbr,
const SCOTCH_Num * const      verttab,
const SCOTCH_Num * const      velotab,
const SCOTCH_Num              edgenbr,
const SCOTCH_Num * const      edgetab,
const SCOTCH_Num * const      permtab,
const SCOTCH_Num * const      peritab,
FILE * restrict const         stream)
{
  FactorStat          statdat;                    /* Structure holding data for performing elimination tree statistics              */
  SCOTCH_Num *        ldadtab;                    /* Array of permuted indices of parent column in factored matrix                  */
  SCOTCH_Num *        lsontab;                    /* Array of permuted indices of first child column in factored matrix             */
  SCOTCH_Num *        lbrotab;                    /* Array of permuted indices of sibling column in factored matrix                 */
  SCOTCH_Num *        fnnztab;                    /* Array of number of row indices, then of weighted factors, in factored matrix   */
  double              fopcsum;                    /* Sum of OPC values for all independent elimination trees in factored matrix     */
  double              heigsum;                    /* Sum of branch heights for all independent elimination trees in factored matrix */
  SCOTCH_Num          fcolnum;                    /* Index of permuted column in factored matrix                                    */
  int                 o;

  if (memAllocGroup ((void **) (void *)
                     &ldadtab, (size_t) (vertnbr * sizeof (SCOTCH_Num)),
                     &lsontab, (size_t) (vertnbr * sizeof (SCOTCH_Num)),
                     &lbrotab, (size_t) (vertnbr * sizeof (SCOTCH_Num)),
                     &fnnztab, (size_t) (vertnbr * sizeof (SCOTCH_Num)), NULL) == NULL) {
    errorPrint ("factorView: out of memory");
    return     (1);
  }
  statdat.velotax = (velotab != NULL) ? (velotab - baseval) : NULL;
  statdat.peritax = peritab - baseval;
  statdat.ldadtax = ldadtab - baseval;
  statdat.lsontax = lsontab - baseval;
  statdat.lbrotax = lbrotab - baseval;
  statdat.fnnztax = fnnztab - baseval;

  if (factorView2 (baseval, vertnbr, verttab - baseval, statdat.velotax,
                   edgetab - baseval, permtab - baseval, peritab - baseval,
                   ldadtab - baseval, lsontab - baseval, lbrotab - baseval, fnnztab - baseval) != 0) {
    errorPrint ("factorView: factored matrix too large");
    memFree    (ldadtab);                         /* Free group leader */
    return     (1);
  }

  statdat.heigmin = SCOTCH_NUMMAX;
  statdat.heigmax =
  statdat.heignbr = 0;
  heigsum         = 0.0L;
  for (fcolnum = 0; fcolnum < vertnbr; fcolnum ++) { /* Get height sum        */
    if (ldadtab[fcolnum] == -1)                   /* If column is a root      */
      factorView3 (&statdat, 1, fcolnum + baseval, &heigsum); /* Scan subtree */
  }
  statdat.heigavg = heigsum / (double) statdat.heignbr;
  statdat.heigdlt = 0.0L;
  statdat.fnnzsum = 0.0L;
  fopcsum         = 0.0L;
  for (fcolnum = 0; fcolnum < vertnbr; fcolnum ++) { /* Get delta        */
    if (ldadtab[fcolnum] == -1)                   /* If column is a root */
      factorView4 (&statdat, 1, fcolnum + baseval, &fopcsum);
  }
  statdat.heigdlt /= (double) statdat.heignbr;

  o = (fprintf (stream, "O\tLeaf=" SCOTCH_NUMSTRING "\nO\tHeight min=" SCOTCH_NUMSTRING "\tmax=" SCOTCH_NUMSTRING "\tavg=%f\tdlt=%f (%5.2f)\n", /* Write tree height statistics */
                (SCOTCH_Num) statdat.heignbr, (SCOTCH_Num) statdat.heigmin, (SCOTCH_Num) statdat.heigmax,
                statdat.heigavg, statdat.heigdlt, ((statdat.heigdlt / statdat.heigavg) * (double) 100.0L)) == EOF);

  o |= (fprintf (stream, "O\tNNZ=%e\nO\tOPC=%e\n",
                 statdat.fnnzsum,
                 fopcsum) == EOF);

  if (o != 0)
    errorPrint ("factorView: bad output");

  memFree (ldadtab);                              /* Free group leader */

  return (o);
}

/* This routine performs the symbolic factorization
** and creates the elimination tree structure. Also,
** it computes the (possibly weighted) sum of factors
** for each permuted column.
*/

static
int
factorView2 (
const SCOTCH_Num              baseval,
const SCOTCH_Num              vertnbr,
const SCOTCH_Num * const      verttax,
const SCOTCH_Num * const      velotax,
const SCOTCH_Num * const      edgetax,
const SCOTCH_Num * const      permtax,
const SCOTCH_Num * const      peritax,
SCOTCH_Num * restrict         ldadtax,
SCOTCH_Num * restrict         lsontax,
SCOTCH_Num * restrict         lbrotax,
SCOTCH_Num * restrict         fnnztax)
{
  SCOTCH_Num * restrict         frowtab;          /* Array for sorting permuted row indices in column being factored */
  SCOTCH_Num * restrict         fnxttab;          /* Index array for chaining permuted row indices in frowtab        */
  SCOTCH_Num ** restrict        facttax;          /* Array of row index arrays for all factored colums to date       */
  SCOTCH_Num                    vertnnd;
  SCOTCH_Num                    pcolnum;

  memSet (lsontax + baseval, ~0, vertnbr * sizeof (SCOTCH_Num)); /* Assume columns have no sons at all */

  if (memAllocGroup ((void **) (void *)
                     &frowtab, (size_t) ((vertnbr + 1) * sizeof (SCOTCH_Num)),
                     &fnxttab, (size_t) ((vertnbr + 1) * sizeof (SCOTCH_Num)),
                     &facttax, (size_t) (vertnbr       * sizeof (SCOTCH_Num *)), NULL) == NULL) {
    errorPrint ("factorView2: out of memory (1)");
    return     (1);
  }
  memSet (facttax, 0, vertnbr * sizeof (SCOTCH_Num *)); /* Set all factored column pointers to NULL */
  facttax -= baseval;

  vertnnd = vertnbr + baseval;
  for (pcolnum = baseval; pcolnum < vertnnd; pcolnum ++) { /* For all columns of the permuted matrix */
    SCOTCH_Num *        fcoltab;
    SCOTCH_Num *        fcolptr;
    SCOTCH_Num          frownbr;
    SCOTCH_Num          frowidx;
    SCOTCH_Num          frowidd;
    SCOTCH_Num          scolnum;
    SCOTCH_Num          icolnum;                  /* Based non-permuted column index */
    SCOTCH_Num          edgenum;
    SCOTCH_Num          edgennd;
    SCOTCH_Num          dcolnum;

    icolnum = peritax[pcolnum];                   /* Get the un-permuted index of the column */

    frownbr    = 1;                               /* Start array of factored terms for column   */
    frowtab[0] = pcolnum;                         /* Add diagonal element as unmoveable starter */
    for (edgenum = verttax[icolnum], edgennd = verttax[icolnum + 1]; edgenum < edgennd; edgenum ++) {
      SCOTCH_Num          prownum;

      prownum = permtax[edgetax[edgenum]];        /* Get permuted row */

      if (prownum >= pcolnum)                     /* Include only permuted row indices of higher rank (lower triangular matrix) */
        frowtab[frownbr ++] = prownum;
    }
    intSort1asc1 (frowtab + 1, frownbr - 1);      /* Sort rows in ascending order */

    frowtab[frownbr ++] = vertnnd;                /* Add trailer                   */
    for (frowidx = 0; frowidx < (frownbr - 1); frowidx ++) /* Create initial links */
      fnxttab[frowidx] = frowidx + 1;
    frowidd = frowidx;                            /* Save index of trailer */

    for (scolnum = lsontax[pcolnum]; scolnum != -1; scolnum = lbrotax[scolnum]) { /* For all son columns in elimination tree */
      const SCOTCH_Num * restrict srowtab;        /* Array of factored indices for son column                                */
      SCOTCH_Num                  srownbr;
      SCOTCH_Num                  srowidx;
      SCOTCH_Num                  frowidx;        /* Index of current factor in singly-linked list        */
      SCOTCH_Num                  foldidx;        /* Index of before-current factor in singly-linked list */
      SCOTCH_Num                  frownum;        /* Value of current factor in singly-linked list        */
      SCOTCH_Num                  fnnzsum;

      srowtab = facttax[scolnum];                 /* Point to array of factors for son column                    */
      srownbr = fnnztax[scolnum];                 /* Get size of array                                           */
      if (velotax != NULL)                        /* If needed, accumulate weight of first extra-diagonal factor */
        fnnzsum = velotax[peritax[srowtab[0]]];

      for (srowidx = 1, frowidx = 0, foldidx = -1, frownum = frowtab[frowidx];
           srowidx < srownbr; srowidx ++) {
        SCOTCH_Num                  srownum;

        srownum = srowtab[srowidx];
        if (velotax != NULL)                      /* If needed, accumulate loads before freeing factored son column */
          fnnzsum += velotax[peritax[srownum]];

        while (frownum < srownum) {               /* While factor to add not in place */
          foldidx = frowidx;                      /* Skip to next position            */
          frowidx = fnxttab[frowidx];
          frownum = frowtab[frowidx];
        }
        if (srownum == frownum)                   /* If factor already in column, skip it */
          continue;

        frowtab[frownbr] = srownum;               /* Add new slot  */
        fnxttab[frownbr] = frowidx;               /* Link new slot */
        fnxttab[foldidx] = frownbr;
        foldidx          = frownbr ++;
      }

      if (velotax != NULL)                        /* If needed, accumulate sum of extra-diagonal factors in column */
        fnnztax[scolnum] = fnnzsum;

      memFree (srowtab);                          /* Free now useless son factored column */
#ifdef SCOTCH_DEBUG_ALL
      facttax[scolnum] = NULL;
#endif /* SCOTCH_DEBUG_ALL */
    }

    frownbr -= 2;                                 /* Remove markers from number of extra-diagonals */
    fnnztax[pcolnum] = frownbr;                   /* Save number of extra-diagonals                */

    if (frownbr <= 0) {                           /* If factored column has no extra-diagonals */
      ldadtax[pcolnum] = -1;                      /* Column has no father                      */
#ifdef SCOTCH_DEBUG_ALL
      lbrotax[pcolnum] = -1;
#endif /* SCOTCH_DEBUG_ALL */
      continue;                                   /* Skip to next column without allocating or linking */
    }

    if ((fcoltab = memAlloc (frownbr * sizeof (SCOTCH_Num))) == NULL) { /* Allocate array for factored column */
      errorPrint ("factorView2: out of memory (2)");
      return     (1);
    }

    for (frowidx = fnxttab[0], fcolptr = fcoltab; frowidx != frowidd; frowidx = fnxttab[frowidx]) /* Fill factored array for column */
      *fcolptr ++ = frowtab[frowidx];

    dcolnum = fcoltab[0];                         /* Get number of father, that it, first extra-diagonal */
    ldadtax[pcolnum] = dcolnum;                   /* Link factored column to the separation tree         */
    lbrotax[pcolnum] = lsontax[dcolnum];
    lsontax[dcolnum] = pcolnum;

    facttax[pcolnum] = fcoltab;                   /* Save factored array */
  }

  memFree (frowtab);                              /* Free group leader */

  return (0);
}

/* This routine traverses the elimination tree from the root(s),
** so as to compute minimum, maximum, and sum of tree branches.
*/

static
void
factorView3 (
FactorStat * restrict const   statptr,
SCOTCH_Num                    levlnum,
SCOTCH_Num                    fcolnum,
double * restrict const       hsumptr)
{
  double                  hsumtmp;

  hsumtmp = 0.0;
  if (statptr->lsontax[fcolnum] != -1) {          /* If node has descendants */
    SCOTCH_Num              csonnum;

    for (csonnum = statptr->lsontax[fcolnum]; csonnum != -1; csonnum = statptr->lbrotax[csonnum])
      factorView3 (statptr, levlnum + 1, csonnum, &hsumtmp);
  }
  else {
    hsumtmp = (double) levlnum;

    statptr->heignbr ++;
    if (levlnum < statptr->heigmin)
      statptr->heigmin = levlnum;
    if (levlnum > statptr->heigmax)
      statptr->heigmax = levlnum;
  }

  *hsumptr += hsumtmp;
}
/* This routine traverses the elimination tree from the root(s),
** so as to compute minimum, maximum, and sum of tree branches.
*/

static
void
factorView4 (
FactorStat * restrict const   statptr,
SCOTCH_Num                    levlnum,
SCOTCH_Num                    fcolnum,
double * restrict const       fopcptr)
{
  SCOTCH_Num              veloval;
  SCOTCH_Num              fnnztmp;
  double                  fopctmp;

  veloval = (statptr->velotax == NULL) ? 1 : statptr->velotax[statptr->peritax[fcolnum]];
  fnnztmp = statptr->fnnztax[fcolnum] + veloval;  /* Get extra-diagonals, plus diagonal */
  fopctmp = (double) (fnnztmp * veloval);
  statptr->fnnzsum += fopctmp;
  fopctmp *= (double) fnnztmp;

  if (statptr->lsontax[fcolnum] != -1) {          /* If node has descendants */
    SCOTCH_Num              csonnum;

    for (csonnum = statptr->lsontax[fcolnum]; csonnum != -1; csonnum = statptr->lbrotax[csonnum])
      factorView4 (statptr, levlnum + 1, csonnum, &fopctmp); /* Accumulate OPC on local sum */
  }
  else
    statptr->heigdlt += fabs ((double) levlnum - statptr->heigavg);

  *fopcptr += fopctmp;                            /* Aggregate local sum at higher level */
}
