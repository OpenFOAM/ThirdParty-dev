/* Copyright 2007-2012,2018,2019,2021,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : metis_graph_part.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Amaury JACQUES (v6.0)                   **/
/**                                                        **/
/**   FUNCTION   : This module is the compatibility        **/
/**                library for the MeTiS partitioning      **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 08 sep 2006     **/
/**                                 to   : 07 jun 2007     **/
/**                # Version 5.1  : from : 06 jun 2009     **/
/**                                 to   : 30 jun 2010     **/
/**                # Version 6.0  : from : 23 dec 2011     **/
/**                                 to   : 19 aug 2019     **/
/**                # Version 6.1  : from : 20 jun 2021     **/
/**                                 to   : 30 dec 2021     **/
/**                # Version 7.0  : from : 21 jan 2023     **/
/**                                 to   : 11 aug 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "metis.h"                                /* Our "metis.h" file */
#include "metis_graph_part.h"

/********************************/
/*                              */
/* These routines compute MeTiS */
/* mapping statistics.          */
/*                              */
/********************************/

/* This routine computes the cut cost of the
** given partition.
** It returns:
** - METIS_OK       : if the metrics could be computed.
** - METIS_ERROR_*  : on error.
*/

int
_SCOTCH_METIS_OutputCut (
const SCOTCH_Num            baseval,
const SCOTCH_Num            vertnnd,
const SCOTCH_Num * const    verttax,
const SCOTCH_Num * const    edgetax,
const SCOTCH_Num * const    edlotax,
const SCOTCH_Num * const    parttax,
SCOTCH_Num * const          csizptr)
{
  SCOTCH_Num          csizsum;

  csizsum = 0;
  if (edlotax == NULL) {                          /* If graph does not have edge weights */
    SCOTCH_Num          vertnum;
    SCOTCH_Num          edgenum;

    for (vertnum = edgenum = baseval; vertnum < vertnnd; vertnum ++) {
      SCOTCH_Num          edgennd;
      SCOTCH_Num          partval;

      partval = parttax[vertnum];
      for (edgennd = verttax[vertnum + 1]; edgenum < edgennd; edgenum ++) {
        if (parttax[edgetax[edgenum]] != partval)
          csizsum ++;
      }
    }
  }
  else {                                          /* Graph has edge weights */
    SCOTCH_Num          vertnum;
    SCOTCH_Num          edgenum;

    for (vertnum = edgenum = baseval; vertnum < vertnnd; vertnum ++) {
      SCOTCH_Num          edgennd;
      SCOTCH_Num          partval;

      partval = parttax[vertnum];
      for (edgennd = verttax[vertnum + 1]; edgenum < edgennd; edgenum ++) {
        SCOTCH_Num          vertend;

        vertend = edgetax[edgenum];
        if (parttax[vertend] != partval)
          csizsum += edlotax[edgenum];
      }
    }
  }

  *csizptr = csizsum / 2;

  return (METIS_OK);
}

/* This routine computes the cost of the
** given partition.
** It returns:
** - METIS_OK       : if the metrics could be computed.
** - METIS_ERROR_*  : on error.
*/

int
_SCOTCH_METIS_OutputVol (
const SCOTCH_Num            baseval,
const SCOTCH_Num            vertnnd,
const SCOTCH_Num * const    verttax,
const SCOTCH_Num * const    edgetax,
const SCOTCH_Num * const    vsiztax,
const SCOTCH_Num            partnbr,
const SCOTCH_Num * const    parttax,
SCOTCH_Num * const          cvolptr)
{
  SCOTCH_Num * restrict nghbtax;
  SCOTCH_Num            vsizval;                  /* Communication volume of current vertex */
  SCOTCH_Num            vertnum;
  SCOTCH_Num            edgenum;
  SCOTCH_Num            cvolsum;

  if ((nghbtax = memAlloc (partnbr * sizeof (SCOTCH_Num))) == NULL) /* Part array is un-based at allocation */
    return (METIS_ERROR_MEMORY);

  memSet (nghbtax, ~0, partnbr * sizeof (SCOTCH_Num));
  nghbtax -= baseval;                             /* Base part array since parts are based in MeTiS */

  vsizval = 1;                                    /* Assume no vertex communication sizes */
  cvolsum = 0;
  for (vertnum = edgenum = baseval; vertnum < vertnnd; vertnum ++) {
    SCOTCH_Num          partval;
    SCOTCH_Num          edgennd;

    partval = parttax[vertnum];
    nghbtax[partval] = vertnum;                   /* Do not count local neighbors in communication volume */
    if (vsiztax != NULL)
      vsizval = vsiztax[vertnum];

    for (edgennd = verttax[vertnum + 1]; edgenum < edgennd; edgenum ++) { /* Based traversal of edge array adjncy */
      SCOTCH_Num          vertend;                /* Based end vertex number                                      */
      SCOTCH_Num          partend;

      vertend = edgetax[edgenum];
      partend = parttax[vertend];
      if (nghbtax[partend] != vertnum) {          /* If first neighbor in this part */
        nghbtax[partend] = vertnum;               /* Set part as accounted for      */
        cvolsum += vsizval;
      }
    }
  }
  *cvolptr = cvolsum;

  memFree (nghbtax + baseval);                    /* Free un-based array */

  return (METIS_OK);
}

/************************************/
/*                                  */
/* These routines are the C API for */
/* MeTiS graph ordering routines.   */
/*                                  */
/************************************/

/* This routine converts an array of doubles
** into an array of proportionate integers.
** It returns:
** - void  : in all cases.
*/

#define EPSILON 1e-6

void
_SCOTCH_METIS_doubleToInt (
const SCOTCH_Num            valunbr,
const double * const        doubtab,
SCOTCH_Num * const          intetab)
{
  double              doubadj;
  SCOTCH_Num          i;

  for (i = 0, doubadj = 1.0; i < valunbr; i ++) {
    double              doubval;
    double              doubtmp;

    doubval  = doubtab[i];
    doubval *= doubadj;                         /* See if renormalization factor works        */
    doubtmp  = doubval - floor (doubval + EPSILON); /* Determine its possible fractional part */
    if (fabs (doubtmp) >= EPSILON) {            /* If a residual fractional part exists       */
      doubtmp = doubadj / doubtmp;              /* Incorporate it in renormalization factor   */
      doubadj = (doubadj * doubtmp) / (double) intGcd ((int) round (doubadj), (int) round (doubtmp));
    }
  }
  for (i = 0; i < valunbr; i ++)
    intetab[i] = (SCOTCH_Num) round (doubtab[i] * doubadj);
}

/* This routine is the interface between MeTiS
** and Scotch. It computes the partition of a
** weighted or unweighted graph.
** It returns:
** - 0   : if the partition could be computed.
** - !0  : on error.
*/

static
int
_SCOTCH_METIS_PartGraph2 (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    ncon,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
const double * const        tpwgts,
SCOTCH_Num * const          part,
SCOTCH_Num                  flagval,
const double * const        kbalval)
{
  SCOTCH_Graph        grafdat;                    /* Scotch graph object to interface with libScotch */
  SCOTCH_Num *        twintab;                    /* Integer array of target weights (if any)        */
  SCOTCH_Strat        stradat;
  SCOTCH_Num          baseval;
  SCOTCH_Num          vertnbr;
  int                 o;

  twintab = NULL;
  if (tpwgts != NULL) {                           /* If weighted part array         */
    double *            twdbtab;                  /* Double array of target weights */
    SCOTCH_Num          i;

    if ((twintab = malloc (*nparts * sizeof (SCOTCH_Num))) == NULL)
      return (METIS_ERROR_MEMORY);
    if ((twdbtab = malloc (*nparts * sizeof (double))) == NULL) {
      free   (twintab);
      return (METIS_ERROR_MEMORY);
    }
    for (i = 0; i < *nparts; i ++)                /* Gather first target constraint; only this one is used */
      twdbtab[i] = tpwgts[i * *ncon];
    _SCOTCH_METIS_doubleToInt (*nparts, twdbtab, twintab); /* Convert balance array to integers */

    free (twdbtab);                               /* Array no longer of use */
  }
  
  SCOTCH_graphInit (&grafdat);

  baseval = *numflag;
  vertnbr = *n;

  o = 1;                                          /* Assume something will go wrong */
  if (SCOTCH_graphBuild (&grafdat, baseval, vertnbr, xadj, xadj + 1, vwgt, NULL,
                         xadj[vertnbr] - baseval, adjncy, adjwgt) == 0) {
    SCOTCH_stratInit          (&stradat);
    SCOTCH_stratGraphMapBuild (&stradat, flagval, *nparts, *kbalval);
#ifdef SCOTCH_DEBUG_ALL
    if (SCOTCH_graphCheck (&grafdat) == 0)        /* TRICK: next instruction called only if graph is consistent */
#endif /* SCOTCH_DEBUG_ALL */
    {
      if (tpwgts == NULL)
        o = SCOTCH_graphPart (&grafdat, *nparts, &stradat, part);
      else {
        SCOTCH_Arch         archdat;

        if (SCOTCH_archInit (&archdat) == 0) {
          if (SCOTCH_archCmpltw (&archdat, *nparts, twintab) == 0) {
            o = SCOTCH_graphMap (&grafdat, &archdat, &stradat, part);
          }
	  SCOTCH_archExit(&archdat);
        }
      }
    }
    SCOTCH_stratExit (&stradat);
  }
  SCOTCH_graphExit (&grafdat);

  if (twintab != NULL)
    free (twintab);

  if (o != 0)
    return (1);

  if (baseval != 0) {                             /* MeTiS part array is based, Scotch is not */
    SCOTCH_Num          vertnum;

    for (vertnum = 0; vertnum < vertnbr; vertnum ++)
      part[vertnum] += baseval;
  }

  return (0);
}

/*
**
*/

static
int
_SCOTCH_METIS_PartGraph (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    ncon,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
const double * const        tpwgts,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          edgecut,
SCOTCH_Num * const          part,
SCOTCH_Num                  flagval,
const double * const        kbalval)
{
  if (_SCOTCH_METIS_PartGraph2 (n, ncon, xadj, adjncy, vwgt, adjwgt, numflag, nparts, tpwgts, part, flagval, kbalval) != 0) {
    *edgecut = -1;                                /* Indicate error */
    return (METIS_ERROR);
  }

  return (_SCOTCH_METIS_OutputCut (*numflag, *n + *numflag, xadj - *numflag, adjncy - *numflag,
                                   (adjwgt != NULL) ? (adjwgt - *numflag) : NULL, part - *numflag, edgecut));
}

/* Scotch does not directly consider communication volume.
** Instead, vertex communication loads are added to the edge
** loads so as to emulate this behavior: heavily weighted
** edges, connected to heavily communicating vertices, will
** be less likely to be cut.
*/

static
int
_SCOTCH_METIS_PartGraph_Volume (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    ncon,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    vsize,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
const double * const        tpwgts,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          volume,
SCOTCH_Num * const          part,
SCOTCH_Num                  flagval,
const double * const        kbalval)
{
  SCOTCH_Num                  baseval;
  SCOTCH_Num                  vertnbr;
  SCOTCH_Num                  vertnum;
  SCOTCH_Num                  edgenum;
  const SCOTCH_Num * restrict edgetax;
  const SCOTCH_Num * restrict vsiztax;

  baseval = *numflag;
  vertnbr = *n;
  edgetax = adjncy - baseval;

  if (vsize == NULL) {                            /* If no communication load data provided */
    if (_SCOTCH_METIS_PartGraph2 (n, ncon, xadj, adjncy, vwgt, NULL, numflag, nparts, tpwgts, part, flagval, kbalval) != 0)
      return (METIS_ERROR);
    vsiztax = NULL;
  }
  else {                                          /* Will have to turn communication volumes into edge loads */
    SCOTCH_Num                  edgenbr;
    SCOTCH_Num * restrict       edlotax;
    int                         o;

    edgenbr = xadj[vertnbr] - baseval;
    if ((edlotax = memAlloc (edgenbr * sizeof (SCOTCH_Num))) == NULL)
      return (METIS_ERROR);
    edlotax -= baseval;                           /* Base access to edlotax */
    vsiztax  = vsize - baseval;

    for (vertnum = 0, edgenum = baseval;          /* Un-based scan of vertex array xadj */
         vertnum < vertnbr; vertnum ++) {
      SCOTCH_Num          vsizval;                /* Communication size of current vertex */
      SCOTCH_Num          edgennd;

      vsizval = vsize[vertnum];
      for (edgennd = xadj[vertnum + 1]; edgenum < edgennd; edgenum ++) { /* Based traversal of edge array adjncy */
        SCOTCH_Num          vertend;              /* Based end vertex number                                     */

        vertend = edgetax[edgenum];
        edlotax[edgenum] = vsizval + vsiztax[vertend];
      }
    }

    o = _SCOTCH_METIS_PartGraph2 (n, ncon, xadj, adjncy, vwgt, edlotax + baseval, numflag, nparts, tpwgts, part,
                                  flagval, kbalval);

    memFree (edlotax + baseval);

    if (o != 0)
      return (METIS_ERROR);
  }

  return (_SCOTCH_METIS_OutputVol (baseval, *n + baseval, xadj - baseval, adjncy - baseval, vsiztax, *nparts, part - baseval, volume));
}

/*
**
*/

int
SCOTCHMETISNAMES (METIS_V3_PartGraphKway) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    wgtflag,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          edgecut,
SCOTCH_Num * const          part)
{
  const SCOTCH_Num *  vwgt2;
  const SCOTCH_Num *  adjwgt2;
  double              kbalval;

  kbalval = 0.01;
  vwgt2   = ((wgtflag == NULL) || ((*wgtflag & 2) != 0)) ? vwgt   : NULL;
  adjwgt2 = ((wgtflag == NULL) || ((*wgtflag & 1) != 0)) ? adjwgt : NULL;

  return (_SCOTCH_METIS_PartGraph (n, 0, xadj, adjncy, vwgt2, adjwgt2,
                                   numflag, nparts, NULL, options, edgecut, part,
                                   SCOTCH_STRATDEFAULT, &kbalval));
}

/*
**
*/

int
SCOTCHMETISNAMES (METIS_V3_PartGraphRecursive) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    wgtflag,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          edgecut,
SCOTCH_Num * const          part)
{
  const SCOTCH_Num *  vwgt2;
  const SCOTCH_Num *  adjwgt2;
  double              kbalval;

  kbalval = 0.01;
  vwgt2   = ((wgtflag == NULL) || ((*wgtflag & 2) != 0)) ? vwgt   : NULL;
  adjwgt2 = ((wgtflag == NULL) || ((*wgtflag & 1) != 0)) ? adjwgt : NULL;

  return (_SCOTCH_METIS_PartGraph (n, 0, xadj, adjncy, vwgt2, adjwgt2,
                                   numflag, nparts, NULL, options, edgecut, part,
                                   SCOTCH_STRATRECURSIVE, &kbalval));
}

/*
**
*/

int
SCOTCHMETISNAMES (METIS_V3_PartGraphVKway) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    vsize,
const SCOTCH_Num * const    wgtflag,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          volume,
SCOTCH_Num * const          part)
{
  const SCOTCH_Num *  vwgt2;
  const SCOTCH_Num *  vsize2;
  double              kbalval;

  kbalval = 0.01;
  vsize2  = ((wgtflag == NULL) || ((*wgtflag & 1) != 0)) ? vsize : NULL;
  vwgt2   = ((wgtflag == NULL) || ((*wgtflag & 2) != 0)) ? vwgt  : NULL;

  return (_SCOTCH_METIS_PartGraph_Volume (n, 0, xadj, adjncy, vwgt2, vsize2,
                                          numflag, nparts, NULL, options, volume, part,
                                          SCOTCH_STRATDEFAULT, &kbalval));
}

/*
**
*/

int
SCOTCHMETISNAMES (METIS_V5_PartGraphKway) (
const SCOTCH_Num * const    nvtxs,
const SCOTCH_Num * const    ncon,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    vsize,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    nparts,
const double * const        tpwgts,
const double * const        ubvec,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          objval,
SCOTCH_Num * const          part)
{
  SCOTCH_Num          baseval;
  double              kbalval;

  kbalval = 0.01;
  baseval = ((options != NULL) && (options != xadj)) ? options[METIS_OPTION_NUMBERING] : 0;

  return ((vsize == NULL)
          ? _SCOTCH_METIS_PartGraph (nvtxs, ncon, xadj, adjncy, vwgt, adjwgt,
                                     &baseval, nparts, tpwgts, options, objval, part,
                                     SCOTCH_STRATDEFAULT, (ubvec != NULL) ? ubvec : &kbalval)
          : _SCOTCH_METIS_PartGraph_Volume (nvtxs, ncon, xadj, adjncy, vwgt, vsize,
                                            &baseval, nparts, tpwgts, options, objval, part,
                                            SCOTCH_STRATDEFAULT, (ubvec != NULL) ? ubvec : &kbalval));
}

/*
**
*/

int
SCOTCHMETISNAMES (METIS_V5_PartGraphRecursive) (
const SCOTCH_Num * const    nvtxs,
const SCOTCH_Num * const    ncon,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    vsize,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    nparts,
const double * const        tpwgts,
const double * const        ubvec,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          objval,
SCOTCH_Num * const          part)
{
  SCOTCH_Num          baseval;

  baseval = ((options != NULL) && (options != xadj)) ? options[METIS_OPTION_NUMBERING] : 0;

  return ((vsize == NULL)
          ? _SCOTCH_METIS_PartGraph (nvtxs, ncon, xadj, adjncy, vwgt, adjwgt,
                                     &baseval, nparts, tpwgts, options, objval, part,
                                     SCOTCH_STRATRECURSIVE, ubvec)
          : _SCOTCH_METIS_PartGraph_Volume (nvtxs, ncon, xadj, adjncy, vwgt, vsize,
                                            &baseval, nparts, tpwgts, options, objval, part,
                                            SCOTCH_STRATRECURSIVE, ubvec));
}

/*******************/
/*                 */
/* MeTiS v3 stubs. */
/*                 */
/*******************/

#if (SCOTCH_METIS_VERSION == 3)

int
SCOTCHMETISNAMEC (METIS_PartGraphKway) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    wgtflag,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          edgecut,
SCOTCH_Num * const          part)
{
  return (SCOTCHMETISNAMES (METIS_V3_PartGraphKway) (n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag,
                                                     nparts, options, edgecut, part));
}

/*
**
*/

int
SCOTCHMETISNAMEC (METIS_PartGraphRecursive) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    wgtflag,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          edgecut,
SCOTCH_Num * const          part)
{
  return (SCOTCHMETISNAMES (METIS_V3_PartGraphRecursive) (n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag,
                                                          nparts, options, edgecut, part));
}

/*
**
*/

int
SCOTCHMETISNAMEC (METIS_PartGraphVKway) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    vsize,
const SCOTCH_Num * const    wgtflag,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    nparts,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          volume,
SCOTCH_Num * const          part)
{
  return (SCOTCHMETISNAMES (METIS_V3_PartGraphVKway) (n, xadj, adjncy, vwgt, vsize, wgtflag, numflag,
                                                      nparts, options, volume, part));
}

#endif /* (SCOTCH_METIS_VERSION == 3) */

/*******************/
/*                 */
/* MeTiS v5 stubs. */
/*                 */
/*******************/

#if (SCOTCH_METIS_VERSION == 5)

int
SCOTCHMETISNAMEC (METIS_PartGraphKway) (
const SCOTCH_Num * const    nvtxs,
const SCOTCH_Num * const    ncon,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    vsize,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    nparts,
const double * const        tpwgts,
const double * const        ubvec,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          objval,
SCOTCH_Num * const          part)
{
  return (SCOTCHMETISNAMES (METIS_V5_PartGraphKway) (nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt,
                                                     nparts, tpwgts, ubvec, options, objval, part));
}

/*
**
*/

int
SCOTCHMETISNAMEC (METIS_PartGraphRecursive) (
const SCOTCH_Num * const    nvtxs,
const SCOTCH_Num * const    ncon,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    vsize,
const SCOTCH_Num * const    adjwgt,
const SCOTCH_Num * const    nparts,
const double * const        tpwgts,
const double * const        ubvec,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          objval,
SCOTCH_Num * const          part)
{
  return (SCOTCHMETISNAMES (METIS_V5_PartGraphRecursive) (nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt,
                                                          nparts, tpwgts, ubvec, options, objval, part));
}

#endif /* (SCOTCH_METIS_VERSION == 5) */
