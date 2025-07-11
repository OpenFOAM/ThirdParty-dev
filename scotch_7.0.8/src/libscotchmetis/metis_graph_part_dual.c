/* Copyright 2020,2021,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : metis_graph_part_dual.c                 **/
/**                                                        **/
/**   AUTHOR     : Marc FUENTES                            **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the compatibility        **/
/**                library for the MeTiS partitioning      **/
/**                routines containing routines relative   **/
/**                to dual graphs                          **/
/**                                                        **/
/**   DATES      : # Version 6.1  : from : 01 sep 2020     **/
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
#include "graph.h"
#include "mesh.h"
#include "scotch.h"
#include "metis.h"                                /* Our "metis.h" file */
#include "metis_graph_dual.h"
#include "metis_graph_part.h"

#define EPSILON                     1.0e-6        /* Rounding percision */

/* This routine computes the partition of the
** dual graph of a mesh.
** It returns:
** - METIS_OK      : if the partition has been successfully computed.
** - METIS_ERROR*  : on error.
*/

int
SCOTCHMETISNAMES (METIS_PartMeshDual) (
const SCOTCH_Num * const    ne,                   /*+ Pointer to number of elements                  +*/
const SCOTCH_Num * const    nn,                   /*+ Pointer to number of nodes                     +*/
const SCOTCH_Num * const    eptr,                 /*+ Element vertex array                           +*/
const SCOTCH_Num * const    eind,                 /*+ Element edge array                             +*/
const SCOTCH_Num * const    vwgt,                 /*+ Element vertex weight array                    +*/
const SCOTCH_Num * const    vsize,                /*+ Element communication volume array             +*/
const SCOTCH_Num * const    ncommon,              /*+ Pointer to number of connectivity common nodes +*/
const SCOTCH_Num * const    nparts,               /*+ Pointer to number of parts to compute          +*/
const double * const        tpwgts,               /*+ Pointer to part weight array                   +*/
const SCOTCH_Num * const    options,              /*+ Optional options array                         +*/
SCOTCH_Num * const          objval,               /*+ Edege cut or communication volume return value +*/
SCOTCH_Num * const          epart,                /*+ Element partition array to be filled           +*/
SCOTCH_Num * const          npart)                /*+ Node partition array to be filled              +*/
{
  SCOTCH_Graph        grafdat;
  SCOTCH_Graph        graftmp;
  SCOTCH_Arch         archdat;
  SCOTCH_Strat        stradat;
  SCOTCH_Num          baseval;                    /* Global base value                   */
  SCOTCH_Num          velmbas;                    /* Base value of mesh element vertices */
  SCOTCH_Num          vnodbas;                    /* Base value of mesh node vertices    */
  SCOTCH_Num          vnodnum;
  SCOTCH_Num          vertnbr;
  SCOTCH_Num *        verttab;
  SCOTCH_Num *        vendtab;
  SCOTCH_Num          edgenbr;
  SCOTCH_Num *        edgetab;
  SCOTCH_Num *        edlotab;
  SCOTCH_Mesh         meshdat;
  int                 o;

  baseval = 0;                                    /* Assume base value is 0 */

  if ((options != NULL) && (options != ne))       /* For the Fortran interface, null arrays are those equal to ne */
    baseval = options[METIS_OPTION_NUMBERING];    /* Get base value                                               */

  if ((tpwgts != NULL) && ((const SCOTCH_Num * const) tpwgts != ne)) {
    SCOTCH_Num * restrict wtgttab;                /* Array of integer loads            */
    double                wtgtsum;                /* Sum of floating-point part weigts */
    SCOTCH_Num            partnum;

    for (partnum = 0, wtgtsum = 0.0; partnum < *nparts; partnum ++)
      wtgtsum += tpwgts[partnum];                 /* Sum floating-point part weights */
    if (fabs (wtgtsum - 1.0) >= EPSILON) {
      SCOTCH_errorPrint ("METIS_PartMeshDual: invalid partition weight sum");
      *objval = METIS_ERROR_INPUT;
      return (METIS_ERROR_INPUT);
    }
    if ((wtgttab = memAlloc (*nparts * sizeof (SCOTCH_Num))) == NULL) {
      SCOTCH_errorPrint ("METIS_PartMeshDual: out of memory (1)");
      *objval = METIS_ERROR_MEMORY;
      return (METIS_ERROR_MEMORY);
    }

    _SCOTCH_METIS_doubleToInt (*nparts, tpwgts, wtgttab); /* Convert array of doubles to array of ints */

    SCOTCH_archInit (&archdat);
    o = SCOTCH_archCmpltw (&archdat, *nparts, wtgttab);
    memFree (wtgttab);

    if (o != 0) {
      SCOTCH_errorPrint ("METIS_PartMeshDual: cannot create weighted architecture");
      SCOTCH_archExit   (&archdat);
      *objval = METIS_ERROR_MEMORY;
      return (METIS_ERROR_MEMORY);
    }
  }
  else {                                          /* Else create unweighted target architecture */
    SCOTCH_archInit  (&archdat);
    SCOTCH_archCmplt (&archdat, *nparts);
  }

  SCOTCH_meshInit (&meshdat);
  if ((o = _SCOTCH_METIS_MeshToDual2 (&meshdat, baseval, *nn, *ne, eptr, eind)) != METIS_OK) {
    SCOTCH_errorPrint ("METIS_PartMeshDual: cannot build dual mesh");
    SCOTCH_archExit   (&archdat);
    *objval = o;                                  /* Error value for the Fortran interface */
    return (o);
  }

  SCOTCH_graphInit (&grafdat);
  o = SCOTCH_meshGraphDual (&meshdat, &grafdat, *ncommon);
  if (o != 0) {
    SCOTCH_errorPrint ("METIS_PartMeshDual: cannot build dual graph");
    SCOTCH_meshExit   (&meshdat);
    SCOTCH_graphExit  (&grafdat);
    SCOTCH_archExit   (&archdat);
    *objval = METIS_ERROR_MEMORY;
    return (METIS_ERROR_MEMORY);
  }

  SCOTCH_graphData (&grafdat, NULL, &vertnbr, &verttab, &vendtab, NULL, NULL, &edgenbr, &edgetab, NULL); /* Get graph topology arrays */

  if ((vsize != NULL) && (vsize != ne)) {
    const SCOTCH_Num * restrict vsiztax;
    Gnum                        vertnum;
    Gnum                        edgenum;
    SCOTCH_Num * restrict       edgetax;
    SCOTCH_Num * restrict       edlotax;

    if ((edlotab = memAlloc (edgenbr * sizeof (SCOTCH_Num))) == NULL) {
      SCOTCH_errorPrint ("METIS_PartMeshDual: out of memory (2)");
      SCOTCH_meshExit   (&meshdat);
      SCOTCH_graphExit  (&grafdat);
      SCOTCH_archExit   (&archdat);
      *objval = METIS_ERROR_MEMORY;
      return (METIS_ERROR_MEMORY);
    }
    edlotax = edlotab - baseval;                  /* Base access to edlotax */
    edgetax = edgetab - baseval;
    vsiztax = vsize   - baseval;

    for (vertnum = 0, edgenum = baseval;          /* Un-based scan of vertex array verttab */
         vertnum < vertnbr; vertnum ++) {
      SCOTCH_Num          vsizval;                /* Communication size of current vertex */
      SCOTCH_Num          edgennd;

      vsizval = vsize[vertnum];
      for (edgennd = vendtab[vertnum]; edgenum < edgennd; edgenum ++) { /* Based traversal of compact edge array */
        SCOTCH_Num          vertend;              /* Based end vertex number                                     */

        vertend = edgetax[edgenum];
        edlotax[edgenum] = vsizval + vsiztax[vertend];
      }
    }
  }
  else                                            /* No edge array */
    edlotab = NULL;

  SCOTCH_graphInit  (&graftmp);
  SCOTCH_graphBuild (&graftmp, baseval, vertnbr, verttab, vendtab, vwgt, NULL, edgenbr, edgetab, edlotab);
  SCOTCH_stratInit  (&stradat);

  if (SCOTCH_graphMap (&graftmp, &archdat, &stradat, epart) == 0) {
    if (baseval != 0) {                           /* MeTiS part array is based, Scotch is not */
      SCOTCH_Num          vertnum;

      for (vertnum = 0; vertnum < vertnbr; vertnum ++)
        epart[vertnum] += baseval;
    }

    if (((vsize != NULL) && (vsize != ne)) ||     /* If computation of communication volume wanted */
        ((options != NULL) && (options != ne) && (options[METIS_OPTION_OBJTYPE] == METIS_OBJTYPE_VOL)))
      o = _SCOTCH_METIS_OutputVol (baseval, vertnbr + baseval, verttab - baseval, edgetab - baseval, vsize - baseval, *nparts, epart - baseval, objval);
    else
      o = _SCOTCH_METIS_OutputCut (baseval, vertnbr + baseval, verttab - baseval, edgetab - baseval,
                                   (edlotab != NULL) ? (edlotab - baseval) : NULL, epart - baseval, objval);
  }
  else
    o = METIS_ERROR;

  SCOTCH_stratExit (&stradat);
  SCOTCH_graphExit (&graftmp);
  if (edlotab != NULL)
    memFree (edlotab);
  SCOTCH_archExit  (&archdat);
  SCOTCH_graphExit (&grafdat);

  if (o != METIS_OK) {
    SCOTCH_errorPrint ("METIS_PartMeshDual: could not partition graph");
    SCOTCH_meshExit   (&meshdat);
    *objval = METIS_ERROR;
    return (METIS_ERROR);
  }

  SCOTCH_meshData (&meshdat, &velmbas, &vnodbas, NULL, NULL, &verttab, &vendtab, NULL, NULL, NULL, NULL, &edgetab, NULL);
  for (vnodnum = 0; vnodnum < *nn; vnodnum ++) {  /* Un-based loop on node vertices */
    SCOTCH_Num          edgenum;

    edgenum = verttab[vnodnum + (vnodbas - baseval)];
    npart[vnodnum] = (edgenum != vendtab[vnodnum + (vnodbas - baseval)]) ? epart[edgetab[edgenum - baseval] - velmbas] : 0; /* If non-isolated node, get part of first neighbor element, else 0 */
  }

  SCOTCH_meshExit (&meshdat);

  return (METIS_OK);
}

/*
**
*/

int
SCOTCHMETISNAMEC (METIS_PartMeshDual) (
const SCOTCH_Num * const    ne,                   /*+ Pointer to number of elements                  +*/
const SCOTCH_Num * const    nn,                   /*+ Pointer to number of nodes                     +*/
const SCOTCH_Num * const    eptr,                 /*+ Element vertex array                           +*/
const SCOTCH_Num * const    eind,                 /*+ Element edge array                             +*/
const SCOTCH_Num * const    vwgt,                 /*+ Element vertex weight array                    +*/
const SCOTCH_Num * const    vsize,                /*+ Element communication volume array             +*/
const SCOTCH_Num * const    ncommon,              /*+ Pointer to number of connectivity common nodes +*/
const SCOTCH_Num * const    nparts,               /*+ Pointer to number of parts to compute          +*/
const double * const        tpwgts,               /*+ Pointer to part weight array                   +*/
const SCOTCH_Num * const    options,              /*+ Optional options array                         +*/
SCOTCH_Num * const          objval,               /*+ Edege cut or communication volume return value +*/
SCOTCH_Num * const          epart,                /*+ Element partition array to be filled           +*/
SCOTCH_Num * const          npart)                /*+ Node partition array to be filled              +*/
{
  return (SCOTCHMETISNAMES (METIS_PartMeshDual) (ne, nn, eptr, eind, vwgt, vsize, ncommon,
                                                 nparts, tpwgts, options, objval, epart, npart));
}
