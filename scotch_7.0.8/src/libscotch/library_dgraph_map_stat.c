/* Copyright 2008-2010,2012,2018,2019,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_map_stat.c               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Clement BARTHELEMY                      **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the distri-  **/
/**                buted mapping routines of the libSCOTCH **/
/**                library.                                **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 10 sep 2024     **/
/**                                 to   : 19 sep 2024     **/
/**                                                        **/
/**   NOTES      : # This code is directly derived from    **/
/**                  the code formerly present in routine  **/
/**                  SCOTCH_dgraphMapView().               **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "context.h"
#include "parser.h"
#include "dgraph.h"
#include "dgraph_halo.h"
#include "arch.h"
#include "dmapping.h"
#include "kdgraph.h"
#include "library_dmapping.h"
#include "library_dgraph_map_stat.h"
#include "ptscotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the mapping routines.            */
/*                                  */
/************************************/

/*+ This routine computes distributed mapping
*** statistics and returns them in the pointer
*** arguments.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphMapStat (
SCOTCH_Dgraph * const         libgrafptr,
const SCOTCH_Dmapping * const libmappptr,
SCOTCH_Num * const            tgtnbrptr,          /* Number of vertices in target architecture                      */
SCOTCH_Num * const            mapnbrptr,          /* Number of target terminal domains used in mapping              */
SCOTCH_Num * const            mapminptr,          /* Minimum load in domain (with respect to domain weight)         */
SCOTCH_Num * const            mapmaxptr,          /* Maximum load in domain (with respect to domain weight)         */
double * const                mapavgptr,          /* Average load in domain (with respect to domain weights)        */
double * const                mapdltptr,          /* Standard deviation of mapping loads                            */
SCOTCH_Num * const            ngbsumptr,          /* Total number of neighbors                                      */
SCOTCH_Num * const            ngbminptr,          /* Minimum number of neighbors                                    */
SCOTCH_Num * const            ngbmaxptr,          /* Maximum number of neighbors                                    */
SCOTCH_Num * const            cdstmaxptr,         /* Maximum distance between domains                               */
SCOTCH_Num                    cdsttab[256],       /* Communication distance histogram array for mapped vertices     */
SCOTCH_Num * const            cmlosumptr,         /* Communication load (sum of edge weights) for mapped vertices   */
SCOTCH_Num * const            cmdisumptr,         /* Communication dilation (sum of distances) for mapped vertices  */
SCOTCH_Num * const            cmexsumptr)         /* Communication expansion (sum of edge weights * distances) etc. */
{
  ArchDom                 dorgdat;                /* Largest domain in architecture          */
  unsigned int * restrict nmskloctab;             /* Local neighbor bitfield                 */
  unsigned int * restrict nmskglbtab;             /* Local neighbor bitfield                 */
  int                     nmskidxnbr;             /* Size of bitfield; int since sent by MPI */
  Gnum * restrict         tgloloctab;             /* Local array of terminal domain loads    */
  Gnum * restrict         tgloglbtab;             /* Global array of terminal domain loads   */
  Gnum * restrict         termgsttax;             /* Terminal domain ghost mapping array     */
  Anum                    tgtnbr;
  Anum                    tgtnum;
  Anum                    mapnbr;
  Gnum                    mapmin;
  Gnum                    mapmax;
  Gnum                    mapsum;
  double                  mapavg;
  double                  mapdlt;
  Gnum                    vertlocnum;
  Gnum                    veloval;
  int                     chekloctab[3];          /* Array of local permissions  */
  int                     chekglbtab[3];          /* Array of global permissions */
  DgraphHaloRequest       requdat;

  Dgraph * const              grafptr = (Dgraph *) CONTEXTOBJECT (libgrafptr);
  const LibDmapping * const   mappptr = (LibDmapping *) libmappptr;
  const Arch * restrict const archptr = &mappptr->m.archdat;
  const Gnum * restrict const veloloctax = grafptr->veloloctax;

  archDomFrst (archptr, &dorgdat);                /* Get architecture domain  */
  tgtnbr = archDomSize (archptr, &dorgdat);       /* Get architecture size    */
  if (tgtnbrptr != NULL)                          /* Set size of architecture */
    *tgtnbrptr = tgtnbr;

  if ((grafptr->vertglbnbr == 0) ||               /* If nothing to do */
      (grafptr->edgeglbnbr == 0)) {
    if (mapnbrptr != NULL)                        /* Fill needed values with empty results */
      *mapnbrptr = 0;
    if (mapminptr != NULL)
      *mapminptr = 0;
    if (mapmaxptr != NULL)
      *mapmaxptr = 0;
    if (mapavgptr != NULL)
      *mapavgptr = 0;
    if (mapdltptr != NULL)
      *mapdltptr = 0;
    if (ngbminptr != NULL)
      *ngbminptr = 0;
    if (ngbmaxptr != NULL)
      *ngbmaxptr = 0;
    if (ngbsumptr != NULL)
      *ngbsumptr = 0;
    if (cdstmaxptr != NULL)
      *cdstmaxptr = 0;
    if (cdsttab != NULL)
      memSet (cdsttab, 0, 256 * sizeof (SCOTCH_Num));
    if (cmlosumptr != NULL)
      *cmlosumptr = 0;
    if (cmdisumptr != NULL)
      *cmdisumptr = 0;
    if (cmexsumptr != NULL)
      *cmexsumptr = 0;

    return (0);
  }

  if (archVar (archptr)) {
    errorPrint (STRINGIFY (SCOTCH_dgraphMapStat) ": not implemented");
    return (1);
  }

  chekloctab[1] =
  chekloctab[2] = LIBDGRAPHMAPSTATNONE;           /* Assume nothing to do */
  if ((cdstmaxptr != NULL) ||
      (cdsttab    != NULL) ||
      (cmlosumptr != NULL) ||
      (cmdisumptr != NULL) ||
      (cmexsumptr != NULL))
    chekloctab[1] = LIBDGRAPHMAPSTATCOMM;
  if ((ngbminptr != NULL) ||
      (ngbmaxptr != NULL) ||
      (ngbsumptr != NULL))
    chekloctab[2] = LIBDGRAPHMAPSTATNGHB;

  if (dgraphGhst (grafptr) != 0) {                /* Compute ghost edge array if not already present */
    errorPrint (STRINGIFY (SCOTCH_dgraphMapStat) ": cannot compute ghost edge array");
    return (1);
  }

  nmskidxnbr    = (tgtnbr + 1 + ((sizeof (int) << 3) - 1)) / (sizeof (int) << 3); /* Size of neighbor subdomain bitfield; TRICK: "+1" to have a "-1" cell for unmapped vertices */
  chekloctab[0] = 0;
  if (memAllocGroup ((void **) (void *)
                     &nmskloctab, (size_t) (nmskidxnbr          * sizeof (unsigned int)),
                     &nmskglbtab, (size_t) (nmskidxnbr          * sizeof (unsigned int)),
                     &tgloloctab, (size_t) ((tgtnbr + 1)        * sizeof (Gnum)), /* TRICK: "+1" to have a "-1" cell for unmapped vertices */
                     &tgloglbtab, (size_t) ((tgtnbr + 1)        * sizeof (Gnum)),
                     &termgsttax, (size_t) (grafptr->vertgstnbr * sizeof (Gnum)), NULL) == NULL) /* TRICK: array not yet based */
    chekloctab[0] = 1;

  if (MPI_Allreduce (chekloctab, chekglbtab, 3, MPI_INT, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint (STRINGIFY (SCOTCH_dgraphMapStat) ": communication error (1)");
    return (1);
  }
  if (chekglbtab[0] != 0) {                       /* If memory allocation error  */
    if (nmskloctab != NULL)                       /* Free group leader if needed */
      memFree (nmskloctab);
    errorPrint (STRINGIFY (SCOTCH_dgraphMapStat) ": out of memory");
    return (1);
  }

  if (dmapTerm (&mappptr->m, grafptr, termgsttax) != 0) { /* Fill array of terminal domain numbers for all local vertices */
    errorPrint (STRINGIFY (SCOTCH_dgraphMapStat) ": cannot build local terminal array");
    memFree    (nmskloctab);                      /* Free group leader */
    return (1);
  }

  dgraphHaloAsync (grafptr, termgsttax, GNUM_MPI, &requdat); /* Exchange array of terminal domain numbers across processes */ 
  termgsttax -= grafptr->baseval;

  memSet (tgloloctab, 0, (tgtnbr + 1) * sizeof (Gnum));
  tgloloctab ++;                                  /* TRICK: trim arrays for "-1" cell */
  tgloglbtab ++;

  veloval = 1;
  for (vertlocnum = grafptr->baseval; vertlocnum < grafptr->vertlocnnd; vertlocnum ++) {
#ifdef SCOTCH_DEBUG_DMAP2
    if ((termgsttax[vertlocnum] < -1) || (termgsttax[vertlocnum] >= tgtnbr)) {
      errorPrint (STRINGIFY (SCOTCH_dgraphMapStat) ": invalid local terminal array");
      tgloloctab[-1] = - (grafptr->veloglbsum + 1); /* Record an error with a value that cannot be compensated for */
      break;
    }
#endif /* SCOTCH_DEBUG_DMAP2 */
    if (veloloctax != NULL)
      veloval = veloloctax[vertlocnum];
    tgloloctab[termgsttax[vertlocnum]] += veloval; /* One more vertex of given weight assigned to this target */
  }

  if (MPI_Allreduce (tgloloctab - 1, tgloglbtab - 1, tgtnbr + 1, GNUM_MPI, MPI_SUM, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint (STRINGIFY (SCOTCH_dgraphMapStat) ": communication error (2)");
    tgloglbtab[-1] = -1;                          /* Set an error similar to that of an invalid array value */
  }
  if (tgloglbtab[-1] < 0) {                       /* If error detected                            */
    dgraphHaloWait (&requdat);                    /* Wait for ghost terminal data to be exchanged */
    memFree        (nmskloctab);                  /* Free group leader                            */
    return (1);
  }

  mapmin = GNUMMAX;
  mapmax = 0;
  mapsum = 0;
  mapnbr = 0;
  for (tgtnum = 0; tgtnum < tgtnbr; tgtnum ++) {  /* For all valid terminal domain numbers (not -1) */
    Gnum                tgloval;

    tgloval = tgloglbtab[tgtnum];
    if (tgloval != 0) {                           /* If target domain has load         */
      ArchDom             domndat;                /* Current domain in architecture    */
      Anum                dwgtval;                /* Weight of terminal domain         */
      Gnum                tgloavg;                /* Average load per target unit load */

      mapnbr ++;                                  /* One more terminal domain used */
      mapsum += tgloval;                          /* Aggregate load mapped onto it */

      archDomTerm (archptr, &domndat, tgtnum);    /* Get weight of terminal domain */
      dwgtval = archDomWght (archptr, &domndat);
      tgloavg = (Gnum) (((double) tgloval / (double) dwgtval) + 0.5L); /* Compute average per domain load */
      if (mapmin > tgloavg)
        mapmin = tgloavg;
      if (mapmax < tgloavg)
        mapmax = tgloavg;
    }
  }
#ifdef SCOTCH_DEBUG_DMAP2
  if ((mapsum + tgloglbtab[-1]) != grafptr->veloglbsum) {
    errorPrint     (STRINGIFY (SCOTCH_dgraphMapStat) ": invalid mapping data");
    dgraphHaloWait (&requdat);                    /* Wait for ghost terminal data to be exchanged */
    memFree        (nmskloctab);                  /* Free group leader                            */
    return (1);
  }
#endif /* SCOTCH_DEBUG_DMAP2 */
  mapavg = (mapnbr == 0) ? 0.0L : ((double) mapsum / (double) mapnbr);
  mapdlt = 0.0L;
  for (tgtnum = 0; tgtnum < tgtnbr; tgtnum ++) {
    if (tgloglbtab[tgtnum] > 0)                   /* If non-empty domain */
      mapdlt += fabs ((double) tgloglbtab[tgtnum] - mapavg);
  }
  mapdlt = (mapnbr != 0) ? mapdlt / ((double) mapnbr * mapavg) : 0.0L;

  if (mapnbrptr != NULL)
    *mapnbrptr = (SCOTCH_Num) mapnbr;
  if (mapminptr != NULL)
    *mapminptr = (SCOTCH_Num) mapmin;
  if (mapmaxptr != NULL)
    *mapmaxptr = (SCOTCH_Num) mapmax;
  if (mapavgptr != NULL)
    *mapavgptr = mapavg;
  if (mapdltptr != NULL)
    *mapdltptr = mapdlt;

  if (dgraphHaloWait (&requdat) != 0) {           /* Wait for ghost terminal data to be exchanged */
    errorPrint (STRINGIFY (SCOTCH_dgraphMapStat) ": cannot complete asynchronous halo exchange");
    memFree    (nmskloctab);                      /* Free group leader */
    return (1);
  }

  if (chekglbtab[1] == LIBDGRAPHMAPSTATCOMM) {    /* If communication data wanted    */
    Gnum                cdstloctab[256 + 3];      /* Local distance histogram array  */
    Gnum                cdstglbtab[256 + 3];      /* Global distance histogram array */
    Gnum                cmlolocsum;               /* Local sum of communication load */
    Gnum                cmdilocsum;               /* Local sum of dilation           */
    Gnum                cmexlocsum;               /* Local sum of expansion          */
    Gnum                vertlocnum;
    Gnum                edloval;

    const Gnum * restrict const vertloctax = grafptr->vertloctax; /* Fast accesses */
    const Gnum * restrict const vendloctax = grafptr->vendloctax;
    const Gnum * restrict const edgegsttax = grafptr->edgegsttax;
    const Gnum * restrict const edloloctax = grafptr->edloloctax;

    memSet (cdstloctab, 0, 256 * sizeof (Gnum));  /* Initialize the data */
    cmlolocsum  =
    cmdilocsum =
    cmexlocsum = 0;

    edloval = 1;
    for (vertlocnum = grafptr->baseval; vertlocnum < grafptr->vertlocnnd; vertlocnum ++) { /* For all local vertices */
      Gnum                termlocnum;
      ArchDom             termdomdat;
      Gnum                edgelocnum;
      Gnum                edgelocnnd;

      termlocnum = termgsttax[vertlocnum];
      if (termlocnum == ~0)                       /* Skip unmapped vertices */
        continue;

      archDomTerm (archptr, &termdomdat, termlocnum);

      for (edgelocnum = vertloctax[vertlocnum], edgelocnnd = vendloctax[vertlocnum];
           edgelocnum < edgelocnnd; edgelocnum ++) {
        ArchDom             termdomend;
        Gnum                termgstend;
        Anum                distval;

        termgstend = termgsttax[edgegsttax[edgelocnum]];
        if (termgstend == ~0)                     /* Skip unmapped end vertices */
          continue;

        distval = 0;
        if (edloloctax != NULL)                   /* Get edge weight if any */
          edloval = edloloctax[edgelocnum];
        if (termgstend != termlocnum) {           /* If not same domain, compute distance */
          archDomTerm (archptr, &termdomend, termgstend);
          distval = archDomDist (archptr, &termdomdat, &termdomend);
        }
        cdstloctab[(distval > 255) ? 255 : distval] += edloval;
        cmlolocsum += edloval;
        cmdilocsum += distval;
        cmexlocsum += distval * edloval;
      }
    }
    cdstloctab[256 + 0] = cmlolocsum;             /* Write local sums for collective summing */
    cdstloctab[256 + 1] = cmdilocsum;
    cdstloctab[256 + 2] = cmexlocsum;

    if (MPI_Allreduce (cdstloctab, cdstglbtab, 256 + 3, GNUM_MPI, MPI_SUM, grafptr->proccomm) != MPI_SUCCESS) {
      errorPrint (STRINGIFY (SCOTCH_dgraphMapStat) ": communication error (3)");
      memFree    (nmskloctab);                    /* Free group leader */
      return (1);
    }

    if (cdsttab != NULL) {
      Anum                cdstnum;

      for (cdstnum = 0; cdstnum < 256; cdstnum ++) /* Arcs have been accounted for twice */
        cdsttab[cdstnum] = cdstglbtab[cdstnum] / 2;
    }
    if (cdstmaxptr != NULL) {
      Anum                cdstmax;

      for (cdstmax = 255; cdstmax != -1; cdstmax --) /* Find longest distance */
        if (cdstglbtab[cdstmax] != 0)
          break;

      *cdstmaxptr = (SCOTCH_Num) cdstmax;
    }
    if (cmlosumptr != NULL)
      *cmlosumptr = (SCOTCH_Num) cdstglbtab[256 + 0] / 2; /* Arcs have been accounted for twice */
    if (cmdisumptr != NULL)
      *cmdisumptr = (SCOTCH_Num) cdstglbtab[256 + 1] / 2;
    if (cmexsumptr != NULL)
      *cmexsumptr = (SCOTCH_Num) cdstglbtab[256 + 2] / 2;
  }

  if (chekglbtab[2] == LIBDGRAPHMAPSTATNGHB) {    /* If neighbor data wanted */
    Anum                ngbsum;
    Anum                ngbmin;
    Anum                ngbmax;
    Anum                tgtnum;

    const Gnum * restrict const vertloctax = grafptr->vertloctax; /* Fast accesses */
    const Gnum * restrict const vendloctax = grafptr->vendloctax;
    const Gnum * restrict const edgegsttax = grafptr->edgegsttax;

    ngbmin = ANUMMAX;
    ngbmax = 0;
    ngbsum = 0;
    for (tgtnum = 0; tgtnum < tgtnbr; tgtnum ++) { /* For all subdomain indices */
      int                 nmskidxnum;
      Gnum                vertlocnum;
      Anum                ngbnbr;

      if (tgloglbtab[tgtnum] <= 0)                /* If empty subdomain, skip it */
        continue;

      memSet (nmskloctab, 0, nmskidxnbr * sizeof (int)); /* Reset neighbor bit mask */

      for (vertlocnum = grafptr->baseval; vertlocnum < grafptr->vertlocnnd; vertlocnum ++) { /* For all local vertices */
        Gnum                termnum;
        Gnum                edgelocnum;
        Gnum                edgelocnnd;

        termnum = termgsttax[vertlocnum];
        if (termnum != tgtnum)                    /* If vertex does not belong to current part or is not mapped, skip it */
          continue;

        for (edgelocnum = vertloctax[vertlocnum], edgelocnnd = vendloctax[vertlocnum];
             edgelocnum < edgelocnnd; edgelocnum ++) {
          Gnum                termend;

          termend = termgsttax[edgegsttax[edgelocnum]];
          if (termend != tgtnum) {                /* If edge is not internal             */
            termend ++;                           /* TRICK: turn unmapped to 0 and so on */
            nmskloctab[termend / (sizeof (int) << 3)] |= 1 << (termend & ((sizeof (int) << 3) - 1)); /* Flag neighbor in bit array */
          }
        }
      }
      nmskloctab[0] &= ~1;                        /* Do not account for unmapped vertices (terminal domain 0 because of "+1") */

      if (MPI_Allreduce (nmskloctab, nmskglbtab, nmskidxnbr, MPI_INT, MPI_BOR, grafptr->proccomm) != MPI_SUCCESS) {
        errorPrint (STRINGIFY (SCOTCH_dgraphMapStat) ": communication error (4)");
        memFree    (nmskloctab);                  /* Free group leader */
        return (1);
      }

      for (nmskidxnum = 0, ngbnbr = 0; nmskidxnum < nmskidxnbr; nmskidxnum ++) {
        unsigned int        nmskbitval;

        for (nmskbitval = nmskglbtab[nmskidxnum]; nmskbitval != 0; nmskbitval >>= 1)
          ngbnbr += nmskbitval & 1;
      }

      ngbsum += ngbnbr;
      if (ngbmin > ngbnbr)
        ngbmin = ngbnbr;
      if (ngbmax < ngbnbr)
        ngbmax = ngbnbr;
    }

    if (ngbminptr != NULL)
      *ngbminptr = (SCOTCH_Num) ngbmin;
    if (ngbmaxptr != NULL)
      *ngbmaxptr = (SCOTCH_Num) ngbmax;
    if (ngbsumptr != NULL)
      *ngbsumptr = (SCOTCH_Num) ngbsum;
  }

  memFree (nmskloctab);                           /* Free group leader */

  return (0);
}
