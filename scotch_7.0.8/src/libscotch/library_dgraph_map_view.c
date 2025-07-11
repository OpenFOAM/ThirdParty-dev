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
/**   NAME       : library_dgraph_map_view.c               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Clement BARTHELEMY                      **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the distri-  **/
/**                buted mapping routines of the libSCOTCH **/
/**                library.                                **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 26 jul 2008     **/
/**                                 to   : 11 aug 2010     **/
/**                # Version 6.0  : from : 29 nov 2012     **/
/**                                 to   : 25 apr 2018     **/
/**                # Version 7.0  : from : 27 aug 2019     **/
/**                                 to   : 19 sep 2024     **/
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
#include "ptscotch.h"

/**************************************/
/*                                    */
/* These routines compute statistics  */
/* on distributed mappings.           */
/*                                    */
/**************************************/

/*+ This routine writes distributed mapping
*** statistics to the given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_dgraphMapView (
SCOTCH_Dgraph * const         libgrafptr,
const SCOTCH_Dmapping * const libmappptr,
FILE * const                  stream)
{
  Anum                          tgtnbr;           /* Number of vertices in target architecture                 */
  Anum                          mapnbr;           /* Number of target terminal domains used in mapping         */
  Gnum                          mapmin;           /* Minimum load in domain (with respect to domain weight)    */
  Gnum                          mapmax;           /* Maximum load in domain (with respect to domain weight)    */
  double                        mapavg;           /* Average load in domain (with respect to domain weights)   */
  double                        mapdlt;           /* Standard deviation of mapping loads                       */
  double                        mapmmy;           /* Ratio of maximum over average                             */
  Anum                          ngbsum;           /* Total number of neighbors                                 */
  Anum                          ngbmin;           /* Minimum number of neighbors                               */
  Anum                          ngbmax;           /* Maximum number of neighbors                               */
  Anum                          cdstnum;          /* Current index in communication distance array             */
  Anum                          cdstmax;          /* Maximum distance reached in histogram array               */
  Gnum                          cdsttab[256];     /* Communication distance histogram array                    */
  Gnum                          cmlosum;          /* Communication load (sum of edge weights)                  */
  Gnum                          cmdisum;          /* Communication dilation (sum of distances)                 */
  Gnum                          cmexsum;          /* Communication expansion (sum of edge weights * distances) */

  const Dgraph * const grafptr = (Dgraph *) CONTEXTOBJECT (libgrafptr);

  if (SCOTCH_dgraphMapStat (libgrafptr, libmappptr, &tgtnbr, &mapnbr, &mapmin, &mapmax, &mapavg, &mapdlt,
                            &ngbsum, &ngbmin, &ngbmax, &cdstmax, cdsttab, &cmlosum, &cmdisum, &cmexsum) != 0) {
    errorPrint (STRINGIFY (SCOTCH_dgraphMapView) ": cannot compute dgraph map stats");
    return (1);
  }

  if (stream == NULL)
    return (0);

  fprintf (stream, "M\tProcessors " GNUMSTRING "/" GNUMSTRING "(%g)\n",
           (Gnum) mapnbr,
           (Gnum) tgtnbr,
           (double) mapnbr / (double) tgtnbr);

  mapmmy = (mapnbr != 0) ? (double) mapmax / (double) mapavg : 0.0L;
  fprintf (stream, "M\tTarget min=" GNUMSTRING "\tmax=" GNUMSTRING "\tavg=%g\tdlt=%g\tmaxavg=%g\n",
           (Gnum) mapmin,
           (Gnum) mapmax,
           mapavg,
           mapdlt,
           mapmmy);

  fprintf (stream, "M\tNeighbors min=" GNUMSTRING "\tmax=" GNUMSTRING "\tsum=" GNUMSTRING "\n",
           (Gnum) ngbmin,
           (Gnum) ngbmax,
           (Gnum) ngbsum);

  fprintf (stream, "M\tCommDilat=%f\t(" GNUMSTRING ")\n",
           (double) cmdisum / (double) (grafptr->edgeglbnbr / 2),
           (Gnum) cmdisum);
  fprintf (stream, "M\tCommExpan=%f\t(" GNUMSTRING ")\n",
           ((cmlosum == 0) ? (double) 0.0L
                           : (double) cmexsum / (double) cmlosum),
           (Gnum) cmexsum);
  fprintf (stream, "M\tCommCutSz=%f\t(" GNUMSTRING ")\n",
           ((cmlosum == 0) ? (double) 0.0L
                           : (double) (cmlosum - cdsttab[0]) / (double) cmlosum),
           (Gnum) (cmlosum - cdsttab[0]));
  fprintf (stream, "M\tCommDelta=%f\n",
           (((double) cmlosum  * (double) cmdisum) == 0.0L)
           ? (double) 0.0L
           : ((double) grafptr->edgeglbnbr) / (double) (cmlosum * 2));

  for (cdstnum = 0; cdstnum <= cdstmax; cdstnum ++) /* Print distance histogram */
    fprintf (stream, "M\tCmlosum[" ANUMSTRING "]=%f\n",
             (Anum) cdstnum,
             (double) cdsttab[cdstnum] / (double) cmlosum);

  return (0);
}
