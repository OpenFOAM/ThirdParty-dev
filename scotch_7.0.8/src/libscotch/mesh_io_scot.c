/* Copyright 2004,2007,2008,2010,2019,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : mesh_io_scot.c                          **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the I/O routines   **/
/**                for handling the Scotch mesh format.    **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 19 jan 2004     **/
/**                                 to   : 19 jan 2004     **/
/**                # Version 5.0  : from : 13 sep 2006     **/
/**                                 to   : 27 feb 2008     **/
/**                # Version 5.1  : from : 11 aug 2010     **/
/**                                 to   : 11 aug 2010     **/
/**                # Version 6.0  : from : 28 apr 2019     **/
/**                                 to   : 28 apr 2019     **/
/**                # Version 7.0  : from : 20 jan 2023     **/
/**                                 to   : 20 jan 2023     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "geom.h"
#include "graph.h"
#include "mesh.h"
#include "mesh_io_scot.h"

/* This routine loads the geometrical mesh
** in the Scotch graph format, and allocates
** the proper structures.
** - 0   : on success.
** - !0  : on error.
*/

int
meshGeomLoadScot (
Mesh * restrict const       meshptr,              /* Graph to load    */
Geom * restrict const       geomptr,              /* Geometry to load */
FILE * const                filesrcptr,           /* Topological data */
FILE * const                filegeoptr,           /* No use           */
const char * const          dataptr)              /* No use           */
{
  double * restrict           coorfiletab;        /* Pointer to geometric data read from file   */
  MeshGeomScotSort * restrict coorsorttab;        /* Pointer to geometric data sorting array    */
  int                         coorsortflag;       /* Flag set if geometric data sorted by label */
  Gnum                        coornbr;            /* Number of geometric coordinates in file    */
  Gnum                        coornum;            /* Number of current coordinate               */
  MeshGeomScotSort * restrict vnodsorttab;        /* Pointer to graph sorting array             */
  int                         vnodsortflag;       /* Flag set if graph data sorted by label     */
  Gnum                        vnodnum;            /* Current graph vertex                       */
  Gnum                        dimnnbr;            /* Dimension of geometry file                 */
  int                         o;

  if (filesrcptr != NULL) {
    if (meshLoad (meshptr, filesrcptr, -1) != 0)
      return (1);
  }

  if (filegeoptr == NULL)
    return (0);

  if ((intLoad (filegeoptr, &dimnnbr) != 1) ||    /* Read type and number of geometry items */
      (intLoad (filegeoptr, &coornbr) != 1) ||
      (dimnnbr < 1)                         ||
      (dimnnbr > 3)) {
    errorPrint ("meshGeomLoadScot: bad input (1)");
    return     (1);
  }
  if ((filesrcptr != NULL) && (meshptr->vnodnbr != coornbr)) {
    errorPrint ("meshGeomLoadScot: inconsistent number of vertices");
    return     (1);
  }

  if (meshptr->vnodnbr == 0)
    return (0);

  if ((geomptr->geomtab == NULL) &&               /* Allocate geometry if necessary */
      ((geomptr->geomtab = (double *) memAlloc (meshptr->vnodnbr * dimnnbr * sizeof (double))) == NULL)) {
    errorPrint ("meshGeomLoadScot: out of memory (1)");
    return     (1);
  }

  if (memAllocGroup ((void **)
                     &coorfiletab, (size_t) (coornbr * dimnnbr * sizeof (double)),
                     &coorsorttab, (size_t) (coornbr           * sizeof (MeshGeomScotSort)),
                     &vnodsorttab, (size_t) (meshptr->vnodnbr  * sizeof (MeshGeomScotSort)), NULL) == NULL) {
    errorPrint ("meshGeomLoadScot: out of memory (2)");
    return     (1);
  }

  o = 0;
  coorsortflag = 1;                               /* Assume geometry data sorted */
  for (coornum = 0; (o == 0) && (coornum < coornbr); coornum ++) {
    Gnum                vlblnum;

    o = 1 - intLoad (filegeoptr, &vlblnum);
    coorsorttab[coornum].labl = vlblnum;
    coorsorttab[coornum].num  = coornum;
    if ((coornum > 0) &&                          /* Check if geometry data sorted */
        (coorsorttab[coornum].labl < coorsorttab[coornum - 1].labl))
      coorsortflag = 0;                           /* Geometry data not sorted */

    o |= 1 - fscanf (filegeoptr, "%lf",           /* Read X coordinate */
                     &coorfiletab[coornum * dimnnbr]);
    if (dimnnbr > 1) {
      o |= 1 - fscanf (filegeoptr, "%lf",         /* Read Y coordinate */
                       &coorfiletab[(coornum * dimnnbr) + 1]);
      if (dimnnbr > 2)
        o |= 1 - fscanf (filegeoptr, "%lf",       /* Read Z coordinate */
                         &coorfiletab[(coornum * dimnnbr) + 2]);
    }
  }
  if (o != 0) {
    errorPrint ("meshGeomLoadScot: bad input (2)");
    memFree    (coorfiletab);                     /* Free group leader */
    return     (1);
  }

  if (coorsortflag != 1)                          /* If geometry data not sorted        */
    intSort2asc1 (coorsorttab, coornbr);          /* Sort sort area by ascending labels */
  for (coornum = 1; coornum < coornbr; coornum ++) { /* Check geometric data integrity  */
    if (coorsorttab[coornum].labl == coorsorttab[coornum - 1].labl) {
      errorPrint ("meshGeomLoadScot: duplicate vertex label");
      memFree    (coorfiletab);                   /* Free group leader */
      return     (1);
    }
  }

  if (meshptr->vlbltax != NULL) {                 /* If graph has vertex labels */
    vnodsortflag = 1;                             /* Assume graph data sorted   */
    for (vnodnum = 0; vnodnum < meshptr->vnodnbr; vnodnum ++) {
      vnodsorttab[vnodnum].labl = meshptr->vlbltax[vnodnum + meshptr->baseval];
      vnodsorttab[vnodnum].num  = vnodnum;
      if ((vnodnum > 0) &&                        /* Check if graph data sorted */
          (vnodsorttab[vnodnum].labl < vnodsorttab[vnodnum - 1].labl))
        vnodsortflag = 0;                         /* Graph data not sorted */
    }
    if (vnodsortflag != 1)                        /* If graph data not sorted             */
      intSort2asc1 (vnodsorttab, meshptr->vnodnbr); /* Sort sort area by ascending labels */
  }
  else {                                          /* Graph does not have vertex labels */
    for (vnodnum = 0; vnodnum < meshptr->vnodnbr; vnodnum ++)
      vnodsorttab[vnodnum].labl =
      vnodsorttab[vnodnum].num  = vnodnum;
  }

  for (coornum = vnodnum = 0; vnodnum < meshptr->vnodnbr; vnodnum ++) { /* For all vertices in graph */
    while ((coornum < coornbr) && (coorsorttab[coornum].labl < vnodsorttab[vnodnum].labl))
      coornum ++;                                 /* Search geometry vertex with same label                           */
    if ((coornum >= coornbr) || (coorsorttab[coornum].labl > vnodsorttab[vnodnum].labl)) { /* If label does not exist */
      errorPrint ("meshGeomLoadScot: vertex geometry data not found (%d)",
                  vnodsorttab[vnodnum].labl);
      memFree    (coorfiletab);                   /* Free group leader */
      return     (1);
    }
    memCpy (&geomptr->geomtab[vnodsorttab[vnodnum].num * dimnnbr], &coorfiletab[coorsorttab[coornum ++].num * dimnnbr], dimnnbr * sizeof (double));
  }

  memFree (coorfiletab);                          /* Free group leader */

  return (0);
}

/* This routine saves the source mesh
** in the Scotch mesh and geometry formats.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
meshGeomSaveScot (
const Mesh * restrict const   meshptr,            /* Mesh to save     */
const Geom * restrict const   geomptr,            /* Geometry to save */
FILE * const                  filesrcptr,         /* Topological data */
FILE * const                  filegeoptr,         /* No use           */
const char * const            dataptr)            /* No use           */
{
  Gnum              vnodnum;
  int               dimnnbr;
  int               o;

  if (filesrcptr != NULL) {
    if (meshSave (meshptr, filesrcptr) != 0)      /* Save mesh structural data */
      return (1);
  }

  dimnnbr = geomptr->dimnnbr;

  o = 0;
  if (geomptr->geomtab != NULL) {                 /* If geometrical data present     */
    o = (fprintf (filegeoptr, GNUMSTRING "\n" GNUMSTRING "\n", /* Output file header */
                  (Gnum) geomptr->dimnnbr,
                  (Gnum) meshptr->vnodnbr) == EOF);

    switch (dimnnbr) {                            /* Output geometry data */
      case 1 :
        for (vnodnum = meshptr->vnodbas; (o == 0) && (vnodnum < meshptr->vnodnnd); vnodnum ++)
          o |= (fprintf (filegeoptr, GNUMSTRING "\t%lf\n",
                         (Gnum) ((meshptr->vlbltax != NULL) ? meshptr->vlbltax[vnodnum] : vnodnum),
                         (double) geomptr->geomtab[(vnodnum - meshptr->vnodbas) * dimnnbr]) == EOF);
        break;
      case 2 :
        for (vnodnum = meshptr->vnodbas; (o == 0) && (vnodnum < meshptr->vnodnnd); vnodnum ++)
          o |= (fprintf (filegeoptr, GNUMSTRING "\t%lf\t%lf\n",
                         (Gnum) ((meshptr->vlbltax != NULL) ? meshptr->vlbltax[vnodnum] : vnodnum),
                         (double) geomptr->geomtab[(vnodnum - meshptr->vnodbas) * dimnnbr],
                         (double) geomptr->geomtab[(vnodnum - meshptr->vnodbas) * dimnnbr + 1]) == EOF);
        break;
      case 3 :
        for (vnodnum = meshptr->vnodbas; (o == 0) && (vnodnum < meshptr->vnodnnd); vnodnum ++)
          o |= (fprintf (filegeoptr, GNUMSTRING "\t%lf\t%lf\t%lf\n",
                         (Gnum) ((meshptr->vlbltax != NULL) ? meshptr->vlbltax[vnodnum] : vnodnum),
                         (double) geomptr->geomtab[(vnodnum - meshptr->vnodbas) * dimnnbr],
                         (double) geomptr->geomtab[(vnodnum - meshptr->vnodbas) * dimnnbr + 1],
                         (double) geomptr->geomtab[(vnodnum - meshptr->vnodbas) * dimnnbr + 2]) == EOF);
        break;
#ifdef SCOTCH_DEBUG_MESH2
      default :
        errorPrint ("meshGeomSaveScot: invalid geometry type");
        return     (1);
#endif /* SCOTCH_DEBUG_MESH2 */
    }

    if (o != 0) {
      errorPrint ("meshGeomSaveScot: bad output");
    }
  }

  return (o);
}
