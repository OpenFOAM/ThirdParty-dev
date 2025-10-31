/* Copyright 2025 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dmesh.h                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the distributed source mesh         **/
/**                structure.                              **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 09 aug 2025     **/
/**                                 to   : 11 aug 2925     **/
/**                                                        **/
/************************************************************/

/*
** The defines.
*/

#define SCOTCH_DMESH_H

/* Mesh flags. */

#define DMESHNONE                   0x0000        /* No options set */

#define DMESHFREEPRIV               0x0001        /* Set if private arrays freed on exit */
#define DMESHFREECOMM               0x0002        /* MPI communicator has to be freed    */
#define DMESHFREETABS               0x0004        /* Set if local arrays freed on exit   */

#define DMESHBITSUSED               0x0007        /* Significant bits for plain distributed mesh routines               */
#define DMESHBITSNOTUSED            0x0008        /* Value above which bits not used by plain distributed mesh routines */

/*
** The type and structure definitions.
*/

/* The mesh basic types, which must be signed. */

#ifndef GNUMMAX                                   /* If mesh.h not included     */
typedef INT                 Gnum;                 /* Vertex or edge number      */
#define GNUMMAX                     (INTVALMAX)   /* Maximum Gnum value         */
#define GNUMMIN                     (-GNUMMAX - 1) /* Minimum signed Gnum value */
#define GNUMSTRING                  INTSTRING     /* String to printf a Gnum    */
#endif /* GNUMMAX */

#ifndef GNUM_MPI
#define GNUM_MPI                    COMM_INT      /* MPI type for Gnum is MPI type for INT */
#endif /* GNUM_MPI */

#ifndef GRAPHPART_MPI
#define GRAPHPART_MPI               COMM_BYTE     /* Raw byte type for graph parts */
#endif /* GRAPHPART_MPI */

/*+ The mesh flag type. +*/

typedef unsigned int DmeshFlag;                  /*+ Mesh property flags +*/

/*+ The vertex part type, in compressed form. From mesh.h +*/

#ifndef SCOTCH_MESH_H
typedef byte MeshPart;
#endif /* SCOTCH_MESH_H */

/* The distributed mesh structure. */

typedef struct Dmesh_ {
  DmeshFlag                 flagval;              /*+ Mesh properties                                             +*/
  Gnum                      baseval;              /*+ Base index for edge/vertex arrays                           +*/
  Gnum                      velmglbnbr;           /*+ Global number of element vertices                           +*/
  Gnum                      velmglbmax;           /*+ Maximum number of local element vertices over all processes +*/
  Gnum                      velmlocnbr;           /*+ Local number of element vertices                            +*/
  Gnum *                    velmloctab;           /*+ Local element vertex beginning index array                  +*/
  Gnum                      eelmglbnbr;           /*+ Global number of element-to-node edges                      +*/
  Gnum                      eelmglbmax;           /*+ Maximum number of local element edges over all processes    +*/
  Gnum                      eelmlocnbr;           /*+ Local number of element vertices                            +*/
  Gnum *                    eelmloctab;           /*+ Edge array holding global neighbor numbers                  +*/
  Gnum                      vnodglbnbr;           /*+ Global number of node vertices                              +*/
  MPI_Comm                  proccomm;             /*+ Mesh communicator                                           +*/
  int                       procglbnbr;           /*+ Number of processes sharing mesh data                       +*/
  int                       proclocnum;           /*+ Number of this process                                      +*/
  Gnum *                    prelvrttab;           /*+ Per-process element vertex distribution array               +*/
} Dmesh;

/*
** The function prototypes.
*/

int                         dmeshInit           (Dmesh * const, MPI_Comm);
void                        dmeshExit           (Dmesh * const);
void                        dmeshFree           (Dmesh * const);
#ifdef SCOTCH_GRAPH_H
int                         dmeshLoad           (Dmesh * const, FILE * const, const Gnum, const GraphLoadFlag);
#endif /* SCOTCH_GRAPH_H */
int                         dmeshBuildAdm       (Dmesh * const, const Gnum, const Gnum, Gnum * const, const Gnum, Gnum * const, const Gnum);
#ifdef SCOTCH_DGRAPH_H
int                         dmeshDgraphDual     (const Dmesh * restrict const, Dgraph * restrict const, const Gnum);
#endif /* SCOTCH_DGRAPH_H */
