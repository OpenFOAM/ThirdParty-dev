/* Copyright 2004,2007,2008,2010-2012,2014,2018,2019,2021,2023-2025 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : gout_c.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a result viewer.                **/
/**                This module contains the main function. **/
/**                                                        **/
/**   DATES      : # Version 2.0  : from : 06 oct 1994     **/
/**                                 to   : 23 dec 1994     **/
/**                # Version 3.0  : from : 14 jul 1995     **/
/**                                 to   : 11 oct 1995     **/
/**                # Version 3.1  : from : 27 mar 1996     **/
/**                                 to   : 03 apr 1996     **/
/**                # Version 3.2  : from : 02 dec 1996     **/
/**                                 to   : 05 jun 1998     **/
/**                # Version 3.3  : from : 29 may 1999     **/
/**                                 to   : 03 jun 1999     **/
/**                # Version 3.4  : from : 03 feb 2000     **/
/**                                 to   : 03 feb 2000     **/
/**                # Version 4.0  : from : 11 dec 2001     **/
/**                                 to   : 08 feb 2004     **/
/**                # Version 5.0  : from : 25 may 2007     **/
/**                                 to   : 25 may 2007     **/
/**                # Version 5.1  : from : 25 oct 2007     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 16 oct 2010     **/
/**                                 to   : 17 apr 2019     **/
/**                # Version 6.1  : from : 04 apr 2021     **/
/**                                 to   : 28 aug 2021     **/
/**                # Version 7.0  : from : 31 aug 2021     **/
/**                                 to   : 13 jul 2025     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes
*/

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "gout_c.h"
#include "gout_o.h"

/*
**  The static and global variables
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
File                        C_fileTab[C_FILENBR] = { /* The file array; public  */
                              { FILEMODER },
                              { FILEMODER },
                              { FILEMODER },
                              { FILEMODEW } };

static unsigned int         C_geoFlag = C_GEOFLAGDEFAULT; /* Geometry flag */

static const char *         C_usageList[] = {     /* Usage */
  "gout [<input source file> [<input geometry file> [<input mapping file> [<output picture file>]]]] <options>",
  "  -g<arguments>       : Geometry parameters :",
  "                          n  : do not read geometry data (matrix display)",
  "                          p  : permute Y and Z geometry dimensions",
  "                          r  : rotate geometry by 90 degrees",
  "  -h                  : Display this help",
  "  -mn                 : Do not read mapping data",
  "  -Oi[{<arguments>}]  : Open Inventor mesh file :",
  "                          c        : color output",
  "                          g        : gray level output",
  "                          r        : remove cut edges",
  "                          v        : view cut edges",
  "  -Om[{<arguments>}]  : PostScript matrix file :",
  "                          e        : EPSF-type output",
  "                          f        : full-page output",
  "  -Op[{<arguments>}]  : PostScript mesh file :",
  "                          c        : color output",
  "                          g        : gray level output",
  "                          e        : EPSF-type output",
  "                          f        : full-page output",
  "                          s        : short clipping (disks excluded)",
  "                          l        : large clipping (disks included)",
  "                          a        : avoid displaying disks",
  "                          d        : display disks",
  "                          r        : remove cut edges",
  "                          v        : view cut edges",
  "                          X=<val>  : maximum x clipping ratio (in [0.0;1.0])",
  "                          x=<val>  : minimum x clipping ratio",
  "                          Y=<val>  : maximum y clipping ratio",
  "                          y=<val>  : minimum y clipping ratio",
  "  -Ot[{<arguments>}]  : Tulip graph file :",
  "                          b        : b/w output",
  "                          c        : color output",
  "                          a        : avoid displaying disks",
  "                          d        : display disks",
  "                          r        : remove cut edges",
  "                          v        : view cut edges",
  "  -Ov[{<arguments>}]  : VTK mesh file :",
  "                          r        : remove cut edges",
  "                          v        : view cut edges",
  "  -V                  : Print program version and copyright",
  "",
  "Default option set is : -Ov{v}",
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
  C_Graph            grafdat;                     /* Source graph   */
  C_Geometry         geo;                         /* Graph geometry */
  C_Mapping          map;                         /* Result mapping */
  int                i, j;

  errorProg ("gout");

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
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'G' :                                /* Geometry parameters */
        case 'g' :
          if ((j = C_geoParse (&argv[i][2])) != 0)
            errorPrint ("main: error in geometry option string '%d'", j);
          break;
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (EXIT_SUCCESS);
        case 'M' :                                /* No-mapping flag */
        case 'm' :
          if (((argv[i][2] != 'N') && (argv[i][2] != 'n')) || (argv[i][3] != '\0'))
            errorPrint ("main: error in mapping option string '%s'", &argv[i][2]);
          C_filenamemapinp = "-";                 /* Default name to avoid opening   */
          C_filepntrmapinp = NULL;                /* NULL file pointer means no file */
          break;
        case 'O' :                                /* Output parameters */
        case 'o' :
          if ((j = outDrawParse (&argv[i][2])) != 0)
            errorPrint ("main: error in output option string (%d)", j);
          break;
        case 'V' :
          fprintf (stderr, "gout, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, SCOTCH_COPYRIGHT_STRING "\n");
          fprintf (stderr, SCOTCH_LICENSE_STRING "\n");
          return  (EXIT_SUCCESS);
        default :
          errorPrint ("main: Unprocessed option '%s'", argv[i]);
      }
    }
  }

  fileBlockOpen (C_fileTab, C_FILENBR);           /* Open all files */

  SCOTCH_graphInit (&grafdat.grafdat);            /* Create graph structure          */
  SCOTCH_graphLoad (&grafdat.grafdat, C_filepntrsrcinp, -1, 3); /* Read source graph */
  SCOTCH_graphData (&grafdat.grafdat, &grafdat.baseval,
                    &grafdat.vertnbr, &grafdat.verttax, &grafdat.vendtax, NULL, &grafdat.vlbltab,
                    &grafdat.edgenbr, &grafdat.edgetax, NULL);
  grafdat.verttax -= grafdat.baseval;
  grafdat.vendtax -= grafdat.baseval;
  grafdat.edgetax -= grafdat.baseval;

  C_geoInit (&geo, &grafdat);                     /* Create geometry structure */
  if (C_geoFlag & C_GEOFLAGUSE)                   /* If geometry is wanted     */
    C_geoLoad (&geo, C_filepntrgeoinp);           /* Read graph geometry       */

  C_mapInit (&map, &grafdat);                     /* Create mapping structure */
  C_mapLoad (&map, C_filepntrmapinp);             /* Read result mapping      */

  outDraw (&grafdat, &geo, &map, C_filepntrdatout); /* Build and write the output */

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  C_mapExit        (&map);                        /* Free data structures */
  C_geoExit        (&geo);
  SCOTCH_graphExit (&grafdat.grafdat);

  return (EXIT_SUCCESS);
}

/***********************************/
/*                                 */
/* These routines handle geometry. */
/*                                 */
/***********************************/

/* This routine parses the source graph
** option string.
** It returns:
** - 0  : if string successfully scanned.
** - 1  : if invalid options
** - 2  : if invalid option arguments.
** - 3  : if syntax error in string.
*/

int
C_geoParse (
const char * const          string)
{
  const char *        charptr;

  for (charptr = string; ; charptr ++) {
    switch (*charptr) {
      case 'N' :                                  /* Do not read geometry data */
      case 'n' :
        C_geoFlag &= ~C_GEOFLAGUSE;
        break;
      case 'P' :                                  /* Permute Y and Z */
      case 'p' :
        C_geoFlag |= C_GEOFLAGPERMUT;
        break;
      case 'R' :                                  /* If want to rotate */
      case 'r' :
        C_geoFlag |= C_GEOFLAGROTATE;
        break;
      case '\0' :
        return (0);
      default   :
        return (1);
    }
  }
}

/* This routine creates a geometry with
** respect to a given source graph.
** It returns:
** - VOID  : in all cases.
*/

void
C_geoInit (
C_Geometry * const          geomptr,
const C_Graph * const       grafptr)
{
  geomptr->grafptr = grafptr;
  geomptr->coortab = NULL;
}

/* This routine deletes a geometry.
** It returns:
** - VOID  : in all cases.
*/

void
C_geoExit (
C_Geometry * const          geomptr)
{
  if (geomptr->coortab != NULL)                   /* If there is a geometry array */
    memFree (geomptr->coortab);                   /* Free it                      */
}

/* This routine loads a mapping.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

/* This is the comparison function used by the
   quicksort algorithm, to sort by increasing
   labels.                                     */

static
int
C_geoLoad2 (
const C_VertSort * const    ver0ptr,
const C_VertSort * const    ver1ptr)
{
  return ((ver0ptr->lablval > ver1ptr->lablval) ?  1 : -1);
}

/** This is the loading routine. **/

int
C_geoLoad (
C_Geometry * const          geomptr,
FILE * const                fileptr)
{
  C_VertSort * restrict lablsorttab;              /* Pointer to graph labels sorting array      */
  SCOTCH_Num            lablvertnum;              /* Number of current vertex in label array    */
  C_GeoVert * restrict  geomcoortab;              /* Pointer to geometric data read from file   */
  C_VertSort *          geomsorttab;              /* Pointer to geometric data sorting array    */
  int                   geomflagval;              /* Flag set if geometric data sorted by label */
  int                   geomtypeval;              /* Type of geometry file                      */
  SCOTCH_Num            geomvertnbr;              /* Number of geometric coordinates in file    */
  SCOTCH_Num            geomvertnum;              /* Number of current vertex in geometry array */
  int                   o;

  const SCOTCH_Num                  baseval = geomptr->grafptr->baseval;
  const SCOTCH_Num                  vertnbr = geomptr->grafptr->vertnbr;
  const SCOTCH_Num * restrict const vlbltab = geomptr->grafptr->vlbltab;

  if ((geomptr->coortab == NULL) &&               /* Allocate geometry if necessary */
      ((geomptr->coortab = (C_GeoVert *) memAlloc (vertnbr * sizeof (C_GeoVert))) == NULL)) {
    errorPrint ("C_geoLoad: out of memory (1)");
    return (1);
  }

  if ((fscanf (fileptr, "%d" SCOTCH_NUMSTRING,    /* Read type and number of geometry items */
               &geomtypeval,
               &geomvertnbr) != 2) ||
      (geomtypeval < 1)           ||
      (geomtypeval > 3)           ||
      (geomvertnbr  < 1)) {
    errorPrint ("C_geoLoad: bad input (1)");
    return (1);
  }

  if (memAllocGroup ((void **) (void *)
                     &geomcoortab, (size_t) (geomvertnbr * sizeof (C_GeoVert)), /* Temporary array for geometry */
                     &geomsorttab, (size_t) (geomvertnbr * sizeof (C_VertSort)),
                     &lablsorttab, (size_t) (vertnbr     * sizeof (C_VertSort)), NULL) == NULL) {
    errorPrint ("C_geoLoad: out of memory (2)");
    return (1);
  }

  o = 0;
  geomflagval = 1;                                /* Assume geometry data sorted */
  switch (geomtypeval) {
    case 1 :                                      /* Load 2D coordinates array */
      for (geomvertnum = 0; (geomvertnum < geomvertnbr) && (o == 0); geomvertnum ++) {
        SCOTCH_Num          lablval;

        if (fscanf (fileptr, SCOTCH_NUMSTRING "%lf",
                    &lablval,
                    &geomcoortab[geomvertnum].x) != 2)
          o = 1;
        geomcoortab[geomvertnum].y =              /* No Y and Z coordinates */
        geomcoortab[geomvertnum].z = 0.0;
        geomsorttab[geomvertnum].lablval = lablval;
        geomsorttab[geomvertnum].vertnum = geomvertnum;

        if (C_geoFlag & C_GEOFLAGROTATE) {        /* Rotate picture if necessary */
          double              t;                  /* Temporary swap variable     */

          t                          = geomcoortab[geomvertnum].y;
          geomcoortab[geomvertnum].y = geomcoortab[geomvertnum].x;
          geomcoortab[geomvertnum].x = - t;
        }

        if ((geomvertnum > 0) &&                  /* Check if geometry data sorted */
            (geomsorttab[geomvertnum].lablval < geomsorttab[geomvertnum - 1].lablval))
          geomflagval = 0;                        /* Geometry data not sorted */
      }
      break;
    case 2 :                                      /* Load 2D coordinates array */
      for (geomvertnum = 0; (geomvertnum < geomvertnbr) && (o == 0); geomvertnum ++) {
        SCOTCH_Num          lablval;

        if (fscanf (fileptr, SCOTCH_NUMSTRING "%lf%lf",
                    &lablval,
                    &geomcoortab[geomvertnum].x,
                    &geomcoortab[geomvertnum].y) != 3)
          o = 1;
        geomcoortab[geomvertnum].z = 0.0;         /* No Z coordinate */
        geomsorttab[geomvertnum].lablval = lablval;
        geomsorttab[geomvertnum].vertnum = geomvertnum;

        if (C_geoFlag & C_GEOFLAGROTATE) {        /* Rotate picture if necessary */
          double              t;                  /* Temporary swap variable     */

          t                          = geomcoortab[geomvertnum].y;
          geomcoortab[geomvertnum].y = geomcoortab[geomvertnum].x;
          geomcoortab[geomvertnum].x = - t;
        }

        if ((geomvertnum > 0) &&                  /* Check if geometry data sorted */
            (geomsorttab[geomvertnum].lablval < geomsorttab[geomvertnum - 1].lablval))
          geomflagval = 0;                        /* Geometry data are not sorted */
      }
      break;
    case 3 :                                      /* Load 3D coordinates array */
      for (geomvertnum = 0; (geomvertnum < geomvertnbr) && (o == 0); geomvertnum ++) {
        SCOTCH_Num          lablval;

        if (fscanf (fileptr, SCOTCH_NUMSTRING "%lf%lf%lf",
                    &lablval,
                    &geomcoortab[geomvertnum].x,
                    &geomcoortab[geomvertnum].y,
                    &geomcoortab[geomvertnum].z) != 4)
          o = 1;
        geomsorttab[geomvertnum].lablval = lablval;
        geomsorttab[geomvertnum].vertnum = geomvertnum;

        if (C_geoFlag & C_GEOFLAGPERMUT) {        /* Rotate picture if necessary */
          double              t;                  /* Temporary swap variable     */

          t                          = geomcoortab[geomvertnum].z;
          geomcoortab[geomvertnum].z = geomcoortab[geomvertnum].y;
          geomcoortab[geomvertnum].y = t;
        }
        if ((geomvertnum > 0) &&                  /* Check if geometry data sorted */
            (geomsorttab[geomvertnum].lablval < geomsorttab[geomvertnum - 1].lablval))
          geomflagval = 0;                        /* Geometry data not sorted */
      }
      break;
    default :
      errorPrint ("C_geoLoad: invalid geometry type (%d)", geomtypeval);
      memFree    (geomcoortab);                   /* Free group leader */
      return (1);
  }
  if (o != 0) {
    errorPrint ("C_geoLoad: bad input (2)");
    memFree    (geomcoortab);                     /* Free group leader */
    return (1);
  }

  if (geomflagval != 1)                           /* If geometry data not sorted        */
    qsort ((char *) geomsorttab, geomvertnbr,     /* Sort sort area by ascending labels */
           sizeof (C_VertSort), (int (*) (const void *, const void *)) C_geoLoad2);
  for (geomvertnum = 1; geomvertnum < geomvertnbr; geomvertnum ++) { /* Check geometric data integrity */
    if (geomsorttab[geomvertnum].lablval == geomsorttab[geomvertnum - 1].lablval) {
      errorPrint ("C_geoLoad: duplicate vertex label");
      memFree    (geomcoortab);                   /* Free group leader */
      return (1);
    }
  }

  if (vlbltab != NULL) {                          /* If graph has vertex labels */
    SCOTCH_Num          lablvertnum;
    int                 lablflagval;              /* Flag set if graph data sorted by label */

    lablflagval = 1;                              /* Assume graph data sorted */
    for (lablvertnum = 0; lablvertnum < vertnbr; lablvertnum ++) {
      lablsorttab[lablvertnum].lablval = vlbltab[lablvertnum];
      lablsorttab[lablvertnum].vertnum = lablvertnum; /* Un-based index for matching */
      if ((lablvertnum > 0) &&                    /* Check if graph data sorted      */
          (lablsorttab[lablvertnum].lablval < lablsorttab[lablvertnum - 1].lablval))
        lablflagval = 0;                          /* Graph data not sorted */
    }
    if (lablflagval != 1)                         /* If graph data not sorted           */
      qsort ((char *) lablsorttab, vertnbr,       /* Sort sort area by ascending labels */
             sizeof (C_VertSort), (int (*) (const void *, const void *)) C_geoLoad2);
  }
  else {                                          /* Graph does not have vertex labels */
    for (lablvertnum = 0; lablvertnum < vertnbr; lablvertnum ++) {
      lablsorttab[lablvertnum].lablval = lablvertnum + baseval;
      lablsorttab[lablvertnum].vertnum = lablvertnum; /* Un-based index for matching */
    }
  }

  for (lablvertnum = geomvertnum = 0; lablvertnum < vertnbr; lablvertnum ++) { /* For all vertices in graph */
    while ((geomvertnum < geomvertnbr) && (geomsorttab[geomvertnum].lablval < lablsorttab[lablvertnum].lablval))
      geomvertnum ++;                             /* Search geometry vertex with same label */
    if ((geomvertnum >= geomvertnbr) ||           /* If label does not exist                */
        (geomsorttab[geomvertnum].lablval > lablsorttab[lablvertnum].lablval)) {
      errorPrint ("C_geoLoad: vertex geometry data not found for label '" SCOTCH_NUMSTRING "'",
                  lablsorttab[lablvertnum].lablval);
      memFree    (geomcoortab);                   /* Free group leader */
      return (1);
    }
    geomptr->coortab[lablsorttab[lablvertnum].vertnum] = geomcoortab[geomsorttab[geomvertnum ++].vertnum];
  }

  memFree (geomcoortab);                          /* Free group leader */

  return (0);
}

/***********************************/
/*                                 */
/* These routines handle mappings. */
/*                                 */
/***********************************/

/* This routine creates a mapping with
** respect to a given source graph.
** It returns:
** - VOID  : in all cases.
*/

void
C_mapInit (
C_Mapping * const           mappptr,
const C_Graph * const       grafptr)
{
  mappptr->grafptr = grafptr;
  mappptr->labltab = NULL;
}

/* This routine deletes a mapping.
** It returns:
** - VOID  : in all cases.
*/

void
C_mapExit (
C_Mapping * const           mappptr)
{
  if (mappptr->labltab != NULL)                   /* If there is a domain array */
    memFree (mappptr->labltab);                   /* Free it                    */
}

/* This routine loads a mapping.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
C_mapLoad (
C_Mapping * const           mappptr,
FILE * const                fileptr)
{
  C_VertSort * restrict lablsorttab;              /* Pointer to graph sorting array            */
  SCOTCH_Num            lablvertnum;              /* Number of current vertex in label array   */
  SCOTCH_Num * restrict mappparttab;              /* Pointer to mapping data read from file    */
  C_VertSort * restrict mappsorttab;              /* Pointer to mapping data sorting array     */
  int                   mappflagval;              /* Flag set if mapping data sorted by label  */
  SCOTCH_Num            mappvertnbr;              /* Number of mapping pairs in file           */
  SCOTCH_Num            mappvertnum;              /* Number of current vertex in mapping array */

  const SCOTCH_Num                  baseval = mappptr->grafptr->baseval;
  const SCOTCH_Num                  vertnbr = mappptr->grafptr->vertnbr;
  const SCOTCH_Num * restrict const vlbltab = mappptr->grafptr->vlbltab;

  if ((mappptr->labltab == NULL) &&               /* Allocate array if necessary */
      ((mappptr->labltab = (SCOTCH_Num *) memAlloc (vertnbr * sizeof (SCOTCH_Num))) == NULL)) {
    errorPrint ("C_mapLoad: out of memory (1)");
    return (1);
  }

  memset (mappptr->labltab, ~0, vertnbr * sizeof (SCOTCH_Num)); /* Pre-initialize mapping */

  if (fileptr == NULL)                            /* If stream is invalid */
    return (0);

  if ((fscanf (fileptr, SCOTCH_NUMSTRING,         /* Read number of mapping pairs */
               &mappvertnbr) != 1) ||
      (mappvertnbr < 1)) {
    errorPrint ("C_mapLoad: bad input (1)");
    return (1);
  }

  if (memAllocGroup ((void **) (void *)
                     &mappparttab, (size_t) (mappvertnbr * sizeof (SCOTCH_Num)),
                     &mappsorttab, (size_t) (mappvertnbr * sizeof (C_VertSort)),
                     &lablsorttab, (size_t) (vertnbr     * sizeof (C_VertSort)), NULL) == NULL) {
    errorPrint ("C_mapLoad: out of memory (2)");
    return (1);
  }

  mappflagval = 1;                                /* Assume mapping data sorted */
  for (mappvertnum = 0; mappvertnum < mappvertnbr; mappvertnum ++) {
    if (fscanf (fileptr, SCOTCH_NUMSTRING SCOTCH_NUMSTRING,
                &mappsorttab[mappvertnum].lablval,
                &mappparttab[mappvertnum]) != 2) {
      errorPrint ("C_mapLoad: bad input (2)");
      memFree    (mappparttab);                   /* Free group leader */
      return (1);
    }
    mappsorttab[mappvertnum].vertnum = mappvertnum;

    if ((mappvertnum > 0) &&                      /* Check if mapping data sorted */
        (mappsorttab[mappvertnum].lablval < mappsorttab[mappvertnum - 1].lablval))
      mappflagval = 0;                            /* Mapping data not sorted */
  }
  if (mappflagval != 1)                           /* If mapping data not sorted         */
      qsort ((char *) mappsorttab, mappvertnbr,   /* Sort sort area by ascending labels */
             sizeof (C_VertSort), (int (*) (const void *, const void *)) C_geoLoad2);
  for (mappvertnum = 1; mappvertnum < mappvertnbr; mappvertnum ++) { /* Check mapping data integrity */
    if (mappsorttab[mappvertnum].lablval == mappsorttab[mappvertnum - 1].lablval) {
      errorPrint ("C_mapLoad: duplicate vertex label");
      memFree    (mappparttab);                   /* Free group leader */
      return (1);
    }
  }

  if (vlbltab != NULL) {                          /* If graph has vertex labels */
    int                 lablflagval;              /* Flag set if graph data sorted by label */

    lablflagval = 1;                              /* Assume graph data sorted */
    for (lablvertnum = 0; lablvertnum < vertnbr; lablvertnum ++) {
      lablsorttab[lablvertnum].lablval = vlbltab[lablvertnum];
      lablsorttab[lablvertnum].vertnum = lablvertnum; /* Un-based index for matching */
      if ((lablvertnum > 0) &&                    /* Check if graph data sorted      */
          (lablsorttab[lablvertnum].lablval < lablsorttab[lablvertnum - 1].lablval))
        lablflagval = 0;                          /* Graph data not sorted */
    }
    if (lablflagval != 1)                         /* If graph data not sorted           */
      qsort ((char *) lablsorttab, vertnbr,       /* Sort sort area by ascending labels */
             sizeof (C_VertSort), (int (*) (const void *, const void *)) C_geoLoad2);
  }
  else {                                          /* Graph does not have vertex labels */
    for (lablvertnum = 0; lablvertnum < vertnbr; lablvertnum ++) {
      lablsorttab[lablvertnum].lablval = lablvertnum + baseval;
      lablsorttab[lablvertnum].vertnum = lablvertnum; /* Un-based index for matching */
    }
  }

  for (lablvertnum = mappvertnum = 0; lablvertnum < vertnbr; lablvertnum ++) { /* For all vertices in graph */
    while ((mappvertnum < mappvertnbr) && (mappsorttab[mappvertnum].lablval < lablsorttab[lablvertnum].lablval))
      mappvertnum ++;                             /* Search geometry vertex with same label */
    if ((mappvertnum >= mappvertnbr) ||           /* If label does not exist                */
        (mappsorttab[mappvertnum].lablval > lablsorttab[lablvertnum].lablval))
      continue;                                   /* This vertex has no related mapping data */
    mappptr->labltab[lablsorttab[lablvertnum].vertnum] = mappparttab[mappsorttab[mappvertnum ++].vertnum];
  }

  memFree (mappparttab);                          /* Free group leader */

  return (0);
}

/**************************************/
/*                                    */
/* The option string parsing routine. */
/*                                    */
/**************************************/

/* This routine parses an option string.
** It returns:
** - 0 : if string successfully scanned.
** - 1 : if invalid code name.
** - 2 : if invalid arguments for the code.
** - 3 : if syntax error in string.
*/

int
C_parse (
const C_ParseCode * const   codetab,              /* Pointer to the code array          */
const C_ParseArg * const    arguptr,              /* Pointer to the code argument array */
int * const                 codeptr,              /* Pointer to the code value to set   */
char * const                strgptr)              /* Pointer to the string to parse     */
{
  int                 codeval;                    /* Code found                       */
  size_t              codesiz;                    /* Code name length                 */
  char                argubuf[128];               /* Buffer for argument scanning     */
  size_t              argusiz;                    /* Length of the current argument   */
  char *              argbptr;                    /* Pointer to beginning of argument */
  char *              argeptr;                    /* Pointer to end of argument       */
  char *              argqptr;                    /* Position of the '=' character    */
  int                 codenum;

  codeval = -1;
  codesiz = 0;                                    /* No code recognized yet               */
  for (codenum = 0; codetab[codenum].nameptr != NULL; codenum ++) { /* For all code names */
    size_t              codetmp;

    codetmp = strlen (codetab[codenum].nameptr);
    if ((strncasecmp (strgptr, codetab[codenum].nameptr, codetmp) == 0) && /* Find longest matching code name */
        (codetmp > codesiz)) {
      codeval = codetab[codenum].codeval;
      codesiz = codetmp;
    }
  }
  if (codesiz == 0)                               /* If no code recognized */
    return (1);                                   /* Return error value    */
  *codeptr = codeval;                             /* Set code value        */

  argbptr = strgptr + codesiz;                    /* Point to the end of the code name */
  if (*argbptr == '{') {                          /* If there are arguments            */
    argbptr ++;                                   /* Point to argument beginning       */
    do {                                          /* For all arguments                 */
      int                 argunum;
      int                 flagval;

      argeptr = strpbrk (argbptr, ",}\0");        /* Search for the argument end  */
      if (*argeptr == '\0')                       /* If there is no end delimiter */
        return (3);                               /* Return syntax error value    */

      argusiz = ((argeptr - argbptr) < 127)       /* Get argument bounded length */
                ? (argeptr - argbptr)
                : 127;
      strncpy (argubuf, argbptr, argusiz);        /* Copy argument to the buffer    */
      argubuf[argusiz] = '\0';                    /* Mark the end of the argument   */
      argqptr = strpbrk (argubuf, "=");           /* Search for the '=' character   */
      if (argqptr != NULL)                        /* If it exists                   */
        *argqptr ++ = '\0';                       /* Turn it into a separating null */

      flagval = 0;                                /* No argument found yet                */
      for (argunum = 0; arguptr[argunum].nameptr != NULL; argunum ++) { /* Scan arguments */
        if ((arguptr[argunum].codeval == codeval) && /* If argument name found            */
            (strcmp (argubuf, arguptr[argunum].nameptr) == 0)) {
          flagval = 1;                            /* Argument found */
          break;
        }
      }
      if (flagval != 1)                           /* If invalid argument     */
        return (2);                               /* Return the proper value */

      if (arguptr[argunum].formptr != NULL) {     /* If there is a value to read    */
        if (argqptr == NULL)                      /* If none has been given however */
          return (2);                             /* Return the error value         */
        if (sscanf (argqptr,                      /* Try to read the argument value */
                    arguptr[argunum].formptr,
                    arguptr[argunum].avalptr) != 1)
          return (2);                             /* Return if error */
      }
      else {                                      /* If no value needed                */
        if (argqptr != NULL)                      /* If there is one however           */
          return (2);                             /* Return the error code             */
        *((char *) arguptr[argunum].avalptr) = argubuf[0]; /* Copy first argument char */
      }

      argbptr = argeptr + 1;                      /* Skip the processed argument         */
    } while (*argeptr != '}');                    /* Loop as long as there are arguments */
  }

  return ((*argbptr == '\0') ? 0 : 3);            /* Check if no extraneous characters */
}
