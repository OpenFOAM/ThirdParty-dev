/* Copyright 2004,2007,2008,2011,2012,2018,2021,2023,2025 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : gout_o.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a result viewer.                **/
/**                This module contains output routines.   **/
/**                                                        **/
/**   DATES      : # Version 2.0  : from : 07 oct 1994     **/
/**                                 to   : 23 dec 1994     **/
/**                # Version 3.0  : from : 14 jul 1995     **/
/**                                 to   : 03 oct 1995     **/
/**                # Version 3.1  : from : 28 mar 1996     **/
/**                                 to   : 03 jun 1996     **/
/**                # Version 3.2  : from : 02 dec 1996     **/
/**                                 to   : 05 jun 1998     **/
/**                # Version 3.3  : from : 29 may 1999     **/
/**                                 to   : 03 jun 1999     **/
/**                # Version 4.0  : from : 11 dec 2001     **/
/**                                 to   : 11 dec 2001     **/
/**                # Version 5.0  : from : 25 may 2007     **/
/**                                 to   : 18 jun 2007     **/
/**                # Version 5.1  : from : 25 oct 2007     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 21 may 2018     **/
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

static O_OutParam           O_outParam = {        /* Parameter structure        */
                              O_OUTTYPEVTKMESH,   /* Default output type        */
                              { 'c', 'v' },       /* OpenInventor mesh defaults */
                              { 'f' },            /* PostScript matrix defaults */
                              { 'f', 'g',         /* PostScript mesh defaults   */
                                'v', 'd',
                                's',
                                { { 0.0, 0.0 } },
                                { { 1.0, 1.0 } } },
                              { 'c', 'v', 'a' },  /* Tulip graph defaults       */
                              { 'v' } };          /* VTK mesh defaults          */

static C_ParseCode          O_outList[] = {       /* Output code list */
                              { O_OUTTYPEINVMESH, "i"  },
                              { O_OUTTYPEPOSMATR, "m"  },
                              { O_OUTTYPEPOSMESH, "p"  },
                              { O_OUTTYPETULMESH, "t"  },
                              { O_OUTTYPEVTKMESH, "v"  },
                              { O_OUTTYPENBR,     NULL } };

static C_ParseArg           O_outArg[] = {        /* Output type argument list */
                              { "c",  O_OUTTYPEINVMESH, NULL,  &O_outParam.InvMesh.coloval   },
                              { "g",  O_OUTTYPEINVMESH, NULL,  &O_outParam.InvMesh.coloval   },
                              { "r",  O_OUTTYPEINVMESH, NULL,  &O_outParam.InvMesh.edgeval   },
                              { "v",  O_OUTTYPEINVMESH, NULL,  &O_outParam.InvMesh.edgeval   },
                              { "e",  O_OUTTYPEPOSMATR, NULL,  &O_outParam.PosMatr.typeval   },
                              { "f",  O_OUTTYPEPOSMATR, NULL,  &O_outParam.PosMatr.typeval   },
                              { "c",  O_OUTTYPEPOSMESH, NULL,  &O_outParam.PosMesh.coloval   },
                              { "g",  O_OUTTYPEPOSMESH, NULL,  &O_outParam.PosMesh.coloval   },
                              { "e",  O_OUTTYPEPOSMESH, NULL,  &O_outParam.PosMesh.typeval   },
                              { "f",  O_OUTTYPEPOSMESH, NULL,  &O_outParam.PosMesh.typeval   },
                              { "l",  O_OUTTYPEPOSMESH, NULL,  &O_outParam.PosMesh.clipval   },
                              { "s",  O_OUTTYPEPOSMESH, NULL,  &O_outParam.PosMesh.clipval   },
                              { "a",  O_OUTTYPEPOSMESH, NULL,  &O_outParam.PosMesh.diskval   },
                              { "d",  O_OUTTYPEPOSMESH, NULL,  &O_outParam.PosMesh.diskval   },
                              { "r",  O_OUTTYPEPOSMESH, NULL,  &O_outParam.PosMesh.edgeval   },
                              { "v",  O_OUTTYPEPOSMESH, NULL,  &O_outParam.PosMesh.edgeval   },
                              { "x",  O_OUTTYPEPOSMESH, "%lf", &O_outParam.PosMesh.pminval.x },
                              { "X",  O_OUTTYPEPOSMESH, "%lf", &O_outParam.PosMesh.pmaxval.x },
                              { "y",  O_OUTTYPEPOSMESH, "%lf", &O_outParam.PosMesh.pminval.y },
                              { "Y",  O_OUTTYPEPOSMESH, "%lf", &O_outParam.PosMesh.pmaxval.y },
                              { "b",  O_OUTTYPETULMESH, NULL,  &O_outParam.TulMesh.coloval   },
                              { "c",  O_OUTTYPETULMESH, NULL,  &O_outParam.TulMesh.coloval   },
                              { "r",  O_OUTTYPETULMESH, NULL,  &O_outParam.TulMesh.edgeval   },
                              { "v",  O_OUTTYPETULMESH, NULL,  &O_outParam.TulMesh.edgeval   },
                              { "a",  O_OUTTYPETULMESH, NULL,  &O_outParam.TulMesh.diskval   },
                              { "d",  O_OUTTYPETULMESH, NULL,  &O_outParam.TulMesh.diskval   },
                              { "r",  O_OUTTYPEVTKMESH, NULL,  &O_outParam.VtkMesh.edgeval   },
                              { "v",  O_OUTTYPEVTKMESH, NULL,  &O_outParam.VtkMesh.edgeval   },
                              { NULL, O_OUTTYPENBR,     "",    NULL                          } };

static double               outcolorcoltab[16][3] = { /* Color list */
                              { 1.00, 0.00, 0.00 }, /* Red          */
                              { 0.00, 1.00, 0.00 }, /* Green        */
                              { 1.00, 1.00, 0.00 }, /* Yellow       */
                              { 0.00, 0.00, 1.00 }, /* Blue         */
                              { 1.00, 0.00, 1.00 }, /* Magenta      */
                              { 0.00, 1.00, 1.00 }, /* Cyan         */
                              { 1.00, 0.50, 0.20 }, /* Orange       */
                              { 0.30, 0.55, 0.00 }, /* Olive        */
                              { 0.72, 0.47, 0.47 }, /* Dark pink    */
                              { 0.33, 0.33, 0.81 }, /* Sea blue     */
                              { 1.00, 0.63, 0.63 }, /* Pink         */
                              { 0.62, 0.44, 0.65 }, /* Violet       */
                              { 0.60, 0.80, 0.70 }, /* Pale green   */
                              { 0.47, 0.20, 0.00 }, /* Brown        */
                              { 0.00, 0.68, 0.68 }, /* Turquoise    */
                              { 0.81, 0.00, 0.40 } }; /* Purple     */

static double               outcolorblwtab[8][3] = { /* Grey list */
                              { 1.00, 1.00, 1.00 },
                              { 0.20, 0.20, 0.20 },
                              { 0.50, 0.50, 0.50 },
                              { 0.80, 0.80, 0.80 },
                              { 0.30, 0.30, 0.30 },
                              { 0.90, 0.90, 0.90 },
                              { 0.40, 0.40, 0.40 },
                              { 0.70, 0.70, 0.70 } };

/****************************************/
/*                                      */
/* This is the color selection routine. */
/*                                      */
/****************************************/

void
outColorBlw (
const SCOTCH_Num            labl,
double                      color[])
{
  if (labl == (-1)) {
    color[0] =
    color[1] =
    color[2] = 1.0L;
  }
  else {
    color[0] = (double) outcolorblwtab[labl % 8][0];
    color[1] = (double) outcolorblwtab[labl % 8][1];
    color[2] = (double) outcolorblwtab[labl % 8][2];
  }
}

void
outColorColor (
const SCOTCH_Num            labl,
double                      color[])
{
  if (labl == (-1)) {
    color[0] =
    color[1] =
    color[2] = 1.0L;
  }
  else {
    color[0] = (double) outcolorcoltab[labl % 16][0];
    color[1] = (double) outcolorcoltab[labl % 16][1];
    color[2] = (double) outcolorcoltab[labl % 16][2];
  }
}

/****************************/
/*                          */
/* The main output routine. */
/*                          */
/****************************/

/* This routine parses the output
** option string.
** It returns:
** - 0  : if string successfully scanned.
** - 1  : if invalid options
** - 2  : if invalid option arguments.
** - 3  : if syntax error in string.
*/

int
outDrawParse (
char * const                strgptr)
{
  return (C_parse (O_outList, O_outArg, (int * const) (void *) &O_outParam.typeval, strgptr));
}

/* This routine is the generic output call.
** It returns:
** - VOID  : in all cases.
*/

void
outDraw (
const C_Graph * const       grafptr,              /* Graph structure */
const C_Geometry * const    geomptr,              /* Graph geometry  */
const C_Mapping * const     mappptr,              /* Result mapping  */
FILE * const                fileptr)              /* Output stream   */
{
  switch (O_outParam.typeval) {
    case O_OUTTYPEINVMESH :                       /* Mesh OpenInventor output type */
      outDrawInvMesh (grafptr, geomptr, mappptr, fileptr);
      break;
    case O_OUTTYPEPOSMATR :                       /* Matrix PostScript output type */
      outDrawPosMatr (grafptr, geomptr, mappptr, fileptr);
      break;
    case O_OUTTYPEPOSMESH :                       /* Mesh PostScript output type */
      outDrawPosMesh (grafptr, geomptr, mappptr, fileptr);
      break;
    case O_OUTTYPETULMESH :                       /* Mesh Tulip output type */
      outDrawTulMesh (grafptr, geomptr, mappptr, fileptr);
      break;
    case O_OUTTYPEVTKMESH :                       /* Mesh Tulip output type */
      outDrawVtkMesh (grafptr, geomptr, mappptr, fileptr);
      break;
    default :
      errorPrint ("outDraw: invalid output method '%d'", O_outParam.typeval);
  }
}

/****************************************/
/*                                      */
/* This is the Inventor output routine. */
/*                                      */
/****************************************/

int
outDrawInvMesh (
const C_Graph * const       grafptr,              /* Graph structure, sorted by vertex index */
const C_Geometry * const    geomptr,              /* Graph geometry, sorted by vertex label  */
const C_Mapping * const     mappptr,              /* Result mapping, sorted by vertex label  */
FILE * const                fileptr)              /* Output stream                           */
{
  void             (* colfptr) (const SCOTCH_Num, double[]); /* Color routine    */
  O_InvMeshPath * restrict  pathtax;              /* Array of path building data */
  SCOTCH_Num * restrict     indxtab;              /* Array of vertex indices     */
  SCOTCH_Num                indxnbr;              /* Number of indices           */
  SCOTCH_Num                indxnum;
  time_t                    timedat;              /* Creation time               */
  double                    colotab[3];           /* Vertex color                */
  SCOTCH_Num                vertnum;

  const SCOTCH_Num                  baseval = grafptr->baseval;
  const SCOTCH_Num                  vertnnd = grafptr->vertnbr + grafptr->baseval;
  const SCOTCH_Num * restrict const verttax = grafptr->verttax;
  const SCOTCH_Num * restrict const vendtax = grafptr->vendtax;
  const SCOTCH_Num * restrict const edgetax = grafptr->edgetax;
  const C_GeoVert * restrict const  coortax = geomptr->coortab - baseval;
  const SCOTCH_Num * restrict const parttax = mappptr->labltab - baseval;

  if (geomptr->coortab == NULL) {
    errorPrint ("outDrawInvMesh: geometry not provided");
    return (1);
  }

  time (&timedat);                                /* Get current time */

  colfptr = (O_outParam.InvMesh.coloval == 'c') ? outColorColor : outColorBlw; /* Select color output routine */

  if (memAllocGroup ((void **) (void *)
                     &indxtab, (size_t) ((grafptr->edgenbr / 2) * 3 * sizeof (SCOTCH_Num)),
                     &pathtax, (size_t) (grafptr->vertnbr           * sizeof (O_InvMeshPath)), NULL) == NULL) {
    errorPrint ("outDrawInvMesh: out of memory");
    return (1);
  }
  pathtax -= baseval;                             /* Base path array */

  for (vertnum = baseval; vertnum < vertnnd; vertnum ++) { /* For all vertices */
    SCOTCH_Num          pathnbr;
    SCOTCH_Num          partval;
    SCOTCH_Num          edgennd;
    SCOTCH_Num          edgenum;

    partval = parttax[vertnum];
    edgenum = verttax[vertnum];
    pathtax[vertnum].edgenum = edgenum;           /* Record start exit index */
    pathnbr = 0;                                  /* No output path yet      */
    for (edgennd = vendtax[vertnum]; edgenum < edgennd; edgenum ++) {
      SCOTCH_Num          vertend;

      vertend = edgetax[edgenum];
      if ((vertend > vertnum) &&                  /* If end vertex can be a next hop */
          ((O_outParam.InvMesh.edgeval != 'r') || /* And this edge can be drawn      */
           (parttax[vertend] == partval)))
        pathnbr ++;                               /* One more path to vertices of higher indices */
    }
    pathtax[vertnum].pathnbr = pathnbr;           /* Record number of output paths */
  }

  indxnbr = 0;                                    /* No indices yet   */
  for (vertnum = baseval; vertnum < vertnnd; ) {  /* For all vertices */
    SCOTCH_Num          verttmp;

    if (pathtax[vertnum].pathnbr == 0) {          /* If no output path for this vertex */
      vertnum ++;                                 /* Skip to next vertex               */
      continue;
    }

    verttmp = vertnum;                            /* Begin with this vertex         */
    indxtab[indxnbr ++] = verttmp;                /* Add it to the current path     */
    do {                                          /* Build path from current vertex */
      SCOTCH_Num          parttmp;
      SCOTCH_Num          vertend;
      SCOTCH_Num          edgennd;
      SCOTCH_Num          edgetmp;

      parttmp = parttax[verttmp];
      edgetmp = pathtax[verttmp].edgenum;
      edgennd = vendtax[verttmp];
      do {                                        /* Search for first valid output vertex */
        vertend = edgetax[edgetmp];
        if ((vertend > verttmp) &&                /* If end vertex can be next hop */
            ((O_outParam.InvMesh.edgeval != 'r') || /* And this edge can be drawn  */
             (parttax[vertend] == parttmp)))
          break;
      } while (++ edgetmp < edgennd);

      pathtax[verttmp].pathnbr --;                /* One less output path remaining    */
      pathtax[verttmp].edgenum = edgetmp + 1;     /* Search from the next position     */
      verttmp = vertend;                          /* Go-on from end vertex             */
      indxtab[indxnbr ++] = verttmp;              /* Add end vertex to current path    */
    } while (pathtax[verttmp].pathnbr > 0);       /* As long as there is an exit route */

    indxtab[indxnbr ++] = -1;                     /* Mark end of current path */
  }

  fprintf (fileptr, "#Inventor V2.0 ascii\n");    /* Write header */
  fprintf (fileptr, "#Title: %s %s %s\n",
           C_filenamesrcinp, C_filenamegeoinp, C_filenamemapinp);
  fprintf (fileptr, "#Creator: Scotch/gout\n");
  fprintf (fileptr, "#CreationDate: %s", ctime (&timedat));

  if (indxnbr == 0)                               /* If nothing to write */
    return (0);

  fprintf (fileptr, "Separator {\n");
  fprintf (fileptr, "  LightModel {\n    model\t\tBASE_COLOR\n  }\n");
  fprintf (fileptr, "  DrawStyle {\n    style\t\tLINES\n  }\n");
  fprintf (fileptr, "  MaterialBinding {\n    value\t\tPER_VERTEX\n  }\n");

  fprintf (fileptr, "  Coordinate3 {\n    point [\n\t%g\t%g\t%g", /* Write vertex coordinates */
           coortax[baseval].x,
           coortax[baseval].y,
           coortax[baseval].z);
  for (vertnum = baseval + 1; vertnum < vertnnd; vertnum ++)
    fprintf (fileptr, ",\n\t%g\t%g\t%g",
             coortax[vertnum].x,
             coortax[vertnum].y,
             coortax[vertnum].z);
  fprintf (fileptr, " ]\n  }\n");

  fprintf (fileptr, "  BaseColor {\n    rgb [");  /* Write color vector */
  for (indxnum = 0; indxnum < indxnbr - 2; indxnum ++) {
    if (indxtab[indxnum] != -1) {
      colfptr (parttax[indxtab[indxnum]], colotab);
      fprintf (fileptr, "\n\t%g\t%g\t%g,",
               (double) colotab[0],
               (double) colotab[1],
               (double) colotab[2]);
    }
  }
  colfptr (parttax[indxtab[indxnbr - 2]], colotab);
  fprintf (fileptr, "\n\t%g\t%g\t%g ]\n  }\n",
           (double) colotab[0],
           (double) colotab[1],
           (double) colotab[2]);

  fprintf (fileptr, "  IndexedLineSet {\n    coordIndex ["); /* Write set of lines */
  for (indxnum = 0; indxnum < indxnbr - 1; indxnum ++) {
    if ((indxnum % 8) == 0)
      fprintf (fileptr, "\n");
    fprintf (fileptr, "\t" SCOTCH_NUMSTRING ",", indxtab[indxnum] - baseval); /* Un-base vertex indices in paths */
  }
  if (((indxnbr - 1) % 8) == 0)
    fprintf (fileptr, "\n");
  fprintf (fileptr, "\t-1 ]\n  }\n");

  fprintf (fileptr, "}\n");                       /* Write end of separator */

  memFree (indxtab);                              /* Free group leader */

  return (0);
}

/*************************************************/
/*                                               */
/* This is the PostScript matrix output routine. */
/*                                               */
/*************************************************/

int
outDrawPosMatr (
const C_Graph * const       grafptr,              /* Graph structure, sorted by vertex index */
const C_Geometry * const    geomptr,              /* Graph geometry, sorted by vertex label  */
const C_Mapping * const     mappptr,              /* Result mapping, sorted by vertex label  */
FILE * const                fileptr)              /* Output stream                           */
{
  int *               nonztab;                    /* Array of non-zero entries               */
  SCOTCH_Num          nonzfrst;                   /* First non-zero entry of area            */
  SCOTCH_Num          nonzlast;                   /* Last non-zero entry of area             */
  double              pictcoorsiz;                /* Number of distinct coordinates          */
  double              pictdispsiz;                /* Size of the matrix display (in inches)  */
  time_t              timedat;                    /* Creation time                           */
  char *              timeptr;                    /* Pointer to string form of creation time */
  SCOTCH_Num          vertnum;

  const SCOTCH_Num                  baseval = grafptr->baseval;
  const SCOTCH_Num                  vertnnd = grafptr->vertnbr + grafptr->baseval;
  const SCOTCH_Num * restrict const verttax = grafptr->verttax;
  const SCOTCH_Num * restrict const vendtax = grafptr->vendtax;
  const SCOTCH_Num * restrict const edgetax = grafptr->edgetax;
  const SCOTCH_Num * restrict const labltax = mappptr->labltab - baseval;

  if ((nonztab = memAlloc ((grafptr->vertnbr + 1) * sizeof (int))) == NULL) {
    errorPrint ("outDrawPosMatr: out of memory");
    return (1);
  }

  time (&timedat);                                /* Get current time */
  timeptr = ctime (&timedat);
  pictcoorsiz = (double) (grafptr->vertnbr + 1);  /* Get matrix size */
  pictdispsiz = MIN (O_PSPICTWIDTH, O_PSPICTHEIGHT);

  if (O_outParam.PosMatr.typeval == 'e') {        /* EPSF-type output */
    fprintf (fileptr, "%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf (fileptr, "%%%%Title: %s %s %s\n",
             C_filenamesrcinp, C_filenamegeoinp, C_filenamemapinp);
    fprintf (fileptr, "%%%%Creator: Scotch/gout\n");
    fprintf (fileptr, "%%%%CreationDate: %s", timeptr);
    fprintf (fileptr, "%%%%BoundingBox: 0 0 %d %d\n",
             (int) (pictdispsiz * O_PSDPI), (int) (pictdispsiz * O_PSDPI));
    fprintf (fileptr, "%%%%Pages: 0\n");
    fprintf (fileptr, "%%%%EndComments\n");
  }
  else {                                          /* Full page output */
    fprintf (fileptr, "%%!PS-Adobe-2.0\n");
    fprintf (fileptr, "%%%%Title: %s %s %s\n",
             C_filenamesrcinp, C_filenamegeoinp, C_filenamemapinp);
    fprintf (fileptr, "%%%%Creator: Scotch/gout\n");
    fprintf (fileptr, "%%%%CreationDate: %s", timeptr);
  }

  fprintf (fileptr, "/p { pop } bind def\n");
  fprintf (fileptr, "/h { 3 1 roll exch 2 copy moveto 2 copy 1 add 5 -3 roll 3 1 roll add exch 2 copy lineto 1 add lineto lineto fill } bind def\n");
  fprintf (fileptr, "/v { 3 copy pop moveto 2 copy add exch pop exch 3 copy pop pop 1 add dup 3 -1 roll lineto exch dup 3 1 roll lineto lineto fill } bind def\n");
  fprintf (fileptr, "/b { 3 copy v 3 copy h pop pop } bind def\n");
  fprintf (fileptr, "/c { 1 3 copy v 3 copy h pop pop } bind def\n");

  fprintf (fileptr, "gsave\n");                   /* Save the context    */
  fprintf (fileptr, "0 setlinecap\n");            /* Perform miter caps  */
  if (O_outParam.PosMatr.typeval == 'f')          /* If full page output */
    fprintf (fileptr, "%d %d translate\n",        /* Center the picture  */
             (int) (O_PSDPI * (O_PSPAGEWIDTH - pictdispsiz)) / 2,
             (int) (O_PSDPI * (O_PSPAGEWIDTH - pictdispsiz)) / 2);
  fprintf (fileptr, "%f %f scale\n",              /* Print scaling factor */
           (double) O_PSDPI * pictdispsiz / pictcoorsiz,
           (double) O_PSDPI * pictdispsiz / pictcoorsiz);
  fprintf (fileptr, "[ 1 0 0 -1 0 %d ] concat\n", /* Reverse Y coordinate */
           (int) (grafptr->vertnbr + 1));
  fprintf (fileptr, "0 setgray newpath\n");       /* Select black color */

  for (vertnum = baseval; vertnum < vertnnd; vertnum ++) {
    SCOTCH_Num          colunum;
    SCOTCH_Num          edgennd;
    SCOTCH_Num          edgenum;

    colunum = (labltax[vertnum] == ~0) ? vertnum : labltax[vertnum];

    fprintf (fileptr, SCOTCH_NUMSTRING "\n",      /* Set column value */
             colunum - baseval);
    memset (nonztab, 0, (colunum - baseval + 2) * sizeof (int));
    for (edgenum = verttax[colunum], edgennd = vendtax[colunum]; edgenum < edgennd; edgenum ++) {
      SCOTCH_Num          coluend;

      coluend = edgetax[edgenum];
      if (coluend < colunum)
        nonztab[coluend - baseval] = 1;
    }
    nonztab[colunum - baseval] = 1;               /* Diagonal is non-zero */
    for (nonzfrst = 0; nonzfrst <= colunum - baseval; nonzfrst ++) {
      if (nonztab[nonzfrst] != 0) {               /* A non-zero has been found */
        for (nonzlast = nonzfrst; nonztab[nonzlast] != 0; nonzlast ++) ;
        if ((nonzlast - nonzfrst) > 1)            /* Draw row block coefficient */
          fprintf (fileptr, SCOTCH_NUMSTRING " " SCOTCH_NUMSTRING " b\n",
                   nonzfrst,
                   (nonzlast - nonzfrst));
        else
          fprintf (fileptr, SCOTCH_NUMSTRING " c\n",
                   nonzfrst);
        nonzfrst = nonzlast - 1;
      }
    }
    fprintf (fileptr, "p ");                      /* Close current column */
  }

  fprintf (fileptr, "\ngrestore\n");              /* Restore context     */
  if (O_outParam.PosMatr.typeval == 'f')          /* If full page output */
    fprintf (fileptr, "showpage\n");              /* Display the page    */

  memFree (nonztab);

  return (0);
}

/***********************************************/
/*                                             */
/* This is the PostScript mesh output routine. */
/*                                             */
/***********************************************/

int
outDrawPosMesh (
const C_Graph * const       grafptr,              /* Graph structure, sorted by vertex index */
const C_Geometry * const    geomptr,              /* Graph geometry, sorted by vertex label  */
const C_Mapping * const     mappptr,              /* Result mapping, sorted by vertex label  */
FILE * const                fileptr)              /* Output stream                           */
{
  SCOTCH_Num *        indxtab;                    /* Array of vertex indices                 */
  SCOTCH_Num          indxnbr;                    /* Number of indices                       */
  SCOTCH_Num          indxnum;
  O_PosMeshPath *     pathtax;                    /* Array of path building data             */
  O_PosMeshVertex *   veextax;                    /* Array of extended 2D data for vertices  */
  O_Point             pictcoormin;                /* Picture minimum and maximum coordinates */
  O_Point             pictcoormax;
  O_Point             pictcoordlt;
  double              pictscalval;                /* Scaling factor                          */
  time_t              timedat;                    /* Creation time                           */
  char *              timeptr;                    /* Pointer to string form of creation time */
  double              colotab[3];                 /* Color values                            */
  SCOTCH_Num          vertnum;

  const SCOTCH_Num                  baseval = grafptr->baseval;
  const SCOTCH_Num                  vertnnd = grafptr->vertnbr + grafptr->baseval;
  const SCOTCH_Num * restrict const verttax = grafptr->verttax;
  const SCOTCH_Num * restrict const vendtax = grafptr->vendtax;
  const SCOTCH_Num * restrict const edgetax = grafptr->edgetax;
  const C_GeoVert * restrict const  coortax = geomptr->coortab - baseval;
  const SCOTCH_Num * restrict const parttax = mappptr->labltab - baseval;

  if (geomptr->coortab == NULL) {
    errorPrint ("outDrawPosMesh: geometry not provided");
    return (1);
  }

  time (&timedat);                                /* Get current time */
  timeptr = ctime (&timedat);

  if (memAllocGroup ((void **) (void *)
                     &indxtab, (size_t) ((grafptr->edgenbr / 2) * 3 * sizeof (SCOTCH_Num)),
                     &pathtax, (size_t) (grafptr->vertnbr           * sizeof (O_PosMeshPath)),
                     &veextax, (size_t) (grafptr->vertnbr           * sizeof (O_PosMeshVertex)), NULL) == NULL) {
    errorPrint ("outDrawPosMesh: out of memory");
    return (1);
  }
  pathtax -= baseval;                             /* Base path and vertex arrays */
  veextax -= baseval;

  for (vertnum = baseval; vertnum < vertnnd; vertnum ++) { /* For all vertices               */
    veextax[vertnum].coorval.x = coortax[vertnum].x + /* Project 3D coordinates into 2D ones */
                                 coortax[vertnum].z * (O_POSMESHISOCOS * O_POSMESHISOREDUC);
    veextax[vertnum].coorval.y = coortax[vertnum].y +
                                 coortax[vertnum].z * (O_POSMESHISOSIN * O_POSMESHISOREDUC);
  }

  pictcoormin.x = pictcoormin.y =  1e30;          /* Pre-set coordinates extrema */
  pictcoormax.x = pictcoormax.y = -1e30;

  if (O_outParam.PosMesh.clipval == 'l') {        /* If clipping encompasses disks */
    for (vertnum = baseval; vertnum < vertnnd; vertnum ++) {
      double              sradmin;                /* Minimum square of circle radius */
      SCOTCH_Num          edgennd;
      SCOTCH_Num          edgenum;

      sradmin = 1e30;                             /* Assume a huge square of radius */
      for (edgenum = verttax[vertnum], edgennd = vendtax[vertnum];
           edgenum < edgennd; edgenum ++) {
        double              sradval;
        SCOTCH_Num          vertend;

        vertend = edgetax[edgenum];
        sradval = (veextax[vertnum].coorval.x - veextax[vertend].coorval.x) *
                  (veextax[vertnum].coorval.x - veextax[vertend].coorval.x) +
                  (veextax[vertnum].coorval.y - veextax[vertend].coorval.y) *
                  (veextax[vertnum].coorval.y - veextax[vertend].coorval.y);
        if (sradval < sradmin)                    /* Get the smallest square of radius */
          sradmin = sradval;
      }
      veextax[vertnum].dradval = sqrt (sradmin) * 0.5; /* Keep the half-distance for radius */

      if ((veextax[vertnum].coorval.x - veextax[vertnum].dradval) < pictcoormin.x) /* Update extrema if necessary */
        pictcoormin.x = veextax[vertnum].coorval.x - veextax[vertnum].dradval;
      if ((veextax[vertnum].coorval.y - veextax[vertnum].dradval) < pictcoormin.y)
        pictcoormin.y = veextax[vertnum].coorval.y - veextax[vertnum].dradval;
      if ((veextax[vertnum].coorval.x + veextax[vertnum].dradval) > pictcoormax.x)
        pictcoormax.x = veextax[vertnum].coorval.x + veextax[vertnum].dradval;
      if ((veextax[vertnum].coorval.y + veextax[vertnum].dradval) > pictcoormax.y)
        pictcoormax.y = veextax[vertnum].coorval.y + veextax[vertnum].dradval;
    }
  }
  else {                                          /* Border clipping                   */
    for (vertnum = baseval; vertnum < vertnnd; vertnum ++) { /* For all vertex indices */
      if (veextax[vertnum].coorval.x < pictcoormin.x) /* Update extrema if necessary   */
        pictcoormin.x = veextax[vertnum].coorval.x;
      if (veextax[vertnum].coorval.y < pictcoormin.y)
        pictcoormin.y = veextax[vertnum].coorval.y;
      if (veextax[vertnum].coorval.x > pictcoormax.x)
        pictcoormax.x = veextax[vertnum].coorval.x;
      if (veextax[vertnum].coorval.y > pictcoormax.y)
        pictcoormax.y = veextax[vertnum].coorval.y;
    }
  }

  pictcoordlt.x = pictcoormax.x - pictcoormin.x;  /* Compute picture extents */
  pictcoordlt.y = pictcoormax.y - pictcoormin.y;
  pictcoormin.x += pictcoordlt.x * O_outParam.PosMesh.pminval.x; /* Resize picture (if necessary) */
  pictcoormin.y += pictcoordlt.y * O_outParam.PosMesh.pminval.y;
  pictcoormax.x -= pictcoordlt.x * (1.0L - O_outParam.PosMesh.pmaxval.x);
  pictcoormax.y -= pictcoordlt.y * (1.0L - O_outParam.PosMesh.pmaxval.y);
  pictcoordlt.x = pictcoormax.x - pictcoormin.x;  /* Recompute picture extents */
  pictcoordlt.y = pictcoormax.y - pictcoormin.y;

  pictscalval = (pictcoordlt.x == 0.0L)           /* Compute scaling factor */
                ? ((pictcoordlt.y == 0.0L)
                   ? 1.0L
                   : (O_PSPICTHEIGHT / pictcoordlt.y))
                : ((pictcoordlt.y == 0.0L)
                   ? (O_PSPICTWIDTH / pictcoordlt.x)
                   : MIN (O_PSPICTWIDTH  / pictcoordlt.x, O_PSPICTHEIGHT / pictcoordlt.y));

  pictcoordlt.x *= pictscalval * O_POSMESHPICTRESOL; /* Rescale extents */
  pictcoordlt.y *= pictscalval * O_POSMESHPICTRESOL;

  for (vertnum = baseval; vertnum < vertnnd; vertnum ++) {
    veextax[vertnum].coorval.x = (veextax[vertnum].coorval.x - pictcoormin.x) * pictscalval * O_POSMESHPICTRESOL; /* Rescale coordinates */
    veextax[vertnum].coorval.y = (veextax[vertnum].coorval.y - pictcoormin.y) * pictscalval * O_POSMESHPICTRESOL;
  }

  if (O_outParam.PosMesh.diskval == 'd') {        /* If disks wanted */
    for (vertnum = baseval; vertnum < vertnnd; vertnum ++) {
      double              sradmin;                /* Minimum square of circle radius */
      SCOTCH_Num          edgennd;
      SCOTCH_Num          edgenum;

      sradmin = 1e30;                             /* Assume huge square of radius */
      for (edgenum = verttax[vertnum], edgennd = vendtax[vertnum]; edgenum < edgennd; edgenum ++) {
        SCOTCH_Num          vertend;
        double              sradval;

        vertend = edgetax[edgenum];
        sradval = (veextax[vertnum].coorval.x - veextax[vertend].coorval.x) *
                  (veextax[vertnum].coorval.x - veextax[vertend].coorval.x) +
                  (veextax[vertnum].coorval.y - veextax[vertend].coorval.y) *
                  (veextax[vertnum].coorval.y - veextax[vertend].coorval.y);
        if (sradval < sradmin)                    /* Get smallest square of radius */
          sradmin = sradval;
      }
      veextax[vertnum].dradval = sqrt (sradmin) * 0.5; /* Keep the half-distance for radius */
      if (veextax[vertnum].dradval < 1.0L)        /* Always get a non-zero radius           */
        veextax[vertnum].dradval = 1.0L;

      veextax[vertnum].visival =                  /* Compute vertex visibility */
        ((veextax[vertnum].coorval.x > - veextax[vertnum].dradval)               &&
         (veextax[vertnum].coorval.x < pictcoordlt.x + veextax[vertnum].dradval) &&
         (veextax[vertnum].coorval.y > - veextax[vertnum].dradval)               &&
         (veextax[vertnum].coorval.y < pictcoordlt.y + veextax[vertnum].dradval)) ? 1 : 0;
    }
  }
  else {                                          /* If disks not wanted */
    for (vertnum = baseval; vertnum < vertnnd; vertnum ++) {
      veextax[vertnum].visival =                  /* Compute vertex visibility */
        ((veextax[vertnum].coorval.x > 0.0L)          &&
         (veextax[vertnum].coorval.x < pictcoordlt.x) &&
         (veextax[vertnum].coorval.y > 0.0L)          &&
         (veextax[vertnum].coorval.y < pictcoordlt.y)) ? 1 : 0;
    }
  }

  for (vertnum = baseval; vertnum < vertnnd; vertnum ++) { /* Prepare to switch to integer coordinates */
    veextax[vertnum].coorval.x += 0.5L;
    veextax[vertnum].coorval.y += 0.5L;
  }
  pictcoordlt.x += 0.5L;
  pictcoordlt.y += 0.5L;

  if (O_outParam.PosMesh.coloval == 'c') {        /* If color output                        */
    for (vertnum = baseval; vertnum < vertnnd; vertnum ++) /* Select color for all vertices */
      veextax[vertnum].coloval = parttax[vertnum] % O_POSMESHCOLNBR;
  }
  else {                                          /* If gray level output                      */
    for (vertnum = baseval; vertnum < vertnnd; vertnum ++) { /* Select color for all vertices  */
      int                 parttmp;                /* Last 8 bits of mapping value              */
      int                 bitsnbr;                /* Number of bits to create color index from */
      int                 coloval;                /* Resulting color index from bit reversal   */

      for (parttmp = (int) (parttax[vertnum] & 255), bitsnbr = 7, coloval = 0; /* Half-tone color */
           bitsnbr > 0;                           /* As long as there are subdivision bits        */
           parttmp >>= 1, coloval <<= 1, bitsnbr --) /* Iterate on each bit                       */
        coloval |= (parttmp & 1);
      veextax[vertnum].coloval = coloval;
    }
  }

  for (vertnum = baseval; vertnum < vertnnd; vertnum ++) { /* For all vertices */
    SCOTCH_Num          pathnbr;
    SCOTCH_Num          partval;
    SCOTCH_Num          edgennd;
    SCOTCH_Num          edgenum;

    partval = parttax[vertnum];
    edgenum = verttax[vertnum];
    pathtax[vertnum].edgenum = edgenum;           /* Record start exit index */
    pathnbr = 0;                                  /* No output path yet      */
    for (edgennd = vendtax[vertnum]; edgenum < edgennd; edgenum ++) {
      SCOTCH_Num          vertend;

      vertend = edgetax[edgenum];
      if ((vertend > vertnum) &&                  /* If end vertex can be a next hop             */
          ((veextax[vertnum].visival | veextax[vertend].visival) != 0) && /* And edge is visible */
          ((O_outParam.PosMesh.edgeval != 'r') || /* And it can be drawn                         */
           (parttax[vertend] == partval)))
        pathnbr ++;                               /* One more path to vertices of higher indices */
    }
    pathtax[vertnum].pathnbr = pathnbr;           /* Record number of output paths */
  }

  indxnbr = 0;                                    /* No indices yet   */
  for (vertnum = baseval; vertnum < vertnnd; ) {  /* For all vertices */
    SCOTCH_Num          verttmp;

    if (pathtax[vertnum].pathnbr == 0) {          /* If no output path for this vertex */
      vertnum ++;                                 /* Skip to next vertex               */
      continue;
    }

    verttmp = vertnum;                            /* Begin with this vertex         */
    indxtab[indxnbr ++] = verttmp;                /* Add it to the current path     */
    do {                                          /* Build path from current vertex */
      SCOTCH_Num          parttmp;
      SCOTCH_Num          vertend;
      SCOTCH_Num          edgennd;
      SCOTCH_Num          edgetmp;

      parttmp = parttax[verttmp];
      edgetmp = pathtax[verttmp].edgenum;
      edgennd = vendtax[verttmp];
      do {                                        /* Search for first valid output vertex */
        vertend = edgetax[edgetmp];
        if ((vertend > verttmp) &&                /* If end vertex can be next hop                 */
            ((veextax[verttmp].visival | veextax[vertend].visival) != 0) && /* And edge is visible */
            ((O_outParam.PosMesh.edgeval != 'r') || /* And it can be drawn                         */
             (parttax[vertend] == parttmp)))
          break;
      } while (++ edgetmp < edgennd);

      pathtax[verttmp].pathnbr --;                /* One less output path remaining    */
      pathtax[verttmp].edgenum = edgetmp + 1;     /* Search from the next position     */
      verttmp = vertend;                          /* Go-on from end vertex             */
      indxtab[indxnbr ++] = verttmp;              /* Add end vertex to current path    */
    } while (pathtax[verttmp].pathnbr > 0);       /* As long as there is an exit route */

    indxtab[indxnbr ++] = -1;                     /* Mark end of current path */
  }

  if (O_outParam.PosMesh.typeval == 'e') {        /* EPSF-type output */
    fprintf (fileptr, "%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf (fileptr, "%%%%Title: %s %s %s\n",
             C_filenamesrcinp, C_filenamegeoinp, C_filenamemapinp);
    fprintf (fileptr, "%%%%Creator: Scotch/gout\n");
    fprintf (fileptr, "%%%%CreationDate: %s", timeptr);
    fprintf (fileptr, "%%%%BoundingBox: 0 0 %d %d\n",
             (int) ((pictcoordlt.x * O_PSDPI) / O_POSMESHPICTRESOL),
             (int) ((pictcoordlt.y * O_PSDPI) / O_POSMESHPICTRESOL));
    fprintf (fileptr, "%%%%Pages: 0\n");
    fprintf (fileptr, "%%%%EndComments\n");
  }
  else {                                          /* Full page output */
    fprintf (fileptr, "%%!PS-Adobe-2.0\n");
    fprintf (fileptr, "%%%%Title: %s %s %s\n",
             C_filenamesrcinp, C_filenamegeoinp, C_filenamemapinp);
    fprintf (fileptr, "%%%%Creator: Scotch/gout\n");
    fprintf (fileptr, "%%%%CreationDate: %s", timeptr);
  }

  fprintf (fileptr, "/A  { 0 360 arc fill } bind def\n"); /* Macro definitions */
  if (O_outParam.PosMesh.coloval == 'c') {        /* If color output           */
    int                 colonum;

    for (colonum = 0; colonum < O_POSMESHCOLNBR; colonum ++) { /* Build color indices */
      outColorColor (colonum, colotab);
      fprintf (fileptr, "/C%c { %g %g %g setrgbcolor } bind def\n",
               ('a' + colonum),
               colotab[0],
               colotab[1],
               colotab[2]);
    }
  }
  fprintf (fileptr, "/G  { 255 div setgray } bind def\n");
  fprintf (fileptr, "/L  { lineto stroke } bind def\n");
  fprintf (fileptr, "/l  { lineto } bind def\n");
  fprintf (fileptr, "/m  { moveto } bind def\n");

  fprintf (fileptr, "gsave\n");                   /* Save context        */
  fprintf (fileptr, "1 setlinecap\n");            /* Use round caps      */
  if (O_outParam.PosMesh.typeval == 'f')          /* If full page output */
    fprintf (fileptr, "%d %d translate\n",        /* Center picture      */
             (int) ((O_PSDPI * (O_PSPAGEWIDTH  * O_POSMESHPICTRESOL - pictcoordlt.x)) /
                    (2 * O_POSMESHPICTRESOL)),
             (int) ((O_PSDPI * (O_PSPAGEHEIGHT * O_POSMESHPICTRESOL - pictcoordlt.y)) /
                    (2 * O_POSMESHPICTRESOL)));
  fprintf (fileptr, "%f %f scale\n",              /* Print scaling factor */
           (double) O_PSDPI / O_POSMESHPICTRESOL,
           (double) O_PSDPI / O_POSMESHPICTRESOL);
  fprintf (fileptr, "newpath 0 0 m %d 0 l %d %d l 0 %d l closepath clip\n", /* Clip picture */
           (int) pictcoordlt.x, (int) pictcoordlt.x,
           (int) pictcoordlt.y, (int) pictcoordlt.y);

  fprintf (fileptr, "0 G\n");                     /* Select black color */
  for (indxnum = 0; indxnum < indxnbr; indxnum ++) {
    fprintf (fileptr, "%d\t%d\tm\n",              /* Set initial point */
             (int) veextax[indxtab[indxnum]].coorval.x,
             (int) veextax[indxtab[indxnum]].coorval.y);
    for (indxnum ++; indxtab[indxnum] != -1; indxnum ++) { /* Build path */
      fprintf (fileptr, "%d\t%d\t%c\n",
               (int) veextax[indxtab[indxnum]].coorval.x,
               (int) veextax[indxtab[indxnum]].coorval.y,
               (indxtab[indxnum + 1] == ~0) ? 'L' : 'l');
    }
  }

  if (O_outParam.PosMesh.diskval == 'd') {        /* If disks wanted */
    int                 coloval;

    coloval = -1;                                 /* No assigned color yet          */
    for (vertnum = baseval; vertnum < vertnnd; vertnum ++) { /* For all vertices    */
      if ((veextax[vertnum].visival != 0) &&      /* If disk is visible             */
          (parttax[vertnum] != -1)) {             /* And is mapped                  */
        if (veextax[vertnum].coloval != coloval) { /* If not same as previous color */
          coloval = veextax[vertnum].coloval;     /* Record the new current color   */

          if (O_outParam.PosMesh.coloval == 'c')  /* Update drawing color */
            fprintf (fileptr, "C%c\n", 'a' + coloval);
          else
            fprintf (fileptr, "%d G\n", coloval);
        }

        fprintf (fileptr, "%d %d %d A\n",         /* Draw the disk */
                (int) veextax[vertnum].coorval.x,
                (int) veextax[vertnum].coorval.y,
                (int) veextax[vertnum].dradval);
      }
    }
  }

  fprintf (fileptr, "grestore\n");                /* Restore the context */
  if (O_outParam.PosMesh.typeval == 'f')          /* If full page output */
    fprintf (fileptr, "showpage\n");              /* Display the page    */

  memFree (indxtab);                              /* Free group leader */

  return (0);
}

/*************************************/
/*                                   */
/* This is the Tulip output routine. */
/*                                   */
/*************************************/

int
outDrawTulMesh (
const C_Graph * const       grafptr,              /* Graph structure, sorted by vertex index */
const C_Geometry * const    geomptr,              /* Graph geometry, sorted by vertex label  */
const C_Mapping * const     mappptr,              /* Result mapping, sorted by vertex label  */
FILE * const                fileptr)              /* Output stream                           */
{
  time_t              timedat;                    /* Creation time */
  char                timestr[64];
  double              colotab[3];                 /* Vertex color  */
  SCOTCH_Num          vertnum;
  SCOTCH_Num          edgeidx;
  char                c;

  const SCOTCH_Num                  baseval = grafptr->baseval;
  const SCOTCH_Num                  vertnnd = grafptr->vertnbr + grafptr->baseval;
  const SCOTCH_Num * restrict const verttax = grafptr->verttax;
  const SCOTCH_Num * restrict const vendtax = grafptr->vendtax;
  const SCOTCH_Num * restrict const edgetax = grafptr->edgetax;
  const C_GeoVert * restrict const  coortax = geomptr->coortab - baseval;
  const SCOTCH_Num * restrict const parttax = mappptr->labltab - baseval;

  if (coortax == NULL) {
    errorPrint ("outDrawTulMesh: geometry not provided");
    return (1);
  }

  time (&timedat);                                /* Get current time */
  strncpy (timestr, ctime (&timedat), 63);
  timestr[63] = '\0';
  timestr[strlen (timestr) - 1] = '\0';

  fprintf (fileptr, "(tlp \"2.0\"\n(author \"Scotch/gout\")\n(date \"%s\")\n(comment \"%s %s %s\")\n", /* Write header */
           timestr,
           C_filenamesrcinp, C_filenamegeoinp, C_filenamemapinp);

  if (grafptr->vertnbr == 0) {                    /* If nothing to write */
    fprintf (fileptr, ")\n");
    return  (0);
  }

  fprintf (fileptr, "(nodes\n");                  /* Write node list */
  for (vertnum = baseval; vertnum < (vertnnd - 1); vertnum ++)
    fprintf (fileptr, SCOTCH_NUMSTRING "%c",
             vertnum,
             (((vertnum - baseval) & 7) == 7) ? '\n' : '\t');
  fprintf (fileptr, SCOTCH_NUMSTRING ")\n",
           vertnum);

  for (vertnum = baseval, edgeidx = 0; vertnum < vertnnd; vertnum ++) {
    SCOTCH_Num          edgennd;
    SCOTCH_Num          edgenum;

    for (edgenum = verttax[vertnum], edgennd = vendtax[vertnum]; edgenum < edgennd; edgenum ++) {
      SCOTCH_Num          vertend;

      vertend = edgetax[edgenum];
      if (vertend <= vertnum)
        continue;

      fprintf (fileptr, "(edge " SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING "\t" SCOTCH_NUMSTRING ")\n",
               edgeidx ++,
               vertnum,
               vertend);
    }
  }

  fprintf (fileptr, "(property 0 layout \"viewLayout\"\n"); /* Write node coordinates */
  c = '\n';
  for (vertnum = baseval; vertnum < vertnnd; vertnum ++) {
    if (vertnum == (vertnnd - 1))
      c = ')';
    fprintf (fileptr, "(node " SCOTCH_NUMSTRING "\t\"(%lf,%lf,%lf)\")%c",
             vertnum,
             (double) coortax[vertnum].x,
             (double) coortax[vertnum].y,
             (double) coortax[vertnum].z,
             c);
  }
  fprintf (fileptr, "\n");

  if (O_outParam.TulMesh.coloval == 'c') {
    fprintf (fileptr, "(property 0 color \"viewColor\"\n(default \"(255,255,255,255)\" \"(0,0,0,0)\")\n"); /* Write node color values */
    c = '\n';
    for (vertnum = baseval; vertnum < vertnnd; vertnum ++) {
      if (vertnum == (vertnnd - 1))
        c = ')';
      outColorColor (parttax[vertnum], colotab);
      fprintf (fileptr, "(node " SCOTCH_NUMSTRING " \"(%d,%d,%d,255)\")%c",
               vertnum,
               (int) (colotab[0] * 255.0),
               (int) (colotab[1] * 255.0),
               (int) (colotab[2] * 255.0),
               c);
    }
    fprintf (fileptr, "\n");
  }

  fprintf (fileptr, "(property 0 size \"viewSize\"\n(default \"(0,0,0)\" \"(0,0,0)\")"); /* Write default node size */
  if (O_outParam.TulMesh.diskval == 'd') {        /* If disks wanted */
    fprintf (fileptr, "\n");
    c = '\n';
    for (vertnum = baseval; vertnum < vertnnd; vertnum ++) {
      SCOTCH_Num          edgenum;
      SCOTCH_Num          edgennd;
      double              sradmin;                /* Minimum square value of disk radius */
      double              dradmin;

      if (vertnum == (vertnnd - 1))
        c = ')';

      sradmin = 1e30;                             /* Huge distance assumed */
      for (edgenum = verttax[vertnum], edgennd = vendtax[vertnum]; edgenum < edgennd; edgenum ++) {
        SCOTCH_Num          vertend;
        double              sradval;              /* Square value of disk radius */

        vertend = edgetax[edgenum];
        sradval = (coortax[vertend].x - coortax[vertnum].x) * (coortax[vertend].x - coortax[vertnum].x) +
                  (coortax[vertend].y - coortax[vertnum].y) * (coortax[vertend].y - coortax[vertnum].y) +
                  (coortax[vertend].z - coortax[vertnum].z) * (coortax[vertend].z - coortax[vertnum].z);
        if (sradval < sradmin)
          sradmin = sradval;
      }
      dradmin = sqrt (sradmin) * (0.5 * O_TULMESHDISKRATIO);
      fprintf (fileptr, "(node " SCOTCH_NUMSTRING " \"(%lf,%lf,%lf)\")%c",
               vertnum,
               dradmin, dradmin, dradmin, c);
    }
    fprintf (fileptr, "\n");
  }
  else
    fprintf (fileptr, ")\n");

  fprintf (fileptr, ")\n");

  return (0);
}

/*****************************************/
/*                                       */
/* This is the VTK output routine, based */
/* on the OpenInventor routine.          */
/*                                       */
/*****************************************/

int
outDrawVtkMesh (
const C_Graph * const       grafptr,              /* Graph structure, sorted by vertex index */
const C_Geometry * const    geomptr,              /* Graph geometry, sorted by vertex label  */
const C_Mapping * const     mappptr,              /* Result mapping, sorted by vertex label  */
FILE * const                fileptr)              /* Output stream                           */
{
  O_VtkMeshPath * restrict  pathtax;              /* Array of path building data */
  SCOTCH_Num                pathnbr;              /* Number of paths created     */
  SCOTCH_Num                pathnum;
  SCOTCH_Num * restrict     indxtab;              /* Array of indices            */
  SCOTCH_Num                indxnbr;              /* Number of indices           */
  SCOTCH_Num                indxnum;
  SCOTCH_Num                vertnum;
  time_t                    timedat;              /* Creation time               */

  const SCOTCH_Num                  baseval = grafptr->baseval;
  const SCOTCH_Num                  vertnnd = grafptr->vertnbr + grafptr->baseval;
  const SCOTCH_Num * restrict const verttax = grafptr->verttax;
  const SCOTCH_Num * restrict const vendtax = grafptr->vendtax;
  const SCOTCH_Num * restrict const edgetax = grafptr->edgetax;
  const C_GeoVert * restrict const  coortax = geomptr->coortab - baseval;
  const SCOTCH_Num * restrict const parttax = mappptr->labltab - baseval;

  if (geomptr->coortab == NULL) {
    errorPrint ("outDrawVtkMesh: geometry not provided");
    return (1);
  }

  time (&timedat);                                /* Get current time */

  if (memAllocGroup ((void **) (void *)
                     &indxtab, (size_t) ((grafptr->edgenbr / 2) * 3 * sizeof (SCOTCH_Num)),
                     &pathtax, (size_t) (grafptr->vertnbr           * sizeof (O_VtkMeshPath)), NULL) == NULL) {
    errorPrint ("outDrawVtkMesh: out of memory");
    return (1);
  }
  pathtax -= baseval;                             /* Base path array */

  for (vertnum = baseval; vertnum < vertnnd; vertnum ++) { /* For all vertices */
    SCOTCH_Num          pathnbr;
    SCOTCH_Num          partval;
    SCOTCH_Num          edgennd;
    SCOTCH_Num          edgenum;

    partval = parttax[vertnum];
    edgenum = verttax[vertnum];
    pathtax[vertnum].edgenum = edgenum;           /* Record start exit index */
    pathnbr = 0;                                  /* No output path yet      */
    for (edgennd = vendtax[vertnum]; edgenum < edgennd; edgenum ++) {
      SCOTCH_Num          vertend;

      vertend = edgetax[edgenum];
      if ((vertend > vertnum) &&                  /* If end vertex can be a next hop */
          ((O_outParam.VtkMesh.edgeval != 'r') || /* And this edge can be drawn      */
           (parttax[vertend] == partval)))
        pathnbr ++;                               /* One more path to vertices of higher indices */
    }
    pathtax[vertnum].pathnbr = pathnbr;           /* Record number of output paths */
  }

  indxnbr = 0;                                    /* No indices yet                        */
  pathnbr = 0;                                    /* No paths yet                          */
  for (vertnum = baseval; vertnum < vertnnd; ) {  /* For all vertices                      */
    SCOTCH_Num          indxtmp;                  /* Index where to place number of points */
    SCOTCH_Num          verttmp;

    if (pathtax[vertnum].pathnbr == 0) {          /* If no output path for this vertex */
      vertnum ++;                                 /* Skip to next vertex               */
      continue;
    }

    indxtmp = indxnbr ++;                         /* Save space for number of points */
    verttmp = vertnum;                            /* Begin with this vertex          */
    indxtab[indxnbr ++] = verttmp;                /* Add it to the current path      */
    do {                                          /* Build path from current vertex  */
      SCOTCH_Num          parttmp;
      SCOTCH_Num          vertend;
      SCOTCH_Num          edgennd;
      SCOTCH_Num          edgetmp;

      parttmp = parttax[verttmp];
      edgetmp = pathtax[verttmp].edgenum;
      edgennd = vendtax[verttmp];
      do {                                        /* Search for first valid output vertex */
        vertend = edgetax[edgetmp];
        if ((vertend > verttmp) &&                /* If end vertex can be next hop */
            ((O_outParam.VtkMesh.edgeval != 'r') || /* And this edge can be drawn  */
             (parttax[vertend] == parttmp)))
          break;
      } while (++ edgetmp < edgennd);

      pathtax[verttmp].pathnbr --;                /* One less output path remaining    */
      pathtax[verttmp].edgenum = edgetmp + 1;     /* Search from the next position     */
      verttmp = vertend;                          /* Go-on from end vertex             */
      indxtab[indxnbr ++] = verttmp;              /* Add end vertex to current path    */
    } while (pathtax[verttmp].pathnbr > 0);       /* As long as there is an exit route */
    indxtab[indxtmp] = indxnbr - indxtmp - 1;     /* Set size of created path          */
    pathnbr ++;                                   /* One more path created             */
  }

  fprintf (fileptr, "# vtk DataFile Version 2.0\n"); /* Write header */
  fprintf (fileptr, "%s %s %s | Created by Scotch/gout | %sASCII\n",
           C_filenamesrcinp, C_filenamegeoinp, C_filenamemapinp,
           ctime (&timedat));

  if (indxnbr == 0)                               /* If nothing to write */
    return (0);

  fprintf (fileptr, "\nDATASET UNSTRUCTURED_GRID\n");

  fprintf (fileptr, "\nPOINTS " SCOTCH_NUMSTRING " float\n",
           grafptr->vertnbr);
  for (vertnum = baseval; vertnum < vertnnd; vertnum ++)
    fprintf (fileptr, "%g\t%g\t%g\n",
             coortax[vertnum].x,
             coortax[vertnum].y,
             coortax[vertnum].z);

  fprintf (fileptr, "\nCELLS " SCOTCH_NUMSTRING " " SCOTCH_NUMSTRING "\n",
           pathnbr,
           indxnbr);
  for (indxnum = 0; indxnum < indxnbr; ) {
    SCOTCH_Num          pontnbr;                  /* Number of points in current segment */    
    SCOTCH_Num          indxtmp;

    pontnbr = indxtab[indxnum];
    for (indxtmp = indxnum + pontnbr; indxnum < indxtmp; )
      fprintf (fileptr, SCOTCH_NUMSTRING "\t", indxtab[indxnum ++]);
    fprintf (fileptr, SCOTCH_NUMSTRING "\n", indxtab[indxnum ++]);
  }

  fprintf (fileptr, "\nCELL_TYPES " SCOTCH_NUMSTRING "\n", /* All cells ate poly-line paths */
           pathnbr);
  for (pathnum = 0; pathnum < pathnbr; pathnum ++)
    fprintf (fileptr, "4\n");

  fprintf (fileptr, "\nPOINT_DATA " SCOTCH_NUMSTRING "\nSCALARS mapValues int\nLOOKUP_TABLE default\n",
           grafptr->vertnbr);
  for (vertnum = baseval; vertnum < vertnnd; vertnum ++) {
    SCOTCH_Num          partval;

    partval = parttax[vertnum];
    fprintf (fileptr, SCOTCH_NUMSTRING "\n",
             (partval == -1) ? 0 : (partval + 1)); /* Unmapped vertices are of VTK color 0 (white) */
  }

  memFree (indxtab);                              /* Free group leader */

  return (0);
}
