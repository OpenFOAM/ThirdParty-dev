/* Copyright 2004,2007,2009-2016,2018,2021,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_arch.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the target   **/
/**                architecture handling routines of the   **/
/**                libSCOTCH library.                      **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 18 aug 1998     **/
/**                                 to   : 18 aug 1998     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to   : 29 mar 1999     **/
/**                # Version 3.4  : from : 01 nov 2001     **/
/**                                 to   : 01 nov 2001     **/
/**                # Version 4.0  : from : 13 jan 2004     **/
/**                                 to   : 13 jan 2004     **/
/**                # Version 5.0  : from : 12 sep 2007     **/
/**                                 to   : 12 sep 2007     **/
/**                # Version 5.1  : from : 05 jun 2009     **/
/**                                 to   : 13 feb 2011     **/
/**                # Version 6.0  : from : 14 feb 2011     **/
/**                                 to   : 05 jun 2018     **/
/**                # Version 6.1  : from : 15 mar 2021     **/
/**                                 to   : 15 mar 2021     **/
/**                # Version 7.0  : from : 21 jan 2023     **/
/**                                 to   : 21 jan 2023     **/
/**                # Version 8.0  : from : 18 feb 2018     **/
/**                                 to   : 18 feb 2018     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "arch_cmplt.h"
#include "arch_cmpltw.h"
#include "arch_hcub.h"
#include "arch_mesh.h"
#include "arch_sub.h"
#include "arch_tleaf.h"
#include "arch_torus.h"
#include "arch_vcmplt.h"
#include "arch_vhcub.h"
#include "scotch.h"

/***************************************/
/*                                     */
/* These routines are the C API for    */
/* the architecture handling routines. */
/*                                     */
/***************************************/

/*+ This routine reserves a memory area
*** of a size sufficient to store a
*** SCOTCH_ArchDom structure.
*** It returns:
*** - !NULL  : if the allocation succeeded.
*** - NULL   : on error.
+*/

SCOTCH_ArchDom *
SCOTCH_archDomAlloc ()
{
  return ((SCOTCH_ArchDom *) memAlloc (sizeof (SCOTCH_ArchDom)));
}

/*+ This routine returns the size, in bytes,
*** of a SCOTCH_ArchDom structure.
*** It returns:
*** - > 0  : in all cases.
+*/

int
SCOTCH_archDomSizeof ()
{
  return (sizeof (SCOTCH_ArchDom));
}

/*+ This routine reserves a memory area
*** of a size sufficient to store a
*** SCOTCH_Arch structure.
*** It returns:
*** - !NULL  : if the allocation succeeded.
*** - NULL   : on error.
+*/

SCOTCH_Arch *
SCOTCH_archAlloc ()
{
  return ((SCOTCH_Arch *) memAlloc (sizeof (SCOTCH_Arch)));
}

/*+ This routine returns the size, in bytes,
*** of a SCOTCH_Arch structure.
*** It returns:
*** - > 0  : in all cases.
+*/

int
SCOTCH_archSizeof ()
{
  return (sizeof (SCOTCH_Arch));
}

/*+ This routine initializes the opaque
*** architecture structure used to handle
*** target architectures in the Scotch library.
*** It returns:
*** - 0   : if the initialization succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_archInit (
SCOTCH_Arch * const         archptr)
{
  if (sizeof (SCOTCH_Num) != sizeof (Anum)) {
    errorPrint (STRINGIFY (SCOTCH_archInit) ": internal error (1)");
    return (1);
  }
  if (sizeof (SCOTCH_Arch) < sizeof (Arch)) {
    errorPrint (STRINGIFY (SCOTCH_archInit) ": internal error (2)");
    return (1);
  }

  return (archInit ((Arch *) archptr));
}

/*+ This routine frees the contents of the
*** given opaque architecture structure.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_archExit (
SCOTCH_Arch * const         archptr)
{
  archExit ((Arch *) archptr);
}

/*+ This routine loads the given opaque
*** architecture structure with the data of
*** the given stream.
*** It returns:
*** - 0   : if the loading succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_archLoad (
SCOTCH_Arch * const         archptr,
FILE * const                stream)
{
  return (archLoad ((Arch *) archptr, stream));
}

/*+ This routine saves the given opaque
*** architecture structure to the given
*** stream.
*** It returns:
*** - 0   : if the saving succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_archSave (
const SCOTCH_Arch * const   archptr,
FILE * const                stream)
{
  return (archSave ((Arch *) archptr, stream));
}

/*+ This routine returns the name of the
*** given target architecture.
*** It returns:
*** - !NULL  : pointer to the name of the
***            target architecture.
+*/

char *
SCOTCH_archName (
const SCOTCH_Arch * const   archptr)
{
  return (archName ((const Arch * const) archptr));
}

/*+ This routine returns the size of the
*** given target architecture.
*** It returns:
*** - !0  : size of the target architecture.
+*/

SCOTCH_Num
SCOTCH_archSize (
const SCOTCH_Arch * const   archptr)
{
  ArchDom             domdat;

  archDomFrst ((Arch *) archptr, &domdat);        /* Get first domain     */
  return (archDomSize ((Arch *) archptr, &domdat)); /* Return domain size */
}

/*+ This routine tells if the given architecture
*** is a variable-sized architecture or not.
*** It returns:
*** - 0  : if the architecture is not variable-sized.
*** - 1  : if the architecture is variable-sized.
+*/

int
SCOTCH_archVar (
const SCOTCH_Arch * const   archptr)
{
  return ((archVar ((Arch *) archptr) != 0) ? 1 : 0);
}

/*+ These routines fill the contents of the given
*** opaque target structure so as to yield target
*** architectures of the given types.
*** It returns:
*** - 0   : if the computation succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_archCmplt (
SCOTCH_Arch * const         archptr,
const SCOTCH_Num            termnbr)
{
  Arch *              tgtarchptr;
  ArchCmplt *         tgtarchdatptr;

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_archCmplt) ": internal error");
    return (1);
  }

  tgtarchptr    = (Arch *) archptr;
  tgtarchdatptr = (ArchCmplt *) (void *) (&tgtarchptr->data);

  tgtarchptr->clasptr    = archClass ("cmplt");
  tgtarchptr->flagval    = tgtarchptr->clasptr->flagval; /* Copy architecture flag */
  tgtarchdatptr->termnbr = (Anum) termnbr;

  return (0);
}

/*
**
*/

int
SCOTCH_archCmpltw (
SCOTCH_Arch * const         archptr,
const SCOTCH_Num            vertnbr,
const SCOTCH_Num * const    velotab)
{
  Arch *              tgtarchptr;

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_archCmpltw) ": internal error");
    return (1);
  }

  tgtarchptr          = (Arch *) archptr;
  tgtarchptr->clasptr = archClass ("cmpltw");
  tgtarchptr->flagval = tgtarchptr->clasptr->flagval; /* Copy architecture flag */

  return (archCmpltwArchBuild ((ArchCmpltw *) (void *) (&tgtarchptr->data), vertnbr, velotab));
}

/*
**
*/

int
SCOTCH_archHcub (
SCOTCH_Arch * const         archptr,
const SCOTCH_Num            dimnnbr)              /*+ Number of dimensions +*/
{
  Arch *              tgtarchptr;
  ArchHcub *          tgtarchdatptr;

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_archHcub) ": internal error");
    return (1);
  }

  tgtarchptr    = (Arch *) archptr;
  tgtarchdatptr = (ArchHcub *) (void *) (&tgtarchptr->data);

  tgtarchptr->clasptr    = archClass ("hcub");
  tgtarchptr->flagval    = tgtarchptr->clasptr->flagval; /* Copy architecture flag */
  tgtarchdatptr->dimnnbr = (Anum) dimnnbr;

  return (0);
}

/*
**
*/

int
SCOTCH_archMesh2 (
SCOTCH_Arch * const         archptr,
const SCOTCH_Num            dimxval,
const SCOTCH_Num            dimyval)
{
  Arch *              tgtarchptr;
  ArchMesh2 *         tgtarchdatptr;

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_archMesh2) ": internal error");
    return (1);
  }

  tgtarchptr    = (Arch *) archptr;
  tgtarchdatptr = (ArchMesh2 *) (void *) (&tgtarchptr->data);

  tgtarchptr->clasptr = archClass ("mesh2D");
  tgtarchptr->flagval = tgtarchptr->clasptr->flagval; /* Copy architecture flag */
  tgtarchdatptr->dimnnbr = 2;
  tgtarchdatptr->c[0]    = (Anum) dimxval;
  tgtarchdatptr->c[1]    = (Anum) dimyval;

  return (0);
}

/*
**
*/

int
SCOTCH_archMesh3 (
SCOTCH_Arch * const         archptr,
const SCOTCH_Num            dimxval,
const SCOTCH_Num            dimyval,
const SCOTCH_Num            dimzval)
{
  Arch *              tgtarchptr;
  ArchMesh3 *         tgtarchdatptr;

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_archMesh3) ": internal error");
    return (1);
  }

  tgtarchptr    = (Arch *) archptr;
  tgtarchdatptr = (ArchMesh3 *) (void *) (&tgtarchptr->data);

  tgtarchptr->clasptr = archClass ("mesh3D");
  tgtarchptr->flagval = tgtarchptr->clasptr->flagval; /* Copy architecture flag */
  tgtarchdatptr->dimnnbr = 3;
  tgtarchdatptr->c[0]    = (Anum) dimxval;
  tgtarchdatptr->c[1]    = (Anum) dimyval;
  tgtarchdatptr->c[2]    = (Anum) dimzval;

  return (0);
}

/*
**
*/

int
SCOTCH_archMeshX (
SCOTCH_Arch * const         archptr,
const SCOTCH_Num            dimnnbr,              /*+ Number of dimensions in architecture +*/
const SCOTCH_Num * const    dimntab)              /*+ Array of dimensions                  +*/
{
  Arch *              tgtarchptr;
  ArchMeshX *         tgtarchdatptr;

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_archMeshX) ": internal error");
    return (1);
  }
  if (dimnnbr > ARCHMESHDIMNMAX) {
    errorPrint (STRINGIFY (SCOTCH_archMeshX) ": too many dimensions");
    return (1);
  }

  tgtarchptr    = (Arch *) archptr;
  tgtarchdatptr = (ArchMeshX *) (void *) (&tgtarchptr->data);

  tgtarchptr->clasptr    = archClass ("meshXD");
  tgtarchptr->flagval    = tgtarchptr->clasptr->flagval; /* Copy architecture flag */
  tgtarchdatptr->dimnnbr = dimnnbr;
  memCpy (tgtarchdatptr->c, dimntab, dimnnbr * sizeof (SCOTCH_Num)); /* Copy dimension array */

  return (0);
}

/*
**
*/

int
SCOTCH_archSub (
SCOTCH_Arch * const         subarchptr,
SCOTCH_Arch * const         orgarchptr,
const SCOTCH_Num            vnumnbr,              /*+ Number of dimensions +*/
const SCOTCH_Num * const    vnumtab)
{
  Arch *              tgtarchptr;
  ArchSub *           tgtarchdatptr;

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_archSub) ": internal error");
    return (2);
  }

  tgtarchptr    = (Arch *) subarchptr;
  tgtarchdatptr = (ArchSub *) (void *) (&tgtarchptr->data);

  tgtarchptr->clasptr = archClass ("sub");
  tgtarchptr->flagval = tgtarchptr->clasptr->flagval; /* Copy architecture flag */

  return (archSubArchBuild (tgtarchdatptr, (Arch *) orgarchptr, (Gnum) vnumnbr, (Gnum *) vnumtab));
}

/*
**
*/

int
SCOTCH_archTleaf (
SCOTCH_Arch * const         archptr,
const SCOTCH_Num            levlnbr,              /*+ Number of levels in architecture            +*/
const SCOTCH_Num * const    sizetab,              /*+ Size array, by increasing level number      +*/
const SCOTCH_Num * const    linktab)              /*+ Link cost array, by increasing level number +*/
{
  Anum                levlnum;
  Anum                sizeval;
  Arch *              tgtarchptr;
  ArchTleaf *         tgtarchdatptr;

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_archTleaf) ": internal error");
    return (1);
  }

  tgtarchptr          = (Arch *) archptr;
  tgtarchdatptr       = (ArchTleaf *) (void *) (&tgtarchptr->data);
  tgtarchptr->clasptr = archClass ("tleaf");
  tgtarchptr->flagval = tgtarchptr->clasptr->flagval; /* Copy architecture flag */

  if ((tgtarchdatptr->sizetab = memAlloc ((levlnbr * 2 + 1) * sizeof (Anum))) == NULL) { /* TRICK: One more slot for linktab[-1] */
    errorPrint (STRINGIFY (SCOTCH_archTleaf) ": out of memory");
    return (1);
  }
  tgtarchdatptr->levlnbr     = (Anum) levlnbr;
  tgtarchdatptr->linktab     = tgtarchdatptr->sizetab + tgtarchdatptr->levlnbr + 1;
  tgtarchdatptr->linktab[-1] = 0;                 /* TRICK: Dummy slot for for level-0 communication */
  tgtarchdatptr->permtab     = NULL;              /* Assume no permutation array                     */

  for (levlnum = 0, sizeval = 1; levlnum < tgtarchdatptr->levlnbr; levlnum ++) {
    tgtarchdatptr->sizetab[levlnum] = sizetab[levlnum];
    tgtarchdatptr->linktab[levlnum] = linktab[levlnum];
    sizeval *= tgtarchdatptr->sizetab[levlnum];
  }
  tgtarchdatptr->termnbr = sizeval;

  return (0);
}

/*
**
*/

int
SCOTCH_archLtleaf (
SCOTCH_Arch * const         archptr,
const SCOTCH_Num            levlnbr,              /*+ Number of levels in architecture            +*/
const SCOTCH_Num * const    sizetab,              /*+ Size array, by increasing level number      +*/
const SCOTCH_Num * const    linktab,              /*+ Link cost array, by increasing level number +*/
const SCOTCH_Num            permnbr,              /*+ Number of permutation indices               +*/
const SCOTCH_Num * const    permtab)              /*+ Permutation array                           +*/
{
  Anum                permnum;
  Arch *              tgtarchptr;
  ArchTleaf *         tgtarchdatptr;

  if (SCOTCH_archTleaf (archptr, levlnbr, sizetab, linktab) != 0) /* Build tree part */
    return (1);
  tgtarchptr          = (Arch *) archptr;
  tgtarchdatptr       = (ArchTleaf *) (void *) (&tgtarchptr->data);
  tgtarchptr->clasptr = archClass ("ltleaf");     /* Override class */

  if ((tgtarchdatptr->permtab = memAlloc (permnbr * 2 * sizeof (Anum))) == NULL) { /* TRICK: space for peritab too */
    errorPrint (STRINGIFY (SCOTCH_archLtleaf) ": out of memory");
    return (1);
  }
  tgtarchdatptr->permnbr = (Anum) permnbr;
  tgtarchdatptr->peritab = tgtarchdatptr->permtab + permnbr;

  for (permnum = 0; permnum < tgtarchdatptr->permnbr; permnum ++)
    tgtarchdatptr->permtab[permnum] = permtab[permnum];
  for (permnum = 0; permnum < tgtarchdatptr->permnbr; permnum ++) /* Build inverse permutation */
    tgtarchdatptr->peritab[tgtarchdatptr->permtab[permnum]] = permnum;

  return (0);
}

/*
**
*/

int
SCOTCH_archTorus2 (
SCOTCH_Arch * const         archptr,
const SCOTCH_Num            dimxval,
const SCOTCH_Num            dimyval)
{
  Arch *              tgtarchptr;
  ArchTorusX *        tgtarchdatptr;

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_archTorus2) ": internal error");
    return (1);
  }

  tgtarchptr    = (Arch *) archptr;
  tgtarchdatptr = (ArchTorusX *) (void *) (&tgtarchptr->data);

  tgtarchptr->clasptr    = archClass ("torus2D");
  tgtarchptr->flagval    = tgtarchptr->clasptr->flagval; /* Copy architecture flag */
  tgtarchdatptr->dimnnbr = 2;
  tgtarchdatptr->c[0]    = (Anum) dimxval;
  tgtarchdatptr->c[1]    = (Anum) dimyval;

  return (0);
}

/*
**
*/

int
SCOTCH_archTorus3 (
SCOTCH_Arch * const         archptr,
const SCOTCH_Num            dimxval,
const SCOTCH_Num            dimyval,
const SCOTCH_Num            dimzval)
{
  Arch *              tgtarchptr;
  ArchTorusX *        tgtarchdatptr;

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_archTorus3) ": internal error");
    return (1);
  }

  tgtarchptr    = (Arch *) archptr;
  tgtarchdatptr = (ArchTorusX *) (void *) (&tgtarchptr->data);

  tgtarchptr->clasptr    = archClass ("torus3D");
  tgtarchptr->flagval    = tgtarchptr->clasptr->flagval; /* Copy architecture flag */
  tgtarchdatptr->dimnnbr = 3;
  tgtarchdatptr->c[0]    = (Anum) dimxval;
  tgtarchdatptr->c[1]    = (Anum) dimyval;
  tgtarchdatptr->c[2]    = (Anum) dimzval;

  return (0);
}

/*
**
*/

int
SCOTCH_archTorusX (
SCOTCH_Arch * const         archptr,
const SCOTCH_Num            dimnnbr,              /*+ Number of dimensions in architecture +*/
const SCOTCH_Num * const    dimntab)              /*+ Array of dimensions                  +*/
{
  Arch *              tgtarchptr;
  ArchTorusX *        tgtarchdatptr;

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_archTorusX) ": internal error");
    return (1);
  }
  if (dimnnbr > ARCHMESHDIMNMAX) {
    errorPrint (STRINGIFY (SCOTCH_archTorusX) ": too many dimensions");
    return (1);
  }

  tgtarchptr    = (Arch *) archptr;
  tgtarchdatptr = (ArchTorusX *) (void *) (&tgtarchptr->data);

  tgtarchptr->clasptr    = archClass ("torusXD");
  tgtarchptr->flagval    = tgtarchptr->clasptr->flagval; /* Copy architecture flag */
  tgtarchdatptr->dimnnbr = dimnnbr;
  memCpy (tgtarchdatptr->c, dimntab, dimnnbr * sizeof (SCOTCH_Num)); /* Copy dimension array */

  return (0);
}

/*
**
*/

int
SCOTCH_archVcmplt (
SCOTCH_Arch * const         archptr)
{
  Arch *              tgtarchptr;

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_archVcmplt) ": internal error");
    return (1);
  }

  tgtarchptr = (Arch *) archptr;

  tgtarchptr->clasptr = archClass ("varcmplt");
  tgtarchptr->flagval = tgtarchptr->clasptr->flagval; /* Copy architecture flag */

  return (0);
}

/*
**
*/

int
SCOTCH_archVhcub (
SCOTCH_Arch * const         archptr)
{
  Arch *              tgtarchptr;

  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_archVhcub) ": internal error");
    return (1);
  }

  tgtarchptr = (Arch *) archptr;

  tgtarchptr->clasptr = archClass ("varhcub");
  tgtarchptr->flagval = tgtarchptr->clasptr->flagval; /* Copy architecture flag */

  return (0);
}
