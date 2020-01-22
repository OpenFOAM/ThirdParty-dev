/* Copyright 2009,2018 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : module.h                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This is the global configuration file   **/
/**                for the ESMUMPS library module.         **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 22 jan 2009     **/
/**                                 to     22 jan 2009     **/
/**                # Version 6.0  : from : 21 may 2018     **/
/**                                 to     21 may 2018     **/
/**                                                        **/
/************************************************************/

#define MODULE_H

/*
** Function renaming.
*/

#if ((! defined SCOTCH_COMMON_EXTERNAL) || (defined SCOTCH_COMMON_RENAME))
#define clockGet                    _SCOTCHclockGet

#define fileNameDistExpand          _SCOTCHfileNameDistExpand 

#define usagePrint                  _SCOTCHusagePrint

#define errorPrint                  SCOTCH_errorPrint
#define errorPrintW                 SCOTCH_errorPrintW
#define errorProg                   SCOTCH_errorProg

#define intLoad                     _SCOTCHintLoad
#define intSave                     _SCOTCHintSave
#define intAscn                     _SCOTCHintAscn
#define intPerm                     _SCOTCHintPerm
#define intRandReset                _SCOTCHintRandReset
#define intRandInit                 _SCOTCHintRandInit
/* #define intRandVal               _SCOTCHintRandVal Already a macro */
#define intSearchDicho              _SCOTCHintSearchDicho
#define intSort1asc1                _SCOTCHintSort1asc1
#define intSort2asc1                _SCOTCHintSort2asc1
#define intSort2asc2                _SCOTCHintSort2asc2
#define intSort3asc1                _SCOTCHintSort3asc1

#define memAllocGroup               _SCOTCHmemAllocGroup
#define memReallocGroup             _SCOTCHmemReallocGroup
#define memOffset                   _SCOTCHmemOffset
#endif /* ((! defined SCOTCH_COMMON_EXTERNAL) || (defined SCOTCH_COMMON_RENAME)) */

#ifndef ESMUMPS_NAME_PREFIX_INTERN
#define ESMUMPS_NAME_PREFIX_INTERN  _ESMUMPS
#define ESMUMPS_NAME_PREFIX_PUBLICFL esmumpsf
#define ESMUMPS_NAME_PREFIX_PUBLICFU ESMUMPSF
#endif /* ESMUMPS_NAME_PREFIX_INTERN */

#ifndef ESMUMPS_NAME_SUFFIX
#define ESMUMPS_NAME_SUFFIX
#endif /* ESMUMPS_NAME_SUFFIX */
#ifndef ESMUMPS_NAME_SUFFIXFL
#define ESMUMPS_NAME_SUFFIXFL       ESMUMPS_NAME_SUFFIX
#define ESMUMPS_NAME_SUFFIXFU       ESMUMPS_NAME_SUFFIX
#endif /* SCOTCH_NAME_SUFFIXFL */

#define ESMUMPS_VOID

#define ESMUMPS_NAME_GLUE2(n,s)     n##s
#define ESMUMPS_NAME_GLUE3(p,n,s)   p##n##s
#define ESMUMPS_NAME_MACRO2(n,s)    ESMUMPS_NAME_GLUE2 (n,s)
#define ESMUMPS_NAME_MACRO3(p,n,s)  ESMUMPS_NAME_GLUE3 (p,n,s)
#define ESMUMPS_NAME_INTERN(f)      ESMUMPS_NAME_MACRO3 (ESMUMPS_NAME_PREFIX_INTERN,f,ESMUMPS_NAME_SUFFIX)
#define ESMUMPS_NAME_PUBLIC(f)      ESMUMPS_NAME_MACRO2 (f,ESMUMPS_NAME_SUFFIX)
#define ESMUMPS_NAME_PUBLICFL(f)    ESMUMPS_NAME_MACRO3 (ESMUMPS_NAME_PREFIX_PUBLICFL,f,ESMUMPS_NAME_SUFFIXFL)
#define ESMUMPS_NAME_PUBLICFU(f)    ESMUMPS_NAME_MACRO3 (ESMUMPS_NAME_PREFIX_PUBLICFU,f,ESMUMPS_NAME_SUFFIXFU)
#define ESMUMPS_FORTRAN(nu,nl,pl,pc) FORTRAN (ESMUMPS_NAME_PUBLICFU(nu),ESMUMPS_NAME_PUBLICFL(nl),pl,pc)

#define dofConstant                 ESMUMPS_NAME_INTERN (dofConstant)
#define dofGraph                    ESMUMPS_NAME_INTERN (dofGraph)
#define dofExit                     ESMUMPS_NAME_INTERN (dofExit)
#define dofInit                     ESMUMPS_NAME_INTERN (dofInit)
#define dofLoad                     ESMUMPS_NAME_INTERN (dofLoad)
#define dofSave                     ESMUMPS_NAME_INTERN (dofSave)

#define graphBuild                  ESMUMPS_NAME_INTERN (graphBuild)
#define graphBuildGraph             ESMUMPS_NAME_INTERN (graphBuildGraph)
#define graphBuildGraph2            ESMUMPS_NAME_INTERN (graphBuildGraph2)

#define orderBase                   ESMUMPS_NAME_INTERN (orderBase)
#define orderCheck                  ESMUMPS_NAME_INTERN (orderCheck)
#define orderExit                   ESMUMPS_NAME_INTERN (orderExit)
#define orderGraph                  ESMUMPS_NAME_INTERN (orderGraph)
#define orderGraphList              ESMUMPS_NAME_INTERN (orderGraphList)
#define orderGraphListStrat         ESMUMPS_NAME_INTERN (orderGraphListStrat)
#define orderGraphStrat             ESMUMPS_NAME_INTERN (orderGraphStrat)
#define orderGrid2                  ESMUMPS_NAME_INTERN (orderGrid2)
#define orderGrid2C                 ESMUMPS_NAME_INTERN (orderGrid2C)
#define orderGrid3                  ESMUMPS_NAME_INTERN (orderGrid3)
#define orderGrid3C                 ESMUMPS_NAME_INTERN (orderGrid3C)
#define orderInit                   ESMUMPS_NAME_INTERN (orderInit)
#define orderLoad                   ESMUMPS_NAME_INTERN (orderLoad)
#define orderMesh                   ESMUMPS_NAME_INTERN (orderMesh)
#define orderMeshList               ESMUMPS_NAME_INTERN (orderMeshList)
#define orderMeshListStrat          ESMUMPS_NAME_INTERN (orderMeshListStrat)
#define orderMeshStrat              ESMUMPS_NAME_INTERN (orderMeshStrat)
#define orderSave                   ESMUMPS_NAME_INTERN (orderSave)

#define symbolCheck                 ESMUMPS_NAME_INTERN (symbolCheck)
#define symbolCost                  ESMUMPS_NAME_INTERN (symbolCost)
#define symbolCosti                 ESMUMPS_NAME_INTERN (symbolCosti)
#define symbolDraw                  ESMUMPS_NAME_INTERN (symbolDraw)
#define symbolDrawColor             ESMUMPS_NAME_INTERN (symbolDrawColor)
#define symbolDrawFunc              ESMUMPS_NAME_INTERN (symbolDrawFunc)
#define symbolExit                  ESMUMPS_NAME_INTERN (symbolExit)
#define symbolInit                  ESMUMPS_NAME_INTERN (symbolInit)
#define symbolLevf                  ESMUMPS_NAME_INTERN (symbolLevf)
#define symbolLoad                  ESMUMPS_NAME_INTERN (symbolLoad)
#define symbolNonzeros              ESMUMPS_NAME_INTERN (symbolNonzeros)
#define symbolRealloc               ESMUMPS_NAME_INTERN (symbolRealloc)
#define symbolSave                  ESMUMPS_NAME_INTERN (symbolSave)
#define symbolTree                  ESMUMPS_NAME_INTERN (symbolTree)

#define symbolKeepAdd               ESMUMPS_NAME_INTERN (symbolKeepAdd)
#define symbolKeepCompute           ESMUMPS_NAME_INTERN (symbolKeepCompute)
#define symbolKeepDel               ESMUMPS_NAME_INTERN (symbolKeepDel)
#define symbolKeepExit              ESMUMPS_NAME_INTERN (symbolKeepExit)
#define symbolKeepHisto             ESMUMPS_NAME_INTERN (symbolKeepHisto)
#define symbolKeepInit              ESMUMPS_NAME_INTERN (symbolKeepInit)
#define symbolKeepPurge             ESMUMPS_NAME_INTERN (symbolKeepPurge)
#define symbolKeepView              ESMUMPS_NAME_INTERN (symbolKeepView)
