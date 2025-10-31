/* Copyright 2020 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : esmumps.h                               **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the symbolic factorization routine. **/
/**                                                        **/
/**   DATES      : # Version 6.1  : from : 05 sep 2020     **/
/**                                 to   : 05 sep 2020     **/
/**                # Version 7.0  : from : 26 aug 2024     **/
/**                                 to   : 26 aug 2024     **/
/**                                                        **/
/**   NOTES      : # This code derives from that of the    **/
/**                  internal "esmumps.h" header file.     **/
/**                                                        **/
/************************************************************/

#ifndef LIB_ESMUMPS_H
#define LIB_ESMUMPS_H

/*
**  The function prototypes.
*/

int                         esmumps             (const SCOTCH_Num n, const SCOTCH_Num iwlen, SCOTCH_Num * const pe, const SCOTCH_Num pfree, SCOTCH_Num * const len, SCOTCH_Num * const iw, SCOTCH_Num * const nv, SCOTCH_Num * const elen, SCOTCH_Num * const last);
int                         esmumpsv            (const SCOTCH_Num n, const SCOTCH_Num iwlen, SCOTCH_Num * const pe, const SCOTCH_Num pfree, SCOTCH_Num * const len, SCOTCH_Num * const iw, SCOTCH_Num * const nv, SCOTCH_Num * const elen, SCOTCH_Num * const last);

int                         esmumps_strat1      (const SCOTCH_Num procnbr, const SCOTCH_Num leafsiz, const int leorval, const SCOTCH_Num cminval, const SCOTCH_Num cmaxval, const double fratval, const int verbval, FILE * const stream, char * const straptr);

#endif /* LIB_ESMUMPS_H */
