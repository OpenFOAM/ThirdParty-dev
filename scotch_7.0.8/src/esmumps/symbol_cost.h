/* Copyright 1998,2020,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : symbol_cost.h                           **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the symbolic matrix cost computing  **/
/**                routine.                                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 14 oct 1998     **/
/**                                 to   : 16 oct 1998     **/
/**                # Version 6.1  : from : 28 aug 2020     **/
/**                                 to   : 28 aug 2020     **/
/**                # Version 7.0  : from : 21 jan 2023     **/
/**                                 to   : 21 jan 2023     **/
/**                                                        **/
/************************************************************/

/*
**  The function prototypes.
*/

#ifdef ESMUMPS_SYMBOL_COST

static void                 symbolCost2         (const SymbolCblk * const cblktax, const SymbolBlok * const bloktax, const Dof * const deofptr, double * const nnzptr, double * const opcptr, const INT cblkmin, const INT cblknbr);

#endif /* ESMUMPS_SYMBOL_COST */
