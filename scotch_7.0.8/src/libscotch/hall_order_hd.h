/* Copyright 2004,2007,2012,2018,2019 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : hall_order_hd.h                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declara-  **/
/**                tions for the block-oriented Halo       **/
/**                Approximate (Multiple) Minimum Degree   **/
/**                ordering routine.                       **/
/**                                                        **/
/**   DATES      : # Version 3.4  : from : 15 may 2001     **/
/**                                 to   : 15 may 2001     **/
/**                # Version 4.0  : from : 10 jan 2003     **/
/**                                 to   : 10 dec 2003     **/
/**                # Version 6.0  : from : 08 mar 2012     **/
/**                                 to   : 07 jun 2018     **/
/**                # Version 6.1  : from : 01 nov 2019     **/
/**                                 to   : 01 nov 2019     **/
/**                                                        **/
/************************************************************/

/*
**  The function prototypes.
*/

void                        hallOrderHdR2Halmd  (const Gnum, const Gnum, const Gnum, const Gnum, Gnum * const, Gnum, Gnum * const, Gnum * const, Gnum * const, Gnum * const, Gnum * const, Gnum * const, Gnum * const, Gnum * const, Gnum * const, Gnum * const);
