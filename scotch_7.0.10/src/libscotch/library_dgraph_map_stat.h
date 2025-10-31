/* Copyright 2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_map_stat.h               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declara-  **/
/**                tions for the library distributed       **/
/**                mapping statistics routine.             **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 19 sep 2024     **/
/**                                 to   : 19 sep 2024     **/
/**                                                        **/
/************************************************************/

/*
** The defines.
*/

#define SCOTCH_LIBRARY_DGRAPH_MAP_STAT_H

/* The statistics category computation flags. */

#define LIBDGRAPHMAPSTATNONE        0x0000        /* No computations to perform       */
#define LIBDGRAPHMAPSTATCOMM        0x0001        /* Compute communication statistics */
#define LIBDGRAPHMAPSTATNGHB        0x0002        /* Compute statistics on neighbors  */
