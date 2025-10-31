/* Copyright 2019,2020 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_hamf_graph.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the original Fortran  **/
/**                HAMF code provided by the MUMPS team.   **/
/**                                                        **/
/**   DATES      : # Version 6.1  : from : 24 nov 2019     **/
/**                                 to   : 10 feb 2020     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "hgraph.h"
#include "hgraph_order_hx.h"
#include "hall_order_hf.h"

#define HGRAPHORDERHFCOMPRAT        1.2L          /*+ Compression ratio +*/

/*
**  The copy-pasted code.
*/

Gnum mfverttab[] = {
  1, 26, 52, 113, 150, 200, 227, 254, 282, 329, 352, 389, 468, 477, 486, 495,
  512, 533, 550, 559, 576, 589, 602, 619, 640, 657, 670, 687, 700, 717, 734, 751,
  768, 781, 802, 815, 832, 845, 862, 879, 892, 905, 922, 939, 964, 981, 998, 1011,
  1032, 1045, 1062, 1075, 1088, 1101, 1114, 1127, 1140, 1157, 1170, 1191, 1208, 1221, 1243, 1261,
  1283, 1301, 1327, 1353, 1371, 1389, 1415, 1437, 1455, 1473, 1499, 1521, 1543, 1565, 1583, 1609,
  1635, 1657, 1675, 1701, 1719, 1741, 1759, 1777, 1795, 1829, 1851, 1877, 1895, 1917, 1935, 1953,
  1971, 1993, 2011, 2041, 2059, 2085, 2107, 2117, 2129, 2150, 2169, 2189, 2197, 2205, 2221, 2230,
  2246, 2264, 2270, 2274, 2283, 2290, 2307, 2311, 2320, 2329, 2338, 2347, 2351, 2377, 2387, 2407,
  2411, 2420, 2432, 2453, 2465, 2470, 2479, 2483, 2501, 2512, 2522, 2528, 2533, 2537, 2541, 2545,
  2554, 2559, 2564, 2570, 2576, 2581, 2587, 2593, 2600, 2609, 2617, 2618, 2619, 2620, 2621, 2623,
  2625, 2628, 2630, 2633, 2635, 2638, 2641, 2643, 2646, 2649, 2651, 2652, 2653
};

Gnum mfvelotab[] = {
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
};

Gnum mfedgetab[] = {
  2, 3, 4, 5, 13, 14, 17, 22, 24, 25, 26, 29, 30, 33, 34, 50,
  52, 65, 67, 73, 77, 106, 105, 104, 103, 1, 4, 5, 6, 12, 14, 15,
  18, 23, 25, 26, 28, 30, 36, 50, 53, 65, 71, 79, 84, 89, 90, 94,
  109, 108, 107, 1, 4, 5, 9, 10, 11, 12, 13, 16, 17, 19, 22, 24,
  26, 29, 33, 34, 41, 42, 43, 48, 49, 50, 52, 57, 59, 60, 61, 62,
  65, 67, 68, 70, 73, 75, 77, 79, 81, 82, 86, 87, 88, 89, 93, 95,
  96, 99, 100, 102, 106, 105, 104, 117, 116, 115, 114, 113, 112, 111, 103, 110,
  1, 2, 3, 5, 11, 12, 14, 17, 18, 23, 24, 25, 26, 30, 48, 50,
  60, 61, 65, 70, 74, 77, 78, 79, 82, 84, 86, 89, 90, 92, 93, 121,
  120, 119, 118, 107, 110, 1, 2, 3, 4, 6, 7, 8, 9, 12, 14, 17,
  24, 25, 26, 30, 31, 32, 34, 36, 37, 39, 44, 50, 52, 53, 54, 58,
  63, 65, 71, 72, 73, 76, 77, 79, 80, 81, 84, 86, 87, 89, 96, 97,
  99, 126, 106, 125, 124, 123, 122, 2, 5, 7, 12, 15, 23, 27, 28, 30,
  36, 37, 38, 53, 54, 71, 76, 84, 89, 90, 94, 101, 129, 109, 128, 127,
  108, 107, 5, 6, 8, 12, 27, 32, 35, 37, 38, 39, 40, 51, 53, 54,
  63, 64, 71, 76, 83, 89, 91, 98, 101, 132, 131, 130, 127, 5, 7, 12,
  31, 32, 35, 39, 40, 54, 58, 63, 64, 72, 76, 80, 83, 89, 97, 98,
  126, 132, 125, 131, 135, 134, 130, 133, 122, 3, 5, 10, 11, 12, 34, 41,
  42, 44, 45, 46, 47, 55, 56, 57, 62, 66, 67, 68, 70, 73, 75, 77,
  80, 81, 86, 87, 88, 89, 93, 96, 97, 99, 100, 102, 106, 105, 125, 104,
  141, 140, 139, 123, 138, 113, 137, 136, 3, 9, 11, 41, 42, 43, 45, 47,
  49, 57, 59, 68, 70, 75, 88, 99, 100, 116, 138, 142, 113, 137, 112, 3,
  4, 9, 10, 12, 17, 19, 41, 42, 43, 46, 48, 49, 56, 57, 59, 60,
  61, 66, 70, 74, 78, 79, 81, 82, 86, 88, 92, 93, 99, 112, 121, 144,
  143, 118, 136, 110, 2, 3, 4, 5, 6, 7, 8, 9, 11, 17, 18, 20,
  23, 25, 27, 30, 32, 36, 37, 38, 39, 42, 44, 46, 53, 54, 55, 56,
  63, 64, 65, 66, 69, 70, 71, 72, 74, 76, 77, 78, 79, 80, 81, 82,
  83, 84, 85, 86, 89, 90, 91, 92, 93, 94, 96, 97, 98, 99, 101, 154,
  125, 153, 131, 152, 134, 151, 150, 149, 148, 129, 147, 146, 145, 120, 144, 118,
  127, 136, 107, 1, 3, 22, 24, 29, 33, 67, 105, 103, 1, 2, 4, 5,
  25, 26, 30, 50, 65, 2, 6, 23, 28, 36, 94, 109, 108, 107, 3, 22,
  29, 59, 62, 67, 75, 95, 105, 117, 116, 115, 114, 113, 112, 111, 103, 1,
  3, 4, 5, 11, 12, 24, 26, 48, 50, 60, 61, 65, 70, 77, 79, 82,
  86, 89, 93, 110, 2, 4, 12, 23, 25, 60, 74, 78, 79, 84, 90, 121,
  120, 119, 118, 107, 110, 3, 11, 48, 49, 59, 61, 70, 112, 110, 12, 55,
  66, 69, 80, 85, 154, 125, 158, 153, 152, 157, 156, 150, 147, 155, 136, 51,
  91, 165, 132, 164, 131, 163, 162, 161, 160, 130, 159, 127, 1, 3, 13, 16,
  24, 29, 33, 67, 95, 105, 114, 111, 103, 2, 4, 6, 12, 15, 18, 25,
  28, 36, 71, 79, 84, 90, 94, 109, 108, 107, 1, 3, 4, 5, 13, 17,
  22, 26, 29, 33, 34, 50, 52, 65, 67, 73, 77, 106, 105, 104, 103, 1,
  2, 4, 5, 12, 14, 18, 23, 26, 30, 50, 65, 79, 84, 89, 90, 107,
  1, 2, 3, 4, 5, 14, 17, 24, 25, 30, 50, 65, 77, 6, 7, 12,
  28, 37, 38, 71, 76, 90, 94, 101, 129, 109, 128, 127, 108, 107, 2, 6,
  15, 23, 27, 36, 94, 129, 109, 128, 127, 108, 107, 1, 3, 13, 16, 22,
  24, 33, 67, 95, 105, 169, 168, 114, 167, 166, 111, 103, 1, 2, 4, 5,
  6, 12, 14, 25, 26, 36, 50, 53, 65, 71, 79, 84, 89, 5, 8, 32,
  40, 58, 64, 72, 97, 126, 132, 125, 131, 135, 134, 130, 133, 122, 5, 7,
  8, 12, 31, 39, 54, 58, 63, 72, 76, 80, 89, 97, 126, 125, 122, 1,
  3, 13, 22, 24, 29, 52, 67, 73, 106, 105, 104, 103, 1, 3, 5, 9,
  24, 44, 50, 52, 58, 73, 77, 87, 96, 97, 99, 126, 106, 125, 124, 123,
  122, 7, 8, 38, 39, 40, 51, 64, 91, 98, 132, 131, 130, 127, 2, 5,
  6, 12, 15, 23, 28, 30, 53, 71, 84, 89, 90, 94, 109, 108, 107, 5,
  6, 7, 12, 27, 38, 53, 54, 71, 76, 89, 101, 127, 6, 7, 12, 27,
  35, 37, 51, 71, 76, 83, 91, 98, 101, 132, 131, 130, 127, 5, 7, 8,
  12, 32, 35, 40, 54, 63, 64, 76, 83, 89, 98, 132, 131, 130, 7, 8,
  31, 35, 39, 64, 98, 126, 132, 131, 135, 130, 133, 3, 9, 10, 11, 42,
  43, 49, 57, 59, 70, 88, 99, 112, 3, 9, 10, 11, 12, 41, 46, 56,
  57, 66, 70, 81, 86, 88, 93, 99, 136, 3, 10, 11, 41, 47, 49, 59,
  68, 70, 75, 88, 116, 138, 142, 113, 137, 112, 5, 9, 12, 34, 45, 55,
  56, 66, 80, 81, 87, 89, 96, 97, 102, 106, 105, 125, 104, 141, 140, 139,
  123, 137, 136, 9, 10, 44, 47, 57, 62, 68, 100, 102, 105, 125, 141, 140,
  139, 138, 113, 137, 9, 11, 12, 42, 48, 56, 66, 74, 81, 92, 93, 121,
  144, 143, 118, 136, 110, 9, 10, 43, 45, 57, 68, 100, 116, 138, 142, 113,
  137, 112, 3, 4, 11, 17, 19, 46, 49, 59, 60, 61, 70, 78, 82, 92,
  112, 121, 144, 143, 118, 136, 110, 3, 10, 11, 19, 41, 43, 48, 59, 61,
  70, 88, 112, 110, 1, 2, 3, 4, 5, 14, 17, 24, 25, 26, 30, 34,
  52, 65, 73, 77, 106, 7, 21, 35, 38, 91, 98, 165, 132, 131, 163, 161,
  130, 127, 1, 3, 5, 24, 33, 34, 50, 67, 73, 77, 106, 105, 104, 2,
  5, 6, 7, 12, 30, 36, 37, 54, 71, 76, 84, 89, 5, 6, 7, 8,
  12, 32, 37, 39, 53, 63, 71, 76, 89, 9, 12, 20, 44, 56, 66, 69,
  80, 81, 154, 125, 152, 136, 9, 11, 12, 42, 44, 46, 55, 66, 80, 81,
  93, 125, 136, 3, 9, 10, 11, 41, 42, 45, 47, 68, 70, 75, 88, 99,
  100, 138, 113, 137, 5, 8, 31, 32, 34, 72, 97, 126, 106, 125, 124, 123,
  122, 3, 10, 11, 16, 19, 41, 43, 48, 49, 61, 68, 70, 75, 88, 117,
  116, 115, 113, 112, 111, 110, 3, 4, 11, 17, 18, 48, 61, 70, 78, 82,
  92, 121, 120, 119, 118, 107, 110, 3, 4, 11, 17, 19, 48, 49, 59, 60,
  70, 82, 112, 110, 3, 9, 16, 45, 67, 75, 95, 99, 100, 102, 105, 169,
  139, 171, 170, 117, 168, 138, 113, 137, 166, 111, 5, 7, 8, 12, 32, 39,
  54, 64, 72, 76, 80, 83, 89, 97, 98, 125, 131, 134, 7, 8, 12, 31,
  35, 39, 40, 63, 72, 76, 80, 83, 98, 126, 132, 125, 131, 135, 134, 130,
  133, 122, 1, 2, 3, 4, 5, 12, 14, 17, 24, 25, 26, 30, 50, 77,
  79, 84, 86, 89, 9, 11, 12, 20, 42, 44, 46, 55, 56, 69, 74, 80,
  81, 85, 92, 93, 154, 125, 153, 152, 150, 147, 145, 144, 118, 136, 1, 3,
  9, 13, 16, 22, 24, 29, 33, 52, 62, 73, 75, 87, 95, 99, 100, 102,
  106, 105, 104, 117, 114, 113, 111, 103, 3, 9, 10, 43, 45, 47, 57, 59,
  75, 88, 99, 100, 116, 138, 142, 113, 137, 112, 12, 20, 55, 66, 80, 83,
  85, 154, 125, 153, 131, 152, 134, 151, 150, 149, 147, 136, 3, 4, 9, 10,
  11, 12, 17, 19, 41, 42, 43, 48, 49, 57, 59, 60, 61, 79, 81, 82,
  86, 88, 93, 99, 112, 110, 2, 5, 6, 7, 12, 23, 27, 30, 36, 37,
  38, 53, 54, 76, 84, 89, 90, 94, 101, 129, 127, 107, 5, 8, 12, 31,
  32, 58, 63, 64, 80, 83, 89, 97, 126, 125, 131, 134, 133, 122, 1, 3,
  5, 9, 24, 33, 34, 50, 52, 67, 77, 87, 96, 99, 102, 106, 105, 104,
  4, 11, 12, 18, 46, 66, 78, 79, 82, 85, 90, 92, 93, 101, 153, 148,
  129, 147, 146, 145, 120, 144, 118, 127, 136, 107, 3, 9, 10, 16, 43, 57,
  59, 62, 67, 68, 88, 95, 99, 100, 102, 105, 117, 116, 115, 113, 112, 111,
  5, 6, 7, 8, 12, 27, 32, 37, 38, 39, 53, 54, 63, 64, 71, 83,
  89, 91, 98, 101, 131, 127, 1, 3, 4, 5, 9, 12, 17, 24, 26, 34,
  50, 52, 65, 73, 79, 81, 86, 87, 89, 96, 99, 106, 4, 11, 12, 18,
  48, 60, 74, 79, 82, 90, 92, 93, 121, 120, 119, 118, 107, 110, 2, 3,
  4, 5, 11, 12, 17, 18, 23, 25, 30, 65, 70, 74, 77, 78, 82, 84,
  86, 89, 90, 92, 93, 120, 118, 107, 5, 8, 9, 12, 20, 32, 44, 55,
  56, 63, 64, 66, 69, 72, 81, 83, 89, 96, 97, 154, 125, 131, 152, 134,
  149, 136, 3, 5, 9, 11, 12, 42, 44, 46, 55, 56, 66, 70, 77, 80,
  86, 89, 93, 96, 97, 99, 125, 136, 3, 4, 11, 12, 17, 48, 60, 61,
  70, 74, 78, 79, 86, 92, 93, 121, 118, 110, 7, 8, 12, 38, 39, 63,
  64, 69, 72, 76, 80, 85, 91, 98, 101, 154, 125, 153, 131, 152, 134, 151,
  150, 149, 148, 127, 2, 4, 5, 6, 12, 18, 23, 25, 30, 36, 53, 65,
  71, 79, 89, 90, 94, 107, 12, 20, 66, 69, 74, 83, 91, 101, 154, 153,
  131, 151, 150, 149, 148, 147, 146, 145, 144, 118, 127, 136, 3, 4, 5, 9,
  11, 12, 17, 42, 65, 70, 77, 79, 81, 82, 89, 93, 96, 99, 3, 5,
  9, 34, 44, 67, 73, 77, 96, 97, 99, 102, 106, 105, 125, 104, 140, 123,
  3, 9, 10, 11, 41, 42, 43, 49, 57, 59, 68, 70, 75, 99, 100, 116,
  113, 112, 2, 3, 4, 5, 6, 7, 8, 9, 12, 17, 25, 30, 32, 36,
  37, 39, 44, 53, 54, 63, 65, 71, 72, 76, 77, 79, 80, 81, 84, 86,
  96, 97, 99, 125, 2, 4, 6, 12, 18, 23, 25, 27, 36, 71, 74, 78,
  79, 84, 94, 101, 129, 146, 120, 118, 127, 107, 7, 12, 21, 35, 38, 51,
  76, 83, 85, 98, 101, 165, 132, 164, 153, 131, 163, 162, 161, 160, 151, 130,
  172, 148, 159, 127, 4, 11, 12, 46, 48, 60, 66, 74, 78, 79, 82, 93,
  121, 144, 143, 118, 136, 110, 3, 4, 9, 11, 12, 17, 42, 46, 56, 66,
  70, 74, 78, 79, 81, 82, 86, 92, 99, 144, 118, 136, 2, 6, 12, 15,
  23, 27, 28, 36, 71, 84, 90, 101, 129, 109, 128, 127, 108, 107, 3, 16,
  22, 29, 62, 67, 75, 105, 169, 170, 117, 168, 114, 167, 113, 166, 111, 103,
  3, 5, 9, 12, 34, 44, 73, 77, 80, 81, 86, 87, 89, 97, 99, 106,
  125, 123, 5, 8, 9, 12, 31, 32, 34, 44, 58, 63, 72, 80, 81, 87,
  89, 96, 126, 106, 125, 124, 123, 122, 7, 8, 12, 35, 38, 39, 40, 51,
  63, 64, 76, 83, 91, 101, 132, 131, 130, 127, 3, 5, 9, 10, 11, 12,
  34, 41, 42, 57, 62, 67, 68, 70, 73, 75, 77, 81, 86, 87, 88, 89,
  93, 96, 100, 102, 106, 105, 104, 113, 3, 9, 10, 45, 47, 57, 62, 67,
  68, 75, 88, 99, 102, 105, 139, 138, 113, 137, 6, 7, 12, 27, 37, 38,
  71, 74, 76, 83, 85, 90, 91, 94, 98, 153, 131, 151, 148, 129, 146, 145,
  120, 118, 127, 107, 3, 9, 44, 45, 62, 67, 73, 75, 87, 99, 100, 106,
  105, 125, 104, 141, 140, 139, 123, 138, 113, 137, 1, 3, 13, 16, 22, 24,
  29, 33, 67, 95, 1, 3, 9, 24, 33, 44, 52, 67, 73, 87, 99, 102,
  1, 3, 9, 13, 16, 22, 24, 29, 33, 44, 45, 52, 62, 67, 73, 75,
  87, 95, 99, 100, 102, 1, 3, 5, 9, 24, 33, 34, 44, 50, 52, 58,
  67, 73, 77, 87, 96, 97, 99, 102, 2, 4, 6, 12, 15, 18, 23, 25,
  27, 28, 36, 60, 71, 74, 78, 79, 84, 90, 94, 101, 2, 6, 15, 23,
  27, 28, 36, 94, 2, 6, 15, 23, 27, 28, 36, 94, 3, 4, 11, 17,
  18, 19, 46, 48, 49, 59, 60, 61, 70, 78, 82, 92, 3, 16, 22, 29,
  59, 62, 67, 75, 95, 3, 10, 11, 16, 19, 41, 43, 47, 48, 49, 59,
  61, 68, 70, 75, 88, 3, 9, 10, 16, 43, 45, 47, 57, 59, 62, 67,
  68, 75, 88, 95, 99, 100, 102, 3, 16, 22, 29, 67, 95, 3, 16, 59,
  75, 3, 10, 16, 43, 47, 59, 68, 75, 88, 3, 16, 59, 62, 67, 75,
  95, 4, 11, 12, 18, 46, 48, 60, 66, 74, 78, 79, 82, 85, 90, 92,
  93, 101, 4, 18, 60, 78, 4, 12, 18, 60, 74, 78, 79, 90, 101, 4,
  11, 18, 46, 48, 60, 78, 82, 92, 5, 8, 31, 32, 34, 58, 64, 72,
  97, 5, 9, 34, 44, 58, 87, 96, 97, 102, 5, 34, 58, 97, 5, 8,
  9, 12, 20, 31, 32, 34, 44, 45, 55, 56, 58, 63, 64, 66, 69, 72,
  80, 81, 83, 87, 89, 96, 97, 102, 5, 8, 31, 32, 34, 40, 58, 64,
  72, 97, 6, 7, 12, 21, 27, 28, 35, 37, 38, 51, 71, 74, 76, 83,
  85, 90, 91, 94, 98, 101, 6, 27, 28, 94, 6, 12, 27, 28, 71, 74,
  90, 94, 101, 7, 8, 21, 31, 35, 38, 39, 40, 51, 64, 91, 98, 7,
  8, 12, 21, 31, 35, 38, 39, 40, 51, 63, 64, 69, 72, 76, 80, 83,
  85, 91, 98, 101, 7, 8, 21, 31, 35, 38, 39, 40, 51, 64, 91, 98,
  8, 31, 40, 64, 72, 8, 12, 31, 63, 64, 69, 72, 80, 83, 8, 31,
  40, 64, 9, 11, 12, 20, 42, 44, 46, 48, 55, 56, 66, 69, 74, 80,
  81, 85, 92, 93, 9, 10, 43, 44, 45, 47, 57, 62, 68, 100, 102, 9,
  10, 43, 45, 47, 57, 62, 68, 100, 102, 9, 44, 45, 62, 100, 102, 9,
  44, 45, 87, 102, 9, 44, 45, 102, 10, 43, 47, 68, 11, 46, 48, 92,
  11, 12, 46, 48, 66, 74, 85, 92, 93, 12, 66, 74, 85, 101, 12, 74,
  85, 90, 101, 12, 20, 66, 69, 74, 85, 12, 74, 83, 85, 91, 101, 12,
  69, 80, 83, 85, 12, 20, 66, 69, 83, 85, 12, 69, 83, 85, 91, 101,
  12, 20, 55, 66, 69, 80, 83, 12, 20, 66, 69, 74, 83, 85, 91, 101,
  12, 20, 55, 66, 69, 80, 83, 85, 20, 20, 20, 20, 21, 91, 21, 91,
  21, 51, 91, 21, 91, 21, 51, 91, 21, 91, 21, 51, 91, 29, 62, 95,
  29, 95, 29, 62, 95, 29, 62, 95, 62, 95, 62, 91
};

Gnum mfvnhdtab[] = {
  22, 49, 101, 144, 194, 221, 250, 273, 317, 346, 382, 448, 475, 486, 492, 503,
  532, 544, 557, 565, 578, 598, 616, 636, 656, 670, 681, 694, 709, 734, 742, 765,
  777, 796, 811, 829, 844, 858, 876, 886, 904, 921, 933, 954, 973, 992, 1005, 1025,
  1043, 1061, 1068, 1085, 1101, 1114, 1123, 1138, 1154, 1164, 1184, 1202, 1219, 1231, 1258, 1274,
  1301, 1317, 1345, 1365, 1378, 1413, 1434, 1449, 1470, 1487, 1514, 1541, 1564, 1577, 1606, 1628,
  1655, 1672, 1690, 1718, 1727, 1759, 1771, 1792, 1828, 1845, 1862, 1889, 1914, 1929, 1942, 1968,
  1987, 2007, 2037, 2054, 2074, 2096
};

void
mfHgraphBuild (
Hgraph *              grafptr)
{
  hgraphInit (grafptr);
  grafptr->s.baseval = 1;
  grafptr->s.vertnbr = 172;
  grafptr->s.vertnnd = 173;
  grafptr->s.verttax = mfverttab - 1;
  grafptr->s.vendtax = grafptr->s.verttax + 1;
  grafptr->s.velotax = mfvelotab - 1;
  grafptr->s.velosum = 516;
  grafptr->s.edgenbr = 2652;
  grafptr->s.edgetax = mfedgetab - 1;
  grafptr->s.edlosum = 2652;
  grafptr->s.degrmax = 79;
  grafptr->vnohnbr = 102;
  grafptr->vnohnnd = 103;
  grafptr->vnhdtax = mfvnhdtab - 1;
  grafptr->enohnbr = 1560;
  grafptr->enlosum = 1560;
}

/*
**  The function prototypes.
*/

void                        mumps_hamf4_ (int *, int *, int *, int *, int64_t *, int64_t *, int64_t *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);

/**********************/
/*                    */
/* The test routines. */
/*                    */
/**********************/

/* Test method 0
*/

int
test_hamf_meth0 (
Hgraph *                    grafptr)              /* Graph to test */
{
  Gnum                vertnum;
  Gnum                vertadj;
  Gnum                vertnew;
  int64_t             edgenew;
  int                 norig;
  int                 n;
  int                 nbelts;
  int                 nbbuck;
  int                 ncmpa;
  int64_t             iwlen;
  int64_t             pfree;
  int *               degtab;
  int *               lentab;
  int *               iwtab;
  int *               nvtab;
  int *               elentab;
  int *               headtab;
  int *               nexttab;
  int *               lasttab;
  int *               wtab;
  int *               wftab;
  int64_t *           petab;

  const Gnum * restrict const verttax = grafptr->s.verttax;
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const velotax = grafptr->s.velotax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;

  n      = grafptr->s.vertnbr;
  norig  = grafptr->s.velosum;
  nbelts = 0;
  nbbuck = norig * 2;
  iwlen  = (int64_t) ((double) grafptr->s.edgenbr * 1.2) + 32;

  petab   = malloc (n * sizeof (int64_t));
  lentab  = malloc (n * sizeof (int));
  nvtab   = malloc (n * sizeof (int));
  elentab = malloc (n * sizeof (int));
  lasttab = malloc (n * sizeof (int));
  degtab  = malloc (n * sizeof (int));
  wftab   = malloc (n * sizeof (int));
  nexttab = malloc (n * sizeof (int));
  wtab    = malloc (n * sizeof (int));
  headtab = malloc ((nbbuck + 2) * sizeof (int));
  iwtab   = malloc (iwlen * sizeof (int));

  int64_t * restrict const  petax   = petab   - 1; /* Base HAMF arrays at base 1 */
  int * restrict const      iwtax   = iwtab   - 1;
  int * restrict const      lentax  = lentab  - 1;
  int * restrict const      nvtax   = nvtab   - 1;
  int * restrict const      elentax = elentab - 1;

  vertadj = 1 - grafptr->s.baseval;
  for (vertnum = grafptr->s.baseval, vertnew = edgenew = 1; /* Process non-halo vertices */
       vertnum < grafptr->vnohnnd; vertnum ++, vertnew ++) {
    int                       degrval;
    int                       edgenum;

    degrval = vendtax[vertnum] - verttax[vertnum];
    petax[vertnew]   = (int64_t) edgenew;
    lentax[vertnew]  = degrval;
    elentax[vertnew] = 0;
    nvtax[vertnew]   = (velotax != NULL) ? velotax[vertnum] : 1;

    for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++, edgenew ++)
      iwtax[edgenew] = edgetax[edgenum] + vertadj;
  }
  for ( ; vertnum < grafptr->s.vertnnd; vertnum ++, vertnew ++) { /* Process halo vertices */
    int                       degrval;
    int                       edgenum;

    degrval = grafptr->s.verttax[vertnum] - grafptr->s.vendtax[vertnum]; /* Negative degree */
    petax[vertnew]   = edgenew;
    lentax[vertnew]  = (degrval != 0) ? degrval : (-1 - grafptr->s.vertnbr);
    elentax[vertnew] = 0;
    nvtax[vertnew]   = (velotax != NULL) ? velotax[vertnum] : 1;

    for (edgenum = grafptr->s.verttax[vertnum];
         edgenum < grafptr->s.vendtax[vertnum]; edgenum ++, edgenew ++)
      iwtax[edgenew] = edgetax[edgenum] + vertadj;
  }

  pfree = edgenew;                             /* Set index to first free area */

  mumps_hamf4_ (&norig, &n, &nbelts, &nbbuck,
                &iwlen, petab, &pfree, lentab, iwtab, nvtab, elentab,
                lasttab, &ncmpa, degtab, wftab, nexttab, wtab, headtab);

  if (ncmpa < 0) {
    errorPrint ("test_hamf_meth0: ordering method internal error");
    return     (1);
  }

  memFree (iwtab);
  memFree (headtab);
  memFree (wtab);
  memFree (nexttab);
  memFree (wftab);
  memFree (degtab);
  memFree (lasttab);
  memFree (elentab);
  memFree (nvtab);
  memFree (lentab);
  memFree (petab);

  return (0);
}

/* Test method 1
*/

int
test_hamf_meth1 (
Hgraph *                    grafptr)              /* Graph to test */
{
  Gnum                n;                          /* Number of nodes to order (with halo or not) */
  Gnum                norig;                      /* Number of nodes in uncompressed graph       */
  Gnum                nbbuck;
  Gnum * restrict     petab;
  Gnum                pfree;
  Gnum * restrict     lentab;
  Gnum                iwlen;
  Gnum * restrict     iwtab;
  Gnum * restrict     nvtab;
  Gnum * restrict     elentab;
  Gnum * restrict     lasttab;
  Gnum * restrict     leaftab;
  Gnum * restrict     secntab;                    /* Array of index to first secondary variable  */
  Gnum * restrict     nexttab;                    /* Array of index of next principal variable   */
  Gnum * restrict     frsttab;
  Gnum * restrict     headtab;                    /* Head array : nbbuck = 2 * n                 */
  Gnum * restrict     cwgttax;                    /* Column weight array                         */
  Gnum                cwgtsiz;
  Gnum                ncmpa;

  n      = grafptr->s.vertnbr;
  norig  = grafptr->s.velosum;
  nbbuck = norig * 2;
  iwlen  = (Gnum) ((double) grafptr->s.edgenbr * HGRAPHORDERHFCOMPRAT) + 32;
  if (iwlen < n)                                  /* TRICK: make sure to be able to re-use array */
    iwlen = n;
  cwgtsiz = (grafptr->s.velotax != NULL) ? n : 0;

  if (memAllocGroup ((void **) (void *)
                     &petab,   (size_t) (n * sizeof (Gnum)),
                     &lentab,  (size_t) (n * sizeof (Gnum)),
                     &nvtab,   (size_t) (n * sizeof (Gnum)),
                     &elentab, (size_t) (n * sizeof (Gnum)),
                     &lasttab, (size_t) (n * sizeof (Gnum)),
                     &leaftab, (size_t) (n * sizeof (Gnum)),
                     &frsttab, (size_t) (n * sizeof (Gnum)),
                     &secntab, (size_t) (n * sizeof (Gnum)),
                     &nexttab, (size_t) (n * sizeof (Gnum)),
                     &cwgttax, (size_t) (cwgtsiz * sizeof (Gnum)), /* Not based yet */
                     &headtab, (size_t) ((nbbuck + 2) * sizeof (Gnum)),
                     &iwtab,   (size_t) (iwlen * sizeof (Gnum)), NULL) == NULL) {
    errorPrint ("test_hamf_meth1: out of memory");
    return     (1);
  }

  hgraphOrderHxFill (grafptr, petab, lentab, iwtab, nvtab, elentab, &pfree);

  hallOrderHfR3Hamdf4 (norig, n, 0, nbbuck, iwlen, petab, pfree,
                       lentab, iwtab, nvtab, elentab, lasttab, &ncmpa,
                       leaftab, secntab, nexttab, frsttab, headtab);
  if (ncmpa < 0) {
    errorPrint ("test_hamf_meth1: ordering method internal error");
    return     (1);
  }

  memFree (petab);                                /* Free group leader */

  return (0);
}

int
test_hamf_graph (
int                         methval)
{
  Hgraph              halgrafdat;
  int                 o;

  mfHgraphBuild (&halgrafdat);
  hgraphCheck (&halgrafdat);

  switch (methval) {
    case 0 :
      o = test_hamf_meth0 (&halgrafdat);
      break;
    case 1 :
      o = test_hamf_meth1 (&halgrafdat);
      break;
    default :
      errorPrint ("test_hamf_graph: invalid test method");
      return     (1);
  }

  return (o);
}

/*********************/
/*                   */
/* The main routine. */
/*                   */
/*********************/

int
main (
int                 argc,
char *              argv[])
{
  int                 methval;

  if (argc > 2) {
    errorPrint ("usage: %s [method]", argv[0]);
    exit (EXIT_FAILURE);
  }

  methval = 0;
  if (argc == 2)
    methval = atoi (argv[1]);

  test_hamf_graph (methval);

  exit (EXIT_SUCCESS);
}
