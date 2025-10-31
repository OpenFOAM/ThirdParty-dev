/* Copyright 2008,2018,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : common_file_compress.h                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the file (de)compression routines.  **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 12 mar 2008     **/
/**                                 to   : 17 mar 2008     **/
/**                # Version 6.0  : from : 08 jul 2018     **/
/**                                 to   : 14 jul 2018     **/
/**                # Version 7.0  : from : 19 jan 2023     **/
/**                                 to   : 06 aug 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/* Buffer size. */

#define FILECOMPRESSDATASIZE        (128 * 1024) /* Size of (de)compressing buffers */

/* Available types of (de)compression. */

typedef enum FileCompressType_ {
  FILECOMPRESSTYPENOTIMPL = -1,                   /* Error code     */
  FILECOMPRESSTYPENONE    = 0,                    /* No compression */
  FILECOMPRESSTYPEBZ2,
  FILECOMPRESSTYPEGZ,
  FILECOMPRESSTYPELZMA
} FileCompressType;

/* (De)compression type slot. */

typedef struct FileCompressTab_ {
  char *                    name;                 /* File extension name  */
  FileCompressType          type;                 /* (De)compression type */
} FileCompressTab;

/*
**  The type and structure definitions.
*/

typedef struct FileCompress_ {
  FileCompressType          typeval;              /*+ Type of (de)compression      +*/
  int                       infdnum;              /*+ Inner file handle (pipe end) +*/
  FILE *                    oustptr;              /*+ Outer stream                 +*/
  byte *                    bufftab;              /*+ Data buffer                  +*/
#ifdef COMMON_PTHREAD_FILE
  pthread_t                 thrdval;              /*+ Spawned thread ID            +*/
#else /* COMMON_PTHREAD_FILE */
  int                       procval;              /*+ Forked process ID            +*/
#endif /* COMMON_PTHREAD_FILE */
} FileCompress;

/*
**  The function prototypes.
*/

#ifdef COMMON_FILE_COMPRESS_BZ2
#ifdef SCOTCH_COMMON_FILE_COMPRESS
static void                 fileCompressBz2     (FileCompress * const  dataptr);
#endif /* SCOTCH_COMMON_FILE_COMPRESS */
#ifdef SCOTCH_COMMON_FILE_DECOMPRESS
static void                 fileDecompressBz2   (FileCompress * const  dataptr);
#endif /* SCOTCH_COMMON_FILE_DECOMPRESS */
#endif /* COMMON_FILE_COMPRESS_Bz2 */
#ifdef COMMON_FILE_COMPRESS_GZ
#ifdef SCOTCH_COMMON_FILE_COMPRESS
static void                 fileCompressGz      (FileCompress * const  dataptr);
#endif /* SCOTCH_COMMON_FILE_COMPRESS */
#ifdef SCOTCH_COMMON_FILE_DECOMPRESS
static void                 fileDecompressGz    (FileCompress * const  dataptr);
#endif /* SCOTCH_COMMON_FILE_DECOMPRESS */
#endif /* COMMON_FILE_COMPRESS_GZ */
#ifdef COMMON_FILE_COMPRESS_LZMA
#ifdef SCOTCH_COMMON_FILE_COMPRESS
static void                 fileCompressLzma    (FileCompress * const  dataptr);
#endif /* SCOTCH_COMMON_FILE_COMPRESS */
#ifdef SCOTCH_COMMON_FILE_DECOMPRESS
static void                 fileDecompressLzma  (FileCompress * const  dataptr);
#endif /* SCOTCH_COMMON_FILE_DECOMPRESS */
#endif /* COMMON_FILE_COMPRESS_LZMA */

int                         fileCompress        (File * const, const int);
void                        fileCompressExit    (File * const);
int                         fileCompressType    (const char * const);
int                         fileDecompress      (File * const, const int);
int                         fileDecompressType  (const char * const);
