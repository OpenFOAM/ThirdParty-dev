EXE		=
LIB		= .a
OBJ		= .o

MAKE		= make
AR		= ar
ARFLAGS		= -ruv
CAT		= cat
CCS		= gcc
CCP		= mpicc
CCD		= gcc
SCOTCH_NAME_SUFFIX=
CFLAGS		= -Wall --coverage -fdiagnostics-show-option `pkg-config --cflags ompi-c` -Drestrict=__restrict -DSCOTCH_GITLAB_SEPARATOR -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_FILE_COMPRESS_BZ2 -DSCOTCH_RENAME
LDFLAGS		= --coverage `pkg-config --libs --static ompi-c` -lz -lbz2 -llzma -lm -lrt -pthread
CP		= cp
FLEX		= flex
LN		= ln
MKDIR		= mkdir
MV		= mv
RANLIB		= ranlib
BISON		= bison
