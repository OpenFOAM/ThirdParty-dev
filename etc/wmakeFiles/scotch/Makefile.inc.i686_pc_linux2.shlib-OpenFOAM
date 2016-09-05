EXE		=
LIB		= .so
OBJ		= .o

MAKE		= make
AR		= gcc
ARFLAGS		= $(WM_CFLAGS) -shared -o
CAT		= cat
CCS		= gcc
CCP		= mpicc
CCD		= mpicc
CFLAGS		= $(WM_CFLAGS) -O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Drestrict=__restrict

ifeq ($(WM_LABEL_OPTION),Int64)
CFLAGS		+= -DINTSIZE64
endif

CLIBFLAGS	= -shared
LDFLAGS		= -Xlinker --no-as-needed $(WM_LDFLAGS) -lz -lm -lrt
CP		= cp
LEX		= flex -Pscotchyy -olex.yy.c
LN		= ln
MKDIR		= mkdir
MV		= mv
RANLIB		= echo
YACC		= bison -pscotchyy -y -b y
