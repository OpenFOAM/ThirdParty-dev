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
/**   NAME       : pthread_win32.h                         **/
/**                                                        **/
/**   AUTHORS    : Clement BARTHELEMY                      **/
/**                Tetsuya MISHIMA                         **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are an implementation of a  **/
/**                (very limited) subset of the pthread    **/
/**                API for win32.                          **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 15 jun 2024     **/
/**                                 to   : 08 aug 2024     **/
/**                                                        **/
/**   NOTES      : # This file is based on earlier work    **/
/**                  by Samuel THIBAULT for the StarPU     **/
/**                  software.                             **/
/**                # No effort is made to translate win32  **/
/**                  error codes to Posix errno's: most    **/
/**                  functions return EINVAL on error.     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

typedef HANDLE pthread_t;
typedef DWORD pthread_attr_t;

typedef unsigned pthread_condattr_t;
typedef CONDITION_VARIABLE pthread_cond_t;

typedef unsigned pthread_mutexattr_t;
typedef struct {
  CRITICAL_SECTION          lock;
  INIT_ONCE                 init;
} pthread_mutex_t;

#define PTHREAD_MUTEX_INITIALIZER { .lock = {0}, .init = INIT_ONCE_STATIC_INIT }

/*
** The pthread compatibility routines.
*/

static inline
int
pthread_create (
pthread_t *             thrdptr,
const pthread_attr_t *  attrptr,
void *               (* funcptr) (void *),
void *                  argptr)
{
  if ((attrptr != NULL) && (*attrptr != 0))
    return (EINVAL);

  *thrdptr = CreateThread (NULL, 0, (LPTHREAD_START_ROUTINE) funcptr, argptr, 0, NULL);
  if (*thrdptr == NULL)
    return (EAGAIN);

  return (0);
}

static inline
pthread_t
pthread_self (void)
{
  return (GetCurrentThread ());
}

static inline
int
pthread_join (
pthread_t                   thrdval,
void **                     resuptr)
{
  switch (WaitForSingleObject (thrdval, INFINITE)) {
    case WAIT_OBJECT_0 :
      break;
    case WAIT_FAILED :                            /* FALL THROUGH */
    default :
      return (EINVAL);
  }
  if (resuptr != NULL) {
    DWORD               resudat;

    if (GetExitCodeThread (thrdval, &resudat))
      *resuptr = (void *) (DWORD_PTR) resudat;
  }

  return (0);
}

static inline
int
pthread_detach (
pthread_t                   thrdval)
{
  if (! CloseHandle (thrdval))
    return (EINVAL);

  return (0);
}

static inline
void
pthread_exit (
void *                      resptr)
{
  ExitThread ((DWORD) (DWORD_PTR) resptr);
}

/*
** The pthread_cond routines.
*/

static inline
int
pthread_cond_init (
pthread_cond_t *            condptr,
const pthread_condattr_t *  attrptr)
{
  if ((attrptr != NULL) && (*attrptr != 0))
    return (EINVAL);
  InitializeConditionVariable (condptr);

  return (0);
}

static inline
int
pthread_cond_wait (
pthread_cond_t *            condptr,
pthread_mutex_t *           muteptr)
{
  if (SleepConditionVariableCS (condptr, &muteptr->lock, INFINITE) != 0)
    return (EINVAL);

  return (0);
}

static inline
int
pthread_cond_broadcast (
pthread_cond_t *            condptr)
{
  WakeAllConditionVariable (condptr);

  return (0);
}

static inline
int
pthread_cond_destroy (
pthread_cond_t *            condptr)              /* Not used */
{
  return (0);
}

/*
** The pthread_mutex routines.
*/

static
BOOL
initCriticalSectionOnce (
INIT_ONCE *                 initOnce,             /* Not used */
void *                      lockptr,
void **                     contptr)              /* Not used */
{
  InitializeCriticalSection (lockptr);

  return (TRUE);
}

static inline
int
pthread_mutex_init (
pthread_mutex_t *           muteptr,
pthread_mutexattr_t *       attrptr)
{
  if ((attrptr != NULL) && (*attrptr != 0))
    return (EINVAL);
  InitOnceInitialize  (&muteptr->init);
  InitOnceExecuteOnce (&muteptr->init, initCriticalSectionOnce, &muteptr->lock, NULL);

  return (0);
}

static inline
int
pthread_mutex_lock (
pthread_mutex_t *           muteptr)
{
  InitOnceExecuteOnce  (&muteptr->init, initCriticalSectionOnce, &muteptr->lock, NULL);
  EnterCriticalSection (&muteptr->lock);

  return (0);
}

static inline
int
pthread_mutex_unlock (
pthread_mutex_t *           muteptr)
{
  InitOnceExecuteOnce  (&muteptr->init, initCriticalSectionOnce, &muteptr->lock, NULL);
  LeaveCriticalSection (&muteptr->lock);

  return (0);
}

static inline
int
pthread_mutex_destroy (
pthread_mutex_t *           muteptr)
{
  DeleteCriticalSection (&muteptr->lock);

  return (0);
}
