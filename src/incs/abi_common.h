/* abi_common.h */

/*
 * Copyright (C) 2008-2012 ABINIT Group (MG)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef _ABINIT_COMMON_H
#define _ABINIT_COMMON_H

/*
 * Language standards requires the existance of pre-defined macros
 * Microsoft Visual C++ does not define __STDC__,
 * Sun Workshop 4.2 supports C94 without setting __STDC_VERSION__ to the proper value
 */

#if defined (__STDC__)
# define PREDEF_STANDARD_C_1989    /** ANSI X3.159-1989 **/
# if defined (__STDC_VERSION__)
#  define PREDEF_STANDARD_C_1990   /** ISO/IEC 9899:1990 **/
#  if (__STDC_VERSION__ >= 199409L)
#   define PREDEF_STANDARD_C_1994  /** ISO/IEC 9899-1:1994 **/
#  endif
#  if (__STDC_VERSION__ >= 199901L)
#   define PREDEF_STANDARD_C_1999  /** ISO/IEC 9899:1999 **/
#  endif
# endif
#endif

/** #define DEBUG_MODE **/

/** define WHEREARG __FILE__,__LINE__**/
#ifdef HAVE_FC_LONG_LINES
# define NEWLINE ;
#else
# define NEWLINE \newline
#endif
/** define WHEREARG NEWLINE __FILE__, NEWLINE __LINE__ **/

/** this does not work with gfort, pgi, **/

#if defined (FC_GNU) || defined(FC_G95) || defined (FC_PGI)
#define QUOTEME(x)     'x'
#else
#define QUOTEME(x)     #x
#endif

#define BYTE_SIZE(array)  PRODUCT(SHAPE(array)) * DBLE(KIND(array))

/** debugging macros so we can pin down message provenance at a glance
#define WHERESTR "[file %s, line %d] "
**/

/*
 * ABI_ basic abinit macros.
 * DBG_ denote macros for debugging. Defined only if abinit is compiled in DEBUG_MODE.
 * MSG_ denote macros for logging.
 * MEM_ for memory profiling and checking.
 * */

#ifdef HAVE_FC_LONG_LINES
#  define ABI_CHECK(expr,str) call assert((expr), str,__FILE__,__LINE__)
#  define ABI_DIE(msg)     call die(msg,__FILE__,__LINE__)
#else
#  define ABI_CHECK(expr,str) call assert((expr), str)
#  define ABI_DIE(msg)     call die(msg)
#endif

#if defined HAVE_FC_LONG_LINES
#  define ABI_CHECK_MPI(ierr,msg)      call check_mpi_ierr(ierr,msg,"PERS",__FILE__,__LINE__)
#  define ABI_CHECK_MPI_PERS(ierr,msg) call check_mpi_ierr(ierr,msg,"PERS",__FILE__,__LINE__)
#else
#  define ABI_CHECK_MPI(ierr,msg) call check_mpi_ierr(ierr,msg,"PERS")
#  define ABI_CHECK_MPI_PERS(ierr,msg) call check_mpi_ierr(ierr,msg,"PERS")
#endif


/* Macros for memory checking and profiling *
 * ERR_ALLOCATE_STAT is defined in the module m_errors.
 * Remember: subtract number of newlines in lineno otherwise the value might be misleading
#define HAVE_MEM_PROFILE 0
*/

/*#define HAVE_MEM_PROFILING*/
#ifdef HAVE_MEM_PROFILING
# ifdef HAVE_FC_LONG_LINES
#  define ABI_ALLOCATE(ARR,SIZE) \
   allocate(ARR SIZE,stat=ABI_ALLOC_STAT) NEWLINE \
   call memocc(ABI_ALLOC_STAT,product(shape(ARR))*kind(ARR),QUOTEME(ARR),ABI_FUNC)
#  define ABI_DEALLOCATE(ARR) \
   ABI_ALLOC_SIZE=-product(shape(ARR))*kind(ARR) NEWLINE \
   deallocate(ARR,stat=ABI_ALLOC_STAT) NEWLINE \
   call memocc(ABI_ALLOC_STAT,ABI_ALLOC_SIZE,QUOTEME(ARR),ABI_FUNC)
# else
#  define ABI_ALLOCATE(ARR,SIZE) \
   allocate(ARR SIZE,stat=ABI_ALLOC_STAT) NEWLINE \
   call memocc(ABI_ALLOC_STAT,product(shape(ARR))*kind(ARR),QUOTEME(ARR),ABI_FUNC)
#  define ABI_DEALLOCATE(ARR) \
   ABI_ALLOC_SIZE=-product(shape(ARR))*kind(ARR) NEWLINE \
   deallocate(ARR,stat=ABI_ALLOC_STAT) NEWLINE \
   call memocc(ABI_ALLOC_STAT,ABI_ALLOC_SIZE,QUOTEME(ARR),ABI_FUNC)
# endif
#else
# define ABI_ALLOCATE(ARR,SIZE) \
  allocate(ARR SIZE,stat=ABI_ALLOC_STAT)
# define ABI_DEALLOCATE(ARR) \
  deallocate(ARR,stat=ABI_ALLOC_STAT)
#endif

#define ABI_POINTER_ALLOCATE(ARR) \
   allocate(ARR)
#define ABI_POINTER_DEALLOCATE(ARR) \
   deallocate(ARR)

#define MEM_ALLOCATE(array, size) \
 allocate(array, stat=ERR_ALLOCATE_STAT) NEWLINE \
 if (ERR_ALLOCATE_STAT /= 0) then        NEWLINE \
  call memerr( QUOTEME(array),int(size*sizeof(kind(array)),8), __FILE__,__LINE__) NEWLINE \
 end if

#define MEM_FREE(array) \
 deallocate(array, stat=ERR_ALLOCATE_STAT) NEWLINE \
 if (ERR_ALLOCATE_STAT /= 0) then        NEWLINE \
  ABI_DIE("deallocation error") NEWLINE \
 end if

/**
* Variadic arguments are supported by several compilers,
* however they have been standardized only in C99.
#define HAVE_VARMACROS
**/

#undef HAVE_VARMACROS

#ifdef HAVE_VARMACROS  /**|| (__STDC_VERSION__ >= 199901L)**/

#define MEM_VARALLOCATE(...) allocate(__VA_ARGS__,stat=ERR_ALLOCATE_STAT) NEWLINE \
 if (ERR_ALLOCATE_STAT /= 0) then NEWLINE \
  ABI_DIE("allocation error") NEWLINE \
 end if

#define MEM_VARFREE(...) deallocate(__VA_ARGS__,stat=ERR_ALLOCATE_STAT) NEWLINE \
 if (ERR_ALLOCATE_STAT /= 0) then NEWLINE \
  ABI_DIE("allocation error") NEWLINE \
 end if

#else

#define MEM_VARALLOCATE allocate
#define MEM_VARFREE deallocate

#endif /** defined(HAVE_VARMACROS) || (__STDC_VERSION__ >= 199901L) **/

#ifdef DEBUG_MODE
#  ifdef HAVE_FC_LONG_LINES
#    define DBG_CHECK(expr,str) call assert((expr), str,__FILE__,__LINE__)
#    define DBG_CHKPT(value) write(std_out,*)__FILE__,":",__LINE__,":",value
/*
* C99 __func__ doesn't work with Fortran code but we might use abilint to define a
* CPP varible, ABI_func, containing the name of the F90 procedure.
#    define DBG_ENTER(mode) call sentinel(1,mode,ABI_func)
#    define DBG_EXIT(mode)  call sentinel(2,mode,ABI_func)
* For the moment we use __FILE__
*/
#    define DBG_ENTER(mode) call sentinel(1,mode,__FILE__,__LINE__)
#    define DBG_EXIT(mode)  call sentinel(2,mode,__FILE__,__LINE__)
#  else

#    define DBG_CHECK(expr,str) call assert((expr),str)
#    define DBG_CHKPT(value) write(std_out,*)value
#    define DBG_ENTER(mode) call sentinel(1,mode)
#    define DBG_EXIT(mode)  call sentinel(2,mode)
#  endif

#ifdef HAVE_VARMACROS  /**|| (__STDC_VERSION__ >= 199901L)**/
#  define DBG_WRITE(...)  write(std_out,*) __VA_ARGS__
#else
  /* FIXME Cannot use this trick with intel since fpp treats ! as the beginning of a comment.
   * defining a new macro works with fpp but gfortran complains about the number of arguments. */
#  define DBG_WRITE !variadic macros not supported
#endif

#else
#  define DBG_CHECK(expr,str)
#  define DBG_CHKPT(value)
#  define DBG_ENTER(mode)
#  define DBG_EXIT(mode)
#  define DBG_WRITE !stripped debugging write
#endif

/* Macro for basic messages (COLLECTIVE and PERSONAL version) */
#ifdef HAVE_FC_LONG_LINES

#  define MSG_COMMENT(msg)      call msg_hndl(msg,"COMMENT","COLL",__FILE__,__LINE__)
#  define MSG_WARNING(msg)      call msg_hndl(msg,"WARNING","COLL",__FILE__,__LINE__)
#  define MSG_ERROR(msg)        call msg_hndl(msg,"ERROR"  ,"COLL",__FILE__,__LINE__)
#  define MSG_BUG(msg)          call msg_hndl(msg,"BUG"    ,"COLL",__FILE__,__LINE__)
#  define MSG_PERS_COMMENT(msg) call msg_hndl(msg,"COMMENT","PERS",__FILE__,__LINE__)
#  define MSG_PERS_WARNING(msg) call msg_hndl(msg,"WARNING","PERS",__FILE__,__LINE__)
#  define MSG_PERS_ERROR(msg)   call msg_hndl(msg,"ERROR"  ,"PERS",__FILE__,__LINE__)
#  define MSG_PERS_BUG(msg)     call msg_hndl(msg,"BUG"    ,"PERS",__FILE__,__LINE__)

#  define ETSF_CHECK_ERROR(lstat,Error_data)   call abietsf_msg_hndl(lstat,Error_data,"COLL",__FILE__,__LINE__)
#  define ETSF_CHECK_MYERROR(lstat,Error_data) call abietsf_msg_hndl(lstat,Error_data,"PERS",__FILE__,__LINE__)
#  define ETSF_WARN(lstat,Error_data) call abietsf_warn(lstat,Error_data,"COLL",__FILE__,__LINE__)

#else
/*
 * Safe macros for emergency cases!
 * Useful if __FILE__ expands to the full path name exceeding
 * the max number of Fortran columns. ISO doesn't define any standard!
 */
#  define MSG_COMMENT(msg)      call msg_hndl(msg,"COMMENT","COLL")
#  define MSG_WARNING(msg)      call msg_hndl(msg,"WARNING","COLL")
#  define MSG_ERROR(msg)        call msg_hndl(msg,"ERROR"  ,"COLL")
#  define MSG_BUG(msg)          call msg_hndl(msg,"BUG"    ,"COLL")
#  define MSG_PERS_COMMENT(msg) call msg_hndl(msg,"COMMENT","PERS")
#  define MSG_PERS_WARNING(msg) call msg_hndl(msg,"WARNING","PERS")
#  define MSG_PERS_ERROR(msg)   call msg_hndl(msg,"ERROR"  ,"PERS")
#  define MSG_PERS_BUG(msg)     call msg_hndl(msg,"BUG"    ,"PERS")

#  define ETSF_CHECK_ERROR(lstat,Error_data)   call abietsf_msg_hndl(lstat,Error_data,"COLL")
#  define ETSF_CHECK_MYERROR(lstat,Error_data) call abietsf_msg_hndl(lstat,Error_data,"PERS")
#  define ETSF_WARN(lstat,Error_data) call abietsf_warn(lstat,Error_data,"COLL")

#endif

/* Does the compiler support => in declarations? */
#ifdef HAVE_FC_NULL
#  define SET2NULL => null()
#else
#  define SET2NULL
#endif


/* Does the compiler support allocatable arrays in datatypes? */
#ifdef HAVE_FC_ALLOCATABLE_DTARRAYS
#  define DTA_ALLOCATABLE_TG  allocatable,target
#  define DTA_ALLOCATABLE     allocatable
#  define DTA_SFREE(arr)      if (allocated(arr)) deallocate(arr)
#  define DTA_NULLIFY(arr)

#else
/* Have to use pointers in datatypes instead of allocatable arrays */
#  define DTA_ALLOCATABLE_TG  pointer
#  define DTA_ALLOCATABLE     pointer
#  define DTA_SFREE(arr)      if (associated(arr)) deallocate(arr)
#  define DTA_NULLIFY(arr)    nullify(arr)
#endif

/* Dummy use of unused arguments to silence compiler warnings */
#define ABI_UNUSED(var) if (.FALSE.) call unused_var(var)

/* #define HAVE_TIMER_ABINIT */
/* #undefine HAVE_TIMER_ABINIT  */

#define ABI_FUNC " "

#if defined HAVE_TIMER_ABINIT && 0
#  define ABI_TIMER_START(key) call timing(1,ABI_FUNC//":"//key)
#  define ABI_TIMER_STOP(key)  call timing(2,ABI_FUNC//":"//key)
#ifdef DEBUG_MODE 
#  define DEV_TIMER_START(key) call timing(1,ABI_FUNC//":"//key)
#  define DEV_TIMER_STOP(key)  call timing(2,ABI_FUNC//":"//key)
#endif

#else
#  define ABI_TIMER_START(key) 
#  define ABI_TIMER_STOP(key) 
#  define DEV_TIMER_START(key) 
#  define DEV_TIMER_STOP(key) 
#endif

/* NetCDF macros */
#define NETCDF_CHECK(nc_call) \
 ncerr = nc_call NEWLINE \
 if ( ncerr /= NF90_NOERR ) then NEWLINE \
  call netcdf_ioerr(ncerr,__FILE__,__LINE__) NEWLINE \
 end if

/* Portable support for /dev/null */
#if defined HAVE_OS_WINDOWS
#define NULL_FILE "NUL"
#else
#define NULL_FILE "/dev/null"
#endif


#endif 
/* _ABINIT_COMMON_H */
