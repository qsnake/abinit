!{\src2tex{textfont=tt}}
!!****f* ABINIT/papi_init
!! NAME
!! papi_init
!!
!! FUNCTION
!! This function initializes papi high level interface , sets up counters
!! to monitor PAPI_FP_OPS and PAPI_TOT_CYC events and starts the counters
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      papif_flops,papif_library_init,papif_perror
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine papi_init()

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'papi_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!Local variables-------------------------------
#ifdef HAVE_TIMER_PAPI
#include "f90papi.h"
#define  C_FLOAT REAL
#define C_LONG_LONG INTEGER*8
#define C_INT INTEGER
character*(PAPI_MAX_STR_LEN) papi_errstr
C_INT :: retval
C_FLOAT :: unused1, unused2, unused4 
C_LONG_LONG :: unused3 
#endif

! *************************************************************************

#ifdef HAVE_TIMER_PAPI

 retval = PAPI_VER_CURRENT
 call PAPIf_library_init(retval)
 if ( retval.NE.PAPI_VER_CURRENT) then
   write(std_out,*) 'Problem library PAPI'
 end if

!First pass. Initializing counter


 call PAPIf_flops(unused1, unused2, unused3, unused4, retval)
 if (retval.NE.PAPI_OK) then
   write(std_out,*) 'Problem to initialize papi high level inteface'
   call papif_perror(retval,papi_errstr,retval)
   write(std_out,*) 'Error code', papi_errstr
 end if ! DEBUG


!call PAPIf_query_event(PAPI_FP_INS, retval) 
!if (retval .NE. PAPI_OK) then
!write(std_out,*) 'Problem query event'
!endif

!inializing papi high level interface , set up counters
!to monitor PAPI_FP_OPS and PAPI_TOT_CYC events and start the counters
!Subsequent calls will read the counters and return total
!real time, total process time, total floting point instructions
!or operations since the start of the mesurement and the Mflop/s rate
!since latests call to PAPI_flops
!call PAPIf_flops(real_time, proc_time, flpops, mflops, retval)
!if (retval.NE.PAPI_OK) then
!write(std_out,*) 'Problem to initialize papi high level inteface'
!endif

#endif

end subroutine papi_init
!!***
