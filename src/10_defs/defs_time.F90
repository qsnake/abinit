!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_time
!! NAME
!! defs_time
!!
!! FUNCTION
!! This module contains accumulators for the timer.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!! Include the name of all routines: better modularity
!!
!! PARENTS
!!    timab,time_accu
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_time

 use defs_basis

#ifdef HAVE_FC_ISO_C_BINDING
 use iso_c_binding
#else
 use m_iso_c_binding
#endif

 implicit none

!mtim determines the maximum number of "timing slots" available
 integer,parameter :: mtim=999

! timeopt is a flag which indicates the suppression or not of the timing.
 integer :: timopt=1

! papiopt is a flag which indicates if there is or not an analysis of speed execution is made. 
! By defaut the analysis is not done 
 integer :: papiopt=0

! initpapiopt is a flag which permits initialisation of overall timing by papi and speed of run
 integer :: initpapiopt=1 

! Number of times that the routine has been called
 integer :: ncount(mtim)

! Store the values 
 real(dp) :: cpu,wall
 integer(C_LONG_LONG) :: flops1
 real(C_FLOAT) :: real_time, proc_time

! Accumulating cpu time (1) and wall to wall time (2) for each "timing slots"
 real(dp)  :: acctim(2,mtim),tzero(2,mtim)

! Accumulating number of floating point operation and cpu time (1) and wall to wall times (2) for each "performance slot"
 real(dp) :: papi_accflops(mtim), papi_acctim(2,mtim)

! Reference value for number of floating point operation and time ( cpu and wall) for each performance slot
 real(dp) :: flops(mtim) , papi_tzero(2,mtim)

! Elapsed time and elapsed number of floating point operation since a reference
 real(dp) :: papi_tottim(2,mtim), papi_totflops(mtim)

end module defs_time
!!***
