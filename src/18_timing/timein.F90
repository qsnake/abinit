!{\src2tex{textfont=tt}}
!!****f* ABINIT/timein
!! NAME
!!  timein
!!
!! FUNCTION
!!  Timing routine. Returns cpu and wall clock time in seconds since some arbitrary start.
!!  For wall clock time, call the F90 intrinsic date_and_time .
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, LSI, MM, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (no inputs)
!!
!! OUTPUT
!!  cpu= cpu time in seconds
!!  wall= wall clock time in seconds
!!
!! NOTES
!!  For CPU time, contains machine-dependent code (choice will be selected
!!  by C preprocessor, see abi_cpu_time).
!!
!! PARENTS
!!      abinit,aim,aim_follow,anaddb,chkexi,cpdrv,drvaim,elphon,first_rec
!!      m_fft_prof,m_fftw3,m_shexc,m_shirley,m_timer,mkifc9,mkphbs,pclock
!!      rdddb9,rsiaf9,rsurf,surf,thm9,timab
!!
!! CHILDREN
!!      date_and_time,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine timein(cpu,wall)

 use m_profiling

 use defs_basis

#ifdef HAVE_FC_ISO_C_BINDING
 use iso_c_binding
#else
 use m_iso_c_binding
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'timein'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing, except_this_one => timein
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: cpu,wall

!Local variables-------------------------------
!scalars
 integer, parameter :: nday(24)=(/31,28,31,30,31,30,31,31,30,31,30,31,&
&                                 31,28,31,30,31,30,31,31,30,31,30,31/)
 integer, save :: month_init,month_now,start=1,year_init
 integer :: months
 character(len=8)   :: date
 character(len=10)  :: time
 character(len=5)   :: zone
 character(len=500) :: message
!arrays
 integer :: values(8)

! *************************************************************************

!CPU time _______________________________________

!It is possible to suppress the call to an external routine, and leave
!the cpu time to 0.0d0, provided the timing of the timer is suppressed in timana.f
!(simply set the loop counter maximum value to 1 in that routine)
 cpu = abi_cpu_time()

!Wallclock time ________________________________

!The following section of code is standard F90, but it is useful only if the intrinsics
!date_and_time is accurate at the 0.01 sec level, which is not the case for a P6 with the pghpf compiler ...
!Year and month initialisation
 if(start==1)then
   start=0
   call date_and_time(date,time,zone,values)
   year_init=values(1)
   month_init=values(2)
 end if

!write(std_out,*)' timein : before date_and_time '

!Uses intrinsic F90 subroutine Date_and_time for
!wall clock (not correct when a change of year happen)
 call date_and_time(date,time,zone,values)

!Compute first the number of seconds from the beginning of the month
 wall=(values(3)*24.0d0+values(5))*3600.0d0+values(6)*60.0d0+values(7)+values(8)*0.001d0

!If the month has changed, compute the number of seconds
!to be added. This fails if the program ran one year !!
 month_now=values(2)
 if(month_now/=month_init)then
   if(year_init+1==values(1))then
     month_now=month_now+12
   end if
   if(month_now<=month_init)then
     write(message, '(a,a,a,a)' ) ch10,&
&     ' timein : BUG -',ch10,&
&     '  Problem with month and year numbers.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   do months=month_init,month_now-1
     wall=wall+86400.0d0*nday(months)
   end do
 end if

!Now take into account bissextile years (I think 2000 is bissextile, but I am not sure ...)
 if(mod(year_init,4)==0 .and. month_init<=2 .and. month_now>2)   wall=wall+3600.0d0
 if(mod(values(1),4)==0 .and. month_init<=14 .and. month_now>14) wall=wall+3600.0d0

!write(std_out,*)' timein : wall at exit ',wall

end subroutine timein
!!***
