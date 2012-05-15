!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_cpu_time
!! NAME
!!  abi_cpu_time
!!
!! FUNCTION
!!  Timing routine. Returns cpu time in seconds since some arbitrary start.
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
!!  cpu_time= cpu time in seconds
!!
!! NOTES
!!  For CPU time, contains machine-dependent code (choice will be selected by c preprocessor).
!!  Note that all supported machines are listed explicitly below; there
!!  is no "else" which covers "other".  The C preprocessor will place
!!  a spurious line of code (see below) into the fortran source unless
!!  preprocessed with -Dflag where flag refers to one of the supported machines.
!!
!!  WARNING: the following list is no more accurate (YP 20060530)
!!
!!  Presently supported flags: "ibm", "hp", "P6", "dec_alpha", "sgi", "vpp", "sun", "mac", "nec", "sr8k"
!!  Previously supported flags:  "ultrix". Might still work !
!!
!!  Calls machine-dependent "mclock" for "ibm" .
!!  Calls ANSI C subroutine "cclock" for "hp" and "sgi".
!!  Calls machine-dependent "etime" for "P6", "mac", "dec_alpha", "sun", "nec" .
!!  Calls machine-dependent "clock" for "vpp"
!!  Calls machine-dependent "xclock" for "sr8k"
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


function abi_cpu_time() result(cpu)

 use defs_basis

#ifdef HAVE_FC_ISO_C_BINDING
 use iso_c_binding
#else
 use m_iso_c_binding
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_cpu_time'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 real(dp) :: cpu

!Local variables-------------------------------
#ifdef HAVE_FC_CPUTIME
 real :: cpu_sp
#elif defined FC_IBM
 integer :: mclock
#elif defined FC_SUN
 real :: tmp(2)
 real :: etime
#elif defined FC_COMPAQ || defined HAVE_OS_MACOSX 
 real :: tmp(2)           !real array only needed by etime
 real(dp) :: etime
#else
 integer :: count_now,count_max,count_rate
#endif

! *************************************************************************

!Machine-dependent timers
#ifdef HAVE_CCLOCK
 call cclock(cpu)

#elif defined HAVE_FC_CPUTIME
!This is the F95 standard subroutine.
 call cpu_time(cpu_sp)
 cpu = cpu_sp

#elif defined FC_IBM
 cpu = mclock()*0.01d0

#elif defined HAVE_OS_MACOSX || defined FC_COMPAQ || defined FC_SUN
 cpu = etime(tmp)

#elif defined FC_FUJITSU
 call clock(cpu,0,2)

#elif defined FC_HITACHI
 call xclock(cpu,5)

#else
!This is the Fortran90 standard subroutine, might not always be sufficiently accurate
 call system_clock(count_now,count_rate,count_max)
 cpu=dble(count_now)/dble(count_rate)
#endif

end function abi_cpu_time
!!***
