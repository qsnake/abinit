!{\src2tex{textfont=tt}}
!!****f* ABINIT/leave_myproc
!! NAME
!! leave_myproc
!!
!! FUNCTION
!! Routine for clean exit of f90 code by one processor
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, NCJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   option=(optional argument, default=0)
!!          0: closes output file, then stops
!!          1: only closes output file
!!          2: only stops
!!
!! OUTPUT
!!  (only writing, then stop)
!!
!! NOTES
!!  By default, it uses "call exit(1)", that is not completely portable.
!!
!! PARENTS
!!      leave_new,leave_test
!!
!! CHILDREN
!!      exit
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine leave_myproc(option)

 use m_profiling

 use defs_basis
#if defined FC_NAG
 use f90_unix
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'leave_myproc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in),optional :: option

!Local variables-------------------------------
!scalars
 integer :: option_
 logical :: testopen

! **********************************************************************

 option_=0;if(present(option)) option_=option

 if (option_==0.or.option_==1) then
   inquire(ab_out,OPENED=testopen)
   if (testopen) close(ab_out)
 end if

 if (option_==0.or.option_==2) then
#if defined FC_NAG
   call exit(-1)
#elif defined HAVE_FC_EXIT
   call exit(1)
#else
   stop 1
#endif
 end if

end subroutine leave_myproc
!!***
