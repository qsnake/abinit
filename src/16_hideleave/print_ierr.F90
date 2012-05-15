!{\src2tex{textfont=tt}}
!!****f* ABINIT/print_ierr
!! NAME
!!  print_ierr
!!
!! FUNCTION
!!  Routine for checking of the error code and clean exit in case this error code is non-zero.
!!  If there is a problem, the value of the error code will be printed, as well as the name of the calling routine,
!!  and the name of the routine that produced the error code.
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2012 ABINIT group (DCA, XG, GMR, NCJ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! ierr=error code from the called routine
!! nrcalled=the name of the routine that has been called
!! nrcalling=the name of the routine that is calling 
!!
!! OUTPUT
!!  (only writing, then stop)
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      leave_new,wrtout,xerror_string
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine print_ierr(ierr,nrcalled,nrcalling)

 use m_profiling

 use defs_basis
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_ierr'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave, except_this_one => print_ierr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: ierr
 character(len=*),intent(in) :: nrcalled,nrcalling

!Local variables-------------------------------
 integer :: ilen,ierror
 character(len=500) :: msg,err_string

! **********************************************************************

 if(ierr/=0)then
   write(msg,'(3a,i6,7a)')' print_ierr : ERROR -',ch10,&
&   ' There is a non-zero return code, ierr=',ierr,ch10,&
&   ' from the routine ',trim(nrcalled),ch10,&
&   ' called by the routine ',trim(nrcalling),'.'
   call wrtout(std_out,msg,"PERS")
   if(nrcalled(1:3)=="MPI")then
     call xerror_string(ierr,err_string,ilen,ierror)
     if(ierror==0)then
       write(std_out,'(3a)')&
&       ' This MPI error code has been associated by the routine MPI_Error_string,',ch10,&
&       ' to the following error message :'
       call wrtout(std_out,msg,"PERS")
       call wrtout(std_out,err_string,"PERS")
     end if
   end if
   call leave_new("PERS")
 end if

end subroutine print_ierr
!!***
