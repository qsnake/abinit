!{\src2tex{textfont=tt}}
!!****f* ABINIT/wrtout_myproc
!! NAME
!!  wrtout_myproc
!!
!! FUNCTION
!!  Do the output for one proc. For parallel or sequential output use wrtout()
!!  instead. Also allows to treat correctly the write operations for Unix (+DOS) and MacOS.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  unit=unit number for writing
!!  message=(character(len=*)) message to be written
!!
!!  mpi_comm= Optional argument
!!            If present, no printing is done
!!            Variables iexit, nwarning and ncomment are
!!            summed over the mpi_comm communicator
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      gstateimg,wrtout
!!
!! CHILDREN
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wrtout_myproc(unit,message,&
&                        mpi_comm) ! optional argument

 use m_profiling

 use defs_basis
 use m_xmpi, only : xsum_mpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrtout_myproc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 character(len=*),intent(in) :: message
 integer,intent(in),optional :: mpi_comm

!Local variables-------------------------------
!scalars
 integer,save :: iexit=0,ncomment=0,nwarning=0
 integer :: ierr,lenmessage,rtnpos
 logical :: print_std_err
 character(len=len(message)) :: messtmp
!arrays
 integer, allocatable :: buf(:)

!******************************************************************

!When I/O are redirected, it is sometimes necessary to reduce counters (saved) values;
!this can be done by passing mpi_comm optional argument to the routine
!In that case, no printing is done.
 if (present(mpi_comm)) then
   ABI_ALLOCATE(buf,(3))
   buf(1)=iexit;buf(2)=ncomment;buf(3)=nwarning
   call xsum_mpi(buf,mpi_comm,ierr)
   iexit=buf(1);ncomment=buf(2);nwarning=buf(3)
   if (iexit/=0) iexit=1
   ABI_DEALLOCATE(buf)
   return
 end if

 print_std_err=(unit==std_out.and.&
& (index(trim(message),'BUG')/=0.or.index(trim(message),'ERROR')/=0))

 if(message/=' ') then
   messtmp=message
   lenmessage=len(message)
!  Here, split the message, according to the char(10)
!  characters (carriage return). This technique is portable accross different OS.
   rtnpos=index(messtmp,ch10)
!  do while(rtnpos/=0)
   do while(rtnpos > 1)
     write(unit, '(a)' ) trim(messtmp(1:rtnpos-1))
     if (print_std_err) write(std_err, '(a)' ) trim(messtmp(1:rtnpos-1))
     messtmp=messtmp(rtnpos+1:lenmessage)
     lenmessage=lenmessage-rtnpos
     rtnpos=index(messtmp,ch10)
   end do
   write(unit, '(a)' ) trim(messtmp)
   if (print_std_err) write(std_err, '(a)' ) trim(messtmp)
 else
   write(unit,*)
   if (print_std_err) write(std_err,*)
 end if

 if( index(trim(message),'BUG') /= 0 )then
   write(unit, '(a)' ) '  Action : contact ABINIT group.'
   if (print_std_err) write(std_err, '(a)' ) '  Action : contact ABINIT group.'
   write(unit,*)
   if (print_std_err) write(std_err,*)
 end if

 if( index(trim(message),'BUG') /= 0   .or. &
& index(trim(message),'Calculation completed') /= 0 )then
   if(nwarning<10000 .and. ncomment<1000)then
     write(unit, '(a,i5,a,i4,a)' ) &
&     '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
   else
     write(unit, '(a,i6,a,i6,a)' ) &
&     '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
   end if
   if(iexit/=0)then
     write(unit, '(a)' ) ' Note : exit requested by the user.'
   end if
 end if

 if( index(trim(message),'Exit') /= 0 )then
   iexit=1
 end if

!Count the number of warnings and comments. Only take into
!account unit 6, in order not to duplicate these numbers.
 if( index(trim(message),'WARNING') /= 0 .and. unit==std_out )then
   nwarning=nwarning+1
 end if
 if( index(trim(message),'COMMENT') /= 0 .and. unit==std_out )then
   ncomment=ncomment+1
 end if

end subroutine wrtout_myproc
!!***
