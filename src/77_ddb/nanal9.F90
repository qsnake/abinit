!{\src2tex{textfont=tt}}
!!****f* ABINIT/nanal9
!!
!! NAME
!! nanal9
!!
!! FUNCTION
!! If plus=0 then substracts the non-analytical part from one dynamical
!!           matrices, with number iqpt.
!! If plus=1 then adds the non-analytical part to the dynamical
!!           matrices, with number iqpt.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (JCC,XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dyew(2,3,natom,3,natom)= Non-analytical part
!! natom= Number of atoms in the unit cell
!! iqpt= Referenced q point for the dynamical matrix
!! nqpt= Number of q points
!! plus= (see above)
!!
!! OUTPUT
!! dynmat(2,3,natom,3,natom,nqpt)
!!  = Dynamical matrices coming from the Derivative Data Base
!!
!! PARENTS
!!      gtdyn9,mkifc9
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine nanal9(dyew,dynmat,iqpt,natom,nqpt,plus)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nanal9'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iqpt,natom,nqpt,plus
!arrays
 real(dp),intent(in) :: dyew(2,3,natom,3,natom)
 real(dp),intent(out) :: dynmat(2,3,natom,3,natom,nqpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,mu,nu
 character(len=500) :: message

! *********************************************************************

 if (plus==0) then

   do ia=1,natom
     do ib=1,natom
       do mu=1,3
         do nu=1,3
!          The following four lines are OK
           dynmat(1,mu,ia,nu,ib,iqpt)=dynmat(1,mu,ia,nu,ib,iqpt)&
&           -dyew(1,mu,ia,nu,ib)
           dynmat(2,mu,ia,nu,ib,iqpt)=dynmat(2,mu,ia,nu,ib,iqpt)&
&           -dyew(2,mu,ia,nu,ib)
!          DEBUG
!          dynmat(1,mu,ia,nu,ib,iqpt)=dynmat(1,mu,ia,nu,ib,iqpt)
!          dynmat(2,mu,ia,nu,ib,iqpt)=dynmat(2,mu,ia,nu,ib,iqpt)
!          ENDDEBUG

         end do
       end do
     end do
   end do

 else if (plus==1) then

!  DEBUG
!  write(std_out,*)' nanal9 : now, only the analytic part '
!  write(std_out,*)' nanal9 : now, only the ewald part '
!  ENDDEBUG
   do ia=1,natom
     do ib=1,natom
       do mu=1,3
         do nu=1,3
!          The following four lines arethe good ones
           dynmat(1,mu,ia,nu,ib,iqpt)=dynmat(1,mu,ia,nu,ib,iqpt)&
&           +dyew(1,mu,ia,nu,ib)
           dynmat(2,mu,ia,nu,ib,iqpt)=dynmat(2,mu,ia,nu,ib,iqpt)&
&           +dyew(2,mu,ia,nu,ib)
!          DEBUG
!          dynmat(1,mu,ia,nu,ib,iqpt)=dyew(1,mu,ia,nu,ib)
!          dynmat(2,mu,ia,nu,ib,iqpt)=dyew(2,mu,ia,nu,ib)
!          dynmat(1,mu,ia,nu,ib,iqpt)=dynmat(1,mu,ia,nu,ib,iqpt)
!          dynmat(2,mu,ia,nu,ib,iqpt)=dynmat(2,mu,ia,nu,ib,iqpt)
!          ENDDEBUG
         end do
       end do
     end do
   end do

 else

   write(message,'(a,a,a,a,a,i4,a)' )&
&   ' nanal9 : BUG -',ch10,&
&   '  The argument "plus" must be equal to 0 or 1.',ch10,&
&   '  The value ',plus,' is not available.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')

 end if

end subroutine nanal9
!!***
