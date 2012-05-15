!{\src2tex{textfont=tt}}
!!****f* ABINIT/q0dy3_calc
!! NAME
!! q0dy3_calc
!!
!! FUNCTION
!! Calculate the q=0 correction term to the dynamical matrix
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG, MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dyew(2,3,natom,3,natom)= dynamical matrix
!!    input, non-corrected, for q=0 if option=1 or 2
!!  natom= number of atom in the unit cell
!!  option= either 1 or 2:
!!     1: use dyew to calculate dyewq0 symmetrical form
!!     2: use dyew to calculate dyewq0 symmetrical form
!!
!! OUTPUT
!!  dyewq0(3,3,natom) = part needed to correct
!!    the dynamical matrix for atom self-interaction.
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Should be used just after each call to ewald3, for both
!! q==0 and the real wavelength.
!!
!! If option=1 or 2, q0dy3_calc uses an Ewald dynamical matrix at q=0,
!! called dyew, to produce a contracted form called dyewq0 :
!! either:
!!   in an unsymmetrical form (if option=1), or
!!   in a symmetrical form (if option=2).
!!
!! The q0dy3_calc should be used in conjunction with the subroutine
!! ewald3 (or ewald9).
!! First, the call of ewald3 with q==0 should be done ,
!!   then the call to q0dy3_calc will produce
!!   the dyewq0 matrix from the (q=0) dyew matrix
!! Second, the call of ewald3 with the real q (either =0 or diff 0)
!!   should be done, then the call to q0dy3_apply
!!   will produce the correct dynamical matrix dyew starting from
!!   the previously calculated dyewq0 and the bare(non-corrected)
!!   dyew matrix
!!
!! PARENTS
!!      hybrid9,mkifc9,respfn
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine q0dy3_calc(natom,dyewq0,dyew,option)

 use m_profiling

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'q0dy3_calc'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,option
!arrays
 real(dp),intent(in) :: dyew(2,3,natom,3,natom)
 real(dp),intent(out) :: dyewq0(3,3,natom)

!Local variables -------------------------
!scalars
 integer :: ia,ib,mu,nu
 character(len=500) :: message

! *********************************************************************

 if(option==1.or.option==2)then
   do mu=1,3
     do nu=1,3
       do ia=1,natom
         dyewq0(mu,nu,ia)=0.0_dp
         do ib=1,natom
           dyewq0(mu,nu,ia)=dyewq0(mu,nu,ia)+dyew(1,mu,ia,nu,ib)
         end do
       end do
     end do
   end do
 else
   write (message, '(3a)') ' q0dy3_calc : ERROR option should be 1 or 2.',ch10,&
&   ' action: correct calling routine'
   MSG_BUG(message)
 end if

 if(option==2)then
   do ia=1,natom
     do mu=1,3
       do nu=mu,3
         dyewq0(mu,nu,ia)=(dyewq0(mu,nu,ia)+dyewq0(nu,mu,ia))/2
         dyewq0(nu,mu,ia)=dyewq0(mu,nu,ia)
       end do
     end do
   end do
 end if

end subroutine q0dy3_calc
!!***
