!{\src2tex{textfont=tt}}
!!****f* ABINIT/ftifc_q2r
!!
!! NAME
!! ftifc_q2r
!!
!! FUNCTION
!!   Generates the Fourier transform of the dynamical matrices
!!   to obtain interatomic forces (real space).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (JCC,XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dynmat(2,3,natom,3,natom,nqpt)
!!  = Dynamical matrices coming from the Derivative Data Base
!! gprim(3,3)= Normalized coordinates in reciprocal space
!! natom= Number of atoms in the unit cell
!! nqpt= Number of q points in the Brillouin zone
!! nrpt= Number of R points in the Big Box
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!           These coordinates are normalized (=> * acell(3)!!)
!! spqpt(3,nqpt)= Reduced coordinates of the q vectors in reciprocal space
!!
!! OUTPUT
!! atmfrc(2,3,natom,3,natom,nrpt)
!!  = Interatomic Forces in real space !!
!!  We used the imaginary part just for debugging !
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      mkifc9,scphon_interpolate_phonon_and_dos
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ftifc_q2r(atmfrc,dynmat,gprim,natom,nqpt,&
&                  nrpt,rpt,spqpt)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ftifc_q2r'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nqpt,nrpt
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),spqpt(3,nqpt)
 real(dp),intent(out) :: atmfrc(2,3,natom,3,natom,nrpt)
 real(dp),intent(in) :: dynmat(2,3,natom,3,natom,nqpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,iqpt,irpt,mu,nu
 real(dp) :: im,kr,re
!arrays
 real(dp) :: kk(3)

! *********************************************************************

!Interatomic Forces from Dynamical Matrices
 atmfrc(:,:,:,:,:,:)=0.0_dp
 do irpt=1,nrpt
   do iqpt=1,nqpt

!    Calculation of the k coordinates in Normalized Reciprocal
!    coordinates
     kk(1)=spqpt(1,iqpt)*gprim(1,1)+spqpt(2,iqpt)*&
&     gprim(1,2)+spqpt(3,iqpt)*gprim(1,3)
     kk(2)=spqpt(1,iqpt)*gprim(2,1)+spqpt(2,iqpt)*&
&     gprim(2,2)+spqpt(3,iqpt)*gprim(2,3)
     kk(3)=spqpt(1,iqpt)*gprim(3,1)+spqpt(2,iqpt)*&
&     gprim(3,2)+spqpt(3,iqpt)*gprim(3,3)

!    Product of k and r
     kr=kk(1)*rpt(1,irpt)+kk(2)*rpt(2,irpt)+kk(3)*&
&     rpt(3,irpt)

!    Get the phase factor
     re=cos(two_pi*kr)
     im=sin(two_pi*kr)

!    Now, big inner loops on atoms and directions
!    The indices are ordered to give better speed
     do ib=1,natom
       do nu=1,3
         do ia=1,natom
           do mu=1,3
!            Real and imaginary part of the interatomic forces
             atmfrc(1,mu,ia,nu,ib,irpt)=atmfrc(1,mu,ia,nu,ib,irpt)&
&             +re*dynmat(1,mu,ia,nu,ib,iqpt)&
&             +im*dynmat(2,mu,ia,nu,ib,iqpt)
!            The imaginary part should be equal to zero !!!!!!
!            atmfrc(2,mu,ia,nu,ib,irpt)=atmfrc(2,mu,ia,nu,ib,irpt)
!            &          +re*dynmat(2,mu,ia,nu,ib,iqpt)
!            &          -im*dynmat(1,mu,ia,nu,ib,iqpt)
           end do
         end do
       end do
     end do

   end do
 end do

!The sum has to be weighted by a normalization factor of 1/nqpt
 atmfrc(:,:,:,:,:,:)=atmfrc(:,:,:,:,:,:)/nqpt


end subroutine ftifc_q2r
!!***
