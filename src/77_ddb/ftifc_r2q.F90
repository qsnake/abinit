!{\src2tex{textfont=tt}}
!!****f* ABINIT/ftifc_r2q
!!
!! NAME
!! ftifc_r2q
!!
!! FUNCTION
!! (r->q):
!!   Generates the Fourier transform of the interatomic forces
!!   to obtain dynamical matrices (reciprocal space).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (JCC,XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! atmfrc(2,3,natom,3,natom,nrpt)
!!  = Interatomic Forces in real space !!
!!  We used the imaginary part just for debugging !
!! gprim(3,3)= Normalized coordinates in reciprocal space
!! natom= Number of atoms in the unit cell
!! nqpt= Number of q points in the Brillouin zone
!! nrpt= Number of R points in the Big Box
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!           These coordinates are normalized (=> * acell(3)!!)
!! spqpt(3,nqpt)= Reduced coordinates of the q vectors in reciprocal space
!! wghatm(natom,natom,nrpt)
!!         = Weights associated to a pair of atoms and to a R vector
!!
!! OUTPUT
!! dynmat(2,3,natom,3,natom,nqpt)
!!  = Dynamical matrices coming from the Derivative Data Base
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      gtdyn9
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ftifc_r2q(atmfrc,dynmat,gprim,natom,nqpt,&
&                  nrpt,rpt,spqpt,wghatm)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ftifc_r2q'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nqpt,nrpt
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),spqpt(3,nqpt)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt)
 real(dp),intent(out) :: dynmat(2,3,natom,3,natom,nqpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,iqpt,irpt,mu,nu
 real(dp) :: facti,factr,im,kr,re
!arrays
 real(dp) :: kk(3)

! *********************************************************************

 dynmat(:,:,:,:,:,:)=0.0_dp
 do iqpt=1,nqpt
   do irpt=1,nrpt

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

!    Get phase factor
     re=cos(two_pi*kr)
     im=sin(two_pi*kr)

!    Inner loop on atoms and directions
     do ib=1,natom
       do ia=1,natom
         if(abs(wghatm(ia,ib,irpt))>1.0d-10)then
           factr=re*wghatm(ia,ib,irpt)
           facti=im*wghatm(ia,ib,irpt)
           do nu=1,3
             do mu=1,3
!              Real and imaginary part of the dynamical matrices
               dynmat(1,mu,ia,nu,ib,iqpt)=dynmat(1,mu,ia,nu,ib,iqpt)&
&               +factr*atmfrc(1,mu,ia,nu,ib,irpt)
!              Atmfrc should be real
!              &       -im*wghatm(ia,ib,irpt)*atmfrc(2,mu,ia,nu,ib,irpt)
               dynmat(2,mu,ia,nu,ib,iqpt)=dynmat(2,mu,ia,nu,ib,iqpt)&
&               +facti*atmfrc(1,mu,ia,nu,ib,irpt)
!              Atmfrc should be real
!              &        +re*wghatm(ia,ib,irpt)*atmfrc(2,mu,ia,nu,ib,irpt)
             end do
           end do
         end if
       end do
     end do
   end do
 end do

end subroutine ftifc_r2q
!!***
