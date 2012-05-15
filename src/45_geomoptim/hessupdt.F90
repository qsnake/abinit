!{\src2tex{textfont=tt}}
!!****f* ABINIT/hessupdt
!! NAME
!! hessupdt
!!
!! FUNCTION
!! Update of the hessian matrix according to the Broyden formula.
!! Could see Numerical Recipes (Fortran), 1986, page 307.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG, JCC).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  iatfix(3,natom)=1 for each atom fixed along specified direction, else 0
!!  natom=number of atoms in unit cell
!!  ndim=size of the hessian and vectors
!!  vin(ndim)=new input vector
!!  vin_prev(ndim)=previous input vector
!!  vout(ndim)=new output vector
!!  vout_prev(ndim)=previous output vector
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  hessin(ndim,ndim)=hessian matrix, updated at output.
!!
!! PARENTS
!!      pred_bfgs,pred_delocint,pred_diisrelax,xfh_recover,xfh_recover_deloc
!!      xfh_recover_new
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine hessupdt(hessin,iatfix,natom,ndim,vin,vin_prev,vout,vout_prev)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hessupdt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ndim
!arrays
 integer,intent(in) :: iatfix(3,natom)
 real(dp),intent(in) :: vin(ndim),vin_prev(ndim),vout(ndim),vout_prev(ndim)
 real(dp),intent(inout) :: hessin(ndim,ndim)

!Local variables-------------------------------
!scalars
 integer :: iatom,idir,ii,jj
 real(dp) :: den1,den2,den3
!arrays
 real(dp) :: bfgs(ndim),din(ndim),dout(ndim),hdelta(ndim)

!***************************************************************************


!write(ab_out,*) 'VECTOR INPUT (vin)'
!do ii=1,ndim,3
!if (ii+2<=ndim)then
!write(ab_out,*) ii,vin(ii:ii+2)
!else
!write(ab_out,*) ii,vin(ii:ndim)
!end if
!end do
!write(ab_out,*) 'VECTOR OUTPUT (vout)'
!do ii=1,ndim,3
!if (ii+2<=ndim)then
!write(ab_out,*) ii,vout(ii:ii+2)
!else
!write(ab_out,*) ii,vout(ii:ndim)
!end if
!end do

!write(ab_out,*) 'VECTOR INPUT (vin_prev)'
!do ii=1,ndim,3
!if (ii+2<=ndim)then
!write(ab_out,*) ii,vin(ii:ii+2)
!else
!write(ab_out,*) ii,vin(ii:ndim)
!end if
!end do
!write(ab_out,*) 'VECTOR OUTPUT (vout_prev)'
!do ii=1,ndim,3
!if (ii+2<=ndim)then
!write(ab_out,*) ii,vout(ii:ii+2)
!else
!write(ab_out,*) ii,vout(ii:ndim)
!end if
!end do

!Difference between new and previous vectors
!Implement fixing of atoms; must discard the change of forces on fixed atoms
 din(:) =vin(:) -vin_prev(:)
 dout(:)=vout(:)-vout_prev(:)
 do iatom=1,natom
   do idir=1,3
     if (iatfix(idir,iatom) == 1) then
       dout(idir+(iatom-1)*3)=0.0_dp
     end if
   end do
 end do

!Compute approximate inverse Hessian times delta fcart
!hdelta=hessin*deltaf
 hdelta(:)=0.0_dp
 do ii=1,ndim
   hdelta(:)=hdelta(:)+hessin(:,ii)*dout(ii)
 end do

!Calculation of dot products for the denominators
 den1=0.0_dp ; den2=0.0_dp
 do ii=1,ndim
   den1=den1+dout(ii)*din(ii)
   den2=den2+dout(ii)*hdelta(ii)
 end do

!DEBUG
!write(std_out,*)' hessupdt : den1,den2',den1,den2
!write(std_out,*)' din ',din
!write(std_out,*)' dout ',dout
!write(std_out,*)' hdelta ',hdelta
!ENDDEBUG

!Denominators are multiplicative
 den1=1.0_dp/den1
 den3=1.0_dp/den2

!Vectors which make a difference between the BROYDEN and
!the DAVIDON scheme.
 bfgs(:)=den1*din(:)-den3*hdelta(:)

!B.F.G.S. updating formula
 do ii=1,ndim
   do jj=1,ndim
     hessin(ii,jj)=hessin(ii,jj) +den1*din(ii)*din(jj) &
&     -den3*hdelta(ii)*hdelta(jj) +den2*bfgs(ii)*bfgs(jj)
   end do
 end do

end subroutine hessupdt
!!***
