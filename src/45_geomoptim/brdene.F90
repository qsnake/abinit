!{\src2tex{textfont=tt}}
!!****f* ABINIT/brdene
!! NAME
!! brdene
!!
!!
!! FUNCTION
!! Update vin according to the Broyden formula, combined
!! with a line minimisation that take into account the total energies.
!! Also transfer vin to vin_prev, vout to vout_prev, and etotal to etotal_prev
!! Could see Numerical Recipes (Fortran), 1986, page 307,
!! as well as Schlegel, J. Comp. Chem. 3, 214 (1982).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  etotal=new total energy (no meaning at output)
!!  hessin(ndim,ndim)=hessian matrix
!!  ndim=size of the hessian and vectors
!!  vout(ndim)=new output vector (no meaning at output)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  etotal_prev=previous total energy; contains input etotal at output
!!  vin(ndim)=new input vector; updated at output
!!  vin_prev(ndim)=previous input vector; contains input vin at output
!!  vout_prev(ndim)=previous output vector; contains input vout at output
!!
!! NOTES
!!
!!
!! PARENTS
!!      pred_bfgs,pred_delocint
!!
!! CHILDREN
!!      findmin
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine brdene(etotal,etotal_prev,hessin,ndim,vin,vin_prev,vout,vout_prev)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'brdene'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndim
 real(dp),intent(in) :: etotal
 real(dp),intent(inout) :: etotal_prev
!arrays
 real(dp),intent(in) :: hessin(ndim,ndim),vout(ndim)
 real(dp),intent(inout) :: vin(ndim),vin_prev(ndim),vout_prev(ndim)

!Local variables-------------------------------
!scalars
 integer :: idim,brd_status
 real(dp) :: d2edv2_1,d2edv2_2,d2edv2_predict,dedv_1,dedv_2,dedv_min
 real(dp) :: dedv_predict,etotal_1,etotal_2,etotal_predict,lambda_1,lambda_2
 real(dp) :: lambda_predict
!arrays
 real(dp),allocatable :: dvin(:),vin_min(:),vout_min(:)

!***************************************************************************

 ABI_ALLOCATE(dvin,(ndim))
 ABI_ALLOCATE(vin_min,(ndim))
 ABI_ALLOCATE(vout_min,(ndim))

 lambda_1=1.0_dp       ; lambda_2=0.0_dp
 etotal_1=etotal      ; etotal_2=etotal_prev
 dvin(:)=vin(:)-vin_prev(:)
 dedv_1=dot_product(vout,dvin)
 dedv_2=dot_product(vout_prev,dvin)
 call findmin(dedv_1,dedv_2,dedv_predict,&
& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,brd_status)

!DEBUG : comes back to usual BFGS !
!lambda_predict=1.0_dp
!dedv_predict=dedv_1
!ENDDEBUG

!Generates vin at the minimum, and an interpolated vout, modified
!to have the right value of dedv_predict, from findmin.
 vin_min(:)=vin_prev(:)+lambda_predict*dvin(:)
 vout_min(:)=vout_prev(:)+lambda_predict*(vout(:)-vout_prev(:))
 dedv_min=dedv_2+lambda_predict*(dedv_1-dedv_2)
!Modify vout_min in order to impose dedv_predict
 vout_min(:)=vout_min(:)+dvin(:)*(dedv_predict-dedv_min)/dot_product(dvin,dvin)

!Previous cartesian coordinates
 etotal_prev=etotal
 vin_prev(:)=vin(:)

!New atomic cartesian coordinates are obtained from vin, hessin and vout
 vin(:)=vin_min(:)
 do idim=1,ndim
   vin(:)=vin(:)-hessin(:,idim)*vout_min(idim)
 end do

!Previous atomic forces
 vout_prev(:)=vout(:)

 ABI_DEALLOCATE(dvin)
 ABI_DEALLOCATE(vin_min)
 ABI_DEALLOCATE(vout_min)

end subroutine brdene





!!***
