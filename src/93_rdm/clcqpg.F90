!{\src2tex{textfont=tt}}
!!****f* ABINIT/clcqpg
!! NAME
!! clcqpg
!!
!! FUNCTION
!! Calculate |q+G| for each q and G
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gvec(3,npwx)=reduced coordinates of the plane waves 
!!  npwx=number of plane waves used for the exchange part
!!  nq=number of q points
!!  qq(3,nq)=coordinates of q points
!!
!! OUTPUT
!!  qpg(npwx,nq)=norm of q+G vector
!!
!! PARENTS
!!      rdm
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine clcqpg(npwx,gvec,gprimd,qq,nq,qpg)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clcqpg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwx,nq
!arrays
 integer,intent(in) :: gvec(3,npwx)
 real(dp),intent(in) :: gprimd(3,3),qq(3,nq)
 real(dp),intent(out) :: qpg(npwx,nq)

!Local variables ------------------------------
!scalars
 integer :: ig,ii,iq
!arrays
 real(dp) :: gmet(3,3),gpq(3)

!************************************************************************

!Compute reciprocal space metrics
 do ii=1,3
   gmet(ii,:)=gprimd(1,ii)*gprimd(1,:)+&
&   gprimd(2,ii)*gprimd(2,:)+&
&   gprimd(3,ii)*gprimd(3,:)
 end do

 do iq=1,nq
   do ig=1,npwx
     gpq(:)=gvec(:,ig)+qq(:,iq)
     qpg(ig,iq)=two_pi*SQRT(dot_product(gpq,MATMUL(gmet,gpq)))
   end do
 end do

end subroutine clcqpg
!!***
