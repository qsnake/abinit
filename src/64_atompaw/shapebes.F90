!{\src2tex{textfont=tt}}
!!****f* ABINIT/shapebes
!! NAME
!! shapebes
!!
!! FUNCTION
!!    Find al and ql parameters for a "Bessel" shape function:
!!    Shape(r)=al1.jl(ql1.r)+al2.jl(ql2.r)
!!      such as Shape(r) and 2 derivatives are zero at r=rc
!!              Intg_0_rc[Shape(r).r^(l+2).dr]=1
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  ll= l quantum number
!!  rc= cut-off radius
!!
!! OUTPUT
!!  al(2)= al coefficients
!!  ql(2)= ql factors
!!
!! PARENTS
!!      psp7in
!!
!! CHILDREN
!!      jbessel,solvbes
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


 subroutine shapebes(al,ql,ll,rc)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shapebes'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: ll
 real(dp) :: rc
!arrays
 real(dp) :: al(2),ql(2)

!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: alpha,beta,det,jbes,jbesp,jbespp,qr
!arrays
 real(dp) :: amat(2,2),bb(2)

! *************************************************************************

 alpha=1._dp;beta=0._dp
 call solvbes(ql,alpha,beta,ll,2)
 ql(1:2)=ql(1:2)/rc

 do ii=1,2
   qr=ql(ii)*rc
   call jbessel(jbes,jbesp,jbespp,ll,1,qr)
   amat(1,ii)=jbesp*ql(ii)
   call jbessel(jbes,jbesp,jbespp,ll+1,0,qr)
   amat(2,ii)=jbes*rc**(ll+2)/ql(ii)  !  Intg_0_rc[jl(qr).r^(l+2).dr]
 end do

 bb(1)=zero;bb(2)=one

 det=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1)
 al(1)=(amat(2,2)*bb(1)-amat(1,2)*bb(2))/det
 al(2)=(amat(1,1)*bb(2)-amat(2,1)*bb(1))/det

end subroutine shapebes
!!***
