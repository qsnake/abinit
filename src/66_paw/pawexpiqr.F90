!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawexpiqr
!! NAME
!! pawexpiqr
!!
!! FUNCTION
!! Compute exp(-i.q.r) for each point of the (fine) rectangular grid
!! around a given atomic site.
!! Used for the determination of phonons at non-zero q wavevector.
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gprimd(3,3)= dimensional primitive translations for reciprocal space
!!  pawfgrtab <type(pawfgrtab_type)>= atomic data given on fine rectangular grid
!!  qphon(3)= wavevector of the phonon
!!  xred(3)= reduced atomic coordinates
!!
!! OUTPUT
!!  pawfgrtab%expiqr(2,nfgd)= exp(i.q.r) around the current atom
!!                            Not allocated if q=0 !
!!
!! PARENTS
!!      pawdij,pawfrnhat,pawgrnl,pawmknhat,respfn
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawexpiqr(gprimd,pawfgrtab,qphon,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawexpiqr'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
!arrays
 real(dp),intent(in) :: gprimd(3,3),qphon(3),xred(3)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab

!Local variables ------------------------------
!scalars
 integer :: ic
 logical :: qne0
 real(dp) :: phase,phase_xred,qx,qy,qz
 character(len=500) :: msg
!arrays

! *************************************************************************

 DBG_ENTER("COLL")

 qne0=(qphon(1)**2+qphon(2)**2+qphon(3)**2>=1.d-15)

!Compatibility tests
 if (pawfgrtab%rfgd_allocated==0.and.qne0) then
   msg='  pawfgrtab()%rfgd array must be allocated  !'
   MSG_BUG(msg)
 end if

!Compute q in cartesian coordinates
 if (qne0) then
   qx=gprimd(1,1)*qphon(1)+gprimd(1,2)*qphon(2)+gprimd(1,3)*qphon(3)
   qy=gprimd(2,1)*qphon(1)+gprimd(2,2)*qphon(2)+gprimd(2,3)*qphon(3)
   qz=gprimd(3,1)*qphon(1)+gprimd(3,2)*qphon(2)+gprimd(3,3)*qphon(3)
   phase_xred=two_pi*(qphon(1)*xred(1)+qphon(2)*xred(2)+qphon(3)*xred(3))
 end if

!Allocate array for exp(i.q.r)
 if (associated(pawfgrtab%expiqr))  then
   ABI_DEALLOCATE(pawfgrtab%expiqr)
 end if
 pawfgrtab%expiqr_allocated=0
 if (qne0) then
   ABI_ALLOCATE(pawfgrtab%expiqr,(2,pawfgrtab%nfgd))
   pawfgrtab%expiqr_allocated=1
 end if

!Compute exp(i.q.r)
 if (qne0) then
   do ic=1,pawfgrtab%nfgd
     phase=two_pi*(qx*pawfgrtab%rfgd(1,ic)+qy*pawfgrtab%rfgd(2,ic) &
&     +qz*pawfgrtab%rfgd(3,ic)) + phase_xred
     pawfgrtab%expiqr(1,ic)=cos(phase)
     pawfgrtab%expiqr(2,ic)=sin(phase)
   end do
 end if

 DBG_EXIT("COLL")

end subroutine pawexpiqr
!!***
