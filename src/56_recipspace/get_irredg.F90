!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_irredg
!! NAME
!! get_irredg
!!
!! FUNCTION
!!  Given a set of reciprocal lattice vectors, find the set of G"s generating the others by symmetry.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MT, VO, AR, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nsym=number of symmetry operations
!!  pinv=-1 if time-reversal can be used, 1 otherwise
!!  npw_k=number of G vectors (for this k-point, as the set of G is k-centered)   
!!  gcurr(3,npw_k)=the list of G vectors
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  symrec(3,3,nsym)=symmetry operations in terms of reciprocal space primitive translations.
!!
!! OUTPUT
!!  nbasek=number of irreducible G vectors found 
!!  cnormk(npw_k)=first nbasek elements are the norm of each irreducible G-vector
!!  gbasek(3,npw_k)=first nbasek elements are the irreducible G vectors
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_gsphere
!!
!! CHILDREN
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine get_irredg(npw_k,nsym,pinv,gprimd,symrec,gcurr,nbasek,gbasek,cnormk)

 use m_profiling

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_irredg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,nsym,pinv
 integer,intent(out) :: nbasek
!arrays
 integer,intent(in) :: gcurr(3,npw_k),symrec(3,3,nsym)
 integer,intent(out) :: gbasek(3,npw_k)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(out) :: cnormk(npw_k)

!Local variables-------------------------------
!scalars
 integer :: ig,irr,isym,jj
 real(dp) :: eps,norm
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: gbas(3),gcur(3),geq(3)
 real(dp) :: gcar(3)

! *************************************************************************

 DBG_ENTER("COLL")

 if (pinv/=1.and.pinv/=-1) then
   write(msg,'(a,i6)')&
&   ' The argument pinv should be -1 or 1, however, pinv =',pinv
   MSG_BUG(msg)
 end if
!
!=== zero irred G vectors found, zeroing output arrays ===
 nbasek=0 ; cnormk(:)=zero ; gbasek(:,:)=0

 do ig=1,npw_k
   gcur(:)=gcurr(:,ig) ; norm=zero
   do jj=1,3
     gcar(jj)=gcur(1)*gprimd(jj,1)+gcur(2)*gprimd(jj,2)+gcur(3)*gprimd(jj,3)
     norm=norm+gcar(jj)**2
   end do
   eps=tol8*norm ; found=.FALSE. ; irr=1
   do while ((.not.found).and.(irr<=nbasek))
     if (ABS(norm-cnormk(irr))<=eps) then
       gbas(:)=gbasek(:,irr)
       isym=1
       do while ((.not.found).and.(isym<=nsym))
         geq(:)=MATMUL(symrec(:,:,isym),gcur)
         found=ALL(geq(:)==gbas(:))
         if (pinv==-1) found=(found.or.ALL(geq==-gbas)) ! For time-reversal
         isym=isym+1
       end do
     end if
     irr=irr+1
   end do
   if (.not.found) then
     nbasek=nbasek+1
     cnormk(nbasek)=norm
     gbasek(:,nbasek)=gcur(:)
   end if
 end do

 DBG_EXIT("COLL")

end subroutine get_irredg
!!***
