!{\src2tex{textfont=tt}}
!!****f* ABINIT/merge_kgirr
!! NAME
!! merge_kgirr
!!
!! FUNCTION
!!  Given a list of irreducible reciprocal vectors associated to different k-centered spheres, 
!!  this subroutine finds the minimal set of G vectors needed to reconstruct the union of the spheres 
!!  through symmetry operations.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! nsym=number of symmetry operations
!! pinv=-1 if time-reversal can be used, 0 otherwise
!! nkpt=number of k-points for k-centered spheres
!! mpw=Max number of G vectors for each k-point
!! sizepw=Max expected number of G vectors found
!! symrec(3,3,nsym)=symmetries in reciprocal space given in reduced coordinates
!! nbasek(nkpt)=number of irred G for each k-point
!! cnormk(mpw,nkpt)=the norm of each k-centered G (only 1:nbase(ik)) is used
!! gbasek(3,mpw,nkpt)
!!
!! OUTPUT
!! nbase=number of irreducible G needed to reconstruct the initial set of spheres
!! gbase(3,sizepw)=irreducible G found in reciprocal coordinates
!! cnorm(sizepw)=Norm of each irred G vector
!! ierr= Exit status, if /=0 the number of G vectors found exceeds sizepw
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

subroutine merge_kgirr(nsym,pinv,nkpt,mpw,sizepw,symrec,nbasek,cnormk,gbasek,nbase,gbase,cnorm,ierr)

 use m_profiling

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'merge_kgirr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpw,nkpt,nsym,pinv,sizepw
 integer,intent(out) :: ierr,nbase
!arrays
 integer,intent(in) :: gbasek(3,mpw,nkpt),nbasek(nkpt),symrec(3,3,nsym)
 integer,intent(out) :: gbase(3,sizepw)
 real(dp),intent(in) :: cnormk(mpw,nkpt)
 real(dp),intent(out) :: cnorm(sizepw)

!Local variables-------------------------------
!scalars
 integer :: ikpt,inb,irgk,isym
 real(dp) :: eps,norm
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: gbas(3),gcur(3),geq(3)

! *************************************************************************

 DBG_ENTER("COLL")

 if (pinv/=1.and.pinv/=-1) then
   write(msg,'(a,i6)')' The argument pinv should be -1 or 1, however, pinv =',pinv
   MSG_BUG(msg)
 end if
!
!=== Start with zero number of G found ===
 nbase=0 ; ierr=0
 do ikpt=1,nkpt
   do irgk=1,nbasek(ikpt)
     gcur(:)=gbasek(:,irgk,ikpt)
     norm=cnormk(irgk,ikpt) ; eps=tol8*norm
     found=.FALSE. ; inb=1
     do while ((.not.found).and.(inb<=nbase))
       if (ABS(norm-cnorm(inb))<=eps) then
         gbas(:)=gbase(:,inb)
         isym=1
         do while ((.not.found).and.(isym<=nsym))
           geq(:)=MATMUL(symrec(:,:,isym),gcur)
           found=ALL(geq(:)==gbas(:))
           if (pinv==-1) found= (found.or.ALL(geq(:)==-gbas(:)) ) ! For time-reversal
           isym=isym+1
         end do
       end if
       inb=inb+1
     end do
     if (.not.found) then
!      === Add to the list ===
       nbase=nbase+1
       if (nbase>sizepw) then
         write(msg,'(2(a,i5),a)')&
&         ' nbase (',nbase,') became greater than sizepw = ',sizepw,' returning ierr=1 '
         MSG_WARNING(msg)
         ierr=1; RETURN
       end if
       cnorm(nbase)=cnormk(irgk,ikpt)
       gbase(:,nbase)=gcur(:)
     end if
   end do
 end do

 DBG_EXIT("COLL")

end subroutine merge_kgirr
!!***
