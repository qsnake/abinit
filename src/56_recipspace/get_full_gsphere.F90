!{\src2tex{textfont=tt}}
!!****f* ABINIT/getfullg
!! NAME
!! getfullg
!!
!! FUNCTION
!!  Reconstruct a G-sphere starting from a set of irreducible lattice vectors
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  pinv=-1 if time-reversal can be used, 1 otherwise
!!  nsym=number of symmetry operations
!!  sizepw=Max expected number of G vectors in the shere
!!  symrec(3,3,nsym)=symmetry operation in reciprocal space
!!  nbase=number of irreducible G vectors
!!  gbase(3,nbase)=irreducible G-vectors
!!  cnorm(nbase)=norm of the irreducible G vectors (supposed not yet sorted) 
!!
!! OUTPUT
!!  maxpw=Number of G vectors found
!!  gbig(3,sizepw)=G vectors in the sphere packed in the first maxpw columns
!!  shlim(nbase)=number of G vectors within each shell
!!  ierr= Exit status, if /=0 the number of G vectors found exceeds sizepw
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  cnorm is a bit redundant since it can be calculated from gbase. However this procedure 
!!  is called by outkss in which cnorm is already calculated and we dont want to do it twice
!!
!! PARENTS
!!      m_gsphere
!!
!! CHILDREN
!!      sort_dp,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getfullg(nbase,nsym,pinv,sizepw,gbase,symrec,cnorm,maxpw,gbig,shlim,ierr)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getfullg'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbase,nsym,pinv,sizepw
 integer,intent(out) :: ierr,maxpw
!arrays
 integer,intent(in) :: gbase(3,nbase),symrec(3,3,nsym)
 integer,intent(out) :: gbig(3,sizepw),shlim(nbase)
 real(dp),intent(inout) :: cnorm(nbase) !sort_dp can change cnorm

!Local variables-------------------------------
!scalars
 integer :: ibase,ig,ilim,ish,isym,itim
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: gcur(3),geq(3)
 integer,allocatable :: gshell(:,:),insort(:),nshell(:)

! *************************************************************************
 
 if (pinv/=1.and.pinv/=-1) then
   write(msg,'(a,i6)')&
&   ' The argument pinv should be -1 or 1, however, pinv =',pinv
   MSG_BUG(msg)
 end if
!
!=== Reorder base g-vectors in order of increasing module ===
 ABI_ALLOCATE(insort,(nbase))
 do ibase=1,nbase
   insort(ibase)=ibase
 end do
 call sort_dp(nbase,cnorm,insort,tol14)
!
!=== Generate all stars of G-vectors ===
!Star of G is the set of all symetrical images of the vector
!gshell contains the symmetrical G at fixed gbase. No need to add an additional dimension
!or initialize to zero the array inside the loop over nbase as we loop over (ish<=nshell(ibase))
 ABI_ALLOCATE(nshell,(nbase))
 ABI_ALLOCATE(gshell,(3,2*nsym))
!
!=== Start with zero number of G vectors found ===
 maxpw=0 ; ierr=0
 do ibase=1,nbase
!  
!  === Loop over all different modules of G ===
!  * Start with zero G vectors found in this star
   nshell(ibase)=0
   gcur(:)=gbase(:,insort(ibase))
!  
!  === Loop over symmetries ===
   do isym=1,nsym
     do itim=pinv,1,2
       geq(:)=itim*MATMUL(symrec(:,:,isym),gcur)
!      
!      * Search for symetric of g and eventually add it:
       found=.FALSE. ; ish=1
       do while ((.not.found).and. (ish<=nshell(ibase)))
         found=ALL(geq(:)==gshell(:,ish))
         ish=ish+1
       end do
       if (.not.found) then
         nshell(ibase)=nshell(ibase)+1
         gshell(:,nshell(ibase))=geq(:)
       end if
     end do
   end do
!  
!  * Was sizepw large enough?
   if ((maxpw+nshell(ibase))>sizepw) then
     write(msg,'(a,i6,2a)')&
&     ' Number of G in sphere exceeds maximum allowed value =',sizepw,ch10,&
&     ' check the value of sizepw in calling routine ' 
     MSG_WARNING(msg)
     ierr=1; RETURN
   end if
!  
!  === Store this shell of Gs in a big array (gbig) ===
   do ig=1,nshell(ibase)
     gbig(:,ig+maxpw)=gshell(:,ig)
   end do
   maxpw=maxpw+nshell(ibase)
 end do ! ibase
!
!=== Compute number of G"s within each shell ===
 ilim=0
 do ibase=1,nbase
   ilim=ilim+nshell(ibase)
   shlim(ibase)=ilim
 end do
!
!=== Print out shell limits ===
 write(msg,'(3a)')&
& ' Shells found:',ch10,&
& ' number of shell    number of G vectors      cut-off energy [Ha] '
 call wrtout(std_out,msg,'COLL')

 do ibase=1,nbase
   write(msg,'(12x,i4,17x,i6,12x,f8.3)')ibase,shlim(ibase),two*pi**2*cnorm(ibase)
   call wrtout(std_out,msg,'COLL')
 end do
 write(msg,'(a)')ch10 
 call wrtout(std_out,msg,'COLL')
 ABI_DEALLOCATE(gshell)
 ABI_DEALLOCATE(insort)
 ABI_DEALLOCATE(nshell)

end subroutine getfullg
!!***
