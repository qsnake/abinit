!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_G_rotation_old
!! NAME
!! setup_G_rotation_old
!!
!! FUNCTION
!! Set up tables indicating rotation of G-vectors
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (GMR, VO, LR, RWG, MT, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gvec(3,npw)=coordinates of plane waves, supposed to be ordered in increasing modulus
!! timrev=if 2, take into account time-reversal, 1 otherwise
!! nsym=number of symmetry operations
!! npw=number of planewaves used
!! symrec(3,3,nsym)=symmetry operations in reciprocal space
!!
!! OUTPUT
!!  grottb(npw,2,nsym)=grottb(G,I,S) contains the index no. of (SI) G in the array gvec where I is the identity or the inversion 
!!  grottbm1(npw,2,nsym)=contains the index no of IS^{-1} G 
!!
!! NOTES: 
!!
!! PARENTS
!!      rdm
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine setup_G_rotation_old(only_one_kpt,nsym,symrec,timrev,npw,gvec,grottb,grottbm1)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_G_rotation_old'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nsym,timrev
 logical,intent(in) :: only_one_kpt
!arrays
 integer,intent(in) :: gvec(3,npw),symrec(3,3,nsym)
 integer,intent(inout) :: grottb(npw,timrev,nsym),grottbm1(npw,timrev,nsym)

!Local variables ------------------------------
!scalars
 integer :: ig1,ig2,isym,itim
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: gbase(3),grot(3)

!************************************************************************

!
!=== Set up G-rotation table ===
!* This loop might be CPU consuming in isolated systems
!Therefore we skip it in case of one single k-point
 if (only_one_kpt) then
!  As only_one_kpt is true, the only symmetry needed is identity
   do ig1=1,npw
     grottb  (ig1,1,1)=ig1
     grottbm1(ig1,1,1)=ig1
!    TODO check if also inversion might enter somewhere!!!
   end do
 else 
!  === Several k-points ===
   do ig1=1,npw
!    ish1=g2sh(ig1) ; ss=shlim(ish1) ; ee=shlim(ish1+1)-1
     gbase(:)=gvec(:,ig1)
     do itim=1,timrev
       do isym=1,nsym
         grot=(3-2*itim)*MATMUL(symrec(:,:,isym),gbase)
         found=.FALSE.
         do ig2=1,npw 
!          do ig2=ss,ee ! Looping on shell of ig1 to have better scalling
           if (ALL(ABS(grot(:)-gvec(:,ig2))==0)) then
             found=.TRUE.
             grottb  (ig1,itim,isym)=ig2
             grottbm1(ig2,itim,isym)=ig1
           end if
         end do 
         if (.not.found) then
           write(msg,'(6a,i5,a,i5,1x,2(3i5,a),a,i3,a,i3)')ch10,&
&           ' setup_G_rotation_old : ERROR-',ch10,&
&           ' G-shell not closed',ch10,&
&           ' Initial G vector ',ig1,'/',npw,gbase(:),' Rotated G vector ',grot(:),ch10,&
           ' Through sym ',isym,' and itim ',itim
           call wrtout(std_out,msg,'COLL') 
           call leave_new('COLL')
         end if
       end do 
     end do 
   end do 
 end if !only_one_kpt

end subroutine setup_G_rotation_old
!!***
