!{\src2tex{textfont=tt}}
!!****f* ABINIT/asria_corr
!! NAME
!! asria_corr
!!
!! FUNCTION
!! Imposition of the Acoustic sum rule on the InterAtomic Forces
!!  or on the dynamical matrix directly from the previously calculated d2asr
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! asr=(0 => no ASR, 1 or 2=> the diagonal element is modified to give the ASR,
!!      5 => impose hermitian solution using lapack call)
!! d2asr=matrix used to store the correction needed to fulfill
!! the acoustic sum rule.
!! mpert =maximum number of ipert
!! natom=number of atom
!!
!! OUTPUT
!! Input/Output:
!! d2cart=matrix of second derivatives of total energy, in cartesian coordinates
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      anaddb,elast9,gath3,instr9,mkphbs,thmeig
!!
!! CHILDREN
!!      wrtout
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine asria_corr(asr,d2asr,d2cart,mpert,natom)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'asria_corr'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: asr,mpert,natom
!arrays
 real(dp),intent(in) :: d2asr(2,3,natom,3,natom)
 real(dp),intent(inout) :: d2cart(2,3,mpert,3,mpert)

!Local variables-------------------------------
!scalars
 integer :: idir1,idir2,ipert1,ipert2
 character(len=500) :: message

! *********************************************************************

 if (asr==0) return

 write(message, '(a)' )&
& ' asria_corr : imposition of the ASR for the interatomic forces.'
 call wrtout(std_out,message,'COLL')
!Apply d2asr
 do ipert2=1,natom
   do idir2=1,3
     do ipert1=1,natom
       do idir1=1,3
         d2cart(:,idir1,ipert1,idir2,ipert2)=&
&         d2cart(:,idir1,ipert1,idir2,ipert2)-&
&         d2asr(:,idir1,ipert1,idir2,ipert2)
       end do
     end do
   end do
 end do

end subroutine asria_corr
!!***
