!{\src2tex{textfont=tt}}
!!****f* ABINIT/getwtk
!! NAME
!! getwtk
!!
!! FUNCTION
!! Routine called by the program optic
!! Presumes kpts are the irreducible ones of a good uniform grid
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (SSharma,MVer,VRecoules,TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  kpt(3,nkpt)=reduced coordinates of k points.
!!  nkpt = number of k points
!!  nsym=Number of symmetry operations.
!!  symrel(3,3,nsym)=symmetry operations
!!
!! OUTPUT
!!  wtk(nkpt)=weight assigned to each k point.
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!      dgemv,wrap2_pmhalf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine getwtk(kpt,nkpt,nsym,symrel,wtk)

 use m_profiling

 use defs_basis
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getwtk'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments -----------------------------------------------
! in
! out
!scalars
 integer,intent(in) :: nkpt,nsym
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: kpt(3,nkpt)
 real(dp),intent(out) :: wtk(nkpt)

!Local variables -----------------------------------------
!scalars
 integer :: ikpt,istar,isym,itim,new,nkpt_tot
 real(dp) :: shift,timsign,tmp
!arrays
 integer :: nstar(nkpt)
 real(dp) :: dkpt(3),kptstar(3,2*nkpt*nsym),rsymrel(3,3,nsym),symkpt(3)
 real(dp) :: tsymkpt(3)

! *************************************************************************

 do isym=1,nsym
   rsymrel(:,:,isym) = dble(symrel(:,:,isym))
 end do

!for each kpt find star and accumulate nkpts
 do ikpt=1,nkpt
   write(std_out,*) ' getwtk : ikpt = ', ikpt
   nstar(ikpt) = 0
   kptstar(:,:) = zero
   do isym=1,nsym

     call dgemv('N',3,3,one,rsymrel(:,:,isym),3,kpt(:,ikpt),1,zero,symkpt,1)

!    is symkpt already in star?
     do itim=0,1
       timsign=one-itim*two
       tsymkpt(:) = timsign*symkpt(:)
       call wrap2_pmhalf(tsymkpt(1),tmp,shift) ;  tsymkpt(1) = tmp
       call wrap2_pmhalf(tsymkpt(2),tmp,shift) ;  tsymkpt(2) = tmp
       call wrap2_pmhalf(tsymkpt(3),tmp,shift) ;  tsymkpt(3) = tmp
       new=1
       do istar=1,nstar(ikpt)
         dkpt(:) = abs(tsymkpt(:)-kptstar(:,istar))
         if ( sum(dkpt) < 1.0d-6) then
           new=0
           exit
         end if
       end do
       if (new==1) then
         nstar(ikpt) = nstar(ikpt)+1
         kptstar(:,nstar(ikpt)) = tsymkpt(:)
       end if
     end do

   end do
!  end do nsym
!  DEBUG
!  write(std_out,*) ' getwtk : nstar = ', nstar(ikpt)
!  write(std_out,*) ' getwtk : star = '
!  write(std_out,*)  kptstar(:,1:nstar(ikpt))
!  ENDDEBUG
 end do
!end do nkpt

 nkpt_tot = sum(nstar)
!write(std_out,*) ' getwtk : nkpt_tot = ', nkpt_tot
 do ikpt=1,nkpt
   wtk(ikpt) = dble(nstar(ikpt))/dble(nkpt_tot)
 end do

end subroutine getwtk
!!***
