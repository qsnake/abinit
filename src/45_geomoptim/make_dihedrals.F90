!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_dihedrals
!! NAME
!! make_dihedrals
!!
!! FUNCTION
!!  (to be completed)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JJ)
!!
!! INPUTS
!!  (to be completed)
!!
!! OUTPUT
!!  (to be completed)
!!
!! PARENTS
!!      make_prim_internals
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine make_dihedrals(badangles,deloc,icenter)

 use m_profiling

 use defs_basis
 use defs_mover
 use m_delocint

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_dihedrals'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icenter
 type(ab_delocint),intent(inout) :: deloc
!arrays
 integer,intent(in) :: badangles(deloc%nang)

!Local variables-------------------------------
!scalars
 integer :: chkdihed,ia1,ia2,ia3,iang,idihed,is1,is2
 integer :: is3,ishift,ja1,ja2,ja3,jang,js1,js2,js3,maxshift
 integer :: minshift
!arrays
 integer,allocatable :: diheds_tmp(:,:,:)

! *************************************************************************

!DEBUG
!write(std_out,*) 'make_dihedrals: enter'
!ENDDEBUG


!tentative first allocation: < 6 dihedrals per angle.
 ABI_ALLOCATE(diheds_tmp,(2,4,6*deloc%nang))

 deloc%ndihed = 0
 diheds_tmp(:,:,:) = 0

 do iang=1,deloc%nang
   if (badangles(iang) == 1) cycle
   ia1 = deloc%angs(1,1,iang)
   is1 = deloc%angs(2,1,iang)
   ia2 = deloc%angs(1,2,iang)
   is2 = deloc%angs(2,2,iang)
   ia3 = deloc%angs(1,3,iang)
   is3 = deloc%angs(2,3,iang)

   do jang=iang+1,deloc%nang
     if (badangles(jang) == 1) cycle
     ja1 = deloc%angs(1,1,jang)
     ja2 = deloc%angs(1,2,jang)
     ja3 = deloc%angs(1,3,jang)
     do ishift=-(icenter-1),(icenter-1)
       js1 = deloc%angs(2,1,jang)+ishift
       js2 = deloc%angs(2,2,jang)+ishift
       js3 = deloc%angs(2,3,jang)+ishift

       chkdihed=0
       if (ia2==ja1 .and. is2==js1) then
         if (ia1==ja2 .and. is1==js2) then
           deloc%ndihed = deloc%ndihed+1
           diheds_tmp(:,1,deloc%ndihed) = (/ia3,is3/)
           diheds_tmp(:,2,deloc%ndihed) = (/ia2,is2/)
           diheds_tmp(:,3,deloc%ndihed) = (/ja2,js2/)
           diheds_tmp(:,4,deloc%ndihed) = (/ja3,js3/)
           chkdihed=1
         else if (ia3==ja2 .and. is3==js2) then
           deloc%ndihed = deloc%ndihed+1
           diheds_tmp(:,1,deloc%ndihed) = (/ia1,is1/)
           diheds_tmp(:,2,deloc%ndihed) = (/ia2,is2/)
           diheds_tmp(:,3,deloc%ndihed) = (/ja2,js2/)
           diheds_tmp(:,4,deloc%ndihed) = (/ja3,js3/)
           chkdihed=1
         end if
       else if (ia2==ja3 .and. is2==js3) then
         if (ia1==ja2 .and. is1==js2) then
           deloc%ndihed = deloc%ndihed+1
           diheds_tmp(:,1,deloc%ndihed) = (/ia3,is3/)
           diheds_tmp(:,2,deloc%ndihed) = (/ia2,is2/)
           diheds_tmp(:,3,deloc%ndihed) = (/ja2,js2/)
           diheds_tmp(:,4,deloc%ndihed) = (/ja1,js1/)
           chkdihed=1
         else if (ia3==ja2 .and. is3==js2) then
           deloc%ndihed = deloc%ndihed+1
           diheds_tmp(:,1,deloc%ndihed) = (/ia1,is1/)
           diheds_tmp(:,2,deloc%ndihed) = (/ia2,is2/)
           diheds_tmp(:,3,deloc%ndihed) = (/ja2,js2/)
           diheds_tmp(:,4,deloc%ndihed) = (/ja1,js1/)
           chkdihed=1
         end if
       end if
       if (deloc%ndihed > 6*deloc%nang) then
         write(std_out,*) 'make_dihedrals : too many dihedrals found > 6*nang'
         stop
       end if
       if (chkdihed == 1) then
         if (   diheds_tmp(1,4,deloc%ndihed) == diheds_tmp(1,1,deloc%ndihed) .and.&
&         diheds_tmp(2,4,deloc%ndihed) == diheds_tmp(2,1,deloc%ndihed) ) then
           write(std_out,*) 'make_dihedrals : Bad dihedral was found: atom1 == atom4. Discarding.'
           diheds_tmp(:,:,deloc%ndihed) = 0
           deloc%ndihed = deloc%ndihed-1
         end if
       end if
     end do
   end do
!  end jang do
 end do
!end iang do

 if (associated(deloc%dihedrals)) nullify(deloc%dihedrals)
 ABI_ALLOCATE(deloc%dihedrals,(2,4,deloc%ndihed))
 do idihed=1,deloc%ndihed
   deloc%dihedrals(:,:,idihed) = diheds_tmp(:,:,idihed)

!  minshift = minval(diheds_tmp(2,:,idihed))
!  if (minshift <= 0) then
!  deloc%dihedrals(2,:,idihed) = deloc%dihedrals(2,:,idihed)+minshift+1
!  end if
!  maxshift = maxval(diheds_tmp(2,:,idihed))
!  if (maxshift > deloc%nrshift) then
!  deloc%dihedrals(2,:,idihed) = deloc%dihedrals(2,:,idihed)-maxshift
!  end if
!  
   minshift = minval(diheds_tmp(2,:,idihed))
   maxshift = maxval(diheds_tmp(2,:,idihed))
   if (minshift <= 0 .or. maxshift > deloc%nrshift) then
     write(std_out,*) ' make_dihedrals : Error : dihedral extends beyond '
     write(std_out,*) '  first neighboring unit cells ! '
     stop
   end if
 end do
 ABI_DEALLOCATE(diheds_tmp)

end subroutine make_dihedrals
!!***
