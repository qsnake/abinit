!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_angles
!! NAME
!! make_angles
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

subroutine make_angles(deloc,icenter,natom)

 use m_profiling

 use defs_basis
 use defs_mover
 use m_delocint

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_angles'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icenter,natom
 type(ab_delocint),intent(inout) :: deloc
!arrays

!Local variables-------------------------------
!scalars
 integer :: ia1,ia2,iang,ibond,is1,is2,ishift,ja1,ja2
 integer :: jbond,js1,js2
!arrays
 integer,allocatable :: angs_tmp(:,:,:)

! *************************************************************************

!DEBUG
!write(std_out,*) 'make_angs: enter'
!ENDDEBUG


!tentative first allocation: < 6 angles per bond.
 ABI_ALLOCATE(angs_tmp,(2,3,72*natom))

 deloc%nang = 0

 do ibond=1, deloc%nbond
   ia1 = deloc%bonds(1,1,ibond)
   is1 = deloc%bonds(2,1,ibond)
   ia2 = deloc%bonds(1,2,ibond)
   is2 = deloc%bonds(2,2,ibond)
   do jbond=ibond+1,deloc%nbond
     ja1 = deloc%bonds(1,1,jbond)
     ja2 = deloc%bonds(1,2,jbond)
     do ishift=-(icenter-1),+(icenter-1)
       js1 = deloc%bonds(2,1,jbond)+ishift
       js2 = deloc%bonds(2,2,jbond)+ishift

       if      (ia1==ja1 .and. is1==js1) then
         deloc%nang = deloc%nang+1
         angs_tmp(:,1,deloc%nang) = (/ia2,is2/)
         angs_tmp(:,2,deloc%nang) = (/ia1,is1/)
         angs_tmp(:,3,deloc%nang) = (/ja2,js2/)

       else if (ia1==ja2 .and. is1==js2) then
         deloc%nang = deloc%nang+1
         angs_tmp(:,1,deloc%nang) = (/ia2,is2/)
         angs_tmp(:,2,deloc%nang) = (/ia1,is1/)
         angs_tmp(:,3,deloc%nang) = (/ja1,js1/)

       else if (ia2==ja2 .and. is2==js2) then
         deloc%nang = deloc%nang+1
         angs_tmp(:,1,deloc%nang) = (/ia1,is1/)
         angs_tmp(:,2,deloc%nang) = (/ia2,is2/)
         angs_tmp(:,3,deloc%nang) = (/ja1,js1/)

       else if (ia2==ja1 .and. is2==js1) then
         deloc%nang = deloc%nang+1
         angs_tmp(:,1,deloc%nang) = (/ia1,is1/)
         angs_tmp(:,2,deloc%nang) = (/ia2,is2/)
         angs_tmp(:,3,deloc%nang) = (/ja2,js2/)

       end if
       if (deloc%nang > 72*natom) then
         write(std_out,*) 'make_angles : too many angles found > 72*natom'
         stop
       end if
     end do
   end do
!  end jbond do
 end do
!end ibond do

 if (associated(deloc%angs)) nullify(deloc%angs)
 ABI_ALLOCATE(deloc%angs,(2,3,deloc%nang))
 do iang=1,deloc%nang
   deloc%angs(:,:,iang) = angs_tmp(:,:,iang)
 end do
 ABI_DEALLOCATE(angs_tmp)

end subroutine make_angles
!!***
