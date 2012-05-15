!{\src2tex{textfont=tt}}
!!****f* ABINIT/sym_map_atoms
!! NAME
!! sym_map_atoms
!!
!! FUNCTION
!! Compute data on how atoms are mapped by the symmetry operations
!!
!! COPYRIGHT
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine sym_map_atoms(atom_map,natom,nsym,symrec,tnons,xred)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sym_map_atoms'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym
!arrays
 integer,intent(out) :: atom_map(4,natom,nsym)
 integer,intent(in) :: symrec(3,3,nsym)
 real(dp),intent(in) :: tnons(3,nsym),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,isym,jatom
!arrays
 real(dp) :: at(3),dr(3),rat(3),rmx(3)

! *************************************************************************
 
 atom_map(:,:,:) = 0
 do isym = 1, nsym
   do iatom = 1, natom
     at(:) = xred(:,iatom) - tnons(:,isym)
     rat(:) = symrec(1,:,isym)*at(1)+&
&     symrec(2,:,isym)*at(2)+&
&     symrec(3,:,isym)*at(3)
     do jatom = 1, natom
       rmx(:) = rat(:) - xred(:,jatom)
       dr(:) = real(nint(rmx(:)))-rmx(:)
       if (maxval(abs(dr(:)))<tol8) then
         atom_map(1,iatom,isym) = jatom
         atom_map(2:4,iatom,isym)=nint(rmx(:))
         exit
       end if
     end do
   end do
 end do

end subroutine sym_map_atoms
!!***
