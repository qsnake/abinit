!{\src2tex{textfont=tt}}
!!****f* ABINIT/print_bonds
!! NAME
!! print_bonds
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JJ)
!!
!! INPUTS
!!  natom=  Number of atoms
!!  ntypat= Number of type of atoms
!!  znucl=  Z number of the atom
!!
!! OUTPUT
!!  bonds= Structure that store all the information about
!!         bonds created by this routine:
!!         nbonds=  Total number of bonds
!!         bondij=  Unitary vector along the bond direction
!!         distij=  Distances between atoms i and j (including shift)
!!         listij= Indices of bonds going from i to j
!!         listji= Indices of bonds going from j to i
!!         indexij= Number of bonds between i and j
!!         indexji= Number of bonds between j and i
!!         tolerance
!!
!! PARENTS
!!      make_angles_new,make_bonds_new
!!
!! CHILDREN
!!
!! SOURCE
  
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine print_bonds(amu,bonds,natom,ntypat,symbol,typat,znucl)

 use m_profiling

use defs_basis
use defs_mover

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_bonds'
 use interfaces_45_geomoptim, except_this_one => print_bonds
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom,ntypat
integer,intent(in) :: typat(natom)
real(dp),intent(in) :: znucl(ntypat)
real(dp) :: amu(ntypat)
character(len=2) :: symbol(ntypat)
type(go_bonds),intent(in) :: bonds

!Local variables ------------------------------
!scalars
integer :: ii,jj,kk

! *********************************************************************

 write(std_out,'(a)') ch10
 write(std_out,'(a,72a,a)') '---BONDS',('-',kk=1,72),ch10
 write(std_out,'(a,i3)') 'Number of atoms:   ',natom
 write(std_out,'(a,i3)') 'Number of bonds:   ',bonds%nbonds
 write(std_out,'(a,f6.3,a,a)') 'Tolerance of bonds: ',bonds%tolerance,' times the sum of covalent radius',ch10

 do ii=1,natom
   write(std_out,'(a,i3)') 'ATOM number:       ',ii
   write(std_out,'(a,f8.3)') '  Z:              ',znucl(typat(ii))
   write(std_out,'(a,f8.3)') '  Weight:         ',amu(typat(ii))
   write(std_out,'(a,a3)') '  Symbol:          ',symbol(typat(ii))
   write(std_out,'(a,i3)') '  Number of bonds: ',bonds%nbondi(ii)

   do jj=1,bonds%nbondi(ii)
     write(std_out,'(a,i3,a,a,i3,a,3f7.3,a,f7.3)') '    [',jj,']',&
&     '    Index of bond: ',bonds%indexi(ii,jj),&
&     '    Unitary vector: ',bonds%bond_vect(:,abs(bonds%indexi(ii,jj))),&
&     '    Bond length: ',bonds%bond_length(abs(bonds%indexi(ii,jj)))
   end do

 end do

 do ii=1,bonds%nbonds

   write(std_out,'(a,i3)') 'BOND Index=',ii
   write(std_out,'(a,3f8.3)') '    Vector',bonds%bond_vect(:,ii)
   write(std_out,'(a,f8.3)')  '    bond Length',bonds%bond_length(ii)

 end do

end subroutine print_bonds
!!***
