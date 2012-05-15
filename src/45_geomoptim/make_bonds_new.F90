!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_bonds_new
!! NAME
!! make_bonds_new
!!
!! FUNCTION
!!  Fill the contents of the bonds structure, that contains
!!  all non redundant bonds that could be generated between
!!  all the atoms in the unitary cell and their adjacent cells
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JJ)
!!
!! INPUTS
!!  natom=  Number of atoms
!!  ntypat= Number of type of atoms
!!  rprimd= Dimensional primitive vectors of the cell
!!  xcart=  Cartesian coordinates of the atoms
!!  znucl=  Z number of the atom
!!
!! OUTPUT
!!  bonds= Structure that store all the information about
!!         bonds created by this routine:
!!         nbonds=  Total number of bonds
!!         nbondi=  Number of bonds for atom i
!!         indexi=  Indeces of bonds for atom i
!!         bond_lenght=  Distances between atoms i and j (including shift)
!!         bond_vect=    Unitary vector for the bond from i to j
!!         tolerance=    The tolerance is multiplied to the
!!                       adition of covalent radius to decide if a bond is created
!!
!! PARENTS
!!      prec_simple
!!
!! CHILDREN
!!      atmdata,print_bonds
!!
!! SOURCE
  
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine make_bonds_new(bonds,natom,ntypat,rprimd,typat,xcart,znucl)

 use m_profiling

use defs_basis
use defs_mover

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_bonds_new'
 use interfaces_32_util
 use interfaces_45_geomoptim, except_this_one => make_bonds_new
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom,ntypat
!arrays
integer,intent(in) :: typat(natom)
real(dp),intent(in) :: znucl(ntypat)
real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
type(go_bonds),intent(inout) :: bonds

!Local variables ------------------------------
!scalars
integer :: ii,jj,kk,ibond,irshift
real(dp) :: rcov1,rcov2
real(dp) :: bl
type(go_bonds) :: bonds_tmp

!arrays
character(len=2) :: symbol(ntypat)
real(dp) :: amu(ntypat)
integer :: shift(3,13) ! Represent all shift vectors that are not equivalent by central symmetry
! For example (1,1,-1) is equivalent to (-1,-1,1)
! It means that bond between atom i in the original cell and atom j in the
! cell with cordinates (1,1,-1) is equivalent to the bond between atom j in 
! the orignal cell and atom i in the cell with coordinates (-1,-1,1)
! The trivial shift (0,0,0) is excluded here
real(dp) :: rcov(ntypat) ! Covalent radius
real(dp) :: rpt(3)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

!write(std_out,*) 'make_bonds 01'
!##########################################################
!### 01. Compute covalent radius

 do ii=1,ntypat
   call atmdata(amu(ii),rcov(ii),symbol(ii),znucl(ii))
 end do

!write(std_out,*) 'make_bonds 02'
!##########################################################
!### 02. Fill the 13 posible shift conecting adjacent cells

 shift(:,:)=reshape( (/ 1,0,0,&
& 0, 1, 0,&
& 0, 0, 1,&
& 1, 1, 0,&
& 1,-1, 0,&
& 0, 1, 1,&
& 0, 1,-1,&
& 1, 0, 1,&
& 1, 0,-1,&
& 1, 1, 1,&
& 1,-1, 1,&
& 1, 1,-1,&
& 1,-1,-1 /), (/ 3, 13 /))

!write(std_out,*) 'make_bonds 03'
!##########################################################
!### 03. Initialize the values of bonds

!The total number of bonds could not be predicted without
!compute all the distances, but the extreme case is linking
!all the atoms within all adjacent cells (natom*natom*13)
!plus the all the bonds inside the original cell (natom*(natom-1))

 bonds_tmp%nbonds=0
 bonds_tmp%tolerance=bonds%tolerance
 ibond=0

 ABI_ALLOCATE(bonds_tmp%bond_vect,(3,natom*natom*14-natom))
 ABI_ALLOCATE(bonds_tmp%bond_length,(natom*natom*14-natom))

!indexi contains the indeces to the bonds
 ABI_ALLOCATE(bonds_tmp%indexi,(natom,natom*natom*14-natom))

 ABI_ALLOCATE(bonds_tmp%nbondi,(natom))

 bonds_tmp%indexi(:,:)=0
 bonds_tmp%nbondi(:)=0
 bonds_tmp%bond_vect(:,:)=0.0
 bonds_tmp%bond_length(:)=0.0

!write(std_out,*) 'make_bonds 04'
!##########################################################
!### 04. Compute the bonds inside the original cell
!### shift=(0,0,0)

 do ii=1,natom
   rcov1 = rcov(typat(ii))

   do jj=ii+1,natom
     rcov2 = rcov(typat(jj))

     bl=bond_length(xcart(:,ii),xcart(:,jj))

     if (bonds_tmp%tolerance*(rcov1+rcov2) > bl) then
!      We have a new bond, nbonds starts from
!      0, so it could be used to index the
!      locations of bondij and distij

!      Increase the number of bonds
       bonds_tmp%nbonds= bonds_tmp%nbonds+1

!      The number of bonds for atoms ii and jj
!      needs to raise by one
       bonds_tmp%nbondi(ii)= bonds_tmp%nbondi(ii)+1
       bonds_tmp%nbondi(jj)= bonds_tmp%nbondi(jj)+1

       bonds_tmp%indexi(ii,bonds_tmp%nbondi(ii))=bonds_tmp%nbonds
!      The value for jj is negative to indicate that
!      the vector is from ii to jj
       bonds_tmp%indexi(jj,bonds_tmp%nbondi(jj))=-bonds_tmp%nbonds

!      The unitary vector is always from ii to jj
       bonds_tmp%bond_vect(:,bonds_tmp%nbonds)=(xcart(:,jj)-xcart(:,ii))/bl
       bonds_tmp%bond_length(bonds_tmp%nbonds)=bl

     end if

   end do !! jj
 end do !! ii

!write(std_out,*) 'make_bonds 05'
!##########################################################
!### 05. Compute the bonds outside the original cell
!###     13 shifts considered

!Bonds between identical atoms but in diferent cells are
!allowed 

 do ii=1,natom
   rcov1 = rcov(typat(ii))
   do jj=1,natom
     rcov2 = rcov(typat(jj))

     do irshift=1,13

       do kk=1,3
         rpt(kk) = xcart(kk,jj)+&
&         shift(1,irshift)*rprimd(kk,1)+ &
&         shift(2,irshift)*rprimd(kk,2)+ &
&         shift(3,irshift)*rprimd(kk,3)
       end do


       bl =bond_length(xcart(:,ii),rpt)

       if (bonds_tmp%tolerance*(rcov1+rcov2) > bl) then

!        We have a new bond, nbonds starts from
!        0, so it could be used to index the
!        locations of bondij and distij

!        Increase the number of bonds
         bonds_tmp%nbonds= bonds_tmp%nbonds+1
         
!        The number of bonds for atoms ii and jj
!        needs to raise by one
         bonds_tmp%nbondi(ii)= bonds_tmp%nbondi(ii)+1
         bonds_tmp%indexi(ii,bonds_tmp%nbondi(ii))=bonds_tmp%nbonds

!        The value for jj is negative to indicate that
!        the vector is from ii to jj
         bonds_tmp%nbondi(jj)= bonds_tmp%nbondi(jj)+1
         bonds_tmp%indexi(jj,bonds_tmp%nbondi(jj))=-bonds_tmp%nbonds

!        The unitary vector is always from ii to jj
         bonds_tmp%bond_vect(:,bonds_tmp%nbonds)=(rpt(:)-xcart(:,ii))/bl
         bonds_tmp%bond_length(bonds_tmp%nbonds)=bl

         if (ii==jj) then
           bonds_tmp%nbonds= bonds_tmp%nbonds+1
         end if

       end if

     end do !! irshift

   end do !! jj
 end do !! ii

 call print_bonds(amu,bonds_tmp,natom,ntypat,symbol,typat,znucl)


!write(std_out,*) 'make_bonds 05'
!##########################################################
!### 05. Deallocate all the arrays inside bonds
!###     allocate them with the right size and fill them 

 if (associated(bonds%bond_vect))then
   ABI_DEALLOCATE(bonds%bond_vect)
   nullify(bonds%bond_vect)
 end if

 if (associated(bonds%bond_length))then
   ABI_DEALLOCATE(bonds%bond_length)
   nullify(bonds%bond_length)
 end if

 if (associated(bonds%nbondi))then
   ABI_DEALLOCATE(bonds%nbondi)
   nullify(bonds%nbondi)
 end if

 if (associated(bonds%indexi))then
   ABI_DEALLOCATE(bonds%indexi)
   nullify(bonds%indexi)
 end if

 bonds%nbonds=bonds_tmp%nbonds

 if (bonds%nbonds>0) then
!  Allocate the arrays with exactly the rigth nbonds
   ABI_ALLOCATE(bonds%bond_vect,(3,bonds%nbonds))
   ABI_ALLOCATE(bonds%bond_length,(bonds%nbonds))
   ABI_ALLOCATE(bonds%indexi,(natom,bonds%nbonds))
   ABI_ALLOCATE(bonds%nbondi,(natom))

!  Fill the values
   bonds%bond_vect(:,1:bonds%nbonds)=bonds_tmp%bond_vect(:,1:bonds%nbonds)
   bonds%bond_length(1:bonds%nbonds)=bonds_tmp%bond_length(1:bonds%nbonds)
   bonds%indexi(:,1:bonds%nbonds)=bonds_tmp%indexi(:,1:bonds%nbonds)
   bonds%nbondi(:)=bonds_tmp%nbondi(:)
 end if


 ABI_DEALLOCATE(bonds_tmp%bond_vect)
 ABI_DEALLOCATE(bonds_tmp%bond_length)
 ABI_DEALLOCATE(bonds_tmp%indexi)
 ABI_DEALLOCATE(bonds_tmp%nbondi)

 nullify(bonds_tmp%bond_vect)
 nullify(bonds_tmp%bond_length)
 nullify(bonds_tmp%indexi)
 nullify(bonds_tmp%nbondi)

end subroutine make_bonds_new
!!***
