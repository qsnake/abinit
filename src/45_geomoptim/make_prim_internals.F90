!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_prim_internals
!! NAME
!! make_prim_internals
!!
!! FUNCTION
!!  Determine the bonds, angles and dihedrals for a starting
!!  geometry, based on covalent radii for the atoms.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! natom  = Number of atoms (dtset%natom)
!! icenter= index of the center of the number of shifts
!! nrshift= dimension of rshift
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! rshift(3,nrshift)=shift in xred that must be done to find all neighbors of
!!                   a given atom within a given number of neighboring shells
!! xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!!   deloc <type(ab_delocint)>=Important variables for
!!   |                           pred_delocint
!!   |
!!   | nang     = Number of angles
!!   | nbond    = Number of bonds
!!   | ncart    = Number of cartesian directions
!!   |             (used for constraints)
!!   | ndihed   = Number of dihedrals
!!   | nrshift  = Dimension of rshift
!!   | ninternal= Number of internal coordinates
!!   |            ninternal=nbond+nang+ndihed+ncart
!!   |
!!   | angs(2,3,nang)  = Indexes to characterize angles
!!   | bonds(2,2,nbond)= For a bond between iatom and jatom
!!   |                   bonds(1,1,nbond) = iatom
!!   |                   bonds(2,1,nbond) = icenter
!!   |                   bonds(1,2,nbond) = jatom
!!   |                   bonds(2,2,nbond) = irshift
!!   | carts(2,ncart)  = Index of total primitive internal,
!!   |                   and atom (carts(2,:))
!!   | dihedrals(2,4,ndihed)= Indexes to characterize dihedrals
!!   |
!!   | rshift(3,nrshift)= Shift in xred that must be done to find
!!   |                    all neighbors of a given atom within a
!!   |                    given number of neighboring shells
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!!   Adds cartesian coordinates if the number of internals with a
!!   given atom is < 4 the chosen coordinate could be optimized
!!   to be less dependent of the internals already incorporated.
!!
!! PARENTS
!!      pred_delocint
!!
!! CHILDREN
!!      make_angles,make_bonds,make_dihedrals
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine make_prim_internals(deloc,icenter,natom,ntypat,rprimd,&
& typat,xcart,znucl)

 use m_profiling

 use defs_basis
 use defs_mover
 use m_delocint

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_prim_internals'
 use interfaces_45_geomoptim, except_this_one => make_prim_internals
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ab_delocint),intent(inout) :: deloc
 integer,intent(in) :: icenter,natom,ntypat
!arrays
 real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: znucl(:) ! znucl(ntypat) or
                                 ! znucl(npsp) ?

!Local variables ------------------------------
! function
!scalars
 integer :: iang,iatom,ibond,icart,idihed,ii
 real(dp) :: spp
!arrays
 integer :: particip_atom(natom)
 integer,allocatable :: badangles(:)
 real(dp) :: rpt1(3),rpt3(3)
!DEBUG
!real(dp) :: rpt2(3)
!ENDDEBUG

!no_abirules

!************************************************************************

!DEBUG
!write(std_out,*) 'make_deloc_internals: enter'
!ENDDEBUG
 particip_atom(:) = 0

 call make_bonds(deloc,natom,ntypat,icenter,rprimd,typat,xcart,znucl)

!DEBUG
!write(std_out,*) size(deloc%bonds,1)
!write(std_out,*) size(deloc%bonds,2)
!write(std_out,*) size(deloc%bonds,3)
!ENDDEBUG

 do ibond=1,deloc%nbond
   write(std_out,'(a,i4,2(2i5,2x))') 'bond ', ibond, deloc%bonds(:,:,ibond)
   particip_atom(deloc%bonds(1,1,ibond)) = particip_atom(deloc%bonds(1,1,ibond))+1
   particip_atom(deloc%bonds(1,2,ibond)) = particip_atom(deloc%bonds(1,2,ibond))+1
 end do

!DEBUG
!return
!ENDDEBUG

 call make_angles(deloc,icenter,natom)
 ABI_ALLOCATE(badangles,(deloc%nang))
 badangles(:) = 0
 do iang=1,deloc%nang
   write(std_out,'(a,i4,3(2i5,2x))') 'angle ', iang, deloc%angs(:,:,iang)
   particip_atom(deloc%angs(1,1,iang)) = particip_atom(deloc%angs(1,1,iang))+1
   particip_atom(deloc%angs(1,2,iang)) = particip_atom(deloc%angs(1,2,iang))+1
   particip_atom(deloc%angs(1,3,iang)) = particip_atom(deloc%angs(1,3,iang))+1

!  DEBUG
!  rpt1(:) = xcart(:,deloc%angs(1,1,iang)) &
!  & + deloc%rshift(1,deloc%angs(2,1,iang))*rprimd(:,1) &
!  & + deloc%rshift(2,deloc%angs(2,1,iang))*rprimd(:,2) &
!  & + deloc%rshift(3,deloc%angs(2,1,iang))*rprimd(:,3)
!  rpt2(:) = xcart(:,deloc%angs(1,2,iang)) &
!  & + deloc%rshift(1,deloc%angs(2,2,iang))*rprimd(:,1) &
!  & + deloc%rshift(2,deloc%angs(2,2,iang))*rprimd(:,2) &
!  & + deloc%rshift(3,deloc%angs(2,2,iang))*rprimd(:,3)
!  rpt3(:) = xcart(:,deloc%angs(1,3,iang)) &
!  & + deloc%rshift(1,deloc%angs(2,3,iang))*rprimd(:,1) &
!  & + deloc%rshift(2,deloc%angs(2,3,iang))*rprimd(:,2) &
!  & + deloc%rshift(3,deloc%angs(2,3,iang))*rprimd(:,3)
!  write(std_out,*) rpt1,rpt2,rpt3,bond_length(rpt1,rpt2),bond_length(rpt2,rpt3)
!  ENDDEBUG

!  check if angles are 180 degrees: discard the dihedrals in that case.
   rpt1(:) = xcart(:,deloc%angs(1,1,iang)) &
&   + deloc%rshift(1,deloc%angs(2,1,iang))*rprimd(:,1) &
&   + deloc%rshift(2,deloc%angs(2,1,iang))*rprimd(:,2) &
&   + deloc%rshift(3,deloc%angs(2,1,iang))*rprimd(:,3) &
&   - xcart(:,deloc%angs(1,2,iang)) &
&   - deloc%rshift(1,deloc%angs(2,2,iang))*rprimd(:,1) &
&   - deloc%rshift(2,deloc%angs(2,2,iang))*rprimd(:,2) &
&   - deloc%rshift(3,deloc%angs(2,2,iang))*rprimd(:,3)

   rpt3(:) = xcart(:,deloc%angs(1,3,iang)) &
&   + deloc%rshift(1,deloc%angs(2,3,iang))*rprimd(:,1) &
&   + deloc%rshift(2,deloc%angs(2,3,iang))*rprimd(:,2) &
&   + deloc%rshift(3,deloc%angs(2,3,iang))*rprimd(:,3) &
&   - xcart(:,deloc%angs(1,2,iang)) &
&   - deloc%rshift(1,deloc%angs(2,2,iang))*rprimd(:,1) &
&   - deloc%rshift(2,deloc%angs(2,2,iang))*rprimd(:,2) &
&   - deloc%rshift(3,deloc%angs(2,2,iang))*rprimd(:,3)
   spp = (rpt1(1)*rpt3(1)+rpt1(2)*rpt3(2)+rpt1(3)*rpt3(3))&
&   / sqrt(rpt1(1)*rpt1(1)+rpt1(2)*rpt1(2)+rpt1(3)*rpt1(3)) &
&   / sqrt(rpt3(1)*rpt3(1)+rpt3(2)*rpt3(2)+rpt3(3)*rpt3(3))
   if (abs(abs(spp) - one) < tol6) then
     write(std_out,*) 'make_prim_internals : an angle is too close to 180 degrees:'
     write(std_out,*) '   will discard dihedrals using it '
     badangles(iang) = 1
   end if
 end do


 call make_dihedrals(badangles,deloc,icenter)
 do idihed=1,deloc%ndihed
   write(std_out,'(a,i4,4(2i5,2x))') 'dihedral ', idihed, deloc%dihedrals(:,:,idihed)
   particip_atom(deloc%dihedrals(1,1,idihed)) = particip_atom(deloc%dihedrals(1,1,idihed))+1
   particip_atom(deloc%dihedrals(1,2,idihed)) = particip_atom(deloc%dihedrals(1,2,idihed))+1
   particip_atom(deloc%dihedrals(1,3,idihed)) = particip_atom(deloc%dihedrals(1,3,idihed))+1
   particip_atom(deloc%dihedrals(1,4,idihed)) = particip_atom(deloc%dihedrals(1,4,idihed))+1

!  DEBUG
!  do ii=1,4
!  write(std_out,'((3E16.6,2x))') xcart(:,deloc%dihedrals(1,ii,idihed)) + &
!  &  deloc%rshift(1,deloc%dihedrals(2,ii,idihed))*rprimd(:,1)   + &
!  &  deloc%rshift(2,deloc%dihedrals(2,ii,idihed))*rprimd(:,2)   + &
!  &  deloc%rshift(2,deloc%dihedrals(2,ii,idihed))*rprimd(:,3)
!  end do
!  ENDDEBUG

 end do

 write(std_out,*) 'make_deloc_internals: nbond,nang,ndihed = ', deloc%nbond,deloc%nang,deloc%ndihed

!Check all atoms participate in at least 4 primitives. Otherwise, we should
!probably add cartesian coordinates to the internal ones.
 deloc%ncart = 0
 do iatom=1,natom
   if (particip_atom(iatom) < 4) then
     write(std_out,*) ' make_prim_internals : Warning : atom ', iatom, &
&     ' does not belong to enough primitives to determine its'
     write(std_out,*) ' position uniquely ! instead : ', particip_atom(iatom)
     write(std_out,*) ' Will add cartesian coordinates to set of internals.'
!    write(std_out,*) ' Not done yet.'
!    stop
     deloc%ncart = deloc%ncart + 4-particip_atom(iatom)
   end if
 end do
 if (associated(deloc%carts)) nullify(deloc%carts)
 ABI_ALLOCATE(deloc%carts ,(2,deloc%ncart))
 icart = 0
 do iatom=1,natom
   if (particip_atom(iatom) < 4) then
!    
!    kind of arbitrary : include first few directions for the atom:
!    x, then y then z
!    
     do ii=1,4-particip_atom(iatom)
       icart = icart+1
       deloc%carts(1,icart) = ii
       deloc%carts(2,icart) = iatom
     end do
   end if
 end do

!ninternal=nbond+nang+ndihed
 deloc%ninternal=deloc%nbond+deloc%nang+deloc%ndihed+deloc%ncart

end subroutine make_prim_internals

!-----------------------------------------------------------------------------
!!***
