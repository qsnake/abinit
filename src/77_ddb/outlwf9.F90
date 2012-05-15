!{\src2tex{textfont=tt}}
!!****f* ABINIT/outlwf9
!!
!! NAME
!! outlwf9
!!
!! FUNCTION
!! Open input file for the ppddb9 code, then
!! echoes the input information.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG,JCC,CL)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell = lattice vector lengths
!! iodyn = fortran record for writing output
!! msym = maximum number of symmetries
!! natom = number of atoms
!! nph1l = number of points in first list
!! nsym = actual number of symmetries
!! ntypat = number of atom types
!! rprim = lattice vectors
!! symrel = symmetry operations (real reduced coordinates)
!! typat = array of atom types
!! xred = reduced atom coordinates
!!
!! OUTPUT
!!  (only writing)
!!
!! SIDE EFFECTS
!! All the other arguments are inputs
!! and are written to the file.
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      mkphbs
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine outlwf9 (acell,iodyn,msym,natom,nph1l,nsym,ntypat,rprim,symrel,typat,xred)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outlwf9'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iodyn,msym,natom,nph1l,nsym,ntypat
!arrays
 integer,intent(in) :: symrel(3,3,msym),typat(natom)
 real(dp),intent(in) :: acell(3),rprim(3,3),xred(3,natom)

!Local variables -------------------------
!Set routine version number here:
!scalars
 integer :: ii

!*********************************************************************

 write(iodyn,'(4i5)') natom,nph1l,nsym,ntypat
 write(iodyn,'(a,a,3f12.8)') 'acell=',ch10,acell(:)
 write(iodyn,'(a,a,3f12.8,a,3f12.8,a,3f12.8)')&
& 'rprim=',ch10,rprim(:,1),ch10,rprim(:,2),ch10,rprim(:,3)
 write(iodyn,'(a,a,i5)') 'natom=',ch10,natom
 write(iodyn,'(a,a,i5)') 'ntypat=',ch10,ntypat
 write(iodyn,'(a)')' xred :'
 do ii=1,natom
   write(iodyn,'(i5,3f12.8)') typat(ii),xred(:,ii)
 end do
 write(iodyn,'(a)')' symrel :'
 do ii=1,nsym
   write(iodyn,'(9i4)') symrel(:,:,ii)
 end do
 write(iodyn,'(a,a,i5)') 'no._of_Q_points=',ch10,nph1l


end subroutine outlwf9
!!***
