!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawpolev
!! NAME
!! pawpolev
!!
!! FUNCTION
!! Compute the PAW term for polarization, named expected value term
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, PH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  natom=number of atoms in cell.
!!  ntypat = number of atom types
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat(natom)=type (integer) for each atom
!!
!! OUTPUT
!! pelev(3)= electronic polarisation. expectation value term (PAW only)
!!
!! PARENTS
!!      berryphase_new
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawpolev(natom,ntypat,pawrhoij,pawtab,pelev,typat)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawpolev'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(out) :: pelev(3)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: iatom,idir,irhoij,ispden,itypat,jrhoij,klmn
 real(dp) :: c1,ro_dlt
!arrays
 integer,dimension(3) :: idirindx = (/4,2,3/)
 real(dp) :: tsec(2)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(560,1,tsec)

!note that when vector r is expanded in real spherical harmonics, the factor
!sqrt(four_pi/three) appears, as in the following
!x = sqrt(four_pi/three)*r*S_{1,1}
!y = sqrt(four_pi/three)*r*S_{1,-1}
!z = sqrt(four_pi/three)*r*S_{1,0}
!
!the moments pawtab()%qijl do not include such normalization factors
!see pawinit.F90 for their definition and computation

 c1=sqrt(four_pi/three)

 pelev=zero
 do idir=1,3
   do iatom=1,natom
     itypat=typat(iatom)
     do ispden=1,pawrhoij(iatom)%nspden
       jrhoij=1
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(irhoij)
         ro_dlt=pawrhoij(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
         pelev(idir)=pelev(idir)+ro_dlt*c1*pawtab(itypat)%qijl(idirindx(idir),klmn)
         jrhoij=jrhoij+pawrhoij(iatom)%cplex
       end do
     end do
   end do
 end do

!write (80,*) 'pelev', pelev

 call timab(560,2,tsec)

 DBG_EXIT("COLL")

end subroutine pawpolev
!!***
