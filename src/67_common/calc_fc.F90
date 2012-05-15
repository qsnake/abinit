!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_fc
!! NAME
!! calc_fc
!!
!! FUNCTION
!! calculation and output of Fermi-contact term at each atomic site
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (SM,VR,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  natom=number of atoms in cell.
!!  ntypat=number of atom types
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  typat(natom)=type (integer) for each atom
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      leave_new,make_fc_paw,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine calc_fc(natom,ntypat,pawrad,pawrhoij,pawtab,psps,typat)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_fc'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_66_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: typat(natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iatom
 character(len=500) :: message
!arrays
 type(nuclear_type) :: fc(natom)

! ************************************************************************

!DEBUG
!write(std_out,*)' calc_fc : enter'
!ENDDEBUG
!Compatibility tests
 if (psps%usepaw /= 1) then
   write (message,'(4a)')' calc_efg : ERROR- ',ch10,&
&   ' usepaw /= 1 but Fermi-contact calculation requires PAW ',ch10
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 do iatom = 1, natom
   ABI_ALLOCATE(fc(iatom)%spden,(pawrhoij(iatom)%nspden))
   fc(iatom)%spden(:) = zero
 end do

 call make_fc_paw(fc,natom,ntypat,pawrhoij,pawrad,pawtab,psps,typat)

 write(message,'(a,a,a)' ) ch10,' Fermi-contact Term Calculation ',ch10
 call wrtout(ab_out,message,'COLL')

 do iatom = 1, natom
   if (pawrhoij(iatom)%nspden == 2) then
     write(message,'(a,i3,a,i3,a,f12.4)') ' Atom ',iatom,', typat ',typat(iatom),': FC total = ',&
&     fc(iatom)%spden(1)+fc(iatom)%spden(2)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,i3,a,i3,a,f12.4)') ' Atom ',iatom,', typat ',typat(iatom),': FC up - down = ',&
&     fc(iatom)%spden(1)-fc(iatom)%spden(2)
     call wrtout(ab_out,message,'COLL')
   else
     write(message,'(a,i3,a,i3,a,f12.4)') ' Atom ',iatom,', typat ',typat(iatom),': FC = ',fc(iatom)%spden(1)
     call wrtout(ab_out,message,'COLL')
   end if
 end do

 write(message,'(3a)')ch10,ch10,ch10
 call wrtout(ab_out,message,'COLL')

 do iatom = 1, natom
   ABI_DEALLOCATE(fc(iatom)%spden)
 end do

!DEBUG
!write(std_out,*)' calc_fc : exit '
!stop
!ENDDEBUG

 end subroutine calc_fc
!!***
