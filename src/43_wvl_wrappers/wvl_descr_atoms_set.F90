!!****f* defs_wvltypes/wvl_descr_atoms_set
!!
!! NAME
!! wvl_descr_atoms_set
!!
!! FUNCTION
!! Defines wvl%atoms% data structure
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=unit cell length scales (bohr)
!! dtset <type(dataset_type)>=all input variables for this dataset
!!
!! OUTPUT
!! wvl <type(wvl_internal_type)>= wavelet type
!!                 | nat      =  number of atoms
!!                 | ntypes   =  number of species
!!                 | alat1    =  acell(1)
!!                 | alat2    =  acell(2)
!!                 | alat3    =  acell(3)
!!                 | iatype   =  types for atoms
!!                 | lfrztyp  =  flag for the movement of atoms.
!!                 | natpol   =  integer related to polarisation at the first step
!!
!! PARENTS
!!      gstate,wvl_memory
!!
!! CHILDREN
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_descr_atoms_set(acell, icoulomb, natom, ntypat, typat, wvl)

 use m_profiling

  use defs_basis
  use defs_wvltypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_descr_atoms_set'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in)                    :: icoulomb, natom, ntypat
  type(wvl_internal_type), intent(inout) :: wvl
  !arrays
  integer, intent(in)                    :: typat(natom)
  real(dp), intent(in)                   :: acell(3)

!Local variables-------------------------------
!scalars
#if defined HAVE_DFT_BIGDFT
  integer :: itype
#endif

! *********************************************************************

#if defined HAVE_DFT_BIGDFT
!We create the atoms_data structure from this dataset
!to be used later in BigDFT routines.
 if (icoulomb == 0) then
   wvl%atoms%geocode = 'P'
 else if (icoulomb == 1) then
   wvl%atoms%geocode = 'F'
 else if (icoulomb == 2) then
   wvl%atoms%geocode = 'S'
 end if
 write(wvl%atoms%units, "(A)") "Bohr"
 wvl%atoms%nat      =  natom
 wvl%atoms%ntypes   =  ntypat
 ABI_ALLOCATE(wvl%atoms%atomnames,(ntypat))
 do itype = 1, ntypat, 1
   write(wvl%atoms%atomnames(itype), "(A,I2)") "At. type", itype
 end do
 wvl%atoms%alat1    =  acell(1)
 wvl%atoms%alat2    =  acell(2)
 wvl%atoms%alat3    =  acell(3)
 ABI_ALLOCATE(wvl%atoms%iatype,(natom))
 wvl%atoms%iatype   = typat
!All atoms are free to move.
 ABI_ALLOCATE(wvl%atoms%ifrztyp,(natom))
 wvl%atoms%ifrztyp  =  0
!No chrage or polarisation on atoms in a first step.
 ABI_ALLOCATE(wvl%atoms%natpol,(natom))
 wvl%atoms%natpol   =  100

 nullify(wvl%atoms%psppar)
 nullify(wvl%atoms%npspcode)
 nullify(wvl%atoms%ixcpsp)
 nullify(wvl%atoms%iasctype)
 nullify(wvl%atoms%nzatom)
 nullify(wvl%atoms%nelpsp)
 nullify(wvl%atoms%amu)
 nullify(wvl%atoms%aocc)
 nullify(wvl%atoms%radii_cf)
 nullify(wvl%atoms%nlcc_ngv)
 nullify(wvl%atoms%nlcc_ngc)
 nullify(wvl%atoms%nlccpar)

#endif  
end subroutine wvl_descr_atoms_set
!!***
