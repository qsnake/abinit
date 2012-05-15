!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_etsf_electrons_put
!! NAME
!! abi_etsf_electrons_put
!!
!!
!! FUNCTION
!!  Output system of electrons to a file, using the ETSF I/O file format.
!!
!! COPYRIGHT
!! Copyright (C) 2006-2012 ABINIT group (Yann Pouillon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | natom  = number of atoms in unit cell
!!   | ntypat = number of types of atoms in unit cell.
!!   | typat(natom) = type integer for each atom in cell
!!   | znucl(ntypat)= real(dp), atomic number of atom type
!!  filapp = character string giving the root to form the name of the ETSF file
!!
!! OUTPUT
!!  Data written in file whose name is filapp//'-etsf.nc'
!!
!! PARENTS
!!      outscfcv,pawmkaewf,sigma
!!
!! CHILDREN
!!      etsf_io_data_write,etsf_io_low_error_to_str,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine abi_etsf_electrons_put(dtset, filapp)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_etsf_electrons_put'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filapp
 type(dataset_type),intent(in) :: dtset

!Local variables-------------------------------
#if defined HAVE_TRIO_ETSF_IO
 type(etsf_groups) :: group_folder
 type(etsf_electrons),target :: electrons_folder
 logical :: lstat
 type(etsf_io_low_error) :: error
 integer, target :: nelect
 real(dp), target :: tsmear
 character(len = etsf_charlen), target :: smearing
 character(len=etsf_io_low_error_len) :: errmess
#endif
!scalars
 character(len=500) :: message
 character(len=fnlen) :: filgeom

! *************************************************************************

#if defined HAVE_TRIO_ETSF_IO
!Initialize filename
 filgeom=trim(filapp)//'-etsf.nc'
 write(message, '(a,a)' ) ' abi_etsf_electrons_put : about to open file ',filgeom
 call wrtout(std_out,message,'COLL')

!Fill-in electrons folder
!FIXME ! define XC and smearing scheme
!MG WARNING, in abinit, unlike ETSF, nelect is real to allow for charging and alchemy!!
 nelect = dtset%nelect
 electrons_folder%number_of_electrons => nelect
 if (dtset%occopt == 3) then
   write(smearing, "(A)") "Fermi-Dirac"
 else if (dtset%occopt == 4) then
   write(smearing, "(A)") "cold smearing of N. Marzari with minimization of the bump"
 else if (dtset%occopt == 5) then
   write(smearing, "(A)") "cold smearing of N. Marzari with monotonic function in the tail"
 else if (dtset%occopt == 6) then
   write(smearing, "(A)") "Methfessel and Paxton"
 else if (dtset%occopt == 7) then
   write(smearing, "(A)") "gaussian"
 else if (dtset%occopt == 8) then
   write(smearing, "(A)") "uniform"
 else
   write(smearing, "(A)") "none"
 end if
 electrons_folder%smearing_scheme => smearing
 tsmear = dtset%tsmear
 electrons_folder%smearing_width => tsmear
 group_folder%electrons => electrons_folder

 call etsf_io_data_write(filgeom, group_folder, lstat, error)
 if (.not. lstat) then
   call etsf_io_low_error_to_str(errmess, error)
   write(message, "(A,A,A,A)") ch10, " abi_etsf_electrons_put: ERROR -", ch10, errmess(1:min(465, len(errmess)))
   call wrtout(std_out, message, 'COLL')
 end if
#endif

end subroutine abi_etsf_electrons_put
!!***
