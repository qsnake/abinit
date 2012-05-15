!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_etsf_geo_put
!! NAME
!! abi_etsf_geo_put
!!
!!
!! FUNCTION
!!  Output system geometry to a file, using the ETSF I/O file format.
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
!!  filapp = character string giving the root to form the name of the GEO file
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | title = a description for the pseudo.
!!
!! OUTPUT
!!  Data written in file whose name is filapp//'-etsf.nc'
!!
!! PARENTS
!!      m_io_kss,outscfcv,pawmkaewf,sigma
!!
!! CHILDREN
!!      atmdata,etsf_io_data_write,etsf_io_low_error_to_str,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine abi_etsf_geo_put(dtset, filapp, psps)

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
#define ABI_FUNC 'abi_etsf_geo_put'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filapp
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays

!Local variables-------------------------------
#if defined HAVE_TRIO_ETSF_IO
 type(etsf_groups) :: group_folder
 type(etsf_geometry),target :: geo_folder
 logical :: lstat
 type(etsf_io_low_error) :: error
 integer, target :: spgroup
 character(len = etsf_io_low_error_len) :: errmess
#endif
!scalars
 integer :: i
 real(dp) :: amu,rcov
 character(len=500) :: message
 character(len=fnlen) :: filgeom
!arrays
 character(len=2),allocatable,target :: symbols(:)
 character(len=80),allocatable,target :: psp_desc(:),symbols_long(:)

! *************************************************************************

#if defined HAVE_TRIO_ETSF_IO
!Initialize filename
 filgeom=trim(filapp)//'-etsf.nc'
 write(message, '(a,a)' ) ' abi_etsf_geo_put : about to open file ',filgeom
 call wrtout(std_out,message,'COLL')

!Set-up atomic symbols
 ABI_ALLOCATE(symbols,(dtset%ntypat))
 ABI_ALLOCATE(symbols_long,(dtset%ntypat))
 ABI_ALLOCATE(psp_desc,(dtset%ntypat))
 do i=1, dtset%ntypat
   call atmdata(amu,rcov,symbols(i),dtset%znucl(i))
   write(symbols_long(i), "(A2,A78)") symbols(i), repeat(char(0), 78)
   write(psp_desc(i), "(A,A)") psps%title(i)(1:min(80, len_trim(psps%title(i)))), &
&   repeat(char(0), max(0, 80 - len_trim(psps%title(i))))
 end do

!Fill-in geometry folder
 if (dtset%spgroup > 0) then
   spgroup = dtset%spgroup
   geo_folder%space_group => spgroup
 end if
!To use rprimd or xred, add it as an argument.
!!! geo_folder%primitive_vectors => rprimd
!!! geo_folder%reduced_symmetry_matrices => dtset%symrel
!!! geo_folder%reduced_symmetry_translations => dtset%tnons
!!! geo_folder%atom_species => dtset%typat
!!! geo_folder%reduced_atom_positions => xred
!!! if (psps%npsp == psps%ntypat) then
!!!   geo_folder%valence_charges => psps%zionpsp
!!! end if
!!! geo_folder%atomic_numbers => dtset%znucl
 geo_folder%atom_species_names => symbols_long
 geo_folder%chemical_symbols => symbols
 geo_folder%pseudopotential_types => psp_desc

 group_folder%geometry => geo_folder

 call etsf_io_data_write(filgeom, group_folder, lstat, error)
 if (.not. lstat) then
   call etsf_io_low_error_to_str(errmess, error)
   write(message, "(A,A,A,A)") ch10, " abi_etsf_geo_put: ERROR -", ch10, errmess(1:min(470, len(errmess)))
   call wrtout(std_out, message, 'COLL')
 end if

!Free memory
 ABI_DEALLOCATE(symbols)
 ABI_DEALLOCATE(symbols_long)
 ABI_DEALLOCATE(psp_desc)
#endif

end subroutine abi_etsf_geo_put
!!***
