!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_crystal_io
!! NAME
!! m_crystal_io
!!
!! FUNCTION
!! Module containing the methods used to do IO on Crystal objects. 
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG, YP, DC)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! TODO 
!! It would be useful to have a set of
!! methods to dump the crystalline structure according to the file format 
!! used by an external graphical software ...
!!
!! PARENTS
!!
!! CHILDREN
!!  
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_crystal_io

 use m_profiling

 use defs_basis
 use m_errors
 use m_crystal
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

 use defs_abitypes,    only : Hdr_type
 use m_numeric_tools,  only : set2unit

 implicit none

 private 
!!***

 public :: init_crystal_from_hdr   ! Initialize the object from the abinit header.
 public :: abi_crystal_put
 ! TODO corresponding reading method is missing, at present use InitCrystalFromHdr

CONTAINS

!!****f* m_crystal_io/init_crystal_from_hdr 
!! NAME
!!  init_crystal_from_hdr
!!
!! FUNCTION
!!  Initializes a crystal_structure data type starting from the abinit header.
!!
!! INPUTS
!!  Hdr<Hdr_type>=the abinit header
!!  timrev ==2 => take advantage of time-reversal symmetry
!!         ==1 ==> do not use time-reversal symmetry 
!!  remove_inv [optional]= if .TRUE. the inversion symmetry is removed from the set of operations
!!  even though it is present in the header
!!
!! OUTPUT
!!  Cryst<crystal_structure>= the data type filled with data reported in the abinit header 
!!
!! TODO
!!  Add information on the use of time-reversal in the Abinit header.
!!
!! PARENTS
!!      m_screening,m_sigma_results,mlwfovlp_qp,mrgscr,rdm,setup_bse
!!      setup_screening,setup_sigma
!!
!! CHILDREN
!!      atmdata,etsf_io_data_write,etsf_io_low_close,etsf_io_low_def_var
!!      etsf_io_low_error_to_str,etsf_io_low_open_modify
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!      wrtout
!!
!! SOURCE

subroutine init_crystal_from_hdr(Cryst,Hdr,timrev,remove_inv)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_crystal_from_hdr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(Hdr_type),intent(in) :: Hdr
 type(crystal_structure),intent(out) :: Cryst 
 integer,intent(in) :: timrev
 logical,optional,intent(in) :: remove_inv

!Local variables-------------------------------
 integer :: space_group
 logical :: rinv,ltest,use_antiferro
! *********************************************************************

 rinv=.FALSE.; if (PRESENT(remove_inv)) rinv=remove_inv
 use_antiferro=(Hdr%nspden==2.and.Hdr%nsppol==1)

 ! === consistency check ===
 ltest = (timrev==1.or.timrev==2)
 ABI_CHECK(ltest,"Wrong value for timrev (1|2)")
 if (use_antiferro) then
  ABI_CHECK(ANY(Hdr%symafm==-1),"Wrong nspden, nsppol, symafm.")
 end if

 space_group=0 !FIXME not known

 call init_crystal(Cryst,space_group,Hdr%natom,Hdr%npsp,Hdr%ntypat,Hdr%nsym,Hdr%rprimd,Hdr%typat,Hdr%xred,&
& Hdr%zionpsp,Hdr%znuclpsp,timrev,use_antiferro,rinv,Hdr%title,&
& Hdr%symrel,Hdr%tnons,Hdr%symafm) ! Optional

end subroutine init_crystal_from_hdr
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/abi_crystal_put
!! NAME
!! abi_crystal_put
!!
!! FUNCTION
!! Output system geometry to a file, using the NETCDF file format and ETSF I/O.
!! Data are taken from the crystal_structure object.
!!
!! INPUTS
!!  Cryst<crystal_structure>=Object defining the unit cell and its symmetries.
!!     %natom=Number of atoms in the unit cell
!!     %ntypat=Number of types of atoms in the unit cell.
!!     %npsp=Number of pseudopotentials.
!!     %space_group=Space group identifier.
!!     %typat(natom)=Type of each atom in the unit cell.
!!     %ziontypat(ntypat)=Charge of the pseudo-ion (No of valence electrons needed to screen exactly the pseudopotential).
!!     %znucl(ntypat)=Nuclear charge for each type of pseudopotential.
!!     %title(npsps)=The content of first line read from the psp file.
!!     %rprimd(3,3)=Real space dimensional primitive translations (bohr)
!!     %xred(3,natom)=Reduced coordinates of atoms.
!!     %symrel=Symmetry operations in real space (in terms of rprimd)
!!     %tnons(3,nsym)=Fractional translations.
!!     %symafm(nsym)=(Anti)ferromagnetic symmetries.
!!  fname=character string giving the root to form the name of the GEO file
!!
!! OUTPUT
!!  Data written in file whose name is fname.
!!
!! NOTES
!!  alchemy not treated, Crystal should be correctly Initialized at the beginning of the run.
!!
!! PARENTS
!!      m_sigma_results
!!
!! CHILDREN
!!      atmdata,etsf_io_data_write,etsf_io_low_close,etsf_io_low_def_var
!!      etsf_io_low_error_to_str,etsf_io_low_open_modify
!!      etsf_io_low_set_write_mode,etsf_io_low_write_dim,etsf_io_low_write_var
!!      wrtout
!!
!! SOURCE
subroutine abi_crystal_put(Cryst,fname)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_crystal_put'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: fname
 type(crystal_structure),target,intent(in) :: Cryst
!arrays

!Local variables-------------------------------
!scalars
 integer :: itypat,ungeo
 integer,target :: space_group
 real(dp) :: amu,rcov
 logical :: lstat,ltest
 character(len=500) :: msg
#if defined HAVE_TRIO_ETSF_IO
 character(len=etsf_io_low_error_len) :: errmess
 character(len=etsf_charlen) :: symmorphic
 !type(ETSF_dims) :: Dims
 type(ETSF_groups) :: Group_folder
 type(ETSF_geometry),target :: Geo_folder
 type(ETSF_io_low_error) :: Error_data
#endif
!arrays
 character(len=2),allocatable,target :: symbols(:)
 character(len=80),allocatable,target :: psp_desc(:),symbols_long(:)

! *************************************************************************

#if defined HAVE_TRIO_ETSF_IO

 !@crystal_structure

 ! === Open the file ===
 write(msg,'(2a)')' abi_etsf_geo_put : about to open file ',TRIM(fname)
 call wrtout(std_out,msg,'COLL')

 call etsf_io_low_open_modify(ungeo,fname,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 ! === Define subset of ETSF-Dimensions and write them on file ===
 ! * Use low level procedures to write dimensions, since 
 !   etsf_io_dims_def correctly complains that some ETSF dims are missing.
 ! * If dimensions already exist, check that definitions are coherent.
 !Dims%number_of_atoms               = Cryst%natom
 !Dims%number_of_atom_species        = Cryst%ntypat
 !Dims%number_of_symmetry_operations = Cryst%nsym
 !call etsf_io_dims_def(ungeo,Dims,lstat,Error_data)
 !if (.not.lstat) goto 1000

 call etsf_io_low_write_dim(ungeo,'number_of_atoms',Cryst%natom,lstat,Error_data=Error_data)  
 if (.not.lstat) goto 1000  
 call etsf_io_low_write_dim(ungeo,'number_of_atoms_species',Cryst%ntypat,lstat,Error_data=Error_data)  
 if (.not.lstat) goto 1000  
 call etsf_io_low_write_dim(ungeo,'number_of_symmetry_operations',Cryst%nsym,lstat,Error_data=Error_data)  
 if (.not.lstat) goto 1000  

 ! === Close the file, due to call to etsf_io_data_write  ===
 call etsf_io_low_close(ungeo,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 ! === Fill-in ETSF geometry folder ===

 ! FIXME alchemy not treated since Cryst should be initialized in invars2
 ltest=(Cryst%npsp==Cryst%ntypat)
 ABI_CHECK(ltest,"alchemy not supported")

 ! === Set-up atomic symbols ===
 ABI_ALLOCATE(symbols,(Cryst%ntypat))
 ABI_ALLOCATE(symbols_long,(Cryst%ntypat))
 ABI_ALLOCATE(psp_desc,(Cryst%ntypat))
 do itypat=1,Cryst%ntypat  
  call atmdata(amu,rcov,symbols(itypat),Cryst%znucl(itypat)) 
  write(symbols_long(itypat),'(a2,a78)') symbols(itypat),REPEAT(CHAR(0),78)
  write(psp_desc(itypat),'(2a)') &
&  Cryst%title(itypat)(1:MIN(80,LEN_TRIM(Cryst%title(itypat)))),REPEAT(CHAR(0),MAX(0,80-LEN_TRIM(Cryst%title(itypat))))
 end do

 space_group=0; if (Cryst%space_group>0) space_group=Cryst%space_group
 Geo_folder%space_group                   => space_group
 Geo_folder%primitive_vectors             => Cryst%rprimd
 Geo_folder%reduced_symmetry_matrices     => Cryst%symrel
 Geo_folder%reduced_symmetry_translations => Cryst%tnons
 Geo_folder%atom_species                  => Cryst%typat
 Geo_folder%reduced_atom_positions        => Cryst%xred
 if (Cryst%npsp==Cryst%ntypat) then
  Geo_folder%valence_charges              => Cryst%ziontypat
 end if
 Geo_folder%atomic_numbers                => Cryst%znucl
 Geo_folder%atom_species_names            => symbols_long
 Geo_folder%chemical_symbols              => symbols
 Geo_folder%pseudopotential_types         => psp_desc

 symmorphic='no'; if (isymmorphic(Cryst)) symmorphic='yes'

 Group_folder%geometry => Geo_folder

 call etsf_io_data_write(fname,Group_folder,lstat,Error_data)
 if (.not.lstat) goto 1000

 ABI_DEALLOCATE(symbols)
 ABI_DEALLOCATE(symbols_long)
 ABI_DEALLOCATE(psp_desc)
 !
 ! === At this point we have an ETSF-compliant file ===
 ! * Add additional stuff for internal use in abinit.
 ! TODO add spinat.

 call etsf_io_low_open_modify(ungeo,fname,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 call etsf_io_low_def_var(ungeo,'symafm',etsf_io_low_integer,&
& (/'number_of_symmetry_operations'/),lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 call etsf_io_low_set_write_mode(ungeo,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 call etsf_io_low_write_var(ungeo,'symafm',Cryst%symafm,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 call etsf_io_low_close(ungeo,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 ! === Handle the error ===
 1000 continue
 if (.not.lstat) then
  call etsf_io_low_error_to_str(errmess,Error_data)
  msg = errmess(1:min(500,len(errmess)))
  MSG_ERROR(msg)
 end if

#else
 msg = 'ETSF-IO support is not activated.'
 MSG_ERROR(msg)
#endif

end subroutine abi_crystal_put
!!***

end MODULE m_crystal_io
