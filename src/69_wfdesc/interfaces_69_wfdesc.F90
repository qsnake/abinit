!!****m* ABINIT/interfaces_69_wfdesc
!! NAME
!! interfaces_69_wfdesc
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/69_wfdesc
!!
!! COPYRIGHT
!! Copyright (C) 2010-2011 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_69_wfdesc

 implicit none

interface
 subroutine classify_bands(Wfd,use_paw_aeur,first_band,last_band,ik_ibz,spin,ngfftf,&  
  &  Cryst,BSt,Pawtab,Pawrad,Pawang,Psps,tolsym,BSym,&  
  &  EDIFF_TOL) ! optional
  use defs_basis
  use m_bands_sym
  use m_crystal
  use defs_datatypes
  use m_wfs
  implicit none
  integer,intent(in) :: first_band
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: last_band
  integer,intent(in) :: spin
  type(bandstructure_type),intent(in) :: BSt
  type(bands_symmetries),intent(out) :: BSym
  type(crystal_structure),intent(in) :: Cryst
  real(dp),intent(in),optional :: EDIFF_TOL
  type(pawang_type),intent(in) :: Pawang
  type(pseudopotential_type),intent(in) :: Psps
  type(wfs_descriptor),intent(inout) :: Wfd
  real(dp),intent(in) :: tolsym
  logical,intent(in) :: use_paw_aeur
  integer,intent(in) :: ngfftf(18)
  type(pawrad_type),intent(inout) :: Pawrad(Cryst%ntypat*Wfd%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 end subroutine classify_bands
end interface

interface
 subroutine rotate_cprj(kpoint,isym,nspinor,nbnds,natom,nsym,typat,indsym,Cprj_in,Cprj_out)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: isym
  integer,intent(in) :: natom
  integer,intent(in) :: nbnds
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsym
  type(cprj_type),intent(in) :: Cprj_in(natom,nspinor*nbnds)
  type(cprj_type),intent(out) :: Cprj_out(natom,nspinor*nbnds)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: kpoint(3)
  integer,intent(in) :: typat(natom)
 end subroutine rotate_cprj
end interface

interface
 function paw_phirotphj(nspinor,natom,typat,zarot_isym,Pawtab,Psps,Cprj_b1,Cprj_b2,conjg_left) result(omat)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nspinor
  type(pseudopotential_type),intent(in) :: Psps
  logical,optional,intent(in) :: conjg_left
  type(cprj_type),intent(in) :: Cprj_b1(natom,nspinor)
  type(cprj_type),intent(in) :: Cprj_b2(natom,nspinor)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
  real(dp) :: omat(2)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: zarot_isym(:,:,:)
 end function paw_phirotphj
end interface

interface
 subroutine outkss(Dtfil,Dtset,ecut,gmet,gprimd,Hdr,&  
  &  kssform,mband,mcg,mcprj,mgfft,mkmem,MPI_enreg,mpsang,mpw,natom,&  
  &  nfft,nkpt,npwarr,nspden,nsppol,nsym,ntypat,occ,Pawtab,Pawfgr,Paw_ij,&  
  &  prtvol,Psps,rprimd,Wffnow,vtrial,xred,cg,usecprj,Cprj,eigen,ierr)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: kssform
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  integer,intent(in) :: usecprj
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(in) :: Dtset
  type(hdr_type),intent(inout) :: Hdr
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pawfgr_type), intent(in) :: Pawfgr
  type(pseudopotential_type),intent(in) :: Psps
  type(wffile_type),intent(inout) :: Wffnow
  real(dp),intent(in) :: ecut
  type(cprj_type),intent(in) :: Cprj(natom,mcprj*usecprj)
  type(paw_ij_type),intent(in) :: Paw_ij(natom*Psps%usepaw)
  type(pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in),target :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: vtrial(nfft,nspden)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine outkss
end interface

interface
 subroutine wfd_mkrho(Wfd,Cryst,Psps,Kmesh,Bstr,ngfftf,nfftf,rhor,&  
  &  optcalc) ! optional arguments
  use m_bz_mesh
  use m_crystal
  use defs_datatypes
  use defs_basis
  use m_wfs
  implicit none
  integer,intent(in) :: nfftf
  integer,intent(in),optional :: optcalc
  type(bandstructure_type),intent(in) :: Bstr
  type(crystal_structure),intent(in) :: Cryst
  type(bz_mesh_type),intent(in) :: Kmesh
  type(pseudopotential_type),intent(in) :: Psps
  type(wfs_descriptor),intent(inout) :: Wfd
  integer,intent(in) :: ngfftf(18)
  real(dp),intent(out) :: rhor(nfftf,Wfd%nspden)
 end subroutine wfd_mkrho
end interface

interface
 subroutine test_charge(nfftf,nelectron_exp,nspden,rhor,ucvol,&  
  &  usepaw,usexcnhat,usefinegrid,compch_sph,compch_fft,omegaplasma)
  use defs_basis
  implicit none
  integer,intent(in) :: nfftf
  integer,intent(in) :: nspden
  integer,intent(in) :: usefinegrid
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  real(dp),intent(in) :: compch_fft
  real(dp),intent(in) :: compch_sph
  real(dp),intent(in) :: nelectron_exp
  real(dp),intent(out) :: omegaplasma
  real(dp),intent(in) :: ucvol
  real(dp),intent(inout) :: rhor(nfftf,nspden)
 end subroutine test_charge
end interface

interface
 subroutine wfd_pawrhoij(Wfd,Cryst,Bst,kptopt,pawrhoij,pawprtvol)
  use m_crystal
  use defs_datatypes
  use m_wfs
  implicit none
  integer,intent(in) :: kptopt
  integer,intent(in) :: pawprtvol
  type(bandstructure_type),intent(in) :: Bst
  type(crystal_structure),intent(in) :: Cryst
  type(wfs_descriptor),intent(inout) :: Wfd
  type(pawrhoij_type),intent(inout) :: pawrhoij(Wfd%natom)
 end subroutine wfd_pawrhoij
end interface

interface
 subroutine wfd_vnlpsi(Wfd,band,ik_ibz,spin,npw_k,Cryst,Psps,GS_hamk,vnl_psi,opaw_psi,Kext)
  use defs_basis
  use m_crystal
  use defs_datatypes
  use m_wfs
  implicit none
  integer,intent(in) :: band
  integer,intent(in) :: ik_ibz
  integer,intent(in) :: npw_k
  integer,intent(in) :: spin
  type(crystal_structure),intent(in) :: Cryst
  type(gs_hamiltonian_type),intent(in) :: GS_hamk
  type(kdata_t),optional,intent(in) :: Kext
  type(pseudopotential_type),intent(in) :: Psps
  type(wfs_descriptor),intent(inout) :: Wfd
  real(dp),intent(out) :: opaw_psi(2,npw_k*Wfd%nspinor*Wfd%usepaw)
  real(dp),intent(out) :: vnl_psi(2,npw_k*Wfd%nspinor)
 end subroutine wfd_vnlpsi
end interface

end module interfaces_69_wfdesc
!!***
