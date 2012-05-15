!!****m* ABINIT/interfaces_62_iowfdenpot
!! NAME
!! interfaces_62_iowfdenpot
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/62_iowfdenpot
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

module interfaces_62_iowfdenpot

 implicit none

interface
 subroutine WffReadEigK(eigen,formeig,headform,ikpt,isppol,mband,mpi_enreg,nband,tim_rwwf,wff)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: nband
  integer,intent(in) :: tim_rwwf
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out) :: eigen((2*mband)**formeig*mband)
 end subroutine WffReadEigK
end interface

interface
 subroutine WffReadSkipK(formeig,headform,ikpt,isppol,mpi_enreg,wff)
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
 end subroutine WffReadSkipK
end interface

interface
 subroutine calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,nunit,ratsph,rhor,rprimd,typat,ucvol,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: nunit
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: ratsph(ntypat)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine calcdensph
end interface

interface
 subroutine initwf(cg,eig_k,formeig,headform,icg,ikpt,ikptsp_old,&  
  &  isppol,mcg,mpi_enreg,&  
  &  nband_k,nkpt,npw,nspinor,occ_k,wff1)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(inout) :: ikptsp_old
  integer,intent(in) :: isppol
  integer,intent(in) :: mcg
  integer,intent(in) :: nband_k
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff1
  real(dp),intent(out) :: cg(2,mcg)
  real(dp),intent(out) :: eig_k((2*nband_k)**formeig*nband_k)
  real(dp),intent(inout) :: occ_k(nband_k)
 end subroutine initwf
end interface

interface
 subroutine ioarr(accessfil,arr,dtset,etotal,fform,fildata,hdr,mpi_enreg,&  
  &  ncplxfft,pawrhoij,rdwr,rdwrpaw,wvl)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: accessfil
  integer,intent(inout) :: fform
  integer,intent(in) :: ncplxfft
  integer,intent(in) :: rdwr
  integer,intent(in) :: rdwrpaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: etotal
  character(len=fnlen),intent(in) :: fildata
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wvl_internal_type), intent(in) :: wvl
  real(dp),intent(inout),target :: arr(ncplxfft,dtset%nspden)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*dtset%usepaw*rdwrpaw)
 end subroutine ioarr
end interface

interface
 subroutine out1dm(fnameabo_app_1dm,natom,nfft,ngfft,nspden,ntypat,&  
  &  rhor,rprimd,typat,ucvol,vtrial,xred,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  character(len=fnlen),intent(in) :: fnameabo_app_1dm
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vtrial(nfft,nspden)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine out1dm
end interface

interface
 subroutine outwant(dtfil,dtset,eig,cg,kg,npwarr,mband,mcg,mpi_enreg,nkpt,nsppol,&  
  mkmem,mpw,wff,prtwant)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer :: mband
  integer :: mcg
  integer :: mkmem
  integer :: mpw
  integer :: nkpt
  integer :: nsppol
  integer :: prtwant
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  real(dp) :: cg(2,mcg)
  real(dp) :: eig(mband*nkpt*nsppol)
  integer :: kg(3,mpw*mkmem)
  integer :: npwarr(nkpt)
 end subroutine outwant
end interface

interface
 subroutine randac(debug,headform1,ikptsp_prev,ikpt,isppol,&  
  &  nband,nkpt,nsppol,wffinp)
  use m_wffile
  implicit none
  integer,intent(in) :: debug
  integer,intent(in) :: headform1
  integer,intent(in) :: ikpt
  integer,intent(inout) :: ikptsp_prev
  integer,intent(in) :: isppol
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(wffile_type),intent(inout) :: wffinp
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine randac
end interface

interface
 subroutine rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,option,unitfile)
  implicit none
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(inout) :: nband_k
  integer,intent(inout) :: npw_k
  integer,intent(inout) :: nspinor
  integer,intent(in) :: option
  integer,intent(in) :: unitfile
 end subroutine rdnpw
end interface

interface
 subroutine read_wfrspa(state,fnametmp_tdwf,eigbnd,iband,isppol,imkmem,ndiel4,ndiel5,ndiel6,wfrspa_extract)
  use defs_basis
  implicit none
  integer,intent(out) :: iband
  integer,intent(in) :: imkmem
  integer,intent(out) :: isppol
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: state
  real(dp),intent(inout) :: eigbnd
  character(len=fnlen),intent(in) :: fnametmp_tdwf
  real(dp),intent(inout) :: wfrspa_extract(ndiel4,ndiel5,ndiel6)
 end subroutine read_wfrspa
end interface

end module interfaces_62_iowfdenpot
!!***
