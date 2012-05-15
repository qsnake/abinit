!!****m* ABINIT/interfaces_59_io_mpi
!! NAME
!! interfaces_59_io_mpi
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/59_io_mpi
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

module interfaces_59_io_mpi

 implicit none

interface
 subroutine chkexi(cpus,filnam,iexit,iout,mpi_enreg,openexit)
  use defs_basis
  use defs_abitypes
  implicit none
  integer          ,intent(out) :: iexit
  integer          ,intent(in) :: iout
  integer          ,intent(in) :: openexit
  real(dp)         ,intent(in) :: cpus
  character(len=fnlen),intent(in) :: filnam
  type(mpi_type)   ,intent(inout) :: mpi_enreg
 end subroutine chkexi
end interface

interface
 subroutine filnam_comm(nfilnam,filnam)
  use defs_basis
  implicit none
  integer, intent(in) :: nfilnam
  character(len=fnlen), intent(inout) :: filnam(nfilnam)
 end subroutine filnam_comm
end interface

interface
 subroutine handle_ncerr(ncerr,message)
  implicit none
  integer         ,intent(in) :: ncerr
  character(len=*),intent(in) :: message
 end subroutine handle_ncerr
end interface

interface
 subroutine hdr_check(fform,fform0,hdr,hdr0,mode_paral,restart,restartpaw)
  use defs_abitypes
  implicit none
  integer,intent(in) :: fform
  integer,intent(in) :: fform0
  integer,intent(out) :: restart
  integer,intent(out) :: restartpaw
  type(hdr_type),intent(in) :: hdr
  type(hdr_type),intent(in) :: hdr0
  character(len=4),intent(in) :: mode_paral
 end subroutine hdr_check
end interface

interface
 subroutine hdr_comm(hdr,master,me,spaceComm)
  use defs_abitypes
  implicit none
  integer, intent(in) :: master
  integer, intent(in) :: me
  integer, intent(in) :: spaceComm
  type(hdr_type),intent(inout) :: hdr
 end subroutine hdr_comm
end interface


!Generic interface of the routines hdr_io
interface hdr_io
 subroutine hdr_io_wfftype(fform,hdr,rdwr,wff)
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(inout) :: fform
  integer,intent(in) :: rdwr
  type(hdr_type),intent(inout) :: hdr
  type(wffile_type),intent(inout) :: wff
 end subroutine hdr_io_wfftype
 subroutine hdr_io_int(fform,hdr,rdwr,unitfi)
  use defs_abitypes
  implicit none
  integer,intent(inout) :: fform
  integer,intent(in) :: rdwr
  integer,intent(in) :: unitfi
  type(hdr_type),intent(inout) :: hdr
 end subroutine hdr_io_int
end interface
!End of the generic interface of hdr_io

interface
 subroutine hdr_io_etsf(fform,hdr,rdwr,unitwff)
  use defs_abitypes
  implicit none
  integer,intent(inout) :: fform
  integer,intent(in) :: rdwr
  integer,intent(in) :: unitwff
  type(hdr_type),intent(inout) :: hdr
 end subroutine hdr_io_etsf
end interface


!Generic interface of the routines hdr_io_netcdf
interface hdr_io_netcdf
 subroutine hdr_io_netcdf_wfftype(fform,hdr,rdwr,wff)
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(inout) :: fform
  integer,intent(in) :: rdwr
  type(hdr_type),intent(inout) :: hdr
  type(wffile_type),intent(inout) :: wff
 end subroutine hdr_io_netcdf_wfftype
 subroutine hdr_io_netcdf_int(fform,hdr,rdwr,unitfi)
  use defs_abitypes
  implicit none
  integer,intent(inout) :: fform
  integer,intent(in) :: rdwr
  integer,intent(in) :: unitfi
  type(hdr_type),intent(inout) :: hdr
 end subroutine hdr_io_netcdf_int
end interface
!End of the generic interface of hdr_io_netcdf


!Generic interface of the routines hdr_skip
interface hdr_skip
 subroutine hdr_skip_int(unitfi,ierr)
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: unitfi
 end subroutine hdr_skip_int
 subroutine hdr_skip_wfftype(wff,ierr)
  use m_wffile
  implicit none
  integer, intent(out) :: ierr
  type(wffile_type),intent(inout) :: wff
 end subroutine hdr_skip_wfftype
end interface
!End of the generic interface of hdr_skip

interface
 subroutine hdr_update(bantot,etot,fermie,hdr,natom,residm,rprimd,occ,mpi_enreg,pawrhoij,usepaw,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: bantot
  integer,intent(in) :: natom
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: etot
  real(dp),intent(in) :: fermie
  type(hdr_type),intent(out) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: residm
  real(dp),intent(in) :: occ(bantot)
  type(pawrhoij_type),intent(in) :: pawrhoij(:)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine hdr_update
end interface

interface
 subroutine hdr_vs_dtset(Hdr,Dtset)
  use defs_abitypes
  implicit none
  type(dataset_type),intent(inout) :: Dtset
  type(hdr_type),intent(in) :: Hdr
 end subroutine hdr_vs_dtset
end interface

interface
 subroutine outxfhist(ab_xfh,natom,option,wff2,ios)
  use defs_mover
  use m_wffile
  implicit none
  integer          ,intent(out) :: ios
  integer          ,intent(in) :: natom
  integer          ,intent(in) :: option
  type(ab_xfh_type),intent(inout) :: ab_xfh
  type(wffile_type),intent(inout) :: wff2
 end subroutine outxfhist
end interface

interface
 subroutine rwwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&  
  &  nband,nband_disk,npw,nspinor,occ,option,optkg,tim_rwwf,wff)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: nband
  integer,intent(inout) :: nband_disk
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: option
  integer,intent(in) :: optkg
  integer,intent(in) :: tim_rwwf
  type(mpi_type), intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(inout),target :: cg(2,mcg)
  real(dp),intent(inout),target :: eigen((2*mband)**formeig*mband)
  integer,intent(inout),target :: kg_k(3,optkg*npw)
  real(dp),intent(inout),target :: occ(mband)
 end subroutine rwwf
end interface

interface
 subroutine readwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&  
  &  nband,nband_disk,npw,nspinor,occ,option,optkg,wff)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: nband
  integer,intent(inout) :: nband_disk
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: option
  integer,intent(in) :: optkg
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(inout),target :: cg(2,mcg)
  real(dp),intent(inout),target :: eigen((2*mband)**formeig*mband)
  integer,intent(inout),target :: kg_k(3,optkg*npw)
  real(dp),intent(inout),target :: occ(mband)
 end subroutine readwf
end interface

interface
 subroutine writewf(cg,eigen,formeig,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&  
  &  nband,nband_disk,npw,nspinor,occ,option,optkg,wff)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: nband
  integer,intent(in) :: nband_disk
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: option
  integer,intent(in) :: optkg
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in),target :: cg(2,mcg)
  real(dp),intent(in),target :: eigen((2*mband)**formeig*mband)
  integer,intent(in),target :: kg_k(3,optkg*npw)
  real(dp),intent(in),target :: occ(mband)
 end subroutine writewf
end interface

interface
 subroutine WffClose(wff,ier)
  use m_wffile
  implicit none
  integer, intent(out) :: ier
  type(wffile_type), intent(inout) :: wff
 end subroutine WffClose
end interface

interface
 subroutine WffDelete(wff,ier)
  use m_wffile
  implicit none
  integer, intent(out) :: ier
  type(wffile_type),intent(inout) :: wff
 end subroutine WffDelete
end interface

interface
 subroutine WffKg(wff,optkg)
  use m_wffile
  implicit none
  integer,intent(in) :: optkg
  type(wffile_type),intent(inout) :: wff
 end subroutine WffKg
end interface

interface
 subroutine WffOffset(wff,sender,spaceComm,ier)
  use m_wffile
  implicit none
  integer          ,intent(out) :: ier
  integer          ,intent(inout) :: sender
  integer          ,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
 end subroutine WffOffset
end interface

interface
 subroutine WffOpen(accesswff,spaceComm,filename,ier,wff,master,me,unwff,&  
  &  spaceComm_mpiio) ! optional argument
  use defs_basis
  use m_wffile
  implicit none
  integer, intent(in) :: accesswff
  integer, intent(out) :: ier
  integer, intent(in) :: master
  integer, intent(in) :: me
  integer, intent(in) :: spaceComm
  integer, intent(in),optional :: spaceComm_mpiio
  integer, intent(in) :: unwff
  character(len=fnlen), intent(in) :: filename
  type(wffile_type), intent(out) :: wff
 end subroutine WffOpen
end interface


!Generic interface of the routines wffreaddatarec
interface wffreaddatarec
 subroutine WffReadDataRec_dp1d(dparray,ierr,ndp,wff)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: ndp
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out) :: dparray(ndp)
 end subroutine WffReadDataRec_dp1d
 subroutine WffReadDataRec_dp2d(dparray,ierr,n1,n2,wff)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out) :: dparray(n1,n2)
 end subroutine WffReadDataRec_dp2d
end interface
!End of the generic interface of wffreaddatarec

interface
 subroutine WffReadNpwRec(ierr,ikpt,isppol,nband_disk,npw,nspinor,wff)
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(out) :: nband_disk
  integer,intent(out) :: npw
  integer,intent(out) :: nspinor
  type(wffile_type),intent(inout) :: wff
 end subroutine WffReadNpwRec
end interface

interface
 subroutine WffReadSkipRec(ierr,nrec,wff)
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: nrec
  type(wffile_type),intent(inout) :: wff
 end subroutine WffReadSkipRec
end interface

interface
 subroutine WffReadWrite_mpio(wff,rdwr,cg,mcg,icg,nband_disk,npwso,npwsotot,depl_mpi_to_seq,ierr)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(in) :: icg
  integer,intent(out) :: ierr
  integer,intent(in) :: mcg
  integer,intent(in) :: nband_disk
  integer,intent(in) :: npwso
  integer,intent(in) :: npwsotot
  integer,intent(in) :: rdwr
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(inout) :: cg(2,mcg)
  integer,intent(in) :: depl_mpi_to_seq(npwso)
 end subroutine WffReadWrite_mpio
end interface


!Generic interface of the routines wffwritedatarec
interface wffwritedatarec
 subroutine WffWriteDataRec_int2d(intarray,ierr,n1,n2,wff)
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: intarray(n1,n2)
 end subroutine WffWriteDataRec_int2d
 subroutine WffWriteDataRec_dp1d(dparray,ierr,ndp,wff)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: ndp
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in) :: dparray(ndp)
 end subroutine WffWriteDataRec_dp1d
 subroutine WffWriteDataRec_dp2d(dparray,ierr,n1,n2,wff)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in) :: dparray(n1,n2)
 end subroutine WffWriteDataRec_dp2d
end interface
!End of the generic interface of wffwritedatarec

interface
 subroutine WffWriteNpwRec(ierr,nband_disk,npw,nspinor,wff,&  
  &  opt_paral) ! optional argument
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: nband_disk
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in),optional :: opt_paral
  type(wffile_type),intent(inout) :: wff
 end subroutine WffWriteNpwRec
end interface

end module interfaces_59_io_mpi
!!***
