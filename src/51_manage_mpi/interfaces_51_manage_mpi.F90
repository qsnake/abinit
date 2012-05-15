!!****m* ABINIT/interfaces_51_manage_mpi
!! NAME
!! interfaces_51_manage_mpi
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/51_manage_mpi
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

module interfaces_51_manage_mpi

 implicit none

interface
 subroutine clnmpi_atom(mpi_enreg)
  use defs_abitypes
  implicit none
  type(mpi_type), intent(inout) :: mpi_enreg
 end subroutine clnmpi_atom
end interface

interface
 subroutine clnmpi_band(mpi_enreg)
  use defs_abitypes
  implicit none
  type(mpi_type), intent(inout) :: mpi_enreg
 end subroutine clnmpi_band
end interface

interface
 subroutine clnmpi_bandfft(mpi_enreg)
  use defs_abitypes
  implicit none
  type(mpi_type), intent(inout) :: mpi_enreg
 end subroutine clnmpi_bandfft
end interface

interface
 subroutine clnmpi_fft(mpi_enreg)
  use defs_abitypes
  implicit none
  type(mpi_type), intent(inout) :: mpi_enreg
 end subroutine clnmpi_fft
end interface

interface
 subroutine clnmpi_grid(mpi_enreg)
  use defs_abitypes
  implicit none
  type(mpi_type), intent(inout) :: mpi_enreg
 end subroutine clnmpi_grid
end interface

interface
 subroutine clnmpi_gs(mpi_enreg)
  use defs_abitypes
  implicit none
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine clnmpi_gs
end interface

interface
 subroutine clnmpi_img(mpi_enreg)
  use defs_abitypes
  implicit none
  type(mpi_type), intent(inout) :: mpi_enreg
 end subroutine clnmpi_img
end interface

interface
 subroutine clnmpi_respfn(mpi_enreg,spaceComm_noparalrf)
  use defs_abitypes
  implicit none
  integer,intent(in) :: spaceComm_noparalrf
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine clnmpi_respfn
end interface

interface
 subroutine cprj_diskinit_r(atind,dimcp,iorder,mkmem,natom,ncpgr,nlmn,nspinor,uncp)
  implicit none
  integer,intent(in) :: dimcp
  integer,intent(in) :: iorder
  integer,intent(in) :: mkmem
  integer,intent(in) :: natom
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nspinor
  integer,intent(in) :: uncp
  integer,intent(in) :: atind(natom)
  integer,intent(in) :: nlmn(dimcp)
 end subroutine cprj_diskinit_r
end interface

interface
 subroutine cprj_diskinit_w(atind,dimcp,iorder,mkmem,natom,ncpgr,nlmn,nspinor,uncp)
  implicit none
  integer,intent(in) :: dimcp
  integer,intent(in) :: iorder
  integer,intent(in) :: mkmem
  integer,intent(in) :: natom
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nspinor
  integer,intent(in) :: uncp
  integer,intent(in) :: atind(natom)
  integer,intent(in) :: nlmn(dimcp)
 end subroutine cprj_diskinit_w
end interface

interface
 subroutine cprj_diskskip(mkmem,ncpgr,nrec,uncp)
  implicit none
  integer,intent(in) :: mkmem
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nrec
  integer,intent(in) :: uncp
 end subroutine cprj_diskskip
end interface

interface
 subroutine cprj_get(atind,cprj_k,cprj,dimcp,iband1,ibg,ikpt,iorder,isppol,mband,&  
  &  mkmem,mpi_enreg,natom,nband,nband_k,nspinor,nsppol,uncp,&  
  &  icpgr,ncpgr) ! optionals arguments
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: dimcp
  integer,intent(in) :: iband1
  integer,intent(in) :: ibg
  integer,intent(in),optional :: icpgr
  integer,intent(in) :: ikpt
  integer,intent(in) :: iorder
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: nband_k
  integer,intent(in),optional :: ncpgr
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: uncp
  type(mpi_type), intent(inout) :: mpi_enreg
  integer,intent(in) :: atind(natom)
  type(cprj_type),intent(in) :: cprj(dimcp,nspinor*mband*mkmem*nsppol)
  type(cprj_type),intent(inout) :: cprj_k(dimcp,nspinor*nband)
 end subroutine cprj_get
end interface

interface
 subroutine cprj_put(atind,cprj_k,cprj,dimcp,iband1,ibg,ikpt,iorder,isppol,mband,&  
  &  mkmem,mpi_enreg,natom,nband,nband_k,nlmn,nspinor,nsppol,spaceComm_band,uncp,&  
  &  to_be_gathered) ! Optional argument
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: dimcp
  integer,intent(in) :: iband1
  integer,intent(in) :: ibg
  integer,intent(in) :: ikpt
  integer,intent(in) :: iorder
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: nband_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: spaceComm_band
  integer,intent(in) :: uncp
  type(mpi_type), intent(inout) :: mpi_enreg
  logical,optional,intent(in) :: to_be_gathered
  integer :: atind(natom)
  type(cprj_type),intent(out) :: cprj(dimcp,nspinor*mband*mkmem*nsppol)
  type(cprj_type),intent(in) :: cprj_k(dimcp,nspinor*nband)
  integer :: nlmn(dimcp)
 end subroutine cprj_put
end interface

interface
 subroutine cprj_exch(natom,n2dim,nlmn,ncpgr,Cprj_send,Cprj_recv,sender,receiver,spaceComm,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n2dim
  integer,intent(in) :: natom
  integer,intent(in) :: ncpgr
  integer,intent(in) :: receiver
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  type(cprj_type),intent(inout) :: Cprj_recv(:,:)
  type(cprj_type),intent(in) :: Cprj_send(:,:)
  integer,intent(in) :: nlmn(natom)
 end subroutine cprj_exch
end interface

interface
 subroutine cprj_mpi_send(natom,n2dim,nlmn,ncpgr,cprj_out,receiver,spaceComm,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n2dim
  integer,intent(in) :: natom
  integer,intent(in) :: ncpgr
  integer,intent(in) :: receiver
  integer,intent(in) :: spaceComm
  type(cprj_type),intent(in) :: cprj_out(:,:)
  integer,intent(in) :: nlmn(natom)
 end subroutine cprj_mpi_send
end interface

interface
 subroutine cprj_mpi_recv(natom,n2dim,nlmn,ncpgr,cprj_in,sender,spaceComm,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n2dim
  integer,intent(in) :: natom
  integer,intent(in) :: ncpgr
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  type(cprj_type),intent(inout) :: cprj_in(:,:)
  integer,intent(in) :: nlmn(natom)
 end subroutine cprj_mpi_recv
end interface

interface
 subroutine cprj_mpi_allgather(cprj_loc,cprj_gat,natom,n2dim,nlmn,ncpgr,nproc,spaceComm,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n2dim
  integer,intent(in) :: natom
  integer,intent(in) :: ncpgr
  integer,intent(in) :: nproc
  integer,intent(in) :: spaceComm
  type(cprj_type),intent(out) :: cprj_gat(:,:)
  type(cprj_type),intent(in) :: cprj_loc(:,:)
  integer,intent(in) :: nlmn(natom)
 end subroutine cprj_mpi_allgather
end interface

interface
 subroutine cprj_bcast(natom,n2dim,nlmn,ncpgr,Cprj,master,spaceComm,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: master
  integer,intent(in) :: n2dim
  integer,intent(in) :: natom
  integer,intent(in) :: ncpgr
  integer,intent(in) :: spaceComm
  type(cprj_type),intent(inout) :: Cprj(natom,n2dim)
  integer,intent(in) :: nlmn(natom)
 end subroutine cprj_bcast
end interface

interface
 subroutine cprj_transpose(cprjin,cprjout,cprj_bandpp,natom,nband,nspinor,spaceComm)
  use defs_datatypes
  implicit none
  integer :: cprj_bandpp
  integer :: natom
  integer :: nband
  integer :: nspinor
  integer :: spaceComm
  type(cprj_type),intent(in) :: cprjin(:,:)
  type(cprj_type),intent(out) :: cprjout(:,:)
 end subroutine cprj_transpose
end interface

interface
 subroutine cprj_transpose_all(cprjin,cprjout,dtfil,mband,mkmem,mpi_enreg,natom,nband,&  
  &  nkpt,nspinor,nsppol,paral_kgb,spaceComm)
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: spaceComm
  type(datafiles_type),intent(in) :: dtfil
  type(mpi_type),intent(inout) :: mpi_enreg
  type(cprj_type),intent(in) :: cprjin(:,:)
  type(cprj_type),intent(out) :: cprjout(:,:)
  integer :: nband(nkpt*nsppol)
 end subroutine cprj_transpose_all
end interface

interface
 subroutine cprj_gather_spin(cprj,cprj_gat,natom,n2size,nspinor,nspinortot,&  
  &  spaceComm_spin,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n2size
  integer,intent(in) :: natom
  integer,intent(in) :: nspinor
  integer,intent(in) :: nspinortot
  integer,intent(in) :: spaceComm_spin
  type(cprj_type),intent(in) :: cprj(:,:)
  type(cprj_type),intent(out) :: cprj_gat(:,:)
 end subroutine cprj_gather_spin
end interface

interface
 subroutine distrb2(mband, nband, nkpt, nsppol, mpi_enreg)
  use defs_abitypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine distrb2
end interface

interface
 subroutine herald(code_name,code_version,iout)
  implicit none
  integer,intent(in) :: iout
  character(len=24),intent(in) :: code_name
  character(len=6),intent(in) :: code_version
 end subroutine herald
end interface

interface
 subroutine initmpi_atom(dtset,mpi_enreg)
  use defs_abitypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_atom
end interface

interface
 subroutine initmpi_band(mpi_enreg,nband,nkpt,nsppol)
  use defs_abitypes
  implicit none
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine initmpi_band
end interface

interface
 subroutine initmpi_fft(dtset,mpi_enreg)
  use defs_abitypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_fft
end interface

interface
 subroutine initmpi_grid(mpi_enreg)
  use defs_abitypes
  implicit none
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_grid
end interface

interface
 subroutine initmpi_gs(dtset,mpi_enreg)
  use defs_abitypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_gs
end interface

interface
 subroutine initmpi_img(dtset,mpi_enreg,option)
  use defs_abitypes
  implicit none
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_img
end interface

interface
 subroutine initmpi_respfn(mpi_enreg,spaceComm,spaceComm_noparalrf)
  use defs_abitypes
  implicit none
  integer,intent(in) :: spaceComm
  integer,intent(out) :: spaceComm_noparalrf
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_respfn
end interface

interface
 subroutine initmpi_seq(mpi_enreg)
  use defs_abitypes
  implicit none
  type(mpi_type),intent(out) :: mpi_enreg
 end subroutine initmpi_seq
end interface

interface
 subroutine initmpi_world(mpi_enreg,nproc)
  use defs_abitypes
  implicit none
  integer, intent(in) :: nproc
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine initmpi_world
end interface

interface
 subroutine leave_test()
  implicit none
 end subroutine leave_test
end interface

interface
 subroutine init_mpi_enreg(mpi_enreg,init_mpi)
  use defs_abitypes
  implicit none
  type(mpi_type),intent(inout) :: MPI_enreg
  logical,optional,intent(in) :: init_mpi
 end subroutine init_mpi_enreg
end interface

interface
 subroutine nullify_mpi_enreg(MPI_enreg)
  use defs_abitypes
  implicit none
  type(mpi_type),intent(inout) :: MPI_enreg
 end subroutine nullify_mpi_enreg
end interface

interface
 subroutine destroy_mpi_enreg(MPI_enreg)
  use defs_abitypes
  implicit none
  type(mpi_type),intent(inout) :: MPI_enreg
 end subroutine destroy_mpi_enreg
end interface

interface
 subroutine copy_mpi_enreg(MPI_enreg1,MPI_enreg2,opt_bandfft)
  use defs_abitypes
  implicit none
  integer :: opt_bandfft
  type(mpi_type),intent(inout) :: MPI_enreg2
  type(mpi_type),intent(inout) :: mpi_enreg1
 end subroutine copy_mpi_enreg
end interface

interface
 subroutine my_indeces(MPI_enreg,ikpt,isppol,nkpt,nsppol,nspinor,npwarr,nband,kindex,bdtot_index,ibg,ikg,ierr)
  use defs_abitypes
  implicit none
  integer,intent(out) :: bdtot_index
  integer,intent(out) :: ibg
  integer,intent(out) :: ierr
  integer,intent(out) :: ikg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(out) :: kindex
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(mpi_type),intent(in) :: MPI_enreg
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
 end subroutine my_indeces
end interface

interface
 subroutine pre_gather(array,array_allgather,n1,n2,n3,n4,mpi_enreg)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: n4
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: array(n1,n2,n4,1)
  real(dp),intent(inout) :: array_allgather(n1,n2,n3,1)
 end subroutine pre_gather
end interface

interface
 subroutine pre_scatter(array,array_allgather,n1,n2,n3,n4,mpi_enreg)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: n4
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(out) :: array(n1,n2,n4,1)
  real(dp),intent(in) :: array_allgather(n1,n2,n3,1)
 end subroutine pre_scatter
end interface

interface
 subroutine pspheads_comm(npsp,pspheads,test_paw)
  use defs_datatypes
  implicit none
  integer,intent(in) :: npsp
  integer,intent(inout) :: test_paw
  type(pspheader_type),intent(inout) :: pspheads(npsp)
 end subroutine pspheads_comm
end interface

interface
 subroutine xcomm_world(mpi_enreg,spaceComm,myrank,mysize)
  use defs_abitypes
  implicit none
  integer,optional,intent(out) :: myrank
  integer,optional,intent(out) :: mysize
  integer,intent(out) :: spaceComm
  type(mpi_type),intent(in) :: mpi_enreg
 end subroutine xcomm_world
end interface

interface
 subroutine xcomm_self(spaceComm)
  implicit none
  integer,intent(out) :: spaceComm
 end subroutine xcomm_self
end interface

interface
 subroutine xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft)
  use defs_abitypes
  implicit none
  integer,intent(out) :: spaceComm
  integer,intent(in),optional :: spaceComm_bandfft
  type(mpi_type),intent(in) :: mpi_enreg
 end subroutine xcomm_init
end interface

interface
 subroutine xmaster_init(mpi_enreg,master)
  use defs_abitypes
  implicit none
  integer,intent(out) :: master
  type(mpi_type),intent(in) :: mpi_enreg
 end subroutine xmaster_init
end interface

interface
 subroutine xmaster_init_fft(mpi_enreg,master)
  use defs_abitypes
  implicit none
  integer,intent(out) :: master
  type(mpi_type),intent(in) :: mpi_enreg
 end subroutine xmaster_init_fft
end interface

interface
 subroutine xme_init(mpi_enreg,me,option_comm)
  use defs_abitypes
  implicit none
  integer,intent(out) :: me
  integer,intent(in),optional :: option_comm
  type(mpi_type),intent(in) :: mpi_enreg
 end subroutine xme_init
end interface

interface
 subroutine xproc_init(mpi_enreg,nproc_max)
  use defs_abitypes
  implicit none
  integer,intent(out) :: nproc_max
  type(mpi_type),intent(in) :: mpi_enreg
 end subroutine xproc_init
end interface

interface
 subroutine xdefineOff(formeig,wff,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)
  use defs_abitypes
  use m_wffile
  implicit none
  integer, intent(in) :: formeig
  integer, intent(in) :: nkpt
  integer, intent(in) :: nspinor
  integer, intent(in) :: nsppol
  type(mpi_type),intent(in) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  integer, intent(in) :: nband(nkpt*nsppol)
  integer, intent(in) :: npwarr(nkpt)
 end subroutine xdefineOff
end interface


!Generic interface of the routines xderiveread
interface xderiveread
 subroutine xderiveRead_int(wff,xval,ierr)
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(out) :: xval
  type(wffile_type),intent(inout) :: wff
 end subroutine xderiveRead_int
 subroutine xderiveRead_int1d(wff,xval,n1,spaceComm,ierr)
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  integer,intent(out) :: xval(:)
 end subroutine xderiveRead_int1d
 subroutine xderiveRead_int2d(wff,xval,n1,n2,spaceComm,ierr)
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  integer,intent(out) :: xval(:,:)
 end subroutine xderiveRead_int2d
 subroutine xderiveRead_dp(wff,xval,ierr)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out) :: xval
 end subroutine xderiveRead_dp
 subroutine xderiveRead_dp1d(wff,xval,n1,spaceComm,ierr)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out) :: xval(:)
 end subroutine xderiveRead_dp1d
 subroutine xderiveRead_dp2d(wff,xval,n1,n2,spaceComm,ierr)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out) :: xval(:,:)
 end subroutine xderiveRead_dp2d
 subroutine xderiveRead_int2d_displ(wff,xval,n1,n2,spaceComm,displace,ierr)
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: displace(:)
  integer,intent(out) :: xval(:,:)
 end subroutine xderiveRead_int2d_displ
 subroutine xderiveRead_dp2d_displ(wff,xval,n1,n2,spaceComm,displace,ierr)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: displace(:)
  real(dp),intent(out) :: xval(:,:)
 end subroutine xderiveRead_dp2d_displ
 subroutine xderiveReadVal_char(wff,xval,n,ierr)
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n
  type(wffile_type),intent(inout) :: wff
  character(len=*),intent(out) :: xval
 end subroutine xderiveReadVal_char
 subroutine xmpi_read_int2d(wff,xval,spaceComm,at_option,ierr)
  use m_wffile
  implicit none
  integer,intent(in) :: at_option
  integer,intent(out) :: ierr
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  integer,intent(out) :: xval(:,:)
 end subroutine xmpi_read_int2d
end interface
!End of the generic interface of xderiveread


!Generic interface of the routines xderivewrite
interface xderivewrite
 subroutine xderiveWrite_int(wff,xval,ierr)
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: xval
  type(wffile_type),intent(inout) :: wff
 end subroutine xderiveWrite_int
 subroutine xderiveWrite_int1d(wff,xval,n1,spaceComm,ierr)
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: xval(:)
 end subroutine xderiveWrite_int1d
 subroutine xderiveWrite_int2d(wff,xval,n1,n2,spaceComm,ierr)
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: xval(:,:)
 end subroutine xderiveWrite_int2d
 subroutine xderiveWrite_dp(wff,xval,ierr)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in) :: xval
 end subroutine xderiveWrite_dp
 subroutine xderiveWrite_dp1d(wff,xval,n1,spaceComm,ierr)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in) :: xval(:)
 end subroutine xderiveWrite_dp1d
 subroutine xderiveWrite_dp2d(wff,xval,n1,n2,spaceComm,ierr)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in) :: xval(:,:)
 end subroutine xderiveWrite_dp2d
 subroutine xderiveWrite_int2d_displ(wff,xval,n1,n2,spaceComm,displace,ierr)
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: displace(:)
  integer,intent(in) :: xval(:,:)
 end subroutine xderiveWrite_int2d_displ
 subroutine xderiveWrite_dp2d_displ(wff,xval,n1,n2,spaceComm,displace,ierr)
  use defs_basis
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: displace(:)
  real(dp),intent(in) :: xval(:,:)
 end subroutine xderiveWrite_dp2d_displ
 subroutine xderiveWrite_char(wff,xval,n,ierr)
  use m_wffile
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n
  type(wffile_type),intent(inout) :: wff
  character(len=*),intent(in) :: xval
 end subroutine xderiveWrite_char
end interface
!End of the generic interface of xderivewrite

end module interfaces_51_manage_mpi
!!***
