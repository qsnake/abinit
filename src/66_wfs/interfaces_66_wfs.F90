!!****m* ABINIT/interfaces_66_wfs
!! NAME
!! interfaces_66_wfs
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/66_wfs
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

module interfaces_66_wfs

 implicit none

interface
 subroutine bestwfs(gcc_block,ghc_block,gscc_block,gscc_calc,&  
  &  gvnlc_block,gvnlc_calc,istwf_k,mpi_enreg,nbdblock,npw_k,nspinor,nvectin,nvectout,wfoptalg)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: gscc_calc
  integer,intent(in) :: gvnlc_calc
  integer,intent(in) :: istwf_k
  integer,intent(in) :: nbdblock
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: nvectin
  integer,intent(in) :: nvectout
  integer,intent(in) :: wfoptalg
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: gcc_block(2,npw_k*nspinor,nbdblock)
  real(dp),intent(inout) :: ghc_block(2,npw_k*nspinor,nbdblock)
  real(dp),intent(inout) :: gscc_block(2,npw_k*nspinor,nbdblock*gscc_calc)
  real(dp),intent(inout) :: gvnlc_block(2,npw_k*nspinor,nbdblock)
 end subroutine bestwfs
end interface

interface
 subroutine envlop(cg,ecut,gmet,icgmod,kg,kpoint,mcg,nband,npw,nspinor)
  use defs_basis
  implicit none
  integer,intent(in) :: icgmod
  integer,intent(in) :: mcg
  integer,intent(in) :: nband
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  real(dp),intent(in) :: ecut
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg(3,npw)
  real(dp),intent(in) :: kpoint(3)
 end subroutine envlop
end interface

interface
 subroutine fxphas(cg,gsc,icg,igsc,istwfk,mcg,mgsc,mpi_enreg,nband_k,npw_k,useoverlap)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: istwfk
  integer,intent(in) :: mcg
  integer,intent(in) :: mgsc
  integer,intent(in) :: nband_k
  integer,intent(in) :: npw_k
  integer,intent(in) :: useoverlap
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(inout) :: gsc(2,mgsc*useoverlap)
 end subroutine fxphas
end interface

interface
 subroutine getdc1(cgq,cprjq,dcwavef,dcwaveprj,ibgq,icgq,istwfk,mcgq,mcprjq,&  
  &  mpi_enreg,natom,nband,npw1,nspinor,optcprj,ortalg,s1cwave0)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: ibgq
  integer,intent(in) :: icgq
  integer,intent(in) :: istwfk
  integer,intent(in) :: mcgq
  integer,intent(in) :: mcprjq
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: npw1
  integer,intent(in) :: nspinor
  integer,intent(in) :: optcprj
  integer,intent(in) :: ortalg
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cgq(2,mcgq)
  type(cprj_type),intent(in) :: cprjq(natom,mcprjq)
  real(dp),intent(out) :: dcwavef(2,npw1*nspinor)
  type(cprj_type),intent(out) :: dcwaveprj(natom,nspinor*optcprj)
  real(dp),intent(in) :: s1cwave0(2,npw1*nspinor)
 end subroutine getdc1
end interface

interface
 subroutine getghc(cpopt,cwavef,cwaveprj,dimffnl,ffnl,filstat,ghc,gsc,gs_ham,&  
  &  gvnlc,kg_k,kinpw,lambda,lmnmax,&  
  &  matblk,mgfft,mpi_enreg,mpsang,mpssoang,&  
  &  natom,ndat,npw,nspinor,ntypat,nvloc,n4,n5,n6,&  
  &  paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,type_calc,vlocal,&  
  &  vxctaulocal) ! optional argument
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: cpopt
  integer,intent(in) :: dimffnl
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: natom
  integer,intent(in) :: ndat
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  integer,intent(in) :: nvloc
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: prtvol
  integer,intent(in) :: sij_opt
  integer,intent(in) :: tim_getghc
  integer,intent(in) :: type_calc
  character(len=fnlen),intent(in) :: filstat
  type(gs_hamiltonian_type),intent(in) :: gs_ham
  real(dp) :: lambda
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cwavef(2,npw*nspinor*ndat)
  type(cprj_type),intent(inout) :: cwaveprj(natom,nspinor*((cpopt+5)/5)*gs_ham%usepaw)
  real(dp),intent(in) :: ffnl(npw,dimffnl,lmnmax,ntypat)
  real(dp),intent(inout) :: ghc(2,npw*nspinor*ndat)
  real(dp),intent(out) :: gsc(2,npw*nspinor*ndat*((sij_opt+1)/2))
  real(dp),intent(inout) :: gvnlc(2,npw*nspinor*ndat)
  integer,intent(in) :: kg_k(3,npw)
  real(dp),intent(in) :: kinpw(npw)
  real(dp),intent(inout) :: ph3d(2,npw,matblk)
  real(dp),intent(inout) :: vlocal(n4,n5,n6,nvloc)
  real(dp),intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)
 end subroutine getghc
end interface

interface
 subroutine getgsc(cg,cprj,dimcprj,dimffnl,ffnl,gs_ham,gsc,ibg,icg,igsc,ikpt,isppol,&  
  &  kg_k,lmnmax,matblk,mcg,mcprj,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,&  
  &  natom,nband,nkpt,npw_k,nspinor,ntypat,ph3d)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: dimffnl
  integer,intent(in) :: ibg
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mgfft
  integer,intent(in) :: mgsc
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: natom
  integer,intent(in) :: nband
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: ntypat
  type(gs_hamiltonian_type),intent(in) :: gs_ham
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,mcg)
  type(cprj_type),intent(in) :: cprj(natom,mcprj)
  integer,intent(in) :: dimcprj(natom)
  real(dp),intent(in) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp),intent(out) :: gsc(2,mgsc)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(inout) :: ph3d(2,npw_k,matblk)
 end subroutine getgsc
end interface

interface
 subroutine listkk(dksqmax,gmet,indkk,kptns1,kptns2,nkpt1,nkpt2,nsym,&  
  &  sppoldbl,symafm,symmat,timrev,use_symrec)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt1
  integer,intent(in) :: nkpt2
  integer,intent(in) :: nsym
  integer,intent(in) :: sppoldbl
  integer,intent(in) :: timrev
  real(dp),intent(out) :: dksqmax
  logical,optional,intent(in) :: use_symrec
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(out) :: indkk(nkpt2*sppoldbl,6)
  real(dp),intent(in) :: kptns1(3,nkpt1)
  real(dp),intent(in) :: kptns2(3,nkpt2)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symmat(3,3,nsym)
 end subroutine listkk
end interface

interface
 subroutine precon(cg,eval,istwf_k,kinpw,mpi_enreg,npw,nspinor,optekin,pcon,vect)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: optekin
  real(dp),intent(in) :: eval
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,npw*nspinor)
  real(dp),intent(in) :: kinpw(npw)
  real(dp),intent(inout) :: pcon(npw)
  real(dp),intent(inout) :: vect(2,npw*nspinor)
 end subroutine precon
end interface

interface
 subroutine precon2(cg,eval,blocksize,iterationnumber,kinpw,&  
  &  mpi_enreg,npw,nspinor,optekin,optpcon,pcon,ghc,vect,vectsize)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: blocksize
  integer,intent(in) :: iterationnumber
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: optekin
  integer,intent(in) :: optpcon
  integer,intent(in) :: vectsize
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(vectsize,blocksize)
  real(dp),intent(in) :: eval(blocksize,blocksize)
  real(dp),intent(in) :: ghc(vectsize,blocksize)
  real(dp),intent(in) :: kinpw(npw)
  real(dp),intent(inout) :: pcon(npw,blocksize)
  real(dp),intent(inout) :: vect(vectsize,blocksize)
 end subroutine precon2
end interface

interface
 subroutine prep_bandfft_tabs(dimffnl,ffnl,gbound,ikpt,kinpw,kpoint,lmnmax,&  
  &  matblk,mgfft,mkmem,mpi_enreg,nkpg,npw_k,ntypat,option,ph3d)
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: dimffnl
  integer, intent(in) :: ikpt
  integer, intent(in) :: lmnmax
  integer, intent(in) :: matblk
  integer, intent(in) :: mgfft
  integer, intent(in) :: mkmem
  integer, intent(in) :: nkpg
  integer, intent(in) :: npw_k
  integer, intent(in) :: ntypat
  integer, intent(in) :: option
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp), intent(in) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  integer, intent(inout) :: gbound(2*mgfft+8,2*option)
  real(dp), intent(in) :: kinpw(npw_k*option)
  real(dp), intent(in) :: kpoint(3)
  real(dp), intent(in) :: ph3d(2,npw_k,matblk)
 end subroutine prep_bandfft_tabs
end interface

interface
 subroutine prep_fourwf(rhoaug,blocksize,cwavef,wfraug,gs_hamk,iblock,ikpt,istwf_k,&  
  &  mgfft,mpi_enreg,nband_k,npw_k,n4,n5,n6,occ_k,paral_kgb,wtk)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer :: blocksize
  integer :: iblock
  integer :: ikpt
  integer :: istwf_k
  integer :: mgfft
  integer :: n4
  integer :: n5
  integer :: n6
  integer :: nband_k
  integer :: npw_k
  integer :: paral_kgb
  type(gs_hamiltonian_type) :: gs_hamk
  type(mpi_type) :: mpi_enreg
  real(dp) :: wtk
  real(dp) :: cwavef(2,npw_k*blocksize)
  real(dp) :: occ_k(nband_k)
  real(dp) :: rhoaug(n4,n5,n6)
  real(dp) :: wfraug(2,n4,n5,n6)
 end subroutine prep_fourwf
end interface

interface
 subroutine prep_getghc(cwavef,dimffnl,dtfil,gs_hamk,gvnlc,gwavef,swavef,ikpt,istwf_k,&  
  &  lambda,lmnmax,matblk,blocksize,mgfft,mpi_enreg,mpsang,mpssoang,natom,npw_k,&  
  &  nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,prtvol,sij_opt,vlocal,&  
  &  vxctaulocal) ! optional argument
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer :: blocksize
  integer :: dimffnl
  integer :: ikpt
  integer :: istwf_k
  integer :: lmnmax
  integer :: matblk
  integer :: mgfft
  integer :: mpsang
  integer :: mpssoang
  integer :: n4
  integer :: n5
  integer :: n6
  integer :: natom
  integer :: npw_k
  integer :: nspinor
  integer :: ntypat
  integer :: nvloc
  integer :: paral_kgb
  integer :: prtvol
  integer :: sij_opt
  type(datafiles_type) :: dtfil
  type(gs_hamiltonian_type) :: gs_hamk
  real(dp) :: lambda
  type(mpi_type) :: mpi_enreg
  real(dp) :: cwavef(2,npw_k*nspinor*blocksize)
  real(dp) :: gvnlc(2,npw_k*nspinor*blocksize)
  real(dp) :: gwavef(2,npw_k*nspinor*blocksize)
  real(dp) :: swavef(2,npw_k*nspinor*blocksize)
  real(dp) :: vlocal(n4,n5,n6,nvloc)
  real(dp), intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)
 end subroutine prep_getghc
end interface

interface
 subroutine prep_index_wavef_bandpp(nproc_band,bandpp,&  
  nspinor,ndatarecv,&  
  recvcounts,rdispls,&  
  index_wavef_band)
  implicit none
  integer,intent(in) :: bandpp
  integer,intent(in) :: ndatarecv
  integer,intent(in) :: nproc_band
  integer,intent(in) :: nspinor
  integer,pointer :: index_wavef_band(:)
  integer,intent(in) :: rdispls(nproc_band)
  integer,intent(in) :: recvcounts(nproc_band)
 end subroutine prep_index_wavef_bandpp
end interface

interface
 subroutine prep_kpgio(accesswff,ecut_eff,exchn2n3d,gmet,istwfk,kg,kptns,fnametmp_kg,mgfft,mkmem,mode_paral,&  
  &  mpi_enreg,mpw,nband,nkpt,npwarr,npwtot,nsppol,unkg)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: accesswff
  integer,intent(in) :: exchn2n3d
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: unkg
  real(dp),intent(in) :: ecut_eff
  character(len=fnlen),intent(in) :: fnametmp_kg
  character(len=4),intent(in) :: mode_paral
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(out) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(out) :: npwarr(nkpt)
  integer,intent(out) :: npwtot(nkpt)
 end subroutine prep_kpgio
end interface

interface
 subroutine prep_nonlop(atindx1,choice,cpopt,cwaveprj,dimenl1,dimenl2,dimffnl,enl,enlout_block,&  
  &  gmet,gprimd,idir,ikpt,indlmn,istwf_k,&  
  &  kpt,lambdablock,lmnmax,matblk,&  
  &  blocksize,mgfft,mpi_enreg,mpsang,mpssoang,&  
  &  natom,nattyp,ngfft,nkpg,nloalg,nnlout,npw_k,&  
  &  nspinor,nspinortot,ntypat,paw_opt,phkxred,ph1d,signs,sij,gsc,&  
  &  tim_nonlop,ucvol,useylm,cwavef,gvnlc,use_gpu_cuda)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: blocksize
  integer,intent(in) :: choice
  integer,intent(in) :: cpopt
  integer,intent(in) :: dimenl1
  integer,intent(in) :: dimenl2
  integer,intent(in) :: dimffnl
  integer,intent(in) :: idir
  integer,intent(in) :: ikpt
  integer,intent(in) :: istwf_k
  integer,intent(in) :: lmnmax
  integer,intent(in) :: matblk
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpssoang
  integer,intent(in) :: natom
  integer,intent(in) :: nkpg
  integer,intent(in) :: nnlout
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: nspinortot
  integer,intent(in) :: ntypat
  integer,intent(in) :: paw_opt
  integer,intent(in) :: signs
  integer :: tim_nonlop
  integer,intent(in),optional :: use_gpu_cuda
  integer,intent(in) :: useylm
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: nloalg(5)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(inout) :: cwavef(2,npw_k*nspinortot*blocksize)
  type(cprj_type),intent(inout),target :: cwaveprj(natom,nspinortot*mpi_enreg%bandpp*((cpopt+5)/5))
  real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2)
  real(dp),intent(out) :: enlout_block(nnlout*blocksize)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: gsc(2,npw_k*nspinor*blocksize*(paw_opt/3))
  real(dp),intent(out) :: gvnlc(2,npw_k*nspinor*blocksize)
  integer,intent(in) :: indlmn(6,lmnmax,ntypat)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(in) :: lambdablock(blocksize)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp),intent(in) :: phkxred(2,natom)
  real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
 end subroutine prep_nonlop
end interface

interface
 subroutine prep_sort_wavef_spin(nproc_band,nspinor,ndatarecv,recvcounts,rdispls,index_wavef)
  implicit none
  integer,intent(in) :: ndatarecv
  integer,intent(in) :: nproc_band
  integer,intent(in) :: nspinor
  integer,pointer :: index_wavef(:)
  integer,intent(in) :: rdispls(nproc_band)
  integer,intent(in) :: recvcounts(nproc_band)
 end subroutine prep_sort_wavef_spin
end interface

interface
 subroutine prep_wavef_sym_do(mpi_enreg,bandpp,nspinor,&  
  ndatarecv,&  
  ndatarecv_tot,ndatasend_sym,tab_proc,&  
  cwavef_alltoall,&  
  sendcounts_sym,sdispls_sym,&  
  recvcounts_sym,rdispls_sym,&  
  ewavef_alltoall_sym,&  
  index_wavef_send)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bandpp
  integer,intent(in) :: ndatarecv
  integer,intent(in) :: ndatarecv_tot
  integer,intent(in) :: ndatasend_sym
  integer,intent(in) :: nspinor
  type(mpi_type),intent(in) :: mpi_enreg
  integer,pointer :: index_wavef_send(:)
  integer,pointer :: rdispls_sym(:)
  integer,pointer :: recvcounts_sym(:)
  integer,pointer :: sdispls_sym(:)
  integer,pointer :: sendcounts_sym(:)
  integer,pointer :: tab_proc(:)
  real(dp),intent(inout) :: cwavef_alltoall(2,ndatarecv*nspinor*bandpp)
  real(dp),pointer :: ewavef_alltoall_sym(:,:)
 end subroutine prep_wavef_sym_do
end interface

interface
 subroutine prep_wavef_sym_undo(mpi_enreg,bandpp,nspinor,&  
  ndatarecv,&  
  ndatarecv_tot,ndatasend_sym,idatarecv0,&  
  gwavef_alltoall,&  
  sendcounts_sym,sdispls_sym,&  
  recvcounts_sym,rdispls_sym,&  
  gwavef_alltoall_sym,&  
  index_wavef_send)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bandpp
  integer,intent(in) :: idatarecv0
  integer,intent(in) :: ndatarecv
  integer,intent(in) :: ndatarecv_tot
  integer,intent(in) :: ndatasend_sym
  integer,intent(in) :: nspinor
  type(mpi_type),intent(in) :: mpi_enreg
  integer,pointer :: index_wavef_send(:)
  integer,pointer :: rdispls_sym(:)
  integer,pointer :: recvcounts_sym(:)
  integer,pointer :: sdispls_sym(:)
  integer,pointer :: sendcounts_sym(:)
  real(dp),intent(inout) :: gwavef_alltoall(2,ndatarecv*nspinor*bandpp)
  real(dp),pointer :: gwavef_alltoall_sym(:,:)
 end subroutine prep_wavef_sym_undo
end interface

interface
 subroutine projbd(cg,direc,iband0,icg,iscg,istwf_k,mcg,mpi_enreg,mscg,nband,&  
  &  npw,nspinor,ortalg,printopt,scg,scprod,scprod_io,tim_projbd,useoverlap)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iband0
  integer,intent(in) :: icg
  integer,intent(in) :: iscg
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mcg
  integer,intent(in) :: mscg
  integer,intent(in) :: nband
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: ortalg
  integer,intent(in) :: printopt
  integer,intent(in) :: scprod_io
  integer,intent(in) :: tim_projbd
  integer,intent(in) :: useoverlap
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(inout) :: direc(2,npw*nspinor)
  real(dp),intent(in) :: scg(2,mscg*useoverlap)
  real(dp),intent(inout) :: scprod(2,nband)
 end subroutine projbd
end interface

interface
 subroutine pw_orthon(icg,igsc,istwf_k,mcg,mgsc,mpi_enreg,nelem,nvec,&  
  &  ortalgo,ovl_vecnm,useoverlap,vecnm)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mcg
  integer,intent(in) :: mgsc
  integer,intent(in) :: nelem
  integer,intent(in) :: nvec
  integer,intent(in) :: ortalgo
  integer,intent(in) :: useoverlap
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: ovl_vecnm(2,mgsc*useoverlap)
  real(dp),intent(inout) :: vecnm(2,mcg)
 end subroutine pw_orthon
end interface

interface
 subroutine sdirot(cg,evec,icg,mcg,ndim,num,npw)
  use defs_basis
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: mcg
  integer,intent(in) :: ndim
  integer,intent(in) :: npw
  integer,intent(in) :: num
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(in) :: evec(2*ndim,num)
 end subroutine sdirot
end interface

interface
 subroutine wfconv(ceksp2,cg1,cg2,debug,ecut1,ecut2,ecut2_eff,&  
  &  eig_k1,eig_k2,exchn2n3d,formeig,gmet1,gmet2,icg1,icg2,&  
  &  ikpt1,ikpt10,ikpt2,indkk,inplace,isppol2,istwfk1,istwfk2,&  
  &  kg1,kg2,kptns1,kptns2,mband1,mband2,mcg1,mcg2,mgfft,mpi_enreg1,mpi_enreg2,mpw1,mpw2,nbd1,nbd2,&  
  &  ngfft,nkpt1,nkpt2,npw1,npw2,nspinor1,nspinor2,nsym,&  
  &  occ_k1,occ_k2,optorth,restart,rprimd2,sppoldbl,symrel,tnons)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: ceksp2
  integer,intent(in) :: debug
  integer,intent(in) :: exchn2n3d
  integer,intent(in) :: formeig
  integer,intent(in) :: icg1
  integer,intent(in) :: icg2
  integer,intent(in) :: ikpt1
  integer,intent(inout) :: ikpt10
  integer,intent(in) :: ikpt2
  integer,intent(in) :: inplace
  integer,intent(in) :: isppol2
  integer,intent(in) :: mband1
  integer,intent(in) :: mband2
  integer,intent(in) :: mcg1
  integer,intent(in) :: mcg2
  integer,intent(in) :: mgfft
  integer,intent(in) :: mpw1
  integer,intent(in) :: mpw2
  integer,intent(in) :: nbd1
  integer,intent(in) :: nbd2
  integer,intent(in) :: nkpt1
  integer,intent(in) :: nkpt2
  integer,intent(inout) :: npw1
  integer,intent(inout) :: npw2
  integer,intent(in) :: nspinor1
  integer,intent(in) :: nspinor2
  integer,intent(in) :: nsym
  integer,intent(in) :: optorth
  integer,intent(in) :: restart
  integer,intent(in) :: sppoldbl
  real(dp),intent(in) :: ecut1
  real(dp),intent(in) :: ecut2
  real(dp),intent(in) :: ecut2_eff
  type(mpi_type),intent(inout) :: mpi_enreg1
  type(mpi_type),intent(inout) :: mpi_enreg2
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: cg1(2,mcg1)
  real(dp),intent(inout) :: cg2(2,mcg2)
  real(dp),intent(inout) :: eig_k1(mband1*(2*mband1)**formeig)
  real(dp),intent(inout) :: eig_k2(mband2*(2*mband2)**formeig)
  real(dp),intent(in) :: gmet1(3,3)
  real(dp),intent(in) :: gmet2(3,3)
  integer,intent(in) :: indkk(nkpt2*sppoldbl,6)
  integer,intent(in) :: istwfk1(nkpt1)
  integer,intent(in) :: istwfk2(nkpt2)
  integer,intent(inout) :: kg1(3,mpw1)
  integer,intent(inout) :: kg2(3,mpw2)
  real(dp),intent(in) :: kptns1(3,nkpt1)
  real(dp),intent(in) :: kptns2(3,nkpt2)
  real(dp),intent(inout) :: occ_k1(mband1)
  real(dp),intent(inout) :: occ_k2(mband2)
  real(dp),intent(in) :: rprimd2(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine wfconv
end interface

interface
 subroutine zprecon3(cg,eval,blocksize,iterationnumber,kinpw,mpi_enreg,npw,nspinor,optekin,optpcon,pcon,ghc,vect,vectsize)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: blocksize
  integer,intent(in) :: iterationnumber
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: optekin
  integer,intent(in) :: optpcon
  integer,intent(in) :: vectsize
  type(mpi_type),intent(inout) :: mpi_enreg
  complex(dpc),intent(in) :: cg(vectsize,blocksize)
  complex(dpc),intent(in) :: eval(blocksize,blocksize)
  complex(dpc),intent(in) :: ghc(vectsize,blocksize)
  real(dp) :: kinpw(npw)
  real(dp),intent(inout) :: pcon(npw,blocksize)
  complex(dpc),intent(inout) :: vect(vectsize,blocksize)
 end subroutine zprecon3
end interface

end module interfaces_66_wfs
!!***
