!!****m* ABINIT/interfaces_79_seqpar_mpi
!! NAME
!! interfaces_79_seqpar_mpi
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/79_seqpar_mpi
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

module interfaces_79_seqpar_mpi

 implicit none

interface
 subroutine inwffil(ask_accurate,cg,dtset,ecut,ecut_eff,eigen,exchn2n3d,&  
  &  formeig,gmet,hdr,ireadwf,istwfk,kg,kptns,localrdwf,mband,&  
  &  mcg,mkmem,mpi_enreg,mpw,nband,ngfft,nkpt,npwarr,&  
  &  nsppol,nsym,occ,optorth,rprimd,symafm,symrel,tnons,unkg,wff1,&  
  &  wffnow,unwff1,unwfnow,wffnm,wft1nm, wvl)
  use defs_basis
  use defs_abitypes
  use m_wffile
  use defs_wvltypes
  implicit none
  integer, intent(in) :: ask_accurate
  integer, intent(in) :: exchn2n3d
  integer, intent(in) :: formeig
  integer, intent(in) :: ireadwf
  integer, intent(in) :: localrdwf
  integer, intent(in) :: mband
  integer, intent(in) :: mcg
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpw
  integer, intent(in) :: nkpt
  integer, intent(in) :: nsppol
  integer, intent(in) :: nsym
  integer, intent(in) :: optorth
  integer, intent(in) :: unkg
  integer, intent(in) :: unwff1
  integer, intent(in) :: unwfnow
  type(dataset_type), intent(in) :: dtset
  real(dp), intent(in) :: ecut
  real(dp), intent(in) :: ecut_eff
  type(hdr_type), intent(inout) :: hdr
  type(mpi_type),intent(inout),target :: mpi_enreg
  type(wffile_type), intent(out) :: wff1
  character(len=fnlen), intent(in) :: wffnm
  type(wffile_type), intent(out) :: wffnow
  character(len=fnlen), intent(in) :: wft1nm
  type(wvl_data), intent(inout) :: wvl
  integer, intent(in) :: ngfft(18)
  real(dp), intent(out),target :: cg(2,mcg)
  real(dp), intent(out),target :: eigen((2*mband)**formeig*mband*nkpt*nsppol)
  real(dp), intent(in) :: gmet(3,3)
  integer, intent(in) :: istwfk(nkpt)
  integer, intent(in) :: kg(3,mpw*mkmem)
  real(dp), intent(in) :: kptns(3,nkpt)
  integer, intent(in),target :: nband(nkpt*nsppol)
  integer, intent(in) :: npwarr(nkpt)
  real(dp), intent(inout) :: occ(mband*nkpt*nsppol)
  real(dp), intent(in) :: rprimd(3,3)
  integer, intent(in) :: symafm(nsym)
  integer, intent(in) :: symrel(3,3,nsym)
  real(dp), intent(in) :: tnons(3,nsym)
 end subroutine inwffil
end interface

interface
 subroutine lobpcgIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&  
  &  kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&  
  &  nband_k,nbdblock,npw_k,ntypat,nvloc,n4,n5,n6,ph3d,prtvol,&  
  &  resid_k,subham,subovl,subvnl,use_subovl,vlocal,&  
  &  vxctaulocal) ! optional argument
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer :: dimffnl
  integer :: icg
  integer :: igsc
  integer :: lmnmax
  integer :: matblk
  integer :: mcg
  integer :: mgfft
  integer :: mgsc
  integer :: mpsang
  integer :: mpssoang
  integer :: n4
  integer :: n5
  integer :: n6
  integer :: natom
  integer :: nband_k
  integer :: nbdblock
  integer :: npw_k
  integer :: ntypat
  integer :: nvloc
  integer :: prtvol
  integer :: use_subovl
  type(datafiles_type) :: dtfil
  type(dataset_type) :: dtset
  type(gs_hamiltonian_type) :: gs_hamk
  type(mpi_type) :: mpi_enreg
  real(dp) :: cg(2,mcg)
  real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp) :: gsc(2,mgsc)
  integer :: kg_k(3,npw_k)
  real(dp) :: kinpw(npw_k)
  real(dp) :: ph3d(2,npw_k,matblk)
  real(dp) :: resid_k(nband_k)
  real(dp) :: subham(nband_k*(nband_k+1))
  real(dp) :: subovl(nband_k*(nband_k+1)*use_subovl)
  real(dp) :: subvnl(nband_k*(nband_k+1)*(1-gs_hamk%usepaw))
  real(dp) :: vlocal(n4,n5,n6,nvloc)
  real(dp), intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)
 end subroutine lobpcgIIwf
end interface

interface
 subroutine lobpcgccIIIwf(dimffnl,dtfil,dtset,ffnl,gs_hamk,iterationnumber,&  
  &  kg_k,kinpw,lmnmax,matblk,mgfft,mpi_enreg,mpsang,&  
  &  mpssoang,natom,npw_k,ntypat,&  
  &  nvloc,n4,n5,n6,pcon,ph3d,prtvol,vlocal,&  
  &  blocksize,bblocksize,vectsize,pflag,&  
  &  blockvectorx,blockvectorbx,blockvectorax,blockvectory,blockvectorby,lambda,&  
  &  blockvectorp,blockvectorbp,blockvectorap,&  
  &  vxctaulocal) ! optional argument
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer :: bblocksize
  integer :: blocksize
  integer :: dimffnl
  integer :: iterationnumber
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
  integer :: ntypat
  integer :: nvloc
  integer :: prtvol
  integer :: vectsize
  type(datafiles_type) :: dtfil
  type(dataset_type) :: dtset
  type(gs_hamiltonian_type) :: gs_hamk
  real(dp) :: lambda
  type(mpi_type) :: mpi_enreg
  complex(dp) :: blockvectorap(vectsize,blocksize)
  complex(dp) :: blockvectorax(vectsize,blocksize)
  complex(dp) :: blockvectorbp(vectsize,blocksize)
  complex(dp) :: blockvectorbx(vectsize,blocksize)
  complex(dp) :: blockvectorby(vectsize,bblocksize)
  complex(dp) :: blockvectorp(vectsize,blocksize)
  complex(dp) :: blockvectorx(vectsize,blocksize)
  complex(dp) :: blockvectory(vectsize,bblocksize)
  real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  integer :: kg_k(3,npw_k)
  real(dp) :: kinpw(npw_k)
  real(dp) :: pcon(npw_k,blocksize)
  logical :: pflag(blocksize)
  real(dp) :: ph3d(2,npw_k,matblk)
  real(dp) :: vlocal(n4,n5,n6,nvloc)
  real(dp), intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)
 end subroutine lobpcgccIIIwf
end interface

interface
 subroutine lobpcgccIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&  
  &  kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&  
  &  nband_k,nbdblock,npw_k,ntypat,nvloc,n4,n5,n6,ph3d,prtvol,&  
  &  resid_k,subham,subovl,subvnl,use_subovl,vlocal,&  
  &  vxctaulocal) ! optional argument
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer :: dimffnl
  integer :: icg
  integer :: igsc
  integer :: lmnmax
  integer :: matblk
  integer :: mcg
  integer :: mgfft
  integer :: mgsc
  integer :: mpsang
  integer :: mpssoang
  integer :: n4
  integer :: n5
  integer :: n6
  integer :: natom
  integer :: nband_k
  integer :: nbdblock
  integer :: npw_k
  integer :: ntypat
  integer :: nvloc
  integer :: prtvol
  integer :: use_subovl
  type(datafiles_type) :: dtfil
  type(dataset_type) :: dtset
  type(gs_hamiltonian_type) :: gs_hamk
  type(mpi_type) :: mpi_enreg
  real(dp) :: cg(2,mcg)
  real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp) :: gsc(2,mgsc)
  integer :: kg_k(3,npw_k)
  real(dp) :: kinpw(npw_k)
  real(dp) :: ph3d(2,npw_k,matblk)
  real(dp) :: resid_k(nband_k)
  real(dp) :: subham(nband_k*(nband_k+1))
  real(dp) :: subovl(nband_k*(nband_k+1)*use_subovl)
  real(dp) :: subvnl(nband_k*(nband_k+1)*(1-gs_hamk%usepaw))
  real(dp) :: vlocal(n4,n5,n6,nvloc)
  real(dp), intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)
 end subroutine lobpcgccIIwf
end interface

interface
 subroutine lobpcgwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,ikpt,&  
  &  kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&  
  &  nband_k,nbdblock,npw_k,ntypat,nvloc,n4,n5,n6,ph3d,prtvol,&  
  &  resid_k,subham,totvnl,vlocal,&  
  &  vxctaulocal) ! optional argument
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer :: dimffnl
  integer :: icg
  integer :: igsc
  integer :: ikpt
  integer :: lmnmax
  integer :: matblk
  integer :: mcg
  integer :: mgfft
  integer :: mgsc
  integer :: mpsang
  integer :: mpssoang
  integer :: n4
  integer :: n5
  integer :: n6
  integer :: natom
  integer :: nband_k
  integer :: nbdblock
  integer :: npw_k
  integer :: ntypat
  integer :: nvloc
  integer :: prtvol
  type(datafiles_type) :: dtfil
  type(dataset_type) :: dtset
  type(gs_hamiltonian_type) :: gs_hamk
  type(mpi_type) :: mpi_enreg
  real(dp) :: cg(2,mcg)
  real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp) :: gsc(2,mgsc)
  integer :: kg_k(3,npw_k)
  real(dp) :: kinpw(npw_k)
  real(dp) :: ph3d(2,npw_k,matblk)
  real(dp) :: resid_k(nband_k)
  real(dp) :: subham(nband_k*(nband_k+1))
  real(dp) :: totvnl((3-gs_hamk%istwf_k)*nband_k*(1-gs_hamk%usepaw),nband_k*(1-gs_hamk%usepaw))
  real(dp) :: vlocal(n4,n5,n6,nvloc)
  real(dp), intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)
 end subroutine lobpcgwf
end interface

interface
 subroutine lobpcgIIIwf(dimffnl,dtfil,dtset,&  
  &  ffnl,gs_hamk,iterationnumber,&  
  &  kg_k,kinpw,lmnmax,matblk,mgfft,mpi_enreg,mpsang,&  
  &  mpssoang,natom,npw_k,ntypat,&  
  &  nvloc,n4,n5,n6,&  
  &  pcon,ph3d,prtvol,vlocal,&  
  &  blocksize,bblocksize,vectsize,pflag,&  
  &  blockvectorx,blockvectorbx,blockvectorax,blockvectorby,lambda,&  
  &  blockvectorp,blockvectorbp,blockvectorap,&  
  &  vxctaulocal) ! optional argument
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer :: bblocksize
  integer :: blocksize
  integer :: dimffnl
  integer :: iterationnumber
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
  integer :: ntypat
  integer :: nvloc
  integer :: prtvol
  integer :: vectsize
  type(datafiles_type) :: dtfil
  type(dataset_type) :: dtset
  type(gs_hamiltonian_type) :: gs_hamk
  type(mpi_type) :: mpi_enreg
  real(dp) :: blockvectorap(vectsize,blocksize)
  real(dp) :: blockvectorax(vectsize,blocksize)
  real(dp) :: blockvectorbp(vectsize,blocksize)
  real(dp) :: blockvectorbx(vectsize,blocksize)
  real(dp) :: blockvectorby(vectsize,bblocksize)
  real(dp) :: blockvectorp(vectsize,blocksize)
  real(dp) :: blockvectorx(vectsize,blocksize)
  real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  integer :: kg_k(3,npw_k)
  real(dp) :: kinpw(npw_k)
  real(dp) :: lambda(blocksize,blocksize)
  real(dp) :: pcon(npw_k,blocksize)
  logical :: pflag(blocksize)
  real(dp) :: ph3d(2,npw_k,matblk)
  real(dp) :: vlocal(n4,n5,n6,nvloc)
  real(dp), intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)
 end subroutine lobpcgIIIwf
end interface

interface
 subroutine mv_3dte(cg,cgindex,cg1,cg3,dtset,dtfil,d3_berry,gmet,&  
  &  i1pert,i3pert,i1dir,i3dir,kneigh,kg_neigh,kptindex,&  
  &  kpt3,mband,mkmem,mkmem_max,mk1mem,mpi_enreg,&  
  &  mpw,mvwtk,natom,nkpt2,nkpt3,nneigh,npwarr,nspinor,&  
  &  nsppol,pwind)
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: i1dir
  integer, intent(in) :: i1pert
  integer, intent(in) :: i3dir
  integer, intent(in) :: i3pert
  integer, intent(in) :: mband
  integer, intent(in) :: mk1mem
  integer, intent(in) :: mkmem
  integer, intent(in) :: mkmem_max
  integer, intent(in) :: mpw
  integer, intent(in) :: natom
  integer, intent(in) :: nkpt2
  integer, intent(in) :: nkpt3
  integer, intent(in) :: nneigh
  integer, intent(in) :: nspinor
  integer, intent(in) :: nsppol
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in) :: dtset
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp), intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp), intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
  real(dp), intent(in) :: cg3(2,mpw*nspinor*mband*mk1mem*nsppol)
  integer, intent(in) :: cgindex(nkpt2,nsppol)
  real(dp), intent(out) :: d3_berry(2,3)
  real(dp), intent(in) :: gmet(3,3)
  integer, intent(in) :: kg_neigh(30,nkpt2,3)
  integer, intent(in) :: kneigh(30,nkpt2)
  real(dp), intent(in) :: kpt3(3,nkpt3)
  integer, intent(in) :: kptindex(2,nkpt3)
  real(dp), intent(in) :: mvwtk(30,nkpt2)
  integer, intent(in) :: npwarr(nkpt2)
  integer, intent(in) :: pwind(mpw,nneigh,mkmem)
 end subroutine mv_3dte
end interface

interface
 subroutine subdiago(cg,eig_k,evec,gsc,icg,igsc,istwf_k,&  
  &  mcg,mgsc,mpi_enreg,nband_k,npw_k,nspinor,paral_kgb,&  
  &  subham,subovl,use_subovl,usepaw)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: icg
  integer,intent(in) :: igsc
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mcg
  integer,intent(in) :: mgsc
  integer,intent(in) :: nband_k
  integer,intent(in) :: npw_k
  integer,intent(in) :: nspinor
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: use_subovl
  integer,intent(in) :: usepaw
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cg(2,mcg)
  real(dp),intent(out) :: eig_k(nband_k)
  real(dp),intent(out) :: evec(2*nband_k,nband_k)
  real(dp),intent(inout) :: gsc(2,mgsc)
  real(dp),intent(inout) :: subham(nband_k*(nband_k+1))
  real(dp),intent(inout) :: subovl(nband_k*(nband_k+1)*use_subovl)
 end subroutine subdiago
end interface

interface
 subroutine tddft(cg,dtfil,dtset,eigen,etotal,gmet,gprimd,gsqcut,&  
  &  kg,kxc,mband,mgfftdiel,mkmem,mpi_enreg,mpw,nfft,ngfftdiel,nkpt,nkxc,&  
  &  npwarr,nspinor,nsppol,occ,ucvol,wffnew)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer, intent(in) :: mband
  integer, intent(in) :: mgfftdiel
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpw
  integer, intent(in) :: nfft
  integer, intent(in) :: nkpt
  integer, intent(in) :: nkxc
  integer, intent(in) :: nspinor
  integer, intent(in) :: nsppol
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in) :: dtset
  real(dp), intent(in) :: etotal
  real(dp), intent(in) :: gsqcut
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp), intent(in) :: ucvol
  type(wffile_type), intent(inout) :: wffnew
  integer, intent(in) :: ngfftdiel(18)
  real(dp), intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp), intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp), intent(in) :: gmet(3,3)
  real(dp), intent(in) :: gprimd(3,3)
  integer, intent(in) :: kg(3,mpw*mkmem)
  real(dp), intent(in) :: kxc(nfft,nkxc)
  integer, intent(in) :: npwarr(nkpt)
  real(dp), intent(in) :: occ(mband*nkpt*nsppol)
 end subroutine tddft
end interface

interface
 subroutine vtorho(afford,atindx,atindx1,cg,compch_fft,cpus,dbl_nnsclo,&  
  &  dielop,dielstrt,dphase,dtefield,dtfil,dtset,&  
  &  eigen,electronpositron,energies,etotal,gbound_diel,&  
  &  gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&  
  &  istep,istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&  
  &  mpsang,natom,nattyp,nfftf,nfftdiel,ngfftdiel,nhat,nkxc,&  
  &  npwarr,npwdiel,nres2,ntypat,nvresid,occ,optforces,&  
  &  optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&  
  &  pwind,pwind_alloc,pwnsfac,resid,residm,rhog,rhor,&  
  &  rmet,rprimd,susmat,symrec,taug,taur,&  
  &  ucvol,wffnew,wffnow,vtrial,wvl,xred,ylm,ylmgr,ylmdiel,&  
  &  vxctau) ! optional argument
  use defs_wvltypes
  use m_paw_dmft
  use m_energies
  use defs_abitypes
  use m_wffile
  use defs_basis
  use m_efield
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer, intent(in) :: afford
  integer, intent(in) :: dbl_nnsclo
  integer, intent(in) :: dielop
  integer, intent(in) :: dielstrt
  integer, intent(in) :: istep
  integer, intent(in) :: istep_mix
  integer, intent(in) :: lmax_diel
  integer, intent(in) :: mcg
  integer, intent(in) :: mgfftdiel
  integer, intent(in) :: mpsang
  integer, intent(in) :: natom
  integer, intent(in) :: nfftdiel
  integer, intent(in) :: nfftf
  integer, intent(in) :: nkxc
  integer, intent(in) :: npwdiel
  integer, intent(in) :: ntypat
  integer, intent(in) :: optforces
  integer, intent(in) :: optres
  integer, intent(in) :: pwind_alloc
  real(dp), intent(out) :: compch_fft
  real(dp), intent(in) :: cpus
  type(efield_type), intent(inout) :: dtefield
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in) :: dtset
  type(electronpositron_type),pointer :: electronpositron
  type(energies_type), intent(inout) :: energies
  real(dp), intent(in) :: etotal
  real(dp), intent(in) :: gsqcut
  type(hdr_type), intent(inout) :: hdr
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp), intent(out) :: nres2
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type), intent(in) :: pawang
  type(pawfgr_type), intent(in) :: pawfgr
  type(pseudopotential_type), intent(in) :: psps
  real(dp), intent(out) :: residm
  real(dp), intent(in) :: ucvol
  type(wffile_type), intent(inout) :: wffnew
  type(wffile_type), intent(inout) :: wffnow
  type(wvl_data), intent(inout) :: wvl
  integer, intent(in) :: ngfftdiel(18)
  integer, intent(in) :: atindx(natom)
  integer, intent(in) :: atindx1(natom)
  real(dp), intent(inout) :: cg(2,mcg)
  real(dp), intent(out) :: dphase(3)
  real(dp), intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer, intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp), intent(in) :: gmet(3,3)
  real(dp), intent(in) :: gprimd(3,3)
  real(dp), intent(out) :: grnl(3*natom)
  integer, intent(in) :: indsym(4,dtset%nsym,natom)
  integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2, &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym), &
  &         2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: kg_diel(3,npwdiel)
  real(dp), intent(inout) :: kxc(nfftf,nkxc)
  integer, intent(in) :: nattyp(ntypat)
  real(dp), intent(out) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(out) :: nvresid(nfftf,dtset%nspden)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_ij_type),intent(in) :: paw_ij(natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp), intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
  real(dp), intent(in) :: ph1ddiel(2,(3*(2*mgfftdiel+1)*natom)*psps%usepaw)
  real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp), intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(inout) :: rhog(2,nfftf)
  real(dp), intent(inout) :: rhor(nfftf,dtset%nspden)
  real(dp), intent(in) :: rmet(3,3)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(out) :: susmat(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden)
  integer, intent(in) :: symrec(3,3,dtset%nsym)
  real(dp), intent(inout) :: taug(2,nfftf*dtset%usekden)
  real(dp), intent(inout) :: taur(nfftf,dtset%nspden*dtset%usekden)
  real(dp), intent(inout) :: vtrial(nfftf,dtset%nspden)
  real(dp), intent(inout),optional :: vxctau(nfftf,dtset%nspden*dtset%usekden,4)
  real(dp), intent(inout) :: xred(3,natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,mpsang*mpsang*psps%useylm)
  real(dp), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine vtorho
end interface

interface
 subroutine vtowfk(cg,cgq,cprj,cpus,dimcprj,dimffnl,dphase_k,dtefield,dtfil,dtset,&  
  &  eig_k,ek_k,ek_k_nd,enl_k,fixed_occ,ffnl,grnl_k,gs_hamk,&  
  &  ibg,icg,ikpt,iscf,isppol,kg_k,kinpw,kpg_k,&  
  &  lmnmax,matblk,mband_cprj,mcg,mcgq,mcprj,mgfft,mkgq,mpi_enreg,mpsang,&  
  &  mpssoang,mpw,natom,nband_k,nkpg,nkpt,nnsclo_now,npw_k,npwarr,&  
  &  ntypat,nvloc,n4,n5,n6,occ_k,optforces,ph3d,prtvol,&  
  &  pwind,pwind_alloc,pwnsfac,pwnsfacq,resid_k,rhoaug,use_dmft,usecprj,vlocal,wtk,zshift,&  
  &  vxctaulocal) ! optional argument
  use defs_basis
  use defs_abitypes
  use m_efield
  use defs_datatypes
  implicit none
  integer, intent(in) :: dimffnl
  integer, intent(in) :: ibg
  integer, intent(in) :: icg
  integer, intent(in) :: ikpt
  integer, intent(in) :: iscf
  integer, intent(in) :: isppol
  integer, intent(in) :: lmnmax
  integer, intent(in) :: matblk
  integer, intent(in) :: mband_cprj
  integer, intent(in) :: mcg
  integer, intent(in) :: mcgq
  integer, intent(in) :: mcprj
  integer, intent(in) :: mgfft
  integer, intent(in) :: mkgq
  integer, intent(in) :: mpsang
  integer, intent(in) :: mpssoang
  integer, intent(in) :: mpw
  integer, intent(in) :: n4
  integer, intent(in) :: n5
  integer, intent(in) :: n6
  integer, intent(in) :: natom
  integer, intent(in) :: nband_k
  integer, intent(in) :: nkpg
  integer, intent(in) :: nkpt
  integer, intent(in) :: nnsclo_now
  integer, intent(in) :: npw_k
  integer, intent(in) :: ntypat
  integer, intent(in) :: nvloc
  integer, intent(in) :: optforces
  integer, intent(in) :: prtvol
  integer, intent(in) :: pwind_alloc
  integer, intent(in) :: use_dmft
  integer, intent(in) :: usecprj
  real(dp), intent(in) :: cpus
  type(efield_type), intent(inout) :: dtefield
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in) :: dtset
  logical,intent(in) :: fixed_occ
  type(gs_hamiltonian_type), intent(inout) :: gs_hamk
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp), intent(in) :: wtk
  real(dp), intent(inout) :: cg(2,mcg)
  real(dp), intent(in) :: cgq(2,mcgq)
  type(cprj_type) :: cprj(natom,mcprj*usecprj)
  integer, intent(in) :: dimcprj(natom*gs_hamk%usepaw)
  real(dp), intent(out) :: dphase_k(3)
  real(dp), intent(out) :: eig_k(nband_k)
  real(dp), intent(out) :: ek_k(nband_k)
  real(dp), intent(out) :: ek_k_nd(nband_k,nband_k*use_dmft)
  real(dp), intent(out) :: enl_k(nband_k*(1-gs_hamk%usepaw))
  real(dp), intent(in) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
  real(dp), intent(out) :: grnl_k(3*natom,nband_k*optforces)
  integer, intent(in) :: kg_k(3,npw_k)
  real(dp), intent(in) :: kinpw(npw_k)
  real(dp), intent(in) :: kpg_k(npw_k,nkpg)
  integer, intent(in) :: npwarr(nkpt)
  real(dp), intent(in) :: occ_k(nband_k)
  real(dp), intent(inout) :: ph3d(2,npw_k,matblk)
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(in) :: pwnsfacq(2,mkgq)
  real(dp), intent(out) :: resid_k(nband_k)
  real(dp), intent(inout) :: rhoaug(n4,n5,n6,nvloc)
  real(dp), intent(inout) :: vlocal(n4,n5,n6,nvloc)
  real(dp), intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)
  real(dp), intent(in) :: zshift(nband_k)
 end subroutine vtowfk
end interface

interface
 subroutine wfsinp(cg,cg_disk,ecut,ecut0,ecut_eff,eigen,exchn2n3d,&  
  &  formeig,gmet,gmet0,headform0,indkk,indkk0,istwfk,&  
  &  istwfk0,kptns,kptns0,localrdwf,mband,mban_dp_rd,&  
  &  mcg,mcg_disk,mkmem,mpi_enreg,mpi_enreg0,mpw,mpw0,nband,nban_dp_rd,&  
  &  ngfft,nkassoc,nkpt,nkpt0,npwarr,npwarr0,nspinor,&  
  &  nspinor0,nsppol,nsppol0,nsym,occ,optorth,prtvol,restart,rprimd,&  
  &  sppoldbl,squeeze,symrel,tnons,wff1,wffnow)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer, intent(in) :: exchn2n3d
  integer, intent(in) :: formeig
  integer, intent(in) :: headform0
  integer, intent(in) :: localrdwf
  integer, intent(in) :: mban_dp_rd
  integer, intent(in) :: mband
  integer, intent(in) :: mcg
  integer, intent(in) :: mcg_disk
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpw
  integer, intent(in) :: mpw0
  integer, intent(in) :: nkassoc
  integer, intent(in) :: nkpt
  integer, intent(in) :: nkpt0
  integer, intent(in) :: nspinor
  integer, intent(in) :: nspinor0
  integer, intent(in) :: nsppol
  integer, intent(in) :: nsppol0
  integer, intent(in) :: nsym
  integer, intent(in) :: optorth
  integer, intent(in) :: prtvol
  integer, intent(in) :: restart
  integer, intent(in) :: sppoldbl
  integer, intent(in) :: squeeze
  real(dp), intent(in) :: ecut
  real(dp), intent(in) :: ecut0
  real(dp), intent(in) :: ecut_eff
  type(mpi_type), intent(inout) :: mpi_enreg
  type(mpi_type), intent(inout) :: mpi_enreg0
  type(wffile_type), intent(inout) :: wff1
  type(wffile_type), intent(inout) :: wffnow
  integer, intent(in) :: ngfft(18)
  real(dp), intent(out) :: cg(2,mcg)
  real(dp), intent(out) :: cg_disk(2,mcg_disk)
  real(dp), intent(out) :: eigen((2*mband)**formeig*mband*nkpt*nsppol)
  real(dp), intent(in) :: gmet(3,3)
  real(dp), intent(in) :: gmet0(3,3)
  integer, intent(in) :: indkk(nkpt*sppoldbl,6)
  integer, intent(in) :: indkk0(nkpt0,nkassoc)
  integer, intent(in) :: istwfk(nkpt)
  integer, intent(in) :: istwfk0(nkpt0)
  real(dp), intent(in) :: kptns(3,nkpt)
  real(dp), intent(in) :: kptns0(3,nkpt0)
  integer, intent(in) :: nban_dp_rd(nkpt0*nsppol0)
  integer, intent(in) :: nband(nkpt*nsppol)
  integer, intent(in) :: npwarr(nkpt)
  integer, intent(in) :: npwarr0(nkpt0)
  real(dp), intent(inout) :: occ(mband*nkpt*nsppol)
  real(dp), intent(in) :: rprimd(3,3)
  integer, intent(in) :: symrel(3,3,nsym)
  real(dp), intent(in) :: tnons(3,nsym)
 end subroutine wfsinp
end interface

interface
 subroutine wvl_vtorho(dtset, energies, irrzon, istep, mpi_enreg,&  
  &  phnons, residm, rhor, rprimd, vtrial, wvl, xred)
  use defs_basis
  use m_energies
  use defs_abitypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: istep
  type(dataset_type), intent(in) :: dtset
  type(energies_type), intent(inout) :: energies
  type(mpi_type), intent(in) :: mpi_enreg
  real(dp), intent(inout) :: residm
  type(wvl_data), intent(inout) :: wvl
  integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2, &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp), intent(inout) :: rhor(dtset%nfft)
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(in) :: vtrial(dtset%nfft * dtset%nspden)
  real(dp), intent(inout) :: xred(3, dtset%natom)
 end subroutine wvl_vtorho
end interface

interface
 subroutine wvl_wfsinp_disk(dtset, hdr0, hdr, mpi_enreg, option,&  
  &  rprimd, wff, wfs, wvl, xred)
  use defs_basis
  use defs_abitypes
  use m_wffile
  use defs_wvltypes
  implicit none
  integer, intent(in) :: option
  type(dataset_type), intent(in) :: dtset
  type(hdr_type), intent(in) :: hdr
  type(hdr_type), intent(in) :: hdr0
  type(mpi_type), intent(in) :: mpi_enreg
  type(wffile_type), intent(in) :: wff
  type(wvl_wf_type), intent(inout) :: wfs
  type(wvl_internal_type), intent(in) :: wvl
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(inout) :: xred(3, dtset%natom)
 end subroutine wvl_wfsinp_disk
end interface

interface
 subroutine wvl_wfsinp_reformat(dtset, mpi_enreg, psps,&  
  &  rprimd, wvl, xred, xred_old)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  implicit none
  type(dataset_type), intent(inout) :: dtset
  type(mpi_type), intent(inout) :: mpi_enreg
  type(pseudopotential_type), intent(in) :: psps
  type(wvl_data), intent(inout) :: wvl
  real(dp), intent(inout) :: rprimd(3,3)
  real(dp), intent(inout) :: xred(3, dtset%natom)
  real(dp), intent(inout) :: xred_old(3, dtset%natom)
 end subroutine wvl_wfsinp_reformat
end interface

interface
 subroutine wvl_wfsinp_scratch(dtset, mpi_enreg, rprimd, wvl, xred)
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  implicit none
  type(dataset_type), intent(in) :: dtset
  type(mpi_type), intent(in) :: mpi_enreg
  type(wvl_data), intent(inout) :: wvl
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(inout) :: xred(3, dtset%natom)
 end subroutine wvl_wfsinp_scratch
end interface

end module interfaces_79_seqpar_mpi
!!***
