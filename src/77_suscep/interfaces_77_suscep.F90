!!****m* ABINIT/interfaces_77_suscep
!! NAME
!! interfaces_77_suscep
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/77_suscep
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

module interfaces_77_suscep

 implicit none

interface
 subroutine acfd_dyson(dtset,freq,gsq,idyson,ig_tiny,igsq_tiny,ikhxc,kg_diel,khxc,&  
  &  ldgapp,mpi_enreg,ndyson,nfft,ngfft,npw_tiny,npwdiel,nspden,&  
  &  option,rcut_coulomb,rhor,rhocut,rprimd,susd_data,suskxcrs,susmat,ucvol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: idyson
  integer,intent(in) :: ikhxc
  integer,intent(in) :: ldgapp
  integer,intent(in) :: ndyson
  integer,intent(in) :: nfft
  integer,intent(in) :: npw_tiny
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: suskxcrs
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: freq
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: rcut_coulomb
  real(dp),intent(in) :: rhocut
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gsq(npwdiel)
  integer,intent(in) :: ig_tiny(npw_tiny,3)
  integer,intent(in) :: igsq_tiny(npw_tiny)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),pointer :: khxc(:,:,:,:,:)
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: susd_data(npwdiel,3)
  real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine acfd_dyson
end interface

interface
 subroutine acfd_intexact(dec,freq,gsq,ikhxc,mband,nband,nkpt,npwdiel,&  
  &  nspden,nsppol,occ,occopt,rcut_coulomb,susmat)
  use defs_basis
  implicit none
  integer,intent(in) :: ikhxc
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  real(dp),intent(out) :: dec
  real(dp),intent(in) :: freq
  real(dp),intent(in) :: rcut_coulomb
  real(dp),intent(in) :: gsq(npwdiel)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine acfd_intexact
end interface

interface
 subroutine dyson_de(ikernel,kernel_diag,kernel_full,npwdiel,nspden,susmat)
  use defs_basis
  implicit none
  integer,intent(in) :: ikernel
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  real(dp),pointer :: kernel_diag(:)
  real(dp),pointer :: kernel_full(:,:,:,:,:)
  real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine dyson_de
end interface

interface
 subroutine dyson_gl(ikernel,kernel_diag,kernel_full,npwdiel,nspden,&  
  &  susd_LDG,susmat)
  use defs_basis
  implicit none
  integer,intent(in) :: ikernel
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  real(dp),pointer :: kernel_diag(:)
  real(dp),pointer :: kernel_full(:,:,:,:,:)
  real(dp),intent(out) :: susd_LDG(npwdiel)
  real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine dyson_gl
end interface

interface
 subroutine dyson_ls(ikernel,kernel_diag,kernel_full,npwdiel,nspden,susmat)
  use defs_basis
  implicit none
  integer,intent(in) :: ikernel
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  real(dp),pointer :: kernel_diag(:)
  real(dp),pointer :: kernel_full(:,:,:,:,:)
  real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine dyson_ls
end interface

interface
 subroutine dyson_sc(kernel_diag,npwdiel,nspden,susd_isc,susmat)
  use defs_basis
  implicit none
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  real(dp),intent(in) :: kernel_diag(npwdiel)
  real(dp),intent(out) :: susd_isc(npwdiel)
  real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 end subroutine dyson_sc
end interface

interface
 subroutine get_g_tiny(gmet,gprimd,gsq_unsorted,ig_tiny,igsq_tiny,index_g,kg,npw,npw_tiny)
  use defs_basis
  implicit none
  integer,intent(in) :: npw
  integer,intent(in) :: npw_tiny
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: gsq_unsorted(npw)
  integer,intent(out) :: ig_tiny(npw_tiny,3)
  integer,intent(out) :: igsq_tiny(npw_tiny)
  integer,intent(out) :: index_g(npw)
  integer,intent(in) :: kg(3,npw)
 end subroutine get_g_tiny
end interface

interface
 subroutine get_susd_null(ig_tiny,igsq_tiny,gsq_input,npwdiel,npw_tiny,&  
  &  sus_gabs,sus_gavg,sus_gdir,sus_input)
  use defs_basis
  implicit none
  integer,intent(in) :: npw_tiny
  integer,intent(in) :: npwdiel
  real(dp),intent(in) :: gsq_input(npwdiel)
  integer,intent(in) :: ig_tiny(npw_tiny,3)
  integer,intent(in) :: igsq_tiny(npw_tiny)
  real(dp),intent(out) :: sus_gabs(npw_tiny)
  real(dp),intent(out) :: sus_gavg(npw_tiny)
  real(dp),intent(out) :: sus_gdir(npw_tiny,3)
  real(dp),intent(in) :: sus_input(npwdiel)
 end subroutine get_susd_null
end interface

interface
 subroutine geteexc_cc(energy,energy_raw,gsq,npwdiel,npw_tiny,rcut_coulomb,susd)
  use defs_basis
  implicit none
  integer,intent(in) :: npw_tiny
  integer,intent(in) :: npwdiel
  real(dp),intent(out) :: energy_raw
  real(dp),intent(in) :: rcut_coulomb
  real(dp),intent(out) :: energy(npw_tiny)
  real(dp),intent(in) :: gsq(npwdiel)
  real(dp),intent(in) :: susd(npwdiel)
 end subroutine geteexc_cc
end interface

interface
 subroutine geteexc_uc(energy,energy_raw,gsq,ig_tiny,npwdiel,npw_tiny,&  
  &  susd)
  use defs_basis
  implicit none
  integer,intent(in) :: npw_tiny
  integer,intent(in) :: npwdiel
  real(dp),intent(out) :: energy_raw
  real(dp),intent(out) :: energy(npw_tiny)
  real(dp),intent(in) :: gsq(npwdiel)
  integer,intent(in) :: ig_tiny(npw_tiny,3)
  real(dp),intent(in) :: susd(npwdiel)
 end subroutine geteexc_uc
end interface

interface
 subroutine getfreqsus(freqs,weights,nfreqs,optfreq,freq1,freq2)
  use defs_basis
  implicit none
  integer,intent(in) :: nfreqs
  integer,intent(in) :: optfreq
  real(dp),intent(in) :: freq1
  real(dp),intent(in) :: freq2
  real(dp),intent(out) :: freqs(nfreqs)
  real(dp),intent(out) :: weights(nfreqs)
 end subroutine getfreqsus
end interface

interface
 subroutine getlambda(idyson,lambda,nlambda,weight)
  use defs_basis
  implicit none
  integer,intent(in) :: idyson
  integer,intent(in) :: nlambda
  real(dp),intent(out) :: lambda(nlambda)
  real(dp),intent(out) :: weight(nlambda)
 end subroutine getlambda
end interface

interface
 subroutine inwffil3(dtset,eigen,hdr,istwfk,mband,mpi_enreg,nband,&  
  &  nkpt,npwarr,nsppol,prtvol,wff1,unwff1,wffnm)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: prtvol
  integer,intent(in) :: unwff1
  type(dataset_type),intent(in) :: dtset
  type(hdr_type),intent(out) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(out) :: wff1
  character(len=fnlen),intent(in) :: wffnm
  real(dp),intent(out) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
 end subroutine inwffil3
end interface

interface
 subroutine k_rpa(gsq,krpa,npw,option,rcut_coulomb)
  use defs_basis
  implicit none
  integer,intent(in) :: npw
  integer,intent(in) :: option
  real(dp),intent(in) :: rcut_coulomb
  real(dp),intent(in) :: gsq(npw)
  real(dp),intent(out) :: krpa(npw)
 end subroutine k_rpa
end interface

interface
 subroutine klocal(ispxc,kg_diel,kxc,kxcg,nfft,ngfft,npwdiel,nspden,option)
  use defs_basis
  implicit none
  integer,intent(in) :: ispxc
  integer,intent(in) :: nfft
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(out) :: kxc(2,npwdiel,nspden,npwdiel,nspden)
  real(dp),intent(in) :: kxcg(2,nfft)
 end subroutine klocal
end interface

interface
 subroutine kxc_alda(dtset,ixc,kxcg,mpi_enreg,nfft,ngfft,nspden,option,rhor,rhocut,rprimd)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: rhocut
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: kxcg(2,nfft,*)
  real(dp),intent(in) :: rhor(nfft,2*nspden-1)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine kxc_alda
end interface

interface
 subroutine kxc_eok(ixceok,kxcg,mpi_enreg,nfft,ngfft,nspden,paral_kgb,rhor,rhocut)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: ixceok
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  type(mpi_type) :: mpi_enreg
  real(dp),intent(in) :: rhocut
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: kxcg(2,nfft,2*nspden-1)
  real(dp),intent(in) :: rhor(nfft,2*nspden-1)
 end subroutine kxc_eok
end interface

interface
 subroutine kxc_pgg(gmet,kg,khxcg,npw,rcut_coulomb,susmat,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: npw
  real(dp),intent(in) :: rcut_coulomb
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg(3,npw)
  real(dp),intent(out) :: khxcg(2,npw,npw)
  real(dp),intent(in) :: susmat(2,npw,npw)
 end subroutine kxc_pgg
end interface

interface
 subroutine prtsusd(gsq,ig_tiny,index_g,npw_tiny,npwdiel,optprt,susd)
  use defs_basis
  implicit none
  integer,intent(in) :: npw_tiny
  integer,intent(in) :: npwdiel
  integer,intent(in) :: optprt
  real(dp),intent(in) :: gsq(npwdiel)
  integer,intent(in) :: ig_tiny(npw_tiny,3)
  integer,intent(in) :: index_g(npwdiel)
  real(dp),intent(in) :: susd(npwdiel)
 end subroutine prtsusd
end interface

interface
 subroutine suscep(dtfil,dtset,iexit,&  
  &  mband,mkmem,mpi_enreg,mpw,natom,nfft,nkpt,&  
  &  nspden,nspinor,nsppol,nsym,occ,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: mband
  integer,intent(inout) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(inout) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine suscep
end interface

interface
 subroutine suscep_dyn(dielar,dtset,&  
  &  eigen,freq,gbound_diel,gprimd,irrzondiel,istwfk,kg,kg_diel,&  
  &  mband,mgfftdiel,mkmem,mpi_enreg,mpw,nband,nband_mx,nfftdiel,nfreq,&  
  &  ngfftdiel,nkpt,npwarr,&  
  &  npwdiel,nspden,nspinor,nsppol,nsym,occ,occopt,phnonsdiel,rprimd,&  
  &  susopt,sus_diag_dyn,susmat_dyn,symafm,symrel,tnons,ucvol,unkg,wff1,wtk)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nband_mx
  integer,intent(in) :: nfftdiel
  integer,intent(in) :: nfreq
  integer,intent(in) :: nkpt
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(inout) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: occopt
  integer,intent(in) :: susopt
  integer,intent(in) :: unkg
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wff1
  integer,intent(in) :: ngfftdiel(18)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: freq(nfreq)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: irrzondiel(nfftdiel**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: phnonsdiel(2,nfftdiel**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: sus_diag_dyn(2,npwdiel,nspden,nfreq)
  real(dp),intent(out) :: susmat_dyn(2,npwdiel,nspden,npwdiel,nspden,nfreq)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine suscep_dyn
end interface

interface
 subroutine suscep_kxc_dyn(cg,dielar,dielop,doccde,dtset,&  
  &  eigen,freq,gbound_diel,gprimd,irrzondiel,istwfk,kg,kg_diel,&  
  &  mband,mgfftdiel,mkmem,mpi_enreg,mpw,nband,nband_mx,nfftdiel,nfreq,&  
  &  ngfftdiel,nkpt,npwarr,&  
  &  npwdiel,nspden,nspinor,nsppol,nsym,occ,occopt,phnonsdiel,prtvol,rprimd,&  
  &  susopt,sus_diag_dyn,susmat_dyn,symafm,symrel,tnons,ucvol,unkg,wff1,wtk)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: dielop
  integer,intent(in) :: mband
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nband_mx
  integer,intent(in) :: nfftdiel
  integer,intent(in) :: nfreq
  integer,intent(in) :: nkpt
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(inout) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: occopt
  integer,intent(in) :: prtvol
  integer,intent(in) :: susopt
  integer,intent(in) :: unkg
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wff1
  integer,intent(in) :: ngfftdiel(18)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: freq(nfreq)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: irrzondiel(nfftdiel**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: phnonsdiel(2,nfftdiel**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: sus_diag_dyn(2,npwdiel,nspden,nfreq)
  real(dp),intent(out) :: susmat_dyn(2,npwdiel,nspden,npwdiel,nspden,nfreq)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine suscep_kxc_dyn
end interface

interface
 subroutine suscep_stat(atindx1,cg,cprj,dielar,dimcprj,doccde,&  
  &  eigen,gbound_diel,gprimd,irrzondiel,istwfk,kg,&  
  &  kg_diel,lmax_diel,&  
  &  mband,mcg,mcprj,mgfftdiel,mkmem,mpi_enreg,mpw,natom,nband,&  
  &  neglect_pawhat,nfftdiel,ngfftdiel,nkpt,npwarr,&  
  &  npwdiel,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&  
  &  pawang,pawtab,phnonsdiel,ph1ddiel,rprimd,&  
  &  susmat,symafm,symrel,tnons,typat,ucvol,unkg,unpaw,usecprj,usepaw,usetimerev,&  
  &  wffnew,wtk,ylmdiel)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: lmax_diel
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mcprj
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: neglect_pawhat
  integer,intent(in) :: nfftdiel
  integer,intent(in) :: nkpt
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: unkg
  integer,intent(in) :: unpaw
  integer,intent(in) :: usecprj
  integer,intent(in) :: usepaw
  integer,intent(in) :: usetimerev
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnew
  integer,intent(in) :: ngfftdiel(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mcg)
  type(cprj_type) :: cprj(natom,mcprj*usecprj)
  real(dp),intent(in) :: dielar(7)
  integer,intent(in) :: dimcprj(natom*usepaw)
  real(dp),intent(in) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: irrzondiel(nfftdiel**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph1ddiel(2,(3*(2*mgfftdiel+1)*natom)*usepaw)
  real(dp),intent(in) :: phnonsdiel(2,nfftdiel**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(ntypat)
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
 end subroutine suscep_stat
end interface

interface
 subroutine susk(atindx1,bdtot_index,cg_mpi,cprj_k,doccde,drhode,eigen,extrap,gbound,&  
  &  gbound_diel,gylmg_diel,icg_mpi,ikpt,isp,istwfk,kg_diel,kg_k_mpi,&  
  &  lmax_diel,mband,mcg,mgfftdiel,mpi_enreg,&  
  &  natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat,ngfftdiel,nkpt,&  
  &  npwdiel,npw_k_mpi,nspden,nspden_eff,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,&  
  &  pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&  
  &  susmat,typat,ucvol,usepaw,wtk)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: bdtot_index
  integer,intent(in) :: extrap
  integer,intent(in),target :: icg_mpi
  integer,intent(in) :: ikpt
  integer,intent(in) :: isp
  integer,intent(in) :: lmax_diel
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: natom
  integer,intent(in) :: nband_k
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: neglect_pawhat
  integer,intent(in) :: nkpt
  integer,intent(in),target :: npw_k_mpi
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: nspden_eff
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: usepaw
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(inout) :: sumdocc
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfftdiel(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in),target :: cg_mpi(2,mcg)
  type(cprj_type) :: cprj_k(natom,nspinor*nband_k*usepaw)
  real(dp),intent(in) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(inout) :: drhode(2,npwdiel,nspden_eff)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(inout) :: gbound(2*mgfftdiel+8,2)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat*usepaw)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in),target :: kg_k_mpi(3,npw_k_mpi)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: occ_deavg(mband)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom*usepaw)
  real(dp),intent(inout) :: rhoextrap(ndiel4,ndiel5,ndiel6,nspinor)
  real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine susk
end interface

interface
 subroutine susk_dyn(bdtot_index,cg,dtset,eigen,extrap,freq,&  
  &  gbound,gbound_diel,icg,ikpt,isp,istwfk,kg_diel,kg_k,&  
  &  mband,mcg,mgfftdiel,mpi_enreg,mpw,&  
  &  nband_k,ndiel4,ndiel5,ndiel6,nfreq,ngfftdiel,nkpt,&  
  &  npwdiel,npw_k,nsppol,occ,occopt,&  
  &  occ_deavg,occ_freq,rhoextrap_dyn,susmat_dyn,susopt,ucvol,wtk)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bdtot_index
  integer,intent(in) :: extrap
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isp
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: mpw
  integer,intent(in) :: nband_k
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: nfreq
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw_k
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: susopt
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfftdiel(18)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: freq(nfreq)
  integer,intent(in) :: gbound(2*mgfftdiel+8,2)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: occ_deavg(mband)
  real(dp),intent(in) :: occ_freq(2,mband,nfreq)
  real(dp),intent(inout) :: rhoextrap_dyn(2,ndiel4,ndiel5,ndiel6,nfreq)
  real(dp),intent(inout) :: susmat_dyn(2,npwdiel,nsppol,npwdiel,nsppol,nfreq)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine susk_dyn
end interface

interface
 subroutine susk_dyn_pgg(bdtot_index,cg,eigen,extrap,freq,&  
  &  gbound,gbound_diel,icg,ikpt,isp,istwfk,kg_diel,kg_k,&  
  &  mband,mcg,mgfftdiel,mpi_enreg,mpw,&  
  &  nband_k,ndiel4,ndiel5,ndiel6,nfreq,ngfftdiel,nkpt,&  
  &  npwdiel,npw_k,nsppol,occ,occopt,&  
  &  occ_deavg,occ_freq,paral_kgb,rhoextrap_dyn,susmat_dyn,susopt,ucvol,wtk)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bdtot_index
  integer,intent(in) :: extrap
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isp
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: mpw
  integer,intent(in) :: nband_k
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: nfreq
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw_k
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: susopt
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfftdiel(18)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: freq(nfreq)
  integer,intent(in) :: gbound(2*mgfftdiel+8,2)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: occ_deavg(mband)
  real(dp),intent(in) :: occ_freq(2,mband,nfreq)
  real(dp),intent(inout) :: rhoextrap_dyn(2,ndiel4,ndiel5,ndiel6,nfreq)
  real(dp),intent(inout) :: susmat_dyn(2,npwdiel,nsppol,npwdiel,nsppol,nfreq)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine susk_dyn_pgg
end interface

interface
 subroutine susk_kxc_dyn(bdtot_index,cg,doccde,drhode,eigen,extrap,freq,&  
  &  gbound,gbound_diel,icg,ikpt,isp,istwfk,kg_diel,kg_k,&  
  &  mband,mcg,mgfftdiel,mkmem,mpi_enreg,mpw,&  
  &  nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,nfreq,ngfftdiel,nkpt,&  
  &  npwdiel,npw_k,nspden,nspinor,nsppol,occ,occopt,&  
  &  occ_deavg,occ_freq,rhoextrap_dyn,sumdocc,susmat_dyn,susopt,ucvol,wtk,kxc)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bdtot_index
  integer,intent(in) :: extrap
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isp
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nband_k
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: nfftdiel
  integer,intent(in) :: nfreq
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw_k
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: susopt
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: sumdocc
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfftdiel(18)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(inout) :: drhode(2,npwdiel,nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: freq(nfreq)
  integer,intent(in) :: gbound(2*mgfftdiel+8,2)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: kxc(:,:,:,:)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: occ_deavg(mband)
  real(dp),intent(in) :: occ_freq(2,mband,nfreq)
  real(dp),intent(inout) :: rhoextrap_dyn(2,ndiel4,ndiel5,ndiel6,nfreq)
  real(dp),intent(inout) :: susmat_dyn(2,npwdiel,nsppol,npwdiel,nsppol,nfreq)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine susk_kxc_dyn
end interface

interface
 subroutine suskmm(atindx1,bdtot_index,cg,cprj_k,doccde,drhode,eigen,extrap,gbound,&  
  &  gbound_diel,gylmg_diel,icg,ikpt,isp,istwfk,kg_diel,kg_k,&  
  &  lmax_diel,mband,mcg,mgfftdiel,mpi_enreg,&  
  &  natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat,ngfftdiel,nkpt,&  
  &  npwdiel,npw_k,nspden,nspden_eff,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,paral_kgb,&  
  &  pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&  
  &  susmat,typat,ucvol,usepaw,wtk)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: bdtot_index
  integer,intent(in) :: extrap
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isp
  integer,intent(in) :: lmax_diel
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: natom
  integer,intent(in) :: nband_k
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: neglect_pawhat
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw_k
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: nspden_eff
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: usepaw
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  real(dp),intent(inout) :: sumdocc
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfftdiel(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mcg)
  type(cprj_type) :: cprj_k(natom,nspinor*nband_k*usepaw)
  real(dp),intent(in) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(inout) :: drhode(2,npwdiel,nspden_eff)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: gbound(2*mgfftdiel+8,2)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat*usepaw)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: occ_deavg(mband)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom*usepaw)
  real(dp),intent(inout) :: rhoextrap(ndiel4,ndiel5,ndiel6,nspinor)
  real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine suskmm
end interface

interface
 subroutine suskmm_dyn(bdtot_index,cg,eigen,extrap,freq,&  
  &  gbound,gbound_diel,icg,ikpt,isp,istwfk,kg_diel,kg_k,&  
  &  mband,mcg,mgfftdiel,mpi_enreg,mpw,&  
  &  nband_k,ndiel4,ndiel5,ndiel6,nfreq,ngfftdiel,nkpt,&  
  &  npwdiel,npw_k,nsppol,occ,occopt,&  
  &  occ_deavg,occ_freq,paral_kgb,rhoextrap_dyn,susmat_dyn,susopt,ucvol,wtk)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bdtot_index
  integer,intent(in) :: extrap
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isp
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: mpw
  integer,intent(in) :: nband_k
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: nfreq
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw_k
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: susopt
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfftdiel(18)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: freq(nfreq)
  integer,intent(in) :: gbound(2*mgfftdiel+8,2)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: occ_deavg(mband)
  real(dp),intent(in) :: occ_freq(2,mband,nfreq)
  real(dp),intent(inout) :: rhoextrap_dyn(2,ndiel4,ndiel5,ndiel6,nfreq)
  real(dp),intent(inout) :: susmat_dyn(2,npwdiel,nsppol,npwdiel,nsppol,nfreq)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine suskmm_dyn
end interface

interface
 subroutine suskmm_kxc_dyn(bdtot_index,cg,doccde,drhode,eigen,extrap,freq,&  
  &  gbound,gbound_diel,icg,ikpt,isp,istwfk,kg_diel,kg_k,&  
  &  mband,mcg,mgfftdiel,mkmem,mpi_enreg,mpw,&  
  &  nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,nfreq,ngfftdiel,nkpt,&  
  &  npwdiel,npw_k,nspden,nspinor,nsppol,occ,occopt,&  
  &  occ_deavg,occ_freq,paral_kgb,rhoextrap_dyn,sumdocc,susmat_dyn,susopt,ucvol,wtk,&  
  &  kxc)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: bdtot_index
  integer,intent(in) :: extrap
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isp
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nband_k
  integer,intent(in) :: ndiel4
  integer,intent(in) :: ndiel5
  integer,intent(in) :: ndiel6
  integer,intent(in) :: nfftdiel
  integer,intent(in) :: nfreq
  integer,intent(in) :: nkpt
  integer,intent(in) :: npw_k
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: susopt
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: sumdocc
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfftdiel(18)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(inout) :: drhode(2,npwdiel,nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: freq(nfreq)
  integer,intent(in) :: gbound(2*mgfftdiel+8,2)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(in) :: kxc(:,:,:,:)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: occ_deavg(mband)
  real(dp),intent(in) :: occ_freq(2,mband,nfreq)
  real(dp),intent(inout) :: rhoextrap_dyn(2,ndiel4,ndiel5,ndiel6,nfreq)
  real(dp),intent(inout) :: susmat_dyn(2,npwdiel,nsppol,npwdiel,nsppol,nfreq)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine suskmm_kxc_dyn
end interface

interface
 subroutine xcacfd(dielar,&  
  &  dtfil,dtset,eigen,freq,gbound_diel,gmet,&  
  &  gprimd,irrzondiel,kg,kg_diel,mband,mgfftdiel,&  
  &  mkmem,mpi_enreg,mpw,nfft,nfftdiel,nfreqsus,ngfft,&  
  &  ngfftdiel,nkpt,npwarr,npwdiel,nspden,nspinor,nsppol,&  
  &  nsym,occ,phnonsdiel,rhor,rprimd,ucvol,wff1,wght_freq)
  use defs_basis
  use defs_abitypes
  use m_wffile
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfftdiel
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftdiel
  integer,intent(in) :: nfreqsus
  integer,intent(in) :: nkpt
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nspden
  integer,intent(inout) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wff1
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftdiel(18)
  real(dp),intent(inout) :: dielar(7)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: freq(nfreqsus)
  integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: irrzondiel(nfftdiel**(1-1/nsym),2+(nspden/4),(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: kg_diel(3,npwdiel)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: phnonsdiel(2,nfftdiel**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(in) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: wght_freq(nfreqsus)
 end subroutine xcacfd
end interface

end module interfaces_77_suscep
!!***
