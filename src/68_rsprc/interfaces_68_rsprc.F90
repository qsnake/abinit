!!****m* ABINIT/interfaces_68_rsprc
!! NAME
!! interfaces_68_rsprc
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/68_rsprc
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

module interfaces_68_rsprc

 implicit none

interface
 subroutine ladielmt(atindx,atindx1,cg,dielmat,dtfil,dtset,&  
  &  eigen,&  
  &  kg,kg_diel,mpi_enreg,&  
  &  nattyp,&  
  &  npwarr,npwdiel,nspinor,&  
  &  occ,&  
  &  ph1d,psps,rhor,rprimd,&  
  &  wffnow,xred,ylm,ladiel)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer, intent(in) :: npwdiel
  integer,intent(inout) :: nspinor
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(in) :: cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp),intent(in) :: dielmat(2,npwdiel,npwdiel)
  real(dp), intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(out) :: ladiel(dtset%nfft)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp), intent(in) :: rhor(dtset%nfft,dtset%nspden)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(in) :: xred(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine ladielmt
end interface

interface
 subroutine lavnl(atindx,atindx1,cg,dtfil,dtset,&  
  &  eigen,&  
  &  kg,lavnlr,mcg,mpi_enreg,&  
  &  nattyp,&  
  &  npwarr,nspinor,&  
  &  occ,&  
  &  ph1d,psps,rhor,rprimd,&  
  &  wffnow,xred,ylm)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: mcg
  integer,intent(in) :: nspinor
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(in) :: cg(2,mcg)
  real(dp), intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  real(dp), intent(out) :: lavnlr(dtset%nfft,dtset%nspden)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp), intent(inout) :: rhor(dtset%nfft,dtset%nspden)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(in) :: xred(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine lavnl
end interface

interface
 subroutine moddiel_csrb(dielar,dtset,gprimd,mpi_enreg,rdiemac,rhor_in)
  use defs_basis
  use defs_abitypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: rdiemac(dtset%nfft,dtset%nspden)
  real(dp),intent(in) :: rhor_in(dtset%nfft,dtset%nspden)
 end subroutine moddiel_csrb
end interface

interface
 subroutine newrho(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,dtn_pc,dtset,etotal,fcart,ffttomix,&  
  &  gmet,grhf,gsqcut,initialized,ispmix,istep,kg_diel,kxc,mgfft,mix,mixtofft,&  
  &  moved_atm_inside,mpi_enreg,nattyp,nfft,nfftmix,ngfft,ngfftmix,nkxc,npawmix,npwdiel,&  
  &  nresid,ntypat,n1xccc,pawrhoij,pawtab,&  
  &  ph1d,psps,rhog,rhor,rprimd,susmat,usepaw,vtrial,wvl,xred)
  use defs_basis
  use m_ab6_mixing
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(inout) :: dbl_nnsclo
  integer,intent(in) :: dielstrt
  integer,intent(in) :: initialized
  integer,intent(in) :: ispmix
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n1xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftmix
  integer,intent(in) :: nkxc
  integer,intent(in) :: npawmix
  integer,intent(in) :: npwdiel
  integer,intent(in) :: ntypat
  integer,intent(in) :: usepaw
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(ab6_mixing_object), intent(inout) :: mix
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_internal_type), intent(in) :: wvl
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftmix(18)
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(inout), target :: dtn_pc(3,dtset%natom)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  integer,intent(in) :: ffttomix(nfft*(1-nfftmix/nfft))
  real(dp), intent(inout) :: gmet(3,3)
  real(dp),intent(in) :: grhf(3,dtset%natom)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  integer,intent(in) :: mixtofft(nfftmix*(1-nfftmix/nfft))
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(inout) :: nresid(nfft,dtset%nspden)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(out) :: rhog(2,nfft)
  real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(in), target :: vtrial(nfft,dtset%nspden)
  real(dp), intent(inout), target :: xred(3,dtset%natom)
 end subroutine newrho
end interface

interface
 subroutine newvtr(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,&  
  &  dtn_pc,dtset,efermi,etotal,fcart,ffttomix,&  
  &  gmet,grhf,gsqcut,&  
  &  initialized,ispmix,&  
  &  istep,&  
  &  kg_diel,kxc,mgfft,mix,mixtofft,&  
  &  moved_atm_inside,mpi_enreg,nattyp,nfft,nfftmix,&  
  &  nhat,nhatgr,nhatgrdim,&  
  &  ngfft,ngfftmix,nkxc,npawmix,npwdiel,&  
  &  nstep,ntypat,n1xccc,optres,optxc,&  
  &  pawrhoij,&  
  &  ph1d,&  
  &  psps,rhor,rprimd,susmat,usepaw,&  
  &  vhartr,vnew_mean,vpsp,vresid,&  
  &  vtrial,vxc,xred,&  
  &  atindx1,cg,deltae,&  
  &  dtfil,eeig,eigen,ek,enl,kg,&  
  &  mcg,nfftf,&  
  &  ngfftf,npwarr,n3xccc,occ,optene,&  
  &  pawfgr,pawtab,&  
  &  resid,rhog,&  
  &  usexcnhat,&  
  &  wffnow,wvl,&  
  &  ylm,xccc3d )
  use defs_abitypes
  use m_wffile
  use defs_basis
  use m_ab6_mixing
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(inout) :: dbl_nnsclo
  integer,intent(in) :: dielstrt
  integer,intent(in) :: initialized
  integer,intent(in) :: ispmix
  integer,intent(in) :: istep
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfft
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftf
  integer,intent(in) :: nfftmix
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: npawmix
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nstep
  integer,intent(in) :: ntypat
  integer,intent(in) :: optene
  integer,intent(in) :: optres
  integer,intent(in) :: optxc
  integer,intent(in) :: usepaw
  integer,intent(in) :: usexcnhat
  real(dp),intent(in) :: deltae
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eeig
  real(dp),intent(in) :: efermi
  real(dp),intent(out) :: ek
  real(dp),intent(out) :: enl
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(ab6_mixing_object),intent(inout) :: mix
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_internal_type), intent(in) :: wvl
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: ngfftmix(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(inout), target :: dtn_pc(3,dtset%natom)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  integer,intent(in) :: ffttomix(nfft*(1-nfftmix/nfft))
  real(dp),intent(inout) :: gmet(3,3)
  real(dp),intent(in) :: grhf(3,dtset%natom)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  integer,intent(in) :: mixtofft(nfftmix*(1-nfftmix/nfft))
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfftf,dtset%nspden*usepaw)
  real(dp),intent(in) :: nhatgr(nfftf,dtset%nspden,3*nhatgrdim)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(inout) :: rhog(2,nfftf)
  real(dp),intent(inout), target :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(in) :: vnew_mean(dtset%nspden)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(inout) :: vresid(nfft,dtset%nspden)
  real(dp),intent(inout) :: vtrial(nfft,dtset%nspden)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout), target :: xred(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine newvtr
end interface

interface
 subroutine prcref(atindx,dielar,dielinv,&  
  &  dielstrt,dtn_pc,dtset,etotal,fcart,ffttomix,gmet,gsqcut,&  
  &  istep,kg_diel,kxc,&  
  &  mgfft,moved_atm_inside,mpi_enreg,&  
  &  nattyp,nfft,nfftprc,ngfft,ngfftprc,nkxc,npawmix,npwdiel,ntypat,n1xccc,&  
  &  optreal,optres,pawrhoij,pawtab,ph1d,psps,rhog,rhoijrespc,rhor,rprimd,&  
  &  susmat,vhartr,vpsp,vresid,vrespc,vxc,wvl,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: dielstrt
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n1xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftprc
  integer,intent(in) :: nkxc
  integer,intent(in) :: npawmix
  integer,intent(in) :: npwdiel
  integer,intent(in) :: ntypat
  integer,intent(in) :: optreal
  integer,intent(in) :: optres
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_internal_type), intent(in) :: wvl
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftprc(18)
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(out) :: dtn_pc(3,dtset%natom)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  integer,intent(in) :: ffttomix(nfft*(1-nfftprc/nfft))
  real(dp),intent(inout) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  integer,intent(in) :: nattyp(ntypat)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: rhoijrespc(npawmix)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(in) :: vresid(nfftprc*optreal,dtset%nspden)
  real(dp),intent(out) :: vrespc(nfftprc*optreal,dtset%nspden)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine prcref
end interface

interface
 subroutine prcref_PMA(atindx,dielar,dielinv,&  
  &  dielstrt,dtn_pc,dtset,fcart,ffttomix,gmet,gsqcut,&  
  &  istep,kg_diel,kxc,lavnlr,&  
  &  mgfft,moved_atm_inside,mpi_enreg,&  
  &  nattyp,nfft,nfftprc,ngfft,ngfftprc,nkxc,npawmix,npwdiel,ntypat,n1xccc,&  
  &  optreal,optres,pawrhoij,ph1d,psps,rhog,rhoijrespc,rhor,rprimd,&  
  &  susmat,vhartr,vpsp,vresid,vrespc,vxc,xred,&  
  &  deltae,efermi,etotal,nfftf,nhat,nhatgr,nhatgrdim,optene,optxc,&  
  &  pawtab,use_lavnlr,usexcnhat,vtrial,wvl)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: dielstrt
  integer,intent(in) :: istep
  integer,intent(in) :: mgfft
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n1xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftf
  integer,intent(in) :: nfftprc
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: npawmix
  integer,intent(in) :: npwdiel
  integer,intent(in) :: ntypat
  integer,intent(in) :: optene
  integer,intent(in) :: optreal
  integer,intent(in) :: optres
  integer,intent(in) :: optxc
  integer,intent(in) :: use_lavnlr
  integer,intent(in) :: usexcnhat
  real(dp),intent(in) :: deltae
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: efermi
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_internal_type), intent(in) :: wvl
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftprc(18)
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(out) :: dtn_pc(3,dtset%natom)
  real(dp),intent(in) :: fcart(3,dtset%natom)
  integer,intent(in) :: ffttomix(nfft*(1-nfftprc/nfft))
  real(dp),intent(inout) :: gmet(3,3)
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(inout) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: lavnlr(dtset%nfft,dtset%nspden*use_lavnlr)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(nfftf,dtset%nspden,3*nhatgrdim)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: rhoijrespc(npawmix)
  real(dp),intent(in) :: rhor(nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(inout) :: vpsp(nfft)
  real(dp),intent(in) :: vresid(nfftprc*optreal,dtset%nspden)
  real(dp),intent(out) :: vrespc(nfftprc*optreal,dtset%nspden)
  real(dp),intent(in) :: vtrial(dtset%nfft,dtset%nspden)
  real(dp),intent(in) :: vxc(nfft,dtset%nspden)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine prcref_PMA
end interface

interface
 subroutine prcrskerker1(dtset,mpi_enreg,nfft,nspden,ngfft,dielar,etotal,gprimd,vresid,vrespc,base)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  type(dataset_type),intent(in) :: dtset
  real(dp) :: etotal
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: base(nfft)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: vresid(nfft,nspden)
  real(dp),intent(out) :: vrespc(nfft,nspden)
 end subroutine prcrskerker1
end interface

interface
 subroutine prcrskerker2(dtset,nfft,nspden,ngfft,dielar,gprimd,rprimd,vresid,vrespc,natom,xred,mpi_enreg,ucvol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: dielar(7)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vresid(nfft,nspden)
  real(dp),intent(out) :: vrespc(nfft,nspden)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine prcrskerker2
end interface

interface
 subroutine prctfvw1(atindx,atindx1,cg,deltae,dtfil,dtset,eeig,&  
  &  efermi,eigen,ek,enl,etotal,fixmom,gsqcut,&  
  &  kg,mband,mcg,mgfft,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nfft,nfftf,ngfftf,&  
  &  nhat,nhatgr,nhatgrdim,&  
  &  nkpt,nkxc,npwarr,nspden,nspinor,nsppol,ntypat,n3xccc,occ,occopt,optene,optxc,&  
  &  pawfgr,&  
  &  ph1d,psps,resid,rhog,rhor,rprimd,&  
  &  usexcnhat,&  
  &  vin_old,vout_unmixed,vpsp,vtrial,&  
  &  wffnow,wvl,xccc3d,xred,ylm)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  use defs_wvltypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftf
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: optene
  integer,intent(in) :: optxc
  integer,intent(in) :: usexcnhat
  real(dp), intent(in) :: deltae
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eeig
  real(dp), intent(in) :: efermi
  real(dp),intent(out) :: ek
  real(dp),intent(out) :: enl
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: fixmom
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_internal_type), intent(in) :: wvl
  integer, intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp), intent(in) :: cg(2,mcg)
  real(dp), intent(in) :: eigen(mband*nkpt*nsppol)
  integer, intent(in) :: kg(3,mpw*mkmem)
  integer, intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfftf,nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(nfftf,nspden,3*nhatgrdim)
  integer, intent(in) :: npwarr(nkpt)
  real(dp), intent(in) :: occ(mband*nkpt*nsppol)
  real(dp), intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp), intent(out) :: resid(mband*nkpt*nsppol)
  real(dp), intent(inout) :: rhog(2,nfftf)
  real(dp), intent(inout) :: rhor(nfftf,nspden)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(inout) :: vin_old(nfftf,nspden)
  real(dp), intent(inout) :: vout_unmixed(nfftf,nspden)
  real(dp), intent(inout) :: vpsp(nfftf)
  real(dp), intent(inout) :: vtrial(nfftf,nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp), intent(in) :: xred(3,natom)
  real(dp), intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 end subroutine prctfvw1
end interface

interface
 subroutine prctfvw2(atindx,atindx1,cg,dtfil,dtset,eeig,&  
  &  efermi,eigen,ek,enl,fixmom,gsqcut,&  
  &  kg,mband,mcg,mgfft,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nfft,nfftf,ngfftf,&  
  &  nhat,nhatgr,nhatgrdim,nkpt,nkxc,&  
  &  npwarr,nspden,nspinor,nsppol,ntypat,n3xccc,occ,occopt,optene,optres,optxc,&  
  &  pawfgr,&  
  &  ph1d,psps,resid,rhog,rhor,rprimd,&  
  &  usexcnhat,vin_old,vpsp,vtrial,&  
  &  wffnow,wvl,xccc3d,xred,ylm)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  use defs_wvltypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftf
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: occopt
  integer,intent(in) :: optene
  integer,intent(in) :: optres
  integer,intent(in) :: optxc
  integer,intent(in) :: usexcnhat
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: eeig
  real(dp), intent(in) :: efermi
  real(dp),intent(out) :: ek
  real(dp),intent(out) :: enl
  real(dp),intent(in) :: fixmom
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_internal_type), intent(in) :: wvl
  integer, intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx(natom)
  integer,intent(in) :: atindx1(natom)
  real(dp), intent(in) :: cg(2,mcg)
  real(dp), intent(in) :: eigen(mband*nkpt*nsppol)
  integer, intent(in) :: kg(3,mpw*mkmem)
  integer, intent(in) :: nattyp(ntypat)
  real(dp),intent(in) :: nhat(nfft,nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(nfft,nspden,3*nhatgrdim)
  integer, intent(in) :: npwarr(nkpt)
  real(dp), intent(in) :: occ(mband*nkpt*nsppol)
  real(dp), intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
  real(dp), intent(out) :: resid(mband*nkpt*nsppol)
  real(dp), intent(inout) :: rhog(2,nfftf)
  real(dp), intent(inout) :: rhor(nfftf,nspden)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(inout) :: vin_old(nfftf,nspden)
  real(dp), intent(inout) :: vpsp(nfftf)
  real(dp), intent(inout) :: vtrial(nfftf,nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp), intent(in) :: xred(3,natom)
  real(dp), intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 end subroutine prctfvw2
end interface

interface
 subroutine prctfw3(deltae,dtset,&  
  &  efermi,etotal,gsqcut,&  
  &  lavnlr,mpi_enreg,&  
  &  nhat,nhatgr,nhatgrdim,&  
  &  nkxc,n3xccc,optene,optxc,&  
  &  psps,rhor_in,rprimd,&  
  &  usexcnhat,&  
  &  vpsp,vresid,vrespc,vtrial,wvl,&  
  &  xccc3d,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nhatgrdim
  integer,intent(in) :: nkxc
  integer,intent(in) :: optene
  integer,intent(in) :: optxc
  integer,intent(in) :: usexcnhat
  real(dp), intent(in) :: deltae
  type(dataset_type),intent(in) :: dtset
  real(dp), intent(in) :: efermi
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_internal_type), intent(in) :: wvl
  real(dp), intent(in) :: lavnlr(dtset%nfft,dtset%nspden)
  real(dp),intent(in) :: nhat(dtset%nfft,dtset%nspden*psps%usepaw)
  real(dp),intent(in) :: nhatgr(dtset%nfft,dtset%nspden,3*nhatgrdim)
  real(dp), intent(in) :: rhor_in(dtset%nfft,dtset%nspden)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(inout) :: vpsp(dtset%nfft)
  real(dp), intent(in) :: vresid(dtset%nfft,dtset%nspden)
  real(dp), intent(out) :: vrespc(dtset%nfft,dtset%nspden)
  real(dp), intent(in) :: vtrial(dtset%nfft,dtset%nspden)
  real(dp),dimension(:),intent(inout) :: xccc3d(n3xccc)
  real(dp), intent(in) :: xred(3,dtset%natom)
 end subroutine prctfw3
end interface

end module interfaces_68_rsprc
!!***
