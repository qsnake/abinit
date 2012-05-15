!!****m* ABINIT/interfaces_77_ddb
!! NAME
!! interfaces_77_ddb
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/77_ddb
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

module interfaces_77_ddb

 implicit none

interface
 subroutine alignph(amu,displ,d2cart,mpert,natom,ntypat,phfrq,typat)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(inout) :: displ(2,3*natom,3*natom)
  real(dp),intent(in) :: phfrq(3*natom)
  integer,intent(in) :: typat(natom)
 end subroutine alignph
end interface

interface
 subroutine anaddb_dtset_clean(anaddb_dtset)
  use defs_abitypes
  implicit none
  type(anaddb_dataset_type), intent(inout) :: anaddb_dtset
 end subroutine anaddb_dtset_clean
end interface

interface
 subroutine anaddb_dtset_nullify(anaddb_dtset)
  use defs_abitypes
  implicit none
  type(anaddb_dataset_type), intent(inout) :: anaddb_dtset
 end subroutine anaddb_dtset_nullify
end interface

interface
 subroutine asrif9(asr,atmfrc,natom,nrpt,rpt,wghatm)
  use defs_basis
  implicit none
  integer,intent(in) :: asr
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine asrif9
end interface

interface
 subroutine axial9(ifccar,vect1,vect2,vect3)
  use defs_basis
  implicit none
  real(dp),intent(in) :: ifccar(3,3)
  real(dp),intent(in) :: vect1(3)
  real(dp),intent(out) :: vect2(3)
  real(dp),intent(out) :: vect3(3)
 end subroutine axial9
end interface

interface
 subroutine bigbx9(brav,choice,mrpt,ngqpt,nqshft,nrpt,rprim,rpt)
  use defs_basis
  implicit none
  integer,intent(in) :: brav
  integer,intent(in) :: choice
  integer,intent(in) :: mrpt
  integer,intent(in) :: nqshft
  integer,intent(out) :: nrpt
  integer,intent(in) :: ngqpt(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(out) :: rpt(3,mrpt)
 end subroutine bigbx9
end interface

interface
 subroutine canat9(brav,natom,rcan,rprim,trans,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: brav
  integer,intent(in) :: natom
  real(dp),intent(out) :: rcan(3,natom)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(out) :: trans(3,natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine canat9
end interface

interface
 subroutine canct9(acell,gprim,ib,index,irpt,natom,nrpt,&  
  &  rcan,rcart,rprim,rpt)
  use defs_basis
  implicit none
  integer,intent(out) :: ib
  integer,intent(in) :: index
  integer,intent(out) :: irpt
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(out) :: rcart(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
 end subroutine canct9
end interface

interface
 subroutine carteig2d(blkflg,blkval,carflg,d2cart,&  
  &  gprimd,iblok,mpert,natom,nblok,rprimd)
  use defs_basis
  implicit none
  integer,intent(in) :: iblok
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  integer,intent(in) :: blkflg(3,mpert,3,mpert,nblok)
  real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
  integer,intent(out) :: carflg(3,mpert,3,mpert)
  real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine carteig2d
end interface

interface
 subroutine carttransf(blkflg,blkval2,carflg,gprimd,iqpt,mband,&  
  &  mpert,msize,natom,nblok,nkpt,rprimd)
  use defs_basis
  implicit none
  integer,intent(in) :: iqpt
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(inout) :: natom
  integer,intent(in) :: nblok
  integer,intent(inout) :: nkpt
  integer,intent(in) :: blkflg(msize,nblok)
  real(dp),intent(inout) :: blkval2(2,msize,mband,nkpt)
  integer,intent(out) :: carflg(3,mpert,3,mpert)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine carttransf
end interface

interface
 subroutine chki8(inti,intt,name)
  implicit none
  integer,intent(in) :: inti
  integer,intent(in) :: intt
  character(len=6),intent(in) :: name
 end subroutine chki8
end interface

interface
 subroutine chkin9(atifc,natifc,natom)
  implicit none
  integer,intent(in) :: natifc
  integer,intent(in) :: natom
  integer,intent(inout) :: atifc(natom)
 end subroutine chkin9
end interface

interface
 subroutine chkr8(reali,realt,name,tol)
  use defs_basis
  implicit none
  character(len=6),intent(in) :: name
  real(dp),intent(in) :: reali
  real(dp),intent(in) :: realt
  real(dp),intent(in) :: tol
 end subroutine chkr8
end interface

interface
 subroutine chkrp9(brav,rprim)
  use defs_basis
  implicit none
  integer,intent(in) :: brav
  real(dp),intent(in) :: rprim(3,3)
 end subroutine chkrp9
end interface

interface
 subroutine cmpar8 (acell,acell8,amu,amu8,dimekb,ecut,ecut8,ekb,ekb8,&  
  &  fullinit,fullinit8,iscf,iscf8,ixc,ixc8,kpt,kpt8,kptnrm,kptnr8,&  
  &  natom,natom8,nband,nband8,ngfft,ngfft8,nkpt,nkpt8,&  
  &  nsppol,nsppo8,nsym,nsym8,ntypat,ntypat8,occ,occ8,&  
  &  occopt,occop8,pawecutdg,pawecutdg8,pawtab,pawtab8,&  
  &  rprim,rprim8,sciss,sciss8,symrel,symre8,&  
  &  tnons,tnons8,tolwfr,tolwf8,typat,typat8,usepaw,wtk,wtk8,xred,xred8,zion,zion8)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: dimekb
  integer,intent(inout) :: fullinit
  integer,intent(in) :: fullinit8
  integer,intent(inout) :: iscf
  integer,intent(in) :: iscf8
  integer,intent(inout) :: ixc
  integer,intent(in) :: ixc8
  integer,intent(inout) :: natom
  integer,intent(in) :: natom8
  integer,intent(inout) :: nkpt
  integer,intent(in) :: nkpt8
  integer,intent(in) :: nsppo8
  integer,intent(inout) :: nsppol
  integer,intent(inout) :: nsym
  integer,intent(in) :: nsym8
  integer,intent(inout) :: ntypat
  integer,intent(in) :: ntypat8
  integer,intent(in) :: occop8
  integer,intent(inout) :: occopt
  integer,intent(in) :: usepaw
  real(dp),intent(inout) :: ecut
  real(dp),intent(in) :: ecut8
  real(dp),intent(in) :: kptnr8
  real(dp),intent(inout) :: kptnrm
  real(dp),intent(inout) :: pawecutdg
  real(dp),intent(in) :: pawecutdg8
  real(dp),intent(inout) :: sciss
  real(dp),intent(in) :: sciss8
  real(dp),intent(in) :: tolwf8
  real(dp),intent(inout) :: tolwfr
  integer,intent(inout) :: nband(*)
  integer,intent(in) :: nband8(*)
  integer,intent(inout) :: ngfft(18)
  integer,intent(in) :: ngfft8(18)
  integer,intent(in) :: symre8(3,3,*)
  integer,intent(inout) :: symrel(3,3,*)
  integer,intent(inout) :: typat(*)
  integer,intent(in) :: typat8(*)
  real(dp),intent(inout) :: acell(3)
  real(dp),intent(in) :: acell8(3)
  real(dp),intent(inout) :: amu(*)
  real(dp),intent(in) :: amu8(*)
  real(dp),intent(inout) :: ekb(dimekb,*)
  real(dp),intent(in) :: ekb8(dimekb,*)
  real(dp),intent(inout) :: kpt(3,*)
  real(dp),intent(in) :: kpt8(3,*)
  real(dp),intent(inout) :: occ(*)
  real(dp),intent(in) :: occ8(*)
  type(pawtab_type),intent(inout) :: pawtab(*)
  type(pawtab_type),intent(in) :: pawtab8(*)
  real(dp),intent(inout) :: rprim(3,3)
  real(dp),intent(in) :: rprim8(3,3)
  real(dp),intent(inout) :: tnons(3,*)
  real(dp),intent(in) :: tnons8(3,*)
  real(dp),intent(inout) :: wtk(*)
  real(dp),intent(in) :: wtk8(*)
  real(dp),intent(inout) :: xred(3,*)
  real(dp),intent(in) :: xred8(3,*)
  real(dp),intent(inout) :: zion(*)
  real(dp),intent(in) :: zion8(*)
 end subroutine cmpar8
end interface

interface
 subroutine complete_gamma(Cryst,nbranch,nsppol,nqptirred,nqpt_full,ep_scalprod,qirredtofull,qpttoqpt,gamma_qpt)
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: ep_scalprod
  integer,intent(in) :: nbranch
  integer,intent(in) :: nqpt_full
  integer,intent(in) :: nqptirred
  integer,intent(in) :: nsppol
  type(crystal_structure),intent(in) :: Cryst
  real(dp), intent(inout) :: gamma_qpt(2,nbranch**2,nsppol,nqpt_full)
  integer,intent(in) :: qirredtofull(nqptirred)
  integer,intent(in) :: qpttoqpt(2,Cryst%nsym,nqpt_full)
 end subroutine complete_gamma
end interface

interface
 subroutine complete_gamma_tr(elph_ds,gamma_qpt_tr,&  
  &  gprimd,indsym,natom,nsym,qpttoqpt,rprimd,&  
  &  symrec,symrel)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  type(elph_type),intent(inout) :: elph_ds
  real(dp), intent(inout) :: gamma_qpt_tr(2,9,elph_ds%nbranch*elph_ds%nbranch, &
  &         elph_ds%nsppol,elph_ds%nqpt_full)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: qpttoqpt(2,nsym,elph_ds%nqpt_full)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine complete_gamma_tr
end interface

interface
 subroutine complete_gkk(elph_ds,gkk_flag,&  
  &  gprimd,indsym,natom,nsym,qpttoqpt,rprimd,&  
  &  symrec,symrel)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(inout) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch, &
  &         elph_ds%k_phon%nkpt,elph_ds%nsppol,elph_ds%nqpt_full)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: qpttoqpt(2,nsym,elph_ds%nqpt_full)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine complete_gkk
end interface

interface
 subroutine completeperts(Cryst,nbranch,nFSband,nkpt,nsppol,gkk_flag,h1_mat_el,h1_mat_el_sq,&  
  &  qpt,symq,qtimrev)
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: nFSband
  integer,intent(in) :: nbranch
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: qtimrev
  type(crystal_structure),intent(in) :: Cryst
  integer,intent(inout) :: gkk_flag(nbranch,nbranch,nkpt,nsppol)
  real(dp),intent(in) :: h1_mat_el(2,nFSband**2,nbranch,nkpt,nsppol)
  real(dp),intent(out) :: h1_mat_el_sq(2,nFSband**2,nbranch**2,nkpt,nsppol)
  real(dp),intent(in) :: qpt(3)
  integer,intent(in) :: symq(4,2,Cryst%nsym)
 end subroutine completeperts
end interface

interface
 subroutine diel9(amu,anaddb_dtset,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&  
  &  iout,lst,mpert,natom,nph2l,ntypat,phfrq,qtol,typat,ucvol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nph2l
  integer,intent(in) :: ntypat
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  real(dp),intent(in) :: qtol
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(out) :: dielt_rlx(3,3)
  real(dp),intent(inout) :: displ(2,3*natom,3*natom)
  real(dp),intent(out) :: epsinf(3,3)
  real(dp),intent(out) :: fact_oscstr(2,3,3*natom)
  real(dp),intent(in) :: lst(nph2l)
  real(dp),intent(in) :: phfrq(3*natom)
  integer,intent(in) :: typat(natom)
 end subroutine diel9
end interface

interface
 subroutine dist9(acell,dist,gprim,natom,nrpt,rcan,rprim,rpt)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  real(dp),intent(in) :: acell(3)
  real(dp),intent(out) :: dist(natom,natom,nrpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
 end subroutine dist9
end interface

interface
 subroutine dtchi(blkval,dchide,dchidt,mpert,natom,ramansr)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: ramansr
  real(dp),intent(in) :: blkval(2,3*mpert*3*mpert*3*mpert)
  real(dp),intent(out) :: dchide(3,3,3)
  real(dp),intent(out) :: dchidt(natom,3,3,3)
 end subroutine dtchi
end interface

interface
 subroutine dtech9(blkval,dielt,iblok,mpert,natom,nblok,zeff)
  use defs_basis
  implicit none
  integer,intent(in) :: iblok
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
  real(dp),intent(out) :: dielt(3,3)
  real(dp),intent(out) :: zeff(3,3,natom)
 end subroutine dtech9
end interface

interface
 subroutine dymfz9(dynmat,natom,nqpt,gprim,option,&  
  &  spqpt,trans)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: option
  real(dp),intent(inout) :: dynmat(2,3,natom,3,natom,nqpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: spqpt(3,nqpt)
  real(dp),intent(in) :: trans(3,natom)
 end subroutine dymfz9
end interface

interface
 subroutine elast9(anaddb_dtset,blkval,d2asr,elast,iblok,iblok_stress,instrain,iout,mpert,&  
  &  natom,nblok,ucvol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iblok
  integer,intent(in) :: iblok_stress
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
  real(dp),intent(inout) :: d2asr(2,3,natom,3,natom)
  real(dp),intent(out) :: elast(6,6)
  real(dp),intent(in) :: instrain(3*natom,6)
 end subroutine elast9
end interface

interface
 subroutine electrooptic(dchide,dieflag,&  
  &  epsinf,fact_oscstr,natom,phfrq,prtmbm,rsus,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: dieflag
  integer,intent(in) :: natom
  integer,intent(in) :: prtmbm
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: dchide(3,3,3)
  real(dp),intent(in) :: epsinf(3,3)
  real(dp),intent(in) :: fact_oscstr(2,3,3*natom)
  real(dp),intent(in) :: phfrq(3*natom)
  real(dp),intent(in) :: rsus(3*natom,3,3)
 end subroutine electrooptic
end interface

interface
 subroutine eli_app_m_1d (delta_1d,lambda_1d,nmatsu,z_1d)
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: z_1d(-nmatsu:nmatsu)
 end subroutine eli_app_m_1d
end interface

interface
 subroutine eli_diag_m_1d (delta_1d,lambda_1d,maxeigval,mustar,nmatsu,tc,z_1d)
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  real(dp),intent(out) :: maxeigval
  real(dp),intent(in) :: mustar
  real(dp),intent(in) :: tc
  real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: z_1d(-nmatsu:nmatsu)
 end subroutine eli_diag_m_1d
end interface

interface
 subroutine eli_lambda_1d (a2f_1d,elph_ds,lambda_1d,nmatsu,tc)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  type(elph_type),intent(in) :: elph_ds
  real(dp),intent(in) :: tc
  real(dp),intent(in) :: a2f_1d(elph_ds%na2f)
  real(dp),intent(out) :: lambda_1d(-nmatsu:nmatsu)
 end subroutine eli_lambda_1d
end interface

interface
 subroutine eli_m_iter_1d (delta_1d,lambda_1d,maxeigval,nmatsu,z_1d)
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  real(dp),intent(out) :: maxeigval
  real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: z_1d(-nmatsu:nmatsu)
 end subroutine eli_m_iter_1d
end interface

interface
 subroutine eli_z_1d (lambda_1d,nmatsu,z_1d)
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
  real(dp),intent(out) :: z_1d(-nmatsu:nmatsu)
 end subroutine eli_z_1d
end interface

interface
 subroutine eliashberg_1d(a2f_1d,elph_ds,mustar)
  use defs_elphon
  use defs_basis
  implicit none
  type(elph_type),intent(in) :: elph_ds
  real(dp),intent(in) :: mustar
  real(dp),intent(in) :: a2f_1d(elph_ds%na2f)
 end subroutine eliashberg_1d
end interface

interface
 subroutine elphon(anaddb_dtset,filnam,acell_in,amu,atmfrc,dielt,dyewq0,gmet,&  
  &  gprim,indsym,mpert,mpi_enreg,natom,nrpt,nsym,ntypat,rcan,rmet,rprim_in,rpt,&  
  &  symrec,symrel,tnons,trans,typat,ucvol,wghatm,xred,zeff)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  type(anaddb_dataset_type) :: anaddb_dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell_in(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(inout) :: dyewq0(3,3,natom)
  character(len=fnlen),intent(in) :: filnam(7)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim_in(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(in) :: trans(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine elphon
end interface

interface
 subroutine ep_fs_weights(ep_b_min, ep_b_max, eigenGS, elphsmear, fermie, gprimd,&  
  &  irredtoGS, kptrlatt, max_occ, minFSband, nband, nFSband, nsppol, telphint, k_obj)
  use defs_elphon
  use defs_basis
  implicit none
  integer, intent(in) :: ep_b_max
  integer, intent(in) :: ep_b_min
  integer, intent(in) :: minFSband
  integer, intent(in) :: nFSband
  integer,intent(in) :: nband
  integer, intent(in) :: nsppol
  integer, intent(in) :: telphint
  real(dp), intent(in) :: elphsmear
  real(dp), intent(in) :: fermie
  type(elph_kgrid_type), intent(inout) :: k_obj
  real(dp), intent(in) :: max_occ
  integer, intent(in) :: kptrlatt(3,3)
  real(dp), intent(in) :: eigenGS(nband,k_obj%nkpt,nsppol)
  real(dp), intent(in) :: gprimd(3,3)
  integer, intent(in) :: irredtoGS(k_obj%nkptirr)
 end subroutine ep_fs_weights
end interface

interface
 subroutine ep_setupqpt (anaddb_dtset,elph_ds,gmet,nsym,qptrlatt,rprimd,symrec,symrel,timrev)
  use defs_elphon
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: nsym
  integer, intent(in) :: timrev
  type(anaddb_dataset_type), intent(in) :: anaddb_dtset
  type(elph_type), intent(inout) :: elph_ds
  integer, intent(out) :: qptrlatt(3,3)
  real(dp), intent(in) :: gmet(3,3)
  real(dp), intent(in) :: rprimd(3,3)
  integer, intent(in) :: symrec(3,3,nsym)
  integer, intent(in) :: symrel(3,3,nsym)
 end subroutine ep_setupqpt
end interface

interface
 subroutine ewald9(acell,dielt,dyew,gmet,gprim,natom,&  
  &  qphon,rmet,rprim,sumg0,ucvol,xred,zeff)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: sumg0
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(out) :: dyew(2,3,natom,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine ewald9
end interface

interface
 subroutine freeze_displ_allmodes(displ, freeze_displ, natom, outfile_radix, phfreq,&  
  &  qphon, rprimd, typat, xcart)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  real(dp), intent(in) :: freeze_displ
  character(len=fnlen),intent(in) :: outfile_radix
  real(dp),intent(in) :: displ(2,3*natom,3*natom)
  real(dp),intent(in) :: phfreq(3*natom)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine freeze_displ_allmodes
end interface

interface
 subroutine ftgam (wghatm,gam_qpt,gam_rpt,gprim,natom,nqpt,nrpt,qtor,rpt,qpt_full)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: nrpt
  integer,intent(in) :: qtor
  real(dp),intent(inout) :: gam_qpt(2,3*natom*3*natom,nqpt)
  real(dp),intent(inout) :: gam_rpt(2,3*natom*3*natom,nrpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qpt_full(3,nqpt)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine ftgam
end interface

interface
 subroutine ftgkk (wghatm,gkk_qpt,gkk_rpt,gkqwrite,gkrwrite,gprim,ikpt_phon0,&  
  &  natom,nkpt_phon,ngkkband,nkpt_used,nqpt,nrpt,nsppol,&  
  &  qtor,rpt,qpt_full,unit_gkk_rpt,unitgkq)
  use defs_basis
  implicit none
  integer,intent(in) :: gkqwrite
  integer,intent(in) :: gkrwrite
  integer,intent(in) :: ikpt_phon0
  integer,intent(in) :: natom
  integer,intent(in) :: ngkkband
  integer,intent(in) :: nkpt_phon
  integer,intent(in) :: nkpt_used
  integer,intent(in) :: nqpt
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: qtor
  integer,intent(in) :: unit_gkk_rpt
  integer,intent(in) :: unitgkq
  real(dp),intent(inout) :: gkk_qpt(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_used,nsppol,nqpt)
  real(dp),intent(inout) :: gkk_rpt(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_used,nsppol,nrpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qpt_full(3,nqpt)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine ftgkk
end interface

interface
 subroutine ftifc_q2r(atmfrc,dynmat,gprim,natom,nqpt,&  
  &  nrpt,rpt,spqpt)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: nrpt
  real(dp),intent(out) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(in) :: dynmat(2,3,natom,3,natom,nqpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: spqpt(3,nqpt)
 end subroutine ftifc_q2r
end interface

interface
 subroutine ftifc_r2q(atmfrc,dynmat,gprim,natom,nqpt,&  
  &  nrpt,rpt,spqpt,wghatm)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: nrpt
  real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(out) :: dynmat(2,3,natom,3,natom,nqpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: spqpt(3,nqpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine ftifc_r2q
end interface

interface
 subroutine fxgkkphase(elph_ds,gkk_flag,h1_mat_el,iqptfull)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptfull
  type(elph_type),intent(in) :: elph_ds
  integer,intent(in) :: gkk_flag(elph_ds%nbranch,elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  real(dp),intent(inout) :: h1_mat_el(2,elph_ds%nFSband,elph_ds%nFSband, &
  &         elph_ds%nbranch,elph_ds%k_phon%nkpt)
 end subroutine fxgkkphase
end interface

interface
 subroutine gam_mult_displ(nbranch, displ_red, gam_bare, gam_now)
  use defs_basis
  implicit none
  integer, intent(in) :: nbranch
  real(dp), intent(in) :: displ_red(2,nbranch,nbranch)
  real(dp), intent(in) :: gam_bare(2,nbranch,nbranch)
  real(dp), intent(out) :: gam_now(2,nbranch,nbranch)
 end subroutine gam_mult_displ
end interface

interface
 subroutine gamma9(gamma,qphon,qphnrm,qtol)
  use defs_basis
  implicit none
  integer,intent(out) :: gamma
  real(dp),intent(in) :: qphnrm
  real(dp),intent(in) :: qtol
  real(dp),intent(in) :: qphon(3)
 end subroutine gamma9
end interface

interface
 subroutine get_all_gkk2(&  
  &  elph_ds,kptirr_phon,kpt_phon,&  
  &  natom,nrpt,&  
  &  phon_ds,rcan,&  
  &  wghatm)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  type(elph_type),intent(inout) :: elph_ds
  type(phon_type),intent(inout) :: phon_ds
  real(dp),intent(in) :: kpt_phon(3,elph_ds%k_phon%nkpt)
  real(dp),intent(in) :: kptirr_phon(3,elph_ds%k_phon%nkptirr)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine get_all_gkk2
end interface

interface
 subroutine get_all_gkq (elph_ds,Cryst,Bst,FSfullpqtofull,nband,n1wf,onegkksize,phon_ds,&  
  &  qpttoqpt,ep_prt_yambo,unitgkk)
  use defs_elphon
  use defs_datatypes
  use m_crystal
  implicit none
  integer,intent(in) :: ep_prt_yambo
  integer,intent(in) :: n1wf
  integer,intent(in) :: nband
  integer,intent(in) :: onegkksize
  integer,intent(in) :: unitgkk
  type(bandstructure_type),intent(in) :: Bst
  type(crystal_structure),intent(in) :: Cryst
  type(elph_type),intent(inout) :: elph_ds
  type(phon_type),intent(inout) :: phon_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  integer,intent(in) :: qpttoqpt(2,Cryst%nsym,elph_ds%nqpt_full)
 end subroutine get_all_gkq
end interface

interface
 subroutine get_all_gkr (elph_ds,gprim,natom,nrpt,onegkksize,rpt,qpt_full,wghatm)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: onegkksize
  type(elph_type),intent(inout) :: elph_ds
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qpt_full(3,elph_ds%nqpt_full)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine get_all_gkr
end interface

interface
 subroutine get_fs_bands(eigenGS,hdr,fermie,ep_b_min,ep_b_max,minFSband,maxFSband,nkptirr)
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: ep_b_max
  integer, intent(in) :: ep_b_min
  integer,intent(out) :: maxFSband
  integer,intent(out) :: minFSband
  integer,intent(out) :: nkptirr
  real(dp),intent(in) :: fermie
  type(hdr_type),intent(in) :: hdr
  real(dp),intent(in) :: eigenGS(hdr%nband(1),hdr%nkpt,hdr%nsppol)
 end subroutine get_fs_bands
end interface

interface
 subroutine get_veloc_tr(elph_ds,mpi_enreg,nband,elph_tr_ds)
  use defs_elphon
  use defs_abitypes
  implicit none
  integer,intent(in) :: nband
  type(elph_type),intent(in) :: elph_ds
  type(elph_tr_type) :: elph_tr_ds
  type(mpi_type), intent(inout) :: mpi_enreg
 end subroutine get_veloc_tr
end interface

interface
 subroutine gtblk9(ddb_blk,iblok,mpert,&  
  &  natom,qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)
  use defs_basis
  use m_ddb_blk
  implicit none
  integer,intent(out) :: iblok
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: rftyp
  type(ddb_blk_type), pointer :: ddb_blk
  real(dp),intent(in) :: qtol
  integer,intent(in) :: rfelfd(4)
  integer,intent(in) :: rfphon(4)
  integer,intent(in) :: rfstrs(4)
  real(dp),intent(inout) :: qphnrm(3)
  real(dp),intent(inout) :: qphon(3,3)
 end subroutine gtblk9
end interface

interface
 subroutine gtdyn9(acell,atmfrc,dielt,dipdip,&  
  &  dyewq0,d2cart,gmet,gprim,mpert,natom,&  
  &  nrpt,qphnrm,qpt,rmet,rprim,rpt,&  
  &  trans,ucvol,wghatm,xred,zeff)
  use defs_basis
  implicit none
  integer,intent(in) :: dipdip
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  real(dp),intent(in) :: qphnrm
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(in) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qpt(3)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: trans(3,natom)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine gtdyn9
end interface

interface
 subroutine hybrid9(acell,asr,atmfrc,dielt,dipdip,dyew,dyewq0,&  
  &  gmet,gprim,iout,natom,nrpt,rcan,rmet,&  
  &  rprim,rpt,ucvol,wghatm,xred,zeff)
  use defs_basis
  implicit none
  integer,intent(in) :: asr
  integer,intent(in) :: dipdip
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(inout) :: dielt(3,3)
  real(dp),intent(inout) :: dyew(2,3,natom,3,natom)
  real(dp),intent(inout) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(inout) :: zeff(3,3,natom)
 end subroutine hybrid9
end interface

interface
 subroutine ifclo9(ifccar,ifcloc,vect1,vect2,vect3)
  use defs_basis
  implicit none
  real(dp),intent(in) :: ifccar(3,3)
  real(dp),intent(out) :: ifcloc(3,3)
  real(dp),intent(in) :: vect1(3)
  real(dp),intent(in) :: vect2(3)
  real(dp),intent(in) :: vect3(3)
 end subroutine ifclo9
end interface

interface
 subroutine init8(dscrpt,filnam,mddb,nddb)
  use defs_basis
  implicit none
  integer,intent(in) :: mddb
  integer,intent(out) :: nddb
  character(len=fnlen),intent(out) :: dscrpt
  character(len=fnlen),intent(out) :: filnam(mddb+1)
 end subroutine init8
end interface

interface
 subroutine init9(filnam)
  use defs_basis
  implicit none
  character(len=fnlen),intent(out) :: filnam(7)
 end subroutine init9
end interface

interface
 subroutine inpphon(displ_cart,pheigval,pheigvec,phfrq,phon_ds,qpt)
  use defs_elphon
  use defs_basis
  implicit none
  type(phon_type),intent(inout) :: phon_ds
  real(dp),intent(out) :: displ_cart(2,3*phon_ds%natom,3*phon_ds%natom)
  real(dp),intent(out) :: pheigval(3*phon_ds%natom)
  real(dp),intent(out) :: pheigvec(2*3*phon_ds%natom*3*phon_ds%natom)
  real(dp),intent(out) :: phfrq(3*phon_ds%natom)
  real(dp),intent(inout) :: qpt(3)
 end subroutine inpphon
end interface

interface
 subroutine instr9(asr,blkval,d2asr,iblok,instrain,iout,mpert,natom,nblok)
  use defs_basis
  implicit none
  integer,intent(in) :: asr
  integer,intent(in) :: iblok
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
  real(dp),intent(in) :: d2asr(2,3,natom,3,natom)
  real(dp),intent(out) :: instrain(3*natom,6)
 end subroutine instr9
end interface

interface
 subroutine integrate_gamma(elph_ds,FSfullpqtofull,nrpt)
  use defs_elphon
  implicit none
  integer, intent(in) :: nrpt
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 end subroutine integrate_gamma
end interface

interface
 subroutine integrate_gamma_alt(elph_ds,elph_tr_ds,Cryst,gprim,kptrlatt,&  
  &  natom,nrpt,nsym,qpttoqpt,rpt,wghatm)
  use defs_elphon
  use defs_basis
  use m_crystal
  implicit none
  integer, intent(in) :: natom
  integer, intent(in) :: nrpt
  integer, intent(in) :: nsym
  type(crystal_structure),intent(in) :: Cryst
  type(elph_type),intent(inout) :: elph_ds
  type(elph_tr_type), intent(inout) :: elph_tr_ds
  integer, intent(in) :: kptrlatt(3,3)
  real(dp), intent(in) :: gprim(3,3)
  integer,intent(in) :: qpttoqpt(2,nsym,elph_ds%nqpt_full)
  real(dp), intent(in) :: rpt(3,nrpt)
  real(dp), intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine integrate_gamma_alt
end interface

interface
 subroutine integrate_gamma_tr(elph_ds,FSfullpqtofull,nrpt,elph_tr_ds)
  use defs_elphon
  implicit none
  integer,intent(in) :: nrpt
  type(elph_type),intent(in) :: elph_ds
  type(elph_tr_type), intent(inout) :: elph_tr_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 end subroutine integrate_gamma_tr
end interface

interface
 subroutine interpolate_gkk (&  
  &  elph_ds,kpt_phon,&  
  &  gprim,natom,&  
  &  nrpt,phon_ds,rpt,&  
  &  wghatm)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  type(elph_type),intent(inout) :: elph_ds
  type(phon_type),intent(inout) :: phon_ds
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: kpt_phon(3,elph_ds%k_phon%nkpt)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine interpolate_gkk
end interface

interface
 subroutine interpolate_phfrq (acell,amu,atmfrc,dielt,dipdip,&  
  &  dyewq0,kptirr_phon,kpt_phon,gmet,&  
  &  gprim,indsym,mpert,msym,natom,&  
  &  nbranch,nFSband,nkpt_phon,nkpt_phonirred,&  
  &  nrpt,nsym,ntypat,phfrq,&  
  &  rmet,rprim,rprimd,rpt,&  
  &  symrel,trans,typat,ucvol,&  
  &  wghatm,xred,zeff)
  use defs_basis
  implicit none
  integer,intent(in) :: dipdip
  integer,intent(in) :: mpert
  integer,intent(in) :: msym
  integer,intent(in) :: nFSband
  integer,intent(in) :: natom
  integer,intent(in) :: nbranch
  integer,intent(in) :: nkpt_phon
  integer,intent(in) :: nkpt_phonirred
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(in) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: kpt_phon(3,nkpt_phon)
  real(dp),intent(in) :: kptirr_phon(3,nkpt_phonirred)
  real(dp),intent(out) :: phfrq(3*natom,nkpt_phon,nkpt_phonirred)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: trans(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine interpolate_phfrq
end interface

interface
 subroutine invars9 (anaddb_dtset,lenstr,natom,qtol,string)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: lenstr
  integer,intent(in) :: natom
  type(anaddb_dataset_type),intent(out) :: anaddb_dtset
  real(dp),intent(in) :: qtol
  character(len=*),intent(in) :: string
 end subroutine invars9
end interface

interface
 subroutine k_neighbors (kpt, kptrlatt,kptrank_t,&  
  &  rel_kpt, kpt_phon_indices)
  use defs_basis
  use m_kptrank
  implicit none
  type(kptrank_type), intent(in) :: kptrank_t
  integer, intent(out) :: kpt_phon_indices(8)
  integer, intent(in) :: kptrlatt(3,3)
  real(dp), intent(in) :: kpt(3)
  real(dp), intent(out) :: rel_kpt(3)
 end subroutine k_neighbors
end interface

interface
 subroutine lin_interpq_gam(gamma_qpt,nbranch,nqbz,nsppol,gam_now,isppol,kptrlatt,qpt)
  use defs_basis
  implicit none
  integer, intent(in) :: isppol
  integer, intent(in) :: nbranch
  integer, intent(in) :: nqbz
  integer, intent(in) :: nsppol
  integer, intent(in) :: kptrlatt(3,3)
  real(dp), intent(out) :: gam_now(2,nbranch**2)
  real(dp),intent(in) :: gamma_qpt(2,nbranch**2,nsppol,nqbz)
  real(dp), intent(in) :: qpt(3)
 end subroutine lin_interpq_gam
end interface

interface
 subroutine mblktyp1(ddbun,dscrpt,filnam,mddb,msym,nddb,vrsddb)
  use defs_basis
  implicit none
  integer,intent(in) :: ddbun
  integer,intent(in) :: mddb
  integer,intent(out) :: msym
  integer,intent(in) :: nddb
  integer,intent(in) :: vrsddb
  character(len=fnlen),intent(in) :: dscrpt
  character(len=fnlen),intent(in) :: filnam(mddb+1)
 end subroutine mblktyp1
end interface

interface
 subroutine mblktyp5 (ddbun,dscrpt,filnam,mddb,msym,nddb,vrsddb)
  use defs_basis
  implicit none
  integer,intent(in) :: ddbun
  integer,intent(in) :: mddb
  integer,intent(out) :: msym
  integer,intent(in) :: nddb
  integer,intent(in) :: vrsddb
  character(len=fnlen),intent(in) :: dscrpt
  character(len=fnlen),intent(in) :: filnam(mddb+1)
 end subroutine mblktyp5
end interface

interface
 subroutine mk_irredpert(indsym,iqptfull,irredpert,&  
  &  natom,nbranch,nqpt,nsym,qpt,qtimrev,symq,symrel)
  implicit none
  integer,intent(in) :: iqptfull
  integer,intent(in) :: natom
  integer,intent(in) :: nbranch
  integer,intent(in) :: nqpt
  integer,intent(in) :: nsym
  integer,intent(in) :: qtimrev
  integer,intent(in) :: qpt(3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(out) :: irredpert(7,nbranch,nbranch,nqpt)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine mk_irredpert
end interface

interface
 subroutine mka2f(Cryst,alter_int_gam,a2f_1d,dos_phon,elph_ds,gprim,kptrlatt,mustar,nrpt,phon_ds,rpt,wghatm)
  use defs_elphon
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: alter_int_gam
  integer,intent(in) :: nrpt
  type(crystal_structure),intent(in) :: Cryst
  type(elph_type),intent(inout) :: elph_ds
  real(dp),intent(in) :: mustar
  type(phon_type),intent(inout) :: phon_ds
  integer, intent(in) :: kptrlatt(3,3)
  real(dp),intent(out) :: a2f_1d(elph_ds%na2f)
  real(dp),intent(out) :: dos_phon(elph_ds%na2f)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(Cryst%natom,Cryst%natom,nrpt)
 end subroutine mka2f
end interface

interface
 subroutine mka2fQgrid(elph_ds,fname)
  use defs_elphon
  use defs_basis
  implicit none
  type(elph_type),intent(in) :: elph_ds
  character(len=fnlen),intent(in) :: fname
 end subroutine mka2fQgrid
end interface

interface
 subroutine mka2f_tr(alter_int_gam,elph_ds,gprim,gprimd,ucvol,natom,nrpt,&  
  &  ntemper,tempermin,temperinc,phon_ds,rpt,wghatm,elph_tr_ds)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: alter_int_gam
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: ntemper
  type(elph_type),intent(inout) :: elph_ds
  type(elph_tr_type) :: elph_tr_ds
  type(phon_type),intent(inout) :: phon_ds
  real(dp),intent(in) :: temperinc
  real(dp),intent(in) :: tempermin
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine mka2f_tr
end interface

interface
 subroutine mkFSkgrid (elph_k, nsym, symrec, timrev)
  use defs_elphon
  implicit none
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  type(elph_kgrid_type),intent(inout) :: elph_k
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine mkFSkgrid
end interface

interface
 subroutine mkfsqgrid(kpt_phon,FStoqpt,nkpt_phon,nFSqpt,tmpFSqpt)
  use defs_basis
  implicit none
  integer,intent(out) :: nFSqpt
  integer,intent(in) :: nkpt_phon
  integer,intent(out) :: FStoqpt(nkpt_phon,nkpt_phon)
  real(dp),intent(in) :: kpt_phon(3,nkpt_phon)
  real(dp),intent(out) :: tmpFSqpt(3,nkpt_phon*nkpt_phon)
 end subroutine mkfsqgrid
end interface

interface
 subroutine mkifc9(acell,amu,anaddb_dtset,&  
  &  ddb_blk,dielt,dyewq0,gmet,gprim,&  
  &  ifc_obj,indsym,iout,mpert,msym,natom,ngqpt_in,&  
  &  nsym,ntypat,rcan,rmet,rprim,&  
  &  symrec,symrel,tcpui,trans,twalli,typat,&  
  &  ucvol,xred,zeff)
  use defs_basis
  use defs_abitypes
  use m_ifc
  use m_ddb_blk
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  type(ddb_blk_type), intent(in) :: ddb_blk
  type(ifc_type), intent(out) :: ifc_obj
  real(dp),intent(in) :: tcpui
  real(dp),intent(in) :: twalli
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngqpt_in(3)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(inout) :: dielt(3,3)
  real(dp),intent(out) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym*natom)
  real(dp),intent(out) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  integer,intent(in) :: symrec(3,3,msym)
  integer,intent(in) :: symrel(3,3,msym)
  real(dp),intent(out) :: trans(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(inout) :: zeff(3,3,natom)
 end subroutine mkifc9
end interface

interface
 subroutine mkph_linwid(Cryst,alter_int_gam,elph_ds,gprim,kptrlatt,nrpt,nqpath,phon_ds,qpath_vertices,rpt,wghatm)
  use defs_elphon
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: alter_int_gam
  integer,intent(in) :: nqpath
  integer,intent(in) :: nrpt
  type(crystal_structure),intent(in) :: Cryst
  type(elph_type),intent(inout) :: elph_ds
  type(phon_type),intent(inout) :: phon_ds
  integer,intent(in) :: kptrlatt(3,3)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qpath_vertices(3,nqpath)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(Cryst%natom,Cryst%natom,nrpt)
 end subroutine mkph_linwid
end interface

interface
 subroutine mkphbs(acell,amu,anaddb_dtset,atmfrc,ddb_blk,&  
  &  d2asr,dielt,dyewq0,outfile_radix,gmet,gprim,indsym,iodyn,&  
  &  mpert,msize,msym,natom,nrpt,nsym,ntypat,&  
  &  qtol,rmet,rprim,rpt,singular,symrel,tcpui,&  
  &  trans,twalli,typat,ucvol,uinvers,vtinvers,wghatm,xred,zeff)
  use defs_basis
  use defs_abitypes
  use m_ddb_blk
  implicit none
  integer,intent(in) :: iodyn
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  type(ddb_blk_type), pointer :: ddb_blk
  character(len=fnlen),intent(in) :: outfile_radix
  real(dp),intent(in) :: qtol
  real(dp),intent(in) :: tcpui
  real(dp),intent(in) :: twalli
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(inout) :: d2asr(2,3,natom,3,natom)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(inout) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(inout) :: singular(1:3*natom*(3*natom-1)/2)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: trans(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: uinvers(1:3*natom*(3*natom-1)/2,1:3*natom*(3*natom-1)/2)
  real(dp),intent(inout) :: vtinvers(1:3*natom*(3*natom-1)/2,1:3*natom*(3*natom-1)/2)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine mkphbs
end interface

interface
 subroutine mkqptequiv(FSfullpqtofull,Cryst,kpt_phon,nkpt_phon,nqpt,qpttoqpt,qpt_full)
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: nkpt_phon
  integer,intent(in) :: nqpt
  type(crystal_structure),intent(in) :: Cryst
  integer,intent(out) :: FSfullpqtofull(nkpt_phon,nqpt)
  real(dp),intent(in) :: kpt_phon(3,nkpt_phon)
  real(dp),intent(in) :: qpt_full(3,nqpt)
  integer,intent(out) :: qpttoqpt(2,Cryst%nsym,nqpt)
 end subroutine mkqptequiv
end interface

interface
 subroutine nanal9(dyew,dynmat,iqpt,natom,nqpt,plus)
  use defs_basis
  implicit none
  integer,intent(in) :: iqpt
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: plus
  real(dp),intent(in) :: dyew(2,3,natom,3,natom)
  real(dp),intent(out) :: dynmat(2,3,natom,3,natom,nqpt)
 end subroutine nanal9
end interface

interface
 subroutine nmsq_gam (accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&  
  &  h1_mat_el_sq,iqptirred)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptirred
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%nkpt,elph_ds%nsppol)
 end subroutine nmsq_gam
end interface

interface
 subroutine nmsq_gam_sumFS(accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&  
  &  h1_mat_el_sq,iqptirred)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptirred
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%nkpt,elph_ds%nsppol)
 end subroutine nmsq_gam_sumFS
end interface

interface
 subroutine nmsq_pure_gkk(accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,&  
  &  h1_mat_el_sq,iqptirred)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptirred
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%nkpt,elph_ds%nsppol)
 end subroutine nmsq_pure_gkk
end interface

interface
 subroutine nmsq_pure_gkk_sumfs(accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,h1_mat_el_sq,iqptirred)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptirred
  type(elph_type),intent(in) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%nkpt,elph_ds%nsppol)
 end subroutine nmsq_pure_gkk_sumfs
end interface

interface
 subroutine normsq_gkq(displ_red,eigvec,elph_ds,FSfullpqtofull,&  
  &  h1_mat_el_sq,iqptirred,phfrq_tmp,qpt_irred,qdata)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptirred
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%nkpt,elph_ds%nsppol)
  real(dp),intent(in) :: phfrq_tmp(elph_ds%nbranch)
  real(dp),intent(out) :: qdata(elph_ds%nbranch,elph_ds%nsppol,3)
  real(dp),intent(in) :: qpt_irred(3,elph_ds%nqptirred)
 end subroutine normsq_gkq
end interface

interface
 subroutine omega_decomp(amu,natom,ntypat,typat,&  
  &  dynmatfl,dynmatsr,dynmatlr,iqpt,nqpt,eigenvec)
  use defs_basis
  implicit none
  integer,intent(in) :: iqpt
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(inout) :: dynmatfl(2,3,natom,3,natom,nqpt)
  real(dp),intent(inout) :: dynmatlr(2,3,natom,3,natom,nqpt)
  real(dp),intent(inout) :: dynmatsr(2,3,natom,3,natom,nqpt)
  real(dp),intent(in) :: eigenvec(2*3*natom*3*natom)
  integer,intent(in) :: typat(natom)
 end subroutine omega_decomp
end interface

interface
 subroutine order_fs_kpts(kptirr,FSirredtoGS,hdr,nkptirr)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nkptirr
  type(hdr_type),intent(in) :: hdr
  integer,intent(out) :: FSirredtoGS(nkptirr)
  real(dp),intent(out) :: kptirr(3,nkptirr)
 end subroutine order_fs_kpts
end interface

interface
 subroutine outelph(elph_ds,enunit,fname)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: enunit
  type(elph_type),intent(in) :: elph_ds
  character(len=fnlen),intent(in) :: fname
 end subroutine outelph
end interface

interface
 subroutine outg2f(deltaene,enemin,enemax,filnam,g2f,g2fsmear,kpnt,mband,nene,nkpt,nqpt,ntetra,telphint,unit_g2f)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nene
  integer,intent(in) :: nkpt
  integer,intent(in) :: nqpt
  integer,intent(in) :: ntetra
  integer,intent(in) :: telphint
  integer,intent(in) :: unit_g2f
  real(dp) :: deltaene
  real(dp) :: enemax
  real(dp) :: enemin
  character(len=fnlen),intent(in) :: filnam
  real(dp) :: g2fsmear
  real(dp) :: g2f(mband,nkpt,nene)
  real(dp) :: kpnt(3,nkpt,nqpt)
 end subroutine outg2f
end interface

interface
 subroutine outlwf9 (acell,iodyn,msym,natom,nph1l,nsym,ntypat,rprim,symrel,typat,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: iodyn
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nph1l
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
  integer,intent(in) :: symrel(3,3,msym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine outlwf9
end interface

interface
 subroutine outphdos(deltaene,dos_phon,enemin,enemax,filnam,g2fsmear,nene,nqpt,ntetra,telphint,unit_phdos)
  use defs_basis
  implicit none
  integer,intent(in) :: nene
  integer,intent(in) :: nqpt
  integer,intent(in) :: ntetra
  integer,intent(in) :: telphint
  integer,intent(in) :: unit_phdos
  real(dp) :: deltaene
  real(dp) :: enemax
  real(dp) :: enemin
  character(len=fnlen),intent(in) :: filnam
  real(dp) :: g2fsmear
  real(dp) :: dos_phon(nene)
 end subroutine outphdos
end interface

interface
 subroutine outvars9 (anaddb_dtset,nunit)
  use defs_abitypes
  implicit none
  integer,intent(in) :: nunit
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
 end subroutine outvars9
end interface

interface
 subroutine piezo9(anaddb_dtset,blkval,dielt_rlx,elast,iblok,instrain,iout,mpert,&  
  &  natom,nblok,piezo,ucvol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iblok
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
  real(dp),intent(in) :: dielt_rlx(3,3)
  real(dp),intent(in) :: elast(6,6)
  real(dp),intent(in) :: instrain(3*natom,6)
  real(dp),intent(out) :: piezo(6,3)
 end subroutine piezo9
end interface

interface
 subroutine prt_gkk_yambo(displ_cart,displ_red,kpt_phon,h1_mat_el,iqpt,&  
  &  natom,nFSband,nkpt_phon,&  
  &  phfrq,qptn)
  use defs_basis
  implicit none
  integer,intent(in) :: iqpt
  integer,intent(in) :: nFSband
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt_phon
  real(dp),intent(in) :: displ_cart(2,3*natom,3*natom)
  real(dp),intent(in) :: displ_red(2,3*natom,3*natom)
  real(dp),intent(in) :: h1_mat_el(2,nFSband*nFSband,3*natom,nkpt_phon,1)
  real(dp),intent(in) :: kpt_phon(3,nkpt_phon)
  real(dp),intent(in) :: phfrq(3*natom)
  real(dp),intent(in) :: qptn(3)
 end subroutine prt_gkk_yambo
end interface

interface
 subroutine prtvsound(eigvec, gmet, natom, phfrq, qphon, ucvol)
  use defs_basis
  implicit none
  integer, intent(in) :: natom
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: eigvec(2,3*natom,3*natom)
  real(dp), intent(in) :: gmet(3,3)
  real(dp), intent(in) :: phfrq(3*natom)
  real(dp), intent(in) :: qphon(3)
 end subroutine prtvsound
end interface

interface
 subroutine ramansus(d2cart,dchide,dchidt,displ,mpert,&  
  &  natom,phfrq,qphon,qphnrm,rsus,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  real(dp),intent(in) :: qphnrm
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(in) :: dchide(3,3,3)
  real(dp),intent(in) :: dchidt(natom,3,3,3)
  real(dp),intent(in) :: displ(2,3*natom,3*natom)
  real(dp),intent(in) :: phfrq(3*natom)
  real(dp),intent(inout) :: qphon(3)
  real(dp),intent(out) :: rsus(3*natom,3,3)
 end subroutine ramansus
end interface

interface
 subroutine rchkGSheader (hdr,natom,nband,unitgkk)
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(out) :: nband
  integer,intent(in) :: unitgkk
  type(hdr_type),intent(out) :: hdr
 end subroutine rchkGSheader
end interface

interface
 subroutine rdddb9(acell,atifc,amu,ddb_blk,&  
  &  ddbun,dimekb,filnam,gmet,gprim,indsym,iout,&  
  &  lmnmax,mband,mpert,msize,msym,&  
  &  natifc,natom,nkpt,nsym,ntypat,&  
  &  occopt,rmet,rprim,symq,symrec,symrel,&  
  &  tnons,typat,ucvol,usepaw,xcart,xred,zion)
  use defs_basis
  use m_ddb_blk
  implicit none
  integer,intent(in) :: ddbun
  integer,intent(in) :: dimekb
  integer,intent(in) :: iout
  integer,intent(in) :: lmnmax
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(in) :: msym
  integer,intent(in) :: natifc
  integer,intent(inout) :: natom
  integer,intent(inout) :: nkpt
  integer,intent(inout) :: nsym
  integer,intent(inout) :: ntypat
  integer,intent(inout) :: occopt
  integer,intent(inout) :: usepaw
  type(ddb_blk_type), pointer :: ddb_blk
  character(len=fnlen),intent(in) :: filnam
  real(dp),intent(out) :: ucvol
  integer,intent(out) :: symq(4,2,*)
  real(dp),intent(out) :: acell(3)
  real(dp),intent(out) :: amu(ntypat)
  integer,intent(inout) :: atifc(natom)
  real(dp),intent(out) :: gmet(3,3)
  real(dp),intent(out) :: gprim(3,3)
  integer,intent(out) :: indsym(4,nsym,natom)
  real(dp),intent(out) :: rmet(3,3)
  real(dp),intent(out) :: rprim(3,3)
  integer,intent(out) :: symrec(3,3,msym)
  integer,intent(out) :: symrel(3,3,msym)
  real(dp),intent(out) :: tnons(3,msym)
  integer,intent(out) :: typat(natom)
  real(dp),intent(out) :: xcart(3,natom)
  real(dp),intent(out) :: xred(3,natom)
  real(dp),intent(out) :: zion(ntypat)
 end subroutine rdddb9
end interface

interface
 subroutine read_el_veloc(mpi_enreg,nband_in,nkpt_in,kpt_in,nsppol_in,elph_tr_ds)
  use defs_elphon
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: nband_in
  integer, intent(in) :: nkpt_in
  integer, intent(in) :: nsppol_in
  type(elph_tr_type), intent(inout) :: elph_tr_ds
  type(mpi_type), intent(in) :: mpi_enreg
  real(dp), intent(in) :: kpt_in(3,nkpt_in)
 end subroutine read_el_veloc
end interface

interface
 subroutine read_gkk(elph_ds,Cryst,Bst,FSfullpqtofull,gkk_flag,n1wf,nband,phon_ds,ep_prt_yambo,unitgkk)
  use defs_elphon
  use defs_datatypes
  use m_crystal
  implicit none
  integer,intent(in) :: ep_prt_yambo
  integer,intent(in) :: n1wf
  integer,intent(in) :: nband
  integer,intent(in) :: unitgkk
  type(bandstructure_type),intent(in) :: Bst
  type(crystal_structure),intent(in) :: Cryst
  type(elph_type),intent(inout) :: elph_ds
  type(phon_type),intent(inout) :: phon_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  integer,intent(out) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch, &
  &         elph_ds%k_phon%nkpt,elph_ds%nsppol,elph_ds%nqptirred)
 end subroutine read_gkk
end interface

interface
 subroutine refineblk(acell,amu,anaddb_dtset,ddb_blk,&  
  &  dielt,gmet,gprim,indsym,iout,&  
  &  mpert,msym,natom,nsym,ntypat,rmet,rprim,&  
  &  symrec,symrel,tcpui,twalli,typat,&  
  &  ucvol,xred,zeff)
  use defs_basis
  use defs_abitypes
  use m_ddb_blk
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  type(ddb_blk_type), pointer :: ddb_blk
  real(dp),intent(in) :: tcpui
  real(dp),intent(in) :: twalli
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(inout) :: dielt(3,3)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym*natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  integer,intent(in) :: symrec(3,3,msym)
  integer,intent(in) :: symrel(3,3,msym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(inout) :: zeff(3,3,natom)
 end subroutine refineblk
end interface

interface
 subroutine relaxpol(blkflg,blkval,etotal,fred,iatfix,indsym,iout,istrfix,&  
  &  mpert,msize,msym,natfix,natom,nstrfix,nsym,ntypat,pel,relaxat,relaxstr,&  
  &  rprimd,strten,symrel,targetpol,typat,ucvol,xcart,xred,zion)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(in) :: msym
  integer,intent(in) :: natfix
  integer,intent(in) :: natom
  integer,intent(in) :: nstrfix
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: relaxat
  integer,intent(in) :: relaxstr
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: istrfix(6)
  integer,intent(in) :: blkflg(msize)
  real(dp),intent(inout) :: blkval(2,msize)
  real(dp),intent(in) :: fred(3,natom)
  integer,intent(in) :: iatfix(natom)
  integer,intent(in) :: indsym(4,msym,natom)
  real(dp),intent(in) :: pel(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: strten(6)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(inout) :: targetpol(3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xcart(3,natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
 end subroutine relaxpol
end interface

interface
 subroutine rsiaf9(acell,atifc,atmfrc,dielt,dipdip,dyewq0,&  
  &  gprim,ifcana,ifcout,iout,natom,nrpt,nsphere,rcan,&  
  &  rifcsph,rprim,rpt,tcpui,twalli,wghatm,zeff)
  use defs_basis
  implicit none
  integer,intent(in) :: dipdip
  integer,intent(in) :: ifcana
  integer,intent(in) :: ifcout
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsphere
  real(dp),intent(in) :: rifcsph
  real(dp),intent(in) :: tcpui
  real(dp),intent(in) :: twalli
  real(dp),intent(in) :: acell(3)
  integer,intent(in) :: atifc(natom)
  real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(in) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(inout) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine rsiaf9
end interface

interface
 subroutine sym_gkk(acell,kphon_full2full,kpt_phon,gkk_flag,gkk_qpt,&  
  &  gprim,indsym,mpert,natom,nbranch,nFSband,nkpt_phon,nqpt,nqptirred,nsym,&  
  &  qptirred,rprim,qpt_full,symrec,symrel,tnons,ucvol,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: nFSband
  integer,intent(in) :: natom
  integer,intent(in) :: nbranch
  integer,intent(in) :: nkpt_phon
  integer,intent(in) :: nqpt
  integer,intent(in) :: nqptirred
  integer,intent(in) :: nsym
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  integer,intent(inout) :: gkk_flag(nbranch,nkpt_phon,nqpt)
  real(dp),intent(inout) :: gkk_qpt(2,nbranch,nFSband,nFSband,nkpt_phon,nqpt)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: kphon_full2full(2,nsym,nkpt_phon)
  real(dp),intent(in) :: kpt_phon(3,nkpt_phon)
  real(dp),intent(in) :: qpt_full(3,nqpt)
  real(dp),intent(in) :: qptirred(3,nqptirred)
  real(dp),intent(in) :: rprim(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine sym_gkk
end interface

interface
 subroutine symdm9(blkflg,blknrm,blkqpt,blktyp,blkval,&  
  &  dynmat,gprim,indsym,mpert,natom,nblok,nqpt,nsym,rfmeth,&  
  &  rprim,spqpt,symrec,symrel)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  integer,intent(in) :: nqpt
  integer,intent(in) :: nsym
  integer,intent(in) :: rfmeth
  integer,intent(in) :: blkflg(3,mpert,3,mpert,nblok)
  real(dp),intent(in) :: blknrm(3,nblok)
  real(dp),intent(in) :: blkqpt(9,nblok)
  integer,intent(in) :: blktyp(nblok)
  real(dp),intent(in) :: blkval(2,3*mpert*3*mpert,nblok)
  real(dp),intent(out) :: dynmat(2,3,natom,3,natom,nqpt)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: spqpt(3,nqpt)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine symdm9
end interface

interface
 subroutine symgamma(elph_ds,kphon_full2full,h1_mat_el,&  
  &  indsym,natom,nsym,symq,symrec)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  type(elph_type),intent(in) :: elph_ds
  real(dp),intent(inout) :: h1_mat_el(2,elph_ds%nFSband,elph_ds%nFSband, &
  &         elph_ds%nbranch,elph_ds%k_phon%nkpt)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: kphon_full2full(2,nsym,elph_ds%k_phon%nkpt)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine symgamma
end interface

interface
 subroutine test_ftgkk(elph_ds,gprim,natom,nrpt,rpt,qpt_full,wghatm)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  type(elph_type),intent(inout) :: elph_ds
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qpt_full(3,elph_ds%nqpt_full)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine test_ftgkk
end interface

interface
 subroutine thm9(acell,amu,anaddb_dtset,atmfrc,dielt,displ,&  
  &  dyewq0,d2cart,eigval,eigvec,gmet,gprim,indsym,iout,mpert,msym,natom,&  
  &  nrpt,nsym,ntypat,outfilename_radix,phfrq,rmet,rprim,rpt,symrec,symrel,tcpui,&  
  &  trans,twalli,typat,ucvol,wghatm,xred,zeff,themflag, udispl, ufreq)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in),optional :: themflag
  integer,intent(in),optional :: udispl
  integer,intent(in),optional :: ufreq
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  character(len=fnlen) :: outfilename_radix
  real(dp),intent(in) :: tcpui
  real(dp),intent(in) :: twalli
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(in) :: dielt(3,3)
  real(dp),intent(out) :: displ(2*3*natom*3*natom)
  real(dp),intent(inout) :: dyewq0(3,3,natom)
  real(dp),intent(out) :: eigval(3*natom)
  real(dp),intent(out) :: eigvec(2,3,natom,3*natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  real(dp),intent(out) :: phfrq(3*natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: trans(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zeff(3,3,natom)
 end subroutine thm9
end interface

interface
 subroutine thmeig(g2fsmear,acell,amu,anaddb_dtset,d2asr,&  
  &  filnam,mband,mpert,msize,natom,nkpt,ntemper,&  
  &  ntypat,rprim,telphint,temperinc,&  
  &  tempermin,thmflag,typat,xred,&  
  &  ddb_blk,ddbun,dimekb,filnam5,iout,&  !new
  &  lmnmax,msym,nblok2,nsym,occopt,symrel,tnons,usepaw,zion,&  
  &  symrec,natifc,gmet,gprim,indsym,rmet,atifc,ucvol,xcart) !new
  use defs_basis
  use defs_abitypes
  use m_ddb_blk
  implicit none
  integer,intent(in) :: ddbun
  integer,intent(in) :: dimekb
  integer,intent(in) :: iout
  integer,intent(in) :: lmnmax
  integer,intent(in) :: mband
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(in) :: msym
  integer,intent(in) :: natifc
  integer,intent(inout) :: natom
  integer,intent(inout) :: nblok2
  integer,intent(inout) :: nkpt
  integer,intent(inout) :: nsym
  integer,intent(in) :: ntemper
  integer,intent(inout) :: ntypat
  integer,intent(inout) :: occopt
  integer,intent(in) :: telphint
  integer,intent(in) :: thmflag
  integer,intent(in) :: usepaw
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  type(ddb_blk_type), pointer :: ddb_blk
  character(len=fnlen),intent(in) :: filnam
  character(len=fnlen),intent(in) :: filnam5
  real(dp),intent(in) :: g2fsmear
  real(dp),intent(in) :: temperinc
  real(dp),intent(in) :: tempermin
  real(dp),intent(out) :: ucvol
  real(dp),intent(inout) :: acell(3)
  real(dp),intent(inout) :: amu(ntypat)
  integer,intent(inout) :: atifc(natom)
  real(dp),intent(inout) :: d2asr(2,3,natom,3,natom)
  real(dp),intent(out) :: gmet(3,3)
  real(dp),intent(out) :: gprim(3,3)
  integer,intent(out) :: indsym(4,nsym,natom)
  real(dp),intent(out) :: rmet(3,3)
  real(dp),intent(inout) :: rprim(3,3)
  integer,intent(out) :: symrec(3,3,msym)
  integer,intent(out) :: symrel(3,3,msym)
  real(dp),intent(out) :: tnons(3,msym)
  integer,intent(inout) :: typat(natom)
  real(dp),intent(out) :: xcart(3,natom)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(out) :: zion(ntypat)
 end subroutine thmeig
end interface

interface
 subroutine wght9(brav,gprim,natom,ngqpt,nqpt,nqshft,nrpt,qshft,rcan,rpt,wghatm)
  use defs_basis
  implicit none
  integer,intent(in) :: brav
  integer,intent(in) :: natom
  integer,intent(in) :: nqpt
  integer,intent(in) :: nqshft
  integer,intent(in) :: nrpt
  integer,intent(inout) :: ngqpt(9)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qshft(3,4)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(out) :: wghatm(natom,natom,nrpt)
 end subroutine wght9
end interface

end module interfaces_77_ddb
!!***
