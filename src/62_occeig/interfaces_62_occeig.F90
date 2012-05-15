!!****m* ABINIT/interfaces_62_occeig
!! NAME
!! interfaces_62_occeig
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/62_occeig
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

module interfaces_62_occeig

 implicit none

interface
 subroutine dos_degeneratewfs(dos_fractions_in,dos_fractions_out,eigen,mband,nband,ndos,nkpt,nsppol)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: ndos
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  real(dp),intent(in) :: dos_fractions_in(nkpt,mband,nsppol,ndos)
  real(dp),intent(out) :: dos_fractions_out(nkpt,mband,nsppol,ndos)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine dos_degeneratewfs
end interface

interface
 subroutine dos_hdr_write(buffer,deltaene,dosdeltae,&  
  &  eigen,enemax,enemin,fermie,mband,nband,nene,&  
  &  nkpt,nsppol,occopt,prtdos,tphysel,tsmear,unitdos)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(out) :: nene
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: prtdos
  integer,intent(in) :: unitdos
  real(dp),intent(in) :: buffer
  real(dp),intent(out) :: deltaene
  real(dp),intent(in) :: dosdeltae
  real(dp),intent(out) :: enemax
  real(dp),intent(out) :: enemin
  real(dp),intent(in) :: fermie
  real(dp),intent(in) :: tphysel
  real(dp),intent(in) :: tsmear
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
 end subroutine dos_hdr_write
end interface

interface
 subroutine gaus_dos(dos_fractions,&  
  &  dos_fractions_paw1,dos_fractions_pawt1,dtset,fermie,eigen,&  
  &  fildata,mbesslang,m_dos_flag,ndosfraction,paw_dos_flag) 
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: m_dos_flag
  integer,intent(in) :: mbesslang
  integer,intent(in) :: ndosfraction
  integer,intent(in) :: paw_dos_flag
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(in) :: fildata
  real(dp),intent(in) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
  real(dp),intent(in) :: dos_fractions_paw1(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*paw_dos_flag)
  real(dp),intent(in) :: dos_fractions_pawt1(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*paw_dos_flag)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 end subroutine gaus_dos
end interface

interface
 subroutine get_dos_1band (dos_fractions,enemin,enemax,&  
  &  integ_dos,nene,nkpt,ndosfraction,&  
  &  partial_dos,tweight,dtweightde)
  use defs_basis
  implicit none
  integer,intent(in) :: ndosfraction
  integer,intent(in) :: nene
  integer,intent(in) :: nkpt
  real(dp),intent(in) :: enemax
  real(dp),intent(in) :: enemin
  real(dp),intent(in) :: dos_fractions(nkpt,ndosfraction)
  real(dp),intent(in) :: dtweightde(nkpt,nene)
  real(dp),intent(out) :: integ_dos(nene,ndosfraction)
  real(dp),intent(out) :: partial_dos(nene,ndosfraction)
  real(dp),intent(in) :: tweight(nkpt,nene)
 end subroutine get_dos_1band
end interface

interface
 subroutine get_dos_1band_m (dos_fractions_m,enemin,enemax,&  
  &  integ_dos_m,nene,nkpt,ndosfraction_m,&  
  &  partial_dos_m,tweight,dtweightde)
  use defs_basis
  implicit none
  integer,intent(in) :: ndosfraction_m
  integer,intent(in) :: nene
  integer,intent(in) :: nkpt
  real(dp),intent(in) :: enemax
  real(dp),intent(in) :: enemin
  real(dp),intent(in) :: dos_fractions_m(nkpt,ndosfraction_m)
  real(dp),intent(in) :: dtweightde(nkpt,nene)
  real(dp),intent(out) :: integ_dos_m(nene,ndosfraction_m)
  real(dp),intent(out) :: partial_dos_m(nene,ndosfraction_m)
  real(dp),intent(in) :: tweight(nkpt,nene)
 end subroutine get_dos_1band_m
end interface

interface
 subroutine get_fsurf_1band(dtset,eigen_in,fermie,klatt,kpt_fullbz,&  
  &  mtetra,nfiner,nkpt_fullbz,tetra_full,tetra_wrap,tolfermi)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mtetra
  integer,intent(in) :: nfiner
  integer,intent(in) :: nkpt_fullbz
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fermie
  real(dp),intent(in) :: tolfermi
  real(dp),intent(in) :: eigen_in(dtset%nkpt)
  real(dp),intent(in) :: klatt(3,3)
  real(dp),intent(in) :: kpt_fullbz(3,nkpt_fullbz)
  integer,intent(in) :: tetra_full(4,2,mtetra)
  integer,intent(in) :: tetra_wrap(3,4,mtetra)
 end subroutine get_fsurf_1band
end interface

interface
 subroutine getnel(doccde,dosdeltae,eigen,entropy,fermie,maxocc,mband,nband,&  
  &  nelect,nkpt,nsppol,occ,occopt,option,tphysel,tsmear,unitdos,wtk)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: option
  integer,intent(in) :: unitdos
  real(dp),intent(in) :: dosdeltae
  real(dp),intent(out) :: entropy
  real(dp),intent(in) :: fermie
  real(dp),intent(in) :: maxocc
  real(dp),intent(out) :: nelect
  real(dp),intent(in) :: tphysel
  real(dp),intent(in) :: tsmear
  real(dp),intent(out) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(out) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine getnel
end interface

interface
 subroutine init_bess_spl(mbess,bessint_delta,mlang,&  
  &  bess_spl,bess_spl_der,x_bess)
  use defs_basis
  implicit none
  integer,intent(in) :: mbess
  integer,intent(in) :: mlang
  real(dp),intent(in) :: bessint_delta
  real(dp),intent(out) :: bess_spl(mbess,mlang)
  real(dp),intent(out) :: bess_spl_der(mbess,mlang)
  real(dp),intent(out) :: x_bess(mbess)
 end subroutine init_bess_spl
end interface

interface
 subroutine init_occ_ent(entfun, limit,&  
  &  nptsdiv2, occfun,occopt,option,smdfun,tphysel,&  
  &  tsmear,tsmearinv, xgrid)
  use defs_basis
  implicit none
  integer,intent(inout) :: nptsdiv2
  integer,intent(in) :: occopt
  integer,intent(in) :: option
  real(dp),intent(out) :: limit
  real(dp),intent(in) :: tphysel
  real(dp),intent(in) :: tsmear
  real(dp),intent(out) :: tsmearinv
  real(dp),intent(inout) :: entfun(-nptsdiv2:nptsdiv2,2)
  real(dp),intent(inout) :: occfun(-nptsdiv2:nptsdiv2,2)
  real(dp),intent(inout) :: smdfun(-nptsdiv2:nptsdiv2,2)
  real(dp),intent(inout) :: xgrid(-nptsdiv2:nptsdiv2)
 end subroutine init_occ_ent
end interface

interface
 subroutine init_ylm_spl(mbessint,bessargmax,bessint_delta,mlang,spl_bessint)
  use defs_basis
  implicit none
  integer,intent(in) :: mbessint
  integer,intent(in) :: mlang
  real(dp),intent(in) :: bessargmax
  real(dp),intent(in) :: bessint_delta
  real(dp),intent(out) :: spl_bessint(mbessint,mlang)
 end subroutine init_ylm_spl
end interface

interface
 subroutine newocc(doccde,eigen,entropy,fermie,fixmom,mband,nband,&  
  &  nelect,nkpt,nspinor,nsppol,occ,occopt,prtvol,stmbias,tphysel,tsmear,wtk)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: occopt
  integer,intent(in) :: prtvol
  real(dp),intent(out) :: entropy
  real(dp),intent(out) :: fermie
  real(dp),intent(in) :: fixmom
  real(dp),intent(in) :: nelect
  real(dp),intent(in) :: stmbias
  real(dp),intent(in) :: tphysel
  real(dp),intent(in) :: tsmear
  real(dp),intent(out) :: doccde(mband*nkpt*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(out) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: wtk(nkpt)
 end subroutine newocc
end interface

interface
 subroutine occeig(doccde_k,doccde_kq,eig0_k,eig0_kq,nband_k,&  
  &  occopt,occ_k,occ_kq,rocceig)
  use defs_basis
  implicit none
  integer,intent(in) :: nband_k
  integer,intent(in) :: occopt
  real(dp),intent(in) :: doccde_k(nband_k)
  real(dp),intent(in) :: doccde_kq(nband_k)
  real(dp),intent(in) :: eig0_k(nband_k)
  real(dp),intent(in) :: eig0_kq(nband_k)
  real(dp),intent(in) :: occ_k(nband_k)
  real(dp),intent(in) :: occ_kq(nband_k)
  real(dp),intent(out) :: rocceig(nband_k,nband_k)
 end subroutine occeig
end interface

interface
 subroutine pareigocc(eigen,formeig,localrdwf,mpi_enreg,mband,nband,nkpt,nsppol,occ,transmit_occ)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: localrdwf
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: transmit_occ
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: eigen(mband*(2*mband)**formeig*nkpt*nsppol)
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
 end subroutine pareigocc
end interface

interface
 subroutine printbxsf(eigen,ewind,fermie,gprimd,kptrlatt,mband,&  
  &  nkptirred,kptirred,nsym,use_afm,symrec,symafm,use_tr,nsppol,shiftk,nshiftk,fname,ierr)
  use defs_basis
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: mband
  integer,intent(in) :: nkptirred
  integer,intent(in) :: nshiftk
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  real(dp),intent(in) :: ewind
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(in) :: fname
  logical,intent(in) :: use_afm
  logical,intent(in) :: use_tr
  integer,intent(in) :: kptrlatt(3,3)
  real(dp),intent(in) :: eigen(mband,nkptirred,nsppol)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kptirred(3,nkptirred)
  real(dp),intent(in) :: shiftk(3,nshiftk)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine printbxsf
end interface

interface
 subroutine prtbltztrp_out (eigen, fermie, fname_radix, kpt,&  
  &  natom, nband, nelec, nkpt, nspinor, nsppol, nsym,&  
  &  rprimd, symrel)
  use defs_basis
  implicit none
  integer, intent(in) :: natom
  integer, intent(in) :: nband
  integer, intent(in) :: nkpt
  integer, intent(in) :: nspinor
  integer, intent(in) :: nsppol
  integer, intent(in) :: nsym
  real(dp), intent(in) :: fermie
  character(len=fnlen), intent(in) :: fname_radix
  real(dp), intent(in) :: nelec
  real(dp), intent(in) :: eigen(nband, nkpt, nsppol)
  real(dp), intent(in) :: kpt(3,nkpt)
  real(dp), intent(in) :: rprimd(3,3)
  integer, intent(in) :: symrel(3,3,nsym)
 end subroutine prtbltztrp_out
end interface

interface
 subroutine recip_ylm (bess_fit,cgcband,iatsph,istwfk,&  
  &  nradint,nradintmax,mlang,mpi_enreg,mpw,natom,natsph,&  
  &  npw_k,ntypat,ph3d,prtsphere,rint,rmax,sum_1ll_1atom,sum_1lm_1atom,typat,ucvol,ylm,znucl)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: istwfk
  integer,intent(in) :: mlang
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: natsph
  integer,intent(in) :: npw_k
  integer,intent(in) :: nradintmax
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtsphere
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: bess_fit(mpw,nradintmax,mlang)
  real(dp),intent(in) :: cgcband(2,npw_k)
  integer,intent(in) :: iatsph(natsph)
  integer,intent(in) :: nradint(natsph)
  real(dp),intent(in) :: ph3d(2,npw_k,natom)
  real(dp),intent(in) :: rint(nradintmax)
  real(dp),intent(in) :: rmax(natom)
  real(dp),intent(out) :: sum_1ll_1atom(mlang,natsph)
  real(dp),intent(out) :: sum_1lm_1atom(mlang*mlang,natsph)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: ylm(mpw,mlang*mlang)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine recip_ylm
end interface

interface
 subroutine sphericaldens(fofg,gnorm,nfft,rmax,sphfofg)
  use defs_basis
  implicit none
  integer,intent(in) :: nfft
  real(dp),intent(in) :: rmax
  real(dp),intent(in) :: fofg(2,nfft)
  real(dp),intent(in) :: gnorm(nfft)
  real(dp),intent(out) :: sphfofg(2,nfft)
 end subroutine sphericaldens
end interface

interface
 subroutine tetrahedron (dos_fractions,dos_fractions_m,dos_fractions_paw1,dos_fractions_pawt1,&  
  &  dtset,fermie,eigen,fildata,mbesslang,m_dos_flag,ndosfraction,paw_dos_flag,rprimd)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: m_dos_flag
  integer,intent(in) :: mbesslang
  integer,intent(in) :: ndosfraction
  integer,intent(in) :: paw_dos_flag
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: fermie
  character(len=fnlen),intent(in) :: fildata
  real(dp),intent(in) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
  real(dp),intent(in) :: dos_fractions_m(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*mbesslang*m_dos_flag)
  real(dp),intent(in) :: dos_fractions_paw1(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*paw_dos_flag)
  real(dp),intent(in) :: dos_fractions_pawt1(dtset%nkpt,dtset%mband, &
  &         dtset%nsppol,ndosfraction*paw_dos_flag)
  real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine tetrahedron
end interface

end module interfaces_62_occeig
!!***
