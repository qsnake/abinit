!!****m* ABINIT/interfaces_68_recursion
!! NAME
!! interfaces_68_recursion
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/68_recursion
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

module interfaces_68_recursion

 implicit none

interface
 subroutine density_rec(an,bn2,rho_out,nrec,&  
  &  fermie,tsmear,rtrotter,&  
  &  dim_trott,tol,inf_ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: dim_trott
  integer,intent(in) :: nrec
  real(dp),intent(in) :: fermie
  real(dp),intent(in) :: inf_ucvol
  real(dp), intent(out) :: rho_out
  real(dp),intent(in) :: rtrotter
  real(dp),intent(in) :: tol
  real(dp),intent(in) :: tsmear
  real(dp),intent(in) :: an(0:nrec)
  real(dp),intent(in) :: bn2(0:nrec)
 end subroutine density_rec
end interface

interface
 subroutine entropyrec(an,bn2,nrec,trotter,ent_out,multce,debug_rec,&  
  &  n_pt_integ,xmax,&  
  &  ent_out1,ent_out2,ent_out3,ent_out4)
  use defs_basis
  implicit none
  integer,intent(in) :: n_pt_integ
  integer,intent(in) :: nrec
  integer,intent(in) :: trotter
  logical,intent(in) :: debug_rec
  real(dp),intent(out) :: ent_out
  real(dp),intent(out) :: ent_out1
  real(dp),intent(out) :: ent_out2
  real(dp),intent(out) :: ent_out3
  real(dp),intent(out) :: ent_out4
  real(dp), intent(in) :: multce
  real(dp), intent(in) :: xmax
  real(dp),intent(in) :: an(0:nrec)
  real(dp),intent(in) :: bn2(0:nrec)
 end subroutine entropyrec
end interface

interface
 subroutine fermisolverec(fermie,rho,a,b2,debug_rec,nb_rec,&  
  &  temperature,trotter,nelect,&  
  &  acc, max_it,&  
  &  long_tranche,mpi_enreg,&  
  &  inf_ucvol,gputopo)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: long_tranche
  integer,intent(in) :: max_it
  integer,intent(in) :: nb_rec
  integer,intent(in) :: trotter
  real(dp),intent(in) :: acc
  logical,intent(in) :: debug_rec
  real(dp), intent(inout) :: fermie
  logical,intent(in) :: gputopo
  real(dp),intent(in) :: inf_ucvol
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: nelect
  real(dp),intent(in) :: temperature
  real(dp), intent(inout) :: a(0:nb_rec,long_tranche)
  real(dp), intent(inout) :: b2(0:nb_rec,long_tranche)
  real(dp), intent(inout) :: rho(long_tranche)
 end subroutine fermisolverec
end interface

interface
 subroutine first_rec(dtset,psps,rset)
  use defs_rectypes
  use defs_abitypes
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(pseudopotential_type),intent(in) :: psps
  type(recursion_type),intent(inout) :: rset
 end subroutine first_rec
end interface

interface
 subroutine getngrec(ngfft,rmet,ngfftrec,nfftrec,recrcut,delta,tronc)
  use defs_basis
  implicit none
  integer,intent(out) :: nfftrec
  real(dp),intent(in) :: delta
  real(dp),intent(in) :: recrcut
  logical,intent(out) :: tronc
  integer,intent(in) :: ngfft(18)
  integer,intent(out) :: ngfftrec(18)
  real(dp),intent(in) :: rmet(3,3)
 end subroutine getngrec
end interface

interface
 subroutine gran_potrec(an,bn2,nrec,trotter,ene_out, mult,&  
  &  debug_rec,n_pt_integ,xmax,&  
  &  ene_out1,ene_out2,ene_out3,ene_out4)
  use defs_basis
  implicit none
  integer,intent(in) :: n_pt_integ
  integer,intent(in) :: nrec
  integer,intent(in) :: trotter
  logical,intent(in) :: debug_rec
  real(dp),intent(out) :: ene_out
  real(dp),intent(out) :: ene_out1
  real(dp),intent(out) :: ene_out2
  real(dp),intent(out) :: ene_out3
  real(dp),intent(out) :: ene_out4
  real(dp), intent(in) :: mult
  real(dp), intent(in) :: xmax
  real(dp), intent(in) :: an(0:nrec)
  real(dp), intent(in) :: bn2(0:nrec)
 end subroutine gran_potrec
end interface

interface
 subroutine green_kernel(ZT_p,inf_rmet,inf_ucvol,mult,mpi_enreg,ngfft,nfft)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  real(dp),intent(in) :: inf_ucvol
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: mult
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: ZT_p(1:2,0:nfft-1)
  real(dp),intent(in) :: inf_rmet(3,3)
 end subroutine green_kernel
end interface

interface
 subroutine nlenergyrec(rset,enl,exppot,ngfft,natom,typat,&  
  &  tsmear,trotter,tol)
  use defs_basis
  use defs_rectypes
  implicit none
  integer , intent(in) :: natom
  integer , intent(in) :: trotter
  real(dp), intent(out) :: enl
  type(recursion_type),intent(in) :: rset
  real(dp), intent(in) :: tol
  real(dp), intent(in) :: tsmear
  integer , intent(in) :: ngfft(18)
  real(dp), intent(in) :: exppot(0:ngfft(1)*ngfft(2)*ngfft(3)-1)
  integer , intent(in) :: typat(natom)
 end subroutine nlenergyrec
end interface

interface
 subroutine pspnl_hgh_rec(psps,temperature,nlrec,debug)
  use defs_basis
  use defs_rectypes
  use defs_datatypes
  implicit none
  logical,intent(in) :: debug
  type(nlpsprec_type),intent(inout) :: nlrec
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: temperature
 end subroutine pspnl_hgh_rec
end interface

interface
 subroutine pspnl_operat_rec(nlrec,metrec,ngfftrec,debug)
  use defs_rectypes
  implicit none
  logical,intent(in) :: debug
  type(metricrec_type),intent(in) :: metrec
  type(nlpsprec_type),intent(inout) :: nlrec
  integer,intent(in) :: ngfftrec(18)
 end subroutine pspnl_operat_rec
end interface

interface
 subroutine recursion(exppot,coordx,coordy,coordz,an,bn2,rho_out,&  
  &  nrec,fermie,tsmear,rtrotter,dim_trott,&  
  &  ZT_p, tol,typat,&  
  &  nlrec,mpi_enreg,&  
  &  nfft,ngfft,metrec,&  
  &  tim_fourdp,natom,projec,tim)
  use defs_basis
  use defs_rectypes
  use defs_abitypes
  implicit none
  integer,intent(in) :: coordx
  integer,intent(in) :: coordy
  integer,intent(in) :: coordz
  integer,intent(in) :: dim_trott
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nrec
  integer,intent(in) :: tim
  integer,intent(in) :: tim_fourdp
  real(dp),intent(in) :: fermie
  type(metricrec_type),intent(in) :: metrec
  type(mpi_type),intent(inout) :: mpi_enreg
  type(nlpsprec_type),intent(in) :: nlrec
  real(dp), intent(out) :: rho_out
  real(dp),intent(in) :: rtrotter
  real(dp),intent(in) :: tol
  real(dp),intent(in) :: tsmear
  integer, intent(in) :: ngfft(18)
  real(dp), intent(in) :: ZT_p(1:2, 0:nfft-1)
  real(dp), intent(out) :: an(0:nrec)
  real(dp), intent(out) :: bn2(0:nrec)
  real(dp), intent(in) :: exppot(0:nfft-1)
  real(dp), intent(in) :: projec(0:,0:,0:,1:,1:)
  integer, intent(in) :: typat(natom)
 end subroutine recursion
end interface

interface
 subroutine recursion_nl(exppot,un,rho_out,rset,ngfft,&  
  &  tsmear,trotter,dim_trott,tol,typat,&  
  &  natom,projec)
  use defs_basis
  use defs_rectypes
  implicit none
  integer,intent(in) :: dim_trott
  integer,intent(in) :: natom
  integer,intent(in) :: trotter
  real(dp), intent(out) :: rho_out
  type(recursion_type),intent(in) :: rset
  real(dp),intent(in) :: tol
  real(dp),intent(in) :: tsmear
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: exppot(0:ngfft(1)*ngfft(2)*ngfft(3)-1)
  real(dp),pointer :: projec(:,:,:,:,:)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: un(0:rset%nfftrec-1)
 end subroutine recursion_nl
end interface

interface
 subroutine vn_nl_rec(vn,natom,typat,ngfftrec,&  
  &  inf_ucvol,nlrec,projec)
  use defs_basis
  use defs_rectypes
  implicit none
  integer,intent(in) :: natom
  real(dp),intent(in) :: inf_ucvol
  type(nlpsprec_type),intent(in) :: nlrec
  integer,intent(in) :: ngfftrec(3)
  real(dp),intent(in) :: projec(0:,0:,0:,1:,1:)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: vn(0:ngfftrec(1)*ngfftrec(2)*ngfftrec(3)-1)
 end subroutine vn_nl_rec
end interface

interface
 subroutine vtorhorec(dtset,&  
  &  ek,enl,entropy,e_eigenvalues,fermie,&  
  &  grnl,initialized,irrzon,nfftf,phnons,&  
  &  rhog, rhor, vtrial,rset,deltastep,rprimd,gprimd)
  use defs_basis
  use defs_rectypes
  use defs_abitypes
  implicit none
  integer,intent(in) :: deltastep
  integer,intent(in) :: initialized
  integer,intent(in) :: nfftf
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: e_eigenvalues
  real(dp),intent(out) :: ek
  real(dp),intent(out) :: enl
  real(dp),intent(out) :: entropy
  real(dp),intent(out) :: fermie
  type(recursion_type),intent(inout) :: rset
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: grnl(3*dtset%natom)
  integer, intent(in) :: irrzon((nfftf)**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp),intent(in) :: phnons(2,(nfftf)**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp),intent(inout) :: rhog(2,nfftf)
  real(dp),intent(out) :: rhor(nfftf,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: vtrial(nfftf,dtset%nspden)
 end subroutine vtorhorec
end interface

end module interfaces_68_recursion
!!***
