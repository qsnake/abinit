!!****m* ABINIT/interfaces_53_spacepar
!! NAME
!! interfaces_53_spacepar
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/53_spacepar
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

module interfaces_53_spacepar

 implicit none

interface
 subroutine dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw,option,vect1,vect2)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: npw
  integer,intent(in) :: option
  real(dp),intent(out) :: doti
  real(dp),intent(out) :: dotr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: vect1(2,npw)
  real(dp),intent(in) :: vect2(2,npw)
 end subroutine dotprod_g
end interface

interface
 subroutine dotprod_v(cplex,dotr,mpi_enreg,nfft,nspden,opt_storage,pot1,pot2)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: opt_storage
  real(dp),intent(out) :: dotr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: pot1(cplex*nfft,nspden)
  real(dp),intent(in) :: pot2(cplex*nfft,nspden)
 end subroutine dotprod_v
end interface

interface
 subroutine dotprod_vn(cplex,dens,dotr,doti,mpi_enreg,nfft,nfftot,nspden,option,pot,ucvol)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  real(dp),intent(out) :: doti
  real(dp),intent(out) :: dotr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: dens(cplex*nfft,nspden)
  real(dp),intent(in) :: pot(cplex*nfft,nspden)
 end subroutine dotprod_vn
end interface

interface
 subroutine matrixelmt_g(ai,ar,diag,istwf_k,mpi_enreg,needimag,npw,nspinor,vect1,vect2)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: needimag
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  real(dp),intent(out) :: ai
  real(dp),intent(out) :: ar
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: diag(npw)
  real(dp),intent(in) :: vect1(2,npw*nspinor)
  real(dp),intent(in) :: vect2(2,npw*nspinor)
 end subroutine matrixelmt_g
end interface

interface
 subroutine mean_fftr(arraysp,meansp,mpi_enreg,nfft,nfftot,nspden)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: arraysp(nfft,nspden)
  real(dp),intent(out) :: meansp(nspden)
 end subroutine mean_fftr
end interface

interface
 subroutine meanvalue_g(ar,diag,filter,istwf_k,mpi_enreg,npw,nspinor,vect,vect1,use_ndo)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: filter
  integer,intent(in) :: istwf_k
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: use_ndo
  real(dp),intent(out) :: ar
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: diag(npw)
  real(dp),intent(in) :: vect(2,npw*nspinor)
  real(dp),intent(in) :: vect1(2,npw*nspinor)
 end subroutine meanvalue_g
end interface

interface
 subroutine multipoles_fftr(arraysp,dipole,mpi_enreg,nfft,ngfft,nspden,rprimd,neworigin)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(3)
  real(dp),intent(in) :: arraysp(nfft,nspden)
  real(dp),intent(out) :: dipole(3,nspden)
  real(dp),intent(in) :: neworigin(3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine multipoles_fftr
end interface

interface
 subroutine multipoles_out(arraysp,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,rprimd,typat,ucvol,xred,ziontypat)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp), intent(in) :: ucvol
  integer,intent(in) :: ngfft(3)
  real(dp),intent(in) :: arraysp(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer, intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: ziontypat(ntypat)
 end subroutine multipoles_out
end interface

interface
 subroutine overlap_g(doti,dotr,mpw,npw_k1,pwind_k,vect1,vect2)
  use defs_basis
  implicit none
  integer,intent(in) :: mpw
  integer,intent(in) :: npw_k1
  real(dp),intent(out) :: doti
  real(dp),intent(out) :: dotr
  integer,intent(in) :: pwind_k(mpw)
  real(dp),intent(in) :: vect1(1:2,0:mpw)
  real(dp),intent(in) :: vect2(1:2,0:mpw)
 end subroutine overlap_g
end interface

interface
 subroutine sqnorm_g(dotr,istwf_k,mpi_enreg,npwsp,vect)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: npwsp
  real(dp),intent(out) :: dotr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: vect(2,npwsp)
 end subroutine sqnorm_g
end interface

interface
 subroutine sqnorm_v(cplex,mpi_enreg,nfft,norm2,nspden,opt_storage,pot)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: opt_storage
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(out) :: norm2
  real(dp),intent(in) :: pot(cplex*nfft,nspden)
 end subroutine sqnorm_v
end interface

end module interfaces_53_spacepar
!!***
