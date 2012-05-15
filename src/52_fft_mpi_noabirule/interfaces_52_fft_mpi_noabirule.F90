!!****m* ABINIT/interfaces_52_fft_mpi_noabirule
!! NAME
!! interfaces_52_fft_mpi_noabirule
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/52_fft_mpi_noabirule
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

module interfaces_52_fft_mpi_noabirule

 implicit none

interface
 subroutine accrho(icplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc,&  
  &  max1,max2,max3,m1,m2,m3,md1,md2proc,md3,mpi_comm,nproc,iproc,paral_kgb,zf,rho,weight)
  use defs_basis
  implicit none
  integer :: icplexwf
  integer :: iproc
  integer :: m1
  integer :: m2
  integer :: m3
  integer :: max1
  integer :: max2
  integer :: max3
  integer :: md1
  integer :: md2proc
  integer :: md3
  integer :: mpi_comm
  integer :: n1
  integer :: n2
  integer :: n3
  integer :: nd1
  integer :: nd2
  integer :: nd3proc
  integer :: ndat
  integer :: nproc
  integer :: paral_kgb
  real(dp), dimension(nd1,nd2,nd3proc) :: rho
  real(dp), dimension(ndat) :: weight
  real(dp), dimension(2,md1,md3,md2proc,ndat) :: zf
 end subroutine accrho
end interface

interface
 subroutine addrho(icplexwf,includelast,nd1,nd2,n2,lot,n1dfft,zw,rhopart,weight)
  use defs_basis
  implicit none
  integer :: icplexwf
  integer :: includelast
  integer :: lot
  integer :: n1dfft
  integer :: n2
  integer :: nd1
  integer :: nd2
  real(dp) :: weight
  real(dp) :: rhopart(nd1,nd2)
  real(dp) :: zw(2,lot,n2)
 end subroutine addrho
end interface

interface
 subroutine applypot(icplexwf,icplex,ndat,n1,n2,n3,nd1,nd2,nd3proc,&  
  &  max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,&  
  &  max1o,max2o,max3o,m1o,m2o,m3o,mpi_comm,nproc,iproc,paral_kgb,pot,zf)
  use defs_basis
  implicit none
  integer :: icplex
  integer :: icplexwf
  integer :: iproc
  integer :: m1i
  integer :: m1o
  integer :: m2i
  integer :: m2o
  integer :: m3i
  integer :: m3o
  integer :: max1i
  integer :: max1o
  integer :: max2i
  integer :: max2o
  integer :: max3i
  integer :: max3o
  integer :: md1
  integer :: md2proc
  integer :: md3
  integer :: mpi_comm
  integer :: n1
  integer :: n2
  integer :: n3
  integer :: nd1
  integer :: nd2
  integer :: nd3proc
  integer :: ndat
  integer :: nproc
  integer :: paral_kgb
  real(kind=dp), dimension(icplex*nd1,nd2,nd3proc) :: pot
  real(kind=dp), dimension(2,md1,md3,md2proc,ndat) :: zf
 end subroutine applypot
end interface

interface
 subroutine back(icplex,mpi_enreg,ndat,n1,n2,n3,nd1,nd2,nd3,nd1eff,nd2proc,nd3proc,option,paral_kgb,zf,zr)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: icplex
  integer :: n1
  integer :: n2
  integer :: n3
  integer :: nd1
  integer :: nd1eff
  integer :: nd2
  integer :: nd2proc
  integer :: nd3
  integer :: nd3proc
  integer :: ndat
  integer :: option
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp), dimension(2,nd1,nd3,nd2proc,ndat) :: zf
  real(dp), dimension(2,nd1eff,nd2,nd3proc,ndat) :: zr
 end subroutine back
end interface

interface
 subroutine back_wf(icplexwf,mpi_enreg,ndat,n1,n2,n3,nd1,nd2,nd3proc,&  
  &  max1,max2,max3,m1,m2,m3,md1,md2proc,md3,nproc,iproc,paral_kgb,zf,zr)
  use defs_basis
  use defs_abitypes
  implicit none
  integer :: icplexwf
  integer :: iproc
  integer :: m1
  integer :: m2
  integer :: m3
  integer :: max1
  integer :: max2
  integer :: max3
  integer :: md1
  integer :: md2proc
  integer :: md3
  integer :: n1
  integer :: n2
  integer :: n3
  integer :: nd1
  integer :: nd2
  integer :: nd3proc
  integer :: ndat
  integer :: nproc
  integer :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp), dimension(2,md1,md3,md2proc,ndat) :: zf
  real(dp), dimension(2,nd1,nd2,nd3proc,ndat) :: zr
 end subroutine back_wf
end interface

interface
 subroutine ctrig(n,trig,after,before,now,isign,ic)
  use defs_basis
  use defs_fftdata
  implicit none
  integer :: ic
  integer :: isign
  integer :: n
  integer :: after(mdata)
  integer :: before(mdata)
  integer :: now(mdata)
  real(dp) :: trig(2,n)
 end subroutine ctrig
end interface

interface
 subroutine fftstp(mm,n1dfft,m,nn,n,zin,zout,trig,after,now,before,isign)
  use defs_basis
  implicit none
  integer :: after
  integer :: before
  integer :: isign
  integer :: m
  integer :: mm
  integer :: n
  integer :: n1dfft
  integer :: nn
  integer :: now
  real(dp) :: trig(2,n)
  real(dp) :: zin(2,mm,m)
  real(dp) :: zout(2,nn,n)
 end subroutine fftstp
end interface

interface
 subroutine fill(nd1,nd3,lot,n1dfft,n3,zf,zw)
  use defs_basis
  implicit none
  integer :: lot
  integer :: n1dfft
  integer :: n3
  integer :: nd1
  integer :: nd3
  real(dp) :: zf(2,nd1,nd3)
  real(dp) :: zw(2,lot,n3)
 end subroutine fill
end interface

interface
 subroutine fill_cent(md1,md3,lot,n1dfft,max3,m3,n3,zf,zw)
  use defs_basis
  implicit none
  integer :: lot
  integer :: m3
  integer :: max3
  integer :: md1
  integer :: md3
  integer :: n1dfft
  integer :: n3
  real(dp) :: zf(2,md1,md3)
  real(dp) :: zw(2,lot,n3)
 end subroutine fill_cent
end interface

interface
 subroutine fill_corn(md1,md3,lot,n1dfft,n3,zf,zw)
  use defs_basis
  implicit none
  integer :: lot
  integer :: md1
  integer :: md3
  integer :: n1dfft
  integer :: n3
  real(dp) :: zf(2,md1,md3)
  real(dp) :: zw(2,lot,n3)
 end subroutine fill_corn
end interface

interface
 subroutine fnorm(ndat,m1,m2,m3,md1,md2,md3,iproc,nproc,zf,sum)
  use defs_basis
  implicit none
  integer :: iproc
  integer :: m1
  integer :: m2
  integer :: m3
  integer :: md1
  integer :: md2
  integer :: md3
  integer :: ndat
  integer :: nproc
  real(dp) :: sum
  real(dp) :: zf(2,md1,md3,md2/nproc,ndat)
 end subroutine fnorm
end interface

interface
 subroutine forw(icplex,mpi_enreg,ndat,n1,n2,n3,nd1,nd2,nd3,nd1eff,nd2proc,nd3proc,option,paral_kgb,zr,zf)
  use defs_basis
  use defs_abitypes
  implicit none
  integer :: icplex
  integer :: n1
  integer :: n2
  integer :: n3
  integer :: nd1
  integer :: nd1eff
  integer :: nd2
  integer :: nd2proc
  integer :: nd3
  integer :: nd3proc
  integer :: ndat
  integer :: option
  integer :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp), dimension(2,nd1,nd3,nd2proc,ndat) :: zf
  real(dp), dimension(2,nd1eff,nd2,nd3proc,ndat) :: zr
 end subroutine forw
end interface

interface
 subroutine forw_wf(icplexwf,mpi_enreg,ndat,n1,n2,n3,nd1,nd2,nd3proc,&  
  &  max1,max2,max3,m1,m2,m3,md1,md2proc,md3,nproc,iproc,paral_kgb,zr,zf)
  use defs_basis
  use defs_abitypes
  implicit none
  integer :: icplexwf
  integer :: iproc
  integer :: m1
  integer :: m2
  integer :: m3
  integer :: max1
  integer :: max2
  integer :: max3
  integer :: md1
  integer :: md2proc
  integer :: md3
  integer :: n1
  integer :: n2
  integer :: n3
  integer :: nd1
  integer :: nd2
  integer :: nd3proc
  integer :: ndat
  integer :: nproc
  integer :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp), dimension(2,md1,md3,md2proc,ndat) :: zf
  real(dp), dimension(2,nd1,nd2,nd3proc,ndat) :: zr
 end subroutine forw_wf
end interface

interface
 subroutine indirect_parallel_Fourier(index,left,mpi_enreg,ngleft,ngright,nleft,nright,paral_kgb,right,sizeindex)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nleft
  integer,intent(in) :: nright
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: sizeindex
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngleft(18)
  integer,intent(in) :: ngright(18)
  integer,intent(in) :: index(sizeindex)
  real(dp),intent(inout) :: left(2,nleft)
  real(dp),intent(in) :: right(2,nright)
 end subroutine indirect_parallel_Fourier
end interface


interface
 subroutine initf(idat,m1,m2,m3,md1,md2,md3,iproc,nproc,zf)
  use defs_basis
  implicit none
  integer :: idat
  integer :: iproc
  integer :: m1
  integer :: m2
  integer :: m3
  integer :: md1
  integer :: md2
  integer :: md3
  integer :: nproc
  real(dp) :: zf(2,0:md1-1,0:md3-1,0:md2/nproc-1)
 end subroutine initf
end interface

interface
 subroutine mpiswitch(j3,n1dfft,Jp2st,J2st,lot,n1,nd2proc,&  
  &  nd3proc,nproc,ioption,zmpi1,zw)
  use defs_basis
  implicit none
  integer :: J2st
  integer :: Jp2st
  integer :: ioption
  integer :: j3
  integer :: lot
  integer :: n1
  integer :: n1dfft
  integer :: nd2proc
  integer :: nd3proc
  integer :: nproc
  real(dp) :: zmpi1(2,n1,nd2proc,nd3proc,nproc)
  real(dp) :: zw(2,lot,n1)
 end subroutine mpiswitch
end interface

interface
 subroutine mpiswitch_cent(j3,n1dfft,Jp2stb,J2stb,lot,max1,md1,m1,n1,md2proc,&  
  &  nd3proc,nproc,ioption,zmpi1,zw,max2,m2,n2)
  use defs_basis
  implicit none
  integer :: J2stb
  integer :: Jp2stb
  integer :: ioption
  integer :: j3
  integer :: lot
  integer :: m1
  integer :: m2
  integer :: max1
  integer :: max2
  integer :: md1
  integer :: md2proc
  integer :: n1
  integer :: n1dfft
  integer :: n2
  integer :: nd3proc
  integer :: nproc
  real(dp) :: zmpi1(2,md1,md2proc,nd3proc,nproc)
  real(dp) :: zw(2,lot,n1)
 end subroutine mpiswitch_cent
end interface

interface
 subroutine mpiswitch_corn(j3,n1dfft,Jp2stb,J2stb,lot,n1,md2,nd3,nproc,zmpi1,zw)
  use defs_basis
  implicit none
  integer :: J2stb
  integer :: Jp2stb
  integer :: j3
  integer :: lot
  integer :: md2
  integer :: n1
  integer :: n1dfft
  integer :: nd3
  integer :: nproc
  real(dp) :: zmpi1(2,n1/2,md2/nproc,nd3/nproc,nproc)
  real(dp) :: zw(2,lot,n1)
 end subroutine mpiswitch_corn
end interface

interface
 subroutine multpot(icplexwf,icplex,includelast,nd1,nd2,n2,lot,n1dfft,pot,zw)
  use defs_basis
  implicit none
  integer :: icplex
  integer :: icplexwf
  integer :: includelast
  integer :: lot
  integer :: n1dfft
  integer :: n2
  integer :: nd1
  integer :: nd2
  real(dp) :: pot(icplex*nd1,nd2)
  real(dp) :: zw(2,lot,n2)
 end subroutine multpot
end interface

interface
 subroutine rnorm(n1,n2,n3,nd1,nd2,nd3,iproc,nproc,rho,sum)
  use defs_basis
  implicit none
  integer :: iproc
  integer :: n1
  integer :: n2
  integer :: n3
  integer :: nd1
  integer :: nd2
  integer :: nd3
  integer :: nproc
  real(dp) :: sum
  real(dp) :: rho(nd1,nd2,nd3/nproc)
 end subroutine rnorm
end interface

interface
 subroutine scramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zw,zmpi2)
  use defs_basis
  implicit none
  integer :: i1
  integer :: j2
  integer :: lot
  integer :: md1
  integer :: md2proc
  integer :: n1dfft
  integer :: n3
  integer :: nnd3
  real(dp) :: zmpi2(2,md1,md2proc,nnd3)
  real(dp) :: zw(2,lot,n3)
 end subroutine scramble
end interface

interface
 subroutine set(n,zmpi)
  use defs_basis
  implicit none
  integer :: n
  real(dp) :: zmpi(n)
 end subroutine set
end interface

interface
 subroutine sg_ctrig(n,trig,aft,bef,now,ris,ic,ind,mfac,mg)
  use defs_basis
  implicit none
  integer,intent(out) :: ic
  integer,intent(in) :: mfac
  integer,intent(in) :: mg
  integer,intent(in) :: n
  real(dp),intent(in) :: ris
  integer,intent(out) :: aft(mfac)
  integer,intent(out) :: bef(mfac)
  integer,intent(out) :: ind(mg)
  integer,intent(out) :: now(mfac)
  real(dp),intent(out) :: trig(2,mg)
 end subroutine sg_ctrig
end interface

interface
 subroutine sg_fftx(fftcache,mfac,mg,nd1,nd2,nd3,n2,n3,z,zbr,&  
  &  trig,aft,now,bef,ris,ind,ic)
  use defs_basis
  implicit none
  integer,intent(in) :: fftcache
  integer,intent(in) :: ic
  integer,intent(in) :: mfac
  integer,intent(in) :: mg
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: nd1
  integer,intent(in) :: nd2
  integer,intent(in) :: nd3
  real(dp),intent(in) :: ris
  integer,intent(in) :: aft(mfac)
  integer,intent(in) :: bef(mfac)
  integer,intent(in) :: ind(mg)
  integer,intent(in) :: now(mfac)
  real(dp),intent(in) :: trig(2,mg)
  real(dp),intent(inout) :: z(2,nd1,nd2,nd3)
  real(dp),intent(inout) :: zbr(2,nd1,nd2,nd3)
 end subroutine sg_fftx
end interface

interface
 subroutine sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3,&  
  &  z,zbr,trig,aft,now,bef,ris,ind,ic)
  use defs_basis
  implicit none
  integer,intent(in) :: fftcache
  integer,intent(in) :: ic
  integer,intent(in) :: mfac
  integer,intent(in) :: mg
  integer,intent(in) :: n1
  integer,intent(in) :: n1i
  integer,intent(in) :: n3
  integer,intent(in) :: n3i
  integer,intent(in) :: nd1
  integer,intent(in) :: nd2
  integer,intent(in) :: nd3
  real(dp),intent(in) :: ris
  integer,intent(in) :: aft(mfac)
  integer,intent(in) :: bef(mfac)
  integer,intent(in) :: ind(mg)
  integer,intent(in) :: now(mfac)
  real(dp),intent(in) :: trig(2,mg)
  real(dp),intent(inout) :: z(2,nd1,nd2,nd3)
  real(dp),intent(inout) :: zbr(2,nd1,nd2,nd3)
 end subroutine sg_ffty
end interface

interface
 subroutine sg_fftz(mfac,mg,nd1,nd2,nd3,n1,n2i,n2,z,zbr,&  
  &  trig,aft,now,bef,ris,ind,ic)
  use defs_basis
  implicit none
  integer,intent(in) :: ic
  integer,intent(in) :: mfac
  integer,intent(in) :: mg
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n2i
  integer,intent(in) :: nd1
  integer,intent(in) :: nd2
  integer,intent(in) :: nd3
  real(dp),intent(in) :: ris
  integer,intent(in) :: aft(mfac)
  integer,intent(in) :: bef(mfac)
  integer,intent(in) :: ind(mg)
  integer,intent(in) :: now(mfac)
  real(dp),intent(in) :: trig(2,mg)
  real(dp),intent(inout) :: z(2,nd1,nd2,nd3)
  real(dp),intent(inout) :: zbr(2,nd1,nd2,nd3)
 end subroutine sg_fftz
end interface

interface
 subroutine slice(mpi_comm,nproc,iproc,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,zf,zr)
  use defs_basis
  implicit none
  integer :: iproc
  integer :: m1
  integer :: m2
  integer :: m3
  integer :: md1
  integer :: md2
  integer :: md3
  integer :: mpi_comm
  integer :: n1
  integer :: n2
  integer :: n3
  integer :: nd1
  integer :: nd2
  integer :: nd3
  integer :: nproc
  real(dp) :: zf
  real(dp) :: zr
 end subroutine slice
end interface

interface
 subroutine switch(n1dfft,n2,lot,n1,lzt,zt,zw)
  use defs_basis
  implicit none
  integer :: lot
  integer :: lzt
  integer :: n1
  integer :: n1dfft
  integer :: n2
  real(dp) :: zt(2,lzt,n1)
  real(dp) :: zw(2,lot,n2)
 end subroutine switch
end interface

interface
 subroutine switch_cent(n1dfft,max2,m2,n2,lot,n1,lzt,zt,zw)
  use defs_basis
  implicit none
  integer :: lot
  integer :: lzt
  integer :: m2
  integer :: max2
  integer :: n1
  integer :: n1dfft
  integer :: n2
  real(dp) :: zt(2,lzt,n1)
  real(dp) :: zw(2,lot,n2)
 end subroutine switch_cent
end interface

interface
 subroutine switch_corn(n1dfft,n2,lot,n1,lzt,zt,zw)
  use defs_basis
  implicit none
  integer :: lot
  integer :: lzt
  integer :: n1
  integer :: n1dfft
  integer :: n2
  real(dp) :: zt(2,lzt,n1)
  real(dp) :: zw(2,lot,n2)
 end subroutine switch_corn
end interface

interface
 subroutine switchreal(includelast,n1dfft,n2,n2eff,lot,n1zt,lzt,zt,zw)
  use defs_basis
  implicit none
  integer :: includelast
  integer :: lot
  integer :: lzt
  integer :: n1dfft
  integer :: n1zt
  integer :: n2
  integer :: n2eff
  real(dp) :: zt(2,lzt,n1zt)
  real(dp) :: zw(2,lot,n2)
 end subroutine switchreal
end interface

interface
 subroutine switchreal_cent(includelast,n1dfft,max2,n2,lot,n1zt,lzt,zt,zw)
  use defs_basis
  implicit none
  integer :: includelast
  integer :: lot
  integer :: lzt
  integer :: max2
  integer :: n1dfft
  integer :: n1zt
  integer :: n2
  real(dp) :: zt(2,lzt,n1zt)
  real(dp) :: zw(2,lot,n2)
 end subroutine switchreal_cent
end interface

interface
 subroutine unfill(nd1,nd3,lot,n1dfft,n3,zw,zf)
  use defs_basis
  implicit none
  integer :: lot
  integer :: n1dfft
  integer :: n3
  integer :: nd1
  integer :: nd3
  real(dp) :: zf(2,nd1,nd3)
  real(dp) :: zw(2,lot,n3)
 end subroutine unfill
end interface

interface
 subroutine unfill_cent(md1,md3,lot,n1dfft,max3,m3,n3,zw,zf)
  use defs_basis
  implicit none
  integer :: lot
  integer :: m3
  integer :: max3
  integer :: md1
  integer :: md3
  integer :: n1dfft
  integer :: n3
  real(dp) :: zf(2,md1,md3)
  real(dp) :: zw(2,lot,n3)
 end subroutine unfill_cent
end interface

interface
 subroutine unfill_corn(md1,md3,lot,n1dfft,n3,zw,zf)
  use defs_basis
  implicit none
  integer :: lot
  integer :: md1
  integer :: md3
  integer :: n1dfft
  integer :: n3
  real(dp) :: zf(2,md1,md3)
  real(dp) :: zw(2,lot,n3)
 end subroutine unfill_corn
end interface

interface
 subroutine unmpiswitch(j3,n1dfft,Jp2st,J2st,lot,n1,nd2proc,nd3proc,nproc,ioption,zw,zmpi1)
  use defs_basis
  implicit none
  integer :: J2st
  integer :: Jp2st
  integer :: ioption
  integer :: j3
  integer :: lot
  integer :: n1
  integer :: n1dfft
  integer :: nd2proc
  integer :: nd3proc
  integer :: nproc
  real(dp) :: zmpi1(2,n1,nd2proc,nd3proc,nproc)
  real(dp) :: zw(2,lot,n1)
 end subroutine unmpiswitch
end interface

interface
 subroutine unmpiswitch_cent(j3,n1dfft,Jp2stf,J2stf,lot,max1,md1,m1,n1,md2proc,nd3proc,nproc,ioption,zw,zmpi1)
  use defs_basis
  implicit none
  integer :: J2stf
  integer :: Jp2stf
  integer :: ioption
  integer :: j3
  integer :: lot
  integer :: m1
  integer :: max1
  integer :: md1
  integer :: md2proc
  integer :: n1
  integer :: n1dfft
  integer :: nd3proc
  integer :: nproc
  real(dp) :: zmpi1(2,md1,md2proc,nd3proc,nproc)
  real(dp) :: zw(2,lot,n1)
 end subroutine unmpiswitch_cent
end interface

interface
 subroutine unmpiswitch_corn(j3,n1dfft,Jp2stf,J2stf,lot,n1,md2,nd3,nproc,zw,zmpi1)
  use defs_basis
  implicit none
  integer :: J2stf
  integer :: Jp2stf
  integer :: j3
  integer :: lot
  integer :: md2
  integer :: n1
  integer :: n1dfft
  integer :: nd3
  integer :: nproc
  real(dp) :: zmpi1(2,n1/2,md2/nproc,nd3/nproc,nproc)
  real(dp) :: zw(2,lot,n1)
 end subroutine unmpiswitch_corn
end interface

interface
 subroutine unscramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zmpi2,zw)
  use defs_basis
  implicit none
  integer :: i1
  integer :: j2
  integer :: lot
  integer :: md1
  integer :: md2proc
  integer :: n1dfft
  integer :: n3
  integer :: nnd3
  real(dp) :: zmpi2(2,md1,md2proc,nnd3)
  real(dp) :: zw(2,lot,n3)
 end subroutine unscramble
end interface

interface
 subroutine unswitch(n1dfft,n2,lot,n1,lzt,zw,zt)
  use defs_basis
  implicit none
  integer :: lot
  integer :: lzt
  integer :: n1
  integer :: n1dfft
  integer :: n2
  real(dp) :: zt(2,lzt,n1)
  real(dp) :: zw(2,lot,n2)
 end subroutine unswitch
end interface

interface
 subroutine unswitch_cent(n1dfft,max2,m2,n2,lot,n1,lzt,zw,zt)
  use defs_basis
  implicit none
  integer :: lot
  integer :: lzt
  integer :: m2
  integer :: max2
  integer :: n1
  integer :: n1dfft
  integer :: n2
  real(dp) :: zt(2,lzt,n1)
  real(dp) :: zw(2,lot,n2)
 end subroutine unswitch_cent
end interface

interface
 subroutine unswitch_corn(n1dfft,n2,lot,n1,lzt,zw,zt)
  use defs_basis
  implicit none
  integer :: lot
  integer :: lzt
  integer :: n1
  integer :: n1dfft
  integer :: n2
  real(dp) :: zt(2,lzt,n1)
  real(dp) :: zw(2,lot,n2)
 end subroutine unswitch_corn
end interface

interface
 subroutine unswitchreal(n1dfft,n2,n2eff,lot,n1zt,lzt,zw,zt)
  use defs_basis
  implicit none
  integer :: lot
  integer :: lzt
  integer :: n1dfft
  integer :: n1zt
  integer :: n2
  integer :: n2eff
  real(dp) :: zt(2,lzt,n1zt)
  real(dp) :: zw(2,lot,n2)
 end subroutine unswitchreal
end interface

interface
 subroutine unswitchreal_cent(n1dfft,max2,n2,lot,n1zt,lzt,zw,zt)
  use defs_basis
  implicit none
  integer :: lot
  integer :: lzt
  integer :: max2
  integer :: n1dfft
  integer :: n1zt
  integer :: n2
  real(dp) :: zt(2,lzt,n1zt)
  real(dp) :: zw(2,lot,n2)
 end subroutine unswitchreal_cent
end interface

interface
 subroutine vgl(idat,scal,m1,m2,m3,md1,md2,md3,iproc,nproc,zf,sum)
  use defs_basis
  implicit none
  integer :: idat
  integer :: iproc
  integer :: m1
  integer :: m2
  integer :: m3
  integer :: md1
  integer :: md2
  integer :: md3
  integer :: nproc
  real(dp) :: scal
  real(dp) :: sum
  real(dp) :: zf(2,0:md1-1,0:md3-1,0:md2/nproc-1)
 end subroutine vgl
end interface

end module interfaces_52_fft_mpi_noabirule
!!***
