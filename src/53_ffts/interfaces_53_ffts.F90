!!****m* ABINIT/interfaces_53_ffts
!! NAME
!! interfaces_53_ffts
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/53_ffts
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

module interfaces_53_ffts

 implicit none

interface
 subroutine ccfft(fftalga,fftcache,inplace,isign,mpi_enreg,normalized,&  
  &  n1,n2,n3,n4,n5,n6,ndat,option,paral_kgb,work1,work2)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: fftalga
  integer,intent(in) :: fftcache
  integer,intent(out) :: inplace
  integer,intent(in) :: isign
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: ndat
  integer,intent(out) :: normalized
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: work1(2,n4*n5*n6*ndat)
  real(dp),intent(out) :: work2(2,n4*n5*n6*ndat)
 end subroutine ccfft
end interface

interface
 subroutine fftpac(ispden,nspden,n1,n2,n3,nd1,nd2,nd3,ngfft,aa,bb,option)
  use defs_basis
  implicit none
  integer,intent(in) :: ispden
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: nd1
  integer,intent(in) :: nd2
  integer,intent(in) :: nd3
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: aa(n1*n2*n3/ngfft(10),nspden)
  real(dp),intent(inout) :: bb(nd1,nd2,nd3)
 end subroutine fftpac
end interface

interface
 subroutine fftw(n1,n2,n3,isign,work1,work2)
  use defs_basis
  implicit none
  integer,intent(in) :: isign
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  real(dp),intent(in) :: work1(:,:,:,:)
  real(dp),intent(out) :: work2(:,:,:,:)
 end subroutine fftw
end interface

interface
 subroutine fftw3_fourdp(cplex,nx,ny,nz,isign,fofg,fofr,fftw_flags)
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,optional,intent(in) :: fftw_flags
  integer,intent(in) :: isign
  integer,intent(in) :: nx
  integer,intent(in) :: ny
  integer,intent(in) :: nz
  real(dp),intent(inout) :: fofg(2,nx*ny*nz)
  real(dp),intent(inout) :: fofr(cplex*nx*ny*nz)
 end subroutine fftw3_fourdp
end interface

interface
 subroutine fourdp_c2c_ip(ff,isign,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: isign
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  complex(dpc),intent(inout) :: ff(nfft)
 end subroutine fourdp_c2c_ip
end interface

interface
 subroutine fourdp_c2c_op(ff,gg,isign,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: isign
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  complex(dpc),intent(in) :: ff(nfft)
  complex(dpc),intent(out) :: gg(nfft)
 end subroutine fourdp_c2c_op
end interface

interface
 subroutine fftw3_fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&  
  &  kg_kin,kg_kout,mgfft,mpi_enreg,ndat,ngfft,npwin,npwout,n4,n5,n6,option,paral_kgb,weight_r,weight_i,&  
  &  use_ndo,fofginb) ! optional Arguments.
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mgfft
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: ndat
  integer,intent(in) :: npwin
  integer,intent(in) :: npwout
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  integer,intent(in),optional :: use_ndo
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: weight_i
  real(dp),intent(in) :: weight_r
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: denpot(cplex*n4,n5,n6)
  real(dp),intent(inout) :: fofgin(2,npwin*ndat)
  real(dp),intent(inout),optional :: fofginb(:,:)
  real(dp),intent(out) :: fofgout(2,npwout*ndat)
  real(dp),intent(inout) :: fofr(2,n4,n5,n6*ndat)
  integer,intent(in) :: gboundin(2*mgfft+8,2)
  integer,intent(in) :: gboundout(2*mgfft+8,2)
  integer,intent(in) :: kg_kin(3,npwin)
  integer,intent(in) :: kg_kout(3,npwout)
 end subroutine fftw3_fourwf
end interface

interface
 subroutine padded_fourwf_cplx(ff,ngfft,nx,ny,nz,ldx,ldy,ldz,mgfft,isign,gbound)
  use defs_basis
  implicit none
  integer,intent(in) :: isign
  integer,intent(in) :: ldx
  integer,intent(in) :: ldy
  integer,intent(in) :: ldz
  integer,intent(in) :: mgfft
  integer,intent(in) :: nx
  integer,intent(in) :: ny
  integer,intent(in) :: nz
  integer,intent(in) :: ngfft(18)
  complex(dpc),intent(inout) :: ff(ldx*ldy*ldz)
  integer,intent(in) :: gbound(2*mgfft+8,2)
 end subroutine padded_fourwf_cplx
end interface

interface
 subroutine fourdp(cplex,fofg,fofr,isign,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: isign
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: fofg(2,nfft)
  real(dp),intent(inout) :: fofr(cplex*nfft)
 end subroutine fourdp
end interface

interface
 subroutine fourdp_6d(cplex,matrix,isign,MPI_enreg,nfft,ngfft,paral_kgb,tim_fourdp)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: isign
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  type(mpi_type),intent(inout) :: MPI_enreg
  integer,intent(in) :: ngfft(18)
  complex(gwpc),intent(inout) :: matrix(nfft,nfft)
 end subroutine fourdp_6d
end interface

interface
 subroutine fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&  
  &  kg_kin,kg_kout,mgfft,mpi_enreg,ndat,ngfft,npwin,npwout,n4,n5,n6,option,&  
  &  paral_kgb,tim_fourwf,weight_r,weight_i,&  
  &  use_gpu_cuda,use_ndo,fofginb) ! Optional arguments
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mgfft
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: ndat
  integer,intent(in) :: npwin
  integer,intent(in) :: npwout
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourwf
  integer,intent(in),optional :: use_gpu_cuda
  integer,intent(in),optional :: use_ndo
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: weight_i
  real(dp),intent(in) :: weight_r
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: denpot(cplex*n4,n5,n6)
  real(dp),intent(inout) :: fofgin(2,npwin*ndat)
  real(dp),intent(inout),optional :: fofginb(:,:)
  real(dp),intent(out) :: fofgout(2,npwout*ndat)
  real(dp),intent(inout) :: fofr(2,n4,n5,n6*ndat)
  integer,intent(in) :: gboundin(2*mgfft+8,2)
  integer,intent(in) :: gboundout(2*mgfft+8,2)
  integer,intent(in) :: kg_kin(3,npwin)
  integer,intent(in) :: kg_kout(3,npwout)
 end subroutine fourwf
end interface

interface
 subroutine indfftrisc(gbound,indpw_k,kg_k,mgfft,ngb,ngfft,npw_k)
  implicit none
  integer,intent(in) :: mgfft
  integer,intent(out) :: ngb
  integer,intent(in) :: npw_k
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: gbound(2*mgfft+4)
  integer,intent(out) :: indpw_k(4,npw_k)
  integer,intent(in) :: kg_k(3,npw_k)
 end subroutine indfftrisc
end interface

interface
 subroutine indgrid(coatofin,fintocoa,nfftc,nfftf,ngfftc,ngfftf)
  implicit none
  integer,intent(in) :: nfftc
  integer,intent(in) :: nfftf
  integer,intent(in) :: ngfftc(18)
  integer,intent(in) :: ngfftf(18)
  integer,intent(out) :: coatofin(nfftc)
  integer,intent(out) :: fintocoa(nfftf)
 end subroutine indgrid
end interface

interface
 subroutine kgindex(indpw_k,kg_k,mask,mpi_enreg,ngfft,npw_k)
  use defs_abitypes
  implicit none
  integer,intent(in) :: npw_k
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  integer,intent(out) :: indpw_k(npw_k)
  integer,intent(in) :: kg_k(3,npw_k)
  logical,intent(out) :: mask(npw_k)
 end subroutine kgindex
end interface

interface
 subroutine sg_fft(fftcache,nd1,nd2,nd3,n1,n2,n3,arr,ftarr,ris)
  use defs_basis
  implicit none
  integer,intent(in) :: fftcache
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: nd1
  integer,intent(in) :: nd2
  integer,intent(in) :: nd3
  real(dp),intent(in) :: ris
  real(dp),intent(inout) :: arr(2,nd1,nd2,nd3)
  real(dp),intent(out) :: ftarr(2,nd1,nd2,nd3)
 end subroutine sg_fft
end interface

interface
 subroutine sg_fftpad(fftcache,mgfft,nd1,nd2,nd3,n1,n2,n3,arr,ftarr,ris,gbound)
  use defs_basis
  implicit none
  integer,intent(in) :: fftcache
  integer,intent(in) :: mgfft
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: nd1
  integer,intent(in) :: nd2
  integer,intent(in) :: nd3
  real(dp),intent(in) :: ris
  real(dp),intent(inout) :: arr(2,nd1,nd2,nd3)
  real(dp),intent(out) :: ftarr(2,nd1,nd2,nd3)
  integer,intent(in) :: gbound(2*mgfft+8,2)
 end subroutine sg_fftpad
end interface

interface
 subroutine sg_fftpx(fftcache,mfac,mg,mgfft,nd1,nd2,nd3,n2,n3,&  
  &  z,zbr,trig,aft,now,bef,ris,ind,ic,gbound)
  use defs_basis
  implicit none
  integer,intent(in) :: fftcache
  integer,intent(in) :: ic
  integer,intent(in) :: mfac
  integer,intent(in) :: mg
  integer,intent(in) :: mgfft
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: nd1
  integer,intent(in) :: nd2
  integer,intent(in) :: nd3
  real(dp),intent(in) :: ris
  integer,intent(in) :: aft(mfac)
  integer,intent(in) :: bef(mfac)
  integer,intent(in) :: gbound(2*mgfft+4)
  integer,intent(in) :: ind(mg)
  integer,intent(in) :: now(mfac)
  real(dp),intent(in) :: trig(2,mg)
  real(dp),intent(inout) :: z(2,nd1,nd2,nd3)
  real(dp),intent(out) :: zbr(2,nd1,nd2,nd3)
 end subroutine sg_fftpx
end interface

interface
 subroutine sg_fftrisc(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,kg_kin,kg_kout,&  
  &  mgfft,ngfft,npwin,npwout,n4,n5,n6,option,weight)
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mgfft
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: npwin
  integer,intent(in) :: npwout
  integer,intent(in) :: option
  real(dp),intent(in) :: weight
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: denpot(cplex*n4,n5,n6)
  real(dp),intent(in) :: fofgin(2,npwin)
  real(dp),intent(out) :: fofgout(2,npwout)
  real(dp),intent(inout) :: fofr(2,n4,n5,n6)
  integer,intent(in) :: gboundin(2*mgfft+8,2)
  integer,intent(in) :: gboundout(2*mgfft+8,2)
  integer,intent(in) :: kg_kin(3,npwin)
  integer,intent(in) :: kg_kout(3,npwout)
 end subroutine sg_fftrisc
end interface

interface
 subroutine sg_fftrisc_2(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,&  
  &  istwf_k,kg_kin,kg_kout,&  
  &  mgfft,ngfft,npwin,npwout,n4,n5,n6,option,weight,&  
  &  luse_ndo,fofgin_p)
  use defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mgfft
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: npwin
  integer,intent(in) :: npwout
  integer,intent(in) :: option
  logical,intent(in),optional :: luse_ndo
  real(dp),intent(in) :: weight
  integer,intent(in) :: ngfft(18)
  real(dp),intent(inout) :: denpot(cplex*n4,n5,n6)
  real(dp),intent(in) :: fofgin(2,npwin)
  real(dp),intent(in),optional :: fofgin_p(:,:)
  real(dp),intent(out) :: fofgout(2,npwout)
  real(dp),intent(inout) :: fofr(2,n4,n5,n6)
  integer,intent(in) :: gboundin(2*mgfft+8,2)
  integer,intent(in) :: gboundout(2*mgfft+8,2)
  integer,intent(in) :: kg_kin(3,npwin)
  integer,intent(in) :: kg_kout(3,npwout)
 end subroutine sg_fftrisc_2
end interface

interface
 subroutine sg_fourwf(cplex,denpot,fftalgc,fofgin,fofgout,fofr,&  
  &  gboundin,gboundout,istwf_k,kg_kin,kg_kout,mgfft,mpi_enreg,n1,n2,n3,&  
  &  npwin,npwout,n4,n5,n6,option,paral_kgb,weight_r,weight_i)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: fftalgc
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mgfft
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: npwin
  integer,intent(in) :: npwout
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: weight_i
  real(dp),intent(in) :: weight_r
  real(dp),intent(inout) :: denpot(cplex*n4,n5,n6)
  real(dp),intent(inout) :: fofgin(2,npwin)
  real(dp),intent(out) :: fofgout(2,npwout)
  real(dp),intent(inout) :: fofr(2,n4,n5,n6)
  integer,intent(in) :: gboundin(2*mgfft+8,2)
  integer,intent(in) :: gboundout(2*mgfft+8,2)
  integer,intent(in) :: kg_kin(3,npwin)
  integer,intent(in) :: kg_kout(3,npwout)
 end subroutine sg_fourwf
end interface

interface
 subroutine sphere(cg,ndat,npw,cfft,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,mpi_enreg,shiftg,symm,xnorm)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iflag
  integer,intent(in) :: istwf_k
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: ndat
  integer,intent(in) :: npw
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: xnorm
  integer,intent(in) :: shiftg(3)
  integer,intent(in) :: symm(3,3)
  real(dp),intent(inout) :: cfft(2,n4,n5,n6*ndat)
  real(dp),intent(inout) :: cg(2,npw*ndat)
  integer,intent(in) :: kg_k(3,npw)
 end subroutine sphere
end interface

interface
 subroutine sphere_fft(cg,ndat,npw,cfft,n1,n2,n3,n4,n5,kg_k,&  
  &  mpi_enreg,nd2proc)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: nd2proc
  integer,intent(in) :: ndat
  integer,intent(in) :: npw
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cfft(2,n4,n5,nd2proc*ndat)
  real(dp),intent(inout) :: cg(2,npw*ndat)
  integer,intent(in) :: kg_k(3,npw)
 end subroutine sphere_fft
end interface

interface
 subroutine sphere_fft1(cg,ndat,npw,cfft,n1,n2,n3,n4,n5,n6,kg_k,&  
  &  mpi_enreg,nd2proc)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: n4
  integer,intent(in) :: n5
  integer,intent(in) :: n6
  integer,intent(in) :: nd2proc
  integer,intent(in) :: ndat
  integer,intent(in) :: npw
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(inout) :: cfft(2,n4,n5,n6)
  real(dp),intent(inout) :: cg(2,npw*ndat)
  integer,intent(in) :: kg_k(3,npw)
 end subroutine sphere_fft1
end interface

interface
 subroutine sphereboundary(gbound,istwf_k,kg_k,mgfft,npw)
  implicit none
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mgfft
  integer,intent(in) :: npw
  integer,intent(out) :: gbound(2*mgfft+8,2)
  integer,intent(in) :: kg_k(3,npw)
 end subroutine sphereboundary
end interface

interface
 subroutine zerosym(array,cplex,mpi_enreg,n1,n2,n3,ig1,ig2,ig3)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,optional,intent(in) :: ig1
  integer,optional,intent(in) :: ig2
  integer,optional,intent(in) :: ig3
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  type(mpi_type) :: mpi_enreg
  real(dp),intent(inout) :: array(cplex,n1*n2*n3)
 end subroutine zerosym
end interface

end module interfaces_53_ffts
!!***
