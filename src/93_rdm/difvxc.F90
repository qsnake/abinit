!{\src2tex{textfont=tt}}
!!****f* ABINIT/difvxc
!! NAME
!! difvxc
!!
!! FUNCTION
!! Given the charge density of a homogenous e-gas, this calculates
!! the differential of the Ceperley Alder exchange correlation potential
!! wrt rho (includes relativistic effects see VXCCA).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! INPUTS
!!  rho=density in real space, at a specific grid point
!!
!! OUTPUT
!!  difvxc=differential of the Ceperley Alder exchange correlation potential wrt rho
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function difvxc(rho)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'difvxc'
 use interfaces_93_rdm, except_this_one => difvxc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: difvxc
 real(dp),intent(in) :: rho

!Local variables ------------------------------
!scalars
 integer,save :: iprint=0

!************************************************************************
 if(rho<1.0e-10) then
   if(iprint/=1) then
     write(std_out,'(/,/,a,e14.6,/,a)')&
&     '***Warning: RHO very small in DIFVXC:',rho,&
&     '  this message will not be printed again'
     iprint=1
   end if
   difvxc=0.0
 else
   difvxc=diffvc(rho)+diffvx(rho)
 end if

 end function difvxc
!!***

!!****f* ABINIT/diffvx
!! NAME
!! diffvx
!!
!! FUNCTION
!! Calculate the differential of the exchange potential wrt Rho for
!! a homogenous e-gas of density RHO. This is done by differenting the
!! relativistic term (see VXCCA).
!!
!! PARENTS
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  rho=density in real space, at a specific grid point
!!
!! OUTPUT
!!  diffvx=differential of the exchange potential wrt Rho
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function diffvx(rho)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'diffvx'
 use interfaces_93_rdm, except_this_one => diffvx
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: diffvx
 real(dp),intent(in) :: rho

!Local variables ------------------------------
!scalars
 real(dp) :: rs

!************************************************************************
 rs=(3.0/(4.0*pi*rho))**third
 diffvx=(0.610887057/(rs*rs))*rel(rho)+vxnr(rho)*difrel(rho)
 diffvx=diffvx*(-4.0*pi/9.0)*(rs**4)

 return
 end function diffvx
!!***


!!****f* ABINIT/diffvc
!! NAME
!! diffvc
!!
!! FUNCTION
!! Given the charge density of a homogenous e-gas, this calculates
!! the differential of the correlation potential wrt rho.
!! It does this by evaluating the differential of the Ceperley-Alder
!! expression given by PERDEW & ZUNGER (PRB,Vol 23,No 10,P.5048,15 May 1981)
!! Gamma,Beta1,Beta2,A,B,C,D are parameters for the fit, in atomic units.
!!
!! PARENTS
!!
!! INPUTS
!!  rho=density in real space, at one point
!!
!! OUTPUT
!!  function diffvc=differential of the correlation potential wrt rho
!!
!! TODO
!!  Use the other abinit routine
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function diffvc(rho)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'diffvc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: diffvc
 real(dp),intent(in) :: rho

!Local variables-------------------------------
!scalars
 real(dp),parameter :: a=0.0311,b=-0.0480,beta1=1.0529,beta2=0.3334,c=0.0020
 real(dp),parameter :: d=-0.0116,gamma=-0.1423
 real(dp) :: rs,stor1,stor2,stor3

! *************************************************************************

!calculate Rs from the density
 rs=(3.0/(4.0*pi*rho))**third
!two different calculations depending if Rs<>1
 if(rs>1.0) then
   stor1=(1.0+beta1*sqrt(rs)+beta2*rs)**(-3.0)
   stor2=-0.41666667*beta1*(rs**(-0.5))-0.5833333*(beta1**2)-0.66666667*beta2
   stor3=-1.75*beta1*beta2*sqrt(rs)-1.3333333*rs*(beta2**2)
   diffvc=gamma*stor1*(stor2+stor3)
 else
   diffvc=a/rs+0.66666667*(c*log(rs)+d)+0.33333333*c
 end if
 diffvc=diffvc*(-4.0*pi/9.0)*(rs**4)
 return
 end function diffvc
!!***


!!****f* ABINIT/difrel
!! NAME
!! difrel
!!
!! FUNCTION
!! Calculate the differential of the relativistic to Vx.
!!
!! PARENTS
!!
!! INPUTS
!!  rho=density in real space, at a specific grid point
!!
!! OUTPUT
!!  function difrel=differential of the relativistic to Vx
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function difrel(rho)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'difrel'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: difrel
 real(dp),intent(in) :: rho

!Local variables-------------------------------
!scalars
 real(dp) :: b,bb,bb1,rs

! *************************************************************************

 rs=(3.0/(4.0*pi*rho))**third
 b=0.0140/rs
 bb=b*b
 bb1=1.0+bb
 difrel=(1.5/(b*bb1))-1.5*log(b+sqrt(bb1))*(1.0+2.0*bb)*(bb1**(-1.5))/bb
 difrel=difrel*(-0.0140)/(rs*rs)

 return
 end function difrel
!!***

!!****f* ABINIT/vxnr
!! NAME
!! vxnr
!!
!! FUNCTION
!! Calculate the exchange potential without any relativistic corrections.
!!
!! INPUTS
!!  rho=density at a specific point
!!
!! OUTPUT
!!  function vxnr=exchange potential
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function vxnr(rho)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vxnr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: vxnr
 real(dp),intent(in) :: rho

!Local variables-------------------------------

! *************************************************************************
 vxnr=-(3*rho/pi)**third

 end function vxnr
!!***


!!****f* ABINIT/vxcca
!! NAME
!! vxcca
!!
!! FUNCTION
!! Calculate the exchange correlation energy using Ceperley-Alder
!! and including relativistic effects (PRB 26 P.4199, 1982)
!! rho should be in e/au**3; Vxc is returned in Hartree
!!
!! INPUTS
!!  rho=density at a specific point
!!
!! OUTPUT
!!  function vxcca=CA potential
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function vxcca(rho)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vxcca'
 use interfaces_93_rdm, except_this_one => vxcca
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: vxcca
 real(dp),intent(in) :: rho

!Local variables ------------------------------

!************************************************************************
 if(rho<1.0e-10) then
   vxcca=zero
 else
   vxcca=vxjas(rho)+vcjas(rho)
 end if

 end function vxcca
!!***


!!****f* ABINIT/vxjas
!! NAME
!! vxjas
!!
!! FUNCTION
!! Calculate the exchange potential with relativistic corrections.
!!
!! INPUTS
!!  rho=density at a specific point
!!
!! OUTPUT
!!  function vxjas=exchange potential
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function vxjas(rho)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vxjas'
 use interfaces_93_rdm, except_this_one => vxjas
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: vxjas
 real(dp),intent(in) :: rho

!Local variables ------------------------------

!************************************************************************
 vxjas=-((3*rho/pi)**0.3333333333)*rel(rho)

 end function vxjas
!!***


!!****f* ABINIT/vcjas
!!
!! NAME
!! vcjas
!!
!! FUNCTION
!! Calculate the correlation potential for a homogenous e-gas,
!! with charge density RHO. Uses the PERDEW & ZUNGER calculations
!! (see diffwc).
!!
!! PARENTS
!!
!! INPUTS
!!  rho=density in real space, at a particular point
!!
!! OUTPUT
!!  vcjas=correlation potential
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function vcjas(rho)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vcjas'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: vcjas
 real(dp),intent(in) :: rho

!Local variables-------------------------------
!scalars
 real(dp),parameter :: a=0.0311,b=-0.0480,beta1=1.0529,beta2=0.3334,c=0.0020
 real(dp),parameter :: d=-0.0116,gamma=-0.1423
 real(dp) :: rs,stor1,stor2

! *************************************************************************

!first calculate Rs from RHO
 rs=(3.0/(4.0*pi*rho))**third
!need to do different calculations depending on Rs<>1
 if(rs>1.0) then
   stor1=(1.0+1.16666667*beta1*sqrt(rs)+1.3333333*beta2*rs)
   stor2=gamma*(1.0+beta1*sqrt(rs)+beta2*rs)**(-2.0)
   vcjas=stor1*stor2
 else
   stor1=a*log(rs)+(b-third*a)
   stor2=0.6666667*c*rs*log(rs)+0.33333333*(2*d-c)*rs
   vcjas=stor1+stor2
 end if

 end function vcjas
!!***


!!****f* ABINIT/rel
!! NAME
!! rel
!!
!! FUNCTION
!! Calculate the relativistic correction to Vx.
!!
!! PARENTS
!!
!! INPUTS
!!  rho=density at one point in real space
!!
!! OUTPUT
!!  rel=relativistic factor
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function rel(rho)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rel'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: rel
 real(dp),intent(in) :: rho

!Local variables-------------------------------
!scalars
 real(dp) :: b,rs

! *************************************************************************

 rs=(3.0/(4.0*pi*rho))**third
 b=0.0140/rs
 rel=-0.5+1.5*log(b+sqrt(1.0+b*b))/(b*sqrt(1.0+b*b))

 end function rel
!!***


!!****f* ABINIT/ckxcldar
!! NAME
!! ckxcldar
!!
!! FUNCTION
!! Calculate exchange-correlation kernel Kxclda in real space
!!
!! INPUTS
!!  nr=number of points of FFT grid
!!  rho2(nr)=density in real space
!!
!! OUTPUT
!!  kxclda=exchange-correlation kernel in real space
!!
!! PARENTS
!!      difvxc
!!
!! CHILDREN
!!      ckxcldar,fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ckxcldar(nr,rho2,kxclda)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ckxcldar'
 use interfaces_93_rdm, except_this_one => ckxcldar
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nr
!arrays
 real(dp),intent(in) :: rho2(nr)
 complex,intent(out) :: kxclda(nr)

!Local variables ------------------------------
!scalars
 integer :: ir

!************************************************************************
 write(std_out,*) 'calculating kxclda'
 do ir=1,nr
   kxclda(ir)=difvxc(rho2(ir))
 end do
 write(std_out,*)

end subroutine ckxcldar
!!***


!!****f* ABINIT/ckxcldag
!! NAME
!! ckxcldag
!!
!! FUNCTION
!! Calculate exchange-correlation kernel Kxclda in G-space
!! put it into Kxclda in igfft form
!! (Needed only if GWGamma calculation - not operational ?!)
!!
!! INPUTS
!!  ngfft1,ngfft2,ngfft3=FFT grid dimensions
!!  nr=number of points of FFT grid
!!  rho2(nr)=density in real space
!!
!! OUTPUT
!!  kxclda(nr)=XC kernel in reciprocal space
!!
!! PARENTS
!!
!! CHILDREN
!!      ckxcldar,fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ckxcldag(ngfft1,ngfft2,ngfft3,nr,paral_kgb,rho2,kxclda)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ckxcldag'
 use interfaces_53_ffts
 use interfaces_93_rdm, except_this_one => ckxcldag
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ngfft1,ngfft2,ngfft3,nr,paral_kgb
!arrays
 real(dp),intent(in) :: rho2(nr)
 complex,intent(out) :: kxclda(nr)

!Local variables-------------------------------
!scalars
 integer :: tim_fourdp
 type(MPI_type) :: mpi_enreg
!arrays
 integer :: ngfft(18)
 real(dp),allocatable :: kxcldag(:,:),kxcldar(:,:)

! *************************************************************************

 ABI_ALLOCATE(kxcldar,(2,nr))
 ABI_ALLOCATE(kxcldag,(2,nr))

 call ckxcldar(nr,rho2,kxclda)
 kxcldar(1,:)=real(kxclda(:))
 kxcldar(2,:)=aimag(kxclda(:))

 write(std_out,*) 'fft kxclda(r)->kxclda(g)'
 write(std_out,*)

!call c3dfft(kxclda,ngfft1a,ngfft2,ngfft3,ngfft1,ngfft2,ngfft3,-1)
 ngfft(1)=ngfft1
 ngfft(2)=ngfft2
 ngfft(3)=ngfft3
 ngfft(4)=2*(ngfft(1)/2)+1
 ngfft(5)=2*(ngfft(2)/2)+1
 ngfft(6)=ngfft(3)
 ngfft(7)=200
 ngfft(8)=256
 ngfft(9)=0
 ngfft(10)=1
 ngfft(11)=0
 ngfft(12)=ngfft2
 ngfft(13)=ngfft3
 ngfft(14)=0
 tim_fourdp=3
 mpi_enreg%nproc_fft=1
 mpi_enreg%me_fft=0
 call fourdp(2,kxcldag,kxcldar,-1,mpi_enreg,nr,ngfft,paral_kgb,tim_fourdp)

!do ig=1, npwvec
!kxcg(ig)=Kxclda(igfft(ig,3,3,3))*fftfct
!end do

!renormalization after fft
!fftfct=1/real(ngfft1*ngfft2*ngfft3)
!kxclda(:)=kxclda(:)*fftfct
 kxclda(:)=kxcldag(1,:)+(0.0,1.0)*kxcldag(2,:)

 ABI_DEALLOCATE(kxcldar)
 ABI_DEALLOCATE(kxcldag)

end subroutine ckxcldag
!!***
