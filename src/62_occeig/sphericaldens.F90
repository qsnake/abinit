!{\src2tex{textfont=tt}}
!!****f* ABINIT/sphericaldens
!! NAME
!! sphericaldens
!!
!! FUNCTION
!! Compute the convolution of a function with 
!! the unity constant function over a sphere of radius rmax .
!! The function is to be given in reciprocal space,
!! the resulting function is also given in reciprocal space.
!! The routine needs the norm of the reciprocal space vectors.
!!
!! The resulting function in reciprocal space can give the
!! integral of the density in any sphere of that radius, centered
!! on any point, by a simple scalar product.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT group (MV,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  fofg(2,nfft)=initial function, in reciprocal space
!!  gnorm(nfft)=norm of the reciprocal space vectors
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  rmax=radius of the sphere
!!
!! OUTPUT
!!  sphfofg(2,nfft)=convoluted function, in reciprocal space
!!
!! PARENTS
!!      dens_in_sph
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine sphericaldens(fofg,gnorm,nfft,rmax,sphfofg)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sphericaldens'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
 real(dp),intent(in) :: rmax
!arrays
 real(dp),intent(in) :: fofg(2,nfft),gnorm(nfft)
 real(dp),intent(out) :: sphfofg(2,nfft)

!Local variables-------------------------------
!scalars
 integer :: ifft
 real(dp) :: factor,int0yy,rmax_2pi,yy
 character(len=500) :: message

! *************************************************************************

 if(nfft<1)then
   write(message,'(a,a,a,a,a,a,i6)') ch10,&
&   ' sphericaldens: BUG -',ch10,&
&   '  The argument nfft should be a positive number,',ch10,&
&   '  however, nfft=',nfft
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 rmax_2pi=two_pi*rmax
 factor=four_pi/(two_pi)**3
 do ifft=1,nfft
   if(abs(gnorm(ifft)) < tol12)then
     sphfofg(1,ifft)=fofg(1,ifft)*four_pi*third*rmax**3
     sphfofg(2,ifft)=fofg(2,ifft)*four_pi*third*rmax**3
   else
     yy=gnorm(ifft)*rmax_2pi
     int0yy=factor*(sin(yy)-yy*cos(yy))/(gnorm(ifft)**3)
     sphfofg(1,ifft)=fofg(1,ifft)*int0yy
     sphfofg(2,ifft)=fofg(2,ifft)*int0yy
   end if
 end do

end subroutine sphericaldens
!!***
