!{\src2tex{textfont=tt}}
!!****f* ABINIT/kxc_eok
!! NAME
!! kxc_eok
!!
!! FUNCTION
!!  Compute the linear (ixceok = 1) or non-linear (ixceok = 2)
!!  energy optimized kernel of Dobson and Wang, in reciprocal
!!  space, on the FFT grid.
!!  [see J. Dobson and J. Wang, Phys. Rev. B 62, 10038 (2000)].
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, MF, XG, GMR, LSI, YMN).
!! This file is distributed under the terms of the
!! GNU General Public License,see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ixceok = 1 linear energy optimized kernel.
!!         = 2 non-linear energy optimized kernel.
!!  mpi_enreg=informations about MPI parallelization
!!  nfft = number of fft grid points.
!!  ngfft(1:3) = integer fft box dimensions, see getng for ngfft(4:8).
!!  nspden = number of spin-density components.
!!  rhor(nfft,nspden) = electron density in real space in electrons/bohr**3
!!   (total in first half and spin-up in second half if nspden = 2).
!!  rhocut = cut-off density for the local kernels (ALDA, EOK),
!!           relative to max(rhor(:,:)).
!! OUTPUT
!!  kxcg(2,nfft,2*nspden-1) = the EOK kernel in reciprocal space, on the FFT grid.
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!!
!! PARENTS
!!      acfd_dyson,xcacfd
!!
!! CHILDREN
!!      fourdp,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine kxc_eok(ixceok,kxcg,mpi_enreg,nfft,ngfft,nspden,paral_kgb,rhor,rhocut)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kxc_eok'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: ixceok,nfft,nspden,paral_kgb
 real(dp),intent(in) :: rhocut
 type(MPI_type) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,2*nspden-1)
 real(dp),intent(out) :: kxcg(2,nfft,2*nspden-1)

!Local variables -------------------------------------------------------
!Maximum value allowed for rs.
!scalars
 integer :: ifft,ikxc,ncut,nkxc,nlop,tim_fourdp
 real(dp),parameter :: rslim=50._dp
 real(dp) :: a2,a3,a4,rho,rhocuttot,rhomin,rs
 character(len=500) :: message
!arrays
 real(dp),allocatable :: kxcr(:,:)

!***********************************************************************

!Check input parameters.

 if (nspden > 1) then
   write (message,'(4a)') ch10,&
&   ' kxc_eok: ERROR - ',ch10,&
&   '  kxc_eok does not work yet for nspden > 1.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Values of a2, a3 and a4 for case 1
 a2 = 0.0_dp
 a3 = 0.0_dp
 a4 = 0.0_dp

 select case (ixceok)
   case (1)
     a2 = -0.51887_dp
     a3 =  4.9359d-03
     a4 = -5.9603d-05
   case (2)
     a2 = -0.50044_dp
     a3 =  4.9653d-03
     a4 = -3.3660d-05
     case default
     write (message,'(4a)') ch10,&
&     ' kxc_eok: ERROR - ',ch10,&
&     '  ixceok /= 1 (linear EOK) or 2 (non-linear EOK).'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
 end select

!Allocate memory.

 nkxc = 2*nspden-1

 ABI_ALLOCATE(kxcr,(nfft,nkxc))

!Calculate the energy optimized kernel in real space.

 nlop = 0

 rhomin = rhocut*maxval(rhor(:,:))

 ncut = 0
 rhocuttot = 0._dp

 do ifft = 1,nfft

   rho = rhor(ifft,1)

   if (rho < rhomin) then
     ncut = ncut+1
     rhocuttot = rhocuttot+rho
     rho = rhomin
   end if

   rs = (3._dp/(4._dp*pi*rho))**(1._dp/3._dp)

   if (rs > rslim) then
     rs = rslim
     nlop = nlop+1
   end if

   kxcr(ifft,1) = a2*rs**2+a3*rs**3+a4*rs**4

 end do

 if (ncut > 0) then
   write (message,'(a,es10.3,3a,i1,a,i6,a,f6.3,3a,f6.3,a)') &
&   ' kxc_eok: WARNING - rhocut = ',rhocut,'.',ch10,&
&   '  For isp = ',1,' the density was cut-off at ',ncut,' (',&
&   100._dp*float(ncut)/float(ifft),'%) grid points.',ch10,&
&   '  These points account for ',100._dp*rhocuttot/sum(rhor(:,1)),'% of the total density.'
   call wrtout(std_out,message,'COLL')
 end if

 if (nlop > 0) then
   write (message,'(3a,f6.2,a,i6,a,f6.3,a)') &
&   ' kxc_eok: WARNING - ',ch10,&
&   '  rs still exceeds ',rslim,' Bohr at ',nlop,' (',&
&   100._dp*float(nlop)/float(ifft),'%) grid points (after cut-off).'
   call wrtout(std_out,message,'COLL')
 end if

!Calculate the Fourier transform of the energy optimized kernel.
 tim_fourdp=0
 do ikxc = 1,nkxc
   call fourdp(1,kxcg(:,:,ikxc),kxcr(:,ikxc),-1,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp)
 end do

!Free memory.

 ABI_DEALLOCATE(kxcr)

end subroutine kxc_eok

!!***
