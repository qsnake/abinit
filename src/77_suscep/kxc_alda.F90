!{\src2tex{textfont=tt}}
!!****f* ABINIT/kxc_alda
!! NAME
!! kxc_alda
!!
!! FUNCTION
!! If option = 1:
!!  Compute the AL(S)DA kernel in reciprocal space, on the FFT grid.
!! If option = 2:
!!  Only computes the up-down channel of the AL(S)DA kernel, on the
!!  FFT grid, for use in the BPG kernel.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, MF, XG, GMR, LSI, YMN).
!! This file is distributed under the terms of the
!! GNU General Public License,see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  ixc = choice of exchange-correlation functional.
!!  mpi_enreg=informations about MPI parallelization
!!  nfft = number of fft grid points.
!!  ngfft(1:3) = integer fft box dimensions, see getng for ngfft(4:8).
!!  nspden = number of spin-density components.
!!  option = 1 compute the AL(S)DA kernel in reciprocal space.
!!         = 2 only computes the up-down channel of the AL(S)DA kernel,
!!             for use in the BPG kernel.
!!  rhor(nfft,nspden) = electron density in real space in electrons/bohr**3
!!   (total in first half and spin-up in second half if nspden = 2).
!!  rhocut = cut-off density for the local kernels (ALDA, EOK),
!!           relative to max(rhor(:,:)).
!!  rprimd(3,3) = dimensional primitive translations for real space in Bohr.
!!
!! OUTPUT
!!  kxcg(2,nfft,*) = the AL(S)DA kernel in reciprocal space, on the FFT grid
!!   (the third dimension is 2*nspden-1 if option = 1, and 1 if option = 2).
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!! Current restrictions are:
!!  a - Spin-polarized case not tested.
!!
!! PARENTS
!!      acfd_dyson
!!
!! CHILDREN
!!      dtsetcopy,dtsetfree,fourdp,leave_new,rhohxc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine kxc_alda(dtset,ixc,kxcg,mpi_enreg,nfft,ngfft,nspden,option,rhor,rhocut,rprimd)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kxc_alda'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_53_abiutil
 use interfaces_53_ffts
 use interfaces_56_xc
!End of the abilint section

 implicit none

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: ixc,nfft,nspden,option
 real(dp),intent(in) :: rhocut
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,2*nspden-1),rprimd(3,3)
 real(dp),intent(out) :: kxcg(2,nfft,*)

!Local variables -------------------------------------------------------
!No improved xc quadrature.
!No core correction.
!Dummy here.
!For debugging purposes (see tests below):
!integer :: i1,i2,i3,k1,n1,n2,n3
!real(dp) :: kx,rho,rhomax,ftest
!scalars
 integer :: ifft,ikxc,isp,n3xccc,ncut,nk3xc,nkxc,optionrhoxc,tim_fourdp
 real(dp),parameter :: gsqcut=1._dp
 real(dp) :: enxc,rhocuttot,rhomin,vxcavg
 character(len=500) :: message
 type(dataset_type) :: dtLocal
!arrays
 real(dp) :: strsxc(6)
 real(dp),allocatable :: dum(:),kxcr(:,:),rhog(:,:),rhorcut(:,:),vhartree(:)
 real(dp),allocatable :: vxc(:,:),xccc3d(:)

!***********************************************************************
!For debugging purposes (see tests below):

!ftest(i1,n1,k1) = 0._dp+1._dp*cos(k1*two_pi*float(i1)/float(n1))

!***********************************************************************

!Check input parameters.

 if (nspden > 2) then
   write (message,'(4a)') ch10,&
&   ' kxc_alda: ERROR - ',ch10,&
&   '  kxc_alda does not work yet for nspden > 2.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Allocate memory.

 ABI_ALLOCATE(rhorcut,(nfft,nspden))
 ABI_ALLOCATE(rhog,(2,nfft))
 ABI_ALLOCATE(vhartree,(nfft))
 ABI_ALLOCATE(vxc,(nfft,nspden))

!Copy the input variables from the current dataset to a temporary one
!to tune some parameters
 call dtsetCopy(dtLocal, dtset)
 dtLocal%intxc = 0
 dtLocal%ixc   = ixc


!to be adjusted for the call to rhohxc
 nk3xc=1

!Cut-off the density.

 rhorcut(:,:) = rhor(:,:)

 do isp = 1,nspden

   rhomin = maxval(rhorcut(:,isp))*rhocut

   ncut = 0
   rhocuttot = 0._dp

   do ifft = 1,nfft
     if (rhorcut(ifft,isp) < rhomin) then
       ncut = ncut+1
       rhocuttot = rhocuttot+rhorcut(ifft,isp)
       rhorcut(ifft,isp) = rhomin
     end if
   end do

   if (ncut > 0) then
     write (message,'(a,es10.3,3a,i1,a,i6,a,f6.3,3a,f6.3,a)') &
&     ' kxc_alda: WARNING - rhocut = ',rhocut,'.',ch10,&
&     '  For isp = ',isp,' the density was cut-off at ',ncut,' (',&
&     100._dp*float(ncut)/float(ifft),'%) grid points.',ch10,&
&     '  These points account for ',100._dp*rhocuttot/sum(rhor(:,isp)),'% of the total density.'
     call wrtout(std_out,message,'COLL')
   end if

 end do

!Calculate the AL(S)DA kernel in real space.

 rhog(:,:) = 0._dp !We do not need the Hartree potential.
 tim_fourdp=0

 if ((option == 1).or.((option == 2).and.(nspden == 2))) then

   nkxc = 2*nspden-1
   n3xccc=0
   ABI_ALLOCATE(kxcr,(nfft,nkxc))
   ABI_ALLOCATE(xccc3d,(n3xccc))

   optionrhoxc = 2 !See rhohxc.f


   call rhohxc(dtLocal,enxc,gsqcut,0,kxcr,mpi_enreg,nfft,ngfft,dum,0,dum,0,nkxc,nk3xc,nspden,n3xccc,&
&   optionrhoxc,rhog,rhorcut,rprimd,strsxc,1,vhartree,vxc,vxcavg,xccc3d)

!  DEBUG
!  fx for tests.
!  write (std_out,'(a)') ' kxc_alda: Using exchange-only kernel for tests.'
!  rhomin = minval(rhor(:,1))
!  rhomax = maxval(rhor(:,1))
!  write (std_out,'(a,es12.5,a,es12.5)') ' kxc_alda: rhomin = ',rhomin,' rhomax = ',rhomax
!  write (std_out,'(a)') ' kxc_alda: loping below 0.2*rhomax.'
!  kx = (3._dp/4._dp)*((3._dp/pi)**(1._dp/3._dp))
!  do ifft = 1,nfft
!  rho = rhor(ifft,1)
!  rho = max(rho,0.2_dp*rhomax)
!  kxcr(ifft,1) = -(4._dp/9._dp)*kx*(rho**(-2._dp/3._dp))
!  write (std_out,'(i4,a,es12.5)') ifft,': ',kxcr(ifft,1)
!  end do
!  write (std_out,'(a,es12.5)') 'kxcrmin: ',minval(kxcr(:,1))
!  write (std_out,'(a,es12.5)') 'kxcrmax: ',maxval(kxcr(:,1))
!  ENDDEBUG

!  DEBUG
!  test kernel.
!  write(std_out,'(a)') ' kxc_alda: Using test kernel for tests.'
!  n1 = ngfft(1) ; n2 = ngfft(2) ; n3 = ngfft(3)
!  do i3 = 0,n3-1
!  do i2 = 0,n2-1
!  do i1 = 0,n1-1
!  ifft = i1+n1*(i2+n2*i3)+1
!  kxcr(ifft,1) = ftest(i1,n1,1)*ftest(i2,n2,2)*ftest(i3,n3,3)
!  end do
!  end do
!  end do
!  ENDDEBUG

!  Calculate the Fourier transform of the AL(S)DA kernel.

   if (option == 1) then
     do ikxc = 1,nkxc
       call fourdp(1,kxcg(:,:,ikxc),kxcr(:,ikxc),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,tim_fourdp)
     end do
   else
     call fourdp(1,kxcg(:,:,1),kxcr(:,2),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,tim_fourdp)
   end if

 else if ((option == 2).and.(nspden == 1)) then

   nkxc = 2
   n3xccc=0
   ABI_ALLOCATE(kxcr,(nfft,nkxc))
   ABI_ALLOCATE(xccc3d,(n3xccc))

   optionrhoxc = -2 !See rhohxc.f

   call rhohxc(dtLocal,enxc,gsqcut,0,kxcr,mpi_enreg,nfft,ngfft,dum,0,dum,0,nkxc,nk3xc,nspden,n3xccc,&
&   optionrhoxc,rhog,rhorcut,rprimd,strsxc,1,vhartree,vxc,vxcavg,xccc3d)

   kxcr(:,2) = 0.5_dp*kxcr(:,2)

!  Calculate the Fourier transform of the up-down channel of the AL(S)DA kernel.

   call fourdp(1,kxcg(:,:,1),kxcr(:,2),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,tim_fourdp)

 else

   write (message,'(4a,i10,a)') ch10,&
&   ' kxc_alda: ERROR - ',ch10,&
&   '  Invalid option = ',option,'.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')

 end if

!DEBUG
!write(std_out,*) ' kxc_alda:  Exc  = ',enxc
!write(std_out,*) ' kxc_alda: <Vxc> = ',vxcavg
!ENDDEBUG

!Free memory.
 call dtsetFree(dtLocal)
 ABI_DEALLOCATE(rhorcut)
 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(vhartree)
 ABI_DEALLOCATE(vxc)
 ABI_DEALLOCATE(kxcr)
 ABI_DEALLOCATE(xccc3d)

end subroutine kxc_alda

!!***
