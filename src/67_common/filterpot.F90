!{\src2tex{textfont=tt}}
!!****f* ABINIT/filterpot
!! NAME
!! filterpot
!!
!! FUNCTION
!! Given a potential, v(r), filter it, so as to let
!! only the active Fourier components, inside gsqcut.
!! When cplex=1, assume q=(0 0 0), and vhartr will be REAL
!! When cplex=2, q must be taken into account, and vhartr will be COMPLEX
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! *Modified code to avoid if statements inside loops to skip G=0.
!!  Replaced if statement on G^2>gsqcut to skip G s outside where
!!  rho(G) should be 0.  Effect is negligible but gsqcut should be
!!  used to be strictly consistent with usage elsewhere in code.
!! *The speed-up is provided by doing a few precomputations outside
!!  the inner loop. One variable size array is needed for this (gq).
!!
!! INPUTS
!!  cplex= if 1, vhartr is REAL, if 2, vhartr is COMPLEX
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!                                 (gsqcut=(boxcut**2)*ecut/(2._dp*(Pi**2))
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nspden=number of spin-components of the potential
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECT
!!  vpot(cplex*nfft,nspden)=potential in real space, either REAL or COMPLEX
!!
!! PARENTS
!!
!! CHILDREN
!!      fourdp,leave_new,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine filterpot(cplex,gmet,gsqcut,nfft,ngfft,nspden,paral_kgb,qphon,vpot)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'filterpot'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,paral_kgb
 real(dp),intent(in) :: gsqcut
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),qphon(3)
 real(dp),intent(inout) :: vpot(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i23,i3,ifft,ig,ii,ii1,ing,ispden,n1,n2,n3,qeq0
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: cutoff,gqg2p3,gqgm12,gqgm13,gqgm23,gs,gs2,gs3
 character(len=500) :: message
 type(MPI_type) :: mpi_enreg
!arrays
 integer :: id(3)
 real(dp) :: tsec(2)
 real(dp),allocatable :: gq(:,:),wkcmpx(:,:),work(:)

! *************************************************************************
!
!DEBUG
!write(std_out,*)' filterpot : enter '
!write(std_out,*)' cplex,nfft,ngfft',cplex,nfft,ngfft
!write(std_out,*)' gsqcut=',gsqcut
!write(std_out,*)' qphon=',qphon
!write(std_out,*)' gmet=',gmet
!ENDDEBUG

!Check that cplex has an allowed value
 if(cplex/=1 .and. cplex/=2)then
   write(message, '(a,a,a,a,i3,a,a)' )ch10,&
&   ' filterpot : BUG -',ch10,&
&   '  From the calling routine, cplex=',cplex,ch10,&
&   '  but the only value allowed are 1 and 2.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!Initialize a few quantities
 cutoff=gsqcut*tolfix
!DEBUG
!cutoff=zero
!ENDDEBUG
!This is to allow q=0
 qeq0=0
 if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) qeq0=1

!If cplex=1 then qphon should be 0 0 0
 if (cplex==1.and. qeq0/=1) then
   write(message, '(a,a,a,a,3e12.4,a,a)' ) ch10,&
&   ' filterpot : BUG -',ch10,&
&   '  cplex=1 but qphon=',qphon,ch10,&
&   '  qphon should be 0 0 0.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 mpi_enreg%me_fft=0
 mpi_enreg%nproc_fft=1

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig+qphon(ii)
   end do
 end do

 ABI_ALLOCATE(wkcmpx,(2,nfft))
 ABI_ALLOCATE(work,(nfft))

!Loop on the spin components
 do ispden=1,nspden

!  Obtain vpot(G) in wkcmpx from input vpot(r)
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(ispden,nfft,vpot,work)
   do ifft=1,cplex*nfft
     work(ifft)=vpot(ifft,ispden)
   end do
!  $OMP END PARALLEL DO

   call timab(82,1,tsec)
   call fourdp(1,wkcmpx,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   call timab(82,2,tsec)

!  Triple loop on each dimension
   do i3=1,n3

!    Precompute some products that do not depend on i2 and i1
     gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
     gqgm23=gq(3,i3)*gmet(2,3)*2
     gqgm13=gq(3,i3)*gmet(1,3)*2

     do i2=1,n2
       gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
       gqgm12=gq(2,i2)*gmet(1,2)*2
       gqg2p3=gqgm13+gqgm12

       i23=n1*((i2-1)+n2*(i3-1))
!      Do the test that eliminates the Gamma point outside
!      of the inner loop
       ii1=1
!      if(i23==0 .and. qeq0==1)then
!      ii1=2
!      wkcmpx(re,1+i23)=0.0_dp
!      wkcmpx(im,1+i23)=0.0_dp
!      end if

!      Final inner loop on the first dimension
!      (note the lower limit)
       do i1=ii1,n1
         gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
         ii=i1+i23
         if(gs>cutoff)then
           wkcmpx(re,ii)=zero
           wkcmpx(im,ii)=zero
         end if
!        End loop on i1
       end do

!      End loop on i2
     end do

!    End loop on i3
   end do

   ABI_DEALLOCATE(gq)

!  DEBUG
!  write(std_out,*)' filterpot : before fourdp'
!  write(std_out,*)' cplex,nfft,ngfft',cplex,nfft,ngfft
!  write(std_out,*)' maxval work1=',maxval(abs(work1))
!  ENDDEBUG

!  Fourier Transform Vpot
   call fourdp(cplex,wkcmpx,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)

!  Obtain filtered vpot(G) from wkcmpx
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(ispden,nfft,vpot,work)
   do ifft=1,cplex*nfft
     vpot(ifft,ispden)=work(ifft)
   end do
!  $OMP END PARALLEL DO

 end do ! ispden=1,nspden

!DEBUG
!write(std_out,*)' filterpot : after fourdp'
!write(std_out,*)' cplex,nfft,ngfft',cplex,nfft,ngfft
!write(std_out,*)' maxval work1=',maxval(abs(work1))
!write(std_out,*)' maxval vhartr=',maxval(abs(vhartr))
!write(std_out,*)' filterpot : set vhartr to zero'
!vpot(:,:)=0.0_dp
!ENDDEBUG

 ABI_DEALLOCATE(wkcmpx)
 ABI_DEALLOCATE(work)

end subroutine filterpot
!!***
