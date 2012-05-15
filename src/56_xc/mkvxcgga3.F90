!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkvxcgga3
!! NAME
!! mkvxcgga3
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! in case of GGA functionals
!! Use the exchange-correlation kernel.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2012 ABINIT group (XG, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gprimd(3,3)=dimensional primitive translations in reciprocal space (bohr^-1)
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see rhohxc.f)
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor1tmp(cplex*nfft,2)=array for first-order electron spin-density
!!   in electrons/bohr**3 (second index corresponds to spin-up and spin-down)
!!
!! OUTPUT
!!  vxc1(cplex*nfft,nspden)=change in exchange-correlation potential
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  For the time being, a rather crude coding, to be optimized ...
!!
!! PARENTS
!!      mkvxc3
!!
!! CHILDREN
!!      xcden,xcpot
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkvxcgga3(cplex,gprimd,kxc,mpi_enreg,nfft,ngfft,nkxc,&
& nspden,paral_kgb,qphon,rhor1tmp,vxc1)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkvxcgga3'
 use interfaces_56_xc, except_this_one => mkvxcgga3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nkxc,nspden,paral_kgb
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),kxc(nfft,nkxc),qphon(3)
 real(dp),intent(in) :: rhor1tmp(cplex*nfft,2)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ir,ishift,ispden,ngrad,nspdentmp,nspgrad
 real(dp) :: coeff_grho_corr,coeff_grho_dn,coeff_grho_up,coeffim_grho_corr
 real(dp) :: coeffim_grho_dn,coeffim_grho_up,gradrho_gradrho1
 real(dp) :: gradrho_gradrho1_dn,gradrho_gradrho1_up,gradrho_gradrho1im
 real(dp) :: gradrho_gradrho1im_dn,gradrho_gradrho1im_up
! character(len=500) :: message
!arrays
 real(dp) :: r0(3),r0_dn(3),r0_up(3),r1(3),r1_dn(3),r1_up(3)
 real(dp) :: r1im(3),r1im_dn(3),r1im_up(3)
 real(dp),allocatable :: dnexcdn(:,:),rho1now(:,:,:),rhortmp(:,:,:)
 real(dp),allocatable :: vxc1tmp(:,:)
 logical :: mgga
!no_abirules
!real(dp) :: ec(7),ex(8)

! *************************************************************************

!DEBUG
!write(std_out,*)' mkggavxc3 : enter '
!write(std_out,*)' mkggavxc3 : cplex,nfft,nkxc,nspden',cplex,nfft,nkxc,nspden
!if(cplex==2)then
!write(std_out,*)' mkvxcgga3 : not yet cplex=2, sorry'
!stop
!end if
!stop
!ENDDEBUG

!metaGGA contributions are not taken into account here
 mgga=.false.

!Treat all cases in a generic way (to be optimized !!)
 nspdentmp=2

!call filterpot(paral_kgb,cplex,gmet,gsqcut,nfft,ngfft,2,qphon,rhor1tmp)

!Compute the gradients of the first-order density
 ishift=0 ; ngrad=2 ; nspdentmp=2
 ABI_ALLOCATE(rho1now,(cplex*nfft,nspdentmp,ngrad*ngrad))
 call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspdentmp,paral_kgb,qphon,rhor1tmp,rho1now)
!Now, rho1now(:,:,1) contains the first-order density, and
!rho1now(:,:,2:4) contains the gradients of the first-order density

!Transfer the ground-state density and its gradient
!to spin-polarized storage
 ABI_ALLOCATE(rhortmp,(nfft,nspdentmp,4))
 if(nspden==1)then
   do ii=1,4
     do ir=1,nfft
       rhortmp(ir,1,ii)=kxc(ir,14+2*ii)*half
       rhortmp(ir,2,ii)=kxc(ir,14+2*ii)*half
     end do
   end do
 else
   do ii=1,4
     do ir=1,nfft
       rhortmp(ir,1,ii)=kxc(ir,15+2*ii)
       rhortmp(ir,2,ii)=kxc(ir,14+2*ii)-kxc(ir,15+2*ii)
     end do
   end do
 end if

!DEBUG
!This debug corresponds to cplex=1
!write(std_out,'(a)') ' mkvxcgga3 :  '
!write(std_out,'(a)') '   ispden   ir     rho1now(ir,ispden,1:4)'
!do ispden=1,2
!do ir=1,nfft
!if(ir<=11 .or. mod(ir,301)==0)then
!write(message,'(2i5,a,4es14.6)')ispden,ir,' ',&
!&    rho1now(ir,ispden,1:4)
!call wrtout(std_out,message,'COLL')
!end if
!end do
!end do
!write(std_out,'(a)') ' mkvxcgga3 :  '
!write(std_out,'(a)') '   ispden   ir     rhortmp(ir,ispden,1:4)'
!do ispden=1,2
!do ir=1,nfft
!if(ir<=11 .or. mod(ir,301)==0)then
!write(message,'(2i5,a,4es14.6)')ispden,ir,' ',&
!&    rhortmp(ir,ispden,1:4)
!call wrtout(std_out,message,'COLL')
!end if
!end do
!end do
!write(std_out,'(a)') ' mkvxcgga3 :  '
!write(std_out,'(a)') '   ispden   ir     kxc(ir,1:4)'
!do ir=1,nfft
!if(ir<=11 .or. mod(ir,301)==0)then
!write(message,'(i5,a,4es14.6)')ir,' ',&
!&   kxc(ir,1:4)
!call wrtout(std_out,message,'COLL')
!end if
!end do
!write(std_out,'(a)') ' mkvxcgga3 :  '
!write(std_out,'(a)') '   ispden   ir     kxc(ir,5:8)'
!do ir=1,nfft
!if(ir<=11 .or. mod(ir,301)==0)then
!write(message,'(i5,a,4es14.6)')ir,' ',&
!&   kxc(ir,5:8)
!call wrtout(std_out,message,'COLL')
!end if
!end do
!ENDDEBUG

!Now, rhortmp(:,:,1) contains the GS density, and
!rhortmp(:,:,2:4) contains the gradients of the GS density
!Now, rho1now(:,:,1) contains the first-order density, and
!rho1now(:,:,2:4) contains the gradients of the first-order density

!DEBUG
!ex(:)=zero ; ec(:)=zero
!ENDDEBUG

!Apply the XC kernel
 nspgrad=5
 ABI_ALLOCATE(dnexcdn,(cplex*nfft,nspgrad))

 if(cplex==1)then  ! Treat real case first

   do ir=1,nfft
     r0_up(:)=rhortmp(ir,1,2:4)   ! grad of spin-up GS rho
     r0_dn(:)=rhortmp(ir,2,2:4)   ! grad of spin-down GS rho
     r0(:)=r0_up(:)+r0_dn(:)      ! grad of GS rho
     r1_up(:)=rho1now(ir,1,2:4)   ! grad of spin-up rho1
     r1_dn(:)=rho1now(ir,2,2:4)   ! grad of spin-down rho1
     r1(:)=r1_up(:)+r1_dn(:)      ! grad of GS rho1
     gradrho_gradrho1_up=r1_up(1)*r0_up(1)+r1_up(2)*r0_up(2)+r1_up(3)*r0_up(3)
     gradrho_gradrho1_dn=r1_dn(1)*r0_dn(1)+r1_dn(2)*r0_dn(2)+r1_dn(3)*r0_dn(3)
     gradrho_gradrho1   =r1(1)*r0(1)+r1(2)*r0(2)+r1(3)*r0(3)
!    DEBUG
!    gradrho_gradrho1=zero
!    ENDDEBUG
     dnexcdn(ir,1)=(kxc(ir,1)+kxc(ir,9))*rho1now(ir,1,1)+&
&     kxc(ir,10)*rho1now(ir,2,1)+&
&     kxc(ir,5)*gradrho_gradrho1_up+&
&     kxc(ir,13)*gradrho_gradrho1
     dnexcdn(ir,2)=(kxc(ir,2)+kxc(ir,11))*rho1now(ir,2,1)+&
&     kxc(ir,10)*rho1now(ir,1,1)+&
&     kxc(ir,6)*gradrho_gradrho1_dn+&
&     kxc(ir,14)*gradrho_gradrho1
     coeff_grho_corr=(kxc(ir,13)*rho1now(ir,1,1)+kxc(ir,14)*rho1now(ir,2,1))+&
&     kxc(ir,15)*gradrho_gradrho1
!    DEBUG
!    dnexcdn(ir,1)=(kxc(ir,1)+kxc(ir,9))*rho1now(ir,1,1)+&
!    &                 kxc(ir,10)*rho1now(ir,2,1)+&
!    &                 kxc(ir,5)*gradrho_gradrho1_up
!    dnexcdn(ir,2)=(kxc(ir,2)+kxc(ir,11))*rho1now(ir,2,1)+&
!    &                 kxc(ir,10)*rho1now(ir,1,1)+&
!    &                 kxc(ir,6)*gradrho_gradrho1_dn
!    coeff_grho_corr=zero
!    coeff_grho_corr=(kxc(ir,13)*rho1now(ir,1,1)+kxc(ir,14)*rho1now(ir,2,1))
!    ENDDEBUG
     coeff_grho_up= kxc(ir,5)*rho1now(ir,1,1)+kxc(ir,7)*gradrho_gradrho1_up
     coeff_grho_dn= kxc(ir,6)*rho1now(ir,2,1)+kxc(ir,8)*gradrho_gradrho1_dn
!    DEBUG
!    ex(1:2)=ex(1:2)+kxc(ir,1:2)*rho1now(ir,1:2,1)*rho1now(ir,1:2,1)
!    rho1rho1(1)=r1_up(1)*r1_up(1)+r1_up(2)*r1_up(2)+r1_up(3)*r1_up(3)
!    rho1rho1(2)=r1_dn(1)*r1_dn(1)+r1_dn(2)*r1_dn(2)+r1_dn(3)*r1_dn(3)
!    ex(3:4)=ex(3:4)+kxc(ir,3:4)*rho1rho1(1:2)
!    ex(5)=ex(5)+two*kxc(ir,5)*gradrho_gradrho1_up*rho1now(ir,1,1)
!    ex(6)=ex(6)+two*kxc(ir,6)*gradrho_gradrho1_dn*rho1now(ir,2,1)
!    ex(7)=ex(7)+kxc(ir,7)*gradrho_gradrho1_up**2
!    ex(8)=ex(8)+kxc(ir,8)*gradrho_gradrho1_dn**2

!    ec(1)=ec(1)+kxc(ir,9)*rho1now(ir,1,1)*rho1now(ir,1,1)
!    ec(2)=ec(2)+two*kxc(ir,10)*rho1now(ir,1,1)*rho1now(ir,2,1)
!    ec(3)=ec(3)+kxc(ir,11)*rho1now(ir,2,1)*rho1now(ir,2,1)
!    rho1rho1(3)=r1(1)*r1(1)+r1(2)*r1(2)+r1(3)*r1(3)
!    ec(4)=ec(4)+kxc(ir,12)*rho1rho1(3)
!    ec(5:6)=ec(5:6)+two*kxc(ir,13:14)*gradrho_gradrho1*rho1now(ir,1:2,1)
!    ec(7)=ec(7)+kxc(ir,15)*gradrho_gradrho1**2
!    ENDDEBUG
!    Reuse the storage in rho1now
     rho1now(ir,1,2:4)=r1_up(:)*(kxc(ir,3)+kxc(ir,12))   &
&     +r1_dn(:)*kxc(ir,12)               &
&     +r0_up(:)*coeff_grho_up            &
&     +(r0_up(:)+r0_dn(:))*coeff_grho_corr
     rho1now(ir,2,2:4)=r1_dn(:)*(kxc(ir,4)+kxc(ir,12))   &
&     +r1_up(:)*kxc(ir,12)               &
&     +r0_dn(:)*coeff_grho_dn            &
&     +(r0_up(:)+r0_dn(:))*coeff_grho_corr
   end do

 else ! if cplex==2

   do ir=1,nfft
     r0_up(:)=rhortmp(ir,1,2:4)   ! grad of spin-up GS rho
     r0_dn(:)=rhortmp(ir,2,2:4)   ! grad of spin-down GS rho
     r0(:)=r0_up(:)+r0_dn(:)      ! grad of GS rho
     r1_up(:)=rho1now(2*ir-1,1,2:4)   ! grad of spin-up rho1
     r1im_up(:)=rho1now(2*ir,1,2:4)   ! grad of spin-up rho1 , im part
     r1_dn(:)=rho1now(2*ir-1,2,2:4)   ! grad of spin-down rho1
     r1im_dn(:)=rho1now(2*ir,2,2:4)   ! grad of spin-down rho1 , im part
     r1(:)=r1_up(:)+r1_dn(:)      ! grad of GS rho1
     r1im(:)=r1im_up(:)+r1im_dn(:)      ! grad of GS rho1, im part
     gradrho_gradrho1_up=r1_up(1)*r0_up(1)+r1_up(2)*r0_up(2)+r1_up(3)*r0_up(3)
     gradrho_gradrho1_dn=r1_dn(1)*r0_dn(1)+r1_dn(2)*r0_dn(2)+r1_dn(3)*r0_dn(3)
     gradrho_gradrho1   =r1(1)*r0(1)+r1(2)*r0(2)+r1(3)*r0(3)
     gradrho_gradrho1im_up=r1im_up(1)*r0_up(1)+r1im_up(2)*r0_up(2)+r1im_up(3)*r0_up(3)
     gradrho_gradrho1im_dn=r1im_dn(1)*r0_dn(1)+r1im_dn(2)*r0_dn(2)+r1im_dn(3)*r0_dn(3)
     gradrho_gradrho1im   =r1im(1)*r0(1)+r1im(2)*r0(2)+r1im(3)*r0(3)
     dnexcdn(2*ir-1,1)=(kxc(ir,1)+kxc(ir,9))*rho1now(2*ir-1,1,1)+&
&     kxc(ir,10)*rho1now(2*ir-1,2,1)+&
&     kxc(ir,5)*gradrho_gradrho1_up+&
&     kxc(ir,13)*gradrho_gradrho1
     dnexcdn(2*ir  ,1)=(kxc(ir,1)+kxc(ir,9))*rho1now(2*ir,1,1)+&
&     kxc(ir,10)*rho1now(2*ir,2,1)+&
&     kxc(ir,5)*gradrho_gradrho1im_up+&
&     kxc(ir,13)*gradrho_gradrho1im
     dnexcdn(2*ir-1,2)=(kxc(ir,2)+kxc(ir,11))*rho1now(2*ir-1,2,1)+&
&     kxc(ir,10)*rho1now(2*ir-1,1,1)+&
&     kxc(ir,6)*gradrho_gradrho1_dn+&
&     kxc(ir,14)*gradrho_gradrho1
     dnexcdn(2*ir  ,2)=(kxc(ir,2)+kxc(ir,11))*rho1now(2*ir,2,1)+&
&     kxc(ir,10)*rho1now(2*ir,1,1)+&
&     kxc(ir,6)*gradrho_gradrho1im_dn+&
&     kxc(ir,14)*gradrho_gradrho1im
     coeff_grho_corr=(kxc(ir,13)*rho1now(2*ir-1,1,1)+kxc(ir,14)*rho1now(2*ir-1,2,1))+&
&     kxc(ir,15)*gradrho_gradrho1
     coeffim_grho_corr=(kxc(ir,13)*rho1now(2*ir,1,1)+kxc(ir,14)*rho1now(2*ir,2,1))+&
&     kxc(ir,15)*gradrho_gradrho1im
     coeff_grho_up= kxc(ir,5)*rho1now(2*ir-1,1,1)+kxc(ir,7)*gradrho_gradrho1_up
     coeffim_grho_up= kxc(ir,5)*rho1now(2*ir,1,1)+kxc(ir,7)*gradrho_gradrho1im_up
     coeff_grho_dn= kxc(ir,6)*rho1now(2*ir-1,2,1)+kxc(ir,8)*gradrho_gradrho1_dn
     coeffim_grho_dn= kxc(ir,6)*rho1now(2*ir,2,1)+kxc(ir,8)*gradrho_gradrho1im_dn
!    Reuse the storage in rho1now
     rho1now(2*ir-1,1,2:4)=r1_up(:)*(kxc(ir,3)+kxc(ir,12))   &
&     +r1_dn(:)*kxc(ir,12)               &
&     +r0_up(:)*coeff_grho_up            &
&     +(r0_up(:)+r0_dn(:))*coeff_grho_corr
     rho1now(2*ir  ,1,2:4)=r1im_up(:)*(kxc(ir,3)+kxc(ir,12))   &
&     +r1im_dn(:)*kxc(ir,12)               &
&     +r0_up(:)*coeffim_grho_up            &
&     +(r0_up(:)+r0_dn(:))*coeffim_grho_corr
     rho1now(2*ir-1,2,2:4)=r1_dn(:)*(kxc(ir,4)+kxc(ir,12))   &
&     +r1_up(:)*kxc(ir,12)               &
&     +r0_dn(:)*coeff_grho_dn            &
&     +(r0_up(:)+r0_dn(:))*coeff_grho_corr
     rho1now(2*ir  ,2,2:4)=r1im_dn(:)*(kxc(ir,4)+kxc(ir,12))   &
&     +r1im_up(:)*kxc(ir,12)               &
&     +r0_dn(:)*coeffim_grho_dn            &
&     +(r0_up(:)+r0_dn(:))*coeffim_grho_corr
   end do

 end if

!Now, dnexcdn(:,1)=d(n.exc)/d(n_up)
!dnexcdn(:,2)=d(n.exc)/d(n_down)
!rho1now(:,:,2:4)=part of vxc1 that comes from FFT

!DEBUG
!This debugging section corresponds to cplex=1
!write(std_out,'(a)') ' mkvxcgga3, before xcpot :  '
!write(std_out,'(a)') '   ispden   ir     dnexcdn '
!do ispden=1,2
!do ir=1,nfft
!if(ir<=11 .or. mod(ir,301)==0)then
!write(message,'(2i5,a,4es14.6)')ispden,ir,' ',&
!&    dnexcdn(ir,ispden)
!call wrtout(std_out,message,'COLL')
!end if
!end do
!end do
!write(std_out,'(a)') '   ispden   ir     rho1now(ir,1:4)'
!do ispden=1,2
!do ir=1,nfft
!if(ir<=11 .or. mod(ir,301)==0)then
!write(message,'(2i5,a,4es14.6)')ispden,ir,' ',&
!&    rho1now(ir,ispden,1:4)
!call wrtout(std_out,message,'COLL')
!end if
!end do
!end do
!ENDDEBUG

 ABI_ALLOCATE(vxc1tmp,(cplex*nfft,nspdentmp))
 vxc1tmp(:,:)=zero
 call xcpot (cplex,dnexcdn,gprimd,ishift,mgga,mpi_enreg,nfft,ngfft,ngrad,nspdentmp,&
& nspgrad,paral_kgb,qphon,rho1now,vxc1tmp)

!Transfer the data from spin-polarized storage
 do ispden=1,nspden
   do ir=1,cplex*nfft
     vxc1(ir,ispden)=vxc1tmp(ir,ispden)
   end do
 end do

!call filterpot(paral_kgb,cplex,gmet,gsqcut,nfft,ngfft,nspden,qphon,vxc1)

!DEBUG
!This debugging section works with cplex=1
!write(std_out,'(a)') ' mkvxcgga3 :  '
!write(std_out,'(a)') '   ispden   ir     rhor1tmp(ir)     vxc1tmp(ir) '
!do ispden=1,2
!do ir=1,nfft
!if(ir<=11 .or. mod(ir,301)==0 .or. abs(vxc1tmp(ir,ispden))>1.d6 )then
!write(message,'(2i5,a,4es14.6)')ispden,ir,' ',&
!&    rhor1tmp(ir,ispden),vxc1tmp(ir,ispden)
!call wrtout(std_out,message,'COLL')
!end if
!end do
!end do
!write(std_out,*)' mkvxcgga3 : exit '
!write(std_out,'(a,4es16.6)' )' ex(1:4)=',ex(1:4)
!write(std_out,'(a,4es16.6)' )' ex(5:8)=',ex(5:8)
!write(std_out,'(a,4es16.6)' )' ec(1:4)=',ec(1:4)
!write(std_out,'(a,4es16.6)' )' ec(5:7)=',ec(5:7)
!write(std_out,*)' sum=',sum(ex(1:8))+sum(ec(1:7))
!stop
!ENDDEBUG

 ABI_DEALLOCATE(dnexcdn)
 ABI_DEALLOCATE(rhortmp)
 ABI_DEALLOCATE(rho1now)
 ABI_DEALLOCATE(vxc1tmp)

end subroutine mkvxcgga3
!!***
