!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkvxc3
!! NAME
!! mkvxc3
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! due to atomic displacement : assemble the first-order
!! density change with the frozen-core density change, then use
!! the exchange-correlation kernel.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2012 ABINIT group (XG, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!     if 2, COMPLEX
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see rhohxc.f)
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used, otherwise, nfft
!!  option=if 0, work only with the XC core-correction,
!!   if 1, treat both density change and XC core correction
!!   if 2, treat only density change
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor1(cplex*nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density, see n3xccc
!!
!! OUTPUT
!!  vxc1(cplex*nfft,nspden)=change in exchange-correlation potential (including
!!   core-correction, if applicable)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      dyxc13,loop3dte,nres2vres,nstdy3,nstpaw3,rhotov3,xc_kernel
!!
!! CHILDREN
!!      leave_new,matr3inv,mkvxcgga3,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkvxc3(cplex,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,n3xccc,option,&
& paral_kgb,qphon,rhor1,rprimd,vxc1,xccc3d1)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkvxc3'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_56_xc, except_this_one => mkvxc3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,n3xccc,nfft,nkxc,nspden,option,paral_kgb
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: kxc(nfft,nkxc),qphon(3)
 real(dp),intent(in) :: rhor1(cplex*nfft,nspden),rprimd(3,3)
 real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
! integer,save :: npass=0
 integer :: ir
 real(dp) :: rho1_dn,rho1_up,rho1im_dn,rho1im_up,rho1re_dn,rho1re_up
 character(len=500) :: message
!arrays
 real(dp) :: gprimd(3,3),tsec(2)
 real(dp),allocatable :: rhor1tmp(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' mkvxc3 : enter '
!if(option==1)stop
!ENDDEBUG

 call timab(181,1,tsec)

 if(nspden/=1 .and. nspden/=2) then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' mkvxc3 : BUG -',ch10,&
&   '  Only for nspden==1 and 2.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Treat first LDA
 if(nkxc/=23)then

!  Case without non-linear core correction
   if(n3xccc==0 .or. option==2)then

     if(option==0)then  ! No straight XC to compute

       vxc1(:,:)=zero

     else               ! XC, without non-linear XC correction

!      Non-spin-polarized
       if(nspden==1)then
         if(cplex==1)then
           do ir=1,nfft
             vxc1(ir,1)=kxc(ir,1)*rhor1(ir,1)
           end do
         else
           do ir=1,nfft
             vxc1(2*ir-1,1)=kxc(ir,1)*rhor1(2*ir-1,1)
             vxc1(2*ir  ,1)=kxc(ir,1)*rhor1(2*ir  ,1)
           end do
         end if ! cplex==1

!        Spin-polarized
       else
         if(cplex==1)then
           do ir=1,nfft
             rho1_dn=rhor1(ir,1)-rhor1(ir,2)
             vxc1(ir,1)=kxc(ir,1)*rhor1(ir,2)+kxc(ir,2)*rho1_dn
             vxc1(ir,2)=kxc(ir,2)*rhor1(ir,2)+kxc(ir,3)*rho1_dn
           end do
         else
           do ir=1,nfft
             rho1re_dn=rhor1(2*ir-1,1)-rhor1(2*ir-1,2)
             rho1im_dn=rhor1(2*ir  ,1)-rhor1(2*ir  ,2)
             vxc1(2*ir-1,1)=kxc(ir,1)*rhor1(2*ir-1,2)+kxc(ir,2)*rho1re_dn
             vxc1(2*ir  ,1)=kxc(ir,1)*rhor1(2*ir  ,2)+kxc(ir,2)*rho1im_dn
             vxc1(2*ir-1,2)=kxc(ir,2)*rhor1(2*ir-1,2)+kxc(ir,3)*rho1re_dn
             vxc1(2*ir  ,2)=kxc(ir,2)*rhor1(2*ir  ,2)+kxc(ir,3)*rho1im_dn
           end do
         end if ! cplex==1
       end if ! nspden==1

     end if ! option==0

!    Treat case with non-linear core correction
   else

     if(option==0)then

       if(nspden==1)then
         if(cplex==1)then
           do ir=1,nfft
             vxc1(ir,1)=kxc(ir,1)*xccc3d1(ir)
           end do
         else
           do ir=1,nfft
             vxc1(2*ir-1,1)=kxc(ir,1)*xccc3d1(2*ir-1)
             vxc1(2*ir  ,1)=kxc(ir,1)*xccc3d1(2*ir  )
           end do
         end if ! cplex==1
       else
         if(cplex==1)then
           do ir=1,nfft
             vxc1(ir,1)=(kxc(ir,1)+kxc(ir,2))*xccc3d1(ir)*half
             vxc1(ir,2)=(kxc(ir,2)+kxc(ir,3))*xccc3d1(ir)*half
           end do
         else
           do ir=1,nfft
             vxc1(2*ir-1,1)=(kxc(ir,1)+kxc(ir,2))*xccc3d1(2*ir-1)*half
             vxc1(2*ir  ,1)=(kxc(ir,1)+kxc(ir,2))*xccc3d1(2*ir  )*half
             vxc1(2*ir-1,2)=(kxc(ir,2)+kxc(ir,3))*xccc3d1(2*ir-1)*half
             vxc1(2*ir  ,2)=(kxc(ir,2)+kxc(ir,3))*xccc3d1(2*ir  )*half
           end do
         end if ! cplex==1
       end if ! nspden==1

     else ! option/=0

       if(nspden==1)then
         if(cplex==1)then
           do ir=1,nfft
             vxc1(ir,1)=kxc(ir,1)*(rhor1(ir,1)+xccc3d1(ir))
           end do
         else
           do ir=1,nfft
             vxc1(2*ir-1,1)=kxc(ir,1)*(rhor1(2*ir-1,1)+xccc3d1(2*ir-1))
             vxc1(2*ir  ,1)=kxc(ir,1)*(rhor1(2*ir  ,1)+xccc3d1(2*ir  ))
           end do
         end if ! cplex==1
       else
         if(cplex==1)then
           do ir=1,nfft
             rho1_dn=rhor1(ir,1)-rhor1(ir,2) + xccc3d1(ir)*half
             rho1_up=rhor1(ir,2)             + xccc3d1(ir)*half
             vxc1(ir,1)=kxc(ir,1)*rho1_up+kxc(ir,2)*rho1_dn
             vxc1(ir,2)=kxc(ir,2)*rho1_up+kxc(ir,3)*rho1_dn
           end do
         else
           do ir=1,nfft
             rho1re_dn=rhor1(2*ir-1,1)-rhor1(2*ir-1,2) + xccc3d1(2*ir-1)*half
             rho1im_dn=rhor1(2*ir  ,1)-rhor1(2*ir  ,2) + xccc3d1(2*ir  )*half
             rho1re_up=rhor1(2*ir-1,2)                 + xccc3d1(2*ir-1)*half
             rho1im_up=rhor1(2*ir  ,2)                 + xccc3d1(2*ir  )*half
             vxc1(2*ir-1,1)=kxc(ir,1)*rho1re_up+kxc(ir,2)*rho1re_dn
             vxc1(2*ir  ,1)=kxc(ir,1)*rho1im_up+kxc(ir,2)*rho1im_dn
             vxc1(2*ir-1,2)=kxc(ir,2)*rho1re_up+kxc(ir,3)*rho1re_dn
             vxc1(2*ir  ,2)=kxc(ir,2)*rho1im_up+kxc(ir,3)*rho1im_dn
           end do
         end if ! cplex==1
       end if ! nspden==1

     end if ! option==0

   end if ! n3xccc==0

!  Treat GGA
 else

   ABI_ALLOCATE(rhor1tmp,(cplex*nfft,2))

!  Generates gprimd
   call matr3inv(rprimd,gprimd)

!  First transfer the data to spin-polarized storage
   if(option==1 .or. option==2)then   ! Treat the density change
     if(nspden==1)then
       do ir=1,cplex*nfft
         rhor1tmp(ir,1)=rhor1(ir,1)*half
         rhor1tmp(ir,2)=rhor1(ir,1)*half
       end do
     else
       do ir=1,cplex*nfft
         rho1_dn=rhor1(ir,1)-rhor1(ir,2)
         rhor1tmp(ir,1)=rhor1(ir,2)
         rhor1tmp(ir,2)=rho1_dn
       end do
     end if ! nspden==1
   else
     do ir=1,cplex*nfft
       rhor1tmp(ir,1)=zero
       rhor1tmp(ir,2)=zero
     end do
   end if
   if( (option==0 .or. option==1) .and. n3xccc/=0)then
     do ir=1,cplex*nfft
       rhor1tmp(ir,1)=rhor1tmp(ir,1)+xccc3d1(ir)*half
       rhor1tmp(ir,2)=rhor1tmp(ir,2)+xccc3d1(ir)*half
     end do
   end if

   call mkvxcgga3(cplex,gprimd,kxc,mpi_enreg,nfft,ngfft,nkxc,&
&   nspden,paral_kgb,qphon,rhor1tmp,vxc1)

   ABI_DEALLOCATE(rhor1tmp)

 end if
 call timab(181,2,tsec)

!DEBUG
!write(std_out,'(a)') ' mkvxc3 :  '
!write(std_out,'(a)') '   ir     rhor1(ir)     vxc1(ir)       kxc(ir) '
!if(cplex==1)then
!do ir=1,nfft
!if(ir<=11 .or. mod(ir,301)==0 )then
!write(message,'(i5,a,4es14.6)')ir,' ',&
!&    rhor1(ir,1),vxc1(ir,1),kxc(ir,1)
!call wrtout(std_out,message,'COLL')
!end if
!end do
!else
!do ir=1,nfft
!if(ir<=11 .or. mod(ir,301)==0 )then
!write(message,'(i5,a,4es14.6)')ir,' ',&
!&    rhor1(2*ir-1,1),vxc1(2*ir-1,1),kxc(ir,1)
!call wrtout(std_out,message,'COLL')
!write(message,'(i5,a,4es14.6)')ir,' ',&
!&    rhor1(2*ir  ,1),vxc1(2*ir  ,1)
!call wrtout(std_out,message,'COLL')
!end if
!end do
!end if
!npass=npass+1
!if(npass==3)stop
!ENDDEBUG

!DEBUG
!write(std_out,*)' mkvxc3 : exit'
!if(option==1)stop
!vxc1(:,:)=zero
!stop
!ENDDEBUG

end subroutine mkvxc3
!!***
