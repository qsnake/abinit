!{\src2tex{textfont=tt}}
!!****f* ABINIT/meanvalue_g
!! NAME
!! meanvalue_g
!!
!! FUNCTION
!!  Compute the mean value of one wavefunction, in reciprocal space,
!!  for an operator that is real, diagonal in reciprocal space:
!!  <wf|op|wf>
!!  For the time being, only spin-independent operators are treated.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT group (XG,BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  diag(npw)=diagonal operator (real, spin-independent !)
!!  filter= if 1, need to filter on the value of diag, that must be less than huge(0.0d0)*1.d-11
!!          otherwise, should be 0
!!  istwf_k=storage mode of the vectors
!!  mpi_enreg=informations about MPI parallelization
!!  npw=number of planewaves of the vector
!!  nspinor=number of spinor components
!!  vect(2,npw*nspinor)=vector
!!  vect1(2,npw*nspinor*use_ndo)=vector1 (=vector in most of the cases)
!!  use_ndo = says if vect=/vect1
!!
!! OUTPUT
!!  ar=mean value
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      energy,forstrnps,vtowfk,vtowfk3
!!
!! CHILDREN
!!      contract_int_ge_val,contract_int_list,leave_new,timab,wrtout,xcomm_init
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine meanvalue_g(ar,diag,filter,istwf_k,mpi_enreg,npw,nspinor,vect,vect1,use_ndo)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'meanvalue_g'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
#if defined DEBUG_CONTRACT
 use interfaces_32_contract
#endif
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: filter,istwf_k,npw,nspinor,use_ndo
 real(dp),intent(out) :: ar
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: diag(npw),vect(2,npw*nspinor)
 real(dp),intent(in) :: vect1(2,npw*nspinor)

!Local variables-------------------------------
!scalars
 integer :: i1,ierr,ipw,old_paral_level,spaceComm
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
!no_abirules
#if defined DEBUG_CONTRACT
 integer :: ii,isp,ispinor
 character(len=11) :: subrnm
#endif

! *************************************************************************

!DEBUG
!write(std_out,*)' meanvalue_g : enter '
!ENDDEBUG

#if defined DEBUG_CONTRACT
 subrnm='meanvalue_g'
 call contract_int_list(subrnm,'filter',filter,(/0,1/),2)
 call contract_int_list(subrnm,'istwf_k',istwf_k,(/ (ii,ii=1,9) /),9)
 call contract_int_ge_val(subrnm,'npw',npw,1)
 call contract_int_list(subrnm,'nspinor',nspinor,(/1,2/),2)
#endif

 if(nspinor==2 .and. istwf_k/=1)then
   write(message,'(a,a,a,a,a,a,i6,a,i6)') ch10,&
&   ' meanvalue_g: BUG -',ch10,&
&   '  When istwf_k/=1, nspinor must be 1,',ch10,&
&   '  however, nspinor=',nspinor,', and istwf_k=',istwf_k
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(use_ndo==1 .and. (nspinor==2 .or. (istwf_k==2 .and.&
& mpi_enreg%me_g0==1))) then
   write(message,'(a,a)') ch10,' use_ndo==1, not tested'
   call wrtout(std_out,message,'COLL')
!  call leave_new('COLL')
 end if

 ar=zero

!Normal storage mode
 if(istwf_k==1)then

!  No filter
   if(filter==0)then

!    $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar) &
!    $OMP&SHARED(vect,diag,npw)
     do ipw=1,npw
!      ar=ar+diag(ipw)*(vect(1,ipw)**2+vect(2,ipw)**2)
       ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
     end do
!    $OMP END PARALLEL DO
     if(nspinor==2)then
!      $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar) &
!      $OMP&SHARED(vect,diag,npw)
       do ipw=1+npw,2*npw
!        ar=ar+diag(ipw-npw)*(vect(1,ipw)**2+vect(2,ipw)**2)
         ar=ar+diag(ipw-npw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
       end do
!      $OMP END PARALLEL DO
     end if

   else ! will filter

!    $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar) &
!    $OMP&SHARED(vect,diag,npw)
     do ipw=1,npw
       if(diag(ipw)<huge(0.0d0)*1.d-11)then
!        ar=ar+diag(ipw)*(vect(1,ipw)**2+vect(2,ipw)**2)
         ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
       end if
     end do
!    $OMP END PARALLEL DO
     if(nspinor==2)then
!      $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar) &
!      $OMP&SHARED(vect,diag,npw)
       do ipw=1+npw,2*npw
         if(diag(ipw-npw)<huge(0.0d0)*1.d-11)then
!          ar=ar+diag(ipw-npw)*(vect(1,ipw)**2+vect(2,ipw)**2)
           ar=ar+diag(ipw-npw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
         end if
       end do
!      $OMP END PARALLEL DO
     end if ! nspinor==2

   end if ! filter==0

 else if(istwf_k>=2)then

!  XG030513 : MPIWF need to know which proc has G=0

   if(filter==0)then

     i1=1
     if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
!      ar=half*diag(1)*vect(1,1)**2 ; i1=2
       ar=half*diag(1)*vect(1,1)*vect1(1,1) ; i1=2
     end if

!    $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar) &
!    $OMP&SHARED(vect,diag,i1,npw)
     do ipw=i1,npw
!      ar=ar+diag(ipw)*(vect(1,ipw)**2+vect(2,ipw)**2)
       ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
     end do
!    $OMP END PARALLEL DO

   else ! filter/=0

     i1=1
     if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
       if(diag(1)<huge(0.0d0)*1.d-11)then
!        ar=half*diag(1)*vect(1,1)**2 ; i1=2
         ar=half*diag(1)*vect(1,1)*vect1(1,1) ; i1=2
       end if
     end if

!    $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar) &
!    $OMP&SHARED(vect,diag,i1,npw)
     do ipw=i1,npw
       if(diag(ipw)<huge(0.0d0)*1.d-11)then
!        ar=ar+diag(ipw)*(vect(1,ipw)**2+vect(2,ipw)**2)
         ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
       end if
     end do
!    $OMP END PARALLEL DO

   end if ! filter==0

   ar=two*ar

 end if ! istwf_k

!XG030513 : MPIWF need to make reduction on ar and ai .
!Init mpi_comm
 if(mpi_enreg%paral_compil_fft==1)then
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart_3d)
   call timab(48,1,tsec)
   call xsum_mpi(ar,spaceComm ,ierr)
   call timab(48,2,tsec)
   mpi_enreg%paral_level=old_paral_level
 end if

!DEBUG
!write(std_out,*)' meanvalue_g : exit'
!stop
!ENDDEBUG

end subroutine meanvalue_g
!!***
