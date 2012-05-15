!{\src2tex{textfont=tt}}
!!****f* ABINIT/matrixelmt_g
!! NAME
!! matrixelmt_g
!!
!! FUNCTION
!!  Compute a matrix element of two wavefunctions, in reciprocal space,
!!  for an operator that is diagonal in reciprocal space:
!!  <wf1|op|wf2>
!!  For the time being, only spin-independent operators are treated.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  diag(npw)=diagonal operator (real, spin-independent !)
!!  istwf_k=storage mode of the vectors
!!  mpi_enreg=informations about MPI parallelization
!!  needimag=0 if the imaginary part is not needed ; 1 if the imaginary part is needed
!!  npw=number of planewaves of the first vector
!!  nspinor=number of spinor components
!!  vect1(2,npw*nspinor)=first vector
!!  vect2(2,npw*nspinor)=second vector
!!
!! OUTPUT
!!  ai=imaginary part of the matrix element
!!  ar=real part of the matrix element
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      vtowfk3
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


subroutine matrixelmt_g(ai,ar,diag,istwf_k,mpi_enreg,needimag,npw,nspinor,vect1,vect2)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'matrixelmt_g'
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
 integer,intent(in) :: istwf_k,needimag,npw,nspinor
 real(dp),intent(out) :: ai,ar
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: diag(npw),vect1(2,npw*nspinor),vect2(2,npw*nspinor)

!Local variables-------------------------------
!scalars
 integer :: i1,ierr,ipw,old_paral_level,spaceComm
 character(len=500) :: message
!arrays
 real(dp) :: buffer2(2),tsec(2)
!no_abirules
#if defined DEBUG_CONTRACT
 integer :: ii,isp,ispinor
 character(len=12) :: subrnm
#endif

! *************************************************************************

!DEBUG
!write(std_out,*)' matrixelmt_g : enter '
!ENDDEBUG

#if defined DEBUG_CONTRACT
 subrnm='matrixelmt_g'
 call contract_int_list(subrnm,'istwf_k',istwf_k,(/ (ii,ii=1,9) /),9)
 call contract_int_list(subrnm,'needimag',needimag,(/0,1/),2)
 call contract_int_ge_val(subrnm,'npw',npw,1)
 call contract_int_list(subrnm,'nspinor',nspinor,(/1,2/),2)
#endif

 if(nspinor==2 .and. istwf_k/=1)then
   write(message,'(a,a,a,a,a,a,i6,a,i6)') ch10,&
&   ' matrixelmt_g: BUG -',ch10,&
&   '  When istwf_k/=1, nspinor must be 1,',ch10,&
&   '  however, nspinor=',nspinor,', and istwf_k=',istwf_k
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 ar=zero
 if(needimag==1)ai=zero

!Normal storage mode
 if(istwf_k==1)then

!  Need only real part
   if(needimag==0)then

     do ipw=1,npw
       ar=ar+diag(ipw)*(vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw))
     end do
     if(nspinor==2)then
       do ipw=1+npw,2*npw
         ar=ar+diag(ipw-npw)*(vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw))
       end do
     end if

   else ! Need also the imaginary part

     do ipw=1,npw
       ar=ar+diag(ipw)*(vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw))
       ai=ai+diag(ipw)*(vect1(1,ipw)*vect2(2,ipw)-vect1(2,ipw)*vect2(1,ipw))
     end do
     if(nspinor==2)then
       do ipw=1+npw,2*npw
         ar=ar+diag(ipw-npw)*(vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw))
         ai=ai+diag(ipw-npw)*(vect1(1,ipw)*vect2(2,ipw)-vect1(2,ipw)*vect2(1,ipw))
       end do
     end if

   end if ! needimag

 else if(istwf_k>=2)then

!  XG030513 : MPIWF need to know which proc has G=0

   i1=1
   if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
     ar=half*diag(1)*vect1(1,1)*vect2(1,1) ; i1=2
   end if

!  Need only real part
   if(needimag==0)then

     do ipw=i1,npw
       ar=ar+diag(ipw)*(vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw))
     end do
     ar=two*ar

   else ! Need also the imaginary part

     do ipw=i1,npw
       ar=ar+diag(ipw)*(vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw))
       ai=ai+diag(ipw)*(vect1(1,ipw)*vect2(2,ipw)-vect1(2,ipw)*vect2(1,ipw))
     end do
     ar=two*ar ; ai=two*ai

   end if ! needimag

 end if ! istwf_k

!XG030513 : MPIWF need to make reduction on ar and ai .
!Init mpi_comm
 if(mpi_enreg%paral_compil_fft==1)then
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   buffer2(1)=ai
   buffer2(2)=ar
   call timab(48,1,tsec)
   call xsum_mpi(buffer2,spaceComm ,ierr)
!  call xsum_mpi(ai,spaceComm ,ierr)
!  call xsum_mpi(ar,spaceComm ,ierr)
   call timab(48,2,tsec)
   ai=buffer2(1)
   ar=buffer2(2)
   mpi_enreg%paral_level=old_paral_level
 end if


!DEBUG
!write(std_out,*)' matrixelmt_g : exit'
!stop
!ENDDEBUG

end subroutine matrixelmt_g
!!***
