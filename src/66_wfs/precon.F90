!{\src2tex{textfont=tt}}
!!****f* ABINIT/precon
!!
!! NAME
!! precon
!!
!! FUNCTION
!! precondition $<G|(H-e_{n,k})|C_{n,k}>$
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MT))
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  $cg(2,npw)=<G|C_{n,k}>$.
!!  $eval=current band eigenvalue=<C_{n,k}|H|C_{n,k}>$.
!!  istwf_k=option parameter that describes the storage of wfs
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!!  mpi_enreg=informations about MPI parallelization
!!  nspinor=number of spinorial components of the wavefunctions
!!  $vect(2,npw)=<G|H|C_{n,k}>$.
!!  npw=number of planewaves at this k point.
!!  optekin= 1 if the kinetic energy used in preconditionning is modified
!!             according to Kresse, Furthmuller, PRB 54, 11169 (1996)
!!           0 otherwise
!!
!!
!! OUTPUT
!!  pcon(npw)=preconditionning matrix
!!  vect(2,npw*nspinor)=<G|(H-eval)|C_{n,k}>*(polynomial ratio)
!!
!! PARENTS
!!      cgwf,cgwf3
!!
!! CHILDREN
!!      timab,wrtout,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine precon(cg,eval,istwf_k,kinpw,mpi_enreg,npw,nspinor,optekin,pcon,vect)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'precon'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: istwf_k,npw,nspinor,optekin
 real(dp),intent(in) :: eval
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: cg(2,npw*nspinor),kinpw(npw)
 real(dp),intent(inout) :: pcon(npw),vect(2,npw*nspinor)

!Local variables-------------------------------
!scalars
 integer :: ierr,ig,igs,ipw1,ispinor,old_paral_level,spaceComm
 real(dp) :: ek0,ek0_inv,fac,poly,xx
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

!DEBUG
!write(std_out,*)' precon: debug, enter.'
!ENDDEBUG

!Compute mean kinetic energy of band
 if(istwf_k==1)then
   ek0=zero
   do ispinor=1,nspinor
     igs=(ispinor-1)*npw
!    $OMP PARALLEL DO PRIVATE(ig) REDUCTION(+:ek0) &
!    $OMP&SHARED(cg,igs,kinpw,npw)
     do ig=1+igs,npw+igs
       if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
         ek0=ek0+kinpw(ig-igs)*(cg(1,ig)**2+cg(2,ig)**2)
       end if
     end do
!    $OMP END PARALLEL DO
   end do

 else if (istwf_k>=2)then
   if (istwf_k==2 .and. mpi_enreg%me_g0 == 1)then
     ek0=zero ; ipw1=2
     if(kinpw(1)<huge(0.0_dp)*1.d-11)ek0=0.5_dp*kinpw(1)*cg(1,1)**2
   else
     ek0=zero ; ipw1=1
   end if
   do ispinor=1,nspinor
     igs=(ispinor-1)*npw
!    $OMP PARALLEL DO PRIVATE(ig) REDUCTION(+:ek0) &
!    $OMP&SHARED(cg,ipw1,kinpw,npw)
     do ig=ipw1+igs,npw+igs
       if(kinpw(ig)<huge(0.0_dp)*1.d-11)then
         ek0=ek0+kinpw(ig)*(cg(1,ig)**2+cg(2,ig)**2)
       end if
     end do
!    $OMP END PARALLEL DO
   end do
   ek0=2.0_dp*ek0
 end if
 old_paral_level= mpi_enreg%paral_level
 mpi_enreg%paral_level=3
 call xcomm_init(mpi_enreg,spaceComm)
 call timab(48,1,tsec)
 call xsum_mpi(ek0,spaceComm,ierr)
 call timab(48,2,tsec)
 mpi_enreg%paral_level= old_paral_level
!DEBUG
!write(std_out,*)' precon : ek0,kinpw(1)=',ek0,kinpw(1)
!ENDDEBUG

 if(ek0<1.0d-10)then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' precon : WARNING -',ch10,&
&   '  The mean kinetic energy of a wavefunction vanishes.',ch10,&
&   '  It is reset to 0.1Ha.'
   call wrtout(std_out,message,'PERS')
   ek0=0.1_dp
 end if

 if (optekin==1) then
   ek0_inv=2.0_dp/(3._dp*ek0)
 else
   ek0_inv=1.0_dp/ek0
 end if
!
!Carry out preconditioning
 do ispinor=1,nspinor
   igs=(ispinor-1)*npw
!  $OMP PARALLEL DO PRIVATE(fac,ig,poly,xx) &
!  $OMP&SHARED(cg,ek0_inv,eval,kinpw,igs,npw,vect)
   do ig=1+igs,npw+igs
     if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
       xx=kinpw(ig-igs)*ek0_inv
!      Teter polynomial ratio
       poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
       fac=poly/(poly+16._dp*xx**4)
       if (optekin==1) fac=two*fac
       pcon(ig-igs)=fac
       vect(1,ig)=( vect(1,ig)-eval*cg(1,ig) )*fac
       vect(2,ig)=( vect(2,ig)-eval*cg(2,ig) )*fac
     else
       pcon(ig-igs)=zero
       vect(1,ig)=zero
       vect(2,ig)=zero
     end if
   end do
!  $OMP END PARALLEL DO
 end do

!DEBUG
!write(std_out,*)' precon: debug, exit.'
!ENDDEBUG

end subroutine precon
!!***
