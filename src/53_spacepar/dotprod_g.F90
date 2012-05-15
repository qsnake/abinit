!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotprod_g
!! NAME
!! dotprod_g
!!
!!
!! FUNCTION
!! Compute dot product of complex vectors vect1 and vect2 (can be the same)
!! Take into account the storage mode of the vectors (istwf_k)
!! If option=1, compute only real part, if option=2 compute also imaginary
!! part. If the number of calls to the dot product scales quadratically
!! with the volume of the system, it is preferable not to
!! call the present routine, but but to write a specially
!! optimized routine, that will avoid many branches related to
!! the existence of 'option' and 'istwf_k'.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  istwf_k=option parameter that describes the storage of wfs
!!  vect1(2,npw)=first vector (one should take its complex conjugate)
!!  vect2(2,npw)=second vector
!!  mpi_enreg=informations about MPI parallelization
!!  npw= (effective) number of planewaves at this k point (including spinorial level)
!!  option= 1 if only real part to be computed,
!!          2 if both real and imaginary.
!!          3 if in case istwf_k==1 must compute real and imaginary parts,
!!               but if  istwf_k >1 must compute only real part
!!
!! OUTPUT
!!  $doti=\Im ( <vect1|vect2> )$ , output only if option=2 and eventually option=3.
!!  $dotr=\Re ( <vect1|vect2> )$
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!! PARENTS
!!      cgwf,cgwf3,corrmetalwf1,eig1fixed,eig2tot,extrapwf,mkresi,nonlop_gpu
!!      nstpaw3,nstwf3,nstwf4,recip_ylm,vtowfk3
!!
!! CHILDREN
!!      contract_int_ge_val,contract_int_list,timab,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw,option,vect1,vect2)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dotprod_g'
 use interfaces_18_timing
#if defined DEBUG_CONTRACT
 use interfaces_32_contract
#endif
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,npw,option
 real(dp),intent(out) :: doti,dotr
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: vect1(2,npw),vect2(2,npw)

!Local variables-------------------------------
!scalars
 integer :: ierr,ipw,old_paral_level,spaceComm
!arrays
 real(dp) :: dotarr(2),tsec(2)
!no_abirules
#if defined DEBUG_CONTRACT
 integer :: ii
 character(len=9) :: subrnm
#endif

! *************************************************************************

#if defined DEBUG_CONTRACT
 subrnm='dotprod_g'
 call contract_int_list(subrnm,'istwf_k',istwf_k,(/ (ii,ii=1,9) /),9)
 call contract_int_ge_val(subrnm,'npw',npw,1)
 call contract_int_list(subrnm,'option',option,(/1,2,3/),3)
#endif

 if(istwf_k==1)then
!  === General k-point
   if(option==1)then
     dotr=zero
!    $OMP PARALLEL DO ORDERED PRIVATE(ipw) REDUCTION(+:dotr) SHARED(vect1,vect2,npw)
     do ipw=1,npw
       dotr=dotr+vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw)
     end do
!    $OMP END PARALLEL DO
   else
     dotr=zero ; doti=zero
!    $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:doti,dotr) SHARED(vect1,vect2,npw)
     do ipw=1,npw
       dotr=dotr+vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw)
       doti=doti-vect1(2,ipw)*vect2(1,ipw)+vect1(1,ipw)*vect2(2,ipw)
     end do
!    $OMP END PARALLEL DO
   end if
 else if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
!  === Gamma k-point
   dotr=half*vect1(1,1)*vect2(1,1)
!  $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:dotr) SHARED(vect1,vect2,npw)
   do ipw=2,npw
     dotr=dotr+vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw)
   end do
!  $OMP END PARALLEL DO
   dotr=two*dotr
   if(option==2) doti=zero
 else
!  === Other TR k-points
   dotr=zero
!  $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:dotr) SHARED(vect1,vect2,npw)
   do ipw=1,npw
     dotr=dotr+vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw)
   end do
!  $OMP END PARALLEL DO
   dotr=two*dotr
   if(option==2) doti=zero
 end if

!Reduction in case of parallelism
 if(mpi_enreg%paral_compil_fft==1)then
   call timab(48,1,tsec)
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   if (option==1.or.istwf_k/=1) then
     call xsum_mpi(dotr,spaceComm,ierr)
     if (mpi_enreg%paral_spin==1) then
       call xsum_mpi(dotr,mpi_enreg%comm_spin,ierr)
     end if
   else
     dotarr(1)=dotr ; dotarr(2)=doti
     call xsum_mpi(dotarr,spaceComm,ierr)
     if (mpi_enreg%paral_spin==1) then
       call xsum_mpi(dotarr, mpi_enreg%comm_spin, ierr)
     end if
     dotr=dotarr(1) ; doti=dotarr(2)
   end if
   mpi_enreg%paral_level=old_paral_level
   call timab(48,2,tsec)
 end if

end subroutine dotprod_g
!!***
