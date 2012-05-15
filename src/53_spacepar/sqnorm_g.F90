!{\src2tex{textfont=tt}}
!!****f* ABINIT/sqnorm_g
!! NAME
!! sqnorm_g
!!
!! FUNCTION
!! Compute the square of the norm of one complex vector vecti, in reciprocal space
!! Take into account the storage mode of the vector (istwf_k)
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
!!  mpi_enreg=informations about MPI parallelization
!!  npwsp= (effective) number of planewaves at this k point.
!!  vect(2,npwsp)=the vector in reciprocal space (npw*nspinor, usually)
!!
!! OUTPUT
!!  dotr= <vect|vect>
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!! PARENTS
!!      cgwf,cgwf3,dens_in_sph,mkresi,vtowfk3
!!
!! CHILDREN
!!      contract_int_ge_val,contract_int_list,timab,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine sqnorm_g(dotr,istwf_k,mpi_enreg,npwsp,vect)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sqnorm_g'
 use interfaces_18_timing
#if defined DEBUG_CONTRACT
 use interfaces_32_contract
#endif
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,npwsp
 real(dp),intent(out) :: dotr
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: vect(2,npwsp)

!Local variables-------------------------------
!scalars
 integer :: ierr,ipw,old_paral_level,spaceComm
!arrays
 real(dp) :: tsec(2)
!no_abirules
#if defined DEBUG_CONTRACT
 integer :: ii
 character(len=8) :: subrnm
#endif

! *************************************************************************

!DEBUG
!write(std_out,*)' sqnorm_g: debug, enter.'
!ENDDEBUG


#if defined DEBUG_CONTRACT
 subrnm='sqnorm_g'
 call contract_int_list(subrnm,'istwf_k',istwf_k,(/ (ii,ii=1,9) /),9)
 call contract_int_ge_val(subrnm,'npwsp',npwsp,1)
#endif

 if(istwf_k==1)then
!  General k-point

   dotr=0.0d0
!  $OMP PARALLEL DO ORDERED PRIVATE(ipw) REDUCTION(+:dotr) &
!  $OMP&SHARED(vect,npwsp)
   do ipw=1,npwsp
     dotr=dotr+vect(1,ipw)**2+vect(2,ipw)**2
   end do
!  $OMP END PARALLEL DO

 else

!  Gamma k-point
   if(istwf_k==2 .and. mpi_enreg%me_g0==1)then

     dotr=half*vect(1,1)**2
!    $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:dotr) &
!    $OMP&SHARED(vect,npwsp)
     do ipw=2,npwsp
       dotr=dotr+vect(1,ipw)**2+vect(2,ipw)**2
     end do
!    $OMP END PARALLEL DO

   else

!    Other TR k-points
     dotr=0.0d0
!    $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:dotr) &
!    $OMP&SHARED(vect,npwsp)
     do ipw=1,npwsp
       dotr=dotr+vect(1,ipw)**2+vect(2,ipw)**2
     end do
!    $OMP END PARALLEL DO

   end if

   dotr=2.0d0*dotr

 end if

!XG030513 : MPIWF reduction on dotr is needed here
!Init mpi_comm
 if(mpi_enreg%paral_compil_fft==1)then
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   call timab(48,1,tsec)
   call xsum_mpi(dotr,spaceComm ,ierr)
   call timab(48,2,tsec)
   mpi_enreg%paral_level=old_paral_level
 end if

!DEBUG
!write(std_out,*)' sqnorm_g: debug, exit.'
!ENDDEBUG

end subroutine sqnorm_g
!!***
