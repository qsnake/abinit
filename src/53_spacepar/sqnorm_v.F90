!{\src2tex{textfont=tt}}
!!****f* ABINIT/sqnorm_v
!! NAME
!! sqnorm_v
!!
!!
!! FUNCTION
!! Compute square of the norm of a potential (integral over FFT grid), to obtain
!! a square residual-like quantity (so the sum of product of values
!! is NOT divided by the number of FFT points, and NOT multiplied by the primitive cell volume).
!! Take into account the spin components of the potentials (nspden),
!! and sum over them.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex=if 1, real space function on FFT grid is REAL, if 2, COMPLEX
!!  mpi_enreg=informations about MPI parallelization
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components
!!  opt_storage: 0, if potential is stored as V^up-up, V^dn-dn, Re[V^up-dn], Im[V^up-dn]
!!               1, if potential is stored as V, B_x, B_y, Bz  (B=magn. field)
!!  pot(cplex*nfft,nspden)=real space potential on FFT grid
!!
!! OUTPUT
!!  norm2= value of the square of the norm
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!! PARENTS
!!      rhotov,rhotov3,vtorho,vtorho3
!!
!! CHILDREN
!!      contract_int_ge_val,contract_int_list,timab,xcomm_world,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine sqnorm_v(cplex,mpi_enreg,nfft,norm2,nspden,opt_storage,pot)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sqnorm_v'
 use interfaces_18_timing
#if defined DEBUG_CONTRACT
 use interfaces_32_contract
#endif
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,opt_storage
 real(dp),intent(out) :: norm2
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: pot(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ierr,ifft,ispden,old_paral_level,spaceComm
 real(dp) :: ar
!arrays
 real(dp) :: tsec(2)
!no_abirules
#if defined DEBUG_CONTRACT
 character(len=8) :: subrnm
#endif

! *************************************************************************

#if defined DEBUG_CONTRACT
 subrnm='sqnorm_v'
!Real or complex inputs are coded
 call contract_int_list(subrnm,'cplex',cplex,(/1,2/),2)
 call contract_int_ge_val(subrnm,'nfft',nfft,1)
 call contract_int_list(subrnm,'nspden',nspden,(/1,2,4/),3)
#endif

 norm2=zero
 do ispden=1,min(nspden,2)
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(cplex,ispden,nfft,pot) REDUCTION(+:norm2)
   do ifft=1,cplex*nfft
     norm2=norm2 + pot(ifft,ispden)**2
   end do
!  $OMP END PARALLEL DO
 end do
 if (nspden==4) then
   ar=zero
   do ispden=3,4
!    $OMP PARALLEL DO PRIVATE(ifft) &
!    $OMP&SHARED(cplex,ispden,nfft,pot) REDUCTION(+:ar)
     do ifft=1,cplex*nfft
       ar=ar + pot(ifft,ispden)**2
     end do
!    $OMP END PARALLEL DO
   end do
   if (opt_storage==0) then
     if (cplex==1) then
       norm2=norm2+two*ar
     else
       norm2=norm2+ar
     end if
   else
     norm2=half*(norm2+ar)
   end if
 end if

!XG030513 : MPIWF reduction (addition) on norm2 is needed here
!Init mpi_comm
 if(mpi_enreg%paral_compil_fft==1)then
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
!  call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
!  The same, use xcomm_init to get the right spaceComm.
   if(mpi_enreg%mode_para=='b')then
     spaceComm=mpi_enreg%comm_fft
     call timab(48,1,tsec)
     call xsum_mpi(norm2,spaceComm ,ierr)
     call timab(48,2,tsec)
   end if
   mpi_enreg%paral_level=old_paral_level
 else if (associated(mpi_enreg%nscatterarr)) then
   call xcomm_world(mpi_enreg,spaceComm)
   call timab(48,1,tsec)
   call xsum_mpi(norm2,spaceComm ,ierr)
   call timab(48,2,tsec)
 end if

end subroutine sqnorm_v
!!***
