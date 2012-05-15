!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotprod_v
!! NAME
!! dotprod_v
!!
!!
!! FUNCTION
!! Compute dot product of two potentials (integral over FFT grid), to obtain
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
!!  cplex=if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  mpi_enreg=informations about MPI parallelization
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components
!!  opt_storage: 0, if potentials are stored as V^up-up, V^dn-dn, Re[V^up-dn], Im[V^up-dn]
!!               1, if potentials are stored as V, B_x, B_y, Bz  (B=magn. field)
!!  pot1(cplex*nfft,nspden)=first real space potential on FFT grid
!!  pot2(cplex*nfft,nspden)=second real space potential on FFT grid
!!
!! OUTPUT
!!  dotr= value of the dot product
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!! PARENTS
!!      dens_in_sph
!!
!! CHILDREN
!!      contract_int_ge_val,contract_int_list,timab,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dotprod_v(cplex,dotr,mpi_enreg,nfft,nspden,opt_storage,pot1,pot2)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dotprod_v'
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
 real(dp),intent(out) :: dotr
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: pot1(cplex*nfft,nspden),pot2(cplex*nfft,nspden)

!Local variables-------------------------------
#if defined DEBUG_CONTRACT
 character(len=9) :: subrnm
#endif
!scalars
 integer :: ierr,ifft,ispden,old_paral_level,spaceComm
 real(dp) :: ar
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

#if defined DEBUG_CONTRACT
 subrnm='dotprod_v'
!Real or complex inputs are coded
 call contract_int_list(subrnm,'cplex',cplex,(/1,2/),2)
 call contract_int_ge_val(subrnm,'nfft',nfft,1)
 call contract_int_list(subrnm,'nspden',nspden,(/1,2,4/),3)
#endif

 dotr=zero
 do ispden=1,min(nspden,2)
!  $OMP PARALLEL DO PRIVATE(ifft) SHARED(ispden,nfft,pot1,pot2) REDUCTION(+:dotr)
   do ifft=1,cplex*nfft
     dotr=dotr + pot1(ifft,ispden)*pot2(ifft,ispden)
   end do
!  $OMP END PARALLEL DO
 end do
 if (nspden==4) then
   ar=zero
   do ispden=3,4
!    $OMP PARALLEL DO PRIVATE(ifft) SHARED(ispden,nfft,pot1,pot2) REDUCTION(+:ar)
     do ifft=1,cplex*nfft
       ar=ar + pot1(ifft,ispden)*pot2(ifft,ispden)
     end do
!    $OMP END PARALLEL DO
   end do
   if (opt_storage==0) then
     if (cplex==1) then
       dotr=dotr+two*ar
     else
       dotr=dotr+ar
     end if
   else
     dotr=half*(dotr+ar)
   end if
 end if

!XG030513 : MPIWF reduction (addition) on dotr is needed here
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

end subroutine dotprod_v
!!***
