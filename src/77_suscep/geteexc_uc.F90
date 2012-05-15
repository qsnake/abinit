!{\src2tex{textfont=tt}}
!!****f* ABINIT/geteexc_uc
!! NAME
!! geteexc_uc
!!
!! FUNCTION
!! For the calculation of the exchange-correlation energy using
!! the fluctuation-dissipation theorem:
!! Calculate the potential energy associated with the susceptibility
!! matrix with the full Coulomb interaction, in reciprocal space.
!! A polynomial extrapolation is used to determine the $\vec G=0$ contributions.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (MF).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ig_tiny(npw_tiny,3)=index of the n-th shortest G vector along the three
!!    reciprocal space directions
!!  gsq(npwdiel)=squares of G vectors
!!  npwdiel=number of plane waves
!!  npw_tiny=considered number of shortest G vectors
!!  susd(npwdiel)=the (real) diagonal of the susceptibility matrix
!!
!! OUTPUT
!!  energy(npw_tiny)=trace of susceptibility matrix times coulomb interaction
!!  energy_raw=dto, but without G=0 contribution
!!
!! PARENTS
!!      acfd_dyson,xcacfd
!!
!! CHILDREN
!!      polyn_coeff
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine geteexc_uc(energy,energy_raw,gsq,ig_tiny,npwdiel,npw_tiny,&
&  susd)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'geteexc_uc'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_tiny,npwdiel
 real(dp),intent(out) :: energy_raw
!arrays
 integer,intent(in) :: ig_tiny(npw_tiny,3)
 real(dp),intent(in) :: gsq(npwdiel),susd(npwdiel)
 real(dp),intent(out) :: energy(npw_tiny)

!Local variables -------------------------
!scalars
 integer :: iorder,ipw,ir
 real(dp) :: trace
!arrays
 real(dp),allocatable :: gsq_data(:),sus_data(:),sus_gavg(:),sus_gdir(:,:)
 real(dp),allocatable :: sus_poly(:)

! *********************************************************************

!Perform allocations
 ABI_ALLOCATE(gsq_data,(npw_tiny))
 ABI_ALLOCATE(sus_data,(npw_tiny))
 ABI_ALLOCATE(sus_gavg,(npw_tiny))
 ABI_ALLOCATE(sus_poly,(npw_tiny))
 ABI_ALLOCATE(sus_gdir,(npw_tiny,3))

!Trace without G=0 term
 trace=0._dp
 do ipw=1,npwdiel
   if(gsq(ipw) > 1.d-12) trace=trace+susd(ipw)/gsq(ipw)
 end do

!Compute the corrections for G=0 term
!Extrapolate along each direction and then average over all directions
 do ir=1,3
   gsq_data(1:npw_tiny)=gsq(ig_tiny(1:npw_tiny,ir))
   sus_data(1:npw_tiny)=susd(ig_tiny(1:npw_tiny,ir))/gsq_data(1:npw_tiny)
   do iorder=npw_tiny,1,-1
     call polyn_coeff(iorder,gsq_data,sus_data,sus_poly)
     sus_data(iorder)=sus_poly(1)
     sus_gdir(iorder,ir)=sus_data(iorder)
   end do
 end do
 sus_gavg(:)=0._dp
 do ir=1,3
   sus_gavg(1:npw_tiny)=sus_gavg(1:npw_tiny)+sus_gdir(1:npw_tiny,ir)
 end do
 sus_gavg(:)=sus_gavg(:)/3._dp

!Add to trace
 energy(1:npw_tiny)=four_pi*(trace+sus_gavg(1:npw_tiny))
 energy_raw=four_pi*trace

!Perform deallocations
 ABI_DEALLOCATE(gsq_data)
 ABI_DEALLOCATE(sus_data)
 ABI_DEALLOCATE(sus_gavg)
 ABI_DEALLOCATE(sus_poly)
 ABI_DEALLOCATE(sus_gdir)

end subroutine geteexc_uc
!!***
