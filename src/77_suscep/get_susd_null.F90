!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_susd_null
!! NAME
!! get_susd_null
!!
!! FUNCTION
!! Calculate $\vec G=0$ diagonal element of the susceptibility matrix times
!! $G^2$ by polynomial extrapolation (of order npw_tiny - 1), from the
!! diagonal elements of $\chi(\vec G,\vec G) G^2$.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (MF).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ig_tiny(npw_tiny,3)=index of the n-th shortest G vectors along the three
!!    reciprocal space directions
!!  igsq_tiny(npw_tiny)=index of the n-th shortest G vectors
!!  gsq_input(npwdiel)=squares of G vectors
!!  npwdiel=number of plane waves
!!  npw_tiny=considered number of shortest G vectors
!!  sus_input(npwdiel)=the (real) diagonal of the susceptibility matrix
!!
!! OUTPUT
!!  sus_gabs(npw_tiny)=extrapolated values from shortest G vectors
!!  sus_gavg(npw_tiny)=arithmetic mean of the directionally extrapolated
!!    values given in sus_gdir array
!!  sus_gdir(npw_tiny,3)=extrapolated values from shortest G vectors along
!!    the three reciprocal space directions
!!
!! TODO
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


subroutine get_susd_null(ig_tiny,igsq_tiny,gsq_input,npwdiel,npw_tiny,&
&  sus_gabs,sus_gavg,sus_gdir,sus_input)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_susd_null'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments-------------------------------------
!scalars
 integer,intent(in) :: npw_tiny,npwdiel
!arrays
 integer,intent(in) :: ig_tiny(npw_tiny,3),igsq_tiny(npw_tiny)
 real(dp),intent(in) :: gsq_input(npwdiel),sus_input(npwdiel)
 real(dp),intent(out) :: sus_gabs(npw_tiny),sus_gavg(npw_tiny)
 real(dp),intent(out) :: sus_gdir(npw_tiny,3)

!Local variables-------------------------------
!scalars
 integer :: iorder,ir
!arrays
 real(dp),allocatable :: gsq_data(:),sus_data(:),sus_poly(:)

! *************************************************************************

!DEBUG
!write(std_out,*) ' %get_susd_null: enter'
!ENDDEBUG

!Perform allocations
 ABI_ALLOCATE(gsq_data,(npw_tiny))
 ABI_ALLOCATE(sus_data,(npw_tiny))
 ABI_ALLOCATE(sus_poly,(npw_tiny))

!Calculate the extrapolated G=0 value of sus_input data
!from shortest G vectors
 gsq_data(1:npw_tiny)=gsq_input(igsq_tiny(1:npw_tiny))
 sus_data(1:npw_tiny)=sus_input(igsq_tiny(1:npw_tiny))/gsq_data(1:npw_tiny)

!DEBUG
!Beautification note: kg_diel has been removed
!write(std_out,*) ' %get_susd_null: show kg_diel, igsq_tiny, gsq_input, sus_input:'
!do iorder=1,npw_tiny
!write(std_out,'(3i4,2x,i3,2x,es14.7,4x,es14.7)') &
!&  kg_diel(1:3,igsq_tiny(iorder)),igsq_tiny(iorder),gsq_input(igsq_tiny(iorder)),&
!&  sus_input(igsq_tiny(iorder))/gsq_input(igsq_tiny(iorder))
!write(std_out,'(a,2(1x,es14.7))') ' gsq_data,sus_data=',gsq_data(iorder), sus_data(iorder)
!end do
!ENDDEBUG

!Perform the polynomial extrapolation to order npw_tiny in g^2
 do iorder=npw_tiny,1,-1
   call polyn_coeff(iorder,gsq_data,sus_data,sus_poly)
   sus_data(iorder)=sus_poly(1)
   sus_gabs(iorder)=sus_data(iorder)
 end do

!from shortest G vectors along each direction
 do ir=1,3
   gsq_data(1:npw_tiny)=gsq_input(ig_tiny(1:npw_tiny,ir))
   sus_data(1:npw_tiny)=sus_input(ig_tiny(1:npw_tiny,ir))/gsq_data(1:npw_tiny)
   do iorder=npw_tiny,1,-1
     call polyn_coeff(iorder,gsq_data,sus_data,sus_poly)
     sus_data(iorder)=sus_poly(1)
     sus_gdir(iorder,ir)=sus_data(iorder)
   end do
 end do

!from averages over all directions
 sus_gavg(:)=0._dp
 do ir=1,3
   sus_gavg(1:npw_tiny)=sus_gavg(1:npw_tiny)+sus_gdir(1:npw_tiny,ir)
 end do
 sus_gavg(:)=sus_gavg(:)/3._dp

!Perform deallocations
 ABI_DEALLOCATE(gsq_data)
 ABI_DEALLOCATE(sus_data)
 ABI_DEALLOCATE(sus_poly)

end subroutine get_susd_null
!!***
