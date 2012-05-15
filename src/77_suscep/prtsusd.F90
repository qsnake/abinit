!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtsusd
!! NAME
!! prtsusd
!!
!! FUNCTION
!! Prints out the real diagonal elements of the type susceptibility matrix
!! times Coulomb interaction, where the $\vec G=0$ contribution exists and
!! can be determined by extrapolation from the values for small $G$ (for
!! the susceptibility, the $\vec G=0$ diagaonal element is zero).
!! Note: Cannot be used for the density, since $|n(\vec G)|^2/G^2$ diverges.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (MF).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gsq(npwdiel)=squares of G vectors
!!  ig_tiny(npw_tiny,3)=plane wave index for shortest G vectors along
!!   principal axes
!!  index_g(npwdiel)=sorting index of plane waves with increasing G^2
!!  npw_tiny=number of shortest G vectors considered
!!  npwdiel=number of plane waves for susceptibility matrix
!!  optprt=1: print all elements
!!        =other: print average over same or similar G vectors
!!  susd(npwdiel)=the real diagonal elements of the susceptibility matrix
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! In development.
!!
!! PARENTS
!!      xcacfd
!!
!! CHILDREN
!!      polyn_coeff
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prtsusd(gsq,ig_tiny,index_g,npw_tiny,npwdiel,optprt,susd)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtsusd'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_tiny,npwdiel,optprt
!arrays
 integer,intent(in) :: ig_tiny(npw_tiny,3),index_g(npwdiel)
 real(dp),intent(in) :: gsq(npwdiel),susd(npwdiel)

!Local variables-------------------------------
!scalars
 integer :: icount,iorder,ipw1,ir
 real(dp) :: gsqthis,scalefactor,susvc
 logical :: tprt
!arrays
 real(dp),allocatable :: gsq_data(:),sus_data(:),sus_gavg(:),sus_gdir(:,:)
 real(dp),allocatable :: sus_poly(:)

! *************************************************************************
!DEBUG
!write(std_out,*) ' %prtsusd: ENTER'
!ENDDEBUG

!Perform allocations
 ABI_ALLOCATE(gsq_data,(npw_tiny))
 ABI_ALLOCATE(sus_data,(npw_tiny))
 ABI_ALLOCATE(sus_gavg,(npw_tiny))
 ABI_ALLOCATE(sus_gdir,(npw_tiny,3))
 ABI_ALLOCATE(sus_poly,(npw_tiny))

!Scaling for e-e interaction energy
!$-1/2\pi$ [fluctuation-dissipation theorem] * $4\pi$ [Coulomb interaction] = -2
 scalefactor=-2._dp ! -1/2pi[fluctuation-dissipation theorem] * 4pi[Coulomb interaction]

!Directional extrapolation followed by spatial average for G=0 term, using polynomial
 do ir=1,3
   gsq_data(1:npw_tiny)=gsq(ig_tiny(1:npw_tiny,ir))
   sus_data(1:npw_tiny)=susd(ig_tiny(1:npw_tiny,ir))/gsq_data(1:npw_tiny)
   do iorder=npw_tiny,1,-1
     call polyn_coeff(iorder,gsq_data,sus_data,sus_poly)
     sus_gdir(iorder,ir)=sus_poly(1)
   end do
 end do
 sus_gavg(:)=0._dp
 do ir=1,3
   sus_gavg(1:npw_tiny)=sus_gavg(1:npw_tiny)+sus_gdir(1:npw_tiny,ir)
 end do
 sus_gavg(:)=sus_gavg(:)/3._dp

!Average over G vectors with same or similar length
 tprt=.false.
 icount=0
 do ipw1=1,npwdiel
   gsqthis=gsq(index_g(ipw1))
   icount=icount+1
   if(ipw1==1) then
     if(gsqthis /= 0._dp) write(std_out,*) ' prtsusd: WARNING: G=0 does not correspond to ipw=0 .'
     susvc=sus_gavg(npw_tiny)
     tprt=.true.
   else if(ipw1==npwdiel) then
     if(abs(gsq(index_g(ipw1-1))-gsqthis) < 1.d-2 ) then
       susvc=susvc+susd(index_g(ipw1))/gsqthis
     else
       susvc=susd(index_g(ipw1))/gsqthis
     end if
     tprt=.true.
   else
     susvc=susvc+susd(index_g(ipw1))/gsqthis
     if( abs(gsq(index_g(ipw1+1))-gsqthis) > 1.d-2 )then
       tprt=.true.
     end if
   end if

!  Print out all values
   if(optprt==1) then
     if(ipw1==1) then
       if(gsqthis /= 0._dp) write(std_out,*) ' prtsusd: WARNING: G=0 does not correspond to ipw=0 .'
       susvc=sus_gavg(npw_tiny)
     else
       susvc=susd(index_g(ipw1))/gsqthis
     end if
     susvc=scalefactor*susvc
     write(tmp_unit,'(2(es14.7,1x))') gsqthis,susvc
!    Print out averaged values only
   else
     if(tprt) then
       susvc=scalefactor*susvc/dble(icount)
       write(tmp_unit,'(2(es14.7,1x),i5)') gsqthis,susvc,icount
!      write(tmp_unit,'(es14.7,1x,es14.7)')   gsqthis,susvc
       icount=0
       susvc=0._dp
       tprt=.false.
     end if
   end if

!  End loop over ipw1
 end do

!Perform deallocations
 ABI_DEALLOCATE(gsq_data)
 ABI_DEALLOCATE(sus_data)
 ABI_DEALLOCATE(sus_gavg)
 ABI_DEALLOCATE(sus_gdir)
 ABI_DEALLOCATE(sus_poly)

end subroutine prtsusd
!!***
