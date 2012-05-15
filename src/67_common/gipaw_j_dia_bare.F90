!{\src2tex{textfont=tt}}
!!****f* ABINIT/gipaw_j_dia_bare
!! NAME
!! gipaw_j_dia_bare
!!
!! FUNCTION
!! compute current vector field due to bare diamagnetic term in GIPAW
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! nfft, ngfft(18): number of grid points and grid details
!! nhat(nfft,nspden): augmentation charge density
!! nspden: number of spin polarization components
!! rhor(nfft,nspden): total charge density (pseudo + nhat)
!! rprimd(3,3): unit cell dimensions
!!
!! OUTPUT
!! jdia(B dir=1..3, j_i=1..3, nfft), current field due to bare diamagnetic term
!!
!! SIDE EFFECTS
!!
!! NOTES
!! The GIPAW current response includes a bare diamagnetic term
!! $-\frac{1}{2c}n(\mathbf{r'})\mathbf{B\wedge r'}$, where the charge density
!! $n$ is the pseudo charge density ($n-\hat{n}$). Here we compute this vector field
!! for each direction of the external B field. Notice that in the GIPAW papers (for 
!! example Yates, Pickard, Mauri, PRB 76, 024401 (2007), this term is absorbed into the
!! perturbation calculation through use of a generalized f-sum rule. We do not do that 
!! here but rather treat it directly.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine gipaw_j_dia_bare(jdia,nfft,ngfft,nhat,nspden,rhor,rprimd)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gipaw_j_dia_bare'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! jdia( (Bvec(1:3), jvec(1:3), (gridpts) )
!scalars
 integer,intent(in) :: nfft,nspden
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: nhat(nfft,nspden),rhor(nfft,nspden),rprimd(3,3)
 real(dp),intent(out) :: jdia(3,3,nfft)

!Local variables-------------------------------
!scalars
 integer :: igfft,igfft1,igfft2,igfft3
!arrays
 real(dp) :: B(3,3),fofr(nfft),rvec(3)
!no_abirules
!

! ************************************************************************

!DEBUG
!write(std_out,*)' gipaw_j_dia_bare : enter'
!ENDDEBUG
!make sure current field is set to zero everywhere
 jdia(:,:,:) = zero

!define the directions of the external B field
 B = zero
 B(1,1) = 1.0; B(2,2) = 1.0; B(3,3) = 1.0;

!construct charge density. Note that $\mathrm{rhor} = \tilde{n} + \hat{n}$, therefore we
!must subtract $\hat{n}$. We also multiply by -1 to convert the electron density to the charge density
!note also that in the nspden > 1 cases, the ispden = 1 component is still the total charge density
 fofr(:) = -(rhor(:,1) - nhat(:,1)) 
 do igfft1 = 1, ngfft(1) ! looping over points in the grid
   do igfft2 = 1, ngfft(2)
     do igfft3 = 1, ngfft(3)
!      here is the index of the grid point
       igfft = (igfft3-1)*ngfft(2)*ngfft(1) + (igfft2-1)*ngfft(1) + igfft1
!      here is the grid point in cartesian coordinates: igfft1 along rprimd(:,1), 
!      igfft2 along rprimd(:,2), and igfft3 along rprimd(:,3)
       rvec(:) = rprimd(:,1)*float(igfft1-1)/float(ngfft(1)) + &
&       rprimd(:,2)*float(igfft2-1)/float(ngfft(2)) + &
&       rprimd(:,3)*float(igfft3-1)/float(ngfft(3))
!      do idir = 1, 3 ! loop over magnetic field directions
!      accumulate -rho*(B X r')/2c into jdia
!      call acrossb(B(idir,:),rvec,Bxr)
!      jdia(idir,:,igfft) = jdia(idir,:,igfft) - fofr(igfft)*Bxr(:)/(2.0*Sp_Lt)
!      end do ! end loop over B field directions
       jdia(1,3,igfft)=fofr(igfft)*rvec(2)/Sp_Lt
       jdia(2,1,igfft)=fofr(igfft)*rvec(3)/Sp_Lt
       jdia(3,2,igfft)=fofr(igfft)*rvec(1)/Sp_Lt
     end do ! end loop over ngfft(3)
   end do ! end loop over ngfft(2)
 end do ! end loop over  ngfft(1)

!DEBUG
!write(std_out,*)' gipaw_j_dia_bare : exit '
!stop
!ENDDEBUG

end subroutine gipaw_j_dia_bare
!!***
