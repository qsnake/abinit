!{\src2tex{textfont=tt}}
!!****f* ABINIT/dyson_sc
!! NAME
!! dyson_sc
!!
!! FUNCTION
!! Compute the diagonal elements of the interacting susceptibility matrix by
!! self-consistently determining the linear density change
!!   $\delta n_{p+1} = \chi_0(\delta v^{\rm ext}+\delta v^{\rm HXC}_{\lambda, p})$,
!! and projecting $\delta n$ on $\delta v^{\rm ext}$, which gives $\chi_{\lambda}$.
!! Here $p$ is the iteration index.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (MF, XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  kernel_diag(npwdiel) = the RPA kernel.
!!  npwdiel = number of planewaves for the susceptibility matrix.
!!  nspden = number of spin-density components.
!!  susmat(2,npwdiel,nspden,npwdiel,nspden) = the Kohn-Sham susceptibility
!!                                            matrix.
!!
!! OUTPUT
!!  susd_isc(npwdiel) = the diagonal elements of the difference between the
!!                      interacting and Kohn-Sham susceptibility matrices.
!!
!! NOTES
!! The linear density change is formally determined by perturbing each
!! component of the external potential and then solving for the induced
!! density change or induced Hartree+XC potential. In this routine only
!! the external potential *change* figures and, for each G component, is
!! of fixed magnitude 'dvext'. Self-consistency is attempted
!! through linear mixing 'bmix' of the induced potential until
!! the density change converges to a fixed tolerance 'btol',
!! and for a maximum number of 'iscmx' iterations.
!! Note that in the RPA case the induced potential for a density change
!! $\delta n(\vec G)$ is simply $4\pi \delta n(\vec G)/G^2$ .
!!
!! WARNINGS
!!  a - Entirely experimental version and meant for the RPA case.
!!  b - Not tested in spin-polarized case.
!!!!
!! PARENTS
!!      acfd_dyson
!!
!! CHILDREN
!!      leave_new,wrtout,zgemv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dyson_sc(kernel_diag,npwdiel,nspden,susd_isc,susmat)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dyson_sc'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: npwdiel,nspden
!arrays
 real(dp),intent(in) :: kernel_diag(npwdiel)
 real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 real(dp),intent(out) :: susd_isc(npwdiel)

!Local variables -------------------------------------------------------
!Settings for the self consistency loop.
!Maximum number of self consistent iterations.
!Mixing factor for linear mixing of density change.
!Relative accuracy in the diagonal density change.
!Strength of the external perturbation.
!scalars
 integer,parameter :: iscmx=10
 integer :: ipw1,ipw2,isc
 real(dp),parameter :: bmix=0.8_dp,btol=1.d-6,dvext=1.d-3
 real(dp) :: bdiff
 character(len=500) :: message
!arrays
 real(dp),parameter :: c1(2)=(/1._dp,0._dp/)
 real(dp),allocatable :: drho_next(:,:),drho_this(:,:),drho_zero(:,:)
 real(dp),allocatable :: susmat_P_kernel(:,:)

!***********************************************************************

!Check input parameters.

 if (nspden > 1) then
   write (message,'(4a)') ch10,&
&   ' dyson_sc: ERROR - ',ch10,&
&   '  dyson_sc does not work yet for nspden > 1.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Allocate memory.

 ABI_ALLOCATE(drho_next,(2,npwdiel))
 ABI_ALLOCATE(drho_this,(2,npwdiel))
 ABI_ALLOCATE(drho_zero,(2,npwdiel))
 ABI_ALLOCATE(susmat_P_kernel,(2,npwdiel))

 drho_zero(:,:) = 0._dp
 drho_next(:,:) = 0._dp

!Loop over external perturbations for each G component.
!(might be restricted to the irreducible wedge of the cell ?)

 do ipw1 = 2,npwdiel

!  Create initial density change.

   drho_zero(1,:) = dvext*(susmat(1,:,1,ipw1,1)-susmat(2,:,1,ipw1,1))
   drho_zero(2,:) = dvext*(susmat(1,:,1,ipw1,1)+susmat(2,:,1,ipw1,1))

!  Initialize iterated density change.

   drho_next(:,:) = drho_zero(:,:)

!  Iterate the density change to self-consistency.

   do isc = 1,iscmx

     drho_this(:,:) = drho_next(:,:)

!    Perform the iteration step.

     if (.true.) then

       do ipw2 = 2,npwdiel
!        (might save some time by multiplying susmat with the kernel in advance).
         susmat_P_kernel(:,ipw2) = drho_this(:,ipw2)*kernel_diag(ipw2)
         drho_next(:,ipw2) = drho_zero(:,ipw2)
       end do

       call ZGEMV('n',npwdiel,npwdiel,c1,susmat,npwdiel,susmat_P_kernel,1,&
&       c1,drho_next,1)

       drho_next(:,1) = 0._dp

     else

       do ipw2 = 2,npwdiel
         susmat_P_kernel(:,ipw2) = drho_this(:,ipw2)*kernel_diag(ipw2)
       end do

       do ipw2 = 2,npwdiel

         drho_next(1,ipw2) = drho_zero(1,ipw2)+sum(( &
&         susmat(1,ipw2,1,:,1)*susmat_P_kernel(1,:)-&
&         susmat(2,ipw2,1,:,1)*susmat_P_kernel(2,:)))

         drho_next(2,ipw2) = drho_zero(2,ipw2)+sum(( &
&         susmat(1,ipw2,1,:,1)*susmat_P_kernel(2,:)+&
&         susmat(2,ipw2,1,:,1)*susmat_P_kernel(1,:)))

       end do

     end if

!    Update the density change by linear mixing.

     drho_next(:,2:) = (1._dp-bmix)*drho_this(:,2:)+bmix*drho_next(:,2:)

!    Terminate if the "diagonal" density change, leading to diagonal of the
!    susceptibility matrix, does not change anymore (quite crude criterion !).

     bdiff = (drho_next(1,ipw1)-drho_this(1,ipw1))&
&     /(drho_next(1,ipw1)+drho_this(1,ipw1))

!    DEBUG
!    write(std_out,*) ' isc = ',isc,' bdiff = ',bdiff
!    ENDDEBUG

     if ((isc > 1).and.(abs(bdiff) < btol)) exit

!    End loop over isc.
   end do

!  Calculate the DIFFERENCE of diagonal components of the interacting and
!  noninteracting susceptibility matrix, these should be real.

   susd_isc(ipw1) = 0.5_dp*(drho_next(1,ipw1)+drho_next(2,ipw1))/dvext &
&   - susmat(1,ipw1,1,ipw1,1)

!  DEBUG
!  check real and imaginary parts, number of iterations.
!  write(std_out,'(a,i6,a,1x,2(1x,e14.6),1x,a,i3)') &
!  &      ' ipw=',     ipw1,&
!  &      ' susd_int=',0.5*(drho_next(1,ipw1)+drho_next(2,ipw1))/dvext,&
!  &                   0.5*(drho_next(2,ipw1)-drho_next(1,ipw1))/dvext,&
!  &      ' isc=',isc
!  ENDDEBUG

!  End loop over ipw1.
 end do

 ABI_DEALLOCATE(drho_next)
 ABI_DEALLOCATE(drho_this)
 ABI_DEALLOCATE(drho_zero)
 ABI_DEALLOCATE(susmat_P_kernel)

!DEBUG
!write(std_out,*) ' %dyson_sc: end and exit'
!ENDDEBUG

end subroutine dyson_sc

!!***
