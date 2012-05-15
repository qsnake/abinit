!{\src2tex{textfont=tt}}
!!****f* ABINIT/dyson_de
!! NAME
!! dyson_de
!!
!! FUNCTION
!! Make a leap-frog step in solving the Dyson equation as a differential equation :
!!  $\chi_{\lambda_{i}}=(\lambda_i-\lambda_{i-1})\chi_{\lambda_{i-1}}
!!                        [d/d\lambda K]_{lambda_{i-1}}\chi_{\lambda_{i-1}}$.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (MF, XG, YMN).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ikernel = 0 if the kernel $(\lambda_i-\lambda_{i-1})[d/d\lambda K]_{lambda_{i-1}}$
!!              is diagonal in reciprocal space,
!!          = 1 otherwise.
!!  kernel_diag = a pointer to a (npwdiel) array containing the real part of
!!                the diagonal elements of the kernel if ikernel == 0.
!!              = an unassociated/unallocated pointer if ikernel = 1.
!!  kernel_full = a pointer to a (2,npwdiel,nspden,npwdiel,nspden) array
!!                containing the full kernel matrix if ikernel = 1.
!!              = an unassociated/unallocated pointer if ikernel = 0.
!!  npwdiel = number of planewaves for the susceptibility matrix.
!!  nspden = number of spin-density components.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  susmat(2,npwdiel,nspden,npwdiel,nspden) = the susceptibility matrix.
!!   on input:  $\chi_{\lambda_{i-1}}$,
!!   on output: $\chi_{\lambda_i}$, updated as described in FUNCTION.
!!
!! NOTES
!! In the RPA case one only needs ikernel = 0.
!!
!! WARNINGS
!!  a - Not tested in spin-polarized case.
!!  b - Tested for the RPA only.
!!
!! TODO
!!  a - Column by column multiplication of susmat*susmat_P_kernel = susmat*kernel*susmat,
!!      to avoid the use of susmat_new, see '!TODO' below.
!!
!! PARENTS
!!      acfd_dyson
!!
!! CHILDREN
!!      zgemm,zhemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dyson_de(ikernel,kernel_diag,kernel_full,npwdiel,nspden,susmat)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dyson_de'
!End of the abilint section

 implicit none

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: ikernel,npwdiel,nspden
!arrays
 real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 real(dp),pointer :: kernel_diag(:),kernel_full(:,:,:,:,:)

!Local variables -------------------------------------------------------
!scalars
 integer :: ipw1,ipw2,isp1,isp2,isp3,ndim
!arrays
 real(dp),parameter :: c0(2)=(/0._dp,0._dp/),c1(2)=(/1._dp,0._dp/)
 real(dp),allocatable :: susmat_P_kernel(:,:,:,:,:),susmat_new(:,:,:,:,:)

!***********************************************************************

!Effective dimension of the susmat matrix, for BLAS subroutines.

 ndim = nspden*npwdiel

!Calculate susmat_P_kernel = susmat*kernel.

 ABI_ALLOCATE(susmat_P_kernel,(2,npwdiel,nspden,npwdiel,nspden))

 if (ikernel == 0) then

!  Case of a diagonal kernel.

   susmat_P_kernel(:,:,:,:,:) = 0._dp

   do isp3 = 1,nspden
     do isp2 = 1,nspden
       do ipw2 = 1,npwdiel
         do isp1 = 1,nspden
           do ipw1 = 1,npwdiel
             susmat_P_kernel(:,ipw1,isp1,ipw2,isp2) = susmat_P_kernel(:,ipw1,isp1,ipw2,isp2)+ &
&             susmat(:,ipw1,isp1,ipw2,isp3)*kernel_diag(ipw2)
           end do
         end do
       end do
     end do
   end do

 else

   if (.true.) then

!    Case of a general kernel: treat susmat as a hermitian matrix.
!    This may be an advantage when used with a sequential BLAS library.

     call ZHEMM('l','l',ndim,ndim,c1,susmat,ndim,&
&     kernel_full,ndim,c0,susmat_P_kernel,ndim)

   else

!    Case of a general kernel: treat susmat as a general matrix.
!    This may be an advantage when used with a parallel BLAS library.

     call ZGEMM('n','n',ndim,ndim,ndim,c1,susmat,ndim,&
&     kernel_full,ndim,c0,susmat_P_kernel,ndim)

   end if

 end if

!Calculate susmat = susmat_P_kernel*susmat = susmat*kernel*susmat.
!TODO : column by column multiplication to avoid the use of susmat_new.

 ABI_ALLOCATE(susmat_new,(2,npwdiel,nspden,npwdiel,nspden))
 susmat_new(:,:,:,:,:) = susmat(:,:,:,:,:)

 if (.true.) then

   call ZHEMM('r','l',ndim,ndim,c1,susmat_new,ndim,&
&   susmat_P_kernel,ndim,c1,susmat,ndim)

 else

   call ZGEMM('n','n',ndim,ndim,ndim,c1,susmat_P_kernel,ndim,&
&   susmat_new,ndim,c1,susmat,ndim)

 end if

 ABI_DEALLOCATE(susmat_new)
 ABI_DEALLOCATE(susmat_P_kernel)

end subroutine dyson_de

!!***
