!{\src2tex{textfont=tt}}
!!****f* ABINIT/dyson_gl
!! NAME
!! dyson_gl
!!
!! FUNCTION
!! Compute the interacting susceptibility matrix from the first order solution
!! of the Dyson equation :
!!  $\chi_{\lambda=1}-\chi_0 = \chi_0[d/d\lambda K]_{\lambda=0}\chi_0$
!! as applied by Lein, Dobson and Gross, J. Comput. Chem. 20, 12 (1999).
!! This routine only returns the real part of the diagonal elements of
!! $\chi_{\lambda=1}-\chi_0$ (summed over spin-density components if appropriate).
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (MF, XG, YMN).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ikernel = 0 if the kernel $[d/d\lambda K]_{\lambda=0}$ is diagonal
!!              in reciprocal space,
!!          = 1 otherwise.
!!  kernel_diag = a pointer to a (npwdiel) array containing the real part of
!!                the diagonal elements of the kernel if ikernel == 0.
!!              = an unassociated/unallocated pointer if ikernel = 1.
!!  kernel_full = a pointer to a (2,npwdiel,nspden,npwdiel,nspden) array
!!                containing the full kernel matrix if ikernel = 1.
!!              = an unassociated/unallocated pointer if ikernel = 0.
!!  npwdiel = number of planewaves for the susceptibility matrix.
!!  nspden = number of spin-density components.
!!  susmat(2,npwdiel,nspden,npwdiel,nspden) = the non-interacting susceptibility
!!                                            matrix $\chi_0$.
!!
!! OUTPUT
!!  susd_LDG(npwdiel) = real part of the diagonal of $\chi_{\lambda=1}-\chi_0$
!!   (summed over spin-density components if appropriate).
!!
!! NOTES
!! In the RPA case one only needs ikernel = 0.
!!
!! WARNINGS
!!  a - Not tested in spin-polarized case.
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


subroutine dyson_gl(ikernel,kernel_diag,kernel_full,npwdiel,nspden,&
&                    susd_LDG,susmat)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dyson_gl'
!End of the abilint section

 implicit none

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: ikernel,npwdiel,nspden
!arrays
 real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 real(dp),intent(out) :: susd_LDG(npwdiel)
 real(dp),pointer :: kernel_diag(:),kernel_full(:,:,:,:,:)

!Local variables -------------------------------------------------------
!scalars
 integer :: ipw1,ipw2,ipw3,isp1,isp2,isp3,ndim
!arrays
 real(dp),parameter :: c0(2)=(/0._dp,0._dp/),c1(2)=(/1._dp,0._dp/)
 real(dp),allocatable :: susmat_P_kernel(:,:,:,:,:)

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

!Calculate the sum over spin-density components of the real part of the
!diagonal of susmat_P_kernel*susmat = susmat*kernel*susmat.

 susd_LDG(:) = 0._dp

 do isp2 = 1,nspden
   do isp1 = 1,nspden
     do ipw1 = 1,npwdiel
       do isp3 = 1,nspden
         do ipw3 = 1,npwdiel
           susd_LDG(ipw1) = susd_LDG(ipw1) &
&           +susmat_P_kernel(1,ipw1,isp1,ipw3,isp3)*susmat(1,ipw3,isp3,ipw1,isp2) &
&           -susmat_P_kernel(2,ipw1,isp1,ipw3,isp3)*susmat(2,ipw3,isp3,ipw1,isp2)
         end do
       end do
     end do
   end do
 end do

 ABI_DEALLOCATE(susmat_P_kernel)

end subroutine dyson_gl

!!***
