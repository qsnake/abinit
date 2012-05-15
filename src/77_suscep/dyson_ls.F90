!{\src2tex{textfont=tt}}
!!****f* ABINIT/dyson_ls
!! NAME
!! dyson_ls
!!
!! FUNCTION
!! Update the interacting susceptibility matrix in a step of the recursive
!! solution of the Dyson equation as a set of linear equations :
!!  $ (1 - \chi_{input}*K_Hxc) \chi_{output} = \chi_{input} $.
!! where (if ikernel=0) : K_Hxc is diagonal
!!       (if ikernel=1) : K_Hxc is a full matrix
!!       (if ikernel=2) : K_Hxc is made of two parts, one diagonal (K_H), and one full matrix (K_xc)
!!                                  to which chi_(input) had already be applied to, in which case
!!                                  kernel_full contains actually chi_0*K_xc .
!! Warning !! At present ikernel=2 implies chi_{input} is equal to chi_0
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (MF, XG, YMN, GO, MG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ikernel = 0 if the kernel $K_{\lambda_i}-K_{\lambda_{i-1}}$ is diagonal
!!              in reciprocal space,
!!          = 1 otherwise.
!!          = 2 if the ALDA is implemented (Onida & Gatti technique)
!!                the product chi_0*Kxc is done in real space, FFT'd, and comes
!!                as an arg of the present routine
!!  kernel_diag = (if ikernel=0) a pointer to a (npwdiel) array containing the real part of
!!                the diagonal elements of the kernel
!!              = (if ikernel=1) an unassociated/unallocated pointer if ikernel = 1.
!!              = (if ikernel=2) a pointer to a (npwdiel) array containing the real part of
!!                the diagonal elements of the Coulomb part of the kernel
!!  kernel_full = (if ikernel=0) an unassociated/unallocated pointer
!!              = (if ikernel=1) a pointer to a (2,npwdiel,nspden,npwdiel,nspden) array
!!                containing the full kernel matrix
!!              = (if ikernel=2) a pointer to a (2,npwdiel,nspden,npwdiel,nspden) array
!!                containing chi_0*Kxc
!!  npwdiel = number of planewaves for the susceptibility matrix.
!!  nspden = number of spin-density components.
!!
!! OUTPUT
!!  (none)
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
!!
!! PARENTS
!!      acfd_dyson
!!
!! CHILDREN
!!      zgemm,zgesv,zhemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dyson_ls(ikernel,kernel_diag,kernel_full,npwdiel,nspden,susmat)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dyson_ls'
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
 integer :: info,ipw1,ipw2,isp1,isp2,isp3,ndim
!arrays
 integer,allocatable :: ipiv(:)
 real(dp),parameter :: c0(2)=(/0._dp,0._dp/),c1m(2)=(/-1._dp,0._dp/)
 real(dp),allocatable :: susmat_P_kernel(:,:,:,:,:)

!***********************************************************************

!DEBUG
!write(std_out,*)' dyson_ls : enter , ikernel=',ikernel
!ENDDEBUG

!Effective dimension of the susmat matrix, for BLAS subroutines.

 ndim = nspden*npwdiel

!Calculate susmat_P_kernel = -susmat*kernel.

 ABI_ALLOCATE(susmat_P_kernel,(2,npwdiel,nspden,npwdiel,nspden))

 if (ikernel == 1 .or. ikernel == 2) then

   if (ikernel == 1) then

     if(.true.)then

!      Case of a general kernel: treat susmat as a hermitian matrix.
!      This may be an advantage when used with a sequential BLAS library.

       call ZHEMM('l','l',ndim,ndim,c1m,susmat,ndim,&
&       kernel_full,ndim,c0,susmat_P_kernel,ndim)

     else

!      Case of a general kernel: treat susmat as a general matrix.
!      This may be an advantage when used with a parallel BLAS library.

       call ZGEMM('n','n',ndim,ndim,ndim,c1m,susmat,ndim,&
&       kernel_full,ndim,c0,susmat_P_kernel,ndim)

     end if

   else if( ikernel == 2) then

     susmat_P_kernel(:,:,:,:,:)=-kernel_full(:,:,:,:,:)

   end if

!  Case of a diagonal kernel.
 else if (ikernel == 0) then

   susmat_P_kernel(:,:,:,:,:) = zero

 end if ! ikernel=0,1, or 2

!Should add the Coulomb diagonal contribution
 if (ikernel == 0 .or. ikernel == 2)then
   do isp3 = 1,nspden
     do isp2 = 1,nspden
       do ipw2 = 1,npwdiel
         do isp1 = 1,nspden
           do ipw1 = 1,npwdiel
             susmat_P_kernel(:,ipw1,isp1,ipw2,isp2) = susmat_P_kernel(:,ipw1,isp1,ipw2,isp2)- &
&             susmat(:,ipw1,isp1,ipw2,isp3)*kernel_diag(ipw2)
           end do
         end do
       end do
     end do
   end do
 end if

!Calculate susmat_P_kernel = 1+susmat_P_kernel = 1-susmat*kernel.
 do isp1 = 1,nspden
   do ipw1 = 1,npwdiel
     susmat_P_kernel(1,ipw1,isp1,ipw1,isp1) = 1._dp+susmat_P_kernel(1,ipw1,isp1,ipw1,isp1)
   end do
 end do

!Calculate susmat = [susmat_P_kernel**-1]*susmat = [(1-susmat*kernel)**-1]*susmat.

 ABI_ALLOCATE(ipiv,(ndim))

 call ZGESV(ndim,ndim,susmat_P_kernel,ndim,ipiv,susmat,ndim,info)

 ABI_DEALLOCATE(susmat_P_kernel)
 ABI_DEALLOCATE(ipiv)

end subroutine dyson_ls

!!***
