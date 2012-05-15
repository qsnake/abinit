!{\src2tex{textfont=tt}}
!!****f* ABINIT/PSolver_kernel
!! NAME
!! PSolver_kernel
!!
!! FUNCTION
!! Build, get or free the kernel matrix used by the Poisson solver to compute the
!! the convolution between 1/r and rho. The kernel is a saved variable. If
!! this routine is called for building while a kernel already exists, it is not
!! recomputed if all parameters (grid step and data size) are unchanged. Otherwise
!! the kernel is freed and recompute again. The build action has a returned variable
!! which is a pointer on the kernel. The get action also returns the kernel, or
!! NULL if none has been associated.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  iaction=0 to free all kernel allocated array,
!!          1 to compute the kernel (parallel case),
!!          2 to get it (parallel case),
!!          3 to compute the kernel (sequential),
!!          4 to get the sequential kernel.
!!  mpi_enreg=MPI-parallelisation information.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! OUTPUT
!!  kernel= associated kernel on build (iaction = 1) and get action (iaction = 2).
!!
!! PARENTS
!!      afterscfloop,gstate,mklocl_realspace,mklocl_wavelets,psolver_hartree
!!      psolver_rhohxc,wvl_vtorho,wvl_wfsinp_scratch
!!
!! CHILDREN
!!      createkernel,leave_new,memocc,wrtout,xcomm_world
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine PSolver_kernel(dtset, iaction, kernel_array, mpi_enreg, rprimd, wvl)
  
  use defs_basis
  use defs_abitypes
  use defs_wvltypes
  use m_profiling
#if defined HAVE_DFT_BIGDFT
  use poisson_solver, only: createKernel
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'PSolver_kernel'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
!End of the abilint section

  implicit none

!Arguments ------------------------------------
  !scalars
  integer,intent(in) :: iaction
  type(dataset_type),intent(in) :: dtset
  type(MPI_type),intent(in) :: mpi_enreg
  type(wvl_internal_type), intent(in) :: wvl
  !arrays
  real(dp), pointer :: kernel_array(:)
  real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------
  !scalars
  !arrays
  integer, save :: kernel_scfOrder 
  integer, save :: icoulomb
  character(len = 1) :: bnd_code
  integer, save :: data_size(3) = (/ -2, -2, -2 /) 
  real(dp), save :: kernel_hgrid(3)   ! Grid step used when creating the kernel.
  real(dp), pointer, save :: kernel(:), kernelseq(:) 

  integer :: comm,me,nproc,i_all,i_stat
  character(len=1) :: geocode
  !scalars
  real(dp) :: hgrid(3)
  character(len=500) :: message
  !arrays
  integer :: current_size(3)

! *************************************************************************

 call xcomm_world(mpi_enreg,comm,myrank=me,mysize=nproc)

 if (dtset%icoulomb == 1) then
   geocode='F'
 else if (dtset%icoulomb == 2) then
   geocode='S'
 end if

!Initialise kernel_array pointer.
 if (maxval(data_size) == -2) then
   nullify(kernel)
   nullify(kernelseq)
 end if

!If iaction == 0, we free the kernel.
 if (iaction == 0) then
   if (associated(kernel)) then
     write(message, "(A)") "Psolver_kernel() : deallocating kernel..."
     call wrtout(std_out, message,'COLL')
     
     i_all=-product(shape(kernel))*kind(kernel)
     ABI_DEALLOCATE(kernel)
     i_stat = ABI_ALLOC_STAT
     call memocc(i_stat,i_all,'kernel',"psolver_kernel")
   end if
   if (associated(kernelseq)) then
     write(message, "(A)") "Psolver_kernel() : deallocating kernelseq..."
     call wrtout(std_out, message,'COLL')

     i_all=-product(shape(kernelseq))*kind(kernelseq)
     ABI_DEALLOCATE(kernelseq)
     i_stat = ABI_ALLOC_STAT
     call memocc(i_stat,i_all,'kernelseq',"psolver_kernel")
   end if
   data_size = (/ -1, -1, -1 /)
   return
 end if


!Action is build or get. We check the sizes before doing anything else.

!Get the size depending on wavelets calculations or not
 if (dtset%usewvl == 0) then
   hgrid(1) = rprimd(1, 1) / dtset%ngfft(1)
   hgrid(2) = rprimd(2, 2) / dtset%ngfft(2)
   hgrid(3) = rprimd(3, 3) / dtset%ngfft(3)

   current_size(:) = dtset%ngfft(1:3)
 else
   hgrid(:) = 0.5d0 * wvl%h(:)
#if defined HAVE_DFT_BIGDFT
   current_size(1:3) = (/ wvl%Glr%d%n1i, wvl%Glr%d%n2i, wvl%Glr%d%n3i /)
#endif
 end if


!If iaction == 2, we get the kernel.
 if (iaction == 2 .or. iaction == 4) then
   if (icoulomb     == dtset%icoulomb  .and. &
&   data_size(1)    == current_size(1) .and. &
&   data_size(2)    == current_size(2) .and. &
&   data_size(3)    == current_size(3) .and. &
&   kernel_hgrid(1) == hgrid(1)        .and. &
&   kernel_hgrid(2) == hgrid(2)        .and. &
&   kernel_hgrid(3) == hgrid(3)        .and. &
&   kernel_scfOrder == dtset%nscforder) then
     if (associated(kernel) .and. iaction == 2) then
       kernel_array => kernel
       return
     end if
     if (associated(kernelseq) .and. iaction == 4) then
       kernel_array => kernelseq
       return
     end if
   end if
!  The kernel must be rebuilt, we continue in the routine.
 end if


!Build action, we do some checks before.

!Check if the box is cubic. Can't work if the grid size is not the same
!in x, y and z directions.
 if (rprimd(1, 2) /= 0.d0 .or. rprimd(1, 3) /= 0.d0 .or. rprimd(2, 3) /= 0.d0) then
   write(message, '(a,a,a,a,9F6.2,a)' ) ch10,&
&   ' Psolver_kernel: BUG -',ch10,&
&   '  box geometry is not a parallelepipede', rprimd, &
&   '  .'
   call wrtout(std_out, message, 'COLL')
   call leave_new('COLL')
 end if

!allocate(rhopot(ngfft(1), ngfft(2), ngfft(3)))

!Compute a new kernel if grid size has changed or if the kernel
!has never been computed.
 if ((iaction == 1 .and. .not. associated(kernel))    .or. &
& (iaction == 3 .and. .not. associated(kernelseq)) .or. &
& icoulomb        /= dtset%icoulomb  .or. &
& data_size(1)    /= current_size(1) .or. &
& data_size(2)    /= current_size(2) .or. &
& data_size(3)    /= current_size(3) .or. &
& kernel_hgrid(1) /= hgrid(1)        .or. &
& kernel_hgrid(2) /= hgrid(2)        .or. &
& kernel_hgrid(3) /= hgrid(3)        .or. &
& kernel_scfOrder /= dtset%nscforder) then
   write(message, "(A,A,A,3I6)") "Psolver_kernel() : building kernel...", ch10, &
&   " | data dimensions:", current_size
   call wrtout(std_out, message, 'COLL')
!  This routine will allocate the pointer kernel with dimensions
!  depending on the number of processors and the datra size.
   if (dtset%icoulomb == 0) then
!    The kernel is built with 'P'eriodic boundary counditions.
     bnd_code = 'P'
   else if (dtset%icoulomb == 1) then
!    The kernel is built with 'F'ree boundary counditions.
     bnd_code = 'F'
   else if (dtset%icoulomb == 2) then
!    The kernel is built with 'S'urface boundary counditions.
     bnd_code = 'S'
   end if
#if defined HAVE_DFT_BIGDFT

   if (iaction == 1 .or. iaction == 2) then
     if (associated(kernel)) then
       i_all=-product(shape(kernel))*kind(kernel)
       ABI_DEALLOCATE(kernel)
       i_stat = ABI_ALLOC_STAT
       call memocc(i_stat,i_all,'kernel',"psolver_kernel")
     end if
     call createKernel(me, nproc, bnd_code, current_size(1), current_size(2), &
&     current_size(3), hgrid(1), hgrid(2), hgrid(3), dtset%nscforder, &
&     kernel, quiet = "no ")
   end if

   if (iaction == 3 .or. iaction == 4) then
     if (associated(kernelseq)) then
       i_all=-product(shape(kernelseq))*kind(kernelseq)
       ABI_DEALLOCATE(kernelseq)
       i_stat = ABI_ALLOC_STAT
       call memocc(i_stat,i_all,'kernelseq',"psolver_kernel")
     end if
     call createKernel(0, 1, bnd_code, current_size(1), current_size(2), &
&     current_size(3), hgrid(1), hgrid(2), hgrid(3), dtset%nscforder, &
&     kernelseq, quiet = "no ")
   end if

#else
   write(message, '(a,a,a,a)' ) ch10,&
&   ' Psolver_kernel: BUG -',ch10,&
&   '  BigDFT is not compile. Use --enable-bigdft during configure.'
   call wrtout(std_out, message, 'COLL')
   call leave_new('COLL')
#endif

!  Storing variables which were used to make the kernel
   icoulomb        = dtset%icoulomb
   data_size(:)    = current_size(:)
   kernel_hgrid(:) = hgrid(:)
   kernel_scfOrder = dtset%nscforder
 end if
 
 if (iaction == 1 .or. iaction == 2) kernel_array => kernel
 if (iaction == 3 .or. iaction == 4) kernel_array => kernelseq

end subroutine PSolver_kernel
!!***
