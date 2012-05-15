!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_grid
!! NAME
!!  initmpi_grid
!!
!! FUNCTION
!!  Initializes the MPI information for the grid:
!!    * 2D if parallization FFT/BAND (!MPI paral_kgb)
!!    * 3D if parallization KPT/FFT/BAND (paral_kgb & MPI)
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2012 ABINIT group
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt.
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! TODO
!!
!! PARENTS
!!      initmpi_fft,invars1
!!
!! CHILDREN
!!      mpi_cart_coords,mpi_cart_create,mpi_cart_sub,mpi_comm_rank
!!      mpi_comm_size,xmpi_abort
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine initmpi_grid(mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi

#if defined HAVE_MPI2 && !defined FC_G95
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_grid'
!End of the abilint section

 implicit none
#if defined HAVE_MPI1 || (defined HAVE_MPI && defined FC_G95)
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(MPI_type),intent(inout) :: mpi_enreg
!Local variables-------------------------------
#if defined HAVE_MPI
 !Variables introduced for MPI version
 integer :: ierr

 !Variables introduced for the bandFFT version
 integer :: np_test,nproc,spacecomm
 logical :: reorder
 logical, allocatable :: periode(:), keepdim(:)
#endif

! *********************************************************************

!DEBUG
!write(std_out,*)' initmpi_grid : enter'
!ENDDEBUG

!Default parameters, for safety of sequential use
 mpi_enreg%me_fft=0
 mpi_enreg%me_kpt=0
 mpi_enreg%me_band=0
 mpi_enreg%me_spin=0
 mpi_enreg%paral_spin=0

#if defined HAVE_MPI

!TEST THE PARAMETERS OF THE 4D GRID
!========================================

 if (mpi_enreg%paral_img==0) then
   nproc=mpi_enreg%nproc
   spacecomm=mpi_enreg%world_comm
 else
   nproc=mpi_enreg%nproc_one_img
   spacecomm=mpi_enreg%comm_one_img
 end if

 if (mpi_enreg%nproc_spin>1) mpi_enreg%paral_spin=1

 if(mpi_enreg%nproc_fft*  &
& mpi_enreg%nproc_band* &
& mpi_enreg%nproc_kpt* &
& mpi_enreg%nproc_spin /= nproc)then

   write(std_out,'(8a,i5,a,i5,a,i5,a,i5,a,i5)') ch10,&
&   ' initmpi_grid : WARNING -',ch10,&
&   '  The number of band*FFT*kpt*spinor processors, npband*npfft*npkpt*npspinor should be',ch10,&
&   '  equal to the total number of processors, nproc.',ch10,&
&   '  However, npband   =',mpi_enreg%nproc_band,&
&   '           npfft    =',mpi_enreg%nproc_fft,&
&   '           npkpt    =',mpi_enreg%nproc_kpt,&
&   '           npspinor =',mpi_enreg%nproc_spin,&
&   '  and nproc=', nproc
!  call leave_new('PERS')
 end if

 write(std_out,*) 'npfft, npband, npspinor and npkpt',&
& mpi_enreg%nproc_fft,mpi_enreg%nproc_band, &
& mpi_enreg%nproc_spin,mpi_enreg%nproc_kpt


!CREATE THE 4D GRID
!==================================================

!Fake values for null communicator
 if (nproc==0) then
   mpi_enreg%nproc_fft    = 0
   mpi_enreg%nproc_band   = 0
   mpi_enreg%nproc_kpt    = 0
   mpi_enreg%nproc_spin   = 0
   mpi_enreg%me_fft       = 0
   mpi_enreg%me_band      = 0
   mpi_enreg%me_kpt       = 0
   mpi_enreg%me_spin      = 0
   mpi_enreg%me_cart_2d   = 0
   mpi_enreg%me_cart_3d   = 0
   mpi_enreg%comm_fft     = MPI_COMM_NULL
   mpi_enreg%comm_band    = MPI_COMM_NULL
   mpi_enreg%comm_kpt     = MPI_COMM_NULL
   mpi_enreg%comm_spin    = MPI_COMM_NULL
   mpi_enreg%comm_bandspin= MPI_COMM_NULL
   mpi_enreg%comm_spinfft = MPI_COMM_NULL
   mpi_enreg%commcart     = MPI_COMM_NULL
   mpi_enreg%commcart_3d  = MPI_COMM_NULL
   mpi_enreg%commcart_4d  = MPI_COMM_NULL
   mpi_enreg%bandpp       = 1
   return
 end if

 mpi_enreg%dimcart=4
 ABI_ALLOCATE(mpi_enreg%sizecart,(mpi_enreg%dimcart))
!valgrind claims this is not deallocated in test v5/72 Can someone knowledgable check?
 ABI_ALLOCATE(periode           ,(mpi_enreg%dimcart))
 ABI_ALLOCATE(mpi_enreg%coords  ,(mpi_enreg%dimcart))

 mpi_enreg%sizecart(1) = mpi_enreg%nproc_fft
 mpi_enreg%sizecart(2) = mpi_enreg%nproc_band
 mpi_enreg%sizecart(3) = mpi_enreg%nproc_kpt
 mpi_enreg%sizecart(4) = mpi_enreg%nproc_spin

 periode(:)=.false.
 reorder   =.false.

!Create the global cartesian 4D- communicator
 call MPI_CART_CREATE(spacecomm,mpi_enreg%dimcart,&
& mpi_enreg%sizecart,periode,reorder,mpi_enreg%commcart_4d,ierr)
 if (ierr /= MPI_SUCCESS ) call xmpi_abort(xmpi_world,ierr)

!Find the index and coordinates of the current processor
 call MPI_COMM_RANK(mpi_enreg%commcart_4d, mpi_enreg%me_cart_4d, ierr)
 call MPI_CART_COORDS(mpi_enreg%commcart_4d, mpi_enreg%me_cart_4d,&
& mpi_enreg%dimcart,mpi_enreg%coords,ierr)
 if (ierr /= MPI_SUCCESS ) call xmpi_abort(xmpi_world,ierr)

 ABI_ALLOCATE(keepdim,(mpi_enreg%dimcart))

!Create the communicator for fft distribution
 keepdim(1)=.true.
 keepdim(2)=.false.
 keepdim(3)=.false.
 keepdim(4)=.false.
 call MPI_CART_SUB(mpi_enreg%commcart_4d, keepdim, mpi_enreg%comm_fft,ierr)
 if (ierr /= MPI_SUCCESS ) call xmpi_abort(xmpi_world,ierr)

!Create the communicator for band distribution
 keepdim(1)=.false.
 keepdim(2)=.true.
 keepdim(3)=.false.
 keepdim(4)=.false.
 call MPI_CART_SUB(mpi_enreg%commcart_4d, keepdim, mpi_enreg%comm_band,ierr)
 if (ierr /= MPI_SUCCESS ) call xmpi_abort(xmpi_world,ierr)

!Create the communicator for kpt distribution
 keepdim(1)=.false.
 keepdim(2)=.false.
 keepdim(3)=.true.
 keepdim(4)=.false.
 call MPI_CART_SUB(mpi_enreg%commcart_4d, keepdim, mpi_enreg%comm_kpt,ierr)
 if (ierr /= MPI_SUCCESS ) call xmpi_abort(xmpi_world,ierr)

!Create the communicator for spinor distribution
 keepdim(1)=.false.
 keepdim(2)=.false.
 keepdim(3)=.false.
 keepdim(4)=.true.
 call MPI_CART_SUB(mpi_enreg%commcart_4d, keepdim, mpi_enreg%comm_spin,ierr)
 if (ierr /= MPI_SUCCESS ) call xmpi_abort(xmpi_world,ierr)

!Create the communicator for band-spinor distribution
 keepdim(1)=.false.
 keepdim(2)=.true.
 keepdim(3)=.false.
 keepdim(4)=.true.
 call MPI_CART_SUB(mpi_enreg%commcart_4d, keepdim, mpi_enreg%comm_bandspin,ierr)
 if (ierr /= MPI_SUCCESS ) call xmpi_abort(xmpi_world,ierr)

!Create the communicator for fft-spinor distribution
 keepdim(1)=.true.
 keepdim(2)=.false.
 keepdim(3)=.false.
 keepdim(4)=.true.
 call MPI_CART_SUB(mpi_enreg%commcart_4d, keepdim, mpi_enreg%comm_spinfft,ierr)
 if (ierr /= MPI_SUCCESS ) call xmpi_abort(xmpi_world,ierr)

!Create the communicator for fft-band distribution
 keepdim(1)=.true.
 keepdim(2)=.true.
 keepdim(3)=.false.
 keepdim(4)=.false.
 call MPI_CART_SUB(mpi_enreg%commcart_4d, keepdim, mpi_enreg%commcart,ierr)
 if (ierr /= MPI_SUCCESS ) call xmpi_abort(xmpi_world,ierr)
 call MPI_COMM_RANK(mpi_enreg%commcart, mpi_enreg%me_cart_2d, ierr)
 if (ierr /= MPI_SUCCESS ) call xmpi_abort(xmpi_world,ierr)

!Create the communicator for fft-band-spinor distribution
 keepdim(1)=.true.
 keepdim(2)=.true.
 keepdim(3)=.false.
 keepdim(4)=.true.
 call MPI_CART_SUB(mpi_enreg%commcart_4d, keepdim, mpi_enreg%commcart_3d,ierr)
 if (ierr /= MPI_SUCCESS ) call xmpi_abort(xmpi_world,ierr)
 call MPI_COMM_RANK(mpi_enreg%commcart_3d, mpi_enreg%me_cart_3d, ierr)
 if (ierr /= MPI_SUCCESS ) call xmpi_abort(xmpi_world,ierr)

!Writings
 call MPI_COMM_SIZE(mpi_enreg%comm_fft,np_test, ierr)
 write(std_out,*) 'mpi_enreg%sizecart(1),np_fft   =' ,mpi_enreg%sizecart(1),np_test
 call MPI_COMM_SIZE(mpi_enreg%comm_band,np_test, ierr)
 write(std_out,*) 'mpi_enreg%sizecart(2),np_band  =',mpi_enreg%sizecart(2),np_test
 call MPI_COMM_SIZE(mpi_enreg%comm_kpt,np_test, ierr)
 write(std_out,*) 'mpi_enreg%sizecart(3),np_kpt   =',mpi_enreg%sizecart(3),np_test
 call MPI_COMM_SIZE(mpi_enreg%comm_spin,np_test, ierr)
 write(std_out,*) 'mpi_enreg%sizecart(4),np_spinor=',mpi_enreg%sizecart(4), np_test

!Define the correspondance with the fft
 mpi_enreg%me_fft  = mpi_enreg%coords(1)
 mpi_enreg%me_band = mpi_enreg%coords(2)
 mpi_enreg%me_kpt  = mpi_enreg%coords(3)
 mpi_enreg%me_spin = mpi_enreg%coords(4)

!WRITE THE PROCESSORS COORDINATES IN THE 3D GRID
!=====================================================
 write(std_out,*) 'in initmpi_grid : me_fft, me_band, me_spin , me_kpt are ',&
& mpi_enreg%me_fft,mpi_enreg%me_band,mpi_enreg%me_spin, mpi_enreg%me_kpt

 ABI_DEALLOCATE(periode)
 ABI_DEALLOCATE(keepdim)

#endif

!DEBUG
!write(std_out,*)' initmpi_grid : exit'
!ENDDEBUG

end subroutine initmpi_grid
!!***




