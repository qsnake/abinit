!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_setngfft
!! NAME
!! wvl_setngfft
!!
!! FUNCTION
!! When wavelets are used, the FFT grid is used to store potentials and
!! density. The size of the grid takes into account the two resolution in wavelet
!! description and also the distribution over processor in the parallel case.
!!
!! The FFT grid is not in strict terms an FFT grid but rather a real space grid.
!! Its dimensions are not directly compatible with FFTs. This is not relevant
!! when using the wavelet part of the code and in the Poisson solver the arrays
!! are extended to match FFT dimensions internally. But for other parts of the
!! code, this must be taken into account.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization (description of the
!!            density and potentials scatterring is allocated and updated).
!!  dtset <type(dataset_type)>=the FFT grid is changed.
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      createdenspotdescriptors,leave_new,wrtout,xcomm_world,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_setngfft(ixc, mgfft, mpi_enreg, natom, nfft, ngfft, nsppol, psps, rprimd, &
     & wvl, wvl_crmult, wvl_frmult, xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi,only : xcomm_size
#if defined HAVE_DFT_BIGDFT
  use Poisson_Solver
  use BigDFT_API, only : createDensPotDescriptors
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_setngfft'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in)                    :: natom, ixc, nsppol
 integer, intent(out)                   :: mgfft, nfft
 real(dp), intent(in)                   :: wvl_crmult, wvl_frmult
 type(MPI_type),intent(inout)           :: mpi_enreg
 type(pseudopotential_type),intent(in)  :: psps
 type(wvl_internal_type), intent(inout) :: wvl
!arrays
 integer, intent(out)                   :: ngfft(13)
 real(dp), intent(inout)                :: rprimd(3,3)
 real(dp), intent(inout)                :: xred(3, natom)

!Local variables-------------------------------
!scalars
 integer :: comm,density_start,me,ngfft3_density,ngfft3_potential
 integer :: nproc,potential_shift
 character(len=500) :: message
 real(dp), allocatable :: xcart(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' wvl_setngfft : enter '
!write(std_out,*)' associated(mpi_enreg%nscatterarr)=',associated(mpi_enreg%nscatterarr)
!stop
!ENDDEBUG

#if defined HAVE_DFT_BIGDFT
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_setngfft : Changing the FFT grid definition.'
 call wrtout(std_out,message,'COLL')

!Change nfft and ngfft
!We use the routine given in the poisson solver part.

 call xcomm_world(mpi_enreg,comm,myrank=me,mysize=nproc)

 ABI_ALLOCATE(xcart,(3, natom))
 call xredxcart(natom, 1, rprimd, xcart, xred)

!number of planes for the density
!mpi_enreg%nscatterarr(jproc, 1) = ngfft3_density
!number of planes for the potential
!mpi_enreg%nscatterarr(jproc, 2) = ngfft3_potential
!starting offset for the potential
!mpi_enreg%nscatterarr(jproc, 3) = density_start + potential_shift - 1
!GGA XC shift between density and potential
!mpi_enreg%nscatterarr(jproc, 4) = potential_shift
 call createDensPotDescriptors(me, nproc, wvl%atoms, wvl%Glr%d, &
& wvl%h(1) / 2., wvl%h(2) / 2., wvl%h(3) / 2., xcart, &
& wvl_crmult, wvl_frmult, psps%gth_params%radii_cf, nsppol, 'D', &
& ixc, "DBL", ngfft3_density, ngfft3_potential, mpi_enreg%ngfft3_ionic, &
& potential_shift, density_start, mpi_enreg%nscatterarr, mpi_enreg%ngatherarr, &
& wvl%rhodsc)

 ABI_DEALLOCATE(xcart)

!Now ngfft will use the density definition (since the potential size
!is always smaller than the density one).
 ngfft(1) = wvl%Glr%d%n1i
 ngfft(2) = wvl%Glr%d%n2i
 ngfft(3) = mpi_enreg%nscatterarr(me, 1)

 nfft = product(ngfft(1:3))
!Set up fft array dimensions ngfft(4,5,6) to avoid cache conflicts
!Code paste from getng()
 ngfft(4) = 2 * (ngfft(1) / 2) + 1
 ngfft(5) = 2 * (ngfft(2) / 2) + 1
 ngfft(6) = ngfft(3)
 if (nproc == 0) then
   ngfft(9)  = 0    ! paral_fft
   ngfft(10) = 1    ! nproc_fft
   ngfft(11) = 0    ! me_fft
   ngfft(12) = 0    ! n2proc
   ngfft(13) = 0    ! n3proc
 else
   ngfft(9)  = 1    ! paral_fft
   ngfft(10) = mpi_enreg%nproc_fft
   ngfft(11) = mpi_enreg%me_fft
   ngfft(12) = ngfft(2)
   ngfft(13) = ngfft(3)
 end if
 write(message, '(a,3I12)' ) &
& '  | ngfft(1:3) is now:    ', ngfft(1:3)
 call wrtout(std_out,message,'COLL')
 write(message, '(a,3I12)' ) &
& '  | ngfft(4:6) is now:    ', ngfft(4:6)
 call wrtout(std_out,message,'COLL')

!Set mgfft
 mgfft= max(ngfft(1), ngfft(2), ngfft(3))

#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_setngfft : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif
end subroutine wvl_setngfft
!!***
