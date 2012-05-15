!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_wfsinp_reformat
!! NAME
!! wvl_wfsinp_reformat
!!
!! FUNCTION
!! This method allocates and initialises wavefunctions with values from disk.
!! See wvl_wfsinp_scratch() or wvl_wfsinp_reformat() from other initialisation
!! routines.
!! 
!! When initialised from scratch or from disk, wvl%wfs%[h]psi comes unallocated
!! and will be allocated inside this routine.
!! When initialised from memory (reformating), wvl%wfs%[h]psi will be reallocated.
!! The projectors are also recomputed.
!!
!! The scalar arrays should be reallocated using dtset%nfft after a call to
!! this routine.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      scfcv_new
!!
!! CHILDREN
!!      copy_old_wavefunctions,deallocate_wfd,first_orthon,leave_new
!!      reformatmywaves,wrtout,wvl_projectors_free,wvl_projectors_set
!!      wvl_setboxgeometry,wvl_setngfft,wvl_wfs_free,wvl_wfs_lr_copy
!!      wvl_wfs_set,xcomm_world,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wvl_wfsinp_reformat(dtset, mpi_enreg, psps,&
     & rprimd, wvl, xred, xred_old)

 use m_profiling

  use defs_basis
  use defs_datatypes
  use defs_abitypes
  use defs_wvltypes
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only : copy_old_wavefunctions, reformatmywaves, first_orthon, deallocate_wfd, wavefunctions_descriptors
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_wfsinp_reformat'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
 use interfaces_43_wvl_wrappers
 use interfaces_51_manage_mpi
 use interfaces_57_iovars
 use interfaces_62_wvl_wfs
!End of the abilint section

  implicit none

!Arguments ------------------------------------
  type(dataset_type), intent(inout)      :: dtset
  type(MPI_type), intent(inout)          :: mpi_enreg
  type(pseudopotential_type), intent(in) :: psps
  type(wvl_data), intent(inout)          :: wvl
  real(dp), intent(inout)                :: rprimd(3,3)
  real(dp), intent(inout)                :: xred_old(3, dtset%natom)
  real(dp), intent(inout)                :: xred(3, dtset%natom)

!Local variables-------------------------------
  integer                  :: nSize_old(3)
  real(dp)                 :: hgrid_old(3)
  real(dp), allocatable    :: xcart(:,:), xcart_old(:,:)
  real(dp), pointer        :: psi_old(:), eigen_old(:)
#if defined HAVE_DFT_BIGDFT
  integer :: comm,me,nproc
  type(wavefunctions_descriptors) :: keys_old
#endif
  character(len=500)       :: message

! *********************************************************************

 write(message, '(a,a)' ) ch10,&
& ' wvl_wfsinp_reformat: reformat the wavefunctions.'
 call wrtout(std_out, message, 'COLL')

#if defined HAVE_DFT_BIGDFT

 call xcomm_world(mpi_enreg,comm,myrank=me,mysize=nproc)

!Convert input xred_old (reduced coordinates) to xcart_old (cartesian)
 ABI_ALLOCATE(xcart_old,(3, dtset%natom))
 call xredxcart(dtset%natom, 1, rprimd, xcart_old, xred_old)

!Copy current to old.
 ABI_ALLOCATE(eigen_old,(wvl%wfs%orbs%norb))
 eigen_old = wvl%wfs%orbs%eval
 hgrid_old = wvl%descr%h
 call copy_old_wavefunctions(nproc, wvl%wfs%orbs, &
& wvl%descr%Glr%d%n1, wvl%descr%Glr%d%n2, wvl%descr%Glr%d%n3, &
& wvl%wfs%Glr%wfd, wvl%wfs%psi, nSize_old(1), nSize_old(2), nSize_old(3), &
& keys_old, psi_old)

!We deallocate the previous projectors.
 call wvl_projectors_free(wvl%projectors)

!Deallocate old wavefunctions
 call wvl_wfs_free(wvl%wfs)

!We change the box geometry.
 call wvl_setBoxGeometry(me, dtset%prtvol, psps%gth_params%radii_cf, rprimd, xred, &
& wvl%descr, dtset%wvl_crmult, dtset%wvl_frmult)
 call wvl_setngfft(dtset%ixc, dtset%mgfft, mpi_enreg, dtset%natom, dtset%nfft, &
& dtset%ngfft, dtset%nsppol, psps, rprimd, wvl%descr, dtset%wvl_crmult, &
& dtset%wvl_frmult, xred)

!We copy the geometry structure.
 call wvl_wfs_lr_copy(wvl%wfs, wvl%descr)
!Reallocate them with new size.
 call wvl_wfs_set(dtset%fixmom, dtset%kpt, me, dtset%natom, sum(dtset%nband), dtset%nkpt, &
& nproc, dtset%nspinor, dtset%nsppol, dtset%nwfshist, dtset%occ_orig, psps, rprimd, &
& wvl%wfs, dtset%wtk, wvl%descr, dtset%wvl_crmult, dtset%wvl_frmult, xred)

!Recopy old eval for precond.
 wvl%wfs%orbs%eval = eigen_old
 ABI_DEALLOCATE(eigen_old)

!We allocate psi.
 ABI_ALLOCATE(wvl%wfs%psi,(wvl%wfs%orbs%npsidim))
 write(message, '(a,a,a,a,I0)' ) ch10, &
& ' wvl_wfsinp_reformat: allocate wavefunctions,', ch10, &
& '  size of the compressed array per proc: ', &
& product(shape(wvl%wfs%psi))
 call wrtout(std_out,message,'COLL')

!Convert input xred (reduced coordinates) to xcart (cartesian)
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 call xredxcart(dtset%natom, 1, rprimd, xcart, xred)

!We transfer the old wavefunctions to the new ones.
 call reformatmywaves(me, wvl%wfs%orbs, wvl%descr%atoms, &
& hgrid_old(1), hgrid_old(2), hgrid_old(3), nSize_old(1), nSize_old(2), &
& nSize_old(3), xcart_old, keys_old, psi_old, wvl%descr%h(1), wvl%descr%h(2), &
& wvl%descr%h(3), wvl%descr%Glr%d%n1, wvl%descr%Glr%d%n2, wvl%descr%Glr%d%n3, xcart, &
& wvl%wfs%Glr%wfd, wvl%wfs%psi)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xcart_old)

!We free the old descriptors and arrays.
 ABI_DEALLOCATE(psi_old)
 call deallocate_wfd(keys_old, "wvl_wfsinp_reformat")

!Reallocate projectors for the new positions.
 call wvl_projectors_set(me, dtset%natom, wvl%projectors, psps, rprimd, &
& wvl%wfs, wvl%descr, dtset%wvl_frmult, xred)

!Orthogonilise new wavefunctions.
 call first_orthon(me, nproc, wvl%wfs%orbs, wvl%wfs%Glr%wfd, wvl%wfs%comms, &
& wvl%wfs%psi, wvl%wfs%hpsi, wvl%wfs%psit, wvl%descr%orthpar)

#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_wfsinp_reformat: BUG -',ch10,&
& '  BigDFT is not compiled. Use --enable-bigdft during configure.'
 call wrtout(std_out, message, 'COLL')
 call leave_new('COLL')
#endif

end subroutine wvl_wfsinp_reformat
!!***
