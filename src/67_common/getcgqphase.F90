!{\src2tex{textfont=tt}}
!!****f* ABINIT/getcgqphase
!! NAME
!! getcgqphase
!!
!! FUNCTION
!! extract phases from wave functions, to cancel contributions to gkk matrix elements
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  timrev = flag for use of time reversal symmetry
!!  cg = input wavefunctions
!!  mcg = dimension of cg = nspinor*mband*mpw*mkmem
!!  cgq = input wavefunctions at k+q
!!  mcgq = dimension of cgq = nspinor*mband*mpw*mkmem
!!  mpi_enreg = datastructure for mpi communication
!!  nkpt_rbz = number of k-points in reduced zone for present q point
!!  npwarr = array of numbers of plane waves for each k-point
!!  npwar1 = array of numbers of plane waves for each k+q point 
!!
!! OUTPUT
!!  phasecg = phase of different wavefunction products <k,n | k+q,n'>
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      smatrix
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"



subroutine getcgqphase(dtset, timrev, cg,  mcg,  cgq, mcgq, mpi_enreg, &
&    nkpt_rbz, npwarr, npwar1, phasecg)

 use m_profiling
 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getcgqphase'
 use interfaces_67_common, except_this_one => getcgqphase
!End of the abilint section

 implicit none

!Arguments -------------------------------
 ! scalars 
 integer, intent(in) :: mcg, mcgq, timrev
 integer, intent(in) :: nkpt_rbz
 type(dataset_type), intent(in) :: dtset

 ! arrays
 integer, intent(in) :: npwarr(nkpt_rbz)
 integer, intent(in) :: npwar1(nkpt_rbz)
 real(dp), intent(in) :: cg(2,mcg)
 real(dp), intent(in) :: cgq(2,mcgq)
 type(MPI_type), intent(in) :: mpi_enreg
 real(dp),intent(out) :: phasecg(2, dtset%mband*dtset%mband*nkpt_rbz&
&     *dtset%nsppol)

!Local variables -------------------------
 ! local vars
 integer :: icg, icgq, isppol, ikpt, ipw
 integer :: istate, iband1, iband2
 integer :: mcg_k, npw_k, npw_q
 integer :: usepaw
 integer :: ddkflag, itrs, job, maxbd, mcg1_k, minbd, shiftbd
 real(dp) :: normsmat

 integer, allocatable :: sflag_k(:)
 integer, allocatable :: pwind_k(:)

 real(dp) :: cg1_dummy(1,1)
 real(dp) :: smat_inv_dummy(1,1,1)
 real(dp) :: smat_k_paw_dummy(1,1,1)
 real(dp) :: dtm_k_dummy(2)
 real(dp), allocatable :: smat_k(:,:,:)
 real(dp), allocatable :: pwnsfac_k(:,:)

! *********************************************************************

 ABI_ALLOCATE(smat_k,(2,dtset%mband,dtset%mband))
 ABI_ALLOCATE(sflag_k,(dtset%mband))

!dummy use of mpi_enreg so abirules stops complaining.
 icg = mpi_enreg%me
!dummy use of timrev so abirules stops complaining.
 icg = timrev

!!MPI data for future use
!call xcomm_init(mpi_enreg,spaceComm_distrb)
!call xme_init(mpi_enreg,me_distrb)
!nproc_distrb=xcomm_size(spaceComm_distrb)

!make trivial association of G vectors: we just want <psi_k| psi_k+q>
!TODO: check this is correct wrt arrangement of kg vectors for k+q
!looks ok : usually made in initberry, from the translations associated
!to the symops, scalar product with the G vectors. The symop is the one
!used to go from the irreducible k to the full zone k. In present context
!we should be using only the reduced zone, and anyhow have the same k-grid
!for the gkk matrix elements and for the cg here...
 ABI_ALLOCATE(pwind_k,(dtset%mpw))
 ABI_ALLOCATE(pwnsfac_k,(4,dtset%mpw))
 do ipw = 1, dtset%mpw
   pwind_k(ipw) = ipw
   pwnsfac_k(1,ipw) = one
   pwnsfac_k(2,ipw) = zero
   pwnsfac_k(3,ipw) = one
   pwnsfac_k(4,ipw) = zero
 end do

!flags for call to smatrix
 usepaw = 0 ! for now
 ddkflag = 0
 itrs = 0
 job = 0
 maxbd = 1
 mcg1_k = 1
 minbd = 1
 shiftbd = 1

!from overlap matrix for each wavefunction, extract phase
 icg = 0
 icgq = 0
 istate = 0

 phasecg = zero
 phasecg(1,:) = one
 do isppol = 1, dtset%nsppol
   do ikpt = 1, nkpt_rbz
     npw_k = npwarr(ikpt)
     npw_q= npwar1(ikpt)
     mcg_k = dtset%mpw*dtset%mband*dtset%nspinor

!    TODO: question: are the k-points in the ibz correctly ordered in cg and cgq? if not the icg below have to be adapted.
!    TODO: if mkmem is less than nkpt_rbz we need to parallelize this, or else read in cg from disk...
     sflag_k = 0 ! make sure all elements are calculated
     smat_k = zero

     call smatrix(cg, cgq, cg1_dummy, ddkflag, dtm_k_dummy, icg, icgq,&
&     itrs, job, maxbd, mcg, mcgq, mcg1_k, minbd,dtset%mpw, dtset%mband,&
&     npw_k, npw_q, dtset%nspinor, pwind_k, pwnsfac_k, sflag_k, shiftbd,&
&     smat_inv_dummy, smat_k, smat_k_paw_dummy, usepaw)

     icg  = icg  + npw_k*dtset%nspinor*dtset%nband(ikpt)
     icgq = icgq + npw_q*dtset%nspinor*dtset%nband(ikpt)


     do iband1 = 1, dtset%nband(ikpt)
       do iband2 = 1, dtset%nband(ikpt)
         istate = istate + 1
!        normalise the overlap matrix element to get just the phase difference phi_k - phi_k+q
         normsmat = sqrt(smat_k(1,iband2, iband1)**2 &
&         + smat_k(2,iband2, iband1)**2)
         if (normsmat > tol12) then
           phasecg(:, istate) = smat_k(:,iband2, iband1) / normsmat
!          NOTE: 21/9/2011 these appear to be always 1, i, or -i, to within 1.e-5 at worst!
         end if
       end do
     end do
   end do
 end do

!eventually do an mpi allreduce over the k-points for phasecg
!

 ABI_DEALLOCATE(sflag_k)
 ABI_DEALLOCATE(smat_k)

end subroutine getcgqphase
!!***
