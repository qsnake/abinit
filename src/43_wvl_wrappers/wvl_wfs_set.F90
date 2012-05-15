!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_wfs_set
!! NAME
!! wvl_wfs_set
!!
!! FUNCTION
!! Compute the access keys for the wavefunctions when the positions
!! of the atoms are given.
!!
!! For memory occupation optimisation reasons, the wavefunctions are not allocated
!! here. See the initialisation routines wvl_wfsinp_disk(), wvl_wfsinp_scratch()
!! and wvl_wfsinp_reformat() to do it. After allocation, use wvl_wfs_free()
!! to deallocate all stuff (descriptors and arrays).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=internal variables used by wavelets, describing
!!   | wvl_internal=desciption of the wavelet box.
!!   | natom=number of atoms.
!!  mpi_enreg=informations about MPI parallelization
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! OUTPUT
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!   | keys=its access keys for compact storage.
!!
!! SIDE EFFECTS
!!  xred(3,natom)=reduced dimensionless atomic coordinates (in fact IN but here
!!                because of INOUT xredxcart() behavior).
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      derf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_wfs_set(fixmom, kpt, me, natom, nband, nkpt, nproc, nspinor, &
     & nsppol, nwfshist, occ, psps, rprimd, wfs, wtk, wvl, wvl_crmult, wvl_frmult, xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only: createWavefunctionsDescriptors, orbitals_descriptors, &
       & orbitals_communicators, allocate_diis_objects, wvl_timing => timing
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_wfs_set'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom, nkpt, nsppol, nspinor, nband, nwfshist,me,nproc
 real(dp), intent(in) :: fixmom, wvl_crmult, wvl_frmult
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_wf_type),intent(out) :: wfs
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 real(dp), intent(in) :: kpt(3,nkpt)
 real(dp), intent(in) :: wtk(nkpt), occ(:)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: idata, norb, norbu, norbd
 real(dp),parameter :: eps_mach=1.d-12
 logical :: parallel
 character(len=500) :: message
!arrays
 real(dp),allocatable :: xcart(:,:)

! *********************************************************************

 parallel = (nproc > 1)

#if defined HAVE_DFT_BIGDFT
!Consistency checks, are all pseudo true GTH pseudo with geometric informations?
 do idata = 1, psps%npsp, 1
   if (.not. psps%gth_params%set(idata)) then
     write(message, '(a,a,a,a,I0,a,a,a)' ) ch10,&
&     ' wvl_wfs_set:  consistency checks failed,', ch10, &
&     '  no GTH parameters found for type number ', idata, '.', ch10, &
&     '  Check your input pseudo files.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   if (.not. psps%gth_params%hasGeometry(idata)) then
     write(message, '(a,a,a,a,a,a)' ) ch10,&
&     ' wvl_wfs_set:  consistency checks failed,', ch10, &
&     '  the given GTH parameters has no geometry informations.', ch10, &
&     '  Upgrade your input pseudo files to GTH with geometric informatoins.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end do

 call wvl_timing(me,'CrtDescriptors','ON')
!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, natom))
 call xredxcart(natom, 1, rprimd, xcart, xred)

!Nullify possibly unset pointers
 nullify(wfs%psi)
 nullify(wfs%hpsi)
 nullify(wfs%psit)

!Static allocations.
 norb = nband
 norbu = 0
 norbd = 0
 if (nsppol == 2) then
   if (fixmom < -real(90, dp)) then
     norbu = min(norb / 2, norb)
   else
     norbu = min(norb / 2 + int(fixmom), norb)
   end if
   norbd = norb - norbu
 else
   norbu = norb
   norbd = 0
 end if

 call orbitals_descriptors(me, nproc,norb,norbu,norbd,nsppol,nspinor, &
      & nkpt,kpt,wtk,wfs%orbs)
 ! We copy occ_orig to wfs%orbs%occup
 wfs%orbs%occup(1:norb * nkpt) = occ(1:norb * nkpt)
 ! We allocate the eigen values storage.
 ABI_ALLOCATE(wfs%orbs%eval,(wfs%orbs%norb * wfs%orbs%nkpts))

 write(message, '(a,a)' ) ch10,&
& ' wvl_init_wfs_type: Create access keys for wavefunctions.'
 call wrtout(std_out,message,'COLL')

 call createWavefunctionsDescriptors(me, wvl%h(1), wvl%h(2), wvl%h(3), &
      & wvl%atoms, xcart, psps%gth_params%radii_cf, &
      & wvl_crmult, wvl_frmult, wfs%Glr)
!The memory is not allocated there for memory occupation optimisation reasons.

 call orbitals_communicators(me,nproc,wfs%Glr,wfs%orbs,wfs%comms)  

 write(message, '(a,2I8)' ) &
& '  | all orbitals have coarse segments, elements:', &
& wfs%Glr%wfd%nseg_c, wfs%Glr%wfd%nvctr_c
 call wrtout(std_out,message,'COLL')
 write(message, '(a,2I8)' ) &
& '  | all orbitals have fine   segments, elements:', &
& wfs%Glr%wfd%nseg_f, 7 * wfs%Glr%wfd%nvctr_f
 call wrtout(std_out,message,'COLL')

!Deallocations
 ABI_DEALLOCATE(xcart)

 call wvl_timing(me,'CrtDescriptors','OF')

 ! allocate arrays necessary for DIIS convergence acceleration
 call allocate_diis_objects(nwfshist,1._dp,&
      & sum(wfs%comms%ncntt(0:nproc-1)), wfs%orbs%nkptsp, wfs%orbs%nspinor, &
      & wfs%orbs%norbd, wfs%diis, "wvl_init_wfs_type")

#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_wfs_set : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif
end subroutine wvl_wfs_set
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/derfcf
!! NAME
!! derfcf
!!
!! FUNCTION
!! Some wrappers for BigDFT which uses different names for the same routines.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!      derf
!!
!! SOURCE
subroutine derfcf(derfc_yy,yy)

 use m_profiling
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'derfcf'
 use interfaces_32_util
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: yy
 real(dp),intent(out) :: derfc_yy

 call derfc(derfc_yy, yy)
end subroutine derfcf
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/derf_ab
!! NAME
!! derf_ab
!!
!! FUNCTION
!! Some wrappers for BigDFT which uses different names for the same routines.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!      derf
!!
!! SOURCE
subroutine derf_ab(derf_yy,yy)

 use m_profiling
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'derf_ab'
 use interfaces_32_util
!End of the abilint section

 implicit none
 real(dp),intent(in) :: yy
 real(dp),intent(out) :: derf_yy

 call derf(derf_yy, yy)
end subroutine derf_ab
!!***
