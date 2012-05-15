!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_tail_corrections
!! NAME
!! wvl_tail_corrections
!!
!! FUNCTION
!! Perform a minimization on the wavefunctions (especially the treatment
!! of the kinetic operator) with exponentialy decreasing functions on
!! boundaries.
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
!!      afterscfloop
!!
!! CHILDREN
!!      calculatetailcorrection,leave_new,wrtout,xallgatherv_mpi,xcomm_world
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_tail_corrections(dtset, energies, etotal, mpi_enreg, psps, &
     & vtrial, wvl, xcart)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_energies, only : energies_type
 use m_xmpi
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only: CalculateTailCorrection
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_tail_corrections'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: etotal
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(energies_type),intent(inout) :: energies
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_data),intent(inout) :: wvl
!arrays
 real(dp),intent(in) :: xcart(3,dtset%natom)
 real(dp),intent(in),target :: vtrial(dtset%nfft)

!Local variables-------------------------------
!scalars
 integer :: ierr,me,nbuf,nproc,nsize,spaceComm,vtrial_shift
 real(dp) :: ekin_sum,epot_sum,eproj_sum
 logical :: parallel
 character(len=500) :: message
!arrays
 integer :: ntails(3)
 real(dp) :: atails(3)
 real(dp),pointer :: vtotal(:)

! *************************************************************************
 
 call xcomm_world(mpi_enreg,spaceComm,myrank=me,mysize=nproc)
 parallel = (nproc > 1)

!Write a message with the total energy before tail corrections.
 etotal = energies%e_kinetic + energies%e_localpsp + energies%e_nonlocalpsp + &
& energies%e_hartree + energies%e_xc - energies%e_vxc + &
& energies%e_ewald + energies%e_corepsp
 write(message,'(a,2x,e19.12)') ' Total energy before tail correction', etotal
 call wrtout(std_out, message, 'COLL')

#if defined HAVE_DFT_BIGDFT
!Calculate kinetic energy correction due to boundary conditions
 nbuf = nint(dtset%tl_radius / dtset%wvl_hgrid)
 ntails = (/ wvl%descr%Glr%d%n1, wvl%descr%Glr%d%n2, wvl%descr%Glr%d%n3 /) + 2 * nbuf
 atails = real(ntails, dp) * dtset%wvl_hgrid
 write(message,'(a,a,i6,a,A,A,3F12.6,A,A,3I12,A)') ch10,&
& ' Tail requires ',nbuf,' additional grid points around cell.', ch10, &
& '  | new acell:', atails, ch10, &
& '  | new box size for wavelets:', ntails, ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

!---reformat potential
 if (parallel) then
   ABI_ALLOCATE(vtotal,(wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i * wvl%descr%Glr%d%n3i))
   nsize = wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i
   vtrial_shift = 1 + nsize * mpi_enreg%nscatterarr(me, 4)
   call xallgatherv_mpi(vtrial(vtrial_shift:dtset%nfft), &
&   nsize * mpi_enreg%nscatterarr(me, 2), &
&   vtotal,  nsize * mpi_enreg%nscatterarr(:,2), &
&   nsize * mpi_enreg%nscatterarr(:,3), spaceComm, ierr)
 else
   vtotal => vtrial
 end if

 call CalculateTailCorrection(me, nproc, wvl%descr%atoms, dtset%tl_radius, &
& wvl%wfs%orbs, wvl%wfs%Glr, wvl%projectors%keys, dtset%tl_nprccg, &
& vtotal, dtset%wvl_hgrid, xcart, psps%gth_params%radii_cf, &
& dtset%wvl_crmult, dtset%wvl_frmult, dtset%nsppol, &
& wvl%projectors%proj, wvl%wfs%psi, .false., ekin_sum, epot_sum, eproj_sum)
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_tail_corrections: BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif

 if (parallel)  then
   ABI_DEALLOCATE(vtotal)
 end if

 energies%e_kinetic = ekin_sum
 energies%e_localpsp = epot_sum - real(2., dp) * energies%e_hartree
 energies%e_nonlocalpsp = eproj_sum
 energies%e_corepsp = real(0., dp)
 etotal = energies%e_kinetic + energies%e_localpsp + energies%e_nonlocalpsp + &
& energies%e_hartree + energies%e_xc - energies%e_vxc + &
& energies%e_ewald + energies%e_corepsp

 write(message,'(a,3(1x,e18.11))') ' ekin_sum,epot_sum,eproj_sum',  & 
 ekin_sum,epot_sum,eproj_sum
 call wrtout(std_out, message, 'COLL')
 write(message,'(a,3(1x,e18.11))') ' ehart,eexcu,vexcu', &
& energies%e_hartree,energies%e_xc,energies%e_vxc
 call wrtout(std_out, message, 'COLL')
 write(message,'(a,2x,e19.12)') ' Total energy with tail correction', etotal
 call wrtout(std_out, message, 'COLL')

!--- End if of tail calculation
end subroutine wvl_tail_corrections
!!***
