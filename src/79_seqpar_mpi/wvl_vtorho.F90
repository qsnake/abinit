!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_vtorho
!! NAME
!! wvl_vtorho
!!
!! FUNCTION
!! Heart of the wavelet resolution, compute new wavefunctions mixed witf previous
!! by computing the gradient of the wavefunctions knowing the external potential.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=input variables.
!!  istep=id of the current iteration (first is 1).
!!  mpi_enreg=informations about MPI parallelization
!!  proj <type(wvl_projector_type)>=projectors informations for wavelets.
!!  vtrial(dtset%nfft)=external potential.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  energies <type(energies_type)>=storage for energies computed here :
!!   | e_kinetic(OUT)=kinetic energy part of total energy
!!   | e_localpsp(OUT)=local pseudopotential part of total energy
!!   | e_nonlocalpsp(OUT)=nonlocal pseudopotential part of total energy
!!  residm=max value for gradient in the minimisation process.
!!  rhor(dtset%nfft)=electron density in r space
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!  xred(3,natom)=reduced dimensionless atomic coordinates (in fact IN but here
!!                because of INOUT xredxcart() behavior).
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      calculate_energy_and_gradient,free_full_potential,full_local_potential
!!      hpsitopsi,leave_new,localhamiltonianapplication
!!      nonlocalhamiltonianapplication,psolver_kernel
!!      synchronizehamiltonianapplication,wrtout,wvl_mkrho,xcomm_world
!!      xredxcart
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wvl_vtorho(dtset, energies, irrzon, istep, mpi_enreg, &
     & phnons, residm, rhor, rprimd, vtrial, wvl, xred)

 use m_profiling

  use defs_basis
  use defs_datatypes
  use defs_abitypes
  use defs_wvltypes
  use m_energies, only : energies_type
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only : LocalHamiltonianApplication, hpsitopsi, locreg_descriptors, &
       & full_local_potential, free_full_potential, calculate_energy_and_gradient, &
       & NonLocalHamiltonianApplication, SynchronizeHamiltonianApplication
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_vtorho'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_62_poisson
 use interfaces_67_common
!End of the abilint section

  implicit none

!Arguments -------------------------------
  type(dataset_type), intent(in)         :: dtset
  type(energies_type), intent(inout)     :: energies
  integer, intent(in)                    :: istep
  type(MPI_type), intent(in)             :: mpi_enreg
  real(dp), intent(inout)                :: residm
  type(wvl_data), intent(inout)          :: wvl
  integer, intent(in)                    :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,&
       & (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp), intent(inout)                :: rhor(dtset%nfft)
  real(dp), intent(in)                   :: rprimd(3, 3)
  real(dp), intent(in)                   :: vtrial(dtset%nfft * dtset%nspden)
  real(dp), intent(inout)                :: xred(3, dtset%natom)
  real(dp), intent(in)                   :: phnons(2,dtset%nfft**(1-1/dtset%nsym),&
       & (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))

!Local variables-------------------------------
  real(dp), save        :: alpha, etotal_local, etotal_min
  integer, save         :: ids
  character(len = 500)  :: message
  real(dp)              :: epot_sum, ekin_sum, eproj_sum
  real(dp)              :: gnrm_zero
  integer               :: comm,me,nproc,nrhodim,i3rho_add
  real(dp), allocatable :: xcart(:, :)
  real(dp), pointer     :: kernelseq(:), potential(:)
  character(len = 1)    :: bndcode
! *********************************************************************

 if (dtset%icoulomb == 0) then
   bndcode = 'P'
 else if (dtset%icoulomb == 1) then
   bndcode = 'F'
 else if (dtset%icoulomb == 2) then
   bndcode = 'S'
 end if

#if defined HAVE_DFT_BIGDFT

 write(message, '(a,a)' ) ch10,&
& ' wvl_vtorho: compute the new density from the trial potential.'
 call wrtout(std_out,message,'COLL')

 call xcomm_world(mpi_enreg,comm,myrank=me,mysize=nproc)

!Initialisation of mixing parameter alpha
 if (istep == 1) then
   alpha        = real(1., dp)
   etotal_min   = real(1.d100, dp)
   etotal_local = real(1.d100, dp)
   ids          = dtset%nwfshist
 end if 

!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 call xredxcart(dtset%natom, 1, rprimd, xcart, xred)

!Define sizes for rho / pot.
 nrhodim   = dtset%nsppol
 i3rho_add = 0
 if (wvl%descr%SIC%approach == 'NK') then
   nrhodim   = 2 * nrhodim
   i3rho_add = wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i * mpi_enreg%nscatterarr(me, 4) + 1
 end if

!allocate the potential in the full box
 call full_local_potential(me, nproc, &
& wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i * mpi_enreg%nscatterarr(me, 2), &
& wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i * wvl%descr%Glr%d%n3i, dtset%nsppol, &
& wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i * mpi_enreg%nscatterarr(me, 1) * nrhodim, &
& i3rho_add, wvl%wfs%orbs%norb, wvl%wfs%orbs%norbp, mpi_enreg%ngatherarr, &
& vtrial, potential)

!Apply vtrial and the projectors to the wavefubctions, computing HPsi.
 call PSolver_kernel(dtset, 4, kernelseq, mpi_enreg, rprimd, wvl%descr)

 call LocalHamiltonianApplication(me, nproc, wvl%descr%atoms, wvl%wfs%orbs, &
& wvl%descr%h(1), wvl%descr%h(2), wvl%descr%h(3), xcart, &
& wvl%wfs%Glr, mpi_enreg%ngatherarr, potential, wvl%wfs%psi, wvl%wfs%hpsi, &
& ekin_sum, epot_sum, energies%e_exactX, energies%e_sicdc, wvl%descr%SIC, &
& wvl%wfs%GPU, pkernel = kernelseq)

 call NonLocalHamiltonianApplication(me, nproc, wvl%descr%atoms, wvl%wfs%orbs, &
& wvl%descr%h(1), wvl%descr%h(2), wvl%descr%h(3), xcart, &
& wvl%projectors%keys, wvl%projectors%proj, wvl%wfs%Glr, wvl%wfs%psi, &
& wvl%wfs%hpsi, eproj_sum)

 call SynchronizeHamiltonianApplication(nproc, wvl%wfs%orbs, wvl%wfs%Glr, wvl%wfs%GPU, &
& wvl%wfs%hpsi, ekin_sum, epot_sum, eproj_sum, energies%e_sicdc, energies%e_exactX)

!deallocate potential
 call free_full_potential(nproc, potential, "wvl_vtorho")

 ABI_DEALLOCATE(xcart)

!WARNING! e_hartree is taken from the previous iteration as e_xc
!Update physical values
 energies%e_kinetic = ekin_sum
 energies%e_localpsp = epot_sum - real(2., dp) * energies%e_hartree
 energies%e_nonlocalpsp = eproj_sum
 energies%e_corepsp = real(0., dp)

!Precondition, minimise (DIIS or steepest descent) and ortho.
!Compute also the norm of the gradient.
 call calculate_energy_and_gradient(istep, me, nproc, wvl%wfs%orbs, wvl%wfs%comms, &
& wvl%wfs%GPU, wvl%wfs%Glr, wvl%descr%h(1), wvl%descr%h(2), wvl%descr%h(3), &
& dtset%wvl_nprccg, 0, ekin_sum, epot_sum, eproj_sum, energies%e_sicdc, &
& energies%e_hartree, energies%e_xc, energies%e_vxc, energies%e_exactX, &
& energies%e_ewald, 0._dp, wvl%wfs%psi, wvl%wfs%psit, wvl%wfs%hpsi, &
& residm, gnrm_zero, wvl%wfs%diis%energy)
 etotal_local = wvl%wfs%diis%energy

 call hpsitopsi(me, nproc, wvl%wfs%orbs, wvl%wfs%Glr, wvl%wfs%comms, istep, &
& wvl%wfs%diis, ids, wvl%wfs%psi, wvl%wfs%psit, wvl%wfs%hpsi, dtset%nsppol, &
& wvl%descr%orthpar)


!Density from wavefunctions.
 call wvl_mkrho(dtset, irrzon, mpi_enreg, phnons, rhor, wvl%wfs, wvl%descr)
 
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_vtorho : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif
end subroutine wvl_vtorho
!!***
