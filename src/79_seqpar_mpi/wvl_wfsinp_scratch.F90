!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_wfsinp_scratch
!! NAME
!! wvl_wfsinp_scratch
!!
!! FUNCTION
!! This method allocates and initialises wavefunctions with values from input guess.
!! See wvl_wfsinp_disk() or wvl_wfsinp_reformat() from other initialisation
!! routines.
!!
!! When initialised from scratch or from disk, wvl%wfs%[h]psi comes unallocated
!! and will be allocated inside this routine.
!! When initialised from memory (reformating), wvl%wfs%[h]psi will be reallocated.
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
!!  hdr0 <type(hdr_type)>=the header of wf, den and pot files (read from restart)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  ireadwf=1 for reading from file, 0 otherwise.
!!  mpi_enreg=informations about MPI parallelization
!!  option=1 for reading a file following ABINIT format, -1 for a BigDFT format.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wff <type(wffile_type)>= structure with informations on wf file.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  wvl <type(wvl_data)>=wavefunctions & projectors informations for wavelets.
!!  xred(3,natom)=reduced dimensionless atomic coordinates (in fact IN but here
!!                because of INOUT xredxcart() behavior).
!!
!! PARENTS
!!      inwffil
!!
!! CHILDREN
!!      calculate_rhocore,createionicpotential,input_wf_diag
!!      kpoints_get_irreductible_zone,leave_new,memocc,psolver_kernel
!!      symmetry_get_n_sym,wrtout,xcomm_world,xredxcart
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wvl_wfsinp_scratch(dtset, mpi_enreg, rprimd, wvl, xred)

  use defs_basis
  use defs_datatypes
  use m_wffile
  use m_profiling
  use m_ab6_kpoints
  use m_ab6_symmetry
  use defs_abitypes
  use defs_wvltypes
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only : createIonicPotential, input_wf_diag, gaussian_basis, input_variables, calculate_rhocore
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_wfsinp_scratch'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_62_poisson
!End of the abilint section

  implicit none

!Arguments -------------------------------
  !scalars
  type(dataset_type), intent(in)        :: dtset
  type(MPI_type), intent(in)            :: mpi_enreg
  type(wvl_data), intent(inout)         :: wvl
  !arrays
  real(dp), intent(in)                  :: rprimd(3, 3)
  real(dp), intent(inout)               :: xred(3, dtset%natom)

!Local variables-------------------------------
  character(len = 500)  :: message
  integer               :: comm,me,nproc, i_stat, nsym, i_all
  integer               :: nvirt, nspin
  real(dp), allocatable :: xcart(:,:)
  real(dp), allocatable :: rhor(:,:), vpsp(:)
  real(dp), pointer     :: kernel(:), kernelseq(:), rhocore(:)
  integer, allocatable  :: irrzon(:,:,:)
  real(dp), allocatable :: phnons(:,:,:)
  character(len = *), parameter :: subname = "wvl_wfsinp_scratch"

#if defined HAVE_DFT_BIGDFT
  type(gaussian_basis) :: Gvirt
  
  ! To be removed, waiting for BigDFT upgrade.
  type(input_variables) :: in
#endif

! *********************************************************************

 write(message, '(a,a)' ) ch10,&
& ' wvl_wfsinp_scratch: wavefunction initialisation.'
 call wrtout(std_out,message,'COLL')

#if defined HAVE_DFT_BIGDFT
 call xcomm_world(mpi_enreg,comm,myrank=me,mysize=nproc)
!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 call xredxcart(dtset%natom, 1, rprimd, xcart, xred)

!We grep the kernel for the Poisson solver.
 call PSolver_kernel(dtset, 2, kernel, mpi_enreg, rprimd, wvl%descr)
!We get a sequential kernel if needed.
 if ((wvl%descr%exctXfac /= 0.0_dp .and. wvl%descr%exctxpar=='OP2P' .or. &
& wvl%descr%SIC%alpha /= 0.0_dp) .and. nproc > 1) then
   call PSolver_kernel(dtset, 4, kernelseq, mpi_enreg, rprimd, wvl%descr)
 else 
   kernelseq => kernel
 end if

!We allocate temporary arrays for rho and vpsp.
!allocate ionic potential
 if (mpi_enreg%ngfft3_ionic > 0) then
   ABI_ALLOCATE(vpsp,(wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i * mpi_enreg%ngfft3_ionic))
 else
   ABI_ALLOCATE(vpsp,(1))
 end if

 call createIonicPotential(wvl%descr%atoms%geocode, &
& me, nproc, wvl%descr%atoms, xcart, wvl%descr%h(1) / 2., &
& wvl%descr%h(2) / 2., wvl%descr%h(3) / 2., &
& real(0, dp), wvl%descr%Glr%d%n1, wvl%descr%Glr%d%n2, wvl%descr%Glr%d%n3, &
& mpi_enreg%ngfft3_ionic, mpi_enreg%nscatterarr(me, 3) + 1, &
& wvl%descr%Glr%d%n1i, wvl%descr%Glr%d%n2i, wvl%descr%Glr%d%n3i, &
& kernel, vpsp, 0.d0, 0, .false.)

 nspin = dtset%nsppol
!Allocate Charge density, Potential in real space
 if (mpi_enreg%nscatterarr(me, 1) > 0) then
   ABI_ALLOCATE(rhor,(wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i * mpi_enreg%nscatterarr(me, 1), nspin))
 else
   ABI_ALLOCATE(rhor,(1, nspin))
 end if

!calculate the irreductible zone, if necessary.
 if (wvl%descr%atoms%symObj >= 0) then
   call symmetry_get_n_sym(wvl%descr%atoms%symObj, nsym, i_stat)
   if (nsym > 1) then
!    Current third dimension is set to 1 always
!    since nspin == nsppol always in BigDFT
     ABI_ALLOCATE(irrzon,(wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i * wvl%descr%Glr%d%n3i,2,1+ndebug))
     i_stat = ABI_ALLOC_STAT
     call memocc(i_stat,irrzon,'irrzon',subname)
     ABI_ALLOCATE(phnons,(2,wvl%descr%Glr%d%n1i * wvl%descr%Glr%d%n2i * wvl%descr%Glr%d%n3i,1+ndebug))
     i_stat = ABI_ALLOC_STAT
     call memocc(i_stat,phnons,'phnons',subname)
     call kpoints_get_irreductible_zone(irrzon, phnons, &
&     wvl%descr%Glr%d%n1i, wvl%descr%Glr%d%n2i, wvl%descr%Glr%d%n3i, &
&     dtset%nsppol, dtset%nsppol, wvl%descr%atoms%symObj, i_stat)
   end if
 end if
 if (.not. allocated(irrzon)) then
!  Allocate anyway to small size other size the bounds check does not pass.
   ABI_ALLOCATE(irrzon,(1,2,1+ndebug))
   i_stat = ABI_ALLOC_STAT
   call memocc(i_stat,irrzon,'irrzon',subname)
   ABI_ALLOCATE(phnons,(2,1,1+ndebug))
   i_stat = ABI_ALLOC_STAT
   call memocc(i_stat,phnons,'phnons',subname)
 end if

!check if non-linear core correction should be applied, and allocate the 
!pointer if it is the case
 call calculate_rhocore(me,wvl%descr%atoms,wvl%descr%Glr%d,xcart,wvl%descr%h(1) / 2., &
& wvl%descr%h(2) / 2., wvl%descr%h(3) / 2., &
& mpi_enreg%nscatterarr(me, 3) - mpi_enreg%nscatterarr(me, 4) + 1, &
& mpi_enreg%nscatterarr(me, 4), mpi_enreg%nscatterarr(me, 1), &
& mpi_enreg%nscatterarr(me, 2),rhocore)

!This routine allocates psi, hpsi and psit.
 nvirt = 0
 in%nspin = dtset%nsppol
 in%exctxpar = wvl%descr%exctxpar
 in%itrpmax = 1
 in%Tel = dtset%tphysel
 in%SIC%approach = wvl%descr%SIC%approach
 in%SIC%ixc = wvl%descr%SIC%ixc
 in%SIC%alpha = wvl%descr%SIC%alpha
 in%SIC%fref = wvl%descr%SIC%fref
 in%orthpar%directDiag = wvl%descr%orthpar%directDiag
 in%orthpar%norbpInguess = wvl%descr%orthpar%norbpInguess
 in%orthpar%bsLow = wvl%descr%orthpar%bsLow
 in%orthpar%bsUp = wvl%descr%orthpar%bsUp
 in%orthpar%methOrtho = wvl%descr%orthpar%methOrtho
 in%orthpar%iguessTol = wvl%descr%orthpar%iguessTol
 call input_wf_diag(me, nproc, &
& wvl%descr%atoms, wvl%descr%rhodsc, wvl%wfs%orbs, nvirt, wvl%wfs%comms, &
& wvl%wfs%Glr, wvl%descr%h(1), wvl%descr%h(2), wvl%descr%h(3), xcart, &
& rhor, rhocore, vpsp, wvl%projectors%keys, wvl%projectors%proj, kernel, kernelseq, &
& dtset%ixc, wvl%wfs%psi, wvl%wfs%hpsi, wvl%wfs%psit, Gvirt, mpi_enreg%nscatterarr, &
& mpi_enreg%ngatherarr, nspin, 0, wvl%descr%atoms%symObj, irrzon, phnons, wvl%wfs%GPU, &
& in)

 i_all=-product(shape(irrzon))*kind(irrzon)
 ABI_DEALLOCATE(irrzon)
 i_stat = ABI_ALLOC_STAT
 call memocc(i_stat,i_all,'irrzon',subname)

 i_all=-product(shape(phnons))*kind(phnons)
 ABI_DEALLOCATE(phnons)
 i_stat = ABI_ALLOC_STAT
 call memocc(i_stat,i_all,'phnons',subname)

 write(message, '(a)' ) &
& '  | wavefunctions have been calculated.'
 call wrtout(std_out,message,'COLL')

 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(vpsp)
 ABI_DEALLOCATE(rhor)

#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_wfsinp_scratch: BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif
end subroutine wvl_wfsinp_scratch
!!***

