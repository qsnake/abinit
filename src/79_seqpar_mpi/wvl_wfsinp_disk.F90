!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_wfsinp_disk
!! NAME
!! wvl_wfsinp_disk
!!
!! FUNCTION
!! This method allocates and initialises wavefunctions with values from disk.
!! See wvl_wfsinp_scratch() or wvl_wfsinp_reformat() from other initialisation
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
!!  mpi_enreg=informations about MPI parallelization
!!  option=1 for reading a file following ABINIT format, -1 for a BigDFT format.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wff <type(wffile_type)>= structure with informations on wf file.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!  xred(3,natom)=reduced dimensionless atomic coordinates (in fact IN but here
!!                because of INOUT xredxcart() behavior).
!!
!! PARENTS
!!      inwffil
!!
!! CHILDREN
!!      first_orthon,leave_new,wrtout,wvl_read,xcomm_world
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wvl_wfsinp_disk(dtset, hdr0, hdr, mpi_enreg, option, &
     & rprimd, wff, wfs, wvl, xred)

 use m_profiling

  use defs_basis
  use defs_datatypes
  use m_wffile

 use defs_abitypes
  use defs_wvltypes
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only : first_orthon
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_wfsinp_disk'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
 use interfaces_62_wvl_wfs
!End of the abilint section

  implicit none

!Arguments -------------------------------
  !scalars
  integer, intent(in)                       :: option
  type(dataset_type), intent(in)            :: dtset
  type(hdr_type), intent(in)                :: hdr0
  type(hdr_type), intent(in)                :: hdr
  type(MPI_type), intent(in)                :: mpi_enreg
  type(wffile_type), intent(in)             :: wff
  type(wvl_wf_type), intent(inout)          :: wfs
  type(wvl_internal_type), intent(in)       :: wvl
  !arrays
  real(dp), intent(in)                      :: rprimd(3, 3)
  real(dp), intent(inout)                   :: xred(3, dtset%natom)

!Local variables-------------------------------
#if defined HAVE_DFT_BIGDFT
  integer :: comm,me,nproc
#endif
  character(len = 500)  :: message

! *********************************************************************

 write(message, '(a,a)' ) ch10,&
& ' wvl_wfsinp_disk: wavefunction initialisation.'
 call wrtout(std_out,message,'COLL')

#if defined HAVE_DFT_BIGDFT

 call xcomm_world(mpi_enreg,comm,myrank=me,mysize=nproc)

!We allocate psi.
 ABI_ALLOCATE(wfs%psi,(wfs%orbs%npsidim))
 write(message, '(a,a,a,a,I0)' ) ch10, &
& ' wvl_wfsinp_disk: allocate wavefunctions,', ch10, &
& '  size of the compressed array per proc: ', &
& product(shape(wfs%psi))
 call wrtout(std_out,message,'COLL')

 call wvl_read(dtset, hdr0, hdr, mpi_enreg, option, rprimd, wff, wfs, wvl, xred)

!We orthogonalise.
 call first_orthon(me, nproc, wfs%orbs, wfs%Glr%wfd, wfs%comms, &
& wfs%psi, wfs%hpsi, wfs%psit, wvl%orthpar)
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_wfs_inp: BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif
end subroutine wvl_wfsinp_disk
!!***
