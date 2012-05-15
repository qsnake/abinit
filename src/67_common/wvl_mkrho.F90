!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_mkrho
!! NAME
!! wvl_mkrho
!!
!! FUNCTION
!! This method is just a wrapper around the BigDFT routine to compute the
!! density from the wavefunctions.
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
!!  mpi_enreg=informations about MPI parallelization
!!  occ(dtset%mband)=occupation numbers.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!!  rhor(dtset%nfft)=electron density in r space
!!
!! SIDE EFFECTS
!!  proj <type(wvl_projector_type)>=projectors informations for wavelets.
!!   | proj(OUT)=computed projectors.
!!
!! PARENTS
!!      afterscfloop,gstate,mkrho,scfcv_new,wvl_vtorho
!!
!! CHILDREN
!!      leave_new,sumrho,wrtout,xcomm_world
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wvl_mkrho(dtset, irrzon, mpi_enreg, phnons, rhor, wfs, wvl)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only : sumrho
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_mkrho'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(wvl_wf_type),intent(inout) :: wfs
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 real(dp),intent(inout) :: rhor(dtset%nfft,dtset%nspden)
 integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,  &
&               (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp), intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3))**(1-1/dtset%nsym),  &
&                                 (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))

!Local variables-------------------------------
!scalars
#if defined HAVE_DFT_BIGDFT
 integer :: comm,me,nproc
#else
 character(len=500) :: message
#endif

! *************************************************************************

#if defined HAVE_DFT_BIGDFT
 call xcomm_world(mpi_enreg,comm,myrank=me,mysize=nproc)

 call sumrho(me, nproc, wfs%orbs, wfs%Glr, &
& wvl%h(1) / 2.d0, wvl%h(2) / 2.d0, wvl%h(3) / 2.d0, &
& wfs%psi, rhor, mpi_enreg%nscatterarr, dtset%nsppol, wfs%GPU, wvl%atoms%symObj, &
& irrzon, phnons, wvl%rhodsc)
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_mkrho : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif
end subroutine wvl_mkrho
!!***
