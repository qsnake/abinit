!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_wfs_free
!!
!! NAME
!! wvl_wfs_free
!!
!! FUNCTION
!! Freeing routine.
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
!!  wfs <type(wvl_wf_type)>=wavefunctions informations in a wavelet basis.
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      deallocate_comms,deallocate_diis_objects,deallocate_lr,deallocate_orbs
!!      leave_new,memocc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_wfs_free(wfs)

 use defs_basis
 use defs_wvltypes
 use m_profiling
#if defined HAVE_DFT_BIGDFT
 use BigDFT_API
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_wfs_free'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

  !Arguments ------------------------------------
  !scalars
  type(wvl_wf_type),intent(inout) :: wfs
  !Local variables -------------------------
  integer :: i_all, i_stat
  character(len = *), parameter :: sub = "wvl_wfs_free"
#ifndef HAVE_DFT_BIGDFT
  character(len=500)       :: message
#endif
  ! *********************************************************************
#if defined HAVE_DFT_BIGDFT
  call deallocate_lr(wfs%Glr, sub)
  call deallocate_orbs(wfs%orbs, sub)
  call deallocate_diis_objects(wfs%diis, sub)
  call deallocate_comms(wfs%comms, sub)
  if (associated(wfs%orbs%eval))  then
    ABI_DEALLOCATE(wfs%orbs%eval)
  end if
#else
  write(message, '(a,a,a,a)' ) ch10,&
       & ' wvl_wfs_free: BigDFT library is not compiled.', ch10, &
       & '   Action, used the flag --enable-bigdft when configuring.'
  call wrtout(std_out,message,'COLL')
  call leave_new('COLL')
#endif
  if (associated(wfs%psi)) then
     i_all=-product(shape(wfs%psi))*kind(wfs%psi)
     ABI_DEALLOCATE(wfs%psi)
     i_stat = ABI_ALLOC_STAT
     call memocc(i_stat,i_all,'wfs%psi',sub)
  end if
  
  if (associated(wfs%hpsi)) then
     i_all=-product(shape(wfs%hpsi))*kind(wfs%hpsi)
     ABI_DEALLOCATE(wfs%hpsi)
     i_stat = ABI_ALLOC_STAT
     call memocc(i_stat,i_all,'wfs%hpsi',sub)
  end if

  if (associated(wfs%psit)) then
     i_all=-product(shape(wfs%psit))*kind(wfs%psit)
     ABI_DEALLOCATE(wfs%psit)
     i_stat = ABI_ALLOC_STAT
     call memocc(i_stat,i_all,'wfs%psit',sub)
  end if
end subroutine wvl_wfs_free
!!***
