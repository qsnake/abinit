
!!****f* ABINIT/wvl_projectors_free
!!
!! NAME
!! wvl_projectors_free
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
!!  proj <type(wvl_projectors_type)>=projectors informations in a wavelet basis.
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      leave_new,memocc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_projectors_free(proj)

 use defs_basis
 use m_profiling
 use defs_wvltypes
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_projectors_free'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

  !Arguments ------------------------------------
  !scalars
  type(wvl_projectors_type),intent(inout) :: proj
  !Local variables -------------------------
  integer :: i_all, i_stat
#ifndef HAVE_DFT_BIGDFT
  character(len=500)       :: message
#endif
  ! *********************************************************************
#if defined HAVE_DFT_BIGDFT
 i_all=-product(shape(proj%keys%nvctr_p))*kind(proj%keys%nvctr_p)
 ABI_DEALLOCATE(proj%keys%nvctr_p)
 i_stat = ABI_ALLOC_STAT
 call memocc(i_stat,i_all,'proj%keys%nvctr_p',"wvl_projectors_free")

 i_all=-product(shape(proj%keys%nseg_p))*kind(proj%keys%nseg_p)
 ABI_DEALLOCATE(proj%keys%nseg_p)
 i_stat = ABI_ALLOC_STAT
 call memocc(i_stat,i_all,'proj%keys%nseg_p',"wvl_projectors_free")

 i_all=-product(shape(proj%keys%keyv_p))*kind(proj%keys%keyv_p)
 ABI_DEALLOCATE(proj%keys%keyv_p)
 i_stat = ABI_ALLOC_STAT
 call memocc(i_stat,i_all,'proj%keys%keyv_p',"wvl_projectors_free")

 i_all=-product(shape(proj%keys%keyg_p))*kind(proj%keys%keyg_p)
 ABI_DEALLOCATE(proj%keys%keyg_p)
 i_stat = ABI_ALLOC_STAT
 call memocc(i_stat,i_all,'proj%keys%keyg_p',"wvl_projectors_free")

 i_all=-product(shape(proj%keys%nboxp_c))*kind(proj%keys%nboxp_c)
 ABI_DEALLOCATE(proj%keys%nboxp_c)
 i_stat = ABI_ALLOC_STAT
 call memocc(i_stat,i_all,'proj%keys%nboxp_c',"wvl_projectors_free")

 i_all=-product(shape(proj%keys%nboxp_f))*kind(proj%keys%nboxp_f)
 ABI_DEALLOCATE(proj%keys%nboxp_f)
 i_stat = ABI_ALLOC_STAT
 call memocc(i_stat,i_all,'proj%keys%nboxp_f',"wvl_projectors_free")
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_projectors_free: BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif
 i_all=-product(shape(proj%proj))*kind(proj%proj)
 ABI_DEALLOCATE(proj%proj)
 i_stat = ABI_ALLOC_STAT
 call memocc(i_stat,i_all,'proj%proj',"wvl_projectors_free")
end subroutine wvl_projectors_free
!!***
