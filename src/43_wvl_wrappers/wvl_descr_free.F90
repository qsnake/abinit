!!****f* defs_wvltypes/wvl_descr_free
!!
!! NAME
!! wvl_descr_free
!!
!! FUNCTION
!! Free the wvl%atoms% datastructure (deallocate or nullify)
!!
!! INPUTS
!! wvl <type(wvl_internal_type)>=internal variables for wavelets
!!
!! OUTPUT
!! wvl <type(wvl_internal_type)>=internal variables for wavelets
!!
!! PARENTS
!!      gstate,wvl_memory
!!
!! CHILDREN
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_descr_free(wvl)

 use m_profiling
  
  use defs_wvltypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_descr_free'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  type(wvl_internal_type), intent(inout) :: wvl
!arrays

!Local variables-------------------------------
!scalars

! *********************************************************************

#if defined HAVE_DFT_BIGDFT

 if (associated(wvl%atoms%atomnames))  then
   ABI_DEALLOCATE(wvl%atoms%atomnames)
 end if
 nullify(wvl%atoms%psppar)
 nullify(wvl%atoms%npspcode)
 nullify(wvl%atoms%ixcpsp)
 if (associated(wvl%atoms%iatype))    then
   ABI_DEALLOCATE(wvl%atoms%iatype)
 end if
 if (associated(wvl%atoms%iasctype))  then
   ABI_DEALLOCATE(wvl%atoms%iasctype)
 end if
 if (associated(wvl%atoms%nzatom))    then
   ABI_DEALLOCATE(wvl%atoms%nzatom)
 end if
 if (associated(wvl%atoms%nelpsp))    then
   ABI_DEALLOCATE(wvl%atoms%nelpsp)
 end if
 if (associated(wvl%atoms%ifrztyp))   then
   ABI_DEALLOCATE(wvl%atoms%ifrztyp)
 end if
 if (associated(wvl%atoms%natpol))    then
   ABI_DEALLOCATE(wvl%atoms%natpol)
 end if
 if (associated(wvl%atoms%aocc))      then
   ABI_DEALLOCATE(wvl%atoms%aocc)
 end if
 if (associated(wvl%atoms%amu))       then
   ABI_DEALLOCATE(wvl%atoms%amu)
 end if
 if (associated(wvl%atoms%nlcc_ngc))  then
   ABI_DEALLOCATE(wvl%atoms%nlcc_ngc)
 end if
 if (associated(wvl%atoms%nlcc_ngv))  then
   ABI_DEALLOCATE(wvl%atoms%nlcc_ngv)
 end if
 if (associated(wvl%atoms%nlccpar))   then
   ABI_DEALLOCATE(wvl%atoms%nlccpar)
 end if
#endif
end subroutine wvl_descr_free
!!***
