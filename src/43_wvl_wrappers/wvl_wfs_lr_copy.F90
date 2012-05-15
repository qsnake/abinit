!!****f* defs_wvltypes/wvl_wfs_lr_copy
!!
!! NAME
!! wvl_wfs_lr_copy
!!
!! FUNCTION
!! Copy the wvl%Glr datastructure geometry part to wfs%Glr.
!!
!! INPUTS
!! wvl <type(wvl_internal_type)> = input localisation region
!!
!! OUTPUT
!! wfs <type(wvl_wf_type)> = output localistaion region
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wvl_wfs_lr_copy(wfs, wvl)

 use m_profiling

  use defs_wvltypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_wfs_lr_copy'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  type(wvl_internal_type), intent(in)  :: wvl
  type(wvl_wf_type), intent(inout)         :: wfs
  !arrays

!Local variables-------------------------------

#if defined HAVE_DFT_BIGDFT
    ! Use global localization region for the moment.
    wfs%Glr%geocode    = wvl%Glr%geocode
    wfs%Glr%hybrid_on  = wvl%Glr%hybrid_on
    wfs%Glr%ns1        = wvl%Glr%ns1
    wfs%Glr%ns2        = wvl%Glr%ns2
    wfs%Glr%ns3        = wvl%Glr%ns3
    wfs%Glr%d%n1       = wvl%Glr%d%n1
    wfs%Glr%d%n2       = wvl%Glr%d%n2
    wfs%Glr%d%n3       = wvl%Glr%d%n3
    wfs%Glr%d%nfl1     = wvl%Glr%d%nfl1
    wfs%Glr%d%nfu1     = wvl%Glr%d%nfu1
    wfs%Glr%d%nfl2     = wvl%Glr%d%nfl2
    wfs%Glr%d%nfu2     = wvl%Glr%d%nfu2
    wfs%Glr%d%nfl3     = wvl%Glr%d%nfl3
    wfs%Glr%d%nfu3     = wvl%Glr%d%nfu3
    wfs%Glr%d%n1i      = wvl%Glr%d%n1i
    wfs%Glr%d%n2i      = wvl%Glr%d%n2i
    wfs%Glr%d%n3i      = wvl%Glr%d%n3i
#endif
  
end subroutine wvl_wfs_lr_copy
!!***
