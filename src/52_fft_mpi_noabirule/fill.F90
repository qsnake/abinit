#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


        subroutine fill(nd1,nd3,lot,n1dfft,n3,zf,zw)

 use m_profiling
        use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fill'
!End of the abilint section

        implicit real(dp) (a-h,o-z)
!Arguments ------------------------------------
        integer :: nd1,nd3,lot,n1dfft,n3
        real(dp) ::zw,zf
        dimension zw(2,lot,n3),zf(2,nd1,nd3)
!Local variables-------------------------------

! *************************************************************************
        do 100,i3=1,n3
        do 100,i1=1,n1dfft
        zw(1,i1,i3)=zf(1,i1,i3)
        zw(2,i1,i3)=zf(2,i1,i3)
100     continue

        return

end subroutine fill


