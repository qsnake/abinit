#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


        subroutine switch(n1dfft,n2,lot,n1,lzt,zt,zw)

 use m_profiling
        use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'switch'
!End of the abilint section

        implicit real(dp) (a-h,o-z)
!Arguments ------------------------------------
        integer :: n1dfft,n2,lot,n1,lzt
        real(dp) :: zt,zw
        dimension zw(2,lot,n2),zt(2,lzt,n1)
!Local variables-------------------------------
! *************************************************************************
        do 200,j=1,n1dfft
        do 100,i=1,n2
        zw(1,j,i)=zt(1,i,j)
        zw(2,j,i)=zt(2,i,j)
100     continue
200     continue
        return

end subroutine switch

