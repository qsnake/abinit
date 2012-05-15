#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


        subroutine fill_cent(md1,md3,lot,n1dfft,max3,m3,n3,zf,zw)

 use m_profiling
        use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fill_cent'
!End of the abilint section

        implicit real(dp) (a-h,o-z)
!Arguments ------------------------------------
        integer:: md1,md3,lot,n1dfft,max3,m3,n3
        real(dp):: zw,zf
        dimension zw(2,lot,n3),zf(2,md1,md3)
!Local variables-------------------------------

! *************************************************************************
!       Here, zero and positive frequencies
        do 90,i3=1,max3+1
        do i1=1,n1dfft
        zw(1,i1,i3)=zf(1,i1,i3)
        zw(2,i1,i3)=zf(2,i1,i3)
        end do
90      continue

!       Fill the center region with zeros
        do 100,i3=max3+2,n3-m3+max3+1
        do i1=1,n1dfft
        zw(1,i1,i3)=zero
        zw(2,i1,i3)=zero
        end do
100     continue

!       Here, negative frequencies
        do 110,i3=max3+2,m3
        do i1=1,n1dfft
        zw(1,i1,i3+n3-m3)=zf(1,i1,i3)
        zw(2,i1,i3+n3-m3)=zf(2,i1,i3)
        end do
110     continue

        return
end subroutine fill_cent
