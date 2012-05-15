#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


        subroutine unswitchreal_cent(n1dfft,max2,n2,lot,n1zt,lzt,zw,zt)

 use m_profiling
        use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unswitchreal_cent'
!End of the abilint section

        implicit real(dp) (a-h,o-z)
!Arguments ------------------------------------
        integer :: n1dfft,max2,n2,lot,n1zt,lzt
        real(dp) :: zw,zt
        dimension zw(2,lot,n2),zt(2,lzt,n1zt)
!Local variables-------------------------------
! *************************************************************************
        do j=1,n1dfft
         zt(1,1,2*j-1)=zw(1,j,1)
         zt(2,1,2*j-1)=zero
         zt(1,1,2*j)  =zw(2,j,1)
         zt(2,1,2*j)  =zero
        end do

        do i=2,max2+1
         do j=1,n1dfft
          zt(1,i,2*j-1)= (zw(1,j,i)+zw(1,j,n2+2-i))*half
          zt(2,i,2*j-1)= (zw(2,j,i)-zw(2,j,n2+2-i))*half
          zt(1,i,2*j)  = (zw(2,j,i)+zw(2,j,n2+2-i))*half
          zt(2,i,2*j)  =-(zw(1,j,i)-zw(1,j,n2+2-i))*half
         end do
        end do



!       Here, zero and positive frequencies
!        do 90,j=1,n1dfft
!        do 90,i=1,max2+1
!        zt(1,i,j)=zw(1,j,i)
!        zt(2,i,j)=zw(2,j,i)
!90      continue

!       Here, negative frequencies
!        if(m2>=max2+2)then
!         do 110,j=1,n1dfft
!         do 110,i=max2+2,m2
!         zt(1,i,j)=zw(1,j,i+n2-m2)
!         zt(2,i,j)=zw(2,j,i+n2-m2)
!110      continue
!        end if

        return

end subroutine unswitchreal_cent
