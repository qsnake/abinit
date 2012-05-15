#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


        subroutine multpot(icplexwf,icplex,includelast,nd1,nd2,n2,lot,n1dfft,pot,zw)

 use m_profiling
        use  defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'multpot'
!End of the abilint section

        implicit real(dp) (a-h,o-z)
 !Arguments ------------------------------------
        integer :: icplexwf,icplex,includelast,nd1,nd2,n2,lot,n1dfft
        real(dp) :: pot,zw
        dimension zw(2,lot,n2),pot(icplex*nd1,nd2)
!Local variables-------------------------------
! *************************************************************************
        if(icplexwf==1)then

         if(icplex==2)then

          write(std_out,*)' multpot : icplexwf=1 and icplex=2, bug '
          stop

         else

!         TO BE SPEEDED UP : should use the same trick as Stefan
          if(includelast==1)then
           do i2=1,n2
            do j=1,n1dfft
             zw(1,j,i2)=zw(1,j,i2)*pot(2*j-1,i2)
             zw(2,j,i2)=zw(2,j,i2)*pot(2*j  ,i2)
            end do
           end do
          else
           do i2=1,n2
            do j=1,n1dfft-1
             zw(1,j,i2)=zw(1,j,i2)*pot(2*j-1,i2)
             zw(2,j,i2)=zw(2,j,i2)*pot(2*j  ,i2)
            end do
             zw(1,n1dfft,i2)=zw(1,n1dfft,i2)*pot(2*n1dfft-1,i2)
           end do
          end if

         end if

        else if(icplexwf==2)then

         if(icplex==1)then

          do 10, i2=1,n2-1,2
          do 10, j=1,n1dfft
          zw(1,j,i2+0)=zw(1,j,i2+0)*pot(j,i2+0)
          zw(2,j,i2+0)=zw(2,j,i2+0)*pot(j,i2+0)
          zw(1,j,i2+1)=zw(1,j,i2+1)*pot(j,i2+1)
          zw(2,j,i2+1)=zw(2,j,i2+1)*pot(j,i2+1)
10        continue
          if(2*(n2/2)/=n2)then
          do 20, j=1,n1dfft
          zw(1,j,n2  )=zw(1,j,n2  )*pot(j,n2  )
          zw(2,j,n2  )=zw(2,j,n2  )*pot(j,n2  )
20        continue
          end if

         else

          do 30, i2=1,n2-1,2
          do 30, j=1,n1dfft
          zw(1,j,i2+0)=zw(1,j,i2+0)*pot(2*j-1,i2+0)&
&                     -zw(2,j,i2+0)*pot(2*j-0,i2+0)
          zw(2,j,i2+0)=zw(2,j,i2+0)*pot(2*j-1,i2+0)&
&                     +zw(1,j,i2+0)*pot(2*j-0,i2+0)
          zw(1,j,i2+1)=zw(1,j,i2+1)*pot(2*j-1,i2+1)&
&                     -zw(2,j,i2+1)*pot(2*j-0,i2+1)
          zw(2,j,i2+1)=zw(2,j,i2+1)*pot(2*j-1,i2+1)&
&                     +zw(1,j,i2+1)*pot(2*j-0,i2+1)
30        continue
          if(2*(n2/2)/=n2)then
          do 40, j=1,n1dfft
          zw(1,j,n2  )=zw(1,j,n2  )*pot(2*j-1,n2  )&
&                     -zw(2,j,n2  )*pot(2*j-0,n2  )
          zw(2,j,n2  )=zw(2,j,n2  )*pot(2*j-1,n2  )&
&                     +zw(1,j,n2  )*pot(2*j-0,n2  )
40        continue
          end if

         end if

        end if

        return

end subroutine multpot
