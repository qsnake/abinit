#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


        subroutine unmpiswitch(j3,n1dfft,Jp2st,J2st,lot,n1,nd2proc,nd3proc,nproc,ioption,zw,zmpi1)

 use m_profiling
        use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unmpiswitch'
!End of the abilint section

        implicit real(dp) (a-h,o-z)

 !Arguments ------------------------------------
        integer :: j3,n1dfft,Jp2st,J2st,lot,n1,nd2proc,nd3proc,nproc,ioption
        real(dp):: zw,zmpi1
       dimension zmpi1(2,n1,nd2proc,nd3proc,nproc),zw(2,lot,n1)
!Local variables-------------------------------
! *************************************************************************
        mfft=0
        if (ioption == 2) then
        do 300,Jp2=Jp2st,nproc
        do 200,J2=J2st,nd2proc
        mfft=mfft+1
        if (mfft.gt.n1dfft) then
        Jp2st=Jp2
        J2st=J2
        return
        end if
        do 100,I1=1,n1
        zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
        zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
100     continue
200     continue
        J2st=1
300     continue
        else
        do 600,Jp2=Jp2st,nproc
        do 500,J2=J2st,nd2proc
        mfft=mfft+1
        if (mfft.gt.n1dfft) then
        Jp2st=Jp2
        J2st=J2
        return
        end if
        ind=(Jp2-1) * nd2proc + J2
        jj2=(ind-1)/nproc +1

        !jjp2=modulo(ind,nproc) +1
        jjp2=modulo(ind-1,nproc)+1

        do 400,I1=1,n1
        zmpi1(1,I1,jj2,j3,jjp2)=zw(1,mfft,I1)
        zmpi1(2,I1,jj2,j3,jjp2)=zw(2,mfft,I1)
400     continue
500     continue
        J2st=1
600     continue
        end if

end subroutine unmpiswitch
