#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


        subroutine unmpiswitch_cent(j3,n1dfft,Jp2stf,J2stf,lot,max1,md1,m1,n1,md2proc,nd3proc,nproc,ioption,zw,zmpi1)

 use m_profiling
        use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unmpiswitch_cent'
!End of the abilint section

        implicit real(dp) (a-h,o-z)
!Arguments ------------------------------------
        integer :: j3,n1dfft,Jp2stf,J2stf,lot,max1,md1,m1,n1,md2proc,nd3proc,nproc
        real(dp):: zw,zmpi1
        dimension zmpi1(2,md1,md2proc,nd3proc,nproc),zw(2,lot,n1)
!Local variables-------------------------------
! *************************************************************************
        mfft=0
        if (ioption == 2) then
        do 300,Jp2=Jp2stf,nproc
        do 200,J2=J2stf,md2proc
        mfft=mfft+1
        if (mfft.gt.n1dfft) then
        Jp2stf=Jp2
        J2stf=J2
        return
        end if

!       Here, zero and positive frequencies
        do 90,I1=1,max1+1
        zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1)
        zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1)
90      continue

!       Here, negative frequencies
        do 110,I1=max1+2,m1
        zmpi1(1,I1,J2,j3,Jp2)=zw(1,mfft,I1+n1-m1)
        zmpi1(2,I1,J2,j3,Jp2)=zw(2,mfft,I1+n1-m1)
110     continue

200     continue
        J2stf=1
300     continue
        else
        do 600,Jp2=Jp2stf,nproc
        do 500,J2=J2stf,md2proc
        mfft=mfft+1
        if (mfft.gt.n1dfft) then
        Jp2stf=Jp2
        J2stf=J2
        return
        end if
        ind=(Jp2-1) * md2proc + J2
        jj2=(ind-1)/nproc +1

        !jjp2=modulo(ind,nproc) +1
        jjp2=modulo(ind-1,nproc)+1

!       Here, zero and positive frequencies
        do 390,I1=1,max1+1
        zmpi1(1,I1,Jj2,j3,Jjp2)=zw(1,mfft,I1)
        zmpi1(2,I1,Jj2,j3,Jjp2)=zw(2,mfft,I1)
390      continue

!       Here, negative frequencies
        do 410,I1=max1+2,m1
        zmpi1(1,I1,Jj2,j3,Jjp2)=zw(1,mfft,I1+n1-m1)
        zmpi1(2,I1,Jj2,j3,Jjp2)=zw(2,mfft,I1+n1-m1)
410     continue

500     continue
        J2stf=1
600     continue
        end if
        end subroutine unmpiswitch_cent
