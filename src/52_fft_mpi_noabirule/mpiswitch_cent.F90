#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


        subroutine mpiswitch_cent(j3,n1dfft,Jp2stb,J2stb,lot,max1,md1,m1,n1,md2proc,&
        &nd3proc,nproc,ioption,zmpi1,zw,max2,m2,n2)

 use m_profiling
        use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mpiswitch_cent'
!End of the abilint section

        implicit real(dp) (a-h,o-z)
!Arguments ------------------------------------
        integer :: j3,n1dfft,Jp2stb,J2stb,lot,max1,md1,m1,n1,md2proc,nd3proc,nproc
        real(dp)  :: zmpi1,zw
        dimension zmpi1(2,md1,md2proc,nd3proc,nproc),zw(2,lot,n1)
!Local variables-------------------------------

! *************************************************************************
        if(.false.)write(std_out,*)m2,max2,n2
        mfft=0
        if (ioption /= 1) then
        do 300,Jp2=Jp2stb,nproc
        do 200,J2=J2stb,md2proc
        mfft=mfft+1
        if (mfft.gt.n1dfft) then
        Jp2stb=Jp2
        J2stb=J2
        return
        end if

!       Here, zero and positive frequencies
!       In zmpi1, they are stored from 1 to max1+1
        do 90,I1=1,max1+1
        zw(1,mfft,I1)=zmpi1(1,I1,J2,j3,Jp2)
        zw(2,mfft,I1)=zmpi1(2,I1,J2,j3,Jp2)
90      continue

!       Fill the center region with zeros
        do 100,I1=max1+2,n1-m1+max1+1
        zw(1,mfft,I1)=0.d0
        zw(2,mfft,I1)=0.d0
100     continue

!       Here, negative frequencies
!       In zmpi1, they are stored from 1 to m1half
        do 110,I1=max1+2,m1
        zw(1,mfft,I1+n1-m1)=zmpi1(1,I1,J2,j3,Jp2)
        zw(2,mfft,I1+n1-m1)=zmpi1(2,I1,J2,j3,Jp2)
110     continue

200     continue
        J2stb=1
300     continue
        else
        do 600,Jp2=Jp2stb,nproc
        do 500,J2=J2stb,md2proc
        mfft=mfft+1
        if (mfft.gt.n1dfft) then
        Jp2stb=Jp2
        J2stb=J2
        return
        end if
        ind=(Jp2-1) * md2proc + J2
        jj2=(ind-1)/nproc +1

        !jjp2=modulo(ind,nproc) +1
        jjp2=modulo(ind-1,nproc)+1

!I gather consecutive I2 indexes in mfft in the modulo case
!       Here, zero and positive frequencies
!       In zmpi1, they are stored from 1 to max1+1
        do 390,I1=1,max1+1
        zw(1,mfft,I1)=zmpi1(1,I1,Jj2,j3,Jjp2)
        zw(2,mfft,I1)=zmpi1(2,I1,Jj2,j3,Jjp2)
390     continue

!       Fill the center region with zeros
        do 400,I1=max1+2,n1-m1+max1+1
        zw(1,mfft,I1)=zero
        zw(2,mfft,I1)=zero
400     continue

!       Here, negative frequencies
!       In zmpi1, they are stored from 1 to m1half
        do 410,I1=max1+2,m1
        zw(1,mfft,I1+n1-m1)=zmpi1(1,I1,Jj2,j3,Jjp2)
        zw(2,mfft,I1+n1-m1)=zmpi1(2,I1,Jj2,j3,Jjp2)
410     continue

500     continue
        J2stb=1
600     continue
        end if
        end subroutine mpiswitch_cent
