#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

        subroutine accrho(icplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc,&
&        max1,max2,max3,m1,m2,m3,md1,md2proc,md3,mpi_comm,nproc,iproc,paral_kgb,zf,rho,weight)

 use m_profiling
! Accumulates the real space density rho from the ndat wavefunctions zf
! by transforming zf into real space and adding all the amplitudes squared
!        INPUT:
!          ZF: input array (note the switch of i2 and i3)
!                real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!                imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!          max1 is positive or zero ; m1 >=max1+1
!          i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!          then, if m1 > max1+1, one has min1=max1-m1+1 and
!          i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!          i2 and i3 have a similar definition of range
!          idat=1,ndat
!          md1,md2,md3: Dimension of ZF
!          md2proc=((md2-1)/nproc)+1 ! maximal number of small box 2nd dim slices for one proc
!        OUTPUT:
!           RHOoutput(i1,i2,i3) = RHOinput(i1,i2,i3) + sum on idat of (FFT(ZF))**2 *weight
!               i1=1,n1 , i2=1,n2 , i3=1,n3
!          mpi_comm: MPI communicator
!          nproc: number of processors used as returned by MPI_COMM_SIZE
!          iproc: [0:nproc-1] number of processor as returned by MPI_COMM_RANK
!           n1,n2,n3: logical dimension of the transform. As transform lengths
!                     most products of the prime factors 2,3,5 are allowed.
!                    The detailed table with allowed transform lengths can
!                    be found in subroutine CTRIG
!           nd1,nd2,nd3: Dimension of RHO
!          nd3proc=((nd3-1)/nproc)+1 ! maximal number of big box 3rd dim slices for one proc
!
!       PERFORMANCE CONSIDERATIONS:
!       The maximum number of processors that can reasonably be used is max(n2/2,n3/2)
!
!       It is very important to find the optimal
!       value of NCACHE. NCACHE determines the size of the work array ZW, that
!       has to fit into cache. It has therefore to be chosen to equal roughly
!        half the size of the physical cache in units of real*8 numbers.
!       The optimal value of ncache can easily be determined by numerical
!       experimentation. A too large value of ncache leads to a dramatic
!       and sudden decrease of performance, a too small value to a to a
!       slow and less dramatic decrease of performance. If NCACHE is set
!       to a value so small, that not even a single one dimensional transform
!       can be done in the workarray zw, the program stops with an error message.
!
!       RESTRICTIONS on USAGE
!  Copyright (C) 2002-2007 Stefan Goedecker, CEA Grenoble
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
        use defs_basis
        use defs_fftdata
#if defined HAVE_MPI2
        use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'accrho'
 use interfaces_52_fft_mpi_noabirule, except_this_one => accrho
!End of the abilint section

        implicit real(dp) (a-h,o-z)
#if defined HAVE_MPI1
 include 'mpif.h'
#endif
!Arguments ------------------------------------
! real space input
        integer            :: icplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc,max1,max2,max3,m1,m2,m3
        integer            :: md1,md2proc,md3,mpi_comm,nproc,iproc,paral_kgb
        REAL(DP), DIMENSION(nd1,nd2,nd3proc) :: rho
! Fourier space output
        REAL(DP), DIMENSION(2,md1,md3,md2proc,ndat) :: zf
! weight for the density accumulation
        REAL(DP), DIMENSION(ndat) :: weight

!Local variables-------------------------------
        REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)   :: rhopart
! work arrays for MPI
        REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: zmpi1
        REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: zmpi2
! work arrays for transpositions
        REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)   :: zt

! cache work array
        REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)   :: zw
! FFT work arrays
        REAL(DP), ALLOCATABLE, DIMENSION(:,:)     :: trig1,trig2,trig3
        INTEGER, ALLOCATABLE, DIMENSION(:)         :: after1,now1,before1, &
             & after2,now2,before2,after3,now3,before3

        integer unused
        unused=0
! *************************************************************************
!!      interface
!!        integer ( kind=4 ) function omp_get_num_threads ( )
!!        end function omp_get_num_threads
!!      end interface
!!      interface
!!        integer ( kind=4 ) function omp_get_thread_num ( )
!!        end function omp_get_thread_num
!!      end interface

        !write(std_out,*)' accrho : enter '

! find cache size that gives optimal performance on machine
        ncache=4*1024
        if (ncache/(4*max(n1,n2,n3)).lt.1) then
                      write(std_out,*) &
&                        ' ncache has to be enlarged to be able to hold at', &
&                        ' least one 1-d FFT of each size even though this will', &
&                        ' reduce the performance for shorter transform lengths'
                       stop
        end if

        lock=0
!$omp parallel  default(private) &
!$omp shared(ndat,n1,n2,n3,nd1,nd2,nd3proc,md1,md2proc,md3,iproc,nproc,ncache,rho,zf,lock,weight)&
!$omp shared(max1,max2,max3,m1,m2,m3,icplexwf)

        iam=0
        npr=1
!!      iam=omp_get_thread_num()
!!      npr=omp_get_num_threads()

!       Effective m1 and m2 (complex-to-complex or real-to-complex)
        n1eff=n1 ; m2eff=m2 ; m1zt=n1
        if(icplexwf==1)then
         n1eff=(n1+1)/2 ; m2eff=m2/2+1 ; m1zt=2*(n1/2+1)
        end if

        lzt=m2eff
        if (mod(m2eff,2).eq.0) lzt=lzt+1
        if (mod(m2eff,4).eq.0) lzt=lzt+1

        nnd3=nd3proc*nproc        ! maximal number of big box 3rd dim slices for all procs

!$omp critical
        ABI_ALLOCATE(trig1,(2,n1))
        ABI_ALLOCATE(after1,(mdata))
        ABI_ALLOCATE(now1,(mdata))
        ABI_ALLOCATE(before1,(mdata))
        ABI_ALLOCATE(trig2,(2,n2))
        ABI_ALLOCATE(after2,(mdata))
        ABI_ALLOCATE(now2,(mdata))
        ABI_ALLOCATE(before2,(mdata))
        ABI_ALLOCATE(trig3,(2,n3))
        ABI_ALLOCATE(after3,(mdata))
        ABI_ALLOCATE(now3,(mdata))
        ABI_ALLOCATE(before3,(mdata))
        ABI_ALLOCATE(zw,(2,ncache/4,2))
        ABI_ALLOCATE(zt,(2,lzt,m1zt))
        ABI_ALLOCATE(zmpi2,(2,md1,md2proc,nnd3))
        if (nproc.gt.1)  then
          ABI_ALLOCATE(zmpi1,(2,md1,md2proc,nnd3))
        end if
        if (npr.gt.1)  then
          ABI_ALLOCATE(rhopart,(nd1,nd2,nd3proc))
        end if
!$omp end critical
        if (npr.gt.1) rhopart(:,:,:)=0.0_dp

        call ctrig(n3,trig3,after3,before3,now3,1,ic3)
        call ctrig(n1,trig1,after1,before1,now1,1,ic1)
        call ctrig(n2,trig2,after2,before2,now2,1,ic2)

!$omp do
        do 12345,idat=1,ndat

! transform along z axis
! input: I1,I3,J2,(Jp2)

        lot=ncache/(4*n3)

        do 3333,j2=1,md2proc
        if (iproc*md2proc+j2.le.m2eff) then

        do 3000,i1=1,m1,lot
        ma=i1
        mb=min(i1+(lot-1),m1)
        n1dfft=mb-ma+1

!  input: I1,I3,J2,(Jp2)
        call fill_cent(md1,md3,lot,n1dfft,max3,m3,n3,zf(1,i1,1,j2,idat),zw(1,1,1))

        inzee=1
        do i=1,ic3
        call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
                    trig3,after3(i),now3(i),before3(i),1)
            inzee=3-inzee
        end do

!  input: I1,i3,J2,(Jp2)
        call scramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zw(1,1,inzee),zmpi2)
!  output: I1,J2,i3,(Jp2)

3000        continue
        end if
3333        continue

! Interprocessor data transposition
! input: I1,J2,j3,jp3,(Jp2)
        if (nproc.gt.1) then
11        continue
!$omp   flush(lock)
        if (mod(lock,npr).ne.iam) goto 11
#if defined HAVE_MPI
 if(paral_kgb == 1) then
        call MPI_ALLTOALL(zmpi2,md1*md2proc*nd3proc, &
                          MPI_double_precision, &
                          zmpi1,md1*md2proc*nd3proc, &
                          MPI_double_precision,mpi_comm,ierr)
 endif
#endif
        lock=lock+1
!$omp   flush(lock)
! output: I1,J2,j3,Jp2,(jp3)
        end if

!DEBUG
!       write(std_out,*)' zmpi2 ='
!       do i3=1,nnd3
!        do i2=1,m2eff
!         do i1=1,md1
!          write(std_out,'(3i4,2es16.6)')i1,i2,i3,zmpi2(1:2,i1,i2,i3)
!         end do
!        end do
!       end do
!       stop
!ENDDEBUG

        do 1212,j3=1,nd3proc
        if (iproc*nd3proc+j3.le.n3) then
        Jp2st=1
        J2st=1

! transform along x axis
        lot=ncache/(4*n1)

        do 1000,j=1,m2eff,lot
        ma=j
        mb=min(j+(lot-1),m2eff)
        n1dfft=mb-ma+1

! input: I1,J2,j3,Jp2,(jp3)
        ioption=0
        if (nproc.eq.1) then
         call mpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot,max1,md1,m1,n1,&
&         md2proc,nd3proc,nproc,ioption,zmpi2,zw(1,1,1),unused, unused,unused)
        else
         call mpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot,max1,md1,m1,n1,&
&         md2proc,nd3proc,nproc,ioption,zmpi1,zw(1,1,1), unused,unused,unused)
        end if
! output: J2,Jp2,I1,j3,(jp3)

! input: I2,I1,j3,(jp3)
        inzee=1
        do i=1,ic1-1
        call fftstp(lot,n1dfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                    trig1,after1(i),now1(i),before1(i),1)
            inzee=3-inzee
        end do

        i=ic1
        call fftstp(lot,n1dfft,n1,lzt,m1zt,zw(1,1,inzee),zt(1,j,1), &
                   trig1,after1(i),now1(i),before1(i),1)

! output: I2,i1,j3,(jp3)

1000        continue
! transform along y axis
        lot=ncache/(4*n2)
        if(icplexwf==1)then
         if(mod(lot,2).ne.0)lot=lot-1 ! needed to introduce jeff
        end if

        do 2000,j=1,n1eff,lot
        ma=j
        mb=min(j+(lot-1),n1eff)
        n1dfft=mb-ma+1
        jeff=j
        includelast=1
        if(icplexwf==1)then
         jeff=2*j-1
         includelast=1
         if(mb==n1eff .and. n1eff*2/=n1)includelast=0
        end if


!  input: I2,i1,j3,(jp3)
        if(icplexwf==2)then
         call switch_cent(n1dfft,max2,m2,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
        else
         call switchreal_cent(includelast,n1dfft,max2,n2,lot,m1zt,lzt,zt(1,1,jeff),zw(1,1,1))
        end if
! output: i1,I2,j3,(jp3)

        inzee=1
        do i=1,ic2
        call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                    trig2,after2(i),now2(i),before2(i),1)
            inzee=3-inzee
        end do

!Accumulate
        if (npr.eq.1) then
        call addrho(icplexwf,includelast,nd1,nd2,n2,lot,n1dfft,&
&        zw(1,1,inzee),rho(jeff,1,j3),weight(idat))
        else
        call addrho(icplexwf,includelast,nd1,nd2,n2,lot,n1dfft,&
&        zw(1,1,inzee),rhopart(jeff,1,j3),weight(idat))
        end if

2000        continue
! output: i1,i2,j3,(jp3)

!DEBUG
!       write(std_out,*)' j3=1, rho ='
!       do i1=1,n1
!        do i2=1,n2
!         write(std_out,'(2i4,2es16.6)')i1,i2,rho(i1,i2,j3)
!        end do
!       end do
!       stop
!ENDDEBUG

        end if
1212        continue

12345   continue
!$omp end do

!$omp critical
! Sum total density from the partial densities per thread
        if (npr.gt.1) then
        do 644, i3=1,nd3proc
        j3=iproc*nd3proc+i3
        if (j3.le.n3) then
        do 643, i2=1,n2
        do 643, i1=1,n1
        rho(i1,i2,i3)=rho(i1,i2,i3)+rhopart(i1,i2,i3)
643        continue
        end if
644        continue
        end if
!$omp end critical

        ABI_DEALLOCATE(trig1)
        ABI_DEALLOCATE(after1)
        ABI_DEALLOCATE(now1)
        ABI_DEALLOCATE(before1)
        ABI_DEALLOCATE(trig2)
        ABI_DEALLOCATE(after2)
        ABI_DEALLOCATE(now2)
        ABI_DEALLOCATE(before2)
        ABI_DEALLOCATE(trig3)
        ABI_DEALLOCATE(after3)
        ABI_DEALLOCATE(now3)
        ABI_DEALLOCATE(before3)
        ABI_DEALLOCATE(zmpi2)
        ABI_DEALLOCATE(zw)
        ABI_DEALLOCATE(zt)
        if (nproc.gt.1)  then
          ABI_DEALLOCATE(zmpi1)
        end if
        if (npr.gt.1)  then
          ABI_DEALLOCATE(rhopart)
        end if
!$omp end parallel

        !write(std_out,*)' accrho : exit '

        return
end subroutine accrho
