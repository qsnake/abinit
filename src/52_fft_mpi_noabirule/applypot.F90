#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

        subroutine applypot(icplexwf,icplex,ndat,n1,n2,n3,nd1,nd2,nd3proc,&
&        max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,&
&        max1o,max2o,max3o,m1o,m2o,m3o,mpi_comm,nproc,iproc,paral_kgb,pot,zf)

 use m_profiling
          use defs_basis
          use defs_fftdata
! Applies the local real space potential to multiple wavefunctions in Fourier space
!          ZF: Wavefunction (input/output) (note the switch of i2 and i3)
!               real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!               imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!          max1 is positive or zero ; m1 >=max1+1
!          i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!          then, if m1 > max1+1, one has min1=max1-m1+1 and
!          i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!          i2 and i3 have a similar definition of range
!          idat=1,ndat
!          md1,md2,md3: Dimension of ZF (input as well as output), distributed on different procs
!          md2proc=((md2-1)/nproc)+1  maximal number of small box 2nd dim slices for one proc
!
!          POT: Potential
!               POT(icplex*i1,i2,i3)
!               icplex=1 or 2 ,  i1=1,n1 , i2=1,n2 , i3=1,n3
!          nd1,nd2,nd3: dimension of pot
!
!          mpi_comm: MPI communicator
!          nproc: number of processors used as returned by MPI_COMM_SIZE
!          iproc: [0:nproc-1] number of processor as returned by MPI_COMM_RANK
!           n1,n2,n3: logical dimension of the transform. As transform lengths
!                     most products of the prime factors 2,3,5 are allowed.
!                    The detailed table with allowed transform lengths can
!                    be found in subroutine CTRIG
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
#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'applypot'
 use interfaces_52_fft_mpi_noabirule, except_this_one => applypot
!End of the abilint section

        implicit real(dp) (a-h,o-z)
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
! real space input
        integer :: icplexwf,icplex,ndat,n1,n2,n3,nd1,nd2,nd3proc
        integer :: paral_kgb,max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3
        integer :: max1o,max2o,max3o,m1o,m2o,m3o,mpi_comm,nproc,iproc
        REAL(KIND=DP), DIMENSION(icplex*nd1,nd2,nd3proc) :: pot
! Fourier space output
        REAL(KIND=DP), DIMENSION(2,md1,md3,md2proc,ndat) :: zf

!Local variables-------------------------------
! work arrays for transpositions
        REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:) :: zt
! work arrays for MPI
        REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: zmpi1
        REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: zmpi2
! cache work array
        REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:) :: zw
! FFT work arrays
        REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:) :: btrig1,btrig2,btrig3
        REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:) :: ftrig1,ftrig2,ftrig3
        INTEGER, ALLOCATABLE, DIMENSION(:) :: after1,now1,before1, &
                          after2,now2,before2,after3,now3,before3

! *************************************************************************
    integer unused
        unused=0

!!      interface
!!        integer ( kind=4 ) function omp_get_num_threads ( )
!!        end function omp_get_num_threads
!!      end interface
!!      interface
!!        integer ( kind=4 ) function omp_get_thread_num ( )
!!        end function omp_get_thread_num
!!      end interface

        !write(std_out,*)' applypot : enter '

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
!$omp shared(icplex,ndat,n1,n2,n3,nd1,nd2,nd3proc,iproc,nproc,ncache,pot,zf,lock)&
!$omp shared(max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3)&
!$omp shared(max1o,max2o,max3o,m1o,m2o,m3o,icplexwf)

        iam=0
        npr=1
!!      iam=omp_get_thread_num()
!!      npr=omp_get_num_threads()
!!      if (nproc.gt.1 .and. npr.ne.2 .and. npr.ne.4)  &
!!          stop 'no communication scheduling provided in applypot'

!       Effective m1 and m2 (complex-to-complex or real-to-complex)
        n1eff=n1 ; m2ieff=m2i ; m2oeff=m2o ; m1zt=n1
        if(icplexwf==1)then
         n1eff=(n1+1)/2 ; m2ieff=m2i/2+1 ; m2oeff=m2o/2+1 ;m1zt=2*(n1/2+1)
        end if

        m2eff=max(m2ieff,m2oeff)
        lzt=m2eff
        if (mod(m2eff,2).eq.0) lzt=lzt+1
        if (mod(m2eff,4).eq.0) lzt=lzt+1

        nnd3=nd3proc*nproc        ! maximal number of big box 3rd dim slices for all procs

!$omp critical
        ABI_ALLOCATE(btrig1,(2,n1))
        ABI_ALLOCATE(ftrig1,(2,n1))
        ABI_ALLOCATE(after1,(mdata))
        ABI_ALLOCATE(now1,(mdata))
        ABI_ALLOCATE(before1,(mdata))
        ABI_ALLOCATE(btrig2,(2,n2))
        ABI_ALLOCATE(ftrig2,(2,n2))
        ABI_ALLOCATE(after2,(mdata))
        ABI_ALLOCATE(now2,(mdata))
        ABI_ALLOCATE(before2,(mdata))
        ABI_ALLOCATE(btrig3,(2,n3))
        ABI_ALLOCATE(ftrig3,(2,n3))
        ABI_ALLOCATE(after3,(mdata))
        ABI_ALLOCATE(now3,(mdata))
        ABI_ALLOCATE(before3,(mdata))
        ABI_ALLOCATE(zw,(2,ncache/4,2))
        ABI_ALLOCATE(zt,(2,lzt,m1zt))
        ABI_ALLOCATE(zmpi2,(2,md1,md2proc,nnd3))
        if (nproc.gt.1)  then
          ABI_ALLOCATE(zmpi1,(2,md1,md2proc,nnd3))
        end if
!$omp end critical

        call ctrig(n3,btrig3,after3,before3,now3,1,ic3)
        call ctrig(n1,btrig1,after1,before1,now1,1,ic1)
        call ctrig(n2,btrig2,after2,before2,now2,1,ic2)
        do  j=1,n1
        ftrig1(1,j)= btrig1(1,j)
        ftrig1(2,j)=-btrig1(2,j)
        end do
        do  j=1,n2
        ftrig2(1,j)= btrig2(1,j)
        ftrig2(2,j)=-btrig2(2,j)
        end do
        do  j=1,n3
        ftrig3(1,j)= btrig3(1,j)
        ftrig3(2,j)=-btrig3(2,j)
        end do

        ndatmod=((ndat+npr-1)/npr)*npr
!        write(std_out,*) 'ndatmod',ndatmod
        do 12345,idat=1,ndatmod
        if (mod(idat-1,npr).eq.iam) then
!        write(std_out,*) 'IN  iam,idat',iam,idat
        if (idat.le.ndat) then

! transform along z axis
! input: I1,I3,J2,(Jp2)

        lot=ncache/(4*n3)

        do 3331,j2=1,md2proc
        if (iproc*md2proc+j2.le.m2ieff) then

        do 3001,i1=1,m1i,lot
        ma=i1
        mb=min(i1+(lot-1),m1i)
        n1dfft=mb-ma+1

!  input: I1,I3,J2,(Jp2)
        call fill_cent(md1,md3,lot,n1dfft,max3i,m3i,n3,zf(1,i1,1,j2,idat),zw(1,1,1))

        inzee=1
        do i=1,ic3
        call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
                    btrig3,after3(i),now3(i),before3(i),1)
            inzee=3-inzee
        end do

!  input: I1,i3,J2,(Jp2)
        call scramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zw(1,1,inzee),zmpi2)
!  output: I1,J2,i3,(Jp2)

3001        continue
        end if
3331        continue
        end if

! Interprocessor data transposition
! input: I1,J2,j3,jp3,(Jp2)
        if (nproc.gt.1) then
11        continue
!$omp   flush(lock)
! communication scheduling
          if ( (npr.eq.4 .and. iam.eq.0 .and. lock.eq.0) .or.  &
               (npr.eq.4 .and. iam.eq.1 .and. lock.eq.1) .or.  &
               (npr.eq.4 .and. iam.eq.2 .and. lock.eq.4) .or.  &
               (npr.eq.4 .and. iam.eq.3 .and. lock.eq.5) .or.  &
               (npr.eq.2 .and. iam.eq.0 .and. lock.eq.0) .or.  &
               (npr.eq.2 .and. iam.eq.1 .and. lock.eq.2) .or.  &
              npr.eq.1                                 ) then
!        write(std_out,'(a,6(x,i3))') 'Applypot 1 ALLTOALL',nproc,iproc,npr,iam,lock,idat
#if defined HAVE_MPI
  if(paral_kgb == 1) then
        if (idat.le.ndat) &
        call MPI_ALLTOALL(zmpi2,md1*md2proc*nd3proc, &
                          MPI_double_precision, &
                          zmpi1,md1*md2proc*nd3proc, &
                          MPI_double_precision,mpi_comm,ierr)
  endif
#endif
        lock=lock+1
          if (npr.eq.4 .and. idat.eq.2) lock=lock+2
          if (npr.eq.2 .and. idat.eq.1) lock=lock+1
          lock=mod(lock,2*npr)
!        write(std_out,'(a,6(x,i3))') 'new lock 1', nproc,iproc,npr,iam,lock
!$omp   flush(lock)
          else
          goto 11
          end if
! output: I1,J2,j3,Jp2,(jp3)
        end if


        if (idat.le.ndat) then
        do 1212,j3=1,nd3proc
        if (iproc*nd3proc+j3.le.n3) then

        Jp2stb=1
        J2stb=1
        Jp2stf=1
        J2stf=1

! transform along x axis
        lot=ncache/(4*n1)

        do 1001,j=1,m2ieff,lot
        ma=j
        mb=min(j+(lot-1),m2ieff)
        n1dfft=mb-ma+1

! input: I1,J2,j3,Jp2,(jp3)
        ioption=0
        if (nproc.eq.1) then
         call mpiswitch_cent(j3,n1dfft,Jp2stb,J2stb,lot,max1i,md1,m1i,n1,&
&         md2proc,nd3proc,nproc,ioption,zmpi2,zw(1,1,1), unused, unused, unused)
        else
         call mpiswitch_cent(j3,n1dfft,Jp2stb,J2stb,lot,max1i,md1,m1i,n1,&
&         md2proc,nd3proc,nproc,ioption,zmpi1,zw(1,1,1), unused, unused, unused)
        end if
! output: J2,Jp2,I1,j3,(jp3)


! input: I2,I1,j3,(jp3)

        inzee=1
        do i=1,ic1-1
        call fftstp(lot,n1dfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                    btrig1,after1(i),now1(i),before1(i),1)
            inzee=3-inzee
        end do

        i=ic1
        call fftstp(lot,n1dfft,n1,lzt,m1zt,zw(1,1,inzee),zt(1,j,1), &
                   btrig1,after1(i),now1(i),before1(i),1)

! output: I2,i1,j3,(jp3)

1001        continue

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
         call switch_cent(n1dfft,max2i,m2i,n2,lot,n1,lzt,zt(1,1,jeff),zw(1,1,1))
        else
         call switchreal_cent(includelast,n1dfft,max2i,n2,lot,m1zt,lzt,zt(1,1,jeff),zw(1,1,1))
        end if
! output: i1,I2,j3,(jp3)

        inzee=1
        do i=1,ic2
        call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                    btrig2,after2(i),now2(i),before2(i),1)
            inzee=3-inzee
        end do
! output: i1,i2,j3,(jp3)

!Multiply with potential in real space
        jx=icplex*(jeff-1)+1
        call multpot(icplexwf,icplex,includelast,nd1,nd2,n2,lot,n1dfft,pot(jx,1,j3),zw(1,1,inzee))

! TRANSFORM BACK IN FOURIER SPACE
! transform along y axis
! input: i1,i2,j3,(jp3)
        do i=1,ic2
        call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                    ftrig2,after2(i),now2(i),before2(i),-1)
            inzee=3-inzee
        end do

!  input: i1,I2,j3,(jp3)
        if(icplexwf==2)then
         call unswitch_cent(n1dfft,max2o,m2o,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,jeff))
        else
         call unswitchreal_cent(n1dfft,max2o,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,jeff))
        end if
! output: I2,i1,j3,(jp3)

2000        continue

! transform along x axis
! input: I2,i1,j3,(jp3)

        lot=ncache/(4*n1)

        do 1002,j=1,m2oeff,lot

        ma=j
        mb=min(j+(lot-1),m2oeff)
        n1dfft=mb-ma+1
        i=1
        call fftstp(lzt,n1dfft,m1zt,lot,n1,zt(1,j,1),zw(1,1,1), &
                   ftrig1,after1(i),now1(i),before1(i),-1)

        inzee=1
        do i=2,ic1
        call fftstp(lot,n1dfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                    ftrig1,after1(i),now1(i),before1(i),-1)
            inzee=3-inzee
        end do

! output: I2,I1,j3,(jp3)

! input: J2,Jp2,I1,j3,(jp3)
        ioption=0
        if (nproc.eq.1) then
         call unmpiswitch_cent(j3,n1dfft,Jp2stf,J2stf,lot,max1o,md1,m1o,n1,&
&         md2proc,nd3proc,nproc,ioption,zw(1,1,inzee),zmpi2)
        else
         call unmpiswitch_cent(j3,n1dfft,Jp2stf,J2stf,lot,max1o,md1,m1o,n1,&
&         md2proc,nd3proc,nproc,ioption,zw(1,1,inzee),zmpi1)
        end if

! output: I1,J2,j3,Jp2,(jp3)

1002        continue

        end if
1212        continue
        end if

! Interprocessor data transposition
! intput: I1,J2,j3,Jp2,(jp3)
        if (nproc.gt.1) then
22      continue
!$omp   flush(lock)
! communication scheduling
        if ( (npr.eq.4 .and. iam.eq.0 .and. lock.eq.6) .or.  &
             (npr.eq.4 .and. iam.eq.1 .and. lock.eq.7) .or.  &
             (npr.eq.4 .and. iam.eq.2 .and. lock.eq.2) .or.  &
             (npr.eq.4 .and. iam.eq.3 .and. lock.eq.3) .or.  &
             (npr.eq.2 .and. iam.eq.0 .and. lock.eq.3) .or.  &
             (npr.eq.2 .and. iam.eq.1 .and. lock.eq.1) .or.  &
              npr.eq.1                                 ) then
!        write(std_out,'(a,6(x,i3))') 'Applypot 2 ALLTOALL',nproc,iproc,npr,iam,lock,idat
#if defined HAVE_MPI
  if(paral_kgb == 1)then
        if (idat.le.ndat) &
        call MPI_ALLTOALL(zmpi1,n1*md2proc*nd3proc, &
                          MPI_double_precision, &
                          zmpi2,n1*md2proc*nd3proc, &
                          MPI_double_precision,mpi_comm,ierr)
  endif
#endif
        lock=lock+1
          if (npr.eq.4 .and. idat.eq.ndatmod-2) lock=lock+2
          if (npr.eq.2 .and. idat.eq.ndatmod-1) lock=lock+1
          lock=mod(lock,2*npr)
!        write(std_out,'(a,6(x,i3))') 'new lock 2', nproc,iproc,npr,iam,lock
!$omp   flush(lock)
          else
          goto 22
          end if
! output: I1,J2,j3,jp3,(Jp2)
        end if

! transform along z axis
! input: I1,J2,i3,(Jp2)

        if (idat.le.ndat) then
        lot=ncache/(4*n3)

        do 3332,j2=1,md2proc
        if (iproc*md2proc+j2.le.m2oeff) then

        do 3002,i1=1,m1o,lot
        ma=i1
        mb=min(i1+(lot-1),m1o)
        n1dfft=mb-ma+1

!   input: I1,J2,i3,(Jp2)
        call unscramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zmpi2,zw(1,1,1))
! output: I1,i3,J2,(Jp2)

        inzee=1
        do i=1,ic3
        call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
        ftrig3,after3(i),now3(i),before3(i),-1)
        inzee=3-inzee
        end do

        call unfill_cent(md1,md3,lot,n1dfft,max3o,m3o,n3,zw(1,1,inzee),zf(1,i1,1,j2,idat))
! output: I1,I3,J2,(Jp2)

3002        continue
        end if
3332    continue
        end if
        end if

!       Complete missing values with complex conjugate
!       Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.

        if(icplexwf==1)then
         do i3=1,m3o
          i3inv=m3o+2-i3
          if(i3==1)i3inv=1
          if(m2oeff>1)then
           do i2=2,m2oeff
            i2inv=m2o+2-i2
            zf(1,1,i3inv,i2inv,idat)= zf(1,1,i3,i2,idat)
            zf(2,1,i3inv,i2inv,idat)=-zf(2,1,i3,i2,idat)
            do i1=2,m1o
             i1inv=m1o+2-i1
             zf(1,i1inv,i3inv,i2inv,idat)= zf(1,i1,i3,i2,idat)
             zf(2,i1inv,i3inv,i2inv,idat)=-zf(2,i1,i3,i2,idat)
            end do
           end do
          end if
         end do
        end if

12345   continue

        ABI_DEALLOCATE(btrig1)
        ABI_DEALLOCATE(ftrig1)
        ABI_DEALLOCATE(after1)
        ABI_DEALLOCATE(now1)
        ABI_DEALLOCATE(before1)
        ABI_DEALLOCATE(btrig2)
        ABI_DEALLOCATE(ftrig2)
        ABI_DEALLOCATE(after2)
        ABI_DEALLOCATE(now2)
        ABI_DEALLOCATE(before2)
        ABI_DEALLOCATE(btrig3)
        ABI_DEALLOCATE(ftrig3)
        ABI_DEALLOCATE(after3)
        ABI_DEALLOCATE(now3)
        ABI_DEALLOCATE(before3)
        ABI_DEALLOCATE(zmpi2)
        ABI_DEALLOCATE(zw)
        ABI_DEALLOCATE(zt)
        if (nproc.gt.1)  then
          ABI_DEALLOCATE(zmpi1)
        end if
!$omp end parallel

        !write(std_out,*)' applypot : exit '

        return
end subroutine applypot
