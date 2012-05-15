#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

        subroutine back(icplex,mpi_enreg,ndat,n1,n2,n3,nd1,nd2,nd3,nd1eff,nd2proc,nd3proc,option,paral_kgb,zf,zr)

 use m_profiling
! Adopt standard convention that isign=1 for backward transform
!        CALCULATES THE DISCRETE FOURIER TRANSFORM ZR(I1,I2,I3)=
!        S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) ZF(j1,j3,j2)
!       in parallel using MPI/OpenMP and BLAS library calls.
!
!       INPUT:
!          ZF: input array (note the switch of i2 and i3)
!               real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!               imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!               i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!       OUTPUT:
!          ZR: output array
!               ZR(1,i1,i2,i3,idat)=real(R(i1,i2,i3,idat))
!               ZR(2,i1,i2,i3,idat)=imag(R(i1,i2,i3,idat))
!               i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!          mpi_enreg%nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!          mpi_enreg%me_fft: [0:mpi_enreg%nproc_fft-1] number of processor as returned by MPI_COMM_RANK
!           n1,n2,n3: logical dimension of the transform. As transform lengths
!                     most products of the prime factors 2,3,5 are allowed.
!                    The detailed table with allowed transform lengths can
!                    be found in subroutine CTRIG
!          nd1,nd2,nd3: Dimension of ZF and ZR
!          nd2proc=((nd2-1)/mpi_enreg%nproc_fft)+1 maximal number of 2nd dim slices
!          nd3proc=((nd3-1)/mpi_enreg%nproc_fft)+1 maximal number of 3rd dim slices
!
!       PERFORMANCE CONSIDERATIONS:
!       The maximum number of processors that can reasonably be used is max(n2,n3)
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
 use defs_abitypes
 use defs_fftdata
#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'back'
 use interfaces_51_manage_mpi
 use interfaces_52_fft_mpi_noabirule, except_this_one => back
!End of the abilint section

        implicit real(dp) (a-h,o-z)
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
! real space input
        integer,intent(in) :: paral_kgb,icplex
        type(MPI_type),intent(inout) :: mpi_enreg
        integer :: ndat,n1,n2,n3,nd1,nd2,nd3,nd1eff,nd2proc,nd3proc,option
        REAL(DP), DIMENSION(2,nd1eff,nd2,nd3proc,ndat) :: zr
! Fourier space output
        REAL(DP), DIMENSION(2,nd1,nd3,nd2proc,ndat) :: zf
!Local variables-------------------------------
        integer :: nproc_fft,me_fft
! work arrays for transpositions
        REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: zt
! work arrays for MPI
        REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: zmpi1
        REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: zmpi2
! cache work array
        REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: zw
! FFT work arrays
        REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: trig1,trig2,trig3
        INTEGER, ALLOCATABLE, DIMENSION(:) :: after1,now1,before1, &
                          after2,now2,before2,after3,now3,before3
        integer :: spaceComm
#if defined HAVE_MPI
        integer :: old_paral_level
#endif
! *************************************************************************

!!      interface
!!        integer ( kind=4 ) function omp_get_num_threads ( )
!!        end function omp_get_num_threads
!!      end interface
!!      interface
!!        integer ( kind=4 ) function omp_get_thread_num ( )
!!        end function omp_get_thread_num
!!      end interface

!DEBUG
!       write(std_out,*)' back : enter '
!ENDDEBUG
        nproc_fft = mpi_enreg%nproc_fft
        me_fft    = mpi_enreg%me_fft

! find cache size that gives optimal performance on machine
        ncache=4*1024

        if (ncache/(4*max(n1,n2,n3)).lt.1) then
                      write(std_out,*) &
&                        ' ncache has to be enlarged to be able to hold at', &
&                        ' least one 1-d FFT of each size even though this will', &
&                        ' reduce the performance for shorter transform lengths'
                       stop
        end if

! check input
        if (nd1.lt.n1) stop 'ERROR:nd1'
        if (nd2.lt.n2) stop 'ERROR:nd2'
        if (nd3.lt.n3) stop 'ERROR:nd3'

        lock=0

! Disabled by MG on Dec  6 2011, omp sections have to be tested, this coding causes a 
! sigfault with nthreads==1. Besides I don't understand why we parallelize over ndata that is usually 1!
! $omp parallel  default(private) &
! $omp shared(ndat,n1,n2,n3,nd1,nd2,nd3,nd2proc,nd3proc,me_fft,nproc_fft,ncache,zr,zf,lock,icplex,paral_kgb)

        iam=0
        npr=1
!!      iam=omp_get_thread_num()
!!      npr=omp_get_num_threads()

!       Effective n1 and n2 (complex-to-complex or real-to-complex)
        n1eff=n1 ; n2eff=n2 ; n1zt=n1
        if(icplex==1)then
         n1eff=(n1+1)/2 ; n2eff=n2/2+1 ; n1zt=2*(n1/2+1)
        end if

        lzt=n2eff
        if (mod(n2eff,2).eq.0) lzt=lzt+1
        if (mod(n2eff,4).eq.0) lzt=lzt+1

        nnd3=nd3proc*nproc_fft        ! maximal number of big box 3rd dim slices for all procs

! $omp critical
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
        ABI_ALLOCATE(zt,(2,lzt,n1zt))
        ABI_ALLOCATE(zmpi2,(2,n1,nd2proc,nnd3))
        if (nproc_fft.gt.1)  then
          ABI_ALLOCATE(zmpi1,(2,n1,nd2proc,nnd3))
        end if
! $omp end critical

        call ctrig(n3,trig3,after3,before3,now3,1,ic3)
        call ctrig(n1,trig1,after1,before1,now1,1,ic1)
        call ctrig(n2,trig2,after2,before2,now2,1,ic2)
!DEBUG
!       write(std_out,'(a,3i4)' )'back,zf n1,n2,n3',n1,n2,n3
!       write(std_out,'(a,3i4)' )'nd1,nd2,nd3proc',nd1,nd2,nd3proc
!      write(std_out,'(a,3i4)' )'m1,m2,m3',m1,m2,m3
!        write(std_out,'(a,3i4)' )'max1,max2,max3',max1,max2,max3
!        write(std_out,'(a,3i4)' )'md1,md2proc,md3',md1,md2proc,md3
!        write(std_out,'(a,3i4)' )'n1eff,m2eff,m1zt',n1eff,m2eff,m1zt
!        do i2=1,nd2proc
!         do i3=1,n3
!          do i1=1,n1
!           write(std_out,'(3i4,2es16.6)')i1,i3,i2,zf(1:2,i1,i3,i2)
!          end do
!        end do
!        end do
!       stop
!ENDDEBUG

! MG this is not very efficient since ndata=1 in many cases
! $omp do
        do 12345,idat=1,ndat

! transform along z axis
! input: I1,I3,J2,(Jp2)


        lot=ncache/(4*n3)

        do 3333,j2=1,nd2proc
        if (me_fft*nd2proc+j2.le.n2eff) then

        do 3000,i1=1,n1,lot
        ma=i1
        mb=min(i1+(lot-1),n1)
        n1dfft=mb-ma+1

!  input: I1,I3,J2,(Jp2)
        call fill(nd1,nd3,lot,n1dfft,n3,zf(1,i1,1,j2,idat),zw(1,1,1))

        inzee=1
        do i=1,ic3
        call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
                    trig3,after3(i),now3(i),before3(i),1)
            inzee=3-inzee
        end do

!  input: I1,i3,J2,(Jp2)
        call scramble(i1,j2,lot,n1dfft,n1,n3,nd2proc,nd3,zw(1,1,inzee),zmpi2)
!  output: I1,J2,i3,(Jp2)

3000        continue
        end if
3333        continue
!DEBUG
!       write(std_out,*)' back zmpi2 ='
!       do i2=1,nd2proc
!        do i1=1,n1
!         do i3=1,n3
!          write(std_out,'(3i4,2es16.6)')i1,i2,i3,zmpi2(1:2,i1,i2,i3)
!         end do
!        end do
!       end do
!ENDDEBUG

! Interprocessor data transposition
! input: I1,J2,j3,jp3,(Jp2)
        if (nproc_fft.gt.1) then
11        continue
! $omp   flush(lock)
        if (mod(lock,npr).ne.iam) goto 11
#if defined HAVE_MPI
  if(paral_kgb == 1) then
        old_paral_level=mpi_enreg%paral_level
        mpi_enreg%paral_level=3
        call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
!DEBUG
! write(std_out,*)' back : will call mpi_alltoall, MPI_double_precision=',MPI_double_precision
!ENDDEBUG
        call MPI_ALLTOALL(zmpi2,2*n1*nd2proc*nd3proc, &
                          MPI_double_precision, &
                          zmpi1,2*n1*nd2proc*nd3proc, &
                          MPI_double_precision,spaceComm,ierr)
!DEBUG
! write(std_out,*)' back : called mpi_alltoall, MPI_double_precision=',MPI_double_precision
!ENDDEBUG

        mpi_enreg%paral_level=old_paral_level
 endif
#endif
        lock=lock+1
! $omp   flush(lock)
! output: I1,J2,j3,Jp2,(jp3)
        end if

        do 1212,j3=1,nd3proc
        if (me_fft*nd3proc+j3.le.n3) then
        Jp2st=1
        J2st=1

!DEBUG
!       if(j3==1)then
!        write(std_out,*)' back : j3==1, zmpi2 '
!        do i2=1,n2eff
!         do i1=1,n1
!          write(std_out,'(2i4,2es16.6)' ) i1,i2,zmpi2(1:2,i1,i2)
!         end do
!        end do
!       end if
!ENDDEBUG


! transform along x axis
        lot=ncache/(4*n1)

        do 1000,j=1,n2eff,lot
        ma=j
        mb=min(j+(lot-1),n2eff)
        n1dfft=mb-ma+1

! input: I1,J2,j3,Jp2,(jp3)
        if (nproc_fft.eq.1) then
         call mpiswitch(j3,n1dfft,Jp2st,J2st,lot,n1,nd2proc,nd3proc,nproc_fft,option,zmpi2,zw(1,1,1))
        else
         call mpiswitch(j3,n1dfft,Jp2st,J2st,lot,n1,nd2proc,nd3proc,nproc_fft,option,zmpi1,zw(1,1,1))
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
        call fftstp(lot,n1dfft,n1,lzt,n1zt,zw(1,1,inzee),zt(1,j,1), &
                   trig1,after1(i),now1(i),before1(i),1)

! output: I2,i1,j3,(jp3)

1000        continue

!DEBUG
!       if(j3==1)then
!        write(std_out,*)' back : j3==1, zt ',j3
!        do i1=1,n1
!         do i2=1,n2
!          write(std_out,'(2i4,2es16.6)' ) i2,i1,zt(1:2,i2,i1)
!         end do
!        end do
!       end if
!ENDDEBUG

! transform along y axis
        lot=ncache/(4*n2)

        do 2000,j=1,n1eff,lot
        ma=j
        mb=min(j+(lot-1),n1eff)
        n1dfft=mb-ma+1
        includelast=1
        if(icplex==1)then
         jeff=2*j-1
         includelast=1
         if(mb==n1eff .and. n1eff*2/=n1)includelast=0
        end if

!  input: I2,i1,j3,(jp3)
        if(icplex==2)then
          call switch(n1dfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
        else
         call switchreal(includelast,n1dfft,n2,n2eff,lot,n1zt,lzt,zt(1,1,jeff),zw(1,1,1))
        end if
! output: i1,I2,j3,(jp3)

!DEBUG
!       if(j3==1)then
!        write(std_out,*)' back : j3==1, zw '
!        do i2=1,n2
!         do i1=1,n1dfft
!          write(std_out,'(2i4,2es16.6)' ) i1,i2,zw(1:2,i1+n1dfft*(i2-1))
!         end do
!        end do
!       end if
!ENDDEBUG

        inzee=1
        do i=1,ic2-1
        call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                    trig2,after2(i),now2(i),before2(i),1)
            inzee=3-inzee
        end do

        i=ic2
        call fftstp(lot,n1dfft,n2,nd1eff,nd2,zw(1,1,inzee),zr(1,j,1,j3,idat), &
                   trig2,after2(i),now2(i),before2(i),1)

2000        continue
! output: i1,i2,j3,(jp3)

!DEBUG
!       if(j3==1)then
!        write(std_out,*)' back : j3==1, zr '
!        do i2=1,n2
!         do i1=1,n1eff
!          write(std_out,'(2i4,2es16.6)' ) i1,i2,zr(1:2,i1,i2)
!         end do
!        end do
!       end if
!       if(icplex==1)stop
!ENDDEBUG

        end if
1212        continue

!DEBUG
!       write(std_out,*)' after j3 loop'
!ENDDEBUG

12345   continue
! $omp end do

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
        if (nproc_fft.gt.1)  then
          ABI_DEALLOCATE(zmpi1)
        end if
! $omp end parallel

!DEBUG
!       write(std_out,*)' back zr: exit '
!       do i3=1,nd3proc
!        do i2=1,n2
!         do i1=1,n1
!          write(std_out,'(3i4,2es16.6)')i1,i2,i3,zr(1:2,i1,i2,i3,1)
!         end do
!        end do
!       end do
!call leave_new("COLL")
!ENDDEBUG
        return
end subroutine back
