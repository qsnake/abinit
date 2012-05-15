#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

        subroutine forw_wf(icplexwf,mpi_enreg,ndat,n1,n2,n3,nd1,nd2,nd3proc,&
&        max1,max2,max3,m1,m2,m3,md1,md2proc,md3,nproc,iproc,paral_kgb,zr,zf)

 use m_profiling
! Does multiple 3-dim backward FFTs from real into Fourier space
! Adopt standard convention that isign=-1 for forward transform
!        CALCULATES THE DISCRETE FOURIERTRANSFORM ZF(I1,I3,I2)=
!        S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) ZR(j1,j2,j3)
!       in parallel using MPI/OpenMP.

!        INPUT:
!          ZR: input array
!               ZR(1,i1,i2,i3,idat)=real(R(i1,i2,i3,idat))
!               ZR(2,i1,i2,i3,idat)=imag(R(i1,i2,i3,idat))
!               i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!          n1,n2,n3: logical dimension of the transform. As transform lengths
!                    most products of the prime factors 2,3,5 are allowed.
!                    The detailed table with allowed transform lengths can
!                    be found in subroutine CTRIG
!          nd1,nd2,nd3: Dimension of ZR
!          nd3proc=((nd3-1)/nproc)+1  maximal number of big box 3rd dim slices for one proc
!        OUTPUT:
!          ZF: output array (note the switch of i2 and i3)
!               real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!               imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!          max1 is positive or zero ; m1 >=max1+1
!          i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!          then, if m1 > max1+1, one has min1=max1-m1+1 and
!          i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!          i2 and i3 have a similar definition of range
!          idat=1,ndat
!          md1,md2,md3: Dimension of ZF
!          md2proc=((md2-1)/nproc)+1  maximal number of small box 2nd dim slices for one proc
!          nproc: number of processors used as returned by MPI_COMM_SIZE
!          iproc: [0:nproc-1] number of processor as returned by MPI_COMM_RANK
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
   use defs_abitypes
   use defs_fftdata
#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'forw_wf'
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_52_fft_mpi_noabirule, except_this_one => forw_wf
!End of the abilint section

        implicit real(dp) (a-h,o-z)
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
! real space input        
        integer :: icplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc
        integer :: max1,max2,max3,m1,m2,m3,md1,md2proc,md3,nproc,iproc,paral_kgb
        REAL(DP), DIMENSION(2,nd1,nd2,nd3proc,ndat) :: zr
! Fourier space output
        REAL(DP), DIMENSION(2,md1,md3,md2proc,ndat) :: zf
!Local variables-------------------------------
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
! *************************************************************************
         integer :: spaceComm
         real(dp) :: tsec(2)
        type(MPI_type),intent(inout) :: mpi_enreg
         integer :: old_paral_level

!!      interface
!!        integer ( kind=4 ) function omp_get_num_threads ( )
!!        end function omp_get_num_threads
!!      end interface
!!      interface
!!        integer ( kind=4 ) function omp_get_thread_num ( )
!!        end function omp_get_thread_num
!!      end interface

        if(mpi_enreg%mode_para=='b') ioption=1
!       call timab(542,1,tsec)

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
!$omp shared(ndat,n1,n2,n3,nd1,nd2,nd3proc,md1,md2proc,md3,iproc,nproc,ncache,zr,zf,lock) &
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
!$omp end critical

        call ctrig(n2,trig2,after2,before2,now2,-1,ic2)
        call ctrig(n1,trig1,after1,before1,now1,-1,ic1)
        call ctrig(n3,trig3,after3,before3,now3,-1,ic3)

!DEBUG
!       write(std_out,'(a,3i4)' )'forw_wf, entre, i1,i2,i3,zr,n1,n2,n3',n1,n2,n3
!       write(std_out,'(a,3i4)' )'nd1,nd2,nd3proc',nd1,nd2,nd3proc
!       write(std_out,'(a,3i4)' )'m1,m2,m3',m1,m2,m3
!       write(std_out,'(a,3i4)' )'max1,max2,max3',max1,max2,max3
!       write(std_out,'(a,3i4)' )'md1,md2proc,md3',md1,md2proc,md3
!       write(std_out,'(a,3i4)' )'n1eff,m2eff,m1zt',n1eff,m2eff,m1zt
!       do i3=1,nd3proc
!        do i2=1,n2
!         do i1=1,n1
!          write(std_out,'(3i4,2es16.6)')i1,i2,i3,zr(1:2,i1,i2,i3)
!         end do
!        end do
!       end do
!       stop
!ENDDEBUG

!$omp do
        do 12345,idat=1,ndat

        do 1212,j3=1,nd3proc
        if (iproc*nd3proc+j3.le.n3) then
        Jp2st=1
        J2st=1

        if(icplexwf==1)then
         n1half=n1/2
         do i2=1,n2
          do i1=1,n1half
           zr(1,i1,i2,j3,idat)=zr(1,2*i1-1,i2,j3,idat)
           zr(2,i1,i2,j3,idat)=zr(1,2*i1  ,i2,j3,idat)
          end do
         end do
!        If odd
         if(n1half*2/=n1)then
          do i2=1,n2
           zr(1,n1eff,i2,j3,idat)=zr(1,n1,i2,j3,idat)
           zr(2,n1eff,i2,j3,idat)=0.0d0
          end do
         end if
        end if

! transform along y axis
! input: i1,i2,j3,(jp3)
        lot=ncache/(4*n2)

        do 2000,j=1,n1eff,lot
        ma=j
        mb=min(j+(lot-1),n1eff)
        n1dfft=mb-ma+1
        i=1
        call fftstp(nd1,n1dfft,nd2,lot,n2,zr(1,j,1,j3,idat),zw(1,1,1), &
                   trig2,after2(i),now2(i),before2(i),-1)

        inzee=1
        do i=2,ic2
        call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
                    trig2,after2(i),now2(i),before2(i),-1)
            inzee=3-inzee
        end do

!  input: i1,I2,j3,(jp3)
        if(icplexwf==2)then
         call unswitch_cent(n1dfft,max2,m2,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,j))
        else
         call unswitchreal_cent(n1dfft,max2,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,2*j-1))
        end if
! output: I2,i1,j3,(jp3)

2000        continue

!DEBUG
!       write(std_out,*)' j3=1, zt ='
!       do i1=1,m1zt
!        do i2=1,lzt
!         write(std_out,'(2i4,2es16.6)')i2,i1,zt(1:2,i2,i1)
!        end do
!       end do
!       stop
!ENDDEBUG

!DEBUG
!       write(std_out,*)' forwf_wf : after unswitch_cent, j3=',j3
!ENDDEBUG

! transform along x axis
! input: I2,i1,j3,(jp3)

        lot=ncache/(4*n1)

        do 1000,j=1,m2eff,lot
        ma=j
        mb=min(j+(lot-1),m2eff)
        n1dfft=mb-ma+1
        i=1
        call fftstp(lzt,n1dfft,m1zt,lot,n1,zt(1,j,1),zw(1,1,1), &
                   trig1,after1(i),now1(i),before1(i),-1)

        inzee=1
        do i=2,ic1
        call fftstp(lot,n1dfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
                    trig1,after1(i),now1(i),before1(i),-1)
            inzee=3-inzee
        end do
! output: I2,I1,j3,(jp3)

! input: J2,Jp2,I1,j3,(jp3)
        if (nproc.eq.1) then
        call unmpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot,max1,md1,m1,n1,&
&         md2proc,nd3proc,nproc,ioption,zw(1,1,inzee),zmpi2)
        else
         call unmpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot,max1,md1,m1,n1,&
&         md2proc,nd3proc,nproc,ioption,zw(1,1,inzee),zmpi1)
! output: I1,J2,j3,Jp2,(jp3)
        end if

1000        continue

        end if
1212        continue

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

!DEBUG
!       write(std_out,*)' forwf_wf : interproc data transp'
!ENDDEBUG


! Interprocessor data transposition
! intput: I1,J2,j3,Jp2,(jp3)
        if (nproc.gt.1) then
11      continue
!$omp   flush(lock)
        if (mod(lock,npr).ne.iam) goto 11
#if defined HAVE_MPI
  if(paral_kgb == 1) then
        old_paral_level=mpi_enreg%paral_level
        mpi_enreg%paral_level=3
        call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
        call timab(544,1,tsec)
        call MPI_ALLTOALL(zmpi1,2*md1*md2proc*nd3proc, &
                          MPI_double_precision, &
                          zmpi2,2*md1*md2proc*nd3proc, &
                          MPI_double_precision,spaceComm,ierr)
        call timab(544,2,tsec)
  endif
#endif
        lock=lock+1
!$omp   flush(lock)
        end if
! output: I1,J2,j3,jp3,(Jp2)

! transform along z axis
! input: I1,J2,i3,(Jp2)

        lot=ncache/(4*n3)

        do 3333,j2=1,md2proc
        if (iproc*md2proc+j2.le.m2eff) then

!DEBUG
!       write(std_out,*)' forwf_wf : before unscramble, j2,md2proc,iproc,m2=',j2,md2proc,iproc,m2
!ENDDEBUG

        do 3000,i1=1,m1,lot
        ma=i1
        mb=min(i1+(lot-1),m1)
        n1dfft=mb-ma+1

!   input: I1,J2,i3,(Jp2)
        call unscramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zmpi2,zw(1,1,1))
! output: I1,i3,J2,(Jp2)

        inzee=1
        do i=1,ic3
        call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
        trig3,after3(i),now3(i),before3(i),-1)
        inzee=3-inzee
        end do

        call unfill_cent(md1,md3,lot,n1dfft,max3,m3,n3,zw(1,1,inzee),zf(1,i1,1,j2,idat))
! output: I1,I3,J2,(Jp2)

3000        continue
        end if
3333    continue

        if(icplexwf==1)then
!        Complete missing values with complex conjugate
!        Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.

         do i3=1,m3
          i3inv=m3+2-i3
          if(i3==1)i3inv=1
          if(m2eff>1)then
           do i2=2,m2eff
            i2inv=m2+2-i2
            zf(1,1,i3inv,i2inv,idat)= zf(1,1,i3,i2,idat)
            zf(2,1,i3inv,i2inv,idat)=-zf(2,1,i3,i2,idat)
            do i1=2,m1
             i1inv=m1+2-i1
             zf(1,i1inv,i3inv,i2inv,idat)= zf(1,i1,i3,i2,idat)
             zf(2,i1inv,i3inv,i2inv,idat)=-zf(2,i1,i3,i2,idat)
            end do
           end do
          end if
         end do

        end if

12345   continue

!DEBUG
!      write(std_out,*)' zf ='
!      do i2=1,md2proc
!       do i3=1,md3
!        do i1=1,md1
!         write(std_out,'(3i4,2es16.6)')i1,i3,i2,zf(1:2,i1,i3,i2)
!        end do
!       end do
!      end do
!      stop
!ENDDEBUG
!DEBUG
!       write(std_out,'(a,3i4)' )'forw_wf exit ,zf n1,n2,n3',n1,n2,n3
!       write(std_out,'(a,3i4)' )'nd1,nd2,nd3proc',nd1,nd2,nd3proc
!       write(std_out,'(a,3i4)' )'m1,m2,m3',m1,m2,m3
!        write(std_out,'(a,3i4)' )'max1,max2,max3',max1,max2,max3
!        write(std_out,'(a,3i4)' )'md1,md2proc,md3',md1,md2proc,md3
!        write(std_out,'(a,3i4)' )'n1eff,m2eff,m1zt',n1eff,m2eff,m1zt
!        do i2=1,md2proc
!         do i3=1,m3
!          do i1=1,m1
!          if(i1 > max1+1) then ;ii1=i1+n1-m1;else; ii1=i1;end if !ii1 is the value in the full
!           if(i3 > max3+1) then ;ii3=i3+n3-m3;else; ii3=i3;end if !ii3 is the value in the full
!          ii2=(i2-1)*mpi_enreg%nproc_fft+mpi_enreg%me_fft+1
!          if(ii2 > max2+1) then;j2=i2-md2proc+n2/mpi_enreg%nproc_fft; else; j2=i2;end if
!          write(std_out,'(6i4,2es16.6)')i1,i3,i2,ii1,ii3,j2,zf(1:2,i1,i3,i2)
!          end do
!         end do
!        end do
!       stop
!ENDDEBUG

!$omp end do
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
!$omp end parallel

!       call timab(542,2,tsec)
        return
end subroutine forw_wf
