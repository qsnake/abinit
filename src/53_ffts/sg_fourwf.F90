!{\src2tex{textfont=tt}}
!!****f* ABINIT/sg_fourwf
!!
!! NAME
!! sg_fourwf
!!
!! FUNCTION
!! Carry out composite Fourier transforms between real and reciprocal (G) space.
!! Wavefunctions, contained in a sphere in reciprocal space,
!! can be FFT to real space. They can also be FFT from real space
!! to a sphere. Also, the density maybe accumulated, and a local
!! potential can be applied.
!!
!! The different options are :
!! - reciprocal to real space and output the result (option=0),
!! - reciprocal to real space and accumulate the density (option=1)
!! - reciprocal to real space, apply the local potential to the wavefunction
!!    in real space and produce the result in reciprocal space (option=2)
!! - real space to reciprocal space (option=3).
!!
!! Schedule of operations
!!(read first the description of the fftalg input variable in abinit_help)
!! - fftalgc=1 : use separate forward and backward transforms
!!     (7/12 savings in execution time);
!! - fftalgc=2 : in case of option=1 and option=2, use routines for composite operation
!!     even faster than 1x1
!!
!! Also for better speed, it uses no F90 construct, except the allocate command
!! and for zeroing arrays.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cplex= if 1 , denpot is real, if 2 , denpot is complex
!!    (cplex=2 only allowed for option=2, and istwf_k=1)
!!    not relevant if option=0 or option=3, so cplex=0 can be used to minimize memory
!! fftalgc=1 or 2 => simple or composite FFT applications
!! fofgin(2,npwin)=holds input wavefunction in G vector basis sphere.
!!                 (intent(in) but the routine sphere can modify it for another iflag)
!! gboundin(2*mgfft+8,2)=sphere boundary info for reciprocal to real space
!! gboundout(2*mgfft+8,2)=sphere boundary info for real to reciprocal space
!! istwf_k=option parameter that describes the storage of wfs
!! kg_kin(3,npwin)=reduced planewave coordinates, input
!! kg_kout(3,npwout)=reduced planewave coordinates, output
!! mgfft=maximum size of 1D FFTs
!! mpi_enreg=informations about MPI parallelization
!! n1,n2,n3=1D FFT sizes
!! npwin=number of elements in fofgin array (for option 0, 1 and 2)
!! npwout=number of elements in fofgout array (for option 2 and 3)
!! n4,n5,n6=dimensions of fofr.
!! option= if 0: do direct FFT
!!         if 1: do direct FFT, then sum the density
!!         if 2: do direct FFT, multiply by the potential, then do reverse FFT
!!         if 3: do reverse FFT only
!! paral_kgb=Flag related to the kpoint-band-fft parallelism
!! weight_r=weight to be used for the accumulation of the density in real space
!!         (needed only when option=1)
!! weight_i=weight to be used for the accumulation of the density in real space
!!         (needed only when option=1)
!!
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! for option==0, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                fofr(2,n4,n5,n6) contains the output Fourier Transform of fofgin;
!!                no use of denpot, fofgout and npwout.
!! for option==1, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                denpot(cplex*n4,n5,n6) contains the input density at input,
!!                and the updated density at output (accumulated);
!!                no use of fofgout and npwout.
!! for option==2, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                denpot(cplex*n4,n5,n6) contains the input local potential;
!!                fofgout(2,npwout) contains the output function;
!! for option==3, fofr(2,n4,n5,n6) contains the input real space wavefunction;
!!                fofgout(2,npwout) contains its output Fourier transform;
!!                no use of fofgin and npwin.
!!
!! PARENTS
!!      fourwf
!!
!! CHILDREN
!!      accrho,applypot,back_wf,forw_wf,sphere,sphere_fft1
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine sg_fourwf(cplex,denpot,fftalgc,fofgin,fofgout,fofr,&
&  gboundin,gboundout,istwf_k,kg_kin,kg_kout,mgfft,mpi_enreg,n1,n2,n3,&
&  npwin,npwout,n4,n5,n6,option,paral_kgb,weight_r,weight_i)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sg_fourwf'
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts, except_this_one => sg_fourwf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,fftalgc,istwf_k,mgfft,n1,n2,n3,n4,n5,n6,npwin
 integer,intent(in) :: npwout,option,paral_kgb
 real(dp),intent(in) :: weight_r,weight_i
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout)
 real(dp),intent(inout) :: denpot(cplex*n4,n5,n6),fofgin(2,npwin)
 real(dp),intent(inout) :: fofr(2,n4,n5,n6)
 real(dp),intent(out) :: fofgout(2,npwout)

!Local variables-------------------------------
!scalars
 integer :: cplexwf,i1,i2,i3,iflag,ig,iproc,isign,m1i,m1o,m2i,m2o,m3i
 integer :: m3o,max1i,max1o,max2i,max2i_plus,max2o,max3i,max3o,md1,md1i,md1o,md2,md2i,md2o,md2proc
 integer :: md3,md3i,md3o,min1i,min1o,min2i,min2i_moins,min2o,min3i,min3o,nd3proc,ndat,nfftot,nproc
! real(dp) :: tsec(2)
 real(dp) :: fim,fre,xnorm
 real(dp),allocatable :: weight_array(:)
 character(len=500) :: message
!arrays
 integer :: nd2proc,shiftg(3),symm(3,3)
 real(dp),allocatable :: work(:,:,:,:)

! *************************************************************************

!the maximum size of weight here is one
 ABI_ALLOCATE(weight_array,(1))

!DEBUG
!write(std_out,*)' sg_fourwf : enter, fftalgc,option= ',fftalgc,option
!if(option==2)fftalgc=1
!write(std_out,*)fofgin(1,1),fofgin(2,1)
!write(std_out,*)fofgin(1,101),fofgin(2,101)
!stop
!ENDDEBUG
!
!call timab(540,1,tsec)
 nfftot=n1*n2*n3

!ifdef DEBUG_MODE
 if (fftalgc<1 .and. fftalgc>2) then
   write(message,'(a,i0,3a)')&
&   '  The input algorithm number fftalgc=',fftalgc,' is not allowed. Must be 1 or 2',ch10,&
&   '  Action : change fftalgc in your input file.'
   MSG_PERS_ERROR(message)
 end if

 if (option<0 .or. option>3) then
   write(message,'(a,i0,3a)')&
&   '  The option number',option,' is not allowed.',ch10,&
&   '  Only option=0, 1, 2 or 3 are allowed presently.'
   MSG_PERS_ERROR(message)
 end if

 if (option==1 .and. cplex/=1) then
   write(message,'(a,i0,a)')&
&   '  With the option number 1, cplex must be 1 but it is cplex=',cplex,'.'
   MSG_PERS_ERROR(message)
 end if

 if ( ALL(cplex/=(/1,2/)) .and. ANY(option==(/1,2/)) ) then
   write(message,'(a,i0,a)')' When option is (1,2) cplex must be 1 or 2, but it is cplex=',cplex,'.'
   MSG_PERS_ERROR(message)
 end if
!endif

 shiftg(:)=0
 symm(:,:)=0
 symm(1,1)=1 ; symm(2,2)=1 ; symm(3,3)=1

 md1i=0 ; md2i=0 ; md3i=0
 md1o=0 ; md2o=0 ; md3o=0
 if(option/=3)then
   max1i=gboundin(2,1) ; min1i=gboundin(1,1)
   max2i=gboundin(4,1) ; min2i=gboundin(3,1)
   if(fftalgc==1 .or. fftalgc==2)then
     if(istwf_k==2 .or. istwf_k==4 .or. istwf_k==6 .or. istwf_k==8)then
       max1i=max(max1i,-min1i)
       min1i=-max1i
     else if(istwf_k==3 .or. istwf_k==5 .or. istwf_k==7 .or. istwf_k==9)then
       max1i=max(max1i,-min1i-1)
       min1i=-max1i-1
     end if
     if(istwf_k>=2 .and. istwf_k<=5)then
       max2i=max(max2i,-min2i)
       min2i=-max2i
     else if(istwf_k>=6 .and. istwf_k<=9)then
       max2i=max(max2i,-min2i-1)
       min2i=-max2i-1
     end if
   end if
   max3i=gboundin(4,2) ; min3i=gboundin(3,2)
   m1i=max1i-min1i+1 ; md1i=2*(m1i/2)+1
   m2i=max2i-min2i+1 ; md2i=2*(m2i/2)+1

   if (mpi_enreg%nproc_fft/=1) then
!    I increase max2i in order to have m2i divisible by nproc_fft
     min2i_moins=(((m2i-1)/mpi_enreg%nproc_fft+1)*mpi_enreg%nproc_fft-m2i)/2
     max2i_plus=((m2i-1)/mpi_enreg%nproc_fft+1)*mpi_enreg%nproc_fft-m2i-min2i_moins
!    max2i=max2i+((m2i-1)/mpi_enreg%nproc_fft+1)*mpi_enreg%nproc_fft-m2i
     max2i=max2i+max2i_plus
     min2i=min2i-min2i_moins
!    careful, to be checked and make sure the max and min are smaller than size of box
     m2i=max2i-min2i+1 ; md2i=2*(m2i/2)+1
   end if

   m3i=max3i-min3i+1 ; md3i=2*(m3i/2)+1
 end if
 if(option==2 .or. option==3)then
   max1o=gboundout(2,1) ; min1o=gboundout(1,1)
   max2o=gboundout(4,1) ; min2o=gboundout(3,1)
   if(fftalgc==1 .or. fftalgc==2)then
     if(istwf_k==2 .or. istwf_k==4 .or. istwf_k==6 .or. istwf_k==8)then
       max1o=max(max1o,-min1o)
       min1o=-max1o
     else if(istwf_k==3 .or. istwf_k==5 .or. istwf_k==7 .or. istwf_k==9)then
       max1o=max(max1o,-min1o-1)
       min1o=-max1o-1
     end if
     if(istwf_k>=2 .and. istwf_k<=5)then
       max2o=max(max2o,-min2o)
       min2o=-max2o
     else if(istwf_k>=6 .and. istwf_k<=9)then
       max2o=max(max2o,-min2o-1)
       min2o=-max2o-1
     end if
   end if
   max3o=gboundout(4,2) ; min3o=gboundout(3,2)
   m1o=max1o-min1o+1 ; md1o=2*(m1o/2)+1
   m2o=max2o-min2o+1 ; md2o=2*(m2o/2)+1
   m3o=max3o-min3o+1 ; md3o=2*(m3o/2)+1
 end if

 md1=max(md1i,md1o)
 md2=max(md2i,md2o)
 md3=max(md3i,md3o)
 md2proc=(m2i-1)/mpi_enreg%nproc_fft+1
 ABI_ALLOCATE(work,(2,md1,md3,md2proc))

 if(option/=3)then
!  Insert fofgin into the small fft box (array work) :
   iflag=2
!  Note the switch of md3 and md2, as they are only
!  needed to dimension work2 inside "sphere"
   if(mpi_enreg%paral_compil_fft==1) then
     if(istwf_k/=1 )then
       write(message,'(a,i0,3a)')&
&       '  The value of istwf_k',istwf_k,' is not allowed.',ch10,&
&       '  Only istwf_k=1 is allowed in FFT parallel mode'
       MSG_BUG(message)
     end if
     nd2proc=((n2-1)/mpi_enreg%nproc_fft) +1
     call sphere_fft1(fofgin,1,npwin,work,m1i,m2i,m3i,md1,md3,md2proc,&
&     kg_kin,mpi_enreg,nd2proc)
   else
     call sphere(fofgin,1,npwin,work,m1i,m2i,m3i,md1,md3,md2,&
&     kg_kin,istwf_k,iflag,mpi_enreg,shiftg,symm,one)
   end if
!  call leave_new("COLL")
!  write(std_out,*) 'i1,i2,i3,work(:,i1,i2,i3)',md1,md3,md2
!  do i3=1,md2proc
!  do i1=1,m1i
!  do i2=1,m3i
!  write(std_out,'(3i4,2e24.10)')i1,i2,i3,work(:,i1,i2,i3)
!  end do;enddo;enddo

 end if
!call leave_new("COLL")
!if(option==2)then
!j1=0
!do i2=1,m2i
!do i3=1,m3i
!do i1=1,m1i
!j1=j1+1
!if(mod(j1,31)==1)write(std_out,'(3i4,2es16.6)' )i1,i2,i3,work(1:2,i1,i3,i2)
!end do
!end do
!end do
!do i3=1,n3
!do i2=1,n2
!do i1=1,n1
!denpot(i1,i2,i3)=1.0d0
!end do
!end do
!end do
!end if

!------------------------------------------------------------------
!Treat non-composite operations

 if(  option==0                                 .or. &
& ((option==1.or.option==2) .and. fftalgc==1) .or. &
& option==3                                         )then

!  DEBUG
!  write(std_out,*)' sg_fourwf : max1i,max2i,max3i=',max1i,max2i,max3i
!  write(std_out,*)' sg_fourwf : min1i,min2i,min3i=',min1i,min2i,min3i
!  write(std_out,*)' sg_fourwf : m1i,m2i,m3i=',m1i,m2i,m3i
!  do i3=1,m3i
!  do i2=1,m2i
!  do i1=1,m1i
!  write(std_out,'(3i4,2es16.6)')i1,i2,i3,work(1:2,i1,i3,i2)
!  end do
!  end do
!  end do
!  stop
!  ENDDEBUG


!  Fourier transform work to fofr (reciprocal to real space)
   if(option/=3)then
     isign=1
     ndat=1 ; nproc=mpi_enreg%nproc_fft ; iproc=mpi_enreg%me_fft
     cplexwf=2
     if(istwf_k==2)cplexwf=1
     call back_wf(cplexwf,mpi_enreg,ndat,n1,n2,n3,n4,n5,(n6-1)/mpi_enreg%nproc_fft+1,&
&     max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,nproc,iproc,paral_kgb,work,fofr)
   end if

!  DEBUG
!  write(std_out,*)' sg_fourwf : after back_wf,i1,i2,i3,fofr'
!  do i3=1,(n3-1)/mpi_enreg%nproc_fft+1
!  do i2=1,n2
!  do i1=1,n1
!  write(std_out,'(3i4,2es16.6)')i1,i2,i3,fofr(1:2,i1,i2,i3)
!  end do
!  end do
!  end do
!  call leave_new("COLL")
!  stop
!  ENDDEBUG
   nd3proc=(n3-1)/mpi_enreg%nproc_fft+1

   if(option==1)then  ! Accumulate density
     do i3=1,nd3proc
       do i2=1,n2
         do i1=1,n1

!          denpot(i1,i2,mpi_enreg%me_fft*nd3proc+i3)=denpot(i1,i2,mpi_enreg%me_fft*nd3proc+i3)+&
!          &            weight_r*(fofr(1,i1,i2,i3)**2+fofr(2,i1,i2,i3)**2)
           denpot(i1,i2,mpi_enreg%me_fft*nd3proc+i3)=denpot(i1,i2,mpi_enreg%me_fft*nd3proc+i3)+&
&           (weight_r*fofr(1,i1,i2,i3)**2+ weight_i*fofr(2,i1,i2,i3)**2)

         end do
       end do
     end do
   end if ! option==1

!  Apply local potential
   if(option==2)then
     if(cplex==1)then
       do i3=1,nd3proc
         do i2=1,n2
           do i1=1,n1
             fofr(1,i1,i2,i3)=denpot(i1,i2,mpi_enreg%me_fft*nd3proc+i3)*fofr(1,i1,i2,i3)
             fofr(2,i1,i2,i3)=denpot(i1,i2,mpi_enreg%me_fft*nd3proc+i3)*fofr(2,i1,i2,i3)
           end do
         end do
       end do
     else if(cplex==2)then
!      MSG_ERROR('cplex==2 not treated in sg_fourwf') ! MT 20110730: why this stop ? seems to work
       do i3=1,(n3-1)/mpi_enreg%nproc_fft+1
         do i2=1,n2
           do i1=1,n1
             fre=fofr(1,i1,i2,i3)
             fim=fofr(2,i1,i2,i3)
             fofr(1,i1,i2,i3)=denpot(2*i1-1,i2,i3)*fre -denpot(2*i1  ,i2,i3)*fim
             fofr(2,i1,i2,i3)=denpot(2*i1-1,i2,i3)*fim +denpot(2*i1  ,i2,i3)*fre
           end do
         end do
       end do
     end if ! cplex
   end if ! option==2

!  Fourier transform fofr to work (real to reciprocal space)
   if(option==2 .or. option==3)then
     isign=-1
     iflag=2
     ndat=1 ; nproc=mpi_enreg%nproc_fft ; iproc=mpi_enreg%me_fft
     cplexwf=2
     if(istwf_k==2)cplexwf=1
     max1o=max1i;max2o=max2i;max3o=max3i;m1o=m1i;m2o=m2i;m3o=m3i
!    DEBUG
!    write(std_out,*) 'before hitting forw_wf, i1,i2,i3,fofr'
!    do i3=1,(n3-1)/mpi_enreg%nproc_fft+1
!    do i2=1,n2
!    do i1=1,n1
!    write(std_out,'(3i4,2es16.6)') i1,i2,i3,fofr(:,i1,i2,i3)
!    end do
!    end do
!    end do

!    ENDDEBUG

     call forw_wf(cplexwf,mpi_enreg,ndat,n1,n2,n3,n4,n5,(n6-1)/mpi_enreg%nproc_fft+1,&
&     max1o,max2o,max3o,m1o,m2o,m3o,md1,md2proc,md3,nproc,iproc,paral_kgb,fofr,work)

!    DEBUG
!    do i3=1,md2proc
!    do i1=1,m1i
!    do i2=1,m3i
!    write(std_out,'(3i4,2e24.10)')i1,i2,i3,work(:,i1,i2,i3)
!    end do;enddo;enddo
!    call leave_new("COLL")
!    ENDDEBUG

   end if

!  ------------------------------------------------------------------
!  Treat composite operations

 else if(fftalgc==2 .and. (option==1 .or. option==2) )then

   if(option==1)then

     ndat=1 ; nproc=1 ; iproc=0
     cplexwf=2
     if(istwf_k==2)cplexwf=1
     weight_array(1)=weight_r
     call accrho(cplexwf,ndat,n1,n2,n3,n4,n5,n6,&
&     max1i,max2i,max3i,m1i,m2i,m3i,md1,md2,md3,mpi_enreg%comm_fft,nproc,iproc,&
&     paral_kgb,work,denpot,weight_array(1:ndat))
   else if(option==2)then

     ndat=1 ; nproc=1 ; iproc=0
     cplexwf=2
     if(istwf_k==2)cplexwf=1
     call applypot(cplexwf,cplex,ndat,n1,n2,n3,n4,n5,n6,&
&     max1i,max2i,max3i,m1i,m2i,m3i,md1,md2,md3,&
&     max1o,max2o,max3o,m1o,m2o,m3o,mpi_enreg%comm_fft,nproc,iproc,&
&     paral_kgb,denpot,work)

   end if

 end if

!End of composite operations
!-----------------------------------------------------------------

!if(option==2)then
!j1=0
!do i2=1,m2o
!do i3=1,m3o
!do i1=1,m1o
!j1=j1+1
!if(mod(j1,31)==1)write(std_out,'(3i4,2es16.6)' )i1,i2,i3,work(1:2,i1,i3,i2)/nfftot
!end do
!end do
!end do
!stop
!end if

 if(option==2 .or. option==3)then

   iflag=-2
   xnorm=1.d0/dble(nfftot)
   if(mpi_enreg%mode_para=='b') then
     do ig=1,npwout
       i1=kg_kout(1,ig); if(i1<0)i1=i1+m1o ; i1=i1+1
       i2=kg_kout(2,ig); if(i2<0)i2=i2+m2o ; i2=i2+1
       i3=kg_kout(3,ig); if(i3<0)i3=i3+m3o ; i3=i3+1
       fofgout(1,ig)=work(1,i1,i3,(i2-1)/mpi_enreg%nproc_fft +1)*xnorm
       fofgout(2,ig)=work(2,i1,i3,(i2-1)/mpi_enreg%nproc_fft +1)*xnorm
!      write(std_out,'(5i3,2e24.12)')i1,i2,(i2-1)/mpi_enreg%nproc_fft +1,i3,ig,fofgout(:,ig)
     end do
   else
     call sphere(fofgout,1,npwout,work,m1o,m2o,m3o,md1,md3,md2,kg_kout,istwf_k,iflag,&
&     mpi_enreg,shiftg,symm,xnorm)
   end if
!  xnorm=1.d0/dble(nfftot)
!  do ig=1,npwout
!  i1=kg_kout(1,ig); if(i1<0)i1=i1+m1o ; i1=i1+1
!  i2=kg_kout(2,ig); if(i2<0)i2=i2+m2o ; i2=i2+1
!  i3=kg_kout(3,ig); if(i3<0)i3=i3+m3o ; i3=i3+1
!  fofgout(1,ig)=work(1,i1,i3,i2)*xnorm
!  fofgout(2,ig)=work(2,i1,i3,i2)*xnorm
!  end do

 end if ! if option==2 or 3

 ABI_DEALLOCATE(work)

 ABI_DEALLOCATE(weight_array)
!call timab(540,2,tsec)

end subroutine sg_fourwf
!!***
