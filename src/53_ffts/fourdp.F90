!{\src2tex{textfont=tt}}
!!****f* ABINIT/fourdp
!! NAME
!! fourdp
!!
!! FUNCTION
!! Conduct Fourier transform of REAL or COMPLEX function f(r)=fofr defined on
!! fft grid in real space, to create complex f(G)=fofg defined on full fft grid 
!! in reciprocal space, in full storage mode, or the reverse operation.
!! For the reverse operation, the final data is divided by nfftot.
!! REAL case when cplex=1, COMPLEX case when cplex=2
!! Usually used for density and potentials.
!!
!! There are two different possibilities :
!!  fftalgb=0 means using the complex-to-complex FFT routine,
!!   irrespective of the value of cplex
!!  fftalgb=1 means using a real-to-complex FFT or a complex-to-complex FFT,
!!   depending on the value of cplex.
!!  The only real-to-complex FFT available is from SGoedecker library.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! isign=sign of Fourier transform exponent: current convention uses
!!  +1 for transforming from G to r 
!!  -1 for transforming from r to G.
!! mpi_enreg=information about MPI parallelization
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! paral_kgb=Flag related to the kpoint-band-fft parallelism
!! tim_fourdp=timing code of the calling routine (can be set to 0 if not attributed)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! fofg(2,nfft)=f(G), complex.
!! fofr(cplex*nfft)=input function f(r) (real or complex)
!!
!! TODO
!! work2 should be a pointer, and allocated inside ccfft, if really needed ...
!!
!! PARENTS
!!      atm2fft,atm2fft3,bethe_salpeter,calc_smeared_density,dens_in_sph
!!      dieltcel,difvxc,dyfro3,energy,fftw3_fourdp,fftwfn,filterpot,forces
!!      fourdp_6d,fresidrsp,ftfvw1,ftfvw2,green_kernel,gstate,hartre,hartre1
!!      hartrestr,initro,jellium,jvec_to_B,kxc_alda,kxc_eok,ladielmt,laplacian
!!      loop3dte,loper3,m_electronpositron,m_fft_prof,m_hidecudarec,m_ppmodel
!!      m_screening,make_efg_el,mklocl_realspace,mklocl_recipspace,moddiel
!!      moddiel_csrb,mrgscr,newrho,newvtr,newvtr3,nonlinear,nres2vres,odamix
!!      pawfrnhat_recipspace,pawmknhat,pawmknhat_psipsi,pawmkrho,prcref
!!      prcref_PMA,prctfvw1,prctfvw2,prctfw3,rdm,recursion,recursion_nl,redgr
!!      respfn,scfcv,scfcv3,screening,setup_positron,sigma,stress,symrhg,tddft
!!      transgrid,vloca3,vlocalstr,xc_kernel,xc_kernel_ADA,xcden,xcpot
!!
!! CHILDREN
!!      back,ccfft,fftw3_fourdp,forw,sg_fft_rc,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fourdp(cplex,fofg,fofr,isign,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

 use m_sgfft,       only : sg_fft_rc

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fourdp'
 use interfaces_18_timing
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts, except_this_one => fourdp
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,isign,nfft,paral_kgb,tim_fourdp
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(inout) :: fofg(2,nfft),fofr(cplex*nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: mfac=11
 integer :: fftalg,fftalga,fftalgb,fftcache,i1,i2,i3,index
 integer :: inplace,ir,max_index,me_fft,n1,n1half1,n1halfm,n2,n2half1,n3,n4
 integer :: n4half1,n5,n5half1,n6,nd2proc,nd3proc,ndat,normalized,nproc_fft
 real(dp) :: ris,xnorm
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: work1(:,:,:,:),work2(:,:,:,:)
 real(dp),allocatable :: workf(:,:,:,:),workr(:,:,:,:)

! *************************************************************************

!DBG_ENTER("COLL")
!write(std_out,*)' fourdp : enter, fftalg,isign=',ngfft(7),isign

!Keep track of timing
 call timab(260+tim_fourdp,1,tsec)

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)
 me_fft=ngfft(11); nproc_fft=ngfft(10)

 fftcache=ngfft(8)
 fftalg  =ngfft(7)
 fftalga =fftalg/100
 fftalgb =mod(fftalg,100)/10

 ris=dble(isign)
 xnorm=1.0d0/dble(n1*n2*n3)
!write(std_out,*)' fourdp :me_fft',me_fft,'nproc_fft',nproc_fft,'nfft',nfft

 if (fftalgb/=0 .and. fftalgb/=1) then
   write(message, '(a,i4,a,a,a,a,a)' )&
&   '  The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&   '  The second digit (fftalg(B)) must be 0 or 1.',ch10,&
&   '  Action : change fftalg in your input file.'
   MSG_PERS_BUG(message)
 end if

 if (fftalgb==1 .and. ALL(fftalga/=(/1,3,4/)) )then
   write(message,'(a,i4,5a)')&
&   '  The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&   '  When fftalg(B) is 1, the allowed values for fftalg(A) are 1 and 4.',ch10,&
&   '  Action : change fftalg in your input file.'
   MSG_PERS_BUG(message)
 end if

 if (n4<n1.or.n5<n2.or.n6<n3) then
   write(message,'(a,3i8,a,3i8)')'  Each of n4,n5,n6=',n4,n5,n6,'must be >= n1, n2, n3 =',n1,n2,n3
   MSG_PERS_BUG(message)
 end if

 if (fftalg==312) then ! call FFTW3 version.
   call fftw3_fourdp(cplex,n1,n2,n3,isign,fofg,fofr)
   call timab(260+tim_fourdp,2,tsec)
   RETURN
 end if

!---------------------------------------------------------
!Here, deal  with the new SG FFT, complex-to-complex case
 if( fftalga==4 .and. (fftalgb==0 .or. cplex==2) )then

   nd2proc=((n2-1)/nproc_fft) +1
   nd3proc=((n6-1)/nproc_fft) +1
   ABI_ALLOCATE(workr,(2,n4,n5,nd3proc))
   ABI_ALLOCATE(workf,(2,n4,n6,nd2proc))
   max_index=0
   if(isign==1)then
     do i3=1,n3
       do i2=1,n2
         if (((i2-1)/(n2/nproc_fft))==me_fft) then
           index=n1*(i2-me_fft*n2/nproc_fft-1+(n2/nproc_fft)*(i3-1))
           if (index > max_index) max_index=index
           do i1=1,n1
             workf(1,i1,i3,i2-me_fft*n2/nproc_fft)=fofg(1,i1+index)
             workf(2,i1,i3,i2-me_fft*n2/nproc_fft)=fofg(2,i1+index)
           end do
         end if
       end do
     end do
     ndat=1  
     call back(2,mpi_enreg,ndat,n1,n2,n3,n4,n5,n6,n4,nd2proc,nd3proc,2,paral_kgb,workf,workr)
     if(cplex==1)then
       do i3=1,n3
         if (((i3-1)/(n3/nproc_fft))==me_fft) then
           do i2=1,n2
             index=n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1))
             do i1=1,n1
               fofr(i1+index)=workr(1,i1,i2,i3-me_fft*n3/nproc_fft)
             end do
           end do
         end if
       end do
     else if(cplex==2)then
       do i3=1,n3
         if (((i3-1)/(n3/nproc_fft))==me_fft) then
           do i2=1,n2
             index=2*n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1))
             do i1=1,n1
               fofr(2*i1-1+index)=workr(1,i1,i2,i3-me_fft*n3/nproc_fft)
               fofr(2*i1  +index)=workr(2,i1,i2,i3-me_fft*n3/nproc_fft)
             end do
           end do
         end if
       end do
     end if

   else if(isign==-1)then

     if(cplex==1)then
       do i3=1,n3
         if (((i3-1)/(n3/nproc_fft))==me_fft) then
           do i2=1,n2
             index=n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1))
!            index=n1*(i2-1+n2*(i3-1))
             do i1=1,n1
               workr(1,i1,i2,i3-me_fft*n3/nproc_fft)=fofr(i1+index)
               workr(2,i1,i2,i3-me_fft*n3/nproc_fft)=zero
             end do
           end do
         end if
       end do
     else if(cplex==2)then
       do i3=1,n3
         if (((i3-1)/(n3/nproc_fft))==me_fft) then
           do i2=1,n2
             index=2*n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1))
             do i1=1,n1
               workr(1,i1,i2,i3-me_fft*n3/nproc_fft)=fofr(2*i1-1+index)
               workr(2,i1,i2,i3-me_fft*n3/nproc_fft)=fofr(2*i1  +index)
             end do
           end do
         end if
       end do
     end if

     ndat=1
     call forw(2,mpi_enreg,ndat,n1,n2,n3,n4,n5,n6,n4,nd2proc,nd3proc,2,paral_kgb,workr,workf)

!    Transfer fft output to the original fft box
     do i2=1,n2
       if (((i2-1)/(n2/nproc_fft))==me_fft) then
         do i3=1,n3
           index=n1*(i2-me_fft*n2/nproc_fft-1+(n2/nproc_fft)*(i3-1))
           do i1=1,n1
             fofg(1,i1+index)=workf(1,i1,i3,i2-me_fft*n2/nproc_fft)*xnorm
             fofg(2,i1+index)=workf(2,i1,i3,i2-me_fft*n2/nproc_fft)*xnorm
           end do
         end do
       end if
     end do
     do ir=1,40
     end do

   end if ! isign

   ABI_DEALLOCATE(workr)
   ABI_DEALLOCATE(workf)

 end if

!---------------------------------------------------------
!Here, deal with the new SG FFT, with real-to-complex
 if(fftalga==4 .and. fftalgb==1 .and. cplex==1)then

   n1half1=n1/2+1 ; n1halfm=(n1+1)/2
   n2half1=n2/2+1
!  n4half1 or n5half1 are the odd integers >= n1half1 or n2half1
   n4half1=(n1half1/2)*2+1
   n5half1=(n2half1/2)*2+1
   ABI_ALLOCATE(workr,(2,n4half1,n5,n6))
   ABI_ALLOCATE(workf,(2,n4,n6,n5half1))
   if(isign==1)then
     do i3=1,n3
       do i2=1,n2half1
         index=n1*(i2-1+n2*(i3-1))
         do i1=1,n1
           workf(1,i1,i3,i2)=fofg(1,i1+index)
           workf(2,i1,i3,i2)=fofg(2,i1+index)
         end do
       end do
     end do

!    DEBUG
!    write(std_out,*)' fourdp : before back '
!    do i2=1,n2half1
!    do i3=1,n3
!    do i1=1,n1
!    index=i1+n1*(i2-1+n2*(i3-1))
!    write(std_out,'(3i4,2es16.6)' ) i1,i3,i2,workf(1:2,i1,i3,i2)
!    end do
!    end do
!    end do
!    ENDDEBUG

     ndat=1;
     nd2proc=((n5-1)/nproc_fft) +1
     nd3proc=((n6-1)/nproc_fft) +1
!    modifier l appel ? n5half1 et n6 ?
     call back(cplex,mpi_enreg,ndat,n1,n2,n3,n4,n5,n6,n4half1,n5half1,n6,2,paral_kgb,workf,workr)

!    DEBUG
!    write(std_out,*)' fourdp : after back '
!    do i3=1,n3
!    do i2=1,n2
!    do i1=1,n1half1
!    index=i1+n1*(i2-1+n2*(i3-1))
!    write(std_out,'(3i4,2es16.6)' ) i1,i2,i3,workr(1:2,i1,i2,i3)
!    end do
!    end do
!    end do
!    stop
!    ENDDEBUG

     do i3=1,n3
       do i2=1,n2
         index=n1*(i2-1+n2*(i3-1))
         do i1=1,n1half1-1
!          copy data
           fofr(2*i1-1+index)=workr(1,i1,i2,i3)
           fofr(2*i1  +index)=workr(2,i1,i2,i3)
         end do
!        If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           fofr(n1+index)=workr(1,n1half1,i2,i3)
         end if
       end do
     end do

!    DEBUG
!    write(std_out,*)' fourdp : stop, isign=1'
!    do i3=1,n3
!    do i2=1,n2
!    do i1=1,n1
!    index=i1+n1*(i2-1+n2*(i3-1))
!    write(std_out,'(3i4,2es16.6)' ) i1,i2,i3,fofr(index)
!    end do
!    end do
!    end do
!    stop
!    ENDDEBUG

   else if(isign==-1)then
     do i3=1,n3
       do i2=1,n2
         index=n1*(i2-1+n2*(i3-1))
         do i1=1,n1half1-1
           workr(1,i1,i2,i3)=fofr(2*i1-1+index)
           workr(2,i1,i2,i3)=fofr(2*i1  +index)
         end do
!        If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           workr(1,n1half1,i2,i3)=fofr(n1+index)
           workr(2,n1half1,i2,i3)=0.0d0
         end if
       end do
     end do
     ndat=1
     call forw(cplex,mpi_enreg,ndat,n1,n2,n3,n4,n5,n6,n4half1,n5half1,n6,2,paral_kgb,workr,workf)
!    Transfer fft output to the original fft box
     do i3=1,n3
       do i2=1,n2half1
         index=n1*(i2-1+n2*(i3-1))
         do i1=1,n1
           fofg(1,i1+index)=workf(1,i1,i3,i2)*xnorm
           fofg(2,i1+index)=workf(2,i1,i3,i2)*xnorm
         end do
       end do
!      Complete missing values with complex conjugate
!      Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.
       if(n2half1>2)then
         do i2=2,n2+1-n2half1
           index=n1*((n2+2-i2)-1)
           if(i3/=1)index=index+n1*n2*((n3+2-i3)-1)
           fofg(1,1+index)= workf(1,1,i3,i2)*xnorm
           fofg(2,1+index)=-workf(2,1,i3,i2)*xnorm
           do i1=2,n1
             fofg(1,n1+2-i1+index)= workf(1,i1,i3,i2)*xnorm
             fofg(2,n1+2-i1+index)=-workf(2,i1,i3,i2)*xnorm
           end do
         end do
       end if
     end do
!    DEBUG
!    write(std_out,*)' fourdp : stop, isign=-1'
!    do i3=1,n3
!    do i2=1,n2
!    do i1=1,n1
!    index=i1+n1*(i2-1+n2*(i3-1))
!    write(std_out,'(3i4,2es16.6)' ) i1,i2,i3,fofg(1:2,index)
!    end do
!    end do
!    end do
!    stop
!    ENDDEBUG

   end if ! isign
   ABI_DEALLOCATE(workr)
   ABI_DEALLOCATE(workf)

 end if

!---------------------------------------------------------
!Here, one calls the complex-to-complex FFT subroutine
 if( (fftalgb==0 .or. cplex==2) .and. fftalga/=4 )then

   ABI_ALLOCATE(work1,(2,n4,n5,n6))
   ABI_ALLOCATE(work2,(2,n4,n5,n6))

   if(isign==1)then

!    Transfer fofg to the expanded fft box
!    $OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofg,n1,n2,n3,work1)
     do i3=1,n3
       do i2=1,n2
         index=n1*(i2-1+n2*(i3-1))
         do i1=1,n1
           work1(1,i1,i2,i3)=fofg(1,i1+index)
           work1(2,i1,i2,i3)=fofg(2,i1+index)
         end do
       end do
     end do
!    $OMP END PARALLEL DO

!    Actual 3D FFT
     call ccfft(fftalga,fftcache,inplace,isign,mpi_enreg,normalized,&
&     n1,n2,n3,n4,n5,n6,1,2,paral_kgb,work1,work2)

!    Take data from expanded box and put it in the original box.
     if(cplex==1)then
!      REAL case

       if(inplace==0)then
!        $OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofr,n1,n2,n3,work2)
         do i3=1,n3
           do i2=1,n2
             index=n1*(i2-1+n2*(i3-1))
             do i1=1,n1
               fofr(i1+index)=work2(1,i1,i2,i3)
             end do
           end do
         end do
!        $OMP END PARALLEL DO
       else if(inplace==1)then
!        $OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofr,n1,n2,n3,work1)
         do i3=1,n3
           do i2=1,n2
             index=n1*(i2-1+n2*(i3-1))
             do i1=1,n1
               fofr(i1+index)=work1(1,i1,i2,i3)
             end do
           end do
         end do
!        $OMP END PARALLEL DO
       end if

     else
!      COMPLEX case

       if(inplace==0)then
!        $OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofr,n1,n2,n3,work2)
         do i3=1,n3
           do i2=1,n2
             index=2*n1*(i2-1+n2*(i3-1))
             do i1=1,n1
               fofr(2*i1-1+index)=work2(1,i1,i2,i3)
               fofr(2*i1  +index)=work2(2,i1,i2,i3)
             end do
           end do
         end do
!        $OMP END PARALLEL DO
       else if(inplace==1)then
!        $OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofr,n1,n2,n3,work1)
         do i3=1,n3
           do i2=1,n2
             index=2*n1*(i2-1+n2*(i3-1))
             do i1=1,n1
               fofr(2*i1-1+index)=work1(1,i1,i2,i3)
               fofr(2*i1  +index)=work1(2,i1,i2,i3)
             end do
           end do
         end do
!        $OMP END PARALLEL DO
       end if

     end if

   else if(isign==-1)then

!    Insert fofr into the augmented fft box
     if(cplex==1)then ! REAL case
!      $OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofr,n1,n2,n3,work1)
       do i3=1,n3
         do i2=1,n2
           index=n1*(i2-1+n2*(i3-1))
           do i1=1,n1
!            copy data
             work1(1,i1,i2,i3)=fofr(i1+index)
             work1(2,i1,i2,i3)=0.0d0
           end do
         end do
       end do
!      $OMP END PARALLEL DO
     else
!      COMPLEX case
!      $OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofr,n1,n2,n3,work1)
       do i3=1,n3
         do i2=1,n2
           index=2*n1*(i2-1+n2*(i3-1))
           do i1=1,n1
!            copy data
             work1(1,i1,i2,i3)=fofr(2*i1-1+index)
             work1(2,i1,i2,i3)=fofr(2*i1  +index)
           end do
         end do
       end do
!      $OMP END PARALLEL DO
     end if ! cplex

!    Actual 3D FFT
     call ccfft(fftalga,fftcache,inplace,isign,mpi_enreg,normalized,&
&     n1,n2,n3,n4,n5,n6,1,2,paral_kgb,work1,work2)

!    Transfer fft output to the original fft box
     if(normalized==0)then
       if(inplace==0)then
!        $OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofg,n1,n2,n3,work2,xnorm)
         do i3=1,n3
           do i2=1,n2
             index=n1*(i2-1+n2*(i3-1))
             do i1=1,n1
               fofg(1,i1+index)=work2(1,i1,i2,i3)*xnorm
               fofg(2,i1+index)=work2(2,i1,i2,i3)*xnorm
             end do
           end do
         end do
!        $OMP END PARALLEL DO
       else if(inplace==1)then
!        $OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofg,n1,n2,n3,work1,xnorm)
         do i3=1,n3
           do i2=1,n2
             index=n1*(i2-1+n2*(i3-1))
             do i1=1,n1
               fofg(1,i1+index)=work1(1,i1,i2,i3)*xnorm
               fofg(2,i1+index)=work1(2,i1,i2,i3)*xnorm
             end do
           end do
         end do
!        $OMP END PARALLEL DO
       end if
     else if(normalized==1)then
       if(inplace==0)then
!        $OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofg,n1,n2,n3,work2)
         do i3=1,n3
           do i2=1,n2
             index=n1*(i2-1+n2*(i3-1))
             do i1=1,n1
               fofg(1,i1+index)=work2(1,i1,i2,i3)
               fofg(2,i1+index)=work2(2,i1,i2,i3)
             end do
           end do
         end do
!        $OMP END PARALLEL DO
       else if(inplace==1)then
!        $OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofg,n1,n2,n3,work1)
         do i3=1,n3
           do i2=1,n2
             index=n1*(i2-1+n2*(i3-1))
             do i1=1,n1
               fofg(1,i1+index)=work1(1,i1,i2,i3)
               fofg(2,i1+index)=work1(2,i1,i2,i3)
             end do
           end do
         end do
!        $OMP END PARALLEL DO
       end if
     end if ! normalized==0 or 1

!    Endif choice of isign
   end if

   ABI_DEALLOCATE(work1)
   ABI_DEALLOCATE(work2)

 end if ! End simple algorithm

!---------------------------------------------------------
!Here sophisticated algorithm based on S. Goedecker routines, only for the REAL case.
!Take advantage of the fact that fofr is real, and that fofg has corresponding symmetry properties.

 if( (fftalgb==1 .and. cplex==1) .and. fftalga/=4 )then
   call sg_fft_rc(cplex,fofg,fofr,isign,nfft,ngfft)
 end if

 call timab(260+tim_fourdp,2,tsec)

!DBG_EXIT("COLL")

end subroutine fourdp
!!***
