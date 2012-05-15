!{\src2tex{textfont=tt}}
!!****f* ABINIT/fourwf
!! NAME
!! fourwf
!!
!! FUNCTION
!! Carry out composite Fourier transforms between real and reciprocal (G) space.
!! Wavefunctions, contained in a sphere in reciprocal space,
!! can be FFT to real space. They can also be FFT from real space
!! to a sphere. Also, the density maybe accumulated, and a local
!! potential can be applied.
!!
!! The different options are :
!! - option=0 --> reciprocal to real space and output the result.
!! - option=1 --> reciprocal to real space and accumulate the density.
!! - option=2 --> reciprocal to real space, apply the local potential to the wavefunction
!!                in real space and produce the result in reciprocal space.
!! - option=3 --> real space to reciprocal space.
!!                NOTE that in this case, fftalg=1x1 MUST be used. This may be changed in the future.
!!
!! The different sections of this routine corresponds to different
!! algorithms, used independently of each others :
!!(read first the description of the fftalg input variable in abinit_help)
!! - fftalg=xx0 : use simple complex-to-complex routines, without zero padding
!!     (rather simple, so can be used to understand how fourwf.f works);
!! - fftalg=1x1 : use S Goedecker routines, with zero padding
!!     (7/12 savings in execution time);
!! - fftalg=1x2 : call even more sophisticated coding also based on S Goedecker routines
!!
!! This routine contains many parts that differ only
!! by small details, in order to treat each case with the better speed.
!! Also for better speed, it uses no F90 construct, except the allocate command
!! and for zeroing arrays.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, FF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cplex= if 1 , denpot is real, if 2 , denpot is complex
!!    (cplex=2 only allowed for option=2, and istwf_k=1)
!!    not relevant if option=0 or option=3, so cplex=0 can be used to minimize memory
!! fofgin(2,npwin)=holds input wavefunction in G vector basis sphere.
!!                 (intent(in) but the routine sphere can modify it for another iflag)
!! gboundin(2*mgfft+8,2)=sphere boundary info for reciprocal to real space
!! gboundout(2*mgfft+8,2)=sphere boundary info for real to reciprocal space
!! istwf_k=option parameter that describes the storage of wfs
!! kg_kin(3,npwin)=reduced planewave coordinates, input
!! kg_kout(3,npwout)=reduced planewave coordinates, output
!! mgfft=maximum size of 1D FFTs
!! mpi_enreg=information about MPI parallelization
!! ndat=number of FFT to do in //
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! npwin=number of elements in fofgin array (for option 0, 1 and 2)
!! npwout=number of elements in fofgout array (for option 2 and 3)
!! n4,n5,n6=ngfft(4),ngfft(5),ngfft(6), dimensions of fofr.
!! option= if 0: do direct FFT
!!         if 1: do direct FFT, then sum the density
!!         if 2: do direct FFT, multiply by the potential, then do reverse FFT
!!         if 3: do reverse FFT only
!! paral_kgb=Flag related to the kpoint-band-fft parallelism
!! tim_fourwf=timing code of the calling routine (can be set to 0 if not attributed)
!! weight_r=weight to be used for the accumulation of the density in real space
!!         (needed only when option=1)
!! weight_i=weight to be used for the accumulation of the density in real space
!!         (needed only when option=1 and (fftalg=4 and fftalgc/=0))
!! fofginb(2,npwin)=holds second input wavefunction in G vector basis sphere.
!!                 (intent(in) but the routine sphere can modify it for another iflag)
!!                 (for non diagonal occupation)
!! use_ndo = use non diagonal occupations.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! for option==0, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                fofr(2,n4,n5,n6) contains the output Fourier Transform of fofgin;
!!                no use of denpot, fofgout and npwout.
!! for option==1, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                denpot(cplex*n4,n5,n6) contains the input density at input,
!!                and the updated density at output (accumulated);
!!                no use of fofgout and npwout.
!! for option==2, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                denpot(cplex*n4,n5,n6) contains the input local potential;
!!                fofgout(2,npwout*ndat) contains the output function;
!! for option==3, fofr(2,n4,n5,n6*ndat) contains the input real space wavefunction;
!!                fofgout(2,npwout*ndat) contains its output Fourier transform;
!!                no use of fofgin and npwin.
!!
!! PARENTS
!!      accrho3,dens_in_sph,getgh1c,getghc,ladielmt,lavnl,m_fft_prof,m_wfs
!!      mkrho,mkrho3,mlwfovlp,overlap_wf,pawmkaewf,pawsushat,prctfvw1,prctfvw2
!!      prep_fourwf,resp3dte,spin_current,susk,susk_dyn,susk_dyn_pgg
!!      susk_kxc_dyn,suskmm,suskmm_dyn,suskmm_kxc_dyn,tddft
!!      vso_realspace_nonlop,vtowfk,wffile,wfkfermi3,wfread
!!
!! CHILDREN
!!      ccfft,fftw3_fourwf,gpu_fourwf,sg_fftpad,sg_fftrisc,sg_fftrisc_2
!!      sg_fourwf,sphere,sphere_fft,timab,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&  kg_kin,kg_kout,mgfft,mpi_enreg,ndat,ngfft,npwin,npwout,n4,n5,n6,option,&
&  paral_kgb,tim_fourwf,weight_r,weight_i,&
&  use_gpu_cuda,use_ndo,fofginb) ! Optional arguments

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fourwf'
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_ffts, except_this_one => fourwf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,mgfft,n4,n5,n6,ndat,npwin,npwout,option,paral_kgb
 integer,intent(in) :: tim_fourwf
 integer,intent(in),optional :: use_gpu_cuda,use_ndo
 real(dp),intent(in) :: weight_r,weight_i
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(dp),intent(inout) :: denpot(cplex*n4,n5,n6),fofgin(2,npwin*ndat)
! real(dp) :: fofginb(2,npwin*ndat)
 real(dp),intent(inout),optional :: fofginb(:,:)
 real(dp),intent(inout) :: fofr(2,n4,n5,n6*ndat)
 real(dp),intent(out) :: fofgout(2,npwout*ndat)

!Local variables-------------------------------
!scalars
 integer :: fftalg,fftalga,fftalgc,fftcache,i1,i2,i3,idat,ier
 integer :: iflag,ig,inplace,isign
 integer :: me_fft,n1,n2,n3,nd2proc,nd3proc
 integer :: nfftot,normalized,nproc_fft,old_paral_level,option_ccfft,spaceComm
 real(dp) :: fim,fre,xnorm
 character(len=500) :: message
 logical :: luse_gpu_cuda,luse_ndo
!arrays
 integer :: shiftg(3),symm(3,3)
 real(dp) :: tsec(2)
 real(dp),allocatable :: work1(:,:,:,:),work2(:,:,:,:),work3(:,:,:,:)
 real(dp),allocatable :: work4(:,:,:,:),work_sum(:,:,:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)cplex,fofgout,fofr,gboundin,gboundout,istwf_k,kg_kout,mgfft,ndat,ngfft,npwin,npwout,n4,n5,n6,option,tim_fourwf,weight_r
!do ipw=1,npwin
!write(std_out,'(4i4,2es16.6)')ipw,kg_kin(1:3,ipw),fofgin(1:2,ipw)
!end do
!stop
!ENDDEBUG

!Accumulate timing
 call timab(840+tim_fourwf,1,tsec)

!Cuda version of fourwf
 luse_gpu_cuda=PRESENT(use_gpu_cuda)
 if (luse_gpu_cuda) luse_gpu_cuda=(luse_gpu_cuda.and.(use_gpu_cuda==1))
 if(luse_gpu_cuda) then
#if defined HAVE_GPU_CUDA
   call gpu_fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&   kg_kin,kg_kout,mgfft,mpi_enreg,ndat,ngfft,npwin,npwout,n4,n5,n6,option,&
&   paral_kgb,tim_fourwf,weight_r,weight_i) !,&
!  &  use_ndo,fofginb)
#endif
   call timab(840+tim_fourwf,2,tsec)
   return
 end if

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3) ; nfftot=n1*n2*n3
 fftcache=ngfft(8)
 fftalg=ngfft(7)
 fftalga=fftalg/100 ; fftalgc=mod(fftalg,10)
 me_fft=ngfft(11) ; nproc_fft=ngfft(10)
!rewind(72); read(72,*) fftalgc

 if(fftalgc<0 .or. fftalgc>2)then
   write(message, '(a,i4,a,a,a,a,a)' )&
&   '  The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&   '  The third digit, fftalg(C), must be 0, 1, or 2',ch10,&
&   '  Action : change fftalg in your input file.'
   MSG_ERROR(message)
 end if

 luse_ndo=.false.
 if(present(use_ndo).and.present(fofginb)) then
   if(use_ndo==1) luse_ndo=.true.
 end if
 if(luse_ndo) then
   if((size(fofginb,2)==0).and.(use_ndo==1)) then
     write(message, '(a,a,a,i4,i5)' )&
&     '  fofginb has a dimension equal to zero and use_ndo==1',ch10,&
&     '  Action : check dimension of fofginb',size(fofginb,2),use_ndo
     MSG_ERROR(message)
   end if
 end if

 if((fftalga/=1.or.fftalgc/=2.or.option/=1).and.luse_ndo)  then
   write(message, '(a,i5,a,a)' )&
&   '  This algorithm for fft is not yet compatible with luse_ndo (dmft)',fftalg,ch10,&
&   '  Action : check fftalg'
   MSG_ERROR(message)
 end if

 if(fftalgc/=0 .and. ALL(fftalga/=(/1,3,4/)) ) then
   write(message, '(a,i4,5a)' )&
&   '  The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&   '  The first digit must be 1,3,4 when the last digit is not 0.',ch10,&
&   '  Action : change fftalg in your input file.'
   MSG_ERROR(message)
 end if

 if( option<0 .or. option>3 )then
   write(message, '(a,i4,a,a,a)' )&
&   '  The option number',option,' is not allowed.',ch10,&
&   '  Only option=0, 1, 2 or 3 are allowed presently.'
   MSG_ERROR(message)
 end if

 if( option==1 .and. cplex/=1 )then
   write(message, '(a,a,a,i4,a)' )&
&   '  With the option number 1, cplex must be 1,',ch10,&
&   '  but it is cplex=',cplex,'.'
   MSG_ERROR(message)
 end if

 if( option==2 .and. (cplex/=1 .and. cplex/=2) )then
   write(message, '(a,a,a,i4,a)' )&
&   '  With the option number 2, cplex must be 1 or 2,',ch10,&
&   '  but it is cplex=',cplex,'.'
   MSG_ERROR(message)
 end if

 shiftg(:)=0
 symm(:,:)=0
 symm(1,1)=1 ; symm(2,2)=1 ; symm(3,3)=1

!Now choose the algorithm
 select case (fftalg)

!  #if defined HAVE_FFT_FFTW3
   case (312)
     if(luse_ndo) then
       call fftw3_fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&       kg_kin,kg_kout,mgfft,mpi_enreg,ndat,ngfft,npwin,npwout,n4,n5,n6,option,paral_kgb,weight_r,weight_i,&
&       use_ndo=use_ndo,fofginb=fofginb) ! optional Arguments. 
     else 
       call fftw3_fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&       kg_kin,kg_kout,mgfft,mpi_enreg,ndat,ngfft,npwin,npwout,n4,n5,n6,option,paral_kgb,weight_r,weight_i)
     end if
!    RETURN   
!    ! XG100720 One should not have a RETURN inside a routine, but on the contrary, exit properly at the end of the routine.
!    e.g. With this RETURN, the timing is corrupted, as timab is not called ...
!    #endif
     case default 

!    ------------------------------------------------------------------
!    Here, use routines that make forwards FFT separately of backwards FFT,
!    in particular, usual 3DFFT library routines, called in ccfft.
     if( fftalgc==0 .or. &
&     (fftalgc==1 .and. fftalga/=4) .or. &
&     (fftalgc==2 .and. fftalga/=4 .and. option==3) )then
       ABI_ALLOCATE(work1,(2,n4,n5,n6*ndat))

       if(option/=3)then ! Insert fofgin into the fft box (array fofr)
         if(fftalga/=4)then
           iflag=1
           call sphere(fofgin,ndat,npwin,fofr,n1,n2,n3,n4,n5,n6,&
&           kg_kin,istwf_k,iflag,mpi_enreg,shiftg,symm,one)
         else if(fftalga==4 .and. fftalgc==0)then
!          write(std_out,*)'fftalga,fftalgc',fftalga,fftalgc
           iflag=2
           ABI_ALLOCATE(work2,(2,n4,n6,n5*ndat))
!          Note the switch of n5 and n6, as they are only
!          needed to dimension work2 inside "sphere"
           nd2proc=((n2-1)/nproc_fft) +1
           nd3proc=((n6-1)/nproc_fft) +1
           ABI_ALLOCATE(work3,(2,n4,n6,nd2proc*ndat))
           ABI_ALLOCATE(work4,(2,n4,n5,nd3proc*ndat))
           if (istwf_k == 1 .and. mpi_enreg%paral_compil_fft==1) then
!            sphere dont need a big array
             work3(:,:,:,:)=0.d0
             call sphere_fft(fofgin,ndat,npwin,work3,n1,n2,n3,n4,n6,&
&             kg_kin,mpi_enreg,nd2proc)
           else
!            sphere needs a big array and communications
             if (nproc_fft == 1 .and. ndat == 1 .and. istwf_k == 1) then
!              dimensions of tab work3 and work2 are identical
!              no need to use work2
               work3(:,:,:,:)=0.d0
!              DEBUG
!              write(std_out,*) 'in fourwf, call sphere',n4,n6,nd2proc,ndat
!              ENDDEBUG
               call sphere(fofgin,ndat,npwin,work3,n1,n2,n3,n4,n6,n5,&
&               kg_kin,istwf_k,iflag,mpi_enreg,shiftg,symm,one)
             else
               work2(:,:,:,:)=0.d0
               call sphere(fofgin,ndat,npwin,work2,n1,n2,n3,n4,n6,n5,&
&               kg_kin,istwf_k,iflag,mpi_enreg,shiftg,symm,one)
               if ( mpi_enreg%paral_compil_fft==1 .and. istwf_k > 1 ) then
                 work3(:,:,:,:)=0.d0
                 old_paral_level=mpi_enreg%paral_level
                 mpi_enreg%paral_level=3
                 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
                 ABI_ALLOCATE(work_sum,(2,n4,n6,n5*ndat))
                 call timab(48,1,tsec)
                 call xsum_mpi(work2,work_sum,2*n4*n6*n5*ndat,spaceComm,ier)
                 call timab(48,2,tsec)
                 mpi_enreg%paral_level=old_paral_level
                 do idat=1,ndat
                   do i2=1,n2
                     if (me_fft==modulo(i2,nproc_fft)) then
                       do i3=1,n3
                         do i1=1,n1
                           work3(1,i1,i3,(i2-1)/nproc_fft +1+nd2proc*(idat-1))=work_sum(1,i1,i3,i2+n5*(idat-1))
                           work3(2,i1,i3,(i2-1)/nproc_fft +1+nd2proc*(idat-1))=work_sum(2,i1,i3,i2+n5*(idat-1))
                         end do
                       end do
                     end if
                   end do
                 end do
                 ABI_DEALLOCATE(work_sum)
               end if
               if (mpi_enreg%paral_compil_fft/=1) then
                 do idat=1,ndat
                   do i2=1,n2
                     do i3=1,n3
                       do i1=1,n1
                         work3(1,i1,i3,i2+nd2proc*(idat-1))=work2(1,i1,i3,i2+n5*(idat-1))
                         work3(2,i1,i3,i2+nd2proc*(idat-1))=work2(2,i1,i3,i2+n5*(idat-1))
                       end do
                     end do
                   end do
                 end do
               end if
             end if
           end if
           if (mpi_enreg%paral_compil_fft==1) then
             option_ccfft=1
           else
             option_ccfft=2
           end if
         end if
         isign=1
!        Fourier transform fofr (reciprocal to real space)
!        The output might be in work1 or fofr, depending on inplace
         if(fftalgc==0)then
           if(fftalga/=4)then ! Call usual 3DFFT library routines 
             call ccfft(fftalga,fftcache,inplace,isign,mpi_enreg,normalized,&
&             n1,n2,n3,n4,n5,n6,ndat,2,paral_kgb,fofr,work1)
           else               ! SG simplest complex-to-complex routine
             call ccfft(fftalga,fftcache,inplace,isign,mpi_enreg,normalized,&
&             n1,n2,n3,n4,n5,n6,ndat,option_ccfft,paral_kgb,work3,work4)
             ABI_DEALLOCATE(work2)
             ABI_DEALLOCATE(work3)
           end if
         else
           inplace=0 ; normalized=0
!          Call SG routine, with zero padding
!          DEBUG
!          write(std_out,*)' fourwf : before sg_fftpad'
!          do i3=1,n3
!          do i2=1,n2
!          do i1=1,n1
!          write(std_out,'(3i4,2es16.6)')i1,i2,i3,fofr(1:2,i1,i2,i3)
!          end do
!          end do
!          end do
!          stop
!          ENDDEBUG
           call sg_fftpad(fftcache,mgfft,n4,n5,n6,n1,n2,n3,fofr,work1,one,gboundin)
!          DEBUG
!          write(std_out,*)' fourwf : after sg_fftpad'
!          do i3=1,n3
!          do i2=1,n2
!          do i1=1,n1
!          write(std_out,'(3i4,2es16.6)')i1,i2,i3,work1(1:2,i1,i2,i3)
!          end do
!          end do
!          end do
!          stop
!          ENDDEBUG

         end if
       end if ! option/=3

       if(option==0 .and. inplace==0)then ! Results are copied to fofr.
!        fofr(:,:,:,:)=0.d0
!        $OMP PARALLEL DO PRIVATE(i1,i2,i3) SHARED(fofr,n1,n2,n3,work1)
         do i3=1,n3
           do i2=1,n2
             do i1=1,n1
               fofr(1,i1,i2,i3)=work1(1,i1,i2,i3)
               fofr(2,i1,i2,i3)=work1(2,i1,i2,i3)
             end do
           end do
         end do
!        $OMP END PARALLEL DO
       end if ! option==0 and inplace==0
!      Note that if option==0 everything is alright already, the output is available in fofr.

       if(option==1)then ! Accumulate density
         if(inplace==0)then
           if ((fftalgc==0) .and. (fftalga==4)) then
             do idat=1,ndat
!              $OMP PARALLEL DO PRIVATE(i1,i2,i3) SHARED(idat,denpot,n1,n2,n3,nproc_fft,me_fft,nd3proc,weight_r,work4)
               do i3=1,n3
                 if (((i3-1)/(n3/nproc_fft))==me_fft) then
                   do i2=1,n2
                     do i1=1,n1
                       denpot(i1,i2,i3)=denpot(i1,i2,i3)+&
&                       weight_r*(work4(1,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))**2+&
&                       work4(2,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))**2)
                     end do
                   end do
                 end if
               end do
!              $OMP END PARALLEL DO
             end do
           else
!            $OMP PARALLEL DO PRIVATE(i1,i2,i3) SHARED(denpot,n1,n2,n3,weight_r,work1)
             do i3=1,n3
               do i2=1,n2
                 do i1=1,n1
                   denpot(i1,i2,i3)=denpot(i1,i2,i3)+&
&                   weight_r*(work1(1,i1,i2,i3)**2+work1(2,i1,i2,i3)**2)
                 end do
               end do
             end do
!            $OMP END PARALLEL DO
           end if
         else if(inplace==1)then
!          $OMP PARALLEL DO PRIVATE(i1,i2,i3) SHARED(denpot,fofr,n1,n2,n3,weight_r)
           do i3=1,n3
             do i2=1,n2
               do i1=1,n1
                 denpot(i1,i2,i3)=denpot(i1,i2,i3)+&
&                 weight_r*(fofr(1,i1,i2,i3)**2+fofr(2,i1,i2,i3)**2)
               end do
             end do
           end do
!          $OMP END PARALLEL DO
         end if
       end if ! option==1

       if(option==2)then ! Apply local potential
         if(cplex==1)then
           if(inplace==0)then
             if ((fftalgc==0) .and. (fftalga==4)) then

!              $OMP PARALLEL DO PRIVATE(i1,i2,i3,idat) SHARED(denpot,fofr,n1,n2,n3,ndat,me_fft,nproc_fft,nd3proc,work1)
               do idat=1,ndat
                 do i3=1,n3
                   if (((i3-1)/(n3/nproc_fft))==me_fft) then
                     do i2=1,n2
                       do i1=1,n1
                         fofr(1,i1,i2,i3+n3*(idat-1))=&
&                         denpot(i1,i2,i3)*work4(1,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))
                         fofr(2,i1,i2,i3+n3*(idat-1))=&
&                         denpot(i1,i2,i3)*work4(2,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))
!                        write(std_out,'(3i3,3e24.10)') i1,i2,i3,denpot(i1,i2,i3),fofr(:,i1,i2,i3+n3*(idat-1))
                       end do
                     end do
                   end if
                 end do
               end do
!              $OMP END PARALLEL DO
             end if
             if ((fftalgc/=0) .or. (fftalga/=4)) then
!              $OMP PARALLEL DO PRIVATE(i1,i2,i3,idat) SHARED(denpot,fofr,n1,n2,n3,ndat,me_fft,nproc_fft,work1)
               do idat=1,ndat
                 do i3=1,n3
                   if (((i3-1)/(n3/nproc_fft))==me_fft) then
                     do i2=1,n2
                       do i1=1,n1
                         fofr(1,i1,i2,i3+n3*(idat-1))=denpot(i1,i2,i3)*work1(1,i1,i2,i3+n3*(idat-1))
                         fofr(2,i1,i2,i3+n3*(idat-1))=denpot(i1,i2,i3)*work1(2,i1,i2,i3+n3*(idat-1))
                       end do
                     end do
                   end if
                 end do
               end do
!              $OMP END PARALLEL DO
             end if
           else if(inplace==1)then
!            $OMP PARALLEL DO PRIVATE(i1,i2,i3,idat) SHARED(denpot,fofr,n1,n2,n3,ndat,me_fft,nproc_fft)
             do idat=1,ndat
               do i3=1,n3
                 if (((i3-1)/(n3/nproc_fft))==me_fft) then
                   do i2=1,n2
                     do i1=1,n1
                       fofr(1,i1,i2,i3+n3*(idat-1))=denpot(i1,i2,i3)*fofr(1,i1,i2,i3+n3*(idat-1))
                       fofr(2,i1,i2,i3+n3*(idat-1))=denpot(i1,i2,i3)*fofr(2,i1,i2,i3+n3*(idat-1))
                     end do
                   end do
                 end if
               end do
             end do
!            $OMP END PARALLEL DO
           end if
         else if(cplex==2)then
           if(inplace==0)then
             if ((fftalgc==0) .and. (fftalga==4)) then
!              $OMP PARALLEL DO PRIVATE(fre,fim,i1,i2,i3,idat) SHARED(denpot,fofr,n1,n2,n3,ndat,me_fft,nproc_fft,work1)
               do idat=1,ndat
                 do i3=1,n3
                   if (((i3-1)/(n3/nproc_fft))==me_fft) then
                     do i2=1,n2
                       do i1=1,n1
                         fre=work4(1,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))
                         fim=work4(2,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))
                         fofr(1,i1,i2,i3+n3*(idat-1))=denpot(2*i1-1,i2,i3)*fre -denpot(2*i1,i2,i3)*fim
                         fofr(2,i1,i2,i3+n3*(idat-1))=denpot(2*i1-1,i2,i3)*fim +denpot(2*i1,i2,i3)*fre
                       end do
                     end do
                   end if
                 end do
               end do
!              $OMP END PARALLEL DO
             end if
             if ((fftalgc/=0) .or. (fftalga/=4)) then
!              $OMP PARALLEL DO PRIVATE(fre,fim,i1,i2,i3,idat) SHARED(denpot,fofr,n1,n2,n3,ndat,me_fft,nproc_fft,work1)
               do idat=1,ndat
                 do i3=1,n3
                   if (((i3-1)/(n3/nproc_fft))==me_fft) then
                     do i2=1,n2
                       do i1=1,n1
                         fre=work1(1,i1,i2,i3+n3*(idat-1))
                         fim=work1(2,i1,i2,i3+n3*(idat-1))
                         fofr(1,i1,i2,i3+n3*(idat-1))=denpot(2*i1-1,i2,i3)*fre -denpot(2*i1,i2,i3)*fim
                         fofr(2,i1,i2,i3+n3*(idat-1))=denpot(2*i1-1,i2,i3)*fim +denpot(2*i1,i2,i3)*fre
                       end do
                     end do
                   end if
                 end do
               end do
!              $OMP END PARALLEL DO
             end if
           else if(inplace==1)then
!            $OMP PARALLEL DO PRIVATE(fre,fim,i1,i2,i3,idat) SHARED(denpot,fofr,n1,n2,n3,ndat,me_fft,nproc_fft)
             do idat=1,ndat
               do i3=1,n3
                 if (((i3-1)/(n3/nproc_fft))==me_fft) then
                   do i2=1,n2
                     do i1=1,n1
                       fre=fofr(1,i1,i2,i3+n3*(idat-1))
                       fim=fofr(2,i1,i2,i3+n3*(idat-1))
                       fofr(1,i1,i2,i3+n3*(idat-1))=denpot(2*i1-1,i2,i3)*fre -denpot(2*i1,i2,i3)*fim
                       fofr(2,i1,i2,i3+n3*(idat-1))=denpot(2*i1-1,i2,i3)*fim +denpot(2*i1,i2,i3)*fre
                     end do
                   end do
                 end if
               end do
             end do
!            $OMP END PARALLEL DO
           end if ! inplace

         end if ! cplex=2

       end if ! option=2

!      The data for option==2 or option==3 is now in fofr.
       if(option==2 .or. option==3)then
         isign=-1

         if(fftalgc==0)then ! Call usual 3DFFT library routines or SG simplest complex-to-complex routine
           if(fftalga==4)then
             ABI_DEALLOCATE(work1)
             ABI_ALLOCATE(work1,(2,n4,n6,n5*ndat))
           end if
           if(option==3 .or. fftalga/=4) then
             call ccfft(fftalga,fftcache,inplace,isign,mpi_enreg,normalized,&
&             n1,n2,n3,n4,n5,n6,ndat,2,paral_kgb,fofr,work1)
           else
!            creation of small arrays
!            !!!        nd3proc=((n5-1)/nproc_fft) +1
             nd3proc=((n6-1)/nproc_fft) +1
             nd2proc=((n2-1)/nproc_fft) +1
             ABI_ALLOCATE(work3,(2,n4,n5,nd3proc*ndat))
             ABI_ALLOCATE(work2,(2,n4,n6,nd2proc*ndat))
!            work3(:,:,:,:)=0.d0
             if (mpi_enreg%paral_compil_fft==1) then
               if (inplace==0) then
                 if (cplex==1) then
                   do idat=1,ndat
                     do i3=1,n3
                       if (((i3-1)/(n3/nproc_fft))==me_fft) then
                         do i2=1,n2
                           do i1=1,n1
                             work3(1,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))=&
&                             denpot(i1,i2,i3)*work4(1,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))
                             work3(2,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))=&
&                             denpot(i1,i2,i3)*work4(2,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))
                           end do
                         end do
                       end if
                     end do
                   end do
                 else
                   do idat=1,ndat
                     do i3=1,n3
                       if (((i3-1)/(n3/nproc_fft))==me_fft) then
                         do i2=1,n2
                           do i1=1,n1
                             fre=work4(1,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))
                             fim=work4(2,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))
                             work3(1,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))=&
&                             denpot(2*i1-1,i2,i3)*fre-denpot(2*i1  ,i2,i3)*fim
                             work3(2,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))=&
&                             denpot(2*i1-1,i2,i3)*fim+denpot(2*i1  ,i2,i3)*fre
                           end do
                         end do
                       end if
                     end do
                   end do
                 end if
               else !inplace /= 0
                 do idat=1,ndat
                   do i3=1,n3
                     if (((i3-1)/(n3/nproc_fft))==me_fft) then
                       do i2=1,n2
                         do i1=1,n1
                           work3(1,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))=fofr(1,i1,i2,i3+n3*(idat-1))
                           work3(2,i1,i2,i3-me_fft*n3/nproc_fft+nd3proc*(idat-1))=fofr(2,i1,i2,i3+n3*(idat-1))
                         end do
                       end do
                     end if
                   end do
                 end do
               end if ! if inplace
               option_ccfft=1
             else !mpi_enreg%paral_compil_fft/=1
               if (nproc_fft /=1 .or. ndat /= 1 ) then
                 do idat=1,ndat
                   do i3=1,n3
                     do i2=1,n2
                       do i1=1,n1
                         work3(1,i1,i2,i3+nd3proc*(idat-1))=fofr(1,i1,i2,i3+n3*(idat-1))
                         work3(2,i1,i2,i3+nd3proc*(idat-1))=fofr(2,i1,i2,i3+n3*(idat-1))
                       end do
                     end do
                   end do
                 end do
                 option_ccfft=2
               end if
             end if
             if (mpi_enreg%paral_compil_fft==1) then
!              DEBUG
!              write(std_out,*) 'before hitting cfft, i1,i2,i3,work3'
!              do i3=1,(n3-1)/mpi_enreg%nproc_fft+1
!              do i2=1,n2
!              do i1=1,n1
!              write(std_out,'(3i4,2es16.6)') i1,i2,i3,work3(:,i1,i2,i3)
!              end do
!              end do
!              end do
!              call leave_new("COLL")
!              ENDDEBUG
               call ccfft(fftalga,fftcache,inplace,isign,mpi_enreg,normalized,&
&               n1,n2,n3,n4,n5,n6,ndat,option_ccfft,paral_kgb,work3,work2)
             else
               if (nproc_fft /=1 .or. ndat /= 1 ) then
                 call ccfft(fftalga,fftcache,inplace,isign,mpi_enreg,normalized,&
&                 n1,n2,n3,n4,n5,n6,ndat,option_ccfft,paral_kgb,work3,work2)
               else
                 call ccfft(fftalga,fftcache,inplace,isign,mpi_enreg,normalized,&
&                 n1,n2,n3,n4,n5,n6,ndat,option_ccfft,paral_kgb,fofr,work1)
               end if
             end if
!            load of work1
             if ((mpi_enreg%paral_compil_fft==1) .and.  ( istwf_k > 1 )) work1(:,:,:,:)=0.d0
             if (mpi_enreg%paral_compil_fft==1) then
               if ( istwf_k > 1 ) then
                 do idat=1,ndat
                   do i2=1,n2
                     if (me_fft==modulo(i2,nproc_fft)) then
                       do i3=1,n3
                         do i1=1,n1
                           work1(1,i1,i3,i2+n5*(idat-1))=&
&                           work2(1,i1,i3,(i2-1)/nproc_fft +1+nd2proc*(idat-1))
                           work1(2,i1,i3,i2+n5*(idat-1))=&
&                           work2(2,i1,i3,(i2-1)/nproc_fft +1+nd2proc*(idat-1))
                         end do
                       end do
                     end if
                   end do
                 end do
               end if
             else
               if (nproc_fft /=1 .or. ndat /= 1 ) then
                 do idat=1,ndat
                   do i2=1,n2
                     do i3=1,n3
                       do i1=1,n1
                         work1(1,i1,i3,i2+n5*(idat-1))=work2(1,i1,i3,i2+nd2proc*(idat-1))
                         work1(2,i1,i3,i2+n5*(idat-1))=work2(2,i1,i3,i2+nd2proc*(idat-1))
!                        write(std_out,'(3i3,2e24.10)') i1,i3,i2+n5*(idat-1),work1(:,i1,i3,i2+n5*(idat-1))
                       end do
                     end do
                   end do
                 end do
               end if
             end if
             ABI_DEALLOCATE(work3)
             if ((mpi_enreg%paral_compil_fft==1) .and.  ( istwf_k > 1 )) then
               old_paral_level=mpi_enreg%paral_level
               mpi_enreg%paral_level=3
               call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
               call timab(48,1,tsec)
               call xsum_mpi(work1,spaceComm,ier)
               call timab(48,2,tsec)
               mpi_enreg%paral_level=old_paral_level
             end if
           end if
         else
           inplace=0 ; normalized=0 ! Call SG routine, with zero padding
           call sg_fftpad(fftcache,mgfft,n4,n5,n6,n1,n2,n3,fofr,work1,-one,gboundout)
         end if
         xnorm=1.d0/dble(nfftot)
         if(normalized==0)then
           if(inplace==0)then
             if(fftalga/=4)then
!              $OMP PARALLEL DO PRIVATE(i1,i2,i3,ig,idat) SHARED(ndat,fofgout,kg_kout,npwout,work1,xnorm)
               do ig=1,npwout
                 i1=kg_kout(1,ig); if(i1<0)i1=i1+n1 ; i1=i1+1
                 i2=kg_kout(2,ig); if(i2<0)i2=i2+n2 ; i2=i2+1
                 i3=kg_kout(3,ig); if(i3<0)i3=i3+n3 ; i3=i3+1
                 do idat=1,ndat
                   fofgout(1,ig)=work1(1,i1,i2,i3)*xnorm
                   fofgout(2,ig)=work1(2,i1,i2,i3)*xnorm
                 end do
               end do
!              $OMP END PARALLEL DO
             else ! if fftalga==4
               if ((mpi_enreg%paral_compil_fft==1) .and. ( istwf_k == 1 )) then
!                $OMP PARALLEL DO PRIVATE(i1,i2,i3,ig,idat) SHARED(fofgout,kg_kout,npwout,ndat,n1,n2,n3,nproc_fft,nd2proc,work2,xnorm)
                 do ig=1,npwout
                   i1=kg_kout(1,ig); if(i1<0)i1=i1+n1 ; i1=i1+1
                   i2=kg_kout(2,ig); if(i2<0)i2=i2+n2 ; i2=i2+1
                   i3=kg_kout(3,ig); if(i3<0)i3=i3+n3 ; i3=i3+1
                   do idat=1,ndat
                     fofgout(1,ig+npwout*(idat-1))=&
&                     work2(1,i1,i3,(i2-1)/nproc_fft +1+nd2proc*(idat-1))*xnorm
                     fofgout(2,ig+npwout*(idat-1))=&
&                     work2(2,i1,i3,(i2-1)/nproc_fft +1+nd2proc*(idat-1))*xnorm
!                    write(std_out,'(4i3,2e24.10)') i1,i2,i3,ig,fofgout(:,ig+npwout*(idat-1))
                   end do
                 end do
!                $OMP END PARALLEL DO
                 ABI_DEALLOCATE(work2)
               else
!                $OMP PARALLEL DO PRIVATE(i1,i2,i3,ig,idat) SHARED(fofgout,kg_kout,npwout,ndat,n1,n2,n3,n5,work1,xnorm)
                 do ig=1,npwout
                   i1=kg_kout(1,ig); if(i1<0)i1=i1+n1 ; i1=i1+1
                   i2=kg_kout(2,ig); if(i2<0)i2=i2+n2 ; i2=i2+1
                   i3=kg_kout(3,ig); if(i3<0)i3=i3+n3 ; i3=i3+1
                   do idat=1,ndat
                     fofgout(1,ig+npwout*(idat-1))=work1(1,i1,i3,i2+n5*(idat-1))*xnorm
                     fofgout(2,ig+npwout*(idat-1))=work1(2,i1,i3,i2+n5*(idat-1))*xnorm
                   end do
                 end do
               end if
!              $OMP END PARALLEL DO
             end if ! fftalga
           else if(inplace==1)then
!            $OMP PARALLEL DO PRIVATE(i1,i2,i3,ig,idat) SHARED(fofgout,fofr,kg_kout,ndat,n1,n2,n3,npwout,xnorm)
             do ig=1,npwout
               i1=kg_kout(1,ig); if(i1<0)i1=i1+n1 ; i1=i1+1
               i2=kg_kout(2,ig); if(i2<0)i2=i2+n2 ; i2=i2+1
               i3=kg_kout(3,ig); if(i3<0)i3=i3+n3 ; i3=i3+1
               do idat=1,ndat
                 fofgout(1,ig+npwout*(idat-1))=fofr(1,i1,i2,i3+n3*(idat-1))*xnorm
                 fofgout(2,ig+npwout*(idat-1))=fofr(2,i1,i2,i3+n3*(idat-1))*xnorm
               end do
             end do
!            $OMP END PARALLEL DO
           end if ! inplace
         else if(normalized==1)then
           if(inplace==0)then
!            $OMP PARALLEL DO PRIVATE(i1,i2,i3,ig,idat) SHARED(fofgout,kg_kout,ndat,n1,n2,n3,npwout,work1)
             do ig=1,npwout
               i1=kg_kout(1,ig); if(i1<0)i1=i1+n1 ; i1=i1+1
               i2=kg_kout(2,ig); if(i2<0)i2=i2+n2 ; i2=i2+1
               i3=kg_kout(3,ig); if(i3<0)i3=i3+n3 ; i3=i3+1
               do idat=1,ndat
                 fofgout(1,ig+npwout*(idat-1))=work1(1,i1,i2,i3+n3*(idat-1))
                 fofgout(2,ig+npwout*(idat-1))=work1(2,i1,i2,i3+n3*(idat-1))
               end do
             end do
!            $OMP END PARALLEL DO
           else if(inplace==1)then
!            $OMP PARALLEL DO PRIVATE(i1,i2,i3,ig,idat) SHARED(fofgout,fofr,n1,n2,n3,kg_kout,npwout)
             do ig=1,npwout
               i1=kg_kout(1,ig); if(i1<0)i1=i1+n1 ; i1=i1+1
               i2=kg_kout(2,ig); if(i2<0)i2=i2+n2 ; i2=i2+1
               i3=kg_kout(3,ig); if(i3<0)i3=i3+n3 ; i3=i3+1
               do idat=1,ndat
                 fofgout(1,ig+npwout*(idat-1))=fofr(1,i1,i2,i3+n3*(idat-1))
                 fofgout(2,ig+npwout*(idat-1))=fofr(2,i1,i2,i3+n3*(idat-1))
               end do
             end do
!            $OMP END PARALLEL DO
           end if ! inplace

         end if ! normalized
       end if ! if option==2 or 3

       if (allocated(work1))  then
         ABI_DEALLOCATE(work1)
       end if
     end if

!    ------------------------------------------------------------------
!    Here, call more sophisticated specialized 3-dimensional fft
!    (zero padding as well as maximize cache reuse) based on S Goedecker routines.
!    Specially tuned for cache architectures.
     if( fftalga==1 .and. fftalgc==2 .and. option/=3)then
       if(luse_ndo) then
         call sg_fftrisc_2(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,&
&         istwf_k,kg_kin,kg_kout,&
&         mgfft,ngfft,npwin,npwout,n4,n5,n6,option,weight_r,luse_ndo=luse_ndo,fofgin_p=fofginb)
       else
!        Note that the arguments are the same as for fourwf
         call sg_fftrisc(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,&
&         istwf_k,kg_kin,kg_kout,&
&         mgfft,ngfft,npwin,npwout,n4,n5,n6,option,weight_r)
       end if
     end if
!    ------------------------------------------------------------------

!    Here, call new FFT from S Goedecker, also sophisticated specialized 3-dimensional fft
!    (zero padding as well as maximize cache reuse)
     if(fftalga==4 .and. fftalgc/=0)then
!      The args are not the same as fourwf, but might be
       call sg_fourwf(cplex,denpot,fftalgc,fofgin,fofgout,fofr,gboundin,gboundout,&
&       istwf_k,kg_kin,kg_kout,&
&       mgfft,mpi_enreg,n1,n2,n3,npwin,npwout,n4,n5,n6,option,paral_kgb,weight_r,weight_i)
     end if

     if (allocated(work4))  then
       ABI_DEALLOCATE(work4)
     end if
     if (allocated(work2))  then
       ABI_DEALLOCATE(work2)
     end if

 end select

!------------------------------------------------------------------

!Accumulate timing
 call timab(840+tim_fourwf,2,tsec)

end subroutine fourwf
!!***
