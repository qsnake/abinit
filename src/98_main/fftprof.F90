!{\src2tex{textfont=tt}}
!!****p* ABINIT/fftprof
!! NAME
!! fftprof
!!
!! FUNCTION
!!  Utility for profiling the FFT libraries supported by abinit.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  (main program)
!!
!! OUTPUT
!!  Timing analysis of the different FFT libraries and algorithms.
!!
!! NOTES
!!  Point-group symmetries are not taken into account in getng during the generation
!!  of the FFT mesh. Therefore the FFT mesh might differ from the one
!!  found by abinit for the same cutoff and Bravais lattice (actually it might be smaller).
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,destroy_fft_prof,destroy_fft_test,destroy_mpi_enreg
!!      fftprof_ncalls_per_test,fftw3_gain_wisdom,fftw3_set_in_place,get_kg
!!      herald,init_fft_test,initmpi_seq,lower,metric,mpi_comm_rank
!!      mpi_comm_size,nullify_fft_prof,nullify_fft_test,nullify_mpi_enreg
!!      print_fft_profs,prof_fourdp,prof_fourwf,prof_rhotwg,prompt,time_fourdp
!!      time_fourdp_cplx,time_fourwf,time_padded_fourwf_cplx,time_rhotwg,wrtout
!!      xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program fftprof

 use defs_basis
 use defs_abitypes
 use m_build_info
 use m_xmpi
 use m_errors
 use m_FFT_prof
#if defined HAVE_MPI2
 use mpi
#endif

 use m_fstrings,   only : lower
 use m_io_tools,   only : prompt
 use m_fftw3,      only : fftw3_set_in_place, fftw3_gain_wisdom
 use m_gsphere,    only : get_kg

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftprof'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments -----------------------------------

!Local variables-------------------------------
!scalars
 integer,parameter :: nsym1=1
 integer :: ndat,ierr,isym,fftcache,it,cplex,ntests,option_fourwf,osc_npw
 integer :: n1,n2,n3,use_wisdom,do_outofplace,ncalls,map2sphere,use_padfft,isign
 integer :: iset,iall,inplace,nsets,avail,max_nthreads,ith,idx
 integer :: necut
 character(len=24) :: codename
 character(len=500) :: header,str_tasks
 real(dp),parameter :: boxcutmin=two
 real(dp) :: ecut,ucvol,osc_ecut
 logical :: test_fourdp,test_fourwf,test_gw,do_bench
 logical :: bench_fourdp,bench_fourwf,bench_rhotwg
 type(MPI_type) :: MPI_enreg
!arrays
 integer,allocatable :: symrel(:,:,:)
 real(dp),parameter :: gamma_point(3)=(/zero,zero,zero/)
 real(dp) :: ecut_arth(2)
 real(dp) :: rprimd(3,3),gmet(3,3),gprimd(3,3),rmet(3,3),kpoint(3)
 type(FFT_test_t),allocatable :: Ftest(:)
 type(FFT_prof_t),allocatable :: Ftprof(:)
 integer,pointer :: osc_gvec(:,:)
 integer,allocatable :: fft_setups(:,:),fourwf_params(:,:)

! *************************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

 call xmpi_init()

 call nullify_mpi_enreg(MPI_enreg)
 call initmpi_seq(MPI_enreg)

#ifdef HAVE_MPI
 MPI_enreg%world_comm=xmpi_world
 MPI_enreg%world_group=MPI_GROUP_NULL
 call MPI_COMM_RANK(xmpi_world,MPI_enreg%me,ierr)
 call MPI_COMM_SIZE(xmpi_world,MPI_enreg%nproc,ierr)
 MPI_enreg%paral_compil=1
 MPI_enreg%paral_compil_respfn=0
 MPI_enreg%paral_level=2
#endif

 codename='FFTPROF'//REPEAT(' ',17)
 call herald(codename,abinit_version,std_out)

 write(std_out,'(a)')" Tool for profiling and testing the FFT libraries used in ABINIT."
 write(std_out,'(a)')" Allowed options are: "
 write(std_out,'(a)')"   fourdp --> Test FFT transforms of density and potentials on the full box."
 write(std_out,'(a)')"   fourwf --> Test FFT transforms of wavefunctions using the zero-pad algorithm."
 write(std_out,'(a)')"   gw_fft --> Test the FFT transforms used in the GW code."
 write(std_out,'(a)')"   all    --> Test all FFT routines."
 write(std_out,'(a)')" "

 call prompt("Enter the string defining the routines to be tested: ",str_tasks)
 call lower(str_tasks)

 iall=INDEX (str_tasks,"all")
 test_fourdp = (iall>0 .or. INDEX(str_tasks,"fourdp")>0 )
 test_fourwf = (iall>0 .or. INDEX(str_tasks,"fourwf")>0 )
 test_gw     = (iall>0 .or. INDEX(str_tasks,"gw_fft")>0 )
 do_bench    = INDEX(str_tasks,"bench")>0

 call prompt("Enter cutoff energy ecut in Hartree: ",ecut)
 call prompt("Enter lattice vectors rprimd in Bohr: ",rprimd)

 call metric(gmet,gprimd,std_out,rmet,rprimd,ucvol)
!
!For the time being symmetries are not used for defining the FFT mesh.
 ABI_ALLOCATE(symrel,(3,3,nsym1))
 do isym=1,nsym1
   symrel(:,:,isym) = RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
 end do

 call prompt("Enter number of calls for each test: ",ncalls)
 ncalls=ABS(ncalls); if (ncalls==0) ncalls=10
 call fftprof_ncalls_per_test(ncalls)

 call prompt("Enter the MAXIMUM number of threads: ",max_nthreads)
#ifndef HAVE_OPENMP
 ABI_CHECK(max_nthreads==1,"nthreads>1 but OMP support is not enabled!")
#endif
!
!List the FFT libraries that will be tested.
!Goedecker FFTs are always available, other libs are optional.
 ndat=1; ntests=3*max_nthreads
!
!First dimension contains (fftalg,fftcache,ndat,nthreads,available).
 ABI_ALLOCATE(fft_setups,(5,ntests))
!
!Default Goedecker library.
 idx=0
 do ith=1,max_nthreads
   fftcache = 16  !fftcache is machine-dependent. Here we use the default value specified in indefo.
   idx = idx + 1
   fft_setups(:,idx) = (/112,fftcache,ndat,ith,1/)
 end do
!
!Goedecker2002
 do ith=1,max_nthreads
   fftcache = 16   !fftcache is machine-dependent. Here we use the default value specified in indefo.
   idx = idx + 1
   fft_setups(:,idx) = (/412,fftcache,ndat,ith,1/)
!  fft_setups(:,idx) = (/412,fftcache,ndat,ith,0/)
 end do
!
!FFTW3.
 avail=0
#ifdef HAVE_FFT_FFTW3
 avail=1
#endif
 do ith=1,max_nthreads
   fftcache = 0 ! Not used.
   idx = idx + 1
   fft_setups(:,idx) = (/312,fftcache,ndat,ith,avail/)
 end do

 kpoint = (/0.3,0.5,0.25/)
 call prompt("Enter the reduced coordinates of the k-point for the wavefunction: ",kpoint)
!
!Compute FFT box.
 ABI_ALLOCATE(Ftest,(ntests))
 ABI_ALLOCATE(Ftprof,(ntests))

 do it=1,ntests ! Needed for crappy compilers that do not support => null() in type declarations.
   call nullify_fft_test(Ftest(it))
   call nullify_fft_prof(Ftprof(it))
 end do

 do it=1,ntests
   call init_FFT_test(Ftest(it),fft_setups(:,it),kpoint,ecut,boxcutmin,rprimd,nsym1,symrel,MPI_enreg)
   if ( ANY(Ftest(it)%ngfft(1:3) /= Ftest(1)%ngfft(1:3)) ) then
     MSG_ERROR("Consistency check assumes equal FFT meshes. Cannot continue")
!    Ftest%results is allocated using nfftot and the consistency btw libs is tested assuming an equal number of FFT-points.
   end if
 end do
!
 use_wisdom = 0
!call prompt(" Do you want to to use FFTW3 wisdom?: ",use_wisdom)
 if (use_wisdom/=0) then
   call wrtout(std_out," Gaining FFTW3 wisdom...","COLL")
   n1 = Ftest(1)%ngfft(1)
   n2 = Ftest(1)%ngfft(2)
   n3 = Ftest(1)%ngfft(3)
   call fftw3_gain_wisdom(n1,n2,n3,ndat)
 end if

 do_outofplace = 0
!call prompt("Do you want to run FFTW3 out-of-place?: ",do_outofplace)
 if (do_outofplace /= 0) then
   call fftw3_set_in_place(.FALSE.)
 else
   call fftw3_set_in_place(.TRUE.)
 end if
!
!=======================
!==== fourdp timing ====
!=======================
 if (test_fourdp) then
   do isign=-1,1,2
     do cplex=1,2
       write(header,'(2(a,i2))')" fourdp with cplex ",cplex,", isign ",isign
       do it=1,ntests
         call time_fourdp(Ftest(it),isign,cplex,Ftprof(it))
       end do
       call print_FFT_profs(Ftprof,header)
       call destroy_FFT_prof(Ftprof)
     end do
   end do
 end if
!
!=======================
!==== fourwf timing ====
!=======================
 if (test_fourwf) then
!  possible combinations of (option, cplex) supported in fourwf.
!  (cplex=2 only allowed for option=2, and istwf_k=1)
   nsets=4; if (Ftest(1)%istwf_k==1) nsets=5
   ABI_ALLOCATE(fourwf_params,(2,nsets))
   fourwf_params(:,1) = (/0,0/)
   fourwf_params(:,2) = (/1,1/)
   fourwf_params(:,3) = (/2,1/)
   fourwf_params(:,4) = (/3,0/)
   if (nsets==5) fourwf_params(:,5) = (/2,2/)

   do iset=1,nsets
     option_fourwf = fourwf_params(1,iset)
     cplex         = fourwf_params(2,iset)
     write(header,'(3(a,i2))')" fourwf with option ",option_fourwf,", cplex ",cplex,", istwf_k ",Ftest(1)%istwf_k

     do it=1,ntests
       call time_fourwf(Ftest(it),cplex,option_fourwf,Ftprof(it))
     end do
     call print_FFT_profs(Ftprof,header)
     call destroy_FFT_prof(Ftprof)
   end do
   ABI_DEALLOCATE(fourwf_params)
 end if
!
!==========================
!==== Test GW routines ====
!==========================
!These routines are used in the GW part, FFTW3 is expected to
!be more efficient as the conversion complex(:) <--> real(2,:) is not needed.
!
 if (test_gw) then
!  
!  ==== fourdp timing with complex arrays ====
   do isign=-1,1,2
     do inplace=0,1
       write(header,'(2(a,i2))')" fourdp_cplx with isign ",isign,", in-place ",inplace
       do it=1,ntests
         call time_fourdp_cplx(Ftest(it),isign,inplace,Ftprof(it))
       end do
       call print_FFT_profs(Ftprof,header)
       call destroy_FFT_prof(Ftprof)
     end do
   end do
!  
!  ==== padded_fourwf_cplx ====
   do isign=-1,1,2
     write(header,'(a,i2)')" padded_fourwf_cplx with isign ",isign
     do it=1,ntests
       call time_padded_fourwf_cplx(Ftest(it),isign,Ftprof(it))
     end do
     call print_FFT_profs(Ftprof,header)
     call destroy_FFT_prof(Ftprof)
   end do
!  
!  ==== rho_tw_g timing ====
   call prompt("Enter ecut in Hartree for GW oscillators: ",osc_ecut)

   nullify(osc_gvec)
   call get_kg(gamma_point,1,osc_ecut,gmet,osc_npw,osc_gvec)
!  TODO should reorder by shells to be consistent with the GW part!
!  Moreover I guess this ordering is more efficient when we have
!  to map the box to the G-sphere!
   map2sphere=1; !map2sphere=0

   do use_padfft=0,1
     write(header,'(2(a,i2))')"rho_tw_g with use_padfft ",use_padfft,", map2sphere ",map2sphere
     do it=1,ntests
       call time_rhotwg(Ftest(it),map2sphere,use_padfft,osc_npw,osc_gvec,Ftprof(it))
     end do
     call print_FFT_profs(Ftprof,header)
     call destroy_FFT_prof(Ftprof)
   end do

   ABI_DEALLOCATE(osc_gvec)
 end if ! test_gw

 if (do_bench) then
   MSG_WARNING("Entering benchmark mode")
   call prompt("Enter ecut0 and ecut_step:",ecut_arth)
   call prompt("Enter necut:",necut)

   bench_fourdp = .TRUE.
!  bench_fourdp = .FALSE.
   if (bench_fourdp) then
     isign=1; cplex=1
     call prof_fourdp(fft_setups,isign,cplex,necut,ecut_arth,boxcutmin,rprimd,nsym1,symrel,MPI_enreg)

     isign=1; cplex=2
     call prof_fourdp(fft_setups,isign,cplex,necut,ecut_arth,boxcutmin,rprimd,nsym1,symrel,MPI_enreg)
   end if

   bench_fourwf = .TRUE. 
!  bench_fourwf = .FALSE.
   if (bench_fourwf) then
     cplex=2; option_fourwf=0
     call prof_fourwf(fft_setups,cplex,option_fourwf,kpoint,necut,ecut_arth,boxcutmin,rprimd,nsym1,symrel,MPI_enreg)
   end if

   bench_rhotwg = .TRUE.
!  bench_rhotwg = .FALSE.
   if (bench_rhotwg) then
     map2sphere=0; use_padfft=0
     call prof_rhotwg(fft_setups,map2sphere,use_padfft,necut,ecut_arth,osc_ecut,boxcutmin,&
&     rprimd,nsym1,symrel,gmet,MPI_enreg)
   end if
 end if
!
!===============================
!=== End of run, free memory ===
!===============================
 call wrtout(std_out,ch10//" Analysis completed.","COLL")

 ABI_DEALLOCATE(fft_setups)
 ABI_DEALLOCATE(symrel)

 call destroy_FFT_test(Ftest)
 ABI_DEALLOCATE(Ftest)
 call destroy_FFT_prof(Ftprof)
 ABI_DEALLOCATE(Ftprof)

 call destroy_mpi_enreg(MPI_enreg)
 call xmpi_end()

 end program fftprof
!!***
