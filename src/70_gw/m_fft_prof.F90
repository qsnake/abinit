!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_FFT_prof
!! NAME
!! m_FFT_prof
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_FFT_prof

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_fftw3
 use m_errors 
#ifdef HAVE_OPENMP
 use omp_lib
#endif

 use m_numeric_tools,  only : imax_loc, arth
 use m_io_tools,       only : get_unit
 use m_fft_mesh,       only : print_ngfft, fftalg_info, eigr, ceigr
 use m_gsphere,        only : get_kg
 use m_oscillators,    only : rho_tw_g

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_12_hide_mpi
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

 private
!!***

!----------------------------------------------------------------------

 integer,private,parameter :: TNAME_LEN=100

!----------------------------------------------------------------------

!!****t* m_fft_mesh/FFT_test_t
!! NAME
!! FFT_test_t
!! 
!! FUNCTION
!! Structure storing the set of paramenters passed to the FFT routines used in 
!! abinit (fourdp|fourwf).
!! 
!! SOURCE

 type,public :: FFT_test_t
   integer :: available=0
   integer :: istwf_k=-1
   integer :: mgfft=-1
   integer :: ndat=-1
   integer :: nfft=-1
   integer :: nthreads=1
   integer :: npw_k=-1
   integer :: npw_kout=-1
   integer :: paral_kgb=-1

   real(dp) :: ecut=zero

   integer :: ngfft(18)=-1

   real(dp) :: kpoint(3)=(/zero,zero,zero/)
   real(dp) :: rprimd(3,3),rmet(3,3)
   real(dp) :: gprimd(3,3),gmet(3,3)

   integer,pointer :: kg_k(:,:)     SET2NULL       
   integer,pointer :: kg_kout(:,:)     SET2NULL       
   integer,pointer :: indpw_k(:)  SET2NULL
  
   type(MPI_type) :: MPI_enreg
 end type FFT_test_t   
!!***

 public :: init_FFT_test
 public :: nullify_FFT_test
 public :: destroy_FFT_test
 public :: print_FFT_test

 interface destroy_FFT_prof
   module procedure destroy_FFT_prof_0D
   module procedure destroy_FFT_prof_1D
 end interface destroy_FFT_prof

 interface destroy_FFT_test
   module procedure destroy_FFT_test_0D
   module procedure destroy_FFT_test_1D
 end interface destroy_FFT_test


!----------------------------------------------------------------------

!!****t* m_fft_mesh/FFT_prof_t
!! NAME
!! FFT_prof_t
!! 
!! FUNCTION
!!  The results of the tests 
!! 
!! SOURCE

 type,public :: FFT_prof_t
   integer :: ncalls
   !% integer :: ndat
   integer :: nthreads 
   real(dp) :: cpu_time
   real(dp) :: wall_time
   character(len=TNAME_LEN) :: test_name   
   complex(dpc),pointer :: results(:)  SET2NULL
 end type FFT_prof_t
!!***

 public :: init_FFT_prof
 public :: nullify_FFT_prof
 public :: destroy_FFT_prof
 public :: print_FFT_profs

!----------------------------------------------------------------------

 public :: fftprof_ncalls_per_test
 !
 ! Timing routines.
 public :: time_fourdp
 public :: time_fourdp_cplx
 public :: time_fourwf
 public :: time_rhotwg
 public :: time_padded_fourwf_cplx
 !
 ! Routines for benchmarks.
 public :: prof_fourdp
 public :: prof_fourwf
 public :: prof_rhotwg

!----------------------------------------------------------------------
 ! Number of calls of each FFT algo, used to have a betters statistics for timing.
 integer,save,private :: NCALLS_FOR_TEST=10

CONTAINS  !====================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/init_FFT_test
!! NAME
!!  init_FFT_test
!!
!! FUNCTION
!!  Creation method for the FFT_test_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof,m_fft_prof
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine init_FFT_test(Ftest,fft_setup,kpoint,ecut,boxcutmin,rprimd,nsym,symrel,MPI_enreg_in)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_FFT_test'
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: nsym
 real(dp),intent(in) :: ecut,boxcutmin
 type(FFT_test_t),intent(inout) :: Ftest
 type(MPI_type),intent(inout) :: MPI_enreg_in
!arrays
 integer,intent(in) :: fft_setup(5),symrel(3,3,nsym)
 real(dp),intent(in) :: kpoint(3),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer,parameter :: option_lob=2
 integer :: fftalg,fftcache,ndat
 real(dp) :: ucvol
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 logical,allocatable :: mask(:)

! *************************************************************************

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 !@FFT_test_t
 Ftest%rprimd = rprimd
 Ftest%rmet   = rmet
 Ftest%gprimd = gprimd
 Ftest%gmet   = gmet
 Ftest%ecut   = ecut

 fftalg   = fft_setup(1)
 fftcache = fft_setup(2)
 ndat     = fft_setup(3)
 ABI_CHECK(ndat==1,"ndat/=1 not supported in fftprof")

 Ftest%nthreads  = fft_setup(4)
 Ftest%available = fft_setup(5)

 Ftest%paral_kgb = 0

 Ftest%kpoint     = kpoint
 Ftest%ndat       = ndat

 Ftest%istwf_k = set_istwfk(kpoint)
                                                             
 call get_kg(Ftest%kpoint,Ftest%istwf_k,ecut,gmet,Ftest%npw_k,Ftest%kg_k) 

 call get_kg(Ftest%kpoint,Ftest%istwf_k,ecut,gmet,Ftest%npw_kout,Ftest%kg_kout) 

 call copy_mpi_enreg(MPI_enreg_in,Ftest%MPI_enreg,opt_bandfft=1)

 Ftest%ngfft(7) = fftalg
 Ftest%ngfft(8) = fftcache

 ! Fill part of ngfft
 call getng(boxcutmin,ecut,gmet,Ftest%MPI_enreg%me_fft,Ftest%mgfft,Ftest%nfft,Ftest%ngfft,Ftest%MPI_enreg%nproc_fft,nsym,&
&  option_lob,Ftest%MPI_enreg%paral_fft,symrel)

 ! Compute the index of each plane wave in the FFT grid.
 ABI_ALLOCATE(Ftest%indpw_k,(Ftest%npw_k))
                                                                                     
 ABI_ALLOCATE(mask,(Ftest%npw_k))
 call kgindex(Ftest%indpw_k,Ftest%kg_k,mask,Ftest%MPI_enreg,Ftest%ngfft,Ftest%npw_k)
 ABI_CHECK(ALL(mask),"FFT parallelism not supported in fftprof")
 ABI_DEALLOCATE(mask)

end subroutine init_FFT_test
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/nullify_fft_test
!! NAME
!!  nullify_fft_test
!!
!! FUNCTION
!!  Nullify all pointers.
!!
!! INPUTS
!!
!! PARENTS
!!      fftprof,m_fft_prof
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine nullify_fft_test(Ftest)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_fft_test'
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 type(FFT_test_t),intent(inout) :: Ftest
!arrays

! *************************************************************************

 ! @FFT_test_t
 nullify(Ftest%kg_k)
 nullify(Ftest%kg_kout)
 nullify(Ftest%indpw_k)

 call nullify_mpi_enreg(Ftest%MPI_enreg)

end subroutine nullify_fft_test
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/destroy_FFT_test_0D
!! NAME
!!  destroy_FFT_test_0D
!!
!! FUNCTION
!!  Destruction method for the FFT_test_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_fft_prof
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine destroy_FFT_test_0D(Ftest)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_FFT_test_0D'
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 type(FFT_test_t),intent(inout) :: Ftest

!Local variables-------------------------------
!scalars

! ********************************************************************* 

 !@FFT_test_t
 Ftest%available=0
 Ftest%istwf_k=-1
 Ftest%mgfft=-1
 Ftest%ndat=-1
 Ftest%nfft=-1
 Ftest%nthreads=1
 Ftest%npw_k=-1
 Ftest%npw_kout=-1
 Ftest%paral_kgb=-1
                                         
 Ftest%ecut=zero
                                         
 Ftest%ngfft=-1
                                         
 Ftest%kpoint =zero
 Ftest%rprimd =zero
 Ftest%rmet   =zero
 Ftest%gprimd =zero
 Ftest%gmet   =zero

 if (associated(Ftest%indpw_k))  then
   ABI_DEALLOCATE(Ftest%indpw_k)
 end if
 if (associated(Ftest%kg_k))     then
   ABI_DEALLOCATE(Ftest%kg_k)
 end if
 if (associated(Ftest%kg_kout))  then
   ABI_DEALLOCATE(Ftest%kg_kout)
 end if

 call destroy_mpi_enreg(Ftest%MPI_enreg)

end subroutine destroy_FFT_test_0D
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/destroy_FFT_test_1D
!! NAME
!!  destroy_FFT_test_1D
!!
!! FUNCTION
!!  Destruction method for the FFT_test_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine destroy_FFT_test_1D(Ftest)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_FFT_test_1D'
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 type(FFT_test_t),intent(inout) :: Ftest(:)

!Local variables-------------------------------
!scalars
 integer :: ii

! ********************************************************************* 

 do ii=LBOUND(Ftest,DIM=1),UBOUND(Ftest,DIM=1)
   call destroy_FFT_test_0D(Ftest(ii))
 end do

end subroutine destroy_FFT_test_1D
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/print_FFT_test
!! NAME
!!  print_FFT_test
!!
!! FUNCTION
!!  Printout of the FFT_test_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine print_FFT_test(Ftest,header,unit,mode_paral,prtvol) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_FFT_test'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral 
 character(len=*),optional,intent(in) :: header
 type(FFT_test_t),intent(in) :: Ftest

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg      
! ********************************************************************* 

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol 
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the FFT test object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 !TODO add additional info
 call print_ngfft(Ftest%ngfft,header="ngfft content",unit=std_out,mode_paral="COLL")

end subroutine print_FFT_test
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/name_of
!! NAME
!!  name_of
!!
!! FUNCTION
!!  Returns a string with info on the test.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function name_of(Ftest) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'name_of'
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 character(len=TNAME_LEN) :: name_of
 type(FFT_test_t),intent(in) :: Ftest

!Local variables-------------------------------
!scalars
 character(len=TNAME_LEN) :: library_name,cplex_mode,padding_mode

! ********************************************************************* 

 call fftalg_info(Ftest%ngfft(7),library_name,cplex_mode,padding_mode)
 !name_of = TRIM(library_name)//"; "//TRIM(cplex_mode)//"; "//TRIM(padding_mode)

 write(name_of,'(i3)')Ftest%ngfft(7)
 name_of = TRIM(library_name)//" ("//TRIM(name_of)//")"

end function name_of
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/init_FFT_prof
!! NAME
!!  init_FFT_prof
!!
!! FUNCTION
!!  Creation method for the FFT_prof_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_fft_prof
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine init_FFT_prof(Ftprof,test_name,nthreads,ncalls,cpu_time,wall_time,results)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_FFT_prof'
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: ncalls,nthreads
 real(dp),intent(in) :: cpu_time,wall_time
 character(len=*),intent(in) :: test_name
 type(FFT_prof_t),intent(out) :: Ftprof
!arrays
 complex(dpc),optional,intent(in) :: results(:)

!Local variables-------------------------------
!scalars

! *************************************************************************

 !@FFT_prof_t
 Ftprof%ncalls    =  ncalls
 Ftprof%nthreads  =  nthreads
 Ftprof%cpu_time  = cpu_time
 Ftprof%wall_time = wall_time
 Ftprof%test_name = test_name

 !nullify(Ftprof%results)
 if (PRESENT(results)) then
   if (associated(Ftprof%results))  then
     ABI_DEALLOCATE(Ftprof%results)
   end if
   ABI_ALLOCATE(Ftprof%results,(SIZE(results)))
   Ftprof%results = results
 end if

end subroutine init_FFT_prof
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/nullify_FFT_prof
!! NAME
!!  nullify_FFT_prof
!!
!! FUNCTION
!!  Nullify all pointers.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine nullify_FFT_prof(Ftprof)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_FFT_prof'
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 type(FFT_prof_t),intent(inout) :: Ftprof
!arrays

! *************************************************************************

 !@FFT_prof_t
 nullify(Ftprof%results)

end subroutine nullify_FFT_prof
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/destroy_FFT_prof_0D
!! NAME
!!  destroy_FFT_prof_0D
!!
!! FUNCTION
!!  Destruction method for the FFT_prof_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_fft_prof
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine destroy_FFT_prof_0D(Ftprof)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_FFT_prof_0D'
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 type(FFT_prof_t),intent(inout) :: Ftprof

! ********************************************************************* 

 !@FFT_prof_t
 Ftprof%ncalls=0
 Ftprof%nthreads=0
 Ftprof%cpu_time=zero
 Ftprof%wall_time=zero
 Ftprof%test_name = "None"
 if (associated(Ftprof%results))  then
   ABI_DEALLOCATE(Ftprof%results)
 end if

end subroutine destroy_FFT_prof_0D
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/destroy_FFT_prof_1D
!! NAME
!!  destroy_FFT_prof_1D
!!
!! FUNCTION
!!  Destruction method for the FFT_prof_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine destroy_FFT_prof_1D(Ftprof)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_FFT_prof_1D'
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 type(FFT_prof_t),intent(inout) :: Ftprof(:)

!Local variables-------------------------------
!scalars
 integer :: ii
! ********************************************************************* 

 !@FFT_prof_t
 do ii=LBOUND(Ftprof,DIM=1),UBOUND(Ftprof,DIM=1)
   call destroy_FFT_prof_0D(Ftprof(ii))
 end do

end subroutine destroy_FFT_prof_1D
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/print_FFT_profs
!! NAME
!!  print_FFT_profs
!!
!! FUNCTION
!!  Printout of the FFT_prof_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine print_FFT_profs(Fprof,header,unit,mode_paral,prtvol) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_FFT_profs'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral 
 character(len=*),optional,intent(in) :: header
 type(FFT_prof_t),intent(in) :: Fprof(:)

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol,ncalls
 integer :: field1_w,ii,ref_lib !ifft
 real(dp) :: mabs_err,mean_err,check_mabs_err,check_mean_err
 character(len=4) :: my_mode
 character(len=500) :: ofmt,hfmt,nafmt
 character(len=500) :: msg      
! ********************************************************************* 

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol 
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg='==== Info on the FFT_prof_t object ===='
 if (PRESENT(header)) msg='==== '//TRIM(ADJUSTL(header))//' ===='

 call wrtout(my_unt,ch10//REPEAT("=",LEN_TRIM(msg)))
 call wrtout(my_unt,msg,my_mode)
 call wrtout(my_unt,REPEAT("=",LEN_TRIM(msg)))

 field1_w=0 ! Width of the field used to print key names.
 do ii=1,SIZE(Fprof)
   field1_w = MAX(field1_w, LEN_TRIM(Fprof(ii)%test_name))
 end do

 if (field1_w==0) RETURN
 field1_w = field1_w + 2 ! To account for ". "
 write(ofmt,*)"(a",field1_w,",2x,2(f7.4,4x),5x,i0,7x,i0,4x,2(es9.2,3x))"

 write(hfmt,*)"(a",field1_w,",2x,a)"
 write(std_out,hfmt)" Library      ","CPU-time   WALL-time   nthreads  ncalls  Max_|Err|   <|Err|>"
 !
 ! Find reference library.
 ref_lib=0
 do ii=1,SIZE(Fprof)
   if (Fprof(ii)%ncalls>0) then
     ref_lib = ii
     EXIT
   end if
 end do
 !ref_lib=3
 !
 ! Write timing analysis and error wrt reference library if available.
 check_mabs_err=zero; check_mean_err=zero
 do ii=1,SIZE(Fprof)
   ncalls = Fprof(ii)%ncalls
   if (ncalls>0) then
     mabs_err = zero; mean_err=zero
     if (ref_lib>0) then
       mabs_err = MAXVAL( ABS(Fprof(ii)%results - Fprof(ref_lib)%results) )
       mean_err = SUM( ABS(Fprof(ii)%results - Fprof(ref_lib)%results) ) / SIZE(Fprof(ref_lib)%results)
       ! Relative error is not a good estimator because some components are close to zero within machine accuracy.
       !mean_err = 100 * MAXVAL( ABS(Fprof(ii)%results - Fprof(1)%results)/ ABS(Fprof(1)%results ))
       !ifft = imax_loc( ABS(Fprof(ii)%results - Fprof(1)%results)/ ABS(Fprof(1)%results) )
       !write(std_out,*) Fprof(ii)%results(ifft),Fprof(1)%results(ifft)
     end if 
     write(std_out,ofmt)&
&      "- "//Fprof(ii)%test_name,Fprof(ii)%cpu_time/ncalls,Fprof(ii)%wall_time/ncalls,Fprof(ii)%nthreads,ncalls,mabs_err,mean_err
     check_mabs_err = MAX(check_mabs_err, mabs_err)
     check_mean_err = MAX(check_mean_err, mean_err)
   else 
     write(nafmt,*)"(a",field1_w,",2x,a)"
     write(std_out,nafmt)"- "//Fprof(ii)%test_name,"   N/A        N/A        N/A       N/A         N/A"
   end if
 end do

 if (ref_lib>0) then
   write(std_out,'(/,2(a,es9.2),2a)')&
&    " Consistency check: MAX(Max_|Err|) = ",check_mabs_err,&
&    ", Max(<|Err|>) = ",check_mean_err,", reference_lib: ",TRIM(Fprof(ref_lib)%test_name) 
 end if
 write(std_out,*)

end subroutine print_FFT_profs
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/time_fourdp
!! NAME
!!  time_fourdp
!!
!! FUNCTION
!!  Profiling of the fourdp routine.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof,m_fft_prof
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine time_fourdp(Ftest,isign,cplex,Ftprof)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'time_fourdp'
 use interfaces_18_timing
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: cplex,isign
 type(FFT_test_t),intent(inout) :: Ftest 
 type(FFT_prof_t),intent(out) :: Ftprof

!Local variables-------------------------------
!scalars
 integer :: icall,i1,i2,i3,n1,n2,n3,ifft
 real(dp) :: cpu_time,wall_time,cpu0,wall0
 real(dp) :: gsq
 character(len=500) :: msg
 character(len=TNAME_LEN) :: test_name
!arrays
 integer :: gg(3)
 integer,parameter :: g0(3)=(/1,2,-1/)
 real(dp),allocatable :: fofg(:,:),fofr(:)
 complex(dpc),allocatable :: results(:)

! ********************************************************************* 

 test_name = name_of(Ftest)
 n1=Ftest%ngfft(1); n2=Ftest%ngfft(2); n3=Ftest%ngfft(3)

 if (Ftest%available==0) then
   call init_FFT_prof(Ftprof,test_name,0,0,zero,zero)
   RETURN
 end if

 ABI_ALLOCATE(fofg,(2,Ftest%nfft))
 ABI_ALLOCATE(fofr,(cplex*Ftest%nfft))
 !
 ! Initialize input data.
 if (isign==1) then ! initialize fofg
   fofg = zero
   ifft=0
   do i3=1,n3
     gg(3)=i3-1; if (i3>1+n3/2) gg(3)=i3-n3-1 ! TODO recheck this
     do i2=1,n2
       gg(2)=i2-1; if (i2>1+n2/2) gg(2)=i2-n2-1
       do i1=1,n1
         gg(1)=i1-1; if (i1>1+n1/2) gg(1)=i1-n1-1
         gsq = two_pi**2 * DOT_PRODUCT(gg,MATMUL(Ftest%gmet,gg))
         ifft=ifft+1
         fofg(1,ifft) = EXP(-gsq)
         fofg(2,ifft) = zero
       end do
     end do
   end do
 
 else ! init fofr 
   if (cplex==2) then
     fofr = eigr(g0,Ftest%nfft,Ftest%ngfft)
   else if (cplex==1) then  
     fofr = REAL(ceigr(g0,Ftest%nfft,Ftest%ngfft))
   else 
     write(msg,'(a,i0)')" Wrong cplex: ",cplex
     MSG_ERROR(msg)
   end if
 end if

 ABI_ALLOCATE(results,(Ftest%nfft))
 results=czero

#ifdef HAVE_OPENMP
 call omp_set_num_threads(Ftest%nthreads)
#endif
 call fftw3_set_nthreads(Ftest%nthreads)

 call timein(cpu0,wall0)

 do icall=1,NCALLS_FOR_TEST
   call fourdp(cplex,fofg,fofr,isign,Ftest%MPI_enreg,Ftest%nfft,Ftest%ngfft,Ftest%paral_kgb,0)
   !
   ! Store results at the first call.
   if (icall==1) then 
     if (isign==1) then
       if (cplex==1) then
         results = CMPLX(fofr, zero)
       else if (cplex==2) then
         results = CMPLX( fofr(1:2*Ftest%nfft:2), fofr(2:2*Ftest%nfft:2) )
       end if
     else if (isign==-1) then
       results = CMPLX(fofg(1,:),fofg(2,:))
     end if
   end if
 end do

 call timein(cpu_time,wall_time)
 cpu_time  = cpu_time-cpu0
 wall_time = wall_time-wall0

 ABI_DEALLOCATE(fofg)
 ABI_DEALLOCATE(fofr)

 call init_FFT_prof(Ftprof,test_name,Ftest%nthreads,NCALLS_FOR_TEST,cpu_time,wall_time,results=results)

 ABI_DEALLOCATE(results)

end subroutine time_fourdp
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/time_fourdp_cplx
!! NAME
!!  time_fourdp_cplx
!!
!! FUNCTION
!!  Profiling of the fourdp_c2c_* routines.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine time_fourdp_cplx(Ftest,isign,inplace,Ftprof)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'time_fourdp_cplx'
 use interfaces_18_timing
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: isign,inplace
 type(FFT_test_t),intent(inout) :: Ftest
 type(FFT_prof_t),intent(out) :: Ftprof

!Local variables-------------------------------
!scalars
 integer :: icall,i1,i2,i3,n1,n2,n3,ifft
 real(dp) :: cpu_time,wall_time,cpu0,wall0,gsq
 character(len=500) :: msg
 character(len=TNAME_LEN) :: test_name
!arrays
 integer,parameter :: g0(3) = (/1,-2,1/)
 integer :: gg(3)
 complex(dpc),allocatable :: ffc(:),ggc(:),results(:)
! ********************************************************************* 

 test_name = name_of(Ftest)

 if (Ftest%available==0) then
   call init_FFT_prof(Ftprof,test_name,0,0,zero,zero)
   RETURN
 end if

 n1=Ftest%ngfft(1); n2=Ftest%ngfft(2); n3=Ftest%ngfft(3)

 ABI_ALLOCATE(ffc,(Ftest%nfft))
 ABI_ALLOCATE(ggc,(Ftest%nfft))
 ABI_ALLOCATE(results,(Ftest%nfft))
 ffc=czero; ggc=czero; results=czero

 if (isign==-1) then
   ffc = ceigr(g0,Ftest%nfft,Ftest%ngfft)
 else if (isign==1) then
   ifft=0
   do i3=1,n3
     gg(3)=i3-1; if (i3>1+n3/2) gg(3)=i3-n3-1 ! TODO recheck this
     do i2=1,n2
       gg(2)=i2-1; if (i2>1+n2/2) gg(2)=i2-n2-1
       do i1=1,n1
         gg(1)=i1-1; if (i1>1+n1/2) gg(1)=i1-n1-1
         gsq = two_pi**2 * DOT_PRODUCT(gg,MATMUL(Ftest%gmet,gg))
         ifft=ifft+1
         ffc(ifft) = EXP(-gsq)
       end do
     end do
   end do

 else 
   write(msg,'(a,i0)')" Wrong isign= ",isign
   MSG_ERROR(msg)
 end if

#ifdef HAVE_OPENMP
 call omp_set_num_threads(Ftest%nthreads)
#endif
 call fftw3_set_nthreads(Ftest%nthreads)

 call timein(cpu0,wall0)

 select case (inplace)
 case (0)
   do icall=1,NCALLS_FOR_TEST
     call fourdp_c2c_op(ffc,ggc,isign,Ftest%mpi_enreg,Ftest%nfft,Ftest%ngfft,Ftest%paral_kgb,0)
     ! Store results at the first call.
     if (icall==1) results = ggc
   end do 
 case (1)
   do icall=1,NCALLS_FOR_TEST
     call fourdp_c2c_ip(ffc,isign,Ftest%mpi_enreg,Ftest%nfft,Ftest%ngfft,Ftest%paral_kgb,0)
     ! Store results at the first call.
     if (icall==1) results = ffc
   end do 
 case default
  write(msg,'(a,i0)')" Wrong value for inplace= ",inplace
  MSG_ERROR(msg)
 end select

 call timein(cpu_time,wall_time)
 cpu_time  = cpu_time-cpu0
 wall_time = wall_time-wall0

 ABI_DEALLOCATE(ffc)
 ABI_DEALLOCATE(ggc)

 call init_FFT_prof(Ftprof,test_name,Ftest%nthreads,NCALLS_FOR_TEST,cpu_time,wall_time,results=results)

 ABI_DEALLOCATE(results)

end subroutine time_fourdp_cplx
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/time_fourwf
!! NAME
!!  time_fourwf
!!
!! FUNCTION
!!  Profiling of the fourwf routine.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof,m_fft_prof
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine time_fourwf(Ftest,cplex,option_fourwf,Ftprof)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'time_fourwf'
 use interfaces_18_timing
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: cplex,option_fourwf
 type(FFT_test_t),intent(inout) :: Ftest
 type(FFT_prof_t),intent(out) :: Ftprof

!Local variables-------------------------------
!scalars
 integer :: n1,n2,n3,n4,n5,n6,npw_out,icall,i1,i2,i3,idx,ipw
 !integer,save :: sunit=100
 real(dp),parameter :: weight_i=one,weight_r=one
 real(dp) :: cpu_time,wall_time,cpu0,wall0,gsq,g0dotr
 logical :: isbuggy
 character(len=500) :: msg
 character(len=TNAME_LEN) :: test_name
!arrays
 integer,parameter :: g0(3)=(/1,-1,2/)
 !integer,parameter :: g0(3)=(/1,0,0/)
 integer :: gg(3)
 integer,allocatable :: gbound_in(:,:),gbound_out(:,:) 
 real(dp),allocatable :: denpot(:,:,:),fofg_in(:,:) 
 real(dp),allocatable :: fofr_4(:,:,:,:),fofg_out(:,:) !,fofr(:)
 complex(dpc),allocatable :: results(:)

! ********************************************************************* 

 test_name = name_of(Ftest)

 ! FIXME option_fourwf==3 makes SG2001 crash!
 ! FIXME SG2001 produces wrong results if option_fourwf==2 and cplex==2 
 isbuggy = (option_fourwf==3 .and. Ftest%ngfft(7)==412) .or. &
&  (option_fourwf==2 .and. cplex==2 .and. Ftest%ngfft(7)==412)

 if (isbuggy.or.Ftest%available==0) then
   call init_FFT_prof(Ftprof,test_name,0,0,zero,zero)
   RETURN
 end if

 npw_out=Ftest%npw_kout

 n1=Ftest%ngfft(1); n2=Ftest%ngfft(2); n3=Ftest%ngfft(3)
 n4=Ftest%ngfft(4); n5=Ftest%ngfft(5); n6=Ftest%ngfft(6)

 ! FFTW3 does not need gbound_in but oh well
 ABI_ALLOCATE(gbound_in,(2*Ftest%mgfft+8,2))
 call sphereboundary(gbound_in,Ftest%istwf_k,Ftest%kg_k,Ftest%mgfft,Ftest%npw_k)

 ABI_ALLOCATE(gbound_out,(2*Ftest%mgfft+8,2))
 call sphereboundary(gbound_out,Ftest%istwf_k,Ftest%kg_kout,Ftest%mgfft,Ftest%npw_kout)

 ABI_ALLOCATE(denpot,(cplex*n4,n5,n6))
 ABI_ALLOCATE(fofg_in,(2,Ftest%npw_k*Ftest%ndat))
 ABI_ALLOCATE(fofg_out,(2,npw_out*Ftest%ndat))
 ABI_ALLOCATE(fofr_4,(2,n4,n5,n6*Ftest%ndat))
 denpot   = zero
 fofg_in  = zero
 fofg_out = zero
 fofr_4   = zero

 ABI_ALLOCATE(results,(Ftest%nfft))
 results=czero

 select case (option_fourwf)
 case (0,1,2)
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
   !!
   do ipw=1,Ftest%npw_k
     gg = Ftest%kg_k(:,ipw)
     gsq = two_pi**2 * DOT_PRODUCT(gg,MATMUL(Ftest%gmet,gg))
     fofg_in(1,ipw) = EXP(-gsq)
     fofg_in(2,ipw) = zero
   end do
   !
   if (option_fourwf==2) then ! Init denpot
     !
     if (cplex==1) then
       do i3=0,n3-1
         do i2=0,n2-1
           do i1=0,n1-1
             g0dotr= two_pi*( g0(1)*(i1/DBLE(n1)) &
&                            +g0(2)*(i2/DBLE(n2)) &
&                            +g0(3)*(i3/DBLE(n3)) )
             denpot(i1+1,i2+1,i3+1)=COS(g0dotr)
           end do
         end do
       end do
     else if (cplex==2) then
       do i3=0,n3-1
         do i2=0,n2-1
           idx=1
           do i1=0,n1-1
             g0dotr= two_pi*( g0(1)*(i1/DBLE(n1)) &
&                            +g0(2)*(i2/DBLE(n2)) &
&                            +g0(3)*(i3/DBLE(n3)) )

             denpot(idx,  i2+1,i3+1)= COS(g0dotr)
             denpot(idx+1,i2+1,i3+1)= SIN(g0dotr)
             idx=idx+2
           end do
         end do
       end do
     end if
   end if

 case (3)
   !! for option==3, fofr(2,n4,n5,n6*ndat) contains the input real space wavefunction;
   !!                fofgout(2,npwout*ndat) contains its output Fourier transform;
   !!                no use of fofgin and npwin.
   idx=0
   do i3=0,n3-1
     do i2=0,n2-1
       do i1=0,n1-1
         g0dotr= two_pi*( g0(1)*(i1/DBLE(n1)) &
                         +g0(2)*(i2/DBLE(n2)) &
                         +g0(3)*(i3/DBLE(n3)) )
         !fofr_4(1,i1+1,i2+1,i3+1)=COS(g0dotr)
         !fofr_4(2,i1+1,i2+1,i3+1)=SIN(g0dotr)
         fofr_4(1,i1+1,i2+1,i3+1)=EXP(-g0dotr**2)
         fofr_4(2,i1+1,i2+1,i3+1)=zero
         idx=idx+1
         results(idx) = CMPLX( fofr_4(1,i1+1,i2+1,i3+1), fofr_4(2,i1+1,i2+1,i3+1))
       end do
     end do
   end do

 case default
   write(msg,'(a,i0)')" Wrong value for option_fourwf: ",option_fourwf
   MSG_ERROR(msg)
 end select 

#ifdef HAVE_OPENMP
 call omp_set_num_threads(Ftest%nthreads)
#endif
 call fftw3_set_nthreads(Ftest%nthreads)

 call timein(cpu0,wall0)

 do icall=1,NCALLS_FOR_TEST

   call fourwf(cplex,denpot,fofg_in,fofg_out,fofr_4,gbound_in,gbound_out,Ftest%istwf_k,&
&    Ftest%kg_k,Ftest%kg_kout,Ftest%mgfft,Ftest%MPI_enreg,Ftest%ndat,Ftest%ngfft,Ftest%npw_k,npw_out,n4,n5,n6,option_fourwf,&
&    Ftest%paral_kgb,0,weight_r,weight_i)
   !
   ! Store results at the first call.
   if (icall==1) then 
     select case (option_fourwf)
     case (0)
       !! for option==0, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
       !!                fofr(2,n4,n5,n6) contains the output Fourier Transform of fofgin;
       !!                no use of denpot, fofgout and npwout.
       idx=0
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             idx=idx+1
             results(idx) = CMPLX(fofr_4(1,i1,i2,i3),fofr_4(2,i1,i2,i3))
           end do
         end do
       end do
     case (1)
       !! for option==1, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
       !!                denpot(cplex*n4,n5,n6) contains the input density at input,
       !!                and the updated density at output (accumulated);
       !!                no use of fofgout and npwout.
       if (cplex==1) then
        idx=0
        do i3=1,n3
          do i2=1,n2
            do i1=1,n1
              idx=idx+1
              results(idx) = CMPLX(denpot(i1,i2,i3),zero)
            end do
          end do
        end do

       else if (cplex==2) then
         idx=0
         do i3=1,n3
           do i2=1,n2
             do i1=1,2*n1,2
               idx=idx+1
               results(idx) = CMPLX(denpot(i1,i2,i3),denpot(i1+1,i2,i3))
             end do
           end do
         end do
       else 
         MSG_ERROR("Wrong cplex")
       end if

     case (2)
       !! for option==2, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
       !!                denpot(cplex*n4,n5,n6) contains the input local potential;
       !!                fofgout(2,npwout*ndat) contains the output function;
       do ipw=1,npw_out
         results(ipw) = CMPLX(fofg_out(1,ipw),fofg_out(2,ipw))
       end do
     case (3)
       !! for option==3, fofr(2,n4,n5,n6*ndat) contains the input real space wavefunction;
       !!                fofgout(2,npwout*ndat) contains its output Fourier transform;
       !!                no use of fofgin and npwin.
       results=czero
       do ipw=1,npw_out
         results(ipw) = CMPLX(fofg_out(1,ipw),fofg_out(2,ipw))
       end do
       !fofr_4=zero
     end select 
   end if
 end do

 call timein(cpu_time,wall_time)
 cpu_time  = cpu_time-cpu0
 wall_time = wall_time-wall0

 ABI_DEALLOCATE(denpot)
 ABI_DEALLOCATE(fofg_in)
 ABI_DEALLOCATE(fofg_out)
 ABI_DEALLOCATE(fofr_4)
 ABI_DEALLOCATE(gbound_in)
 ABI_DEALLOCATE(gbound_out)

 call init_FFT_prof(Ftprof,test_name,Ftest%nthreads,NCALLS_FOR_TEST,cpu_time,wall_time,results=results)

!BEGINDEBUG
 !if (option_fourwf==2.and.cplex==1) then
 !if (option_fourwf==3) then
 !  sunit=sunit+1
 !  write(sunit,*)"option_fourwf",option_fourwf," ngfft ",Ftest%ngfft
 !  do idx=1,SIZE(results)
 !    write(sunit,*)results(idx)
 !  end do
 !end if
!ENDDEBUG

 ABI_DEALLOCATE(results)

end subroutine time_fourwf
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/fftprof_ncalls_per_test
!! NAME
!!  fftprof_ncalls_per_test
!!
!! FUNCTION
!!  Helper function used to set the number of calls to  be used in each time_* routine.
!!
!! INPUTS
!!   nc=Number of calls to be used.
!!
!! SIDE EFFECTS 
!!  NCALLS_FOR_TEST = nc
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine fftprof_ncalls_per_test(ncalls)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftprof_ncalls_per_test'
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: ncalls

! ********************************************************************* 

 NCALLS_FOR_TEST = ncalls

end subroutine fftprof_ncalls_per_test
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/time_rhotwg
!! NAME
!!  time_rhotwg
!!
!! FUNCTION
!!  Profiling of the rho_tw_g  routine.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof,m_fft_prof
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine time_rhotwg(Ftest,map2sphere,use_padfft,osc_npw,osc_gvec,Ftprof)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'time_rhotwg'
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: map2sphere,use_padfft,osc_npw
 type(FFT_test_t),intent(inout) :: Ftest
 type(FFT_prof_t),intent(out) :: Ftprof
!arrays 
 integer,intent(in) :: osc_gvec(3,osc_npw)

!Local variables-------------------------------
!scalars
 integer,parameter :: paral_kgb=0,nspinor=1,dim_rtwg=1,istwfk1=1
 !integer,save :: rtw_unt=300
 integer :: icall,istat,ifft,itim1,itim2,nfft,ipw
 integer :: n1,n2,n3,n4,n5,n6,cplex,i1,i2,i3
 real(dp) :: cpu_time,wall_time,cpu0,wall0
 complex(dpc) :: ktabp1=cone,ktabp2=cone
 character(len=TNAME_LEN) :: test_name
 logical :: not_implemented
 type(MPI_type) :: MPI_enreg_seq
!arrays
 !integer,parameter :: g1(3)=(/1,1,2/),g2(3)=(/-1,-2,-1/)
 integer,parameter :: g1(3)=(/-1,0,0/),g2(3)=(/1,0,0/)
 integer,allocatable :: gbound(:,:) 
 integer,allocatable :: ktabr1(:),ktabr2(:)
 integer,allocatable :: igfftg0(:)
 real(dp),parameter :: spinrot1(4)=(/one,zero,zero,one/),spinrot2(4)=(/one,zero,zero,one/)
 real(dp),allocatable :: denpot(:,:,:),fofg_in(:,:) 
 real(dp),allocatable :: fofr_4(:,:,:,:),fofg_out(:,:) 
 logical,allocatable :: mask(:)
 complex(dpc),allocatable :: results(:)
 complex(gwpc),allocatable :: rhotwg(:)
 complex(gwpc),allocatable :: wfn1(:),wfn2(:)
 complex(dpc),allocatable :: uu(:),usk(:)

! ********************************************************************* 

 test_name = name_of(Ftest)

 nfft = Ftest%nfft
 n1=Ftest%ngfft(1); n2=Ftest%ngfft(2); n3=Ftest%ngfft(3)
 n4=Ftest%ngfft(4); n5=Ftest%ngfft(5); n6=Ftest%ngfft(6)

 ! TODO: zero-pad not available with SG2001 routines.
 not_implemented = (use_padfft==1.and.Ftest%ngfft(7) == 412)

 if (Ftest%available==0.or.not_implemented) then
   call init_FFT_prof(Ftprof,test_name,0,0,zero,zero)
   RETURN
 end if

 call initmpi_seq(MPI_enreg_seq)

 itim1=1; itim2=1
 ABI_ALLOCATE(ktabr1,(nfft))
 ABI_ALLOCATE(ktabr2,(nfft))

 do ifft=1,nfft
   ktabr1(ifft)= ifft 
   ktabr2(ifft)= ifft 
 end do

 ABI_ALLOCATE(igfftg0,(osc_npw*map2sphere))
 ABI_ALLOCATE(mask,(osc_npw))

 call kgindex(igfftg0,osc_gvec,mask,MPI_enreg_seq,Ftest%ngfft,osc_npw)
 ABI_CHECK(ALL(mask)," FFT parallelism not supported")
 ABI_DEALLOCATE(mask)

 ABI_ALLOCATE(gbound,(2*Ftest%mgfft+8,2*use_padfft))
 if (use_padfft==1) then  
   call sphereboundary(gbound,istwfk1,osc_gvec,Ftest%mgfft,osc_npw)
 end if

 ABI_ALLOCATE(wfn1,(nfft*nspinor))
 ABI_ALLOCATE(wfn2,(nfft*nspinor))

 wfn1 = ceigr(g1,nfft,Ftest%ngfft)
 wfn2 = ceigr(g2,nfft,Ftest%ngfft)

 ABI_ALLOCATE(rhotwg,(osc_npw*dim_rtwg))
 ABI_ALLOCATE(results,(osc_npw*dim_rtwg))

#ifdef HAVE_OPENMP
 call omp_set_num_threads(Ftest%nthreads)
#endif
 call fftw3_set_nthreads(Ftest%nthreads)

 call timein(cpu0,wall0)

 do icall=1,NCALLS_FOR_TEST

   if (.FALSE. .and. use_padfft==1 .and. Ftest%ngfft(7)/=312) then
   !if (use_padfft==1 .and. Ftest%ngfft(7)/=312) then
     MSG_WARNING("Calling fourwf instead of rho_tw_g")
     ABI_ALLOCATE(uu,(nfft))
     ABI_ALLOCATE(usk,(nfft))
     uu  = wfn1(ktabr1)*ktabp1; if (itim1==1) uu  = CONJG(uu)
     usk = wfn2(ktabr2)*ktabp2; if (itim2==2) usk = CONJG(usk)
     uu  = uu * usk

     cplex=1
     ABI_ALLOCATE(denpot,(cplex*n4,n5,n6))
     ABI_ALLOCATE(fofg_in,(2,Ftest%npw_k*Ftest%ndat))
     ABI_ALLOCATE(fofg_out,(2,osc_npw*Ftest%ndat))
     ABI_ALLOCATE(fofr_4,(2,n4,n5,n6*Ftest%ndat))
     denpot   = zero
     fofg_in  = zero
     fofg_out = zero
     fofr_4   = zero

     ifft=0
     do i3=1,n3
       do i2=1,n2
         do i1=1,n1
           ifft=ifft+1
           fofr_4(1,i1,i2,i3)= REAL(uu(ifft))
           fofr_4(2,i1,i2,i3)= AIMAG(uu(ifft))
         end do
       end do
     end do

     ABI_DEALLOCATE(uu)
     ABI_DEALLOCATE(usk)

     call fourwf(cplex,denpot,fofg_in,fofg_out,fofr_4,gbound,gbound,istwfk1,&
&      Ftest%kg_k,osc_gvec,Ftest%mgfft,Ftest%MPI_enreg,Ftest%ndat,Ftest%ngfft,Ftest%npw_k,osc_npw,n4,n5,n6,3,&
&      Ftest%paral_kgb,0,one,one)

     ! Store results at the first call.
     if (icall==1) then
       results=czero
       do ipw=1,osc_npw
         results(ipw) = CMPLX(fofg_out(1,ipw),fofg_out(2,ipw))
       end do
     end if

     ABI_DEALLOCATE(denpot)
     ABI_DEALLOCATE(fofg_in)
     ABI_DEALLOCATE(fofg_out)
     ABI_DEALLOCATE(fofr_4)

   else
     call rho_tw_g(paral_kgb,nspinor,osc_npw,nfft,Ftest%ngfft,map2sphere,use_padfft,igfftg0,gbound,&
&      wfn1,itim1,ktabr1,ktabp1,spinrot1,&
&      wfn2,itim2,ktabr2,ktabp2,spinrot2,&
&      dim_rtwg,rhotwg,0,MPI_enreg_seq)
     ! Store results at the first call.
     if (icall==1) results = rhotwg
   end if
 end do

 call timein(cpu_time,wall_time)
 cpu_time  = cpu_time-cpu0
 wall_time = wall_time-wall0

 ABI_DEALLOCATE(ktabr1)
 ABI_DEALLOCATE(ktabr2)
 ABI_DEALLOCATE(igfftg0)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(gbound)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(wfn1)
 ABI_DEALLOCATE(wfn2)

 call init_FFT_prof(Ftprof,test_name,Ftest%nthreads,NCALLS_FOR_TEST,cpu_time,wall_time,results)

!BEGINDEBUG
! rtw_unt=rtw_unt+1
! write(rtw_unt,*)" rho_tw_g with ngfft ",Ftest%ngfft
! do ifft=1,SIZE(results)
!   write(rtw_unt,*)results(ifft)
! end do
!ENDDEBUG

 ABI_DEALLOCATE(results)

end subroutine time_rhotwg
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/time_padded_fourwf_cplx
!! NAME
!!  time_padded_fourwf_cplx
!!
!! FUNCTION
!!  Profiling of the padded_fourwf_cplx routine.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      destroy_fft_prof,destroy_fft_test,get_kg,init_fft_test,nullify_fft_test
!!      time_rhotwg
!!
!! SOURCE

subroutine time_padded_fourwf_cplx(Ftest,isign,Ftprof)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'time_padded_fourwf_cplx'
 use interfaces_18_timing
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: isign
 type(FFT_test_t),intent(inout) :: Ftest
 type(FFT_prof_t),intent(out) :: Ftprof

!Local variables-------------------------------
!scalars
 integer,parameter :: istwfk1=1
 integer :: icall,i1,i2,i3,n1,n2,n3,n4,n5,n6,ifft,npw_k,idx,ipw
 real(dp) :: cpu_time,wall_time,cpu0,wall0,gsq,g0dotr
 logical :: not_implemented
 character(len=500) :: msg
 character(len=TNAME_LEN) :: test_name
!arrays
 integer,parameter :: g0(3) = (/1,-2,1/)
 integer,pointer :: kg_k(:,:)
 integer :: gg(3)
 integer,allocatable :: gbound(:,:)
 complex(dpc),allocatable :: ffc(:),ffc_save(:),results(:),cwork(:)
! ********************************************************************* 

 test_name = name_of(Ftest)

 ! TODO: zero-pad not available with SG2001 routines.
 not_implemented = (Ftest%ngfft(7) == 412)

 if (Ftest%available==0.or.not_implemented) then
   call init_FFT_prof(Ftprof,test_name,0,0,zero,zero)
   RETURN
 end if

 n1=Ftest%ngfft(1); n2=Ftest%ngfft(2); n3=Ftest%ngfft(3)
 n4=Ftest%ngfft(4); n5=Ftest%ngfft(5); n6=Ftest%ngfft(6)

 ABI_ALLOCATE(ffc,(n4*n5*n6))
 ABI_ALLOCATE(ffc_save,(n4*n5*n6))
 ABI_ALLOCATE(results,(Ftest%nfft))
 ffc=czero; ffc_save=czero; results=czero

 if (isign==-1) then
   do i3=0,n3-1
     do i2=0,n2-1
       do i1=0,n1-1
         g0dotr= two_pi*( g0(1)*(i1/DBLE(n1)) &
&                        +g0(2)*(i2/DBLE(n2)) &
&                        +g0(3)*(i3/DBLE(n3)) )
         ifft= 1 + i1 + i2*n4 + i3*n4*n5
         ffc(ifft)=DCMPLX(DCOS(g0dotr),DSIN(g0dotr))
       end do
     end do
   end do

 else if (isign==1) then
   do i3=1,n3
     gg(3)=i3-1; if (i3>1+n3/2) gg(3)=i3-n3-1 ! TODO recheck this
     do i2=1,n2
       gg(2)=i2-1; if (i2>1+n2/2) gg(2)=i2-n2-1
       do i1=1,n1
         gg(1)=i1-1; if (i1>1+n1/2) gg(1)=i1-n1-1
         gsq = two_pi**2 * DOT_PRODUCT(gg,MATMUL(Ftest%gmet,gg))
         ifft= i1 + (i2-1)*n4 + (i3-1)*n4*n5
         if (half*gsq <= Ftest%ecut) then
           ffc(ifft) = EXP(-gsq)
         else 
           ffc(ifft) = czero
         end if
       end do
     end do
   end do

 else 
   write(msg,'(a,i0)')" Wrong isign= ",isign
   MSG_ERROR(msg)
 end if
 ffc_save = ffc
 !
 ! Do not use istwfk triks here
 call get_kg(Ftest%kpoint,istwfk1,Ftest%ecut,Ftest%gmet,npw_k,kg_k) 

 ABI_ALLOCATE(gbound,(2*Ftest%mgfft+8,2))
 call sphereboundary(gbound,istwfk1,kg_k,Ftest%mgfft,npw_k)
 ABI_DEALLOCATE(kg_k)

#ifdef HAVE_OPENMP
 call omp_set_num_threads(Ftest%nthreads)
#endif
 call fftw3_set_nthreads(Ftest%nthreads)

 call timein(cpu0,wall0)

 do icall=1,NCALLS_FOR_TEST
   call padded_fourwf_cplx(ffc,Ftest%ngfft,n1,n2,n3,n4,n5,n6,Ftest%mgfft,isign,gbound)
   ! Store results at the first call.
   if (icall==1) then 
     if (isign==-1) then 
        ABI_ALLOCATE(cwork,(n1*n2*n3))
        ifft=0
        do i3=1,n3
          do i2=1,n2
            do i1=1,n1
              ifft = ifft+1 
              idx = i1 + (i2-1)*n4 + (i3-1)*n4*n5
              cwork(ifft) = ffc(idx)
            end do
          end do
        end do

        do ipw=1,npw_k
          ifft = Ftest%indpw_k(ipw)
          results(ipw) = cwork(ifft)
        end do
        ABI_DEALLOCATE(cwork)
     else if (isign==+1)  then
       results = ffc(1:Ftest%nfft)
     end if
   end if
   ffc = ffc_save ! Restore input to avoid numerical exceptions.
 end do

 call timein(cpu_time,wall_time)
 cpu_time  = cpu_time-cpu0
 wall_time = wall_time-wall0

 ABI_DEALLOCATE(gbound)
 ABI_DEALLOCATE(ffc)
 ABI_DEALLOCATE(ffc_save)

 call init_FFT_prof(Ftprof,test_name,Ftest%nthreads,NCALLS_FOR_TEST,cpu_time,wall_time,results=results)

 ABI_DEALLOCATE(results)

end subroutine time_padded_fourwf_cplx
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/prof_fourdp
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine prof_fourdp(fft_setups,isign,cplex,necut,ecut_arth,boxcutmin,rprimd,nsym,symrel,MPI_enreg_in)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prof_fourdp'
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: nsym,isign,cplex,necut
 real(dp),intent(in) :: boxcutmin
 type(MPI_type),intent(inout) :: MPI_enreg_in
!arrays
 integer,intent(in) :: fft_setups(:,:),symrel(3,3,nsym)
 real(dp),intent(in) :: ecut_arth(2)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: iec,nsetups,set,funt
 type(FFT_test_t) :: Ftest
 type(FFT_prof_t) :: Ftprof
 character(len=500) :: msg,frm
 character(len=fnlen) :: fname
!arrays
 integer :: ngfft_ecut(18,necut)
 real(dp),parameter :: k_gamma(3)=zero
 real(dp) :: ecut_list(necut)
 real(dp),allocatable :: prof_res(:,:,:) 

! ********************************************************************* 

 nsetups = SIZE(fft_setups,DIM=2)
 ecut_list = arth(ecut_arth(1),ecut_arth(2),necut)
 !
 ! Open file and write header with info.
 funt = get_unit()
 write(fname,'(2(a,i1))')"PROF_fourdp_cplex",cplex,"_isign",isign
 open(file=fname,unit=funt)

 write(msg,'(2(a,i0))')"Benchmark: routine = fourdp, cplex =",cplex,", isign=",isign
 write(std_out,*)" Running "//TRIM(msg)

 write(funt,'(a)')"# "//TRIM(msg)
 do set=1,nsetups
   write(funt,'(a,5(a,i0))') "#",&
&    "  fftalg = "   ,fft_setups(1,set),&
&    ", fftcache = " ,fft_setups(2,set),&
&    ", ndat = "     ,fft_setups(3,set),&
&    ", nthreads = " ,fft_setups(4,set),&
&    ", available = ",fft_setups(5,set)
 end do

 ABI_ALLOCATE(prof_res,(2,necut,nsetups))

 do set=1,nsetups
   !
   do iec=1,necut
     call nullify_fft_test(Ftest)
     call init_FFT_test(Ftest,fft_setups(:,set),k_gamma,ecut_list(iec),boxcutmin,rprimd,nsym,symrel,MPI_enreg_in)

     call time_fourdp(Ftest,isign,cplex,Ftprof)

     prof_res(1,iec,set) = Ftprof%cpu_time /Ftprof%ncalls 
     prof_res(2,iec,set) = Ftprof%wall_time/Ftprof%ncalls
     !
     ! Save FFT divisions.
     if (set==1) ngfft_ecut(:,iec) = Ftest%ngfft

     call destroy_FFT_prof(Ftprof)
     call destroy_FFT_test(Ftest)
   end do
   !
 end do
 !
 ! Write the wall-time as a function of ecut.
 write(frm,*)"(f7.1,3i4,",nsetups,"(f7.4))"
 do iec=1,necut
   write(funt,frm) ecut_list(iec),ngfft_ecut(4:6,iec),prof_res(2,iec,:)
 end do

 close(funt)
 ABI_DEALLOCATE(prof_res)

end subroutine prof_fourdp
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/prof_fourwf
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine prof_fourwf(fft_setups,cplex,option,kpoint,necut,ecut_arth,boxcutmin,rprimd,nsym,symrel,MPI_enreg_in)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prof_fourwf'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: nsym,cplex,necut,option
 real(dp),intent(in) :: boxcutmin
 type(MPI_type),intent(inout) :: MPI_enreg_in
!arrays
 integer,intent(in) :: fft_setups(:,:),symrel(3,3,nsym)
 real(dp),intent(in) :: ecut_arth(2)
 real(dp),intent(in) :: kpoint(3),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: iec,nsetups,set,funt,istwf_k   
 type(FFT_test_t) :: Ftest
 type(FFT_prof_t) :: Ftprof
 character(len=500) :: msg,frm
 character(len=fnlen) :: fname
!arrays
 integer :: ngfft_ecut(18,necut)
 real(dp) :: ecut_list(necut)
 real(dp),allocatable :: prof_res(:,:,:) 

! ********************************************************************* 

 nsetups = SIZE(fft_setups,DIM=2)
 ecut_list = arth(ecut_arth(1),ecut_arth(2),necut)
 istwf_k = set_istwfk(kpoint)
 !
 ! Open file and write header with info.
 funt = get_unit()
 write(fname,'(3(a,i1))')"PROF_fourwf_cplex",cplex,"_option",option,"_istwfk",istwf_k
 open(file=fname,unit=funt)

 write(msg,'(3(a,i1))')"Benchmark: routine = fourwf, cplex = ",cplex,", option= ",option,", istwfk= ",istwf_k
 write(std_out,*)" Running "//TRIM(msg)

 write(funt,'(a)')"# "//TRIM(msg)
 do set=1,nsetups
   write(funt,'(a,5(a,i0))') "#",&
&    "  fftalg = "   ,fft_setups(1,set),&
&    ", fftcache = " ,fft_setups(2,set),&
&    ", ndat = "     ,fft_setups(3,set),&
&    ", nthreads = " ,fft_setups(4,set),&
&    ", available = ",fft_setups(5,set)
 end do

 ABI_ALLOCATE(prof_res,(2,necut,nsetups))

 do set=1,nsetups
   !
   do iec=1,necut
     call nullify_fft_test(Ftest)
     call init_FFT_test(Ftest,fft_setups(:,set),kpoint,ecut_list(iec),boxcutmin,rprimd,nsym,symrel,MPI_enreg_in)

     call time_fourwf(Ftest,cplex,option,Ftprof)

     prof_res(1,iec,set) = Ftprof%cpu_time /Ftprof%ncalls 
     prof_res(2,iec,set) = Ftprof%wall_time/Ftprof%ncalls
     !
     ! Save FFT divisions.
     if (set==1) ngfft_ecut(:,iec) = Ftest%ngfft

     call destroy_FFT_prof(Ftprof)
     call destroy_FFT_test(Ftest)
   end do
   !
 end do
 !
 ! Write the wall-time as a function of ecut.
 write(frm,*)"(f7.1,3i4,",nsetups,"(f7.4))"
 do iec=1,necut
   write(funt,frm) ecut_list(iec),ngfft_ecut(4:6,iec),prof_res(2,iec,:)
 end do

 close(funt)
 ABI_DEALLOCATE(prof_res)

end subroutine prof_fourwf
!!***

!----------------------------------------------------------------------

!!****f* m_FFT_prof/prof_rhotwg
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine prof_rhotwg(fft_setups,map2sphere,use_padfft,necut,ecut_arth,osc_ecut,boxcutmin,&
&  rprimd,nsym,symrel,gmet,MPI_enreg_in)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prof_rhotwg'
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: nsym,necut,map2sphere,use_padfft
 real(dp),intent(in) :: boxcutmin,osc_ecut
 type(MPI_type),intent(inout) :: MPI_enreg_in
!arrays
 integer,intent(in) :: fft_setups(:,:),symrel(3,3,nsym)
 real(dp),intent(in) :: ecut_arth(2)
 real(dp),intent(in) :: rprimd(3,3),gmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: iec,nsetups,set,funt,osc_npw
 type(FFT_test_t) :: Ftest
 type(FFT_prof_t) :: Ftprof
 character(len=500) :: msg,frm
 character(len=fnlen) :: fname
!arrays
 integer,pointer :: osc_gvec(:,:)
 integer :: ngfft_ecut(18,necut)
 real(dp),parameter :: k_gamma(3)=zero
 real(dp) :: ecut_list(necut)
 real(dp),allocatable :: prof_res(:,:,:) 

! ********************************************************************* 

 nsetups = SIZE(fft_setups,DIM=2)
 ecut_list = arth(ecut_arth(1),ecut_arth(2),necut)
 !
 ! Open file and write header with info.
 funt = get_unit()
 write(fname,'(2(a,i1))')"PROF_rhotwg_map2sphere",map2sphere,"_use_padfft",use_padfft
 open(file=fname,unit=funt)

 write(msg,'(2(a,i0),a,f5.1)')&
&  "Benchmark: routine = rho_tw_g, map2sphere = ",map2sphere,", use_padfft = ",use_padfft,", osc_ecut = ",osc_ecut
 write(std_out,*)" Running "//TRIM(msg)

 write(funt,'(a)')"# "//TRIM(msg)
 do set=1,nsetups
   write(funt,'(a,5(a,i0))') "#",&
&    "  fftalg = "   ,fft_setups(1,set),&
&    ", fftcache = " ,fft_setups(2,set),&
&    ", ndat = "     ,fft_setups(3,set),&
&    ", nthreads = " ,fft_setups(4,set),&
&    ", available = ",fft_setups(5,set)
 end do

 nullify(osc_gvec)
 call get_kg((/zero,zero,zero/),1,osc_ecut,gmet,osc_npw,osc_gvec)
!  TODO should reorder by shells to be consistent with the GW part!
!  Moreover I guess this ordering is more efficient when we have
!  to map the box to the G-sphere!

 ABI_ALLOCATE(prof_res,(2,necut,nsetups))

 do set=1,nsetups
   !
   do iec=1,necut
     call nullify_fft_test(Ftest)
     call init_FFT_test(Ftest,fft_setups(:,set),k_gamma,ecut_list(iec),boxcutmin,rprimd,nsym,symrel,MPI_enreg_in)

     call time_rhotwg(Ftest,map2sphere,use_padfft,osc_npw,osc_gvec,Ftprof)

     prof_res(1,iec,set) = Ftprof%cpu_time /Ftprof%ncalls 
     prof_res(2,iec,set) = Ftprof%wall_time/Ftprof%ncalls
     !
     ! Save FFT divisions.
     if (set==1) ngfft_ecut(:,iec) = Ftest%ngfft

     call destroy_FFT_prof(Ftprof)
     call destroy_FFT_test(Ftest)
   end do
   !
 end do
 !
 ! Write the wall-time as a function of ecut.
 write(frm,*)"(f7.1,3i4,",nsetups,"(f7.4))"
 do iec=1,necut
   write(funt,frm) ecut_list(iec),ngfft_ecut(4:6,iec),prof_res(2,iec,:)
 end do

 close(funt)
 ABI_DEALLOCATE(prof_res)
 ABI_DEALLOCATE(osc_gvec)

end subroutine prof_rhotwg
!!***

!----------------------------------------------------------------------

END MODULE m_FFT_prof
!!***
