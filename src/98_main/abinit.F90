!{\src2tex{textfont=tt}}
!!****p* ABINIT/abinit
!! NAME
!! abinit
!!
!! FUNCTION
!! Main routine for conducting Density-Functional Theory calculations or Many-Body Perturbation Theory calculations.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MKV, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! The new user is strongly adviced to read the
!! latest version of the file ~abinit/doc/users/new_user_guide.html
!! before trying to modify or even use the code.
!! Even experienced users of the code should also be careful in coding,
!! please read the latest version of the file ~abinit/doc/developers/rules_coding
!!
!! The present main routine drives the following operations :
!!
!! 1) Eventually initialize MPI
!! 2) Initialize overall timing of run
!! 3) Print greeting for interactive user and
!!    Read names of files (input, output, rootinput, rootoutput, roottemporaries),
!!    create the name of the status file, initialize the status subroutine.
!! 4) Open output file and print herald at top of output and log files
!! 5) Read the input file, and store the information in a long string of characters
!! 6) Take ndtset from the input string, then allocate
!!    the arrays whose dimensions depends only on ndtset
!! 7) Continue to analyze the input string, and allocate the remaining arrays.
!!    Also modulate the timing according to timopt.
!! 8) Finish to read the "file" file completely,
!!    and also initialize pspheads (the pseudopotential header information)
!! 9) Provide defaults for the variables that have not yet been initialized.
!! 10) Perform some global initialization, depending on the value of
!! pseudopotentials, parallelism variables, or macro input variables
!! 11) Call the main input routine.
!! 12) Echo input data to output file and log file
!! 13) Perform additional checks on input data
!!  At this stage, all the information from the "files" file and "input" file
!!  have been read and checked.
!! ___________________________________________
!! 14) Perform main calculation  (call driver)
!! -------------------------------------------
!!
!! 15) Give final echo of coordinates, etc.
!! 16) Timing analysis
!! 17) Bibliographical recommendations
!! 18) Delete the status file, and, for build-in tests,
!!       analyse the correctness of results
!! 19) Write the final timing, close the output file, and write a final line
!!       to the log file
!! 20) Eventual cleaning of MPI run
!!
!! TODO
!!  Reduce the number of tasks in this main routine ! Create subroutine to make the full init of input variables,
!!  including the echo ...
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!      ab6_invars_free,ab6_invars_get_abinit_vars,ab6_invars_load
!!      ab6_invars_set_flags,ab6_invars_set_mpi,abi_io_redirect,chkinp,chkvars
!!      date_and_time,destroy_mpi_enreg,destroy_results_out,driver,dump_config
!!      dump_cpp_options,dump_optim,gather_results_out,herald,init_mpi_enreg
!!      init_results_out,iofn1,leave_new,leave_test,mpi_allreduce,out_acknowl
!!      outvars,outxml_finalise,outxml_open,papi_init,parsefile,print_kinds
!!      setdevice_cuda,status,testfi,timab,timana,timein,unsetdevice_cuda
!!      wrtout,xmpi_end
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program abinit

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_parameters
 use m_ab6_invars
 use m_build_info
 use m_cppopts_dumper
 use m_optim_dumper
 use m_profiling
 use m_results_out
 use m_xmpi
 use m_errors
#if defined HAVE_MPI2
 use mpi
#endif
#ifdef HAVE_GPU_CUDA
 use m_initcuda, only: setdevice_cuda,unsetdevice_cuda
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abinit'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_57_iovars
 use interfaces_95_drive
#if defined HAVE_GPU_CUDA
#endif
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments -----------------------------------
!Local variables-------------------------------
!
!===============================================================================
!  abinit_version designate overall code version
!  mpw=maximum number of planewaves in basis sphere
!  unkg,unkgq,unkg1,unkg,unwff1,unwff2,unwffgs,unwfkq,
!  unwft1,unwft2,unylm,unylm1,unpaw,unpaw1 and unpawq...
!!  define input and output unit numbers.
!   These unit numbers are transferred
!   down to the adequate routines.
!   Other unit numbers (ab_in,ab_out,std_out,tmp_unit)
!   have been defined in defs_basis.f .
!  The array filnam is used for the name of input and output files,
!  and roots for generic input, output or temporary files.
!  Pseudopotential file names are set in iofn2, and are contained in pspheads.
!  The name filstat will be needed beyond gstate to check
!  the appearance of the "exit" flag, to make a hasty exit, as well as
!  in order to output the status of the computation.
!==============================================================================
! Declarations
! Define input and output unit numbers (do not forget, unit 5 and 6
! are standard input and output)
! Also, unit number 21, 22 and 23 are used in nstdy3, for the 3 dot
! wavefunctions. Others unit numbers will be used in the case
! of the variational and 2n+1 expressions.
! In defs_basis, one defines :
!  std_in=5, ab_in=5, std_out=6, ab_out=7, tmp_unit=9, tmp_unit2=10
! integer,parameter :: unddb=16,unddk=15,unkg1=19,unkg=17,unkgq=18
! integer,parameter :: unpaw=26,unpaw1=27,unpawq=28,unpos=30
! integer,parameter :: unwff1=1,unwff2=2,unwffgs=3,unwfkq=4,unwft1=11
! integer,parameter :: unwft2=12,unwftgs=13,unwftkq=14,unylm=24,unylm1=25
! Define "level of the routine", for debugging purposes
 integer,parameter :: level=1
 integer :: choice,dmatpuflag,ierr,iexit,ii,iounit,istatr,istatshft,dtsetsId
 integer :: lenstr
 integer :: mu,natom,ndtset
 integer :: ndtset_alloc,nfft,nkpt,npsp
 integer :: nsppol,papiopt,prtvol,timopt,use_gpu_cuda
 integer,allocatable :: nband(:),npwtot(:)
 real(dp) :: cpui,etotal,walli
 real(dp) :: strten(6),tsec(2)
 real(dp),allocatable :: fred(:,:),xred(:,:)
 character(len=24) :: codename
 character(len=5000) :: message
 character(len=strlen) :: string
 character(len=fnlen) :: filstat
 character(len=fnlen) :: filnam(5)
 type(dataset_type),pointer  :: dtsets(:)
 type(MPI_type) :: mpi_enreg
 type(pspheader_type),pointer :: pspheads(:)
 type(results_out_type),allocatable,target :: results_out(:)
 type(results_out_type),pointer :: results_out_all(:)
 type(ab_maxvals) :: mxvals
 logical :: test_img,use_results_all,xml_output=.false.
 integer :: values(8)
 character(len=5) :: strzone
 character(len=8) :: strdat
 character(len=10) :: strtime
#if defined HAVE_MPI
 real(dp) :: tsec_s(2)
#endif

!******************************************************************
!
!0) Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

!1) Eventually initialize MPI
!Pay attention: it may be initialzed again in finddistrproc
 call init_mpi_enreg(mpi_enreg)

!Init global entities defined in the abinit modules.
!call abi_init_globals()

!2) Initialize overall timing of run:
#ifdef HAVE_TIMER_PAPI
 call papi_init()
#endif

 call timein(cpui,walli)

 call timab(1,0,tsec)

!Start to accumulate time for the entire run. The end of accumulation is in timana.f
 call timab(1,1,tsec)

!3) Print greeting for interactive user,
!read names of files (input, output, rootinput, rootoutput, roottemporaries),
!create the name of the status file, initialize the status subroutine.

 call timab(41,3,tsec)

 call iofn1(filnam,filstat,mpi_enreg)

 iexit=99
 call status(0,filstat,iexit,level,'first status  ')

!4) Open output file and print herald at top of output and log files

 if(mpi_enreg%me==0)then

   open (unit=ab_out,file=filnam(2),form='formatted',status='new')
   rewind (unit=ab_out)
   codename='ABINIT'//repeat(' ',18)
   call herald(codename,abinit_version,ab_out)
   call herald(codename,abinit_version,std_out)
   call dump_config()
   call dump_optim()
   call dump_cpp_options()
!  Write names of files
   write(message, '(a,a,a,a,a,a,a,a,a,a,a,a)' )&
&   '- input  file    -> ',trim(filnam(1)),ch10,&
&   '- output file    -> ',trim(filnam(2)),ch10,&
&   '- root for input  files -> ',trim(filnam(3)),ch10,&
&   '- root for output files -> ',trim(filnam(4)),ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!5) Read the file, stringify it and return the number of datasets.
 call status(0,filstat,iexit,level,'call parsefile')
 call parsefile(filnam(1), lenstr, ndtset, string)

!6~11) Call the parser from the parser module.
 call ab6_invars_set_flags(.true., .true., status_file = filstat, &
& timab_tsec = tsec)

 call ab6_invars_set_mpi(mpi_enreg)

 call ab6_invars_load(dtsetsId, string, lenstr, ndtset, .true., .true.)

 call timab(44,1,tsec)

 if(mpi_enreg%me>=0) then
   call ab6_invars_get_abinit_vars(dtsetsId, dtsets, pspheads,&
&   mxvals, istatr, istatshft, papiopt, timopt, dmatpuflag)
   ndtset_alloc = size(dtsets) - 1
   npsp = size(pspheads)

!  12) Echo input data to output file and log file

   call status(0,filstat,iexit,level,'call outvars(1')

!  For evolving variables, and results
   ABI_ALLOCATE(results_out,(0:ndtset_alloc))

!  Initialize results_out datastructure
   call init_results_out(dtsets,1,1,&
&   mxvals%mxnatom,&
&   mxvals%mxmband_upper,&
&   mxvals%mxnkpt,&
&   mxvals%mxnsppol,&
&   results_out)

!  Gather contributions to results_out from images of the cell, if needed
   test_img=(mxvals%mxnimage/=1.and.maxval(dtsets(:)%npimage)>1)
   use_results_all=.false.
   if (test_img) then
     use_results_all=(mpi_enreg%me==0)
     if (use_results_all)  then
       ABI_ALLOCATE(results_out_all,(0:ndtset_alloc))
     end if
     call gather_results_out(dtsets,results_out,results_out_all,use_results_all,&
&     allgather=.false.,master=0)
   else
     results_out_all => results_out
   end if

   if(mpi_enreg%me==0) then

!    Echo input to output file on unit ab_out, and to log file on unit 06 :
     choice=1
     do ii=1,2
       if(ii==1)iounit=ab_out
       if(ii==2)iounit=std_out

       call outvars (choice,dmatpuflag,dtsets,trim(filnam(4)),&
&       iounit,istatr,istatshft,mpi_enreg, mxvals,&
&       ndtset,ndtset_alloc,npsp,results_out_all,timopt)
     end do

     if (dtsets(1)%prtxml == 1) then
       call outxml_open(trim(filnam(4)))
       call date_and_time(strdat,strtime,strzone,values)
       xml_output = .true.
     else
       xml_output = .false.
     end if

   end if ! End of me==0 section

!  Clean memory
   if (test_img.and.mpi_enreg%me==0) then
     call destroy_results_out(results_out_all)
     ABI_DEALLOCATE(results_out_all)
   end if

!  This synchronization is not strictly needed, but without it,
!  there are problems with Tv1#93 in parallel, PGI compiler, on Intel/PC
   call leave_test()

   call timab(44,2,tsec)

!  13) Perform additional checks on input data

   call timab(45,3,tsec)
   call status(0,filstat,iexit,level,'call chkinp   ')


   call chkinp(dtsets,ab_out,mpi_enreg,ndtset,ndtset_alloc,npsp,pspheads)

!  Check whether the string only contains valid keywords
   call chkvars(string)

!  At this stage, all the information from the "files" file and "input" file
!  have been read and checked.

!  14) Perform main calculation
!  The timing is done in gstate

   call status(0,filstat,iexit,level,'call driver   ')

   prtvol=dtsets(1)%prtvol
   if(prtvol==-level)then
     write(message,'(a1,a,a1,a,i1,a)') ch10,' abinit : before driver ',&
&     ch10,'  prtvol=-',level,', debugging mode => stop '
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  Printout of kinds and precisions.
   if (mpi_enreg%me==0) call print_kinds(std_out)

!  Eventually activate GPU
   use_gpu_cuda=0
#if defined HAVE_GPU_CUDA
   do ii=1,ndtset_alloc
     if (dtsets(ii)%use_gpu_cuda==1) use_gpu_cuda=1
   end do
   call setdevice_cuda(use_gpu_cuda)
#endif

   call timab(45,2,tsec)

   call driver(abinit_version,cpui,dtsets,filnam,filstat,&
&   mpi_enreg,ndtset,ndtset_alloc,npsp,pspheads,results_out)

   call status(0,filstat,iexit,level,'after driver  ')

!  15) Give final echo of coordinates, etc.
   call timab(46,1,tsec)

   write(message,'(a,a,a,62a,80a)') ch10,&
&   '== END DATASET(S) ',('=',mu=1,62),ch10,('=',mu=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

!  Gather contributions to results_out from images of the cell, if needed
   if (test_img) then
     if (use_results_all)  then
       ABI_ALLOCATE(results_out_all,(0:ndtset_alloc))
     end if
     call gather_results_out(dtsets,results_out,results_out_all,use_results_all,&
&     allgather=.false.,master=0)
   end if

   if(mpi_enreg%me==0) then
!    Echo input to output file on unit ab_out, and to log file on unit std_out.
     choice=2
     do ii=1,2
       if(ii==1)iounit=ab_out
       if(ii==2)iounit=std_out
       write(iounit,*)' '
       call outvars (choice,dmatpuflag,dtsets,trim(filnam(4)),&
&       iounit,istatr,istatshft,mpi_enreg,mxvals,&
&       ndtset,ndtset_alloc,npsp,results_out_all,timopt)
       if(ii==2)write(std_out,*)' '
     end do
   end if ! mpi_enreg%me==0

!  Clean memory
   if (test_img.and.mpi_enreg%me==0) then
     call destroy_results_out(results_out_all)
     ABI_DEALLOCATE(results_out_all)
   else
     nullify(results_out_all)
   end if

!  In prevision of the next two calls, some variables need to be transfered.
!  They concern the case ndtset<2, and nimage=1 so take first value.
   natom=dtsets(1)%natom ; nkpt=dtsets(1)%nkpt ; nsppol=dtsets(1)%nsppol
   nfft=dtsets(1)%nfft   ; etotal=results_out(1)%etotal(1)
   ABI_ALLOCATE(nband,(nkpt*nsppol))
   ABI_ALLOCATE(npwtot,(nkpt))
   ABI_ALLOCATE(fred,(3,natom))
   ABI_ALLOCATE(xred,(3,natom))
   fred(:,:)  =results_out(1)%fred(:,1:natom,1)
   nband(:)   =dtsets(1)%nband(1:nkpt*nsppol)
   npwtot(:)  =results_out(1)%npwtot(1:nkpt,1)
   strten(:)  =results_out(1)%strten(:,1)
   xred(:,:)  =results_out(1)%xred(:,1:natom,1)

   call timab(46,2,tsec)

!  16) Timing analysis

   if(timopt/=0)then
     call status(0,filstat,iexit,level,'call timana   ')
     call timana (mpi_enreg, natom, nband, ndtset, nfft, nkpt, npwtot, nsppol, timopt, papiopt)
   else
#if defined HAVE_MPI
     if(mpi_enreg%me==0)then ! This is for the automatic tests
       write(ab_out,'(5a)')ch10,ch10,'- Timing analysis has been suppressed with timopt=0',ch10,ch10
     end if
#endif
   end if

!  17) Bibliographical recommendations
   if(mpi_enreg%me==0) then
     do ii=1,2
       if(ii==1)iounit=ab_out
       if(ii==2)iounit=std_out
       call out_acknowl(dtsets,iounit,ndtset_alloc,npsp,pspheads)
     end do
   end if ! mpi_enreg%me==0


!  18) Delete the status file, and, for build-in tests,
!  analyse the correctness of results.

   if(ndtset==0)then
     call testfi(etotal,filnam,filstat,fred,natom,strten,xred)
   else
     open (tmp_unit,file=filstat,form='formatted',status='unknown')
     close (tmp_unit,status='delete')
!    $call delete_file(filstat)
   end if

!  One should have here the explicit deallocation of all arrays
   call destroy_results_out(results_out)

   call ab6_invars_free(dtsetsId)

!  deallocate(dtsets)
!  deallocate(pspheads)
   ABI_DEALLOCATE(fred)
   ABI_DEALLOCATE(nband)
   ABI_DEALLOCATE(npwtot)
   ABI_DEALLOCATE(results_out)
   ABI_DEALLOCATE(xred)

!  19) Write the final timing, close the output file, and write a final line to the log file

   call timein(tsec(1),tsec(2))
   tsec(1)=tsec(1)-cpui
   tsec(2)=tsec(2)-walli

#if defined HAVE_MPI
   write(std_out, '(a,i4,a,f13.1,a,f13.1)' )&
&   ' Proc.',mpi_enreg%me,' individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
   if(mpi_enreg%me==0)then
     write(ab_out, '(a,a,a,i4,a,f13.1,a,f13.1)' ) &
&     '-',ch10,'- Proc.',mpi_enreg%me,' individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
   end if
   call MPI_ALLREDUCE(tsec,tsec_s,2,MPI_DOUBLE_PRECISION,MPI_SUM,xmpi_world,ierr)
   tsec=tsec_s
#endif

   write(message, '(a,80a,a,a,a)' ) ch10,('=',mu=1,80),ch10,ch10,' Calculation completed.'
   call wrtout(ab_out,message,'COLL')

   if(mpi_enreg%me==0)then
     write(ab_out, '(a,f13.1,a,f13.1)' ) &
&     '+Overall time at end (sec) : cpu=',tsec(1),'  wall=',tsec(2)
     close (unit=ab_out)

     if (xml_output) then
       call outxml_finalise(tsec, values)
     end if

     write(message, '(a,a)' ) ch10,' Calculation completed.'
     call wrtout(std_out,message,'COLL')
   end if
 end if
!20) Eventual cleaning of MPI (and/or GPU) run
 call destroy_mpi_enreg(MPI_enreg)
 call xmpi_end()
#if defined HAVE_GPU_CUDA
 call unsetdevice_cuda(use_gpu_cuda)
#endif

 end program abinit
!!***
