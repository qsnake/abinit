!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_timer
!! NAME
!! m_timer
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_timer

 use m_profiling

 use defs_basis
 use m_build_info
 use m_optim_dumper
 use m_cppopts_dumper
 !%use m_errors TODO move m_timer to a lower level
 !%use m_io_tools,   only : get_unit

 implicit none

 private

 integer,parameter,public :: TIMER_KEY_LENGTH = 100

!----------------------------------------------------------------------

! private variables.
 integer,parameter :: TIMER_SIZE_DEFAULT = 2000
 integer,save      :: TIMER_SIZE         = TIMER_SIZE_DEFAULT
 integer,save      :: NCOLLISIONS        = 0
 real(dp),save     :: TIMER_CPU_TIME0 
 real(dp),save     :: TIMER_WALL_TIME0

#ifdef DEBUG_MODE
 logical,parameter :: errors_are_fatal=.TRUE.
#else
 logical,parameter :: errors_are_fatal=.FALSE.
#endif
!!***

!!****t* m_timer/entry_type
!! NAME
!! entry_type
!! 
!! FUNCTION
!!  Structure storing timing data associated to a section of the code.
!! 
!! SOURCE

 type,private :: entry_type
   character(len=TIMER_KEY_LENGTH) :: key=""
   ! The name of the entry.

   integer :: ncalls=0
   ! The number of time the section of code has been executed.

   integer(i1b) :: status=0
   ! 0 --> entry is not defined.
   ! 1 --> timing analysis is in execution.
   ! 2 --> timing analysis is completed.

   real(dp) :: cpu_time_tot=zero
   ! The total CPU time of the entry.

   real(dp) :: wall_time_tot=zero
   ! The total wall time of the entry.

   real(dp) :: tzero(2)
   ! Temporary array used to store the initial CPU and the WALL time.
 end type entry_type

 type(entry_type),save,allocatable :: Timer(:)

 public :: init_timer        ! Initialize the timer
 public :: destroy_timer     ! Destroy the timer.
 public :: timing            ! Time a section of the code.
 public :: print_timer       ! Printout of the timer.

CONTAINS  !===========================================================
!!***

!!****f* m_timer/init_timer
!! NAME
!!  init_timer
!!
!! FUNCTION
!!  Allocate the Timer. Size is given by TIMER_SIZE.
!!
!! PARENTS
!!      driver,m_timer
!!
!! CHILDREN
!!      dump_config,dump_cpp_options,dump_optim,leave_new,timein,wrtout
!!
!! SOURCE

subroutine init_timer()

!Local variables-------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_timer'
 use interfaces_18_timing
!End of the abilint section

 integer :: istat
 
! *************************************************************************

 if (allocated(Timer))  then
   ABI_DEALLOCATE(Timer)
 end if
 ABI_ALLOCATE(Timer,(TIMER_SIZE))
 istat = ABI_ALLOC_STAT
 if (istat/=0) stop "Out-of-memory while allocating Timer"

 NCOLLISIONS = 0
 call timein(TIMER_CPU_TIME0,TIMER_WALL_TIME0)

end subroutine init_timer
!!***

!----------------------------------------------------------------------

!!****f* m_timer/destroy_timer
!! NAME
!!  destroy_timer
!!
!! FUNCTION
!!  Deallocate the Timer, reset NCOLLISIONS and restore TIMER_SIZE to the default value.
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      dump_config,dump_cpp_options,dump_optim,leave_new,timein,wrtout
!!
!! SOURCE

subroutine destroy_timer()

! *************************************************************************

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_timer'
!End of the abilint section

 TIMER_SIZE = TIMER_SIZE_DEFAULT
 ABI_DEALLOCATE(Timer)

end subroutine destroy_timer
!!***

!----------------------------------------------------------------------

!!****f* m_timer/hash
!! NAME
!!  hash
!!
!! FUNCTION
!!  Hashing function used to map the keyword of the timer onto the array index.
!!
!! INPUTS
!!  string=The keyword associated to the timed section.
!!
!! OUTPUT
!!  hash=The hash value corresponding to string.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function hash(string)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hash'
!End of the abilint section

 integer :: hash
 character(len=*),intent(in) :: string

!Local variables-------------------------------
!scalars
 integer :: ii

! *************************************************************************

 hash = 0 
 do ii=1,LEN_TRIM(string)
   hash = hash + ICHAR(string(ii:ii)) ! TODO not efficient at all! First coding.
 end do

 hash = MOD(hash,TIMER_SIZE)+1

end function hash
!!***

!----------------------------------------------------------------------

!!****f* m_timer/timing
!! NAME
!!  timing
!!
!! FUNCTION
!!  Timing routine. 
!!
!! INPUTS
!!  key=The keyword with the name of the section.
!!  topt= 0 to start the timer, 1 to stop
!!
!! SIDE EFFECTS
!!  Timer
!!
!! PARENTS
!!
!! CHILDREN
!!      dump_config,dump_cpp_options,dump_optim,leave_new,timein,wrtout
!!
!! SOURCE

subroutine timing(topt,key)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'timing'
 use interfaces_18_timing
!End of the abilint section

 integer,intent(in) :: topt 
 character(len=*),intent(in) :: key 

!Local variables-------------------------------
!scalars
 integer :: ihs
 real(dp) :: cpu,wall
 !character(len=500) :: msg
 
! *************************************************************************
 
 ihs = set_index_from_key(key)

 Timer(ihs)%key = key

 call timein(cpu,wall)

 SELECT CASE (topt)

 CASE (1) ! start the timer for this entry.
   Timer(ihs)%ncalls = Timer(ihs)%ncalls +1

   Timer(ihs)%tzero(1)=cpu
   Timer(ihs)%tzero(2)=wall

   if (ALL(Timer(ihs)%status /= (/0,2/))) then
    write(std_out,'(a,i0)')" Status of entry "//TRIM(key)//" is ",Timer(ihs)%status
    write(std_out,*)"Likely misplaced or missing call to ABI_TIMER_START or ABI_TIMER_STOP"
    if (errors_are_fatal) stop
   end if

   Timer(ihs)%status = 1

 CASE (2) ! accumulate time for this entry.

   Timer(ihs)%cpu_time_tot  = Timer(ihs)%cpu_time_tot  + cpu  -Timer(ihs)%tzero(1)
   Timer(ihs)%wall_time_tot = Timer(ihs)%wall_time_tot + wall -Timer(ihs)%tzero(2)

   if (ALL(Timer(ihs)%status /= (/0,1/))) then
     write(std_out,'(a,i0)')" Status of entry "//TRIM(key)//" is ",Timer(ihs)%status
     write(std_out,*)"Likely misplaced or missing calls to ABI_TIMER_START or ABI_TIMER_STOP"
     if (errors_are_fatal) stop
   end if

   Timer(ihs)%status = 2

 CASE DEFAULT
   write(std_out,'(a,i0,2a)')" Wrong value for topt: ",topt," key= ",TRIM(key)
   stop
 END SELECT

end subroutine timing
!!***

!----------------------------------------------------------------------

!!****f* m_timer/enlarge_timer
!! NAME
!!  enlarge_timer
!!
!! FUNCTION
!!  Reallocate Timer with larger dimension while preserving previous data.
!!
!! INPUTS
!!  new_size=The new size of the timer
!!
!! SIDE EFFECTS 
!!  Timer reallocated with new size. TIMER_SIZE is changed accordingly.
!!
!! PARENTS
!!      m_timer
!!
!! CHILDREN
!!      dump_config,dump_cpp_options,dump_optim,leave_new,timein,wrtout
!!
!! SOURCE

subroutine enlarge_timer(new_size)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'enlarge_timer'
!End of the abilint section

 integer,intent(in) :: new_size

!Local variables-------------------------------
!scalars
 integer :: old_size,istat
!arrays
 type(entry_type),allocatable :: Tmp_Timer(:)
 
! *************************************************************************
 
 old_size = SIZE(Timer)

 if (new_size<=old_size) stop "New_size too short"

 ABI_ALLOCATE(Tmp_Timer,(new_size))
 istat = ABI_ALLOC_STAT
 if (istat/=0) stop "Out-of-memory in Tmp_Timer"

 Tmp_Timer(1:old_size) = Timer

 TIMER_SIZE = new_size
 call init_timer()

 Timer(1:old_size) = Tmp_Timer
 ABI_DEALLOCATE(Tmp_Timer)

end subroutine enlarge_timer
!!***

!----------------------------------------------------------------------

!!****f* m_timer/set_index_from_key
!! NAME
!!  set_index_from_key
!!
!! FUNCTION
!!  Set the index associated to the keyword. If a collision occurs, key is 
!!  given the first subsequent index that is not being used by any another entry.
!!  If all entries are occupied, Timer is re-allocated with larger size and the first
!!  new index is used.
!!
!! INPUTS
!!  key=The keyword associated to the timed entry.
!!
!! SIDE EFFECTS
!!  Timer. See description above.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function set_index_from_key(key)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'set_index_from_key'
!End of the abilint section

 integer :: set_index_from_key
 character(len=*),intent(in) :: key 

!Local variables-------------------------------
!scalars
 integer :: ihs
 logical :: collide
 !character(len=500) :: msg
 
! *************************************************************************

 ihs = hash(key)
                                                             
 collide = (Timer(ihs)%key /= key .and. Timer(ihs)%key/="")
 !collide = (Timer(ihs)%key /= key .and. Timer(ihs)%status/=0)
                                                             
 do while (collide .and. ihs<TIMER_SIZE)
   write(std_out,*)TRIM(key)," collides with ",TRIM(Timer(ihs)%key)
   NCOLLISIONS = NCOLLISIONS +1
   ihs = ihs+1
   collide = (Timer(ihs)%key /= key .and. Timer(ihs)%key/="")
   !collide = (Timer(ihs)%key /= key .and. Timer(ihs)%status/=0)
 end do
 !
 ! Enlarge the buffer if all slots are occupied.
 if (ihs==TIMER_SIZE.and.collide) then
   call enlarge_timer(TIMER_SIZE + TIMER_SIZE_DEFAULT)
   ihs = ihs+1
 end if

 set_index_from_key = ihs

end function set_index_from_key
!!***

!----------------------------------------------------------------------

!!****f* m_timer/get_index_from_key
!! NAME
!!  get_index_from_key
!!
!! FUNCTION
!!  Returns the index associated to the keyword without changing the time. Returns 0 if not found.
!!
!! INPUTS
!!  key=The keyword associated to the timed entry.
!!
!! OUTPUT
!!  The index of key. 0 if key is not found.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function get_index_from_key(key)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_index_from_key'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: get_index_from_key
 character(len=*),intent(in) :: key 

!Local variables-------------------------------
!scalars
 integer :: ihs
 logical :: found
! *************************************************************************

 ihs = hash(key)

 found = .FALSE.
 do while (.not.found .and. ihs<TIMER_SIZE)
   found = (Timer(ihs)%key == key)
   ihs = ihs+1
 end do

 if (found) then  
   get_index_from_key = ihs-1
 else
   get_index_from_key = 0
 end if

end function get_index_from_key
!!***

!----------------------------------------------------------------------

!!****f* m_timer/print_timer
!! NAME
!!  print_timer
!!
!! FUNCTION
!!
!! INPUTS
!!  fname=File name for printing. Use STDOUT to print on the standard output.
!!  [keysel]=List of keys to be printed. By default all entries are printed. 
!!  [prtvol]=Integer flag governing verbosity output (defaults to 0)
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      dump_config,dump_cpp_options,dump_optim,leave_new,timein,wrtout
!!
!! SOURCE

subroutine print_timer(fname,keysel,prtvol)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_timer'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
!End of the abilint section

 integer,optional,intent(in) :: prtvol
 character(len=fnlen),intent(in) :: fname
 character(len=TIMER_KEY_LENGTH),optional,intent(in) :: keysel(:)

!Local variables-------------------------------
!scalars
 integer :: ii,my_prtvol,field1_w,ikey,unt,ios
 real(dp) :: cpu_now,wall_now
 character(len=500) :: ofmt 
 
! *************************************************************************

 my_prtvol=0; if (PRESENT(prtvol)) my_prtvol=prtvol

 if (fname=="STDOUT") then
   unt = std_out
 else
   !unt = get_unit()
   unt = tmp_unit
   open(file=fname,unit=unt,form="formatted",status="replace",iostat=ios)
   if (ios/=0) then
     call wrtout(std_out,'Error Opening file '//TRIM(fname),"COLL")
     call leave_new("COLL")
   end if
 end if
 !

 if (my_prtvol>0) then

   ! Dump build info.
   write(unt,'(a)')"<config>"
   call dump_config(unt)
   write(unt,'(a)')"</config>"
  
   write(unt,'(a)')"<optim>"
   call dump_optim(unt)
   write(unt,'(a)')"</optim>"
  
   write(unt,'(a)')"<cpp_options>"
   call dump_cpp_options(unt)
   write(unt,'(a)')"</cpp_options>"
   !
   ! Dump timer info.
   call timein(cpu_now,wall_now)
   write(unt,'(a)')"<timer_info>"
   write(unt,*)" timer_cpu_time = ",cpu_now - TIMER_CPU_TIME0
   write(unt,*)" timer_wall_time = ",wall_now - TIMER_WALL_TIME0
   write(unt,*)" timer_size = ",TIMER_SIZE
   write(unt,*)" ncollisions = ",NCOLLISIONS
   write(unt,'(a)')"</timer_info>"

 endif
 !
 ! Dump timing results.
 !
 field1_w=0 ! Width of the field used to print key names.
 do ii=1,SIZE(Timer)
   if (Timer(ii)%status/=0) field1_w = MAX(field1_w, LEN_TRIM(Timer(ii)%key))
 end do

 if (field1_w==0) RETURN

 write(unt,'(a)')"<timer_data>"
 write(unt,*)" Entry name     CPU-time    WALL-time   ncalls"

 write(ofmt,*)"(a",field1_w,",2x,2(f12.3,2x),3x,i0)"
 if (.not.PRESENT(keysel)) then ! write entire set of keys.
   do ii=1,SIZE(Timer)
     if (Timer(ii)%status==0) CYCLE
     write(unt,ofmt)TRIM(Timer(ii)%key),Timer(ii)%cpu_time_tot,Timer(ii)%wall_time_tot,Timer(ii)%ncalls
   end do

 else  ! write only keysel keys.
   do ikey=1,SIZE(keysel)
     ii = get_index_from_key(keysel(ikey))
     if (ii>0 .and. Timer(ii)%status/=0) then
       write(unt,ofmt)TRIM(Timer(ii)%key),Timer(ii)%cpu_time_tot,Timer(ii)%wall_time_tot,Timer(ii)%ncalls
     end if
   end do
 end if

 write(unt,'(a)')"</timer_data>"

 if (unt/=std_out) close(unt)

end subroutine print_timer
!!***

!----------------------------------------------------------------------

END MODULE m_timer
!!***
