!{\src2tex{textfont=tt}}
!!****f* ABINIT/herald
!! NAME
!!  herald
!!
!! FUNCTION
!!  Prints out a message to unit iout giving info about current
!!  code, version of code, platform, and starting date.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, LSI, MM, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  code_name= code name
!!  code_version= code version
!!  iout=unit number for output
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      abinit,aim,anaddb,cut3d,fftprof,mrgddb,mrggkk,mrgscr,newsp,optic,ujdet
!!
!! CHILDREN
!!      date_and_time,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine herald(code_name,code_version,iout)

 use m_profiling

 use defs_basis
 use m_build_info
 use m_build_info_fake

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'herald'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: iout
 character(len=24),intent(in) :: code_name
 character(len=6),intent(in) :: code_version

!Local variables-------------------------------
 integer :: day,dd,ja,jy,jm,jdn,mm,mm_rel,yyyy,yyyy_rel
 integer :: values(8)
 character(len=5) :: strzone
 character(len=8) :: strdat
 character(len=10) :: strtime
 character(len=500) :: message
!no_abirules
 character(len=3) :: daynam(7)=(/'Mon','Tue','Wed','Thu','Fri','Sat','Sun'/)
 character(len=3), parameter :: monnam(12)=(/'Jan','Feb','Mar','Apr','May','Jun',&
&                                       'Jul','Aug','Sep','Oct','Nov','Dec'/)

! *************************************************************************

!RELEASE TIME FROM ABIRULES
 yyyy_rel=2012
 mm_rel=3
!END OF RELEASE TIME

!The technique used hereafter is the only one that we have found to obtain
!perfect transferability across platforms and OS.
 write(iout, '(/,a,a,a,a,a)' ) &
& '.Version ',code_version,' of ',trim(code_name),' '
#if defined HAVE_MPI
 write(iout, '(a,a,a,/)' ) '.(MPI version, prepared for a ',&
& build_target,' computer) '
#else
 write(iout, '(a,a,a,/)' ) '.(sequential version, prepared for a ',&
& build_target,' computer) '
#endif

!GNU GPL license
 write(iout, '(a,/,a,a,a,/,a,/,a,/,a,/)' ) &
& '.Copyright (C) 1998-2012 ABINIT group . ',&
& ' ',trim(code_name),' comes with ABSOLUTELY NO WARRANTY.',&
& ' It is free software, and you are welcome to redistribute it',&
& ' under certain conditions (GNU General Public License,',&
& ' see ~abinit/COPYING or http://www.gnu.org/copyleft/gpl.txt).'

 if(trim(code_name)=='OPTIC')then
   write(iout, '(a,a,a,/,a,/,a,/,a,/,a,/)' ) &
&   ' ',trim(code_name),' has originally been developed by',&
&   ' Sangeeta Sharma and incorporated in ABINIT with the help of M. Verstraete.',&
&   ' Please refer to : ',&
&   ' S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl, Phys. Rev. B 67, 165332 (2003), and',&
&   ' S. Sharma and C. Ambrosch-Draxl, Physica Scripta T 109 (2004).'
 end if

 write(iout, '(a,/,a,/,a,/,a,/,a)' ) &
& ' ABINIT is a project of the Universite Catholique de Louvain,',&
& ' Corning Inc. and other collaborators, see ~abinit/doc/developers/contributors.txt .',&
& ' Please read ~abinit/doc/users/acknowledgments.html for suggested',&
& ' acknowledgments of the ABINIT effort.',&
& ' For more information, see http://www.abinit.org .'

!Get year, month and day
 call date_and_time(strdat,strtime,strzone,values)
 yyyy=values(1)
 mm=values(2)
 dd=values(3)

!Get day of the week
 if (mm.gt.2) then
   jy=yyyy
   jm=mm+1
 else
   jy=yyyy-1
   jm=mm+13
 end if
 jdn=int(365.25d0*jy)+int(30.6001d0*jm)+dd+1720995
 ja=int(0.01d0*jy)
 jdn=jdn+2-ja+int(quarter*ja)
 day=mod(jdn,7)+1

!Print date in nice format (* new format *)
 write(iout, '(/,a,a,1x,i2,1x,a,1x,i4,a,/,a,i2,a,i2,a)' ) &
& '.Starting date : ',daynam(day),dd,monnam(mm),yyyy,'.','- ( at ',values(5),'h',values(6),' )'

!Print date in nice format (* old format *)
!write(iout, '(a,a,1x,i2,1x,a,1x,i4,a)' ) &
!& '.Starting date : ',daynam(day),dd,monnam(mm),yyyy,'.'

 write(iout,*)' '

!Impose a maximal life cycle of 3 years
 if(yyyy>yyyy_rel+3 .or. (yyyy==yyyy_rel+3 .and. mm>mm_rel) ) then 
   write(message, '(8a,i4,5a)' ) ch10,&
&   ' herald :  WARNING -',ch10,&
&   '  The starting date is more than 3 years after the initial release',ch10,&
&   '  of this version of ABINIT, namely ',monnam(mm_rel),' ',yyyy_rel,'.',ch10,&
&   '  This version of ABINIT is not supported anymore.',ch10,&
&   '  Action : please, switch to a more recent version of ABINIT.'
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')

!  Gives a warning beyond 2 years
 else if(yyyy>yyyy_rel+2 .or. (yyyy==yyyy_rel+2 .and. mm>mm_rel) ) then
   write(message, '(8a,i4,6a)' ) ch10,&
&   ' herald :  WARNING -',ch10,&
&   '  The starting date is more than 2 years after the initial release',ch10,&
&   '  of this version of ABINIT, namely ',monnam(mm_rel),' ',yyyy_rel,'.',ch10,&
&   '  Note that the use beyond 3 years after the release will not be supported.',ch10,&
&   '  Action : please, switch to a more recent version of ABINIT.',ch10
   call wrtout(std_out,message,'COLL')
!  Writing to iout might make the tests fail
!  call wrtout(iout,message,'COLL')
 end if

end subroutine herald
!!***
