!{\src2tex{textfont=tt}}
!!****f* ABINIT/rdnpw
!! NAME
!! rdnpw
!!
!! FUNCTION
!! Read the line that contains npw from a kg file (option=0)
!! or wf file (option=1 or option==2).
!! Then, skip the next line.
!! Also performs some checks, related to npw_k and nband_k.
!! The arguments ikpt and isppol are only needed for the error message
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ikpt=index of current k point (needed for error message)
!!  isppol=spin polarization currently treated (needed for error message)
!!  nband_k=number of bands at this k-point (input if option==0 or option==1)
!!  npw_k=number of plane waves at this k-point (input if option==0 or option==1)
!!  option=if 0, read a kg file; if 1 or 2, read a wf file
!!   this is important for the content of the line and error messages
!!   if option=2, no checking
!!  unitfile=unit of the file to be read
!!
!! OUTPUT
!!  if option==2, npw_k, nspinor and nband_k are output
!!
!! SIDE EFFECTS
!!  nspinor=number of spinorial components of the wavefunctions (input if option==0 or option==1)
!!
!! NOTES
!!
!! PARENTS
!!      ctocprj,dyfnl3,eltfrkin3,eltfrnl3,energy,forstrnps,ladielmt,lavnl,mkrho
!!      mkrho3,newkpt,nselt3,nstdy3,nstpaw3,optics_paw,optics_vloc,pawmkaewf
!!      prctfvw1,prctfvw2,rhofermi3,suscep_dyn,suscep_kxc_dyn,suscep_stat,tddft
!!      vtorho,vtorho3
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,option,unitfile)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rdnpw'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikpt,isppol,option,unitfile
 integer,intent(inout) :: nband_k,npw_k,nspinor

!Local variables-------------------------------
!scalars
 integer :: ierr,nband_disk,npw_disk,nspinor_disk
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,*) ' rdnpw : enter, debug, unitfile= ',unitfile
!ENDDEBUG

 nband_disk=0
 if(option==0) read(unitfile,IOSTAT=ierr)npw_disk
 if(option==1 .or. option==2) read(unitfile,IOSTAT=ierr)npw_disk,nspinor_disk,nband_disk

 if(ierr/=0)then
   write(message, '(4a,i5,a,i5)' ) ch10,&
&   ' rdnpw : BUG -',ch10,&
&   '  Reading npw record of disk file unit',unitfile,', gives iostat=',ierr
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

!DEBUG
!write(std_out,*) ' rdnpw : option,npw_disk,nspinor_disk,nband_disk=',option,npw_disk,nspinor_disk,nband_disk
!ENDDEBUG

 if(option==0)then

!  Check agreement with npw_k
   if (npw_k/=npw_disk) then
     write(message, '(4a,i5,a,i2,3a,i6,a,i6,2a,i6,a)' ) ch10,&
&     ' rdnpw : BUG -',ch10,&
&     '  At k point number',ikpt,', with spin polarization',isppol,',',ch10,&
&     '  for kg disk file unit',unitfile,', the value of npw_disk=',npw_disk,ch10,&
&     '  disagrees with argument npw_k=',npw_k,'. IO problem.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if

 else if(option==1)then

!  Check agreement with npw_k
   if (npw_k/=npw_disk) then
     write(message, '(a,a,a,a,i5,a,i2,a,a,a,i6,a,i6,a,a,i6,a)' ) ch10,&
&     ' rdnpw : BUG -',ch10,&
&     '  At k point number',ikpt,', with spin polarization',isppol,',',ch10,&
&     '  for wf file unit',unitfile,', the value of npw_disk=',npw_disk,ch10,&
&     '  disagrees with argument npw_k=',npw_k,'. IO problem.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if

!  Check agreement with nspinor
   if (nspinor/=nspinor_disk) then
     write(message, '(a,a,a,a,i5,a,i2,a,a,a,i6,a,i6,a,a,i6,a)' ) ch10,&
&     ' rdnpw : BUG -',ch10,&
&     '  At k point number',ikpt,', with spin polarization',isppol,',',ch10,&
&     '  for wf file unit',unitfile,', the value of nspinor_disk=',nspinor_disk,ch10,&
&     '  disagrees with argument nspinor=',nspinor,'. IO problem.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if

!  Check agreement with nband_k
   if (nband_k/=nband_disk) then
     write(message, '(4a,i5,a,i2,3a,i6,a,i6,2a,i6,a)' ) ch10,&
&     ' rdnpw : BUG -',ch10,&
&     '  At k point number',ikpt,', with spin polarization',isppol,',',ch10,&
&     '  for wf file unit',unitfile,', the value of nband_disk=',nband_disk,ch10,&
&     '  disagrees with argument nband_k=',nband_k,'. IO problem.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if

 else if(option==2)then

   npw_k=npw_disk
   nspinor=nspinor_disk
   nband_k=nband_disk

 else

   write(message, '(a,a,a,a,i5)' ) ch10,&
&   ' rdnpw : BUG -',ch10,&
&   '  Only option=0, 1 or 2 are allowed, while option=',option
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')

 end if

!Skip the next line (k+G)
 read(unitfile,IOSTAT=ierr)

 if(ierr/=0)then
   write(message, '(4a,i5,a,i5)' ) ch10,&
&   ' rdnpw : BUG -',ch10,&
&   '  Reading next record of disk file unit',unitfile,', gives iostat=',ierr
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

end subroutine rdnpw
!!***
