!{\src2tex{textfont=tt}}
!!****f* ABINIT/rwwan
!! NAME
!! rwwan
!!
!! FUNCTION
!! Reads (or write) the eigenvectors corresponding to the Wannier (interpolated) bands
!! from (to) file (depending on the input variable irwfl).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! irwfl= identifier for read/write status = -1/+1
!! iwf= identifier for the Wannier file
!! lambda(nqpt,nwnn,maxqsize,2)= interpolation coefficients of the Wannier bands
!! maxqsize= maximum size of the global window
!! natom= number of atoms per unit cell
!! nqpt= number of q points in the whole BZ
!! nwnn=no of wannier functions to be generated
!! wannvect(nqpt,nwnn,natom,3,2)= eigenvector of the Wannier states (if to write)
!!
!! OUTPUT
!! wannvect(nqpt,nwnn,natom,3,2)= eigenvector of the Wannier states (if to read)
!!
!! SIDE EFFECTS
!! qpoint(nqpt,3)= coordinates of the q points
!!
!! PARENTS
!!      lwf,wanvec
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine rwwan(irwfl,iwf,natom,nqpt,nwnn,qpoint,wannvect)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rwwan'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: irwfl,iwf,natom,nqpt,nwnn
!arrays
 real(dp),intent(inout) :: qpoint(nqpt,3),wannvect(nqpt,nwnn,natom,3,2)

!Local variables-------------------------------
!scalars
 integer :: iiqpt,iiwann,iqpt,iwann,tao
 character(len=500) :: message

!******************************************************************
!BEGIN EXECUTABLE SECTION

!DEBUG
!write(std_out,*) ' rwwan : enter'
!ENDDEBUG

!write the Wannier information
 if (irwfl==1) then
   do iqpt=1,nqpt                        ! loop over Q points
     write(iwf,'(i4,3f10.5)') iqpt,qpoint(iqpt,:)
     do iwann=1,nwnn                         ! loop over Wannier states
       write(iwf,'(i4)') iwann
       do tao=1,natom                     ! loop over atoms
         write(iwf,'(6f12.7)') wannvect(iqpt,iwann,tao,1,1),wannvect(iqpt,iwann,tao,2,1),wannvect(iqpt,iwann,tao,3,1),&
&         wannvect(iqpt,iwann,tao,1,2),wannvect(iqpt,iwann,tao,2,2),wannvect(iqpt,iwann,tao,3,2)
       end do                              !atoms
     end do                                !Wannier states
   end do                                 !q points

!  read the Wannier information
 elseif (irwfl==-1) then
   do iqpt=1,nqpt                        ! loop over Q points
     read(iwf,'(i4,3f10.5)') iiqpt,qpoint(iqpt,:)
     do iwann=1,nwnn                         ! loop over Wannier states
       read(iwf,'(i4)') iiwann
       do tao=1,natom                     ! loop over atoms
         read(iwf,'(6f12.7)') wannvect(iqpt,iwann,tao,1,1),wannvect(iqpt,iwann,tao,2,1),wannvect(iqpt,iwann,tao,3,1),&
&         wannvect(iqpt,iwann,tao,1,2),wannvect(iqpt,iwann,tao,2,2),wannvect(iqpt,iwann,tao,3,2)
       end do                              !atoms
     end do                                !Wannier states
   end do                                 !q points

!  it is not possible => BUG
 else
   write(message, '(a,a,a,i4,a,a,a)' )&
&   ' wanwec : BUG -',ch10,&
&   '  irwfl is ',irwfl,' and it should be +1 or -1 .',ch10,&
&   '  Action : contact the ABINIT team.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if



!ENDDEBUG
!write(std_out,*) ' rwwan : exit'
!ENDDEBUG

end subroutine rwwan
!!***
