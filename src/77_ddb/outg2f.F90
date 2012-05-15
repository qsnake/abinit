!{\src2tex{textfont=tt}}
!!****f* ABINIT/outg2f
!! NAME
!! outg2f
!!
!! FUNCTION
!!  Output g2f function to file. FIXME: Paul, please explain what g2f is.
!!  Probably a variant on the Eliashberg spectral function a2F
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2012 ABINIT group (PB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  only write
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      thmeig
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine outg2f(deltaene,enemin,enemax,filnam,g2f,g2fsmear,kpnt,mband,nene,nkpt,nqpt,ntetra,telphint,unit_g2f)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outg2f'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nene,nkpt,nqpt,ntetra,telphint,unit_g2f
 character(len=fnlen),intent(in) :: filnam
 real(dp) :: deltaene,enemin,enemax,g2fsmear
!arrays
 real(dp) :: g2f(mband,nkpt,nene),kpnt(3,nkpt,nqpt)

!Local variables-------------------------------
!scalars
 integer :: iband,ikpt,iomega,iost
 real(dp) :: omega
 character(len=fnlen) :: outfile
 character(len=500) :: message
!arrays

! *************************************************************************

!output the g2f
 outfile = trim(filnam) // '_G2F'
 open (unit=unit_g2f,file=outfile,status='unknown',iostat=iost)
 if (iost /= 0) then
   write (message,'(3a)')' thmeig : ERROR- opening file ',trim(outfile),' as new'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 write(std_out,*) ' g2f function'
 write (unit_g2f,'(a)') '#'
 write (unit_g2f,'(a)') '# ABINIT package : g2f file'
 write (unit_g2f,'(a)') '#'
 write (unit_g2f,'(a,I10)') '#     number of qpoints integrated over : ', nqpt
 write (unit_g2f,'(a,I10)') '#     number of energy points : ', nene
 write (unit_g2f,'(a,E16.6,a,E16.6,a)') '#       between omega_min = ', enemin, &
& ' Ha and omega_max = ', enemax, ' Ha'
 if(telphint==1)then
   write (unit_g2f,'(a,E16.6)') '#   and the smearing width for gaussians is ', g2fsmear
   write (unit_g2f,'(a)') '#'
 end if
 if(telphint==0)then
   write (unit_g2f,'(a,I10)') '#   number of tetrahedrons', ntetra
   write (unit_g2f,'(a)') '#'
 end if

!Write only the a2f function for the first K point
!ikpt=1
 do ikpt=1,nkpt
   write(unit_g2f,'(a,3es16.8)')' Kpt :', kpnt(:,ikpt,1)
   do iband=1,mband
     write(unit_g2f,*) 'band :', iband
     omega = enemin
     do iomega=1,nene
       write (unit_g2f,*) omega*Ha_eV*1000, g2f(iband, ikpt,iomega)
       omega=omega+deltaene
     end do
   end do
 end do

 close (unit=unit_g2f)

end subroutine outg2f
!!***



