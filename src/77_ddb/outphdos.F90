!{\src2tex{textfont=tt}}
!!****f* ABINIT/outphdos
!! NAME
!! outphdos
!!
!! FUNCTION
!!  Print out phonon density of states
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!   deltaene = step on energy/frequency grid, in Hartree
!!   dos_phon = phonon DOS calculated on a grid
!!   enemin = minimal frequency
!!   enemax = maximal frequency
!!   filnam = file name for output to disk
!!   g2fsmear = smearing width
!!   nene = number of points on energy axis
!!   nqpt = number of q-points
!!   ntetra = number of tetrahedra, if tetrahedron interpolation is used
!!   telphint = flag for el-phonon interpolation method (to indicate Gaussian or tetrahedron integration)
!!   unit_phdos = unit for phonon DOS output
!!
!!
!! OUTPUT
!!  only write
!!
!! SIDE EFFECTS
!!
!! NOTES
!!   FIXME
!!   overcomplete inputs. Eliminate unit_phdos (just filnam) and deltaene (gotten from max-min/nene)
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


subroutine outphdos(deltaene,dos_phon,enemin,enemax,filnam,g2fsmear,nene,nqpt,ntetra,telphint,unit_phdos)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outphdos'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nene,nqpt,ntetra,telphint,unit_phdos
 character(len=fnlen),intent(in) :: filnam
 real(dp) :: deltaene,enemin,enemax,g2fsmear
!arrays
 real(dp) :: dos_phon(nene)

!Local variables-------------------------------
!scalars
 integer :: iomega,iost,step10
 real(dp) :: dos_effective,omega
 character(len=fnlen) :: outfile
 character(len=500) :: message
!arrays

! *************************************************************************

 outfile = trim(filnam) // '_PDS'
 write(message, '(3a)')ch10,&
& ' Will write phonon DOS in file ',trim(outfile)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message, '(4a)')ch10,&
& ' For checking purposes, write ten values in the present file.',ch10,&
& '       Index    Energy (in Ha)      DOS '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 open (unit=unit_phdos,file=outfile,status='replace',iostat=iost)
 if (iost /= 0) then
   write (message,'(3a)')' thmeig : ERROR- opening file ',trim(outfile),' as new'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 write (unit_phdos,'(a)') '#'
 write (unit_phdos,'(a)') '# ABINIT package : phonon DOS file'
 write (unit_phdos,'(a)') '#'
 write (unit_phdos,'(a,i10)') '#   Number of Qpoints integrated over : ', nqpt
 write (unit_phdos,'(a,i10)') '#   Number of energy points : ', nene
 write (unit_phdos,'(a,es16.6,a,es16.6,a)') '#   between omega_min = ', enemin, &
& ' Ha and omega_max = ', enemax, ' Ha'
 if(telphint==1)then
   write (unit_phdos,'(a,es16.6)') '#   The smearing width for gaussians is ', g2fsmear
 end if
 if(telphint==0)then
   write (unit_phdos,'(a,i10)') '#   Number of tetrahedrons', ntetra
 end if
 write (unit_phdos,'(a)') '#'
 write (unit_phdos,'(a)') '#      Index    Energy (in Ha)      DOS '

 omega = enemin
 do iomega=1,nene
   dos_effective=dos_phon(iomega)
   if(abs(dos_effective)<tol16)then
     dos_effective=zero
   end if
   step10=nene/10
   if(mod(iomega,step10)==1)write (std_out,'(i10,es18.6,es18.6)')iomega, omega, dos_effective
   if(mod(iomega,step10)==1)write (ab_out,'(i10,es18.6,es18.6)')iomega, omega, dos_effective
   write (unit_phdos, '(i10,es18.6,es18.6)')iomega, omega, dos_effective
   omega=omega+deltaene
 end do
 
 close (unit=unit_phdos)

end subroutine outphdos
!!***
