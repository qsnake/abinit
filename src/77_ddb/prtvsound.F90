!!****f* ABINIT/prtvsound
!!
!! NAME
!! prtvsound
!!
!! FUNCTION
!!  From the frequencies for acoustic modes at small q, estimate speed of sound and Debye temperature
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! eigvec(2,3*natom,3*natom) = phonon eigenvectors at present q-point
!! gmet(3,3) = metric tensor in reciprocal space.
!! natom = number of atoms in the unit cell
!! phfrq(3*natom) = phonon frequencies at present q-point
!! qphon(3) = phonon q-point
!! ucvol = unit cell volume
!!
!! OUTPUT
!!
!! NOTES
!! 
!! PARENTS
!!      mkphbs
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prtvsound(eigvec, gmet, natom, phfrq, qphon, ucvol)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtvsound'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------

 integer, intent(in) :: natom
 real(dp), intent(in) :: ucvol
 real(dp), intent(in) :: gmet(3,3)
 real(dp), intent(in) :: qphon(3)
 real(dp), intent(in) :: phfrq(3*natom)
 real(dp), intent(in) :: eigvec(2,3*natom,3*natom)

!Local variables -------------------------

 integer :: imode, iatom, isacoustic
 character(len=500) :: message
 real(dp) :: qnormcart, speedofsound, tdebye
 real(dp) :: qtmp(3)

! *********************************************************************

 do imode = 1, 3*natom
   
!  check if this mode is acoustic like: scalar product of all displacement vectors are collinear
   isacoustic = 1
   do iatom = 2, natom
     if (sum(eigvec(1,(iatom-1)*3+1:(iatom-1)*3+3, imode)*eigvec(1,1:3, imode)) + &
&     sum(eigvec(2,(iatom-1)*3+1:(iatom-1)*3+3, imode)*eigvec(2,1:3, imode)) < zero) isacoustic = 0
   end do
   if (isacoustic == 0) cycle

   write (message, '(a,I6,a,3F12.4)') ' Found acoustic mode ', imode, &
&   ' for |q| in red coord < 0.25 ; q = ', qphon
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   
   qtmp = matmul(gmet, qphon)
   qnormcart = two * pi * sqrt(sum(qphon*qtmp))
   speedofsound = phfrq(imode) / qnormcart

!  from phonon frequency, estimate speed of sound by linear interpolation from Gamma
   write (message, '(2a,a,E20.10,a,a,F20.5)') &
&   ' Speed of sound for this q and mode:',ch10,&
&   '   in atomic units: ', speedofsound, ch10,&
&   '   in SI units m/s: ', speedofsound * Bohr_Ang * 1.d-10 / Time_Sec
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')


!  also estimate partial Debye temperature, = energy if this band went to zone edge
   tdebye = speedofsound * pi * (six / pi / ucvol)**(third)
   write (message, '(2a,a,E20.10,a,a,F20.5)') &
&   ' Partial Debye temperature for this q and mode:',ch10,&
&   '   in atomic units: ', tdebye, ch10,&
&   '   in SI units K  : ', tdebye * Ha_K
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   write (message, '(a)') ''
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 end do

end subroutine prtvsound
!!***
