!!****f* ABINIT/scphon_free_energy
!! NAME
!! scphon_free_energy
!!
!! FUNCTION
!! Calculate the phonon Free energy, from the input Density of States.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! istep= number of the present iteration of the SC phonon calculations, for
!!   printing purposes
!! phonon_dos=array containing the phonon DOS
!! scphon_temp= phononic temperature in Ha
!!
!! OUTPUT
!! free_energy= value of the Free energy
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      scphon
!!
!! CHILDREN
!!      simpson_int,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon_free_energy(free_energy,istep,t_phonon_dos,scphon_temp)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_phdos

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_free_energy'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep
 real(dp),intent(in) :: scphon_temp
 real(dp),intent(out) :: free_energy
 type(phonon_dos_type),intent(in) :: t_phonon_dos

!Local variables-------------------------------
!scalars
 integer :: ifreq
 real(dp) :: freq,inv_temp
 character(len=500) :: message
!arrays
 real(dp),allocatable :: free_energy_int(:),free_energy_integrand(:)

! *************************************************************************

 inv_temp = one / scphon_temp

!integrate free energy from DOS
 ABI_ALLOCATE(free_energy_integrand,(t_phonon_dos%nomega))
 ABI_ALLOCATE(free_energy_int,(t_phonon_dos%nomega))
 free_energy_integrand=zero
 do ifreq=1,t_phonon_dos%nomega
   freq=t_phonon_dos%omega_min+dble(ifreq-1)*&
&   (t_phonon_dos%omega_max-t_phonon_dos%omega_min)/dble(t_phonon_dos%nomega-1)
!  0 point energy
   free_energy_integrand(ifreq)=t_phonon_dos%phdos(ifreq)*freq*half
!  rest of the free energy
   if ((freq*inv_temp) < 30._dp .and. (freq*inv_temp) > tol12) then
     free_energy_integrand(ifreq)=free_energy_integrand(ifreq) &
&     + t_phonon_dos%phdos(ifreq)*scphon_temp*log(one-exp(-freq*inv_temp))
!    else if((freq*inv_temp) < tol12) then
!    free_energy_integrand(ifreq)=free_energy_integrand(ifreq) &
!    &     + t_phonon_dos%phdos(ifreq)*scphon_temp*
   end if
 end do
 free_energy_int=zero
 call simpson_int(t_phonon_dos%nomega,t_phonon_dos%omega_step,free_energy_integrand,free_energy_int)
 free_energy = free_energy_int(t_phonon_dos%nomega)


 write (message,'(a,I6,a,E10.3,a)') ' Free energy at iteration ', istep, &
& ' and temperature T= ',scphon_temp/kb_HaK, ' K is:'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')
 write (message,'(a,i6,E20.8)') 'FREEEN ',  istep, free_energy
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 ABI_DEALLOCATE(free_energy_int)
 ABI_DEALLOCATE(free_energy_integrand)

end subroutine scphon_free_energy
!!***


