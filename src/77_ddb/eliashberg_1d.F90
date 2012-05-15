!{\src2tex{textfont=tt}}
!!****f* ABINIT/eliashberg_1d
!!
!! NAME
!! eliashberg_1d
!!
!! FUNCTION
!!  Solve the Eliashberg equations in the isotropic case
!!   First the linearized case, which allows the estimation of Tc
!!   then the full case which gives the gap as a function of temperature.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  a2f_1d = alpha^2F function averaged over the FS (only energy dependence)
!!  elph_ds = datastructure with phonon matrix elements
!!  gkk2 = gkk2 matrix elements on full FS grid for each phonon mode
!!  natom = number of atoms
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      eli_diag_m_1d,eli_lambda_1d,eli_z_1d,wrtout
!!
!! NOTES
!!  na2f = number of frequency points for alpha^2F function
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine eliashberg_1d(a2f_1d,elph_ds,mustar)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_elphon
 use m_io_tools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eliashberg_1d'
 use interfaces_14_hidewrite
 use interfaces_77_ddb, except_this_one => eliashberg_1d
!End of the abilint section

 implicit none

!Arguments ------------------------------------
  ! needed for phonon interpolation
!scalars
 real(dp),intent(in) :: mustar
 type(elph_type),intent(in) :: elph_ds
!arrays
 real(dp),intent(in) :: a2f_1d(elph_ds%na2f)

!Local variables-------------------------------
  ! for diagonalization of gammma matrix
  ! output variables for gtdyn9+phfrq3
!scalars
 integer :: iiter,imatsu
 integer :: maxiter,nmatsu,unit_del,unit_lam,unit_z
 real(dp) :: maxeigval,omega_cutoff
 real(dp) :: tc
 character(len=fnlen) :: fname
!arrays
 real(dp),allocatable :: delta_1d(:),lambda_1d(:),mm_1d(:,:),z_1d(:)

! *********************************************************************

 call wrtout(std_out,'Solving the 1-D Eliashberg equation (isotropic case)',"COLL")

 if (elph_ds%nsppol /= 1) then
   write(std_out,*) 'eliashberg_1d is not coded for nsppol > 1 yet'
   return
 end if

!maximum number of iterations to find T_c
 maxiter=30

!Fix nmatsu. Should add test in iiter loop to check if
!omega_cutoff is respected
 nmatsu = 50
!write(std_out,*) ' eliashberg_1d : nmatsu = ', nmatsu

 ABI_ALLOCATE(lambda_1d,(-nmatsu:nmatsu))
 ABI_ALLOCATE(z_1d,(-nmatsu:nmatsu))
 ABI_ALLOCATE(delta_1d,(-nmatsu:nmatsu))
 ABI_ALLOCATE(mm_1d,(-nmatsu:nmatsu,-nmatsu:nmatsu))

 unit_lam=get_unit()
 fname=trim(elph_ds%elph_base_name) // "_LAM"
 open (UNIT=unit_lam,FILE=fname,STATUS='REPLACE')
 unit_z=get_unit()
 fname=trim(elph_ds%elph_base_name) // "_Z"
 open (UNIT=unit_z,FILE=fname,STATUS='REPLACE')
 unit_del=get_unit()
 fname=trim(elph_ds%elph_base_name) // "_DEL"
 open (UNIT=unit_del,FILE=fname,STATUS='REPLACE')

!
!1) use linearized Eliashberg equation to find Tc
!! \sum_j \mathbf{M}_{ij} \Delta_j = \zeta \cdot \Delta_i $ $i,j = 1 .. n_{\mathrm{Matsubara}}$
!$\zeta = 1$ gives T$_c$ $\beta = \frac{1}{\mathrm{T}}$ $\omega_i = (2 i + 1) \pi \mathrm{T}$
!! \mathbf{M}_{ij} = \frac{\pi}{\beta} \frac{\lambda (\omega_i - \omega_j)}{Z (\omega_i)}$
!! Z (\omega_i) = 1 + \frac{\pi}{\beta \omega_i} \sum_j \lambda(\omega_i - \omega_j) \mathrm{sgn}(\omega_j)$
!

!initial guess for T$_c$ in Hartree (1Ha =3.067e5 K)
 tc = 0.0001
!
!big iterative loop
!
 do iiter=1,maxiter

   omega_cutoff = (two*nmatsu+one) * pi * tc

!  
!  calculate array of lambda values
!  
   call eli_lambda_1d (a2f_1d,elph_ds,lambda_1d,nmatsu,tc)
   write (unit_lam,'(a)') '#'
   write (unit_lam,'(a)') '# ABINIT package : lambda file'
   write (unit_lam,'(a)') '#'
   write (unit_lam,'(a,I10,a)') '# lambda_1d array containing 2*', nmatsu, '+1 Matsubara frequency points'
   write (unit_lam,'(a,E16.6,a,E16.6)') '#  from ', -omega_cutoff, ' to ', omega_cutoff
   write (unit_lam,'(a)') '#  lambda_1d is the frequency dependent coupling constant '
   write (unit_lam,'(a)') '#  in the Eliashberg equations '
   write (unit_lam,'(a)') '#'
   do imatsu=-nmatsu,nmatsu
     write (unit_lam,*) imatsu,lambda_1d(imatsu)
   end do
   write (unit_lam,*)

!  
!  calculate array of z values
!  
   call eli_z_1d (lambda_1d,nmatsu,z_1d)
   write (unit_z,'(a)') '#'
   write (unit_z,'(a)') '# ABINIT package : Z file'
   write (unit_z,'(a)') '#'
   write (unit_z,'(a,I10,a)') '# z_1d array containing 2*', nmatsu, '+1 Matsubara frequency points'
   write (unit_z,'(a,E16.6,a,E16.6)') '# from ', -omega_cutoff, ' to ', omega_cutoff
   write (unit_z,'(a)') '# z_1d is the renormalization factor in the Eliashberg equations'
   write (unit_z,'(a)') '#'
   do imatsu=-nmatsu,nmatsu
     write (unit_z,*) imatsu,z_1d(imatsu)
   end do

!  !
!  ! apply M matrix until a maximal eigenvalue is found.
!  !
!  call eli_m_iter_1d (delta_1d,lambda_1d,maxeigval,nmatsu,z_1d)

!  
!  diagonalize M brute forcefully
!  
   call eli_diag_m_1d(delta_1d,lambda_1d,maxeigval,mustar,nmatsu,tc,z_1d)

   write (unit_del,'(a)') '#'
   write (unit_del,'(a)') '# eliashberg_1d : delta_1d = '
   write (unit_del,'(a)') '#'
   write (unit_del,'(a,i6,a)') '# delta_1d array containing 2*', nmatsu, '+1 Matsubara frequency points'
   write (unit_z,'(a,E16.6,a,E16.6)') '# from ', -omega_cutoff, ' to ', omega_cutoff
   write (unit_z,'(a)') '# delta_1d is the gap function in the Eliashberg equations'
   write (unit_z,'(a)') '#'
   do imatsu=-nmatsu,nmatsu
     write (unit_del,*) imatsu,delta_1d(imatsu)
   end do
   write (unit_del,*)

!  if eigenvalue is < 1 increase T
!  else if eigenvalue is > 1 decrease T
!  if eigenvalue ~= 1 stop
!  
   if (abs(maxeigval-one) < tol8) then
     write(std_out,*) 'Eliashberg Tc found = ', tc, ' (Ha) = ', tc/kb_HaK, ' (K)'
     exit
   else if (maxeigval > 0.001_dp) then
     tc = tc * maxeigval
   else
     write(std_out,*) 'maxeigval is very small'
     tc = tc * 1000.0_dp
   end if


 end do
!end iiter do

 if (abs(maxeigval-one) > tol8) then
   write(std_out,*) 'eliashberg_1d : Tc not converged. ', maxeigval, ' /= 1'
   write(std_out,*) 'Eliashberg Tc nonetheless = ', tc, ' (Ha) = ', tc/kb_HaK, ' (K)'
 end if

 ABI_DEALLOCATE(lambda_1d)
 ABI_DEALLOCATE(z_1d)
 ABI_DEALLOCATE(delta_1d)
 ABI_DEALLOCATE(mm_1d)

 close (UNIT=unit_z)
 close (UNIT=unit_lam)
 close (UNIT=unit_del)

 write(std_out,*) ' eliashberg_1d : end '


end subroutine eliashberg_1d
!!***
