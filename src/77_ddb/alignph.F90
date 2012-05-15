!{\src2tex{textfont=tt}}
!!****f* ABINIT/alignph
!!
!! NAME
!! alignph
!!
!! FUNCTION
!! Construct linear combinations of the phonon eigendisplacements
!! of degenerate modes in order to align the mode effective charges
!! along the axes of the cartesian frame.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! displ(2,3*natom,3*natom)=
!! the displacements of atoms in cartesian coordinates.
!! The first index means either the real or the imaginary part,
!! The second index runs on the direction and the atoms displaced
!! The third index runs on the modes.
!! d2cart(2,3,mpert,3,mpert)=
!!  dynamical matrix, effective charges, dielectric tensor,....
!!  all in cartesian coordinates
!! mpert =maximum number of ipert
!! natom=number of atoms in unit cell
!! ntypat=number of types of atoms
!! phfrq(3*natom)=phonon frequencies (square root of the dynamical
!!  matrix eigenvalues, except if these are negative, and in this
!!  case, give minus the square root of the absolute value
!!  of the matrix eigenvalues). Hartree units.
!! typat(natom)=integer label of each type of atom (1,2,...)
!!
!! OUTPUT
!! displ(2,3*natom,3*natom)=
!! the displacements of atoms in cartesian coordinates.
!! The eigendisplacements of degenerate modes have been aligned along
!! the cartesian axes.
!!
!! PARENTS
!!      diel9
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine alignph(amu,displ,d2cart,mpert,natom,ntypat,phfrq,typat)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'alignph'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,ntypat
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: amu(ntypat),d2cart(2,3,mpert,3,mpert),phfrq(3*natom)
 real(dp),intent(inout) :: displ(2,3*natom,3*natom)

!Local variables -------------------------
!scalars
 integer :: i1,idir1,idir2,ii,imode,ipert1,jj,kk
 real(dp) :: c1,c2,c3,c4,c5,c6,dtm,mod_,theta
!arrays
 integer,allocatable :: deg(:)
 real(dp) :: coeff(3)
 real(dp),allocatable :: modez(:,:,:),oscstr(:,:,:),vec(:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)'alignph : enter'
!stop
!ENDDEBUG

 ABI_ALLOCATE(deg,(3*natom))
 deg(:) = 0

!Find degenerate modes

 imode = 0
 do while (imode < 3*natom)
   imode = imode + 1
   if (imode == 3*natom) then
     deg(imode) = 1
   else if (abs(phfrq(imode) - phfrq(imode+1)) > tol6) then
     deg(imode) = 1
   else
     deg(imode) = 2
     if (imode < 3*natom - 1) then
       if (abs(phfrq(imode) - phfrq(imode+2)) < tol6) then
         deg(imode) = 3
         imode = imode + 1
       end if
     end if
     imode = imode + 1
   end if
 end do

!Get the oscillator strength and mode effective charge for each mode
!In case of a degenerate mode align the vectors along
!the axes of the cartesian frame

 ABI_ALLOCATE(oscstr,(2,3,3*natom))
 ABI_ALLOCATE(modez,(2,3,3*natom))
 ABI_ALLOCATE(vec,(3,3*natom))

 do ii=1,2
   do imode=1,3*natom
     do idir2=1,3
       oscstr(ii,idir2,imode)=0.0_dp
       modez(ii,idir2,imode)=0.0_dp
       do idir1=1,3
         do ipert1=1,natom
           i1=idir1+(ipert1-1)*3
           oscstr(ii,idir2,imode)=oscstr(ii,idir2,imode)+&
&           displ(ii,i1,imode)*&
&           d2cart(1,idir1,ipert1,idir2,natom+2)
           modez(ii,idir2,imode)=modez(ii,idir2,imode)+&
&           displ(ii,i1,imode)*&
&           d2cart(1,idir1,ipert1,idir2,natom+2)*&
&           sqrt(amu(typat(ipert1))*amu_emass)
         end do
       end do
     end do
   end do
 end do

 imode = 1
 do while (imode <= 3*natom)

   if (deg(imode) == 2) then

     if (abs(modez(1,1,imode)) > tol8) then
       theta = atan(-1._dp*modez(1,1,imode+1)/modez(1,1,imode))
       vec(1,:) = displ(1,:,imode)
       vec(2,:) = displ(1,:,imode+1)
       displ(1,:,imode) = cos(theta)*vec(1,:) - sin(theta)*vec(2,:)
       displ(1,:,imode+1) = sin(theta)*vec(1,:) + cos(theta)*vec(2,:)
     end if

   else if (deg(imode) == 3) then

     vec(1,:) = displ(1,:,imode)
     vec(2,:) = displ(1,:,imode+1)
     vec(3,:) = displ(1,:,imode+2)
     do ii = 1,3
       coeff(:) = 0._dp
       if (ii == 1) then
         jj = 2 ; kk = 3
       else if (ii == 2) then
         jj = 1 ; kk = 3
       else
         jj = 1 ; kk = 2
       end if
       coeff(ii) = 1._dp
       c1 = modez(1,jj,imode+ii-1)
       c2 = modez(1,jj,imode+jj-1)
       c3 = modez(1,jj,imode+kk-1)
       c4 = modez(1,kk,imode+ii-1)
       c5 = modez(1,kk,imode+jj-1)
       c6 = modez(1,kk,imode+kk-1)
       dtm = c2*c6 - c3*c5
       if (abs(dtm) > tol8) then
         coeff(jj) = (c3*c4 - c1*c6)/dtm
         coeff(kk) = (c1*c5 - c2*c4)/dtm
       end if
       mod_ = sqrt(1._dp + coeff(jj)*coeff(jj) + coeff(kk)*coeff(kk))
       coeff(:) = coeff(:)/mod_
       displ(1,:,imode+ii-1) = coeff(1)*vec(1,:) + coeff(2)*vec(2,:) + &
&       coeff(3)*vec(3,:)
     end do

   end if

   imode = imode + deg(imode)

 end do



 ABI_DEALLOCATE(deg)
 ABI_DEALLOCATE(oscstr)
 ABI_DEALLOCATE(modez)
 ABI_DEALLOCATE(vec)

end subroutine alignph
!!***
