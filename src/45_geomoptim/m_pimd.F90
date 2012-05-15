!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pimd
!! NAME
!!  m_pimd
!!
!! FUNCTION
!!  This module provides several routines and datatypes for the
!!  Path-Integral Molecular Dynamics (PIMD) implementation.
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (GG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_pimd

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors
 use m_random_zbq

 implicit none

 private

!public procedures
 public :: pimd_print
 public :: pimd_is_restart
 public :: pimd_temperature
 public :: pimd_initvel
 public :: pimd_langevin_random
 public :: pimd_langevin_random_bar
 public :: pimd_langevin_random_init
 public :: pimd_energies
 public :: pimd_forces
 public :: pimd_langevin_forces
 public :: pimd_nosehoover_forces
 public :: pimd_stresses
 public :: pimd_diff_stress
 public :: pimd_predict_taylor
 public :: pimd_predict_verlet
 public :: pimd_nosehoover_propagate
 public :: pimd_coord_transform
 public :: pimd_force_transform
 public :: pimd_mass_spring
!!***

!!****t* m_pimd/pimd_type
!! NAME
!! pimd_type
!!
!! FUNCTION
!! Datatype with the variables required to perform PIMD
!!
!! NOTES
!!
!! SOURCE

 type,public :: pimd_type
! Scalars
  integer,pointer :: irandom
  integer,pointer :: natcon
  integer,pointer :: natom
  integer,pointer :: natfix
  integer,pointer :: natfixx
  integer,pointer :: natfixy
  integer,pointer :: natfixz
  integer,pointer :: nconeq
  integer,pointer :: nnos
  integer,pointer :: ntypat
  integer,pointer :: optcell
  integer,pointer :: pitransform
  real(dp),pointer :: vis
  real(dp),pointer :: bmass
  real(dp),pointer :: dtion
  real(dp),pointer :: friction
! Arrays
  integer,pointer :: iatcon(:)
  integer,pointer :: iatfix(:)
  integer,pointer :: iatfixx(:)
  integer,pointer :: iatfixy(:)
  integer,pointer :: iatfixz(:)
  integer,pointer :: typat(:)
  real(dp),pointer :: amu(:)
  real(dp),pointer :: mdtemp(:)
  real(dp),pointer :: pimass(:)
  real(dp),pointer :: qmass(:)
  real(dp),pointer :: strtarget(:)
  real(dp),pointer :: wtatcon(:)
 end type pimd_type
!!***

CONTAINS

!===========================================================
!!***

!!****f* m_pimd/pimd_is_restart
!! NAME
!!  pimd_is_restart
!!
!! FUNCTION
!!  Determine whether this is a PIMD restart or not:
!!  test on value of velocities and corresponding temperature
!!
!! INPUTS
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  vel(3,natom,nimage)=velocities for each image of the cell
!!
!! OUTPUT
!!  pimd_is_restart=1 if temperature is not zero
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function pimd_is_restart(mass,vel)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_is_restart'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: pimd_is_restart
!arrays
 real(dp),intent(in) :: mass(:,:),vel(:,:,:)
!Local variables-------------------------------
!scalars
 real(dp),parameter :: zero_temp=tol7
!arrays

!************************************************************************

 pimd_is_restart=0
 if (pimd_temperature(mass,vel)>zero_temp) pimd_is_restart=1

end function pimd_is_restart
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_temperature
!! NAME
!!  pimd_temperature
!!
!! FUNCTION
!!  Compute temperature from velocities and masses
!!
!! INPUTS
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  vel(3,natom,nimage)=velocities for each image of the cell
!!
!! OUTPUT
!!  pimd_temperature=temperature (from all images of the cell)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function pimd_temperature(mass,vel)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_temperature'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: pimd_temperature
!arrays
 real(dp),intent(in) :: mass(:,:),vel(:,:,:)
!Local variables-------------------------------
!scalars
 integer :: iatom,idir,iimage,imass,natom,natom_mass,ndir,nimage,nmass
 real(dp) :: v2
 character(len=500) :: msg
!arrays

!************************************************************************

 ndir=size(vel,1);natom=size(vel,2);nimage=size(vel,3)
 natom_mass=size(mass,1);nmass=size(mass,2)
 if (ndir/=3.or.natom<=0.or.nimage<=0) then
   msg='Wrong sizes for vel array !'
   MSG_BUG(msg)
 end if
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=nimage)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

 v2=zero
 do iimage=1,nimage
   imass=min(nmass,iimage)
   do iatom=1,natom
     do idir=1,3
       v2=v2+vel(idir,iatom,iimage)*vel(idir,iatom,iimage)*mass(iatom,imass)
     end do
   end do
 end do
 pimd_temperature=v2/(dble(3*natom*nimage)*kb_HaK)

end function pimd_temperature
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_print
!! NAME
!!  pimd_print
!!
!! FUNCTION
!!  Print out results related to PIMD (for given time step)
!!
!! INPUTS
!!  eharm=harmonic energy
!!  eharm_virial=harmonic energy from virial
!!  epot=potential energy
!!  forces(3,natom,trotter)=forces on atoms in each cell
!!  irestart=1 if this is a restart
!!  itimimage=index of time step
!!  natom=number of atoms
!!  prtstress=flag for stress tensor printing
!!  prtvolimg=printing volume
!!  stress(3,3)=stress tensor
!!  temperature1,temperature2=temperatures at t and t+dt
!!  trotter=Trotter number
!!  vel(3,natom,trotter)=velocities of atoms in in each cell
!!  xcart(3,natom,trotter)=cartesian coordinates of atoms in in each cell
!!  xred(3,natom,trotter)=reduced coordinates of atoms in in each cell
!!
!! OUTPUT
!!  -- only printing --
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_npt
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_print(eharm,eharm_virial,epot,forces,irestart,itimimage,natom,prtstress,prtvolimg,&
&                     stress,temperature1,temperature2,trotter,vel,xcart,xred)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: irestart,itimimage,natom,prtstress,prtvolimg,trotter
 real(dp),intent(in) :: eharm,eharm_virial,epot,temperature1,temperature2
!arrays
 real(dp),intent(in) :: forces(3,natom,trotter),stress(3,3),vel(3,natom,trotter)
 real(dp),intent(in) :: xcart(3,natom,trotter),xred(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage
 character(len=500) :: msg
!arrays
 real(dp) :: forcetot(3)
 real(dp),allocatable :: centroid(:,:),qudeloc(:)

!************************************************************************

!Temperature
 if(itimimage==1)then
   if(irestart==0) then
     write(msg,'(2a)') ch10,' This is a PIMD calculation from scratch'
   else if (irestart==1) then
     write(msg,'(2a)') ch10,' This is a RESTART calculation'
   end if
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a,f12.5,a)') &
&   ' In the initial configuration, the temperature is ',temperature1,' K'
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if
 write(msg,'(2a,i5,a,f12.5,a)') ch10,&
&  ' At PIMD time step ',itimimage,', the temperature is',temperature2,' K'
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

!Energies
 write(msg,'(4a,2(f18.9,a,a,a),f18.9,a)') ch10,&
   ' Energy components:',ch10, &
&  '   Harmonic  energy          =',eharm ,' Ha',ch10, &
&  '   Harmonic  energy (virial) =',eharm_virial,' Ha',ch10, &
&  '   Potential energy          =',epot  ,' Ha'
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

!Stress tensor and pressure
 write(msg,'(2a,3(2a,3f18.9))') ch10,&
&   ' Stress tensor from virial theorem (Ha/Bohr^3):',ch10, &
&   '   ',stress(1,1),stress(1,2),stress(1,3),ch10, &
&   '   ',stress(2,1),stress(2,2),stress(2,3),ch10, &
&   '   ',stress(3,1),stress(3,2),stress(3,3)
 if (prtstress==1) then
   call wrtout(ab_out,msg,'COLL')
 end if
 call wrtout(std_out,msg,'COLL')
 write(msg,'(a,f18.9,a)') ' Pressure=', &
&  -third*(stress(1,1)+stress(2,2)+stress(3,3))*HaBohr3_GPa,' GPa'
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

!Total force
 forcetot=zero
 do iimage=1,trotter
   do iatom=1,natom
     do ii=1,3
       forcetot(ii)=forcetot(ii)+forces(ii,iatom,iimage)
     end do
   end do
 end do
 write(msg,'(2a,3f18.10)') ch10,' Total force=',forcetot(1:3)
 call wrtout(std_out,msg,'COLL')

!Positions
 write(msg,'(2a)') ch10,' Old atomic positions:'
 call wrtout(std_out,msg,'COLL')
 do iimage=1,trotter
   select case(iimage)
   case(1)
    write(msg,'(a)') ' xred'
   case(2,3,4,5,6,7,8,9)
     write(msg,'(a,i1,a)') ' xred_',iimage,'img'
   case default
     write(msg,'(a,i2,a)') ' xred_',iimage,'img'
   end select
   call wrtout(std_out,msg,'COLL')
   do iatom=1,natom
     write(msg,'(3f18.10)') xred(1:3,iatom,iimage)
     call wrtout(std_out,msg,'COLL')
   end do
 end do

!Velocities
 write(msg,'(2a)') ch10,' Velocities:'
 call wrtout(std_out,msg,'COLL')
 do iimage=1,trotter
   select case(iimage)
   case(1)
     write(msg,'(a)') ' vel'
   case(2,3,4,5,6,7,8,9)
     write(msg,'(a,i1,a)') ' vel_',iimage,'img'
   case default
     write(msg,'(a,i2,a)') ' vel_',iimage,'img'
   end select
   call wrtout(std_out,msg,'COLL')
   do iatom=1,natom
     write(msg,'(3f18.10)') vel(1:3,iatom,iimage)
     call wrtout(std_out,msg,'COLL')
   end do
 end do

!Centroids and wave-packet spatial spreads
 ABI_ALLOCATE(centroid,(3,natom))
 ABI_ALLOCATE(qudeloc,(natom))
 centroid=zero;qudeloc=zero
 do iimage=1,trotter
   do iatom=1,natom
     do ii=1,3
       centroid(ii,iatom)=centroid(ii,iatom)+xcart(ii,iatom,iimage)
     end do
   end do
 end do
 centroid=centroid/dble(trotter)
 do iimage=1,trotter
   do iatom=1,natom
     do ii=1,3
       qudeloc(iatom)=qudeloc(iatom)+((xcart(ii,iatom,iimage)-centroid(ii,iatom))**2)
     end do
   end do
 end do
 qudeloc=sqrt(qudeloc/dble(trotter))
 write(msg,'(4a)') ch10,' Centroids and wave-packet spatial spreads (cart. coord.):',ch10,&
&  ' iat        centroid_x        centroid_y        centroid_y    spatial_spread'
 call wrtout(std_out,msg,'COLL')
 do iatom=1,natom
   write(msg,'(i4,4f18.10)') iatom,centroid(1:3,iatom),qudeloc(iatom)
   call wrtout(std_out,msg,'COLL')
 end do
 ABI_DEALLOCATE(centroid)
 ABI_DEALLOCATE(qudeloc)

!Fake statement
 return;ii=prtvolimg

end subroutine pimd_print
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_initvel
!! NAME
!!  pimd_initvel
!!
!! FUNCTION
!!  Initialize velocities for PIMD with a gaussian distribution
!!  fixing the center of mass
!!
!! INPUTS
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  temperature=temperature used to define velocities
!!  trotter=Trotter number
!!
!! OUTPUT
!!  vel(3,natom,trotter)=velocities of atoms in in each cell
!!
!! SIDE EFFECTS
!!  iseed=seed for random number generator
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_npt
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_initvel(iseed,mass,natom,temperature,trotter,vel)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_initvel'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,trotter
 integer,intent(inout) :: iseed
 real(dp),intent(in) :: temperature
!arrays
 real(dp),intent(in) :: mass(:,:)
 real(dp),intent(out) :: vel(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,natom_mass,nmass
 real(dp) :: mtot,rescale_vel
 character(len=500) :: msg
!arrays
 real(dp) :: mvini(3)

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

 mtot=sum(mass(1:natom,1:nmass))
 if (nmass==1) mtot=mtot*trotter

 do iimage=1,trotter
   imass=min(nmass,iimage)
   do iatom=1,natom
     do ii=1,3
       vel(ii,iatom,iimage)=sqrt(kb_HaK*temperature/mass(iatom,imass))*cos(two_pi*uniformrandom(iseed))
       vel(ii,iatom,iimage)=vel(ii,iatom,iimage)*sqrt(-two*log(uniformrandom(iseed)))
     end do
   end do
 end do

!Make sure that the (sum of m_i v_i) at step zero is zero
 mvini=zero
 do iimage=1,trotter
   imass=min(nmass,iimage)
   do iatom=1,natom
     do ii=1,3
      mvini(ii)=mvini(ii)+mass(iatom,imass)*vel(ii,iatom,iimage)
     end do
   end do
 end do
 do iimage=1,trotter
   do iatom=1,natom
     do ii=1,3
       vel(ii,iatom,iimage)=vel(ii,iatom,iimage)-(mvini(ii)/mtot)
     end do
   end do
 end do
!at this point (sum of m_i v_i) at step zero is zero

!Now rescale the velocities to give the exact temperature
 rescale_vel=sqrt(temperature/pimd_temperature(mass,vel))
 vel(:,:,:)=vel(:,:,:)*rescale_vel

end subroutine pimd_initvel
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_langevin_random
!! NAME
!!  pimd_langevin_random
!!
!! FUNCTION
!!  Generate a set of random numbers to be used for PIMD Langevin algorithm
!!
!! INPUTS
!!  irandom=option for random number generator:
!!          1:uniform random routine provided within Abinit package
!!          2:Fortran 90 random number generator
!!          3:ZBQLU01 non deterministic random number generator
!!  langev(natom)=Langevin factors
!!  mass(natom)=masses of atoms
!!  natom=number of atoms
!!  trotter=Trotter number
!!  zeroforce=flag; if 1 keep sum of forces equal to zero
!!
!! OUTPUT
!!  alea(3,natom,trotter)=set of random numbers
!!
!! SIDE EFFECTS
!!  iseed=seed for random number generator (used only if irandom=1)
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_langevin_random(alea,irandom,iseed,langev,mass,natom,trotter,zeroforce)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_langevin_random'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: irandom,natom,trotter,zeroforce
 integer,intent(inout) :: iseed
!arrays
 real(dp),intent(in) :: langev(natom),mass(natom)
 real(dp),intent(out) :: alea(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage
 real(dp) :: mtot,r1,r2
!arrays
 real(dp) :: total(3)

!************************************************************************

!Draw random numbers
 do iimage=1,trotter
   do iatom=1,natom
     do ii=1,3
       select case(irandom)
       case(1)
         r1=uniformrandom(iseed)
         r2=uniformrandom(iseed)
       case(2)
         call random_number(r1)
         call random_number(r2)
       case(3)
         r1=ZBQLU01(zero)
         r2=ZBQLU01(zero)
       end select
       alea(ii,iatom,iimage)= cos(two_pi*r1)*sqrt(-log(r2)*two)
     end do
   end do
 end do

!Make sure that the sum of random forces is zero
 if(zeroforce==1)then
   total=zero
   mtot=sum(mass(1:natom))*trotter
   do iimage=1,trotter
     do iatom=1,natom
       do ii=1,3
         total(ii)=total(ii)+langev(iatom)*alea(ii,iatom,iimage)
       end do
     end do
   end do
   do iimage=1,trotter
     do iatom=1,natom
       do ii=1,3
         alea(ii,iatom,iimage)= alea(ii,iatom,iimage)-(total(ii)*mass(iatom))/(langev(iatom)*mtot)
       end do
     end do
   end do
!  now random forces have been rescaled so that their sum is zero
 end if

end subroutine pimd_langevin_random
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_langevin_random_bar
!! NAME
!!  pimd_langevin_random_bar
!!
!! FUNCTION
!!  Generate a set of random numbers to be used for the barostat of PIMD Langevin algorithm
!!
!! INPUTS
!!  irandom=option for random number generator:
!!          1:uniform random routine provided within Abinit package
!!          2:Fortran 90 random number generator
!!          3:ZBQLU01 non deterministic random number generator
!!
!! OUTPUT
!!  alea_bar(3,3)=set of random numbers
!!
!! SIDE EFFECTS
!!  iseed=seed for random number generator (used only if irandom=1)
!!
!! PARENTS
!!      pimd_langevin_npt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_langevin_random_bar(alea_bar,irandom,iseed)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_langevin_random_bar'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: irandom
 integer,intent(inout) :: iseed
!arrays
 real(dp),intent(out) :: alea_bar(3,3)
!Local variables-------------------------------
!scalars
 integer :: ii,jj
 real(dp) :: r1,r2
!arrays

!************************************************************************

!Draw random numbers
 do ii=1,3
   do jj=1,3
     select case(irandom)
     case(1)
       r1=uniformrandom(iseed)
       r2=uniformrandom(iseed)
     case(2)
       call random_number(r1)
       call random_number(r2)
     case(3)
       r1=ZBQLU01(zero)
       r2=ZBQLU01(zero)
     end select
     alea_bar(ii,jj)= cos(two_pi*r1)*sqrt(-log(r2)*two)
   end do
 end do

!Symmetrize
 alea_bar(1,2)=half*(alea_bar(1,2)+alea_bar(2,1))
 alea_bar(1,3)=half*(alea_bar(1,3)+alea_bar(3,1))
 alea_bar(2,3)=half*(alea_bar(2,3)+alea_bar(3,2))
 alea_bar(2,1)=alea_bar(1,2)
 alea_bar(3,1)=alea_bar(1,3)
 alea_bar(3,2)=alea_bar(2,3)

end subroutine pimd_langevin_random_bar
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_langevin_random_init
!! NAME
!!  pimd_langevin_random_init
!!
!! FUNCTION
!!  Initialize random number generator to be used for PIMD Langevin algorithm
!!
!! INPUTS
!!  irandom=option for random number generator:
!!          1:uniform random routine provided within Abinit package
!!          2:Fortran 90 random number generator
!!          3:ZBQLU01 non deterministic random number generator
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  iseed=seed for random number generator (used only if irandom=1)
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_langevin_random_init(irandom,iseed)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_langevin_random_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: irandom
 integer,intent(inout) :: iseed
!arrays
!Local variables-------------------------------
!scalars
!arrays

!************************************************************************

 if (irandom==3) then
   call ZBQLINI(0)
 end if
 
!Fake statement
 return;if (.false.) iseed=zero

end subroutine pimd_langevin_random_init
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_energies
!! NAME
!!  pimd_energies
!!
!! FUNCTION
!!  In the case od PIMD, compute the several contribution to total energy
!!
!! INPUTS
!!  etotal_img(trotter)= energy (from DFT) for each cell
!!  forces(3,natom,trotter)=forces (from DFT) on atoms in each cell
!!  natom=number of atoms
!!  spring(natom,spring_dim)=spring constants (spring_dim=1 or trotter)
!!  trotter=Trotter number
!!  xcart(3,natom,trotter)=cartesian coordinates of atoms in each cell at t
!!
!! OUTPUT
!!  eharm       =harmonic energy
!!  eharm_virial=harmonic energy from virial estimator
!!  epot        =potential energy
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_npt
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_energies(eharm,eharm_virial,epot,etotal_img,forces,natom,spring,trotter,xcart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_energies'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,trotter
 real(dp),intent(out) :: eharm,eharm_virial,epot
!arrays
 real(dp),intent(in) :: etotal_img(trotter),forces(3,natom,trotter)
 real(dp),intent(in) :: xcart(3,natom,trotter)
 real(dp),intent(in) :: spring(:,:)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,iimagep,ispring,natom_spring,nspring
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: centroid(:,:)

!************************************************************************

 natom_spring=size(spring,1);nspring=size(spring,2)
 if (natom/=natom_spring.or.(nspring/=1.and.nspring/=trotter)) then
   msg='Wrong dimensions for array spring !'
   MSG_BUG(msg)
 end if

!Compute the centroid
 ABI_ALLOCATE(centroid,(3,natom))
 centroid=zero
 do iimage=1,trotter
   do iatom=1,natom
     do ii=1,3
       centroid(ii,iatom)=centroid(ii,iatom)+xcart(ii,iatom,iimage)
     end do
   end do
 end do
 centroid=centroid/dble(trotter)

!Potential energy
 epot=sum(etotal_img(1:trotter))/dble(trotter)

!Harmonic energy
 eharm=zero
 do iimage=1,trotter
   ispring=min(nspring,iimage)
   iimagep=iimage+1;if(iimage==trotter)iimagep=1
   do iatom=1,natom
     do ii=1,3
       eharm=eharm+half*spring(iatom,ispring)*((xcart(ii,iatom,iimagep)-xcart(ii,iatom,iimage))**2)
     end do
   end do
 end do

!Harmonic energy from virial estimator
 eharm_virial=zero
 do iimage=1,trotter
   do iatom=1,natom
     do ii=1,3
       eharm_virial=eharm_virial-(xcart(ii,iatom,iimage)-centroid(ii,iatom)) &
    &              *forces(ii,iatom,iimage)
     end do
   end do
 end do
 eharm_virial=eharm_virial/dble(two*trotter)

 ABI_DEALLOCATE(centroid)

end subroutine pimd_energies
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_forces
!! NAME
!!  pimd_forces
!!
!! FUNCTION
!!  Modify forces in order to take into account PIMD contribution
!!
!! INPUTS
!!  natom=number of atoms
!!  spring(natom,spring_dim)=spring constants (spring_dim=1 or trotter)
!!  transform=coordinate transformation:
!!            0: no tranformation
!!            1: normal mode transformation
!!            2: staging transformation
!!  trotter=Trotter number
!!  xcart(3,natom,trotter)=cartesian coordinates of atoms in each cell
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  forces(3,natom,trotter)=
!!    at input:  forces from electronic calculation
!!    at output: forces from electronic calculation + quantum spring contribution
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_forces(forces,natom,spring,transform,trotter,xcart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_forces'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,transform,trotter
!arrays
 real(dp),intent(in) :: xcart(3,natom,trotter)
 real(dp),intent(in) :: spring(:,:)
 real(dp),intent(inout) :: forces(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,iimagem,iimagep,ispring,natom_spring,nspring
 character(len=500) :: msg
!arrays

!************************************************************************

 natom_spring=size(spring,1);nspring=size(spring,2)
 if (natom/=natom_spring.or.(nspring/=1.and.nspring/=trotter)) then
   msg='Wrong dimensions for array spring !'
   MSG_BUG(msg)
 end if

 if (transform==0) then
   do iimage=1,trotter
     ispring=min(nspring,iimage)
     iimagep=iimage+1; iimagem=iimage-1
     if(iimage==trotter) iimagep=1
     if(iimage==1)       iimagem=trotter
     do iatom=1,natom
       do ii=1,3
         forces(ii,iatom,iimage)= &
&               forces(ii,iatom,iimage)/dble(trotter) &
&             - spring(iatom,ispring)*(two*xcart(ii,iatom,iimage)-xcart(ii,iatom,iimagem) &
&                                                                -xcart(ii,iatom,iimagep))
       end do
     end do
   end do

 else
   do iimage=1,trotter
     ispring=min(nspring,iimage)
     iimagep=iimage+1; iimagem=iimage-1
     if(iimage==trotter) iimagep=1
     if(iimage==1)       iimagem=trotter
     do iatom=1,natom
       do ii=1,3
         forces(ii,iatom,iimage)= &
&               forces(ii,iatom,iimage)/dble(trotter) &
&             - spring(iatom,ispring)*xcart(ii,iatom,iimage)
       end do
     end do
   end do

 end if

end subroutine pimd_forces
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_langevin_forces
!! NAME
!!  pimd_langevin_forces
!!
!! FUNCTION
!!  Compute Langevin contribution to PIMD forces
!!
!! INPUTS
!!  alea(3,natom,trotter)=set of random numbers
!!  forces(3,natom,trotter)=forces without Langevin contribution
!!  friction=friction factor
!!  langev(natom)=Langevin factors
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  trotter=Trotter number
!!  vel(3,natom,trotter)=velocities of atoms in each cell
!!
!! OUTPUT
!!  forces_langevin(3,natom,trotter)=forces including Langevin contribution
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_langevin_forces(alea,forces,forces_langevin,friction,&
&                               langev,mass,natom,trotter,vel)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_langevin_forces'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,trotter
 real(dp),intent(in) :: friction
!arrays
 real(dp),intent(in) :: alea(3,natom,trotter),forces(3,natom,trotter)
 real(dp),intent(in) :: langev(natom),vel(3,natom,trotter)
 real(dp),intent(in) :: mass(:,:)
 real(dp),intent(out) :: forces_langevin(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,natom_mass,nmass
 character(len=500) :: msg
!arrays

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

 do iimage=1,trotter
   imass=min(nmass,iimage)
   do iatom=1,natom
     do ii=1,3
       forces_langevin(ii,iatom,iimage)=forces(ii,iatom,iimage) &
&                    + langev(iatom)*alea(ii,iatom,iimage) &
&                    - friction*mass(iatom,imass)*vel(ii,iatom,iimage)
     end do
   end do
 end do

end subroutine pimd_langevin_forces
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_noseehoover_forces
!! NAME
!!  pimd_nosehoover_forces
!!
!! FUNCTION
!!  Compute Nose-Hoover contribution to PIMD forces
!!  by adding friction force of thermostat number one
!!
!! INPUTS
!!  dzeta(3,natom,trotter,nnos)=variables of thermostats, in (atomic time unit)^(-1)
!!                              used only when a coordinate transformation is applied (transfom/=0)
!!  forces(3,natom,trotter)=forces without Nose-Hoover contribution
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  nnos=number of thermostats
!!  transform=coordinate transformation:
!!            0: no tranformation
!!            1: normal mode transformation
!!            2: staging transformation
!!  trotter=Trotter number
!!  vel(3,natom,trotter)=velocities of atoms in each cell
!!  zeroforce=if 1, forces are constrained to be zero
!!  zeta(nnos)=variables of thermostats, in (atomic time unit)^(-1)
!!             used only when no coordinate transformation is applied (transfom==0)
!!
!! OUTPUT
!!  forces_nosehoover(3,natom,trotter)=forces including thermostat contribution
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_nosehoover_forces(dzeta,forces,forces_nosehoover,mass,natom,&
&                                 nnos,transform,trotter,vel,zeroforce,zeta)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_nosehoover_forces'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nnos,transform,trotter,zeroforce
!arrays
 real(dp),intent(in) :: dzeta(3,natom,trotter,nnos),forces(3,natom,trotter)
 real(dp),intent(in) :: vel(3,natom,trotter),zeta(nnos)
 real(dp),intent(in) :: mass(:,:)
 real(dp),intent(out) :: forces_nosehoover(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,natom_mass,nmass
 character(len=500) :: msg
!arrays
 real(dp) :: total(3)

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

 if (transform==0) then
   do iimage=1,trotter
     imass=min(nmass,iimage)
     do iatom=1,natom
       do ii=1,3
         forces_nosehoover(ii,iatom,iimage)=forces(ii,iatom,iimage) &
&            - mass(iatom,imass)*zeta(1)*vel(ii,iatom,iimage)
       end do
     end do
   end do

 else if (transform==1) then

! TO BE IMPLEMENTED

 else if (transform==2) then
   do iimage=1,trotter
     imass=min(nmass,iimage)
     do iatom=1,natom
       do ii=1,3
         forces_nosehoover(ii,iatom,iimage)=forces(ii,iatom,iimage) &
&            - mass(iatom,imass)*dzeta(ii,iatom,iimage,1)*vel(ii,iatom,iimage)
       end do
     end do
   end do

 end if

!Make sure that the sum of random forces is zero
 if(zeroforce==1)then
   total=zero
   do iimage=1,trotter
     do iatom=1,natom
       do ii=1,3
         total(ii)=total(ii)+forces_nosehoover(ii,iatom,iimage)
       end do
     end do
   end do
   total=total/dble(natom*trotter)
   do iimage=1,trotter
     do iatom=1,natom
       do ii=1,3
         forces_nosehoover(ii,iatom,iimage)=forces_nosehoover(ii,iatom,iimage)-total(ii)
       end do
     end do
   end do
 end if

end subroutine pimd_nosehoover_forces
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_stresses
!! NAME
!!  pimd_stresses
!!
!! FUNCTION
!!  In the case od PIMD, compute the stress tensor from virial theorem
!!
!! INPUTS
!!  forces_pimd(3,natom,trotter)=PIMD forces on atoms in each cell
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  trotter=Trotter number
!!  vel(3,natom,trotter)=velocities of atoms in each cell
!!  volume=volume of each cell (common to all cells)
!!  xcart(3,natom,trotter)=cartesian coordinates of atoms in each cell at t
!!
!! OUTPUT
!!  stress_pimd(3,3)=stress tensor for PIMD
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_npt
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_stresses(forces_pimd,mass,natom,stress_pimd,trotter,vel,volume,xcart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_stresses'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,trotter
 real(dp),intent(in) :: volume
!arrays
 real(dp),intent(in) :: forces_pimd(3,natom,trotter),vel(3,natom,trotter),xcart(3,natom,trotter)
 real(dp),intent(in) :: mass(:,:)
 real(dp),intent(out) :: stress_pimd(3,3)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,jj,natom_mass,nmass
 character(len=500) :: msg
!arrays

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

 stress_pimd=zero

!1-Kinetic part
 do iimage=1,trotter
  imass=min(nmass,iimage)
  do iatom=1,natom
    do ii=1,3
      do jj=1,3
        stress_pimd(ii,jj)=stress_pimd(ii,jj)+ &
&          mass(iatom,imass)*vel(ii,iatom,iimage)*vel(jj,iatom,iimage)
      end do
    end do
  end do
 end do

!2-Potential part
 do iimage=1,trotter
  do iatom=1,natom
    do ii=1,3
      do jj=1,3
        stress_pimd(ii,jj)=stress_pimd(ii,jj)+ &
&          xcart(ii,iatom,iimage)*forces_pimd(jj,iatom,iimage)
      end do
    end do
  end do
 end do
 stress_pimd=stress_pimd/(three*volume)

end subroutine pimd_stresses
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_diff_stress
!! NAME
!!  pimd_diff_stress
!!
!! FUNCTION
!!  Compute the difference between the stress tensor and the stress target
!!
!! INPUTS
!!  stress_pimd(3,3)=stress tensor for PIMD
!!  stress_target(6)=stress target
!!
!! OUTPUT
!!  pimd_diff_stress(3,3)=difference between stresses and stress target
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function pimd_diff_stress(stress_pimd,stress_target)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_diff_stress'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 real(dp),intent(in) :: stress_pimd(3,3),stress_target(6)
 real(dp) :: pimd_diff_stress(3,3)
!Local variables-------------------------------
!scalars
!arrays

!************************************************************************

 pimd_diff_stress(1,1)=stress_pimd(1,1)-stress_target(1)
 pimd_diff_stress(2,2)=stress_pimd(2,2)-stress_target(2)
 pimd_diff_stress(3,3)=stress_pimd(3,3)-stress_target(3)
 pimd_diff_stress(2,3)=stress_pimd(2,3)-stress_target(4)
 pimd_diff_stress(3,2)=stress_pimd(3,2)-stress_target(4)
 pimd_diff_stress(1,3)=stress_pimd(1,3)-stress_target(5)
 pimd_diff_stress(3,1)=stress_pimd(3,1)-stress_target(5)
 pimd_diff_stress(1,2)=stress_pimd(1,2)-stress_target(6)
 pimd_diff_stress(2,1)=stress_pimd(2,1)-stress_target(6)

end function pimd_diff_stress
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_predict_taylor
!! NAME
!!  pimd_predict_taylor
!!
!! FUNCTION
!!  Predict new atomic positions using a Taylor algorithm (first time step) - for PIMD
!!
!! INPUTS
!!  dtion=time step
!!  forces(3,natom,trotter)=PIMD forces on atoms in each cell
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  trotter=Trotter number
!!  vel(3,natom,trotter)=velocities of atoms in each cell
!!  xcart(3,natom,trotter)=cartesian coordinates of atoms in each cell at t
!!
!! OUTPUT
!!  xcart_next(3,natom,trotter)=cartesian coordinates of atoms in each cell at t+dt
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_predict_taylor(dtion,forces,mass,natom,trotter,vel,xcart,xcart_next)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_predict_taylor'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,trotter
 real(dp),intent(in) :: dtion
!arrays
 real(dp),intent(in) :: forces(3,natom,trotter),vel(3,natom,trotter),xcart(3,natom,trotter)
 real(dp),intent(in) :: mass(:,:)
 real(dp),intent(out) :: xcart_next(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,natom_mass,nmass
 character(len=500) :: msg
!arrays

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

 do iimage=1,trotter
   imass=min(nmass,iimage)
   do iatom=1,natom
     do ii=1,3
       xcart_next(ii,iatom,iimage)=xcart(ii,iatom,iimage) &
&             + half*dtion*dtion*forces(ii,iatom,iimage)/mass(iatom,imass) &
&             + dtion*vel(ii,iatom,iimage)
     end do
   end do
 end do

end subroutine pimd_predict_taylor
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_predict_verlet
!! NAME
!!  pimd_predict_verlet
!!
!! FUNCTION
!!  Predict new atomic positions using a Verlet algorithm - for PIMD
!!
!! INPUTS
!!  dtion=time step
!!  forces(3,natom,trotter)=PIMD forces on atoms in each cell
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  trotter=Trotter number
!!  xcart(3,natom,trotter)=cartesian coordinates of atoms in each cell at t
!!  xcart_prev(3,natom,trotter)=cartesian coordinates of atoms in each cell at t-dt
!!
!! OUTPUT
!!  xcart_next(3,natom,trotter)=cartesian coordinates of atoms in each cell at t+dt
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_langevin_nvt,pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_predict_verlet(dtion,forces,mass,natom,trotter,xcart,xcart_next,xcart_prev)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_predict_verlet'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,trotter
 real(dp),intent(in) :: dtion
!arrays
 real(dp),intent(in) :: forces(3,natom,trotter)
 real(dp),intent(in) :: xcart(3,natom,trotter),xcart_prev(3,natom,trotter)
 real(dp),intent(in) :: mass(:,:)
 real(dp),intent(out) :: xcart_next(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,natom_mass,nmass
 character(len=500) :: msg
!arrays

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if

 do iimage=1,trotter
   imass=min(nmass,iimage)
   do iatom=1,natom
     do ii=1,3
       xcart_next(ii,iatom,iimage)= &
&         two*xcart(ii,iatom,iimage) &
&       - xcart_prev(ii,iatom,iimage) &
&       + dtion*dtion*forces(ii,iatom,iimage)/mass(iatom,imass)
     end do
   end do
 end do

end subroutine pimd_predict_verlet
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_nosehoover_propagate
!! NAME
!!  pimd_nosehoover_propagate
!!
!! FUNCTION
!!  Propagate thermostat variables (Nose-Hoover algorithm) - for PIMD
!!
!! INPUTS
!!  dtion=time step
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  natom=number of atoms
!!  nnos=number of thermostats
!!  qmass(nnos)=masses of thermostats
!!  temperature=temperature
!!  transform=coordinate transformation:
!!            0: no tranformation
!!            1: normal mode transformation
!!            2: staging transformation
!!  trotter=Trotter number
!!  vel(3,natom,trotter)=velocities of atoms in each cell
!!
!! OUTPUT
!!  dzeta(3,natom,trotter,nnos)=variables of thermostats, in (atomic time unit)^(-1)
!!                              used only when a coordinate transformation is applied (transfom/=0)
!!  zeta(nnos)=variables of thermostats, in (atomic time unit)^(-1)
!!             used only when no coordinate transformation is applied (transfom==0)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_nosehoover_propagate(dtion,dzeta,mass,natom,nnos,qmass,temperature,&
&                                    transform,trotter,vel,zeta)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_nosehoover_propagate'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nnos,transform,trotter
 real(dp),intent(in) :: dtion,temperature
!arrays
 real(dp),intent(in) :: qmass(nnos),vel(3,natom,trotter)
 real(dp),intent(in) :: mass(:,:)
 real(dp),intent(inout) :: dzeta(3,natom,trotter,nnos),zeta(nnos)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,imass,inos,natom_mass,nmass
 real(dp) :: gg,ec,kt
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: thermforces(:,:,:,:)

!************************************************************************

 natom_mass=size(mass,1);nmass=size(mass,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if
 if (nnos<3) then
   msg='Not available for nnos<3 !'
   MSG_BUG(msg)
 end if 

 kt=temperature*kb_HaK

!=== No transformation ===================================================
 if (transform==0) then
   gg=dble(3*natom*trotter)

   ec=zero
   do iimage=1,trotter
     imass=min(nmass,iimage)
     do iatom=1,natom
       do ii=1,3
         ec=ec+mass(iatom,imass)*vel(ii,iatom,iimage)**2
       end do
     end do
   end do

   do inos=1,nnos
     do iimage=1,trotter
       do iatom=1,natom
         do ii=1,3
           if (inos==1) then
             zeta(1)=zeta(1)+dtion*((ec-gg*kt)/qmass(1)-zeta(1)*zeta(2))
           else if(inos==nnos)then
             zeta(nnos)=zeta(nnos)+dtion*(qmass(nnos-1)*zeta(nnos-1)*zeta(nnos-1)-kt)/qmass(nnos)
           else
             zeta(inos)=zeta(inos)+dtion*((qmass(inos-1)*zeta(inos-1)*zeta(inos-1)-kt)/qmass(inos) &
&                                        -zeta(inos)*zeta(inos+1))
           end if
         end do
       end do
     end do
   end do

!=== Normal mode transformation ==========================================
 else if (transform==1) then

!  TO BE IMPLEMENTED

!=== Staging transformation ==============================================
 else if (transform==2) then

   ABI_ALLOCATE(thermforces,(3,natom,trotter,nnos))
   thermforces=zero
   do inos=1,nnos
     do iimage=1,trotter
       imass=min(nmass,iimage)
       do iatom=1,natom
         do ii=1,3
           if (inos==1) then
             thermforces(ii,iatom,iimage,inos)=mass(iatom,imass)*vel(ii,iatom,iimage)**2-kt
           else if (inos==nnos-1) then
             thermforces(ii,iatom,iimage,inos)=qmass(nnos-1) &
&                 *(dzeta(ii,iatom,iimage,nnos-2)**2+dzeta(ii,iatom,iimage,nnos)**2)-two*kt
           else
             thermforces(ii,iatom,iimage,inos)=qmass(inos)*dzeta(ii,iatom,iimage,inos-1)**2-kt
           end if
         end do
       end do
     end do
   end do

   do inos=1,nnos
     do iimage=1,trotter
       do iatom=1,natom
         do ii=1,3
           if (inos==nnos) then
             dzeta(ii,iatom,iimage,nnos)=dzeta(ii,iatom,iimage,nnos) &
&                    + dtion*((thermforces(ii,iatom,iimage,nnos)/qmass(nnos)) &
&                            -dzeta(ii,iatom,iimage,nnos)*dzeta(ii,iatom,iimage,nnos-1))
           else if (inos==nnos-1) then
             dzeta(ii,iatom,iimage,nnos-1)=dzeta(ii,iatom,iimage,nnos-1) &
&                    + dtion*((thermforces(ii,iatom,iimage,nnos-1)/qmass(nnos-1)) &
&                            -dzeta(ii,iatom,iimage,nnos)*dzeta(ii,iatom,iimage,nnos-1))
           else
             dzeta(ii,iatom,iimage,inos)=dzeta(ii,iatom,iimage,inos) &
&                    + dtion*((thermforces(ii,iatom,iimage,inos)/qmass(inos)) &
&                            -dzeta(ii,iatom,iimage,inos)*dzeta(ii,iatom,iimage,inos+1))
           end if
         end do
       end do
     end do
   end do

   ABI_DEALLOCATE(thermforces)

 end if ! transform

end subroutine pimd_nosehoover_propagate
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_coord_transform
!! NAME
!!  pimd_coord_transform
!!
!! FUNCTION
!!  Apply a coordinate transformation on a given vector field
!!  (defined for each atom in in cell) - for PIMD
!!  Possible choices for the transformation:
!!    0: no transformation
!!    1: normal mode tranformation
!!    2: staging transformation
!!
!! INPUTS
!!  ioption=option given the direction of the transformation
!!     +1: from primitive coordinates to transformed coordinates
!!     -1: from transformed coordinates to primitive coordinates
!!  natom=number of atoms
!!  transform=coordinate transformation:
!!            0: no tranformation
!!            1: normal mode transformation
!!            2: staging transformation
!!  trotter=Trotter number
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  array(3,natom,trotter)=array to be transformed
!!
!! PARENTS
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_coord_transform(array,ioption,natom,transform,trotter)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_coord_transform'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ioption,natom,transform,trotter
!arrays
 real(dp),intent(inout) :: array(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage,iimagem,iimagep,iimage2
!arrays
 real(dp),allocatable :: array_temp(:,:,:)

!************************************************************************

!=== No transformation ===================================================
 if (transform==0) then

   return

!=== Normal mode transformation ==========================================
 else if (transform==1) then

!  TO BE IMPLEMENTED
!  ---------------- From primitive to transformed coordinates ------------
   if (ioption==+1) then
!  ---------------- From transformed to primitive coordinates ------------
   else if (ioption==-1) then
   end if ! ioption

!=== Staging transformation ==============================================
 else if (transform==2) then

!  ---------------- From primitive to transformed coordinates ------------
   if (ioption==+1) then

     ABI_ALLOCATE(array_temp,(3,natom,trotter))
     array_temp=zero
     do iimage=1,trotter
       iimagep=iimage+1;if(iimage==trotter) iimagep=1
       iimagem=iimage-1;if(iimage==1) iimagem=trotter
       do iatom=1,natom
         do ii=1,3
           array_temp(ii,iatom,iimage)=(dble(iimagem)*array(ii,iatom,iimagep)+array(ii,iatom,1))/dble(iimage)
         end do
       end do
     end do
     if (trotter>1) then
       do iimage=2,trotter
         do iatom=1,natom
           do ii=1,3
             array(ii,iatom,iimage)=array(ii,iatom,iimage)-array_temp(ii,iatom,iimage)
           end do
         end do
       end do
     end if
     ABI_DEALLOCATE(array_temp)

!  ---------------- From transformed to primitive coordinates ------------
   else if (ioption==-1) then

     ABI_ALLOCATE(array_temp,(3,natom,trotter))
     array_temp=zero
     do iimage=1,trotter
       do iatom=1,natom
         do ii=1,3
           array_temp(ii,iatom,iimage)=array(ii,iatom,iimage)
         end do
       end do
     end do
     if (trotter>1) then
       do iimage=2,trotter
         do iatom=1,natom
           do ii=1,3
             do iimage2=iimage,trotter
               array_temp(ii,iatom,iimage)=array_temp(ii,iatom,iimage) &
&                     +array(ii,iatom,iimage2)*(dble(iimage-1))/(dble(iimage2-1))
             end do
           end do
         end do
       end do
     end if
     array(:,:,:)=array_temp(:,:,:)
     ABI_DEALLOCATE(array_temp)

   end if ! ioption

 end if ! transform

end subroutine pimd_coord_transform
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_force_transform
!! NAME
!!  pimd_force_transform
!!
!! FUNCTION
!!  Apply a coordinate transformation on forces (defined for each atom in in cell) - for PIMD
!!  Possible choices for the transformation:
!!    0: no transformation
!!    1: normal mode tranformation
!!    2: staging transformation
!!
!! INPUTS
!!  ioption=option given the direction of the transformation
!!     +1: from primitive coordinates to transformed coordinates
!!     -1: from transformed coordinates to primitive coordinates
!!  natom=number of atoms
!!  transform=coordinate transformation:
!!            0: no tranformation
!!            1: normal mode transformation
!!            2: staging transformation
!!  trotter=Trotter number
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  forces(3,natom,trotter)=array containing forces
!!
!! NOTES
!!  Back transformation (ioption=-1) not implemented !
!!
!! PARENTS
!!      pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_force_transform(forces,ioption,natom,transform,trotter)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_force_transform'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ioption,natom,transform,trotter
!arrays
 real(dp),intent(inout) :: forces(3,natom,trotter)
!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iimage
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: forces_temp(:,:,:)

!************************************************************************

 if (ioption==-1) then
   msg='Back transformation not implemented !'
   MSG_BUG(msg)
 end if

!=== No transformation ===================================================
 if (transform==0) then

   return

!=== Normal mode transformation ==========================================
 else if (transform==1) then

!  TO BE IMPLEMENTED

!=== Staging transformation ==============================================
 else if (transform==2) then

   ABI_ALLOCATE(forces_temp,(3,natom,trotter))
   forces_temp=zero
   do iimage=1,trotter
     do iatom=1,natom
       do ii=1,3
         forces_temp(ii,iatom,1)=forces_temp(ii,iatom,1)+forces(ii,iatom,iimage)
       end do
     end do
   end do
   if (trotter>1) then
     do iimage=2,trotter
       do iatom=1,natom
         do ii=1,3
           forces_temp(ii,iatom,iimage)=forces(ii,iatom,iimage) &
&               +forces_temp(ii,iatom,iimage-1)*(dble(iimage-2)/dble(iimage-1))
         end do
       end do
     end do
   end if
   forces=forces_temp
   ABI_DEALLOCATE(forces_temp)

 end if ! transform

end subroutine pimd_force_transform
!!***

!----------------------------------------------------------------------

!!****f* m_pimd/pimd_mass_spring
!! NAME
!!  pimd_mass_spring
!!
!! FUNCTION
!!  Compute masses and spring constants for PIMD. Eventually apply a coordinate transformation.
!!  Possible choices for the transformation:
!!    0: no transformation
!!    1: normal mode tranformation
!!    2: staging transformation
!!
!! INPUTS
!!  inertmass(natom)=inertial masses of atoms
!!  kt=kT constant
!!  natom=number of atoms
!!  quantummass(natom)=quantum masses of atoms
!!  transform=coordinate transformation:
!!            0: no tranformation
!!            1: normal mode transformation
!!            2: staging transformation
!!  trotter=Trotter number
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mass(natom,mass_dim)=masses of atoms (mass_dim=1 or trotter)
!!  spring(natom,mass_dim)=spring constants of atoms (mass_dim=1 or trotter)
!!
!! NOTES
!!  Back transformation (ioption=-1) not implemented !
!!
!! PARENTS
!!      pimd_langevin_npt,pimd_langevin_nvt,pimd_nosehoover_nvt
!!
!! CHILDREN
!!
!! SOURCE

subroutine pimd_mass_spring(inertmass,kt,mass,natom,quantummass,spring,transform,trotter)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_mass_spring'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,transform,trotter
 real(dp),intent(in) :: kt
!arrays
 real(dp),intent(in) :: inertmass(natom),quantummass(natom)
 real(dp),intent(out) :: mass(:,:),spring(:,:)
!Local variables-------------------------------
!scalars
 integer :: iimage,natom_mass,natom_spring,nmass,nspring
 character(len=500) :: msg
!arrays

!************************************************************************

 natom_mass  =size(mass  ,1);nmass  =size(mass  ,2)
 natom_spring=size(spring,1);nspring=size(spring,2)
 if (natom/=natom_mass.or.(nmass/=1.and.nmass/=trotter)) then
   msg='Wrong dimensions for array mass !'
   MSG_BUG(msg)
 end if
 if (natom/=natom_spring.or.(nspring/=1.and.nspring/=trotter)) then
   msg='Wrong dimensions for array spring !'
   MSG_BUG(msg)
 end if

!=== No transformation ===================================================
 if (transform==0) then

   mass(1:natom,1)=inertmass(1:natom)
   if (nmass>1) then
     do iimage=2,trotter
       mass(1:natom,iimage)=inertmass(1:natom)
     end do
   end if

   spring(1:natom,1)=quantummass(1:natom)*dble(trotter)*kt*kt
   if (nspring>1) then
     do iimage=2,trotter
       spring(1:natom,iimage)=quantummass(1:natom)*dble(trotter)*kt*kt
     end do
   end if

!=== Normal mode transformation ==========================================
 else if (transform==1) then

!    TO BE IMPLEMENTED

!=== Staging transformation ==============================================
 else if (transform==2) then

   mass(1:natom,1)=inertmass(1:natom)
   if (nmass>1) then
     do iimage=2,trotter
       mass(1:natom,iimage)=inertmass(1:natom)*dble(iimage)/dble(iimage-1)
     end do
   end if

   spring(1:natom,1)=mass(1:natom,1)*dble(trotter)*kt*kt
   if (nspring>1) then
     do iimage=2,trotter
       spring(1:natom,iimage)=mass(1:natom,iimage)*dble(trotter)*kt*kt
     end do
   end if
 end if

end subroutine pimd_mass_spring
!!***

END MODULE m_pimd
!!***
