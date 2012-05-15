!{\src2tex{textfont=tt}}
!!****f* ABINIT/pimd_nosehoover_nvt
!! NAME
!! pimd_nosehoover_nvt
!!
!! FUNCTION
!! Predicts new positions in Path Integral Molecular Dynamics using Nose-Hoover in the NVT ensemble.
!! Given the positions at time t and t-dtion, an estimation of the velocities at time t,
!! the forces and an estimation of the stress at time t, and an estimation of the cell at time t,
!! computes in the Path Integral Molecular Dynamics framework the new positions at time t+dtion,
!! computes self-consistently the velocities, the stress and the cell at time t and produces
!! an estimation of the velocities, stress and new cell at time t+dtion
!! No change of acell and rprim at present.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (GG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  etotal(trotter)=electronic total energy for all images
!!  itimimage=number of the current time for image propagation (itimimage+1 is to be predicted here)
!!  natom=dimension of vel_timimage and xred_timimage
!!  pimd_param=datastructure that contains all the parameters necessary to Path-Integral MD
!!  prtvolimg=printing volume
!!  rprimd(3,3)=dimensionless unit cell vectors (common to all images)
!!  trotter=Trotter number (total number of images)
!!  volume=voume of unit cell (common to all images)
!!  xred(3,natom,trotter)=reduced coordinates of atoms for all images at time t (present time step)
!!  xred_prev(3,natom,trotter)=reduced coordinates of atoms for all images at time t-dt (previous time step)
!!
!! OUTPUT
!!  xred_next(3,natom,trotter)=reduced coordinates of atoms for all images at time t+dt (next time step)
!!
!! SIDE EFFECTS
!!  forces(3,natom,trotter)=forces over atoms for all images
!!    at input,  electronic forces
!!    at output, electronic forces + quantum spring contribution
!!  vel(3,natom,trotter)=velocies of atoms for all images
!!    at input,  values at time t
!!    at output, values at time t+dt
!!
!! NOTES
!!  Thermization by Nose-Hoover chains according to
!!  Martyna, Klein, Tuckerman, J. Chem. Phys. 97, 2635 (1992)
!!  Tuckerman, Marx, Klein, Parrinello, J. Chem. Phys. 104, 5579 (1996)

!!
!! PARENTS
!!      predict_pimd
!!
!! CHILDREN
!!      pimd_coord_transform,pimd_energies,pimd_force_transform,pimd_forces
!!      pimd_initvel,pimd_mass_spring,pimd_nosehoover_forces
!!      pimd_nosehoover_propagate,pimd_predict_taylor,pimd_predict_verlet
!!      pimd_print,pimd_stresses,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pimd_nosehoover_nvt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&
&                              rprimd,trotter,vel,volume,xred,xred_next,xred_prev)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_pimd

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_nosehoover_nvt'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage,natom,prtvolimg,trotter
 real(dp),intent(in) :: volume
 type(pimd_type),intent(in) :: pimd_param
!arrays
 real(dp),intent(in) :: etotal(trotter),rprimd(3,3)
 real(dp),intent(in),target :: xred(3,natom,trotter),xred_prev(3,natom,trotter)
 real(dp),intent(out) :: xred_next(3,natom,trotter)
 real(dp),intent(inout) :: forces(3,natom,trotter),vel(3,natom,trotter)

!Local variables-------------------------------
!Options
 real(dp),parameter :: tolerance=tol7 ! SCF tolerance
!scalars
 integer :: idum=-5
 integer :: iimage,irestart,ndof,nnos,pitransform,prtstress
 real(dp) :: dtion,eharm,eharm2,epot,initemp,kt
 real(dp) :: temperature1,temperature2,temp2_prev,thermtemp,tol
 character(len=500) :: msg
!arrays
 real(dp) :: stress_pimd(3,3)
 real(dp),allocatable :: dzeta(:,:,:,:),forces_orig(:,:,:),forces_pimd(:,:,:)
 real(dp),allocatable :: inertmass(:),masseff(:,:),qmass(:),quantummass(:),springeff(:,:)
 real(dp),allocatable :: xcart(:,:,:),xcart_next(:,:,:),xcart_prev(:,:,:),zeta(:)
 real(dp),pointer :: xred_tmp(:,:)

! *************************************************************************

!############# Initializations ###########################

!Allocation of local arrays
 ABI_ALLOCATE(xcart,(3,natom,trotter))
 ABI_ALLOCATE(xcart_prev,(3,natom,trotter))
 ABI_ALLOCATE(xcart_next,(3,natom,trotter))
 ABI_ALLOCATE(forces_orig,(3,natom,trotter))
 ABI_ALLOCATE(forces_pimd,(3,natom,trotter))
 ABI_ALLOCATE(inertmass,(natom))
 ABI_ALLOCATE(quantummass,(natom))

!Fill in the local variables
 ndof=3*natom*trotter
 quantummass(1:natom)=pimd_param%amu   (pimd_param%typat(1:natom))*amu_emass
 inertmass  (1:natom)=pimd_param%pimass(pimd_param%typat(1:natom))*amu_emass
 initemp=pimd_param%mdtemp(1);thermtemp=pimd_param%mdtemp(2)
 dtion=pimd_param%dtion;pitransform=pimd_param%pitransform
 kt=thermtemp*kb_HaK
 forces_orig=forces

!Allocation/initialization of local variables used for Nose-Hoover chains
!Associated variables:
!nnos = number of thermostats
!dzeta(3,natom,trotter,nnos) = variables of thermostats, in (atomic time unit)^(-1)
!qmass(nnos) = masses of thermostats
!specific to PIMD: pitransform = coordinate transformation (0:no; 1:normal mode; 2:staging)
 nnos=pimd_param%nnos
 ABI_ALLOCATE(qmass,(nnos))
 ABI_ALLOCATE(zeta,(nnos))
 ABI_ALLOCATE(dzeta,(3,natom,trotter,nnos))
 qmass(1:nnos)=pimd_param%qmass(1:nnos)
 zeta=zero;dzeta=zero

!Effective masses and spring constants (according to pitransform)
 ABI_ALLOCATE(masseff,(natom,trotter))
 ABI_ALLOCATE(springeff,(natom,trotter))
 call pimd_mass_spring(inertmass,kt,masseff,natom,quantummass,springeff,pitransform,trotter)

!Recommended value of Nose mass
 write(msg,'(2a,f9.5,3a)') ch10,&
& ' Recommended value of Nose mass is',one/(dble(trotter)*kt),'(atomic units)',ch10,&
& '(see Tuckerman et al, J. Chem. Phys. 104, 5579 (1996))'
 call wrtout(std_out,msg,'COLL')

!Compute cartesian coordinates
 do iimage=1,trotter
   xred_tmp => xred(:,:,iimage)      ! Needed because xred is inten(inout) in xredxcart
   call xredxcart(natom,1,rprimd,xcart     (:,:,iimage),xred_tmp)
   xred_tmp => xred_prev(:,:,iimage) ! Needed because xred_tmp is inten(inout) in xredxcart
   call xredxcart(natom,1,rprimd,xcart_prev(:,:,iimage),xred_tmp)
 end do

!Transform the coordinates and forces (according to pitransform)
 call pimd_coord_transform(xcart     ,+1,natom,pitransform,trotter)
 call pimd_coord_transform(xcart_prev,+1,natom,pitransform,trotter)
 call pimd_force_transform(forces    ,+1,natom,pitransform,trotter)

!Determine if it is a restart or not
!If this is a calculation from scratch,generate random distribution of velocities
 irestart=1;if (itimimage==1) irestart=pimd_is_restart(masseff,vel)
 if (irestart==0) then
   call pimd_initvel(idum,masseff,natom,initemp,trotter,vel)
 end if

!Compute temperature at t
 temperature1=pimd_temperature(masseff,vel)

!Transform the velocities (according to pitransform)
 call pimd_coord_transform(vel,+1,natom,pitransform,trotter)


!################## Images evolution #####################

!Compute PIMD and 1st thermostat contributions to forces
 call pimd_forces(forces,natom,springeff,pitransform,trotter,xcart)
 call pimd_nosehoover_forces(dzeta,forces,forces_pimd,masseff,natom,&
& nnos,pitransform,trotter,vel,1,zeta)

!Compute atomic positions at t+dt
 if (itimimage<=1) then

!  === 1st time step: single Taylor algorithm
!  Predict positions
   call pimd_predict_taylor(dtion,forces_pimd,masseff,natom,trotter,&
&   vel,xcart,xcart_next)
!  Estimate the velocities at t+dt/2
   vel=(xcart_next-xcart)/dtion
!  Propagate the thermostat variables
   call pimd_nosehoover_propagate(dtion,dzeta,masseff,natom,nnos,qmass,&
&   thermtemp,pitransform,trotter,vel,zeta)
 else

!  === Other time steps: Verlet algorithm + SC cycle
!  Predict positions
   call pimd_predict_verlet(dtion,forces_pimd,masseff,natom,trotter,&
&   xcart,xcart_next,xcart_prev)
!  Propagate the thermostat variables
   call pimd_nosehoover_propagate(dtion,dzeta,masseff,natom,nnos,qmass,&
&   thermtemp,pitransform,trotter,vel,zeta)
!  Self-consistent loop
   temperature2=pimd_temperature(masseff,vel)
   temp2_prev=temperature2; tol=one
   do while (tol>tolerance)
!    Recompute a (better) estimation of the velocity at time step t
     vel = (xcart_next - xcart_prev) / (two*dtion)
     temperature2=pimd_temperature(masseff,vel)
!    Reestimate the force
     call pimd_nosehoover_forces(dzeta,forces,forces_pimd,masseff,natom,&
&     nnos,pitransform,trotter,vel,1,zeta)
!    Compute new positions
     call pimd_predict_verlet(dtion,forces_pimd,masseff,natom,trotter,&
&     xcart,xcart_next,xcart_prev)
!    Propagate the thermostat variables
     call pimd_nosehoover_propagate(dtion,dzeta,masseff,natom,nnos,qmass,&
&     thermtemp,pitransform,trotter,vel,zeta)
!    Compute variation of temperature (to check convergence of SC loop)
     tol=dabs(temperature2-temp2_prev)/dabs(temp2_prev)
     temp2_prev=temperature2
   end do ! End self-consistent loop

 end if ! itimimage==1

!Come back to primitive coordinates and velocities
 call pimd_coord_transform(xcart     ,-1,natom,pitransform,trotter)
 call pimd_coord_transform(xcart_prev,-1,natom,pitransform,trotter)
 call pimd_coord_transform(vel       ,-1,natom,pitransform,trotter)

!Compute new temperature
 temperature2=pimd_temperature(masseff,vel)

!Compute contributions to energy
 call pimd_energies(eharm,eharm2,epot,etotal,forces_orig,natom,springeff,trotter,xcart)

!Compute stress tensor at t from virial theorem
 call pimd_stresses(forces_pimd,masseff,natom,stress_pimd,trotter,vel,volume,xcart)

!If possible, estimate the velocities at t+dt
 if (itimimage>1) then
   vel = (three*xcart_next - four*xcart + xcart_prev)/(two * dtion)
 end if


!############# Final operations ############################

!Print messages
 prtstress=1;if (prtvolimg>=2) prtstress=0
 call pimd_print(eharm,eharm2,epot,forces_pimd,irestart,itimimage,natom,prtstress,prtvolimg,&
& stress_pimd,temperature1,temperature2,trotter,vel,xcart,xred)

!Come back to reduced coordinates
 do iimage=1,trotter
   call xredxcart(natom,-1,rprimd,xcart_next(:,:,iimage),xred_next(:,:,iimage))
 end do

!Free memory
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xcart_prev)
 ABI_DEALLOCATE(xcart_next)
 ABI_DEALLOCATE(forces_orig)
 ABI_DEALLOCATE(forces_pimd)
 ABI_DEALLOCATE(inertmass)
 ABI_DEALLOCATE(quantummass)
 ABI_DEALLOCATE(masseff)
 ABI_DEALLOCATE(springeff)
 ABI_DEALLOCATE(qmass)
 ABI_DEALLOCATE(dzeta)
 ABI_DEALLOCATE(zeta)

end subroutine pimd_nosehoover_nvt
!!***


