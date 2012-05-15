!{\src2tex{textfont=tt}}
!!****f* ABINIT/pimd_nosehoover_npt
!! NAME
!! pimd_nosehoover_npt
!!
!! FUNCTION
!! Predicts new positions in Path Integral Molecular Dynamics using Nose-Hoover in the NPT ensemble.
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
!!  rprimd(3,3)=dimensionless unit cell vectors (common to all images) at time t (present time step)
!!  rprimd_prev(3,3)=dimensionless unit cell vectors (common to all images) at time t-dt (previous time step)
!!  trotter=Trotter number (total number of images)
!!  volume=voume of unit cell (common to all images)
!!  xred(3,natom,trotter)=reduced coordinates of atoms for all images at time t (present time step)
!!  xred_prev(3,natom,trotter)=reduced coordinates of atoms for all images at time t-dt (previous time step)
!!
!! OUTPUT
!!  rprimd_next(3,3)=dimensionless unit cell vectors (common to all images) at time t+dt (next time step)
!!  xred_next(3,natom,trotter)=reduced coordinates of atoms for all images at time t+dt (next time step)
!!
!! SIDE EFFECTS
!!  forces(3,natom,trotter)=forces over atoms for all images
!!    at input,  electronic forces
!!    at output, electronic forces + Langevin contribution
!!  vel(3,natom,trotter)=velocies of atoms for all images
!!    at input,  values at time t
!!    at output, values at time t+dt
!!
!! PARENTS
!!      predict_pimd
!!
!! CHILDREN
!!      pimd_energies,pimd_initvel,pimd_print,pimd_stresses,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pimd_nosehoover_npt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&
&          rprimd,rprimd_next,rprimd_prev,trotter,vel,volume,xred,xred_next,xred_prev)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_pimd

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_nosehoover_npt'
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimimage,natom,prtvolimg,trotter
 real(dp),intent(in) :: volume
 type(pimd_type),intent(in) :: pimd_param
!arrays
 real(dp),intent(in) :: etotal(trotter),rprimd(3,3),rprimd_prev(3,3)
 real(dp),intent(in),target :: xred(3,natom,trotter),xred_prev(3,natom,trotter)
 real(dp),intent(out) :: rprimd_next(3,3),xred_next(3,natom,trotter)
 real(dp),intent(inout) :: forces(3,natom,trotter),vel(3,natom,trotter)

!Local variables-------------------------------
!Options
 real(dp),parameter :: tolerance=tol7 ! SCF tolerance
!scalars
 integer :: idum=-5
 integer :: iimage,irestart,ndof,nnos,pitransform,prtstress
 real(dp) :: dtion,eharm,eharm2,epot,initemp,kt,temperature1,temperature2,thermtemp
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

!Compute cartesian coordinates
 do iimage=1,trotter
   xred_tmp => xred(:,:,iimage)      ! Needed because xred is inten(inout) in xredxcart
   call xredxcart(natom,1,rprimd,xcart     (:,:,iimage),xred_tmp)
   xred_tmp => xred_prev(:,:,iimage) ! Needed because xred_tmp is inten(inout) in xredxcart
   call xredxcart(natom,1,rprimd,xcart_prev(:,:,iimage),xred_tmp)
 end do

!Determine if it is a restart or not
!If this is a calculation from scratch,generate random distribution of velocities
 irestart=1;if (itimimage==1) irestart=pimd_is_restart(masseff,vel)
 if (irestart==0) then
   call pimd_initvel(idum,masseff,natom,initemp,trotter,vel)
 end if

!Compute temperature at t
 temperature1=pimd_temperature(masseff,vel)

!################## Images evolution #####################

!This is temporary
 xcart_next=zero
 rprimd_next=rprimd_prev
 temperature2=pimd_temperature(masseff,vel)

!Compute contributions to energy
 call pimd_energies(eharm,eharm2,epot,etotal,forces_orig,natom,springeff,trotter,xcart)

!Compute stress tensor at t from virial theorem
 call pimd_stresses(forces_pimd,masseff,natom,stress_pimd,trotter,vel,volume,xcart)


!############# Final operations ############################

!Print messages
 prtstress=1
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
 
end subroutine pimd_nosehoover_npt
!!***


