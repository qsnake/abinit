!{\src2tex{textfont=tt}}
!!****f* ABINIT/pimd_langevin_npt
!! NAME
!! pimd_langevin_npt
!!
!! FUNCTION
!! Predicts new positions in Path Integral Molecular Dynamics using Langevin thermostat in the NPT ensemble.
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
!!    at output, electronic forces + quantum spring contribution
!!  vel(3,natom,trotter)=velocies of atoms for all images
!!    at input,  values at time t
!!    at output, values at time t+dt
!!
!! NOTES
!!  Here follows PIMD in the NPT ensemble within the Langevin barostat algorithm
!1  of Quigley and Probert: J. Chem. Phys. 120, 11432 (2004)
!!  and Comput. Phys. Comm. 169, 322 (2005)

!!
!! PARENTS
!!      predict_pimd
!!
!! CHILDREN
!!      matr3inv,pimd_energies,pimd_forces,pimd_initvel,pimd_langevin_forces
!!      pimd_langevin_random,pimd_langevin_random_bar,pimd_langevin_random_init
!!      pimd_mass_spring,pimd_predict_taylor,pimd_print,pimd_stresses,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pimd_langevin_npt(etotal,forces,itimimage,natom,pimd_param,prtvolimg,&
&          rprimd,rprimd_next,rprimd_prev,trotter,vel,volume,xred,xred_next,xred_prev)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_pimd
 use m_random_zbq

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pimd_langevin_npt'
 use interfaces_32_util
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
!        Option for the Langevin algorithm correction
 integer,parameter :: ilangevin=0
!        The following option forces the total of forces to be zero
!         It prevents the translation of the center of mass
!         If it is zero, no constraint on mass center is applied
 integer,parameter :: zeroforce=1
!        Tolerance for the SC cycle
 real(dp),parameter :: tolerance=tol7

!scalars
 integer :: idum=-5
 integer :: iatom,ii,iimage,irestart,jj,ndof,prtstress
 real(dp) :: dtion,eharm,eharm2,epot,friction,frictionbar,initemp,kt,scalebar
 real(dp) :: temperature1,temperature2,temp2_prev,thermtemp,tol,tracepg,wg
!arrays
 real, parameter :: identity(3,3)=reshape((/(one,(zero,ii=1,3),jj=1,2),one/),(/3,3/))
 real(dp) :: aleabar(3,3),ddh(3,3),diffstress(3,3),dstrhh(3,3),fg(3,3),invrprimd(3,3)
 real(dp) :: langev_bar(3,3),pg(3,3),pgdh(3,3),stress_pimd(3,3),strtarget(6),tmp(3,3)
 real(dp),allocatable :: alea(:,:,:),forces_orig(:,:,:),forces_pimd(:,:,:),forces_pimd_red(:,:)
 real(dp),allocatable :: fsup(:,:),hxredpoint(:,:,:),inertmass(:),langev(:)
 real(dp),allocatable ::  mass(:,:),quantummass(:),spring(:,:)
 real(dp),allocatable :: xcart(:,:,:),xcart_next(:,:,:),xcart_prev(:,:,:)
 real(dp),allocatable :: xredpoint(:,:,:)
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
 dtion=pimd_param%dtion
 kt=thermtemp*kb_HaK
 friction=pimd_param%vis
 wg=pimd_param%bmass
 strtarget(:)=pimd_param%strtarget(:) ! imposed stress tensor
 frictionbar=pimd_param%friction      ! friction coeff of barostat
 scalebar=sqrt(two*frictionbar*wg*kt/dtion)
 forces_orig=forces

!Masses and spring constants
 ABI_ALLOCATE(mass,(natom,1))
 ABI_ALLOCATE(spring,(natom,1))
 call pimd_mass_spring(inertmass,kt,mass,natom,quantummass,spring,0,trotter)

!Correction Langevin algorithm bias
 ABI_ALLOCATE(alea,(3,natom,trotter))
 ABI_ALLOCATE(langev,(natom))
 langev(:)=sqrt(two*friction*inertmass(:)*kt/dtion)
 if(ilangevin==1)then
   langev(:)=langev(:)*sqrt(one-(friction*dtion/(two*inertmass(:))))
 end if

!Random number generator initialization
 if(itimimage<=1) then
   call pimd_langevin_random_init(pimd_param%irandom,idum)
 end if

!Compute cartesian coordinates
 do iimage=1,trotter
   xred_tmp => xred(:,:,iimage)      ! Needed because xred is inten(inout) in xredxcart
   call xredxcart(natom,1,rprimd,xcart     (:,:,iimage),xred_tmp)
   xred_tmp => xred_prev(:,:,iimage) ! Needed because xred_tmp is inten(inout) in xredxcart
   call xredxcart(natom,1,rprimd,xcart_prev(:,:,iimage),xred_tmp)
 end do

!Determine if it is a restart or not
!If this is a calculation from scratch,generate random distribution of velocities
 irestart=1;if (itimimage==1) irestart=pimd_is_restart(mass,vel)
 if (irestart==0) then
   call pimd_initvel(idum,mass,natom,initemp,trotter,vel)
 end if

 if (itimimage<=1) then

!  ========================= FIRST TIME STEP =======================================

!  Compute temperature at t
   temperature1=pimd_temperature(mass,vel)

!  Generate random numbers
   call pimd_langevin_random(alea,pimd_param%irandom,idum,langev,mass,natom,trotter,zeroforce)

!  Compute PIMD and Langevin contributions to forces
   call pimd_forces(forces,natom,spring,0,trotter,xcart)
   call pimd_langevin_forces(alea,forces,forces_pimd,friction,langev,mass,natom,trotter,vel)

!  Predict positions
   call pimd_predict_taylor(dtion,forces_pimd,mass,natom,trotter,&
&   vel,xcart,xcart_next)

!  Compute stress tensor at t from virial theorem
   call pimd_stresses(forces_pimd,mass,natom,stress_pimd,trotter,vel,volume,xcart)

!  Estimate the velocities at t+dt/2
   vel=(xcart_next-xcart)/dtion

!  Compute new temperature
   temperature2=pimd_temperature(mass,vel)

!  rprimd at time t+dt
   rprimd_next=rprimd

 else

!  ========================= OTHER TIME STEPS ======================================

!  Additional allocations
   ABI_ALLOCATE(hxredpoint,(3,natom,trotter))
   ABI_ALLOCATE(xredpoint,(3,natom,trotter))
   ABI_ALLOCATE(forces_pimd_red,(3,natom))
   ABI_ALLOCATE(fsup,(3,natom))

!  Initialization of Verlet algorithm
   if (itimimage==2) then
     ddh(:,:)=zero;pg(:,:)=zero
   end if
   tracepg=pg(1,1)+pg(2,2)+pg(3,3)

!  Momenta hxredpoint = H ds/dt
   if (itimimage==2) then
     hxredpoint=vel
   else
     do iimage=1,trotter
       hxredpoint(:,:,iimage)=matmul(rprimd,vel(:,:,iimage))
     end do
   end if

!  Compute temperature at t
   temperature1=pimd_temperature(mass,hxredpoint)

!  Generate random numbers
   call pimd_langevin_random(alea,pimd_param%irandom,idum,langev,mass,natom,trotter,zeroforce)

!  Generate random numbers for the barostat
   call pimd_langevin_random_bar(aleabar,pimd_param%irandom,idum)
   langev_bar=matmul(aleabar,rprimd)*scalebar

!  Compute PIMD and Langevin contributions to forces
   call pimd_forces(forces,natom,spring,0,trotter,xcart)
   call pimd_langevin_forces(alea,forces,forces_pimd,friction,langev,mass,natom,trotter,vel)

!  Compute inverse of matrix rprimd
   call matr3inv(rprimd,tmp)
   invrprimd=tmp
   invrprimd(1,2)=tmp(2,1);invrprimd(2,1)=tmp(1,2)
   invrprimd(1,3)=tmp(3,1);invrprimd(3,1)=tmp(1,3)
   invrprimd(2,3)=tmp(3,2);invrprimd(3,2)=tmp(2,3)

!  Compute difference between instantaneous stress and imposed stress (barostat)
   call pimd_stresses(forces_pimd,mass,natom,stress_pimd,trotter,vel,volume,xcart)
   diffstress=pimd_diff_stress(stress_pimd,strtarget)

!  Compute "force" on supercell vectors
   dstrhh=matmul(diffstress,rprimd)
   pgdh=dstrhh+matmul(pg,ddh/wg)
   fg(:,:)=volume*diffstress(:,:)+temperature1*kb_HaK*identity(:,:)+wg*pgdh(:,:)

!  Evolve the supercell (fist estimation)
   rprimd_next=two*rprimd-rprimd_prev+(fg-wg*frictionbar*ddh+langev_bar)*dtion*dtion/wg

!  Evolve atomic positions (first estimation)
   tmp=matmul(invrprimd,ddh)
   do iimage=1,trotter
     call xredxcart(natom,-1,rprimd,hxredpoint(:,:,iimage),xredpoint(:,:,iimage))
     fsup(:,:)=matmul(tmp,xredpoint(:,:,iimage))
     call xredxcart(natom,-1,rprimd,forces_pimd(:,:,iimage),forces_pimd_red)
     do iatom=1,natom
       xred_next(:,iatom,iimage)= &
&       two*xred(:,iatom,iimage) - xred_prev(:,iatom,iimage) &
&       +(forces_pimd_red(:,iatom)-two*inertmass(iatom)*fsup(:,iatom) &
&       -tracepg*inertmass(iatom)*xredpoint(:,iatom,iimage)/(wg*dble(ndof))) &
&       *dtion*dtion/inertmass(iatom)
     end do
   end do

!  Self-consistent loop
   temperature2=pimd_temperature(mass,vel)
   temp2_prev=temperature2; tol=one
   do while (tol>tolerance)
!    Reestimate dH/dt at t
     ddh(:,:)=half*(rprimd_next(:,:)-rprimd_prev(:,:))/dtion
!    Reestimate the scaled velocities at t
     do iimage=1,trotter
       xredpoint(:,:,iimage)=(xred_next(:,:,iimage)-xred_prev(:,:,iimage))/dtion
       call xredxcart(natom,1,rprimd,hxredpoint(:,:,iimage),xredpoint(:,:,iimage))
     end do
!    Reestimate the forces
     call pimd_langevin_forces(alea,forces,forces_pimd,friction,langev,mass,natom,trotter,hxredpoint)
!    Compute variation of temperature (to check convergence of SC loop)
     temperature2=pimd_temperature(mass,xredpoint)
     tol=dabs(temperature2-temp2_prev)/dabs(temp2_prev)
     temp2_prev=temperature2
!    Recompute the new temperature
     temperature2=pimd_temperature(mass,hxredpoint)
!    Recompute pg
     pg=matmul(ddh,invrprimd)*wg
     tracepg=pg(1,1)+pg(2,2)+pg(3,3)
!    Recompute difference between instantaneous stress and imposed stress (barostat)
     call pimd_stresses(forces,mass,natom,stress_pimd,trotter,hxredpoint,volume,xcart)
     diffstress=pimd_diff_stress(stress_pimd,strtarget)
!    Recompute "force" on supercell vectors
     dstrhh=matmul(diffstress,rprimd)
     pgdh=dstrhh+matmul(pg,ddh/wg)
     fg(:,:)=volume*diffstress(:,:)+temperature2*kb_HaK*identity(:,:)+wg*pgdh(:,:)
!    Evolve the supercell (first estimation)
     rprimd_next=two*rprimd-rprimd_prev+(fg-wg*frictionbar*ddh+langev_bar)*dtion*dtion/wg
!    Evolve atomic positions (better estimation):
     tmp=matmul(invrprimd,ddh)
     do iimage=1,trotter
       call xredxcart(natom,-1,rprimd,hxredpoint(:,:,iimage),xredpoint(:,:,iimage))
       fsup(:,:)=matmul(tmp,xredpoint(:,:,iimage))
       call xredxcart(natom,-1,rprimd,forces_pimd(:,:,iimage),forces_pimd_red)
       do iatom=1,natom
         xred_next(:,iatom,iimage)= &
&         two*xred(:,iatom,iimage) - xred_prev(:,iatom,iimage) &
&         +(forces_pimd_red(:,iatom)-two*inertmass(iatom)*fsup(:,iatom) &
&         -tracepg*inertmass(iatom)*xredpoint(:,iatom,iimage)/(wg*dble(ndof))) &
&         *dtion*dtion/inertmass(iatom)
       end do
     end do
   end do ! End self-consistent loop

!  Computation of true temperature from true velocities
   do iimage=1,trotter
     vel(:,:,iimage)=hxredpoint(:,:,iimage)+matmul(ddh,xred(:,:,iimage))
   end do
   temperature2=pimd_temperature(mass,vel)

!  Computation of the real stress tensor at t
   call pimd_stresses(forces_pimd,mass,natom,stress_pimd,trotter,vel,volume,xcart)

!  Estimate the velocities at t+dt
   vel = (three*xred_next - four*xred + xred_prev)/(two * dtion)

!  Deallocations (Verlet algo)
   ABI_DEALLOCATE(xredpoint)
   ABI_DEALLOCATE(hxredpoint)
   ABI_DEALLOCATE(forces_pimd_red)
   ABI_DEALLOCATE(fsup)

 end if ! itimimage==1

!############# Final operations ############################

!Compute contributions to energy
 call pimd_energies(eharm,eharm2,epot,etotal,forces_orig,natom,spring,trotter,xcart)

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
 ABI_DEALLOCATE(mass)
 ABI_DEALLOCATE(spring)
 ABI_DEALLOCATE(alea)
 ABI_DEALLOCATE(langev)

end subroutine pimd_langevin_npt
!!***
