!{\src2tex{textfont=tt}}
!!****f* ABINIT/mklocl_wavelets
!!
!! NAME
!! mklocl_wavelets
!!
!! FUNCTION
!! Compute the ionic local potential when the pseudo-potentials are GTH, using
!! the special decomposition of these pseudo. The resulting potential is computed with
!! free boundary conditions. It gives the same result than mklocl_realspace for the
!! GTH pseudo only with a different way to compute the potential.
!!
!! Optionally compute :
!!  option=1 : local ionic potential throughout unit cell
!!  option=2 : contribution of local ionic potential to E gradient wrt xred
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=internal variables used by wavelets, describing
!!   | wvl=desciption of the wavelet box.
!!   | natom=number of atoms.
!!  mpi_enreg=informations about MPI parallelization
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! OUTPUT
!!  (if option==1) vpsp(dtset%nfft)=the potential resulting from the ionic
!!                 density of charge.
!!  (if option==2) grtn(3,natom)=grads of Etot wrt tn. These gradients are in
!!                 reduced coordinates. Multiply them by rprimd to get
!!                 gradients in cartesian coordinates.
!!
!! SIDE EFFECTS
!!  xred(3,natom)=reduced dimensionless atomic coordinates (in fact IN but here
!!                because of INOUT xredxcart() behavior).
!!
!! PARENTS
!!      mklocl
!!
!! CHILDREN
!!      createionicpotential,leave_new,local_forces,psolver_hartree
!!      psolver_kernel,wrtout,xcomm_world,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mklocl_wavelets(dtset, grtn, mpi_enreg, option, rhor, rprimd, &
     & vpsp, wvl, xcart)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only : createIonicPotential, local_forces
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mklocl_wavelets'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
 use interfaces_62_poisson
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 real(dp),intent(in) :: rhor(dtset%nfft,dtset%nspden),rprimd(3,3)
 real(dp),intent(inout) :: grtn(3,dtset%natom),vpsp(dtset%nfft)
 real(dp),intent(inout) :: xcart(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: i,i1,i2,i3,ia,ierr,igeo,me,nproc,shift,spaceComm
 real(dp) :: energy
 character(len=500) :: message
 character(len = 1) :: bndcode
!arrays
 real(dp) :: epot(3)
 real(dp),allocatable :: gxyz(:,:),vhartr(:)
 real(dp),pointer :: kernel(:)
 character(len=20) :: atomnames(100)

! *********************************************************************

 if (dtset%icoulomb == 0) then
   bndcode = 'P'
 else if (dtset%icoulomb == 1) then
   bndcode = 'F'
 else if (dtset%icoulomb == 2) then
   bndcode = 'S'
 end if


#if defined HAVE_DFT_BIGDFT

 call xcomm_world(mpi_enreg,spaceComm,myrank=me,mysize=nproc)

!We get the kernel for the Poisson solver (used to go from the ionic
!charge to the potential or to compute the Hartree potential).
!If the kernel is uncomputed, it does it now.
 call PSolver_kernel(dtset, 2, kernel, mpi_enreg, rprimd, wvl)

 shift = 1 + wvl%Glr%d%n1i * wvl%Glr%d%n2i * &
& mpi_enreg%nscatterarr(me, 4)
 if (option == 1) then
   write(message, '(a,a)' ) ch10,&
&   ' mklocl_wavelets: Create local potential from ions.'
   call wrtout(std_out,message,'COLL')

!  Call the BigDFT routine...
   call createIonicPotential(bndcode, me, nproc, wvl%atoms, &
&   xcart, wvl%h(1) / 2., wvl%h(2) / 2., wvl%h(3) / 2., &
&   real(0, dp), wvl%Glr%d%n1, wvl%Glr%d%n2, wvl%Glr%d%n3, &
&   mpi_enreg%ngfft3_ionic, mpi_enreg%nscatterarr(me, 3) + 1, &
&   wvl%Glr%d%n1i, wvl%Glr%d%n2i, wvl%Glr%d%n3i, &
&   kernel, vpsp(shift), 0.d0, 0, .false.)

   if (maxval(dtset%efield) > zero) then
     write(message, '(a,a)' ) ch10,&
&     ' mklocl_wavelets: Add the electric field.'
     call wrtout(std_out,message,'COLL')

!    We add here the electric field since in BigDFT, the field must be on x...
     epot(:) = real(0.25, dp) * dtset%efield(:) * dtset%wvl_hgrid
     do i3 = 1, mpi_enreg%ngfft3_ionic, 1
       ia = (i3 - 1) * wvl%Glr%d%n1i * wvl%Glr%d%n2i
       do i2 = -14, 2 * wvl%Glr%d%n2 + 16, 1
         i = ia + (i2 + 14) * wvl%Glr%d%n1i
         do i1 = -14, 2 * wvl%Glr%d%n1 + 16, 1
           i = i + 1
           vpsp(shift + i) = vpsp(shift + i) + &
&           epot(1) * real(i1 - wvl%Glr%d%n1, dp) + &
&           epot(2) * real(i2 - wvl%Glr%d%n2, dp) + &
&           epot(3) * real(i3 - wvl%Glr%d%n3, dp)
         end do
       end do
     end do
   end if

 else if (option == 2) then
!  Dummy arguments
   atomnames(:) = "'Unknown name'"

!  Compute forces
   write(message, '(a)' ) 'mklocl_wavelets: compute local forces.'
   call wrtout(std_out,message,'COLL')

!  Compute Hartree's potential from rhor.
   ABI_ALLOCATE(vhartr,(dtset%nfft))
   call PSolver_hartree(dtset, energy, mpi_enreg, rhor, rprimd, vhartr, wvl)

!  Allocate temporary array for forces.
   ABI_ALLOCATE(gxyz,(3, dtset%natom))

!  calculate local part of the forces grtn (BigDFT routine)
   call local_forces(me, wvl%atoms, xcart, &
&   wvl%h(1) / 2., wvl%h(2) / 2., wvl%h(3) / 2., &
&   wvl%Glr%d%n1, wvl%Glr%d%n2, wvl%Glr%d%n3, &
&   mpi_enreg%nscatterarr(me, 2), mpi_enreg%nscatterarr(me, 3) + 1, &
&   wvl%Glr%d%n1i, wvl%Glr%d%n2i, wvl%Glr%d%n3i, rhor(shift, 1), &
&   vhartr(shift), gxyz)
   ABI_DEALLOCATE(vhartr)

   if (nproc > 1) then
     call xsum_mpi(gxyz, spaceComm, ierr)
   end if

!  Forces should be in reduced coordinates.
   do ia = 1, dtset%natom, 1
     do igeo = 1, 3, 1
       grtn(igeo, ia) = - rprimd(1, igeo) * gxyz(1, ia) - &
&       rprimd(2, igeo) * gxyz(2, ia) - &
&       rprimd(3, igeo) * gxyz(3, ia)
     end do
   end do

!  Deallocate local variables
   ABI_DEALLOCATE(gxyz)
 else ! option switch
   write(message, '(a,a,a,a)' ) ch10,&
&   ' mklocl_wavelets : internal error, option should be 1 or 2.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' mklocl_wavelets : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif
end subroutine mklocl_wavelets
!!***
