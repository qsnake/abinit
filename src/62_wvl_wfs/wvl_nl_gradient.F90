!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_nl_gradient
!! NAME
!! wvl_nl_gradient
!!
!! FUNCTION
!! Compute the non local part of the wavefunction gradient.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      afterscfloop,vtorho
!!
!! CHILDREN
!!      leave_new,nonlocal_forces,wrtout,xcomm_world,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_nl_gradient(grnl, mpi_enreg, natom, rprimd, wvl, xcart)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only: nonlocal_forces
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_nl_gradient'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 type(MPI_type),intent(inout) :: mpi_enreg
 type(wvl_data),intent(inout) :: wvl
!arrays
 real(dp),intent(in) :: xcart(3,natom),rprimd(3,3)
 real(dp),intent(inout) :: grnl(3,natom)

!Local variables-------------------------------
!scalars
 integer :: ia,ierr,igeo,me,nproc,spaceComm
 character(len=500) :: message
!arrays
 real(dp),allocatable :: gxyz(:,:)

! *************************************************************************

!Compute forces
 write(message, '(a,a)' ) 'wvl_nl_gradient(): compute non-local part to gradient.'
 call wrtout(std_out,message,'COLL')

!Nullify output arrays.
 grnl(:, :) = zero
 
 ABI_ALLOCATE(gxyz,(3, natom))
 gxyz(:,:) = zero
!Add the nonlocal part of the forces to grtn (BigDFT routine)
#if defined HAVE_DFT_BIGDFT
 call xcomm_world(mpi_enreg,spaceComm,myrank=me,mysize=nproc)
 call nonlocal_forces(me, wvl%descr%Glr%d%n1, wvl%descr%Glr%d%n2, wvl%descr%Glr%d%n3, &
& wvl%descr%h(1), wvl%descr%h(2), wvl%descr%h(3), wvl%descr%atoms, &
& xcart, wvl%wfs%orbs, wvl%projectors%keys, wvl%projectors%proj, wvl%wfs%Glr%wfd, &
& wvl%wfs%psi, gxyz, .true.)
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_nl_gradient: BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,message,'COLL')
 call leave_new('COLL')
#endif

 if (nproc > 1) then
   call xsum_mpi(gxyz, spaceComm, ierr)
 end if

!Forces should be in reduced coordinates.
 do ia = 1, natom, 1
   do igeo = 1, 3, 1
     grnl(igeo, ia) = - rprimd(1, igeo) * gxyz(1, ia) - &
&     rprimd(2, igeo) * gxyz(2, ia) - &
&     rprimd(3, igeo) * gxyz(3, ia)
   end do
 end do
 ABI_DEALLOCATE(gxyz)

end subroutine wvl_nl_gradient
!!***
