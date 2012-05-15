!{\src2tex{textfont=tt}}
!!****f* ABINIT/psps_init_global
!! NAME
!! psps_init_global
!!
!! FUNCTION
!! Allocate and initialise all part of psps structure that are independent
!! of a given dataset.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! npsp=the number of read pseudo files.
!! pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! psps=<type pseudopotential_type>the pseudopotentials description
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      psp2params_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psps_init_global(mtypalch, npsp, psps, pspheads)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psps_init_global'
 use interfaces_65_psp, except_this_one => psps_init_global
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mtypalch,npsp
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
!scalars
 integer :: ii, mpsang, n1xccc

! *************************************************************************

!Allocation of some arrays independent of the dataset
 ABI_ALLOCATE(psps%filpsp,(npsp))
 ABI_ALLOCATE(psps%pspcod,(npsp))
 ABI_ALLOCATE(psps%pspdat,(npsp))
 ABI_ALLOCATE(psps%pspso,(npsp))
 ABI_ALLOCATE(psps%pspxc,(npsp))
 ABI_ALLOCATE(psps%title,(npsp))
 ABI_ALLOCATE(psps%zionpsp,(npsp))
 ABI_ALLOCATE(psps%znuclpsp,(npsp))
 call psp2params_init(psps%gth_params, npsp)

 psps%filpsp(1:npsp)=pspheads(1:npsp)%filpsp
 psps%pspcod(1:npsp)=pspheads(1:npsp)%pspcod
 psps%pspdat(1:npsp)=pspheads(1:npsp)%pspdat
 psps%pspso(1:npsp)=pspheads(1:npsp)%pspso
 psps%pspxc(1:npsp)=pspheads(1:npsp)%pspxc
 psps%title(1:npsp)=pspheads(1:npsp)%title
 psps%zionpsp(1:npsp)=pspheads(1:npsp)%zionpsp
 psps%znuclpsp(1:npsp)=pspheads(1:npsp)%znuclpsp

!Set values independant from dtset
 psps%npsp   = npsp
!Note that mpsang is the max of 1+lmax, with minimal value 1 (even for local psps, at present)
!mpsang=max(maxval(pspheads(1:npsp)%lmax)+1,1) ! might not work with HP compiler
!n1xccc=maxval(pspheads(1:npsp)%xccc)
 mpsang=1
 n1xccc=pspheads(1)%xccc
 do ii=1,psps%npsp
   mpsang=max(pspheads(ii)%lmax+1,mpsang)
   n1xccc=max(pspheads(ii)%xccc,n1xccc)
 end do
 psps%mpsang = mpsang
 psps%n1xccc = n1xccc
!Determine here whether the calculation is PAW
 psps%usepaw  =0
 if (pspheads(1)%pspcod==7) psps%usepaw=1  ! If paw, all pspcod necessarily are 7 (see iofn2)
 psps%mtypalch = mtypalch
end subroutine psps_init_global
!!***
