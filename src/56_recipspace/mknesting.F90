!{\src2tex{textfont=tt}}
!!****f* ABINIT/mknesting
!! NAME
!! mknesting
!!
!! FUNCTION
!!  Calculate the nesting factor over the dense k-grid,
!!  interpolate the values along a given q path
!!  and write the data on file in the X-Y format or
!!  in the XCrysden format (XSF)
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nkpt = number of k points
!!  kpt(3,nkpt) = k points
!!  nkx, nky, nkz = number of k-point along each direction
!!  nband = number of bands to be considered in the calculation
!!  weight(nband,nkpt) =  integration weights for each k-point and band
!!  nqpath = number of points requested along the trajectory
!!  qpath_vertices = vertices of the reciprocal space trajectory
!!  base_name = prefix of the output file
!!  gprimd(3,3) dimensional reciprocal lattice vectors
!!  gmet = metric in reciprocal space
!!  prtnest = flags governing the format of the output file
!! OUTPUT
!!  only write to file
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      elphon,outscfcv
!!
!! CHILDREN
!!      bfactor,destroy_kptrank,mkkptrank,outnesting,wrap2_zero_one,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h" 

subroutine mknesting(nkpt,kpt,kptrlatt,nband,weight,nqpath,&
& qpath_vertices,nqptfull,qptfull,base_name,gprimd,gmet,prtnest,qptrlatt,&
& nsym, symrec)

 use m_profiling

 use defs_basis
 use m_errors
 use m_kptrank

 !use m_crystal,   only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mknesting'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_56_recipspace, except_this_one => mknesting
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nband,nkpt,nqpath,prtnest
 integer, intent(in) :: nqptfull
 integer, intent(in), optional :: nsym
 character(len=fnlen),intent(in) :: base_name
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 integer,intent(in),optional :: symrec(3,3,*)
 real(dp),intent(in) :: gprimd(3,3),kpt(3,nkpt)
 real(dp),intent(in) :: qptfull(3,nqptfull)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: qpath_vertices(3,nqpath)
 real(dp),intent(in) :: weight(nband,nkpt)
 integer,intent(in)  :: qptrlatt(3,3)

!Local variables-------------------------------
!scalars
 integer :: ikpt,jkpt,kindex,maxrank
 integer :: ik1, ik2, ik3, nkptfull
 real(dp) :: res
 character(len=500) :: message
 type(kptrank_type) :: kptrank_t
!arrays
 integer,allocatable :: kptrank(:),ktable(:)
 character(len=fnlen) :: tmpname
 real(dp) :: tmpkpt(3)
 real(dp),allocatable :: nestfactor(:),nestordered(:)
 real(dp), allocatable :: kptfull(:,:)

! *************************************************************************

 if (  kptrlatt(1,2) /= 0 .or. kptrlatt(1,3) /= 0 .or. kptrlatt(2,1) /= 0       &
& .or. kptrlatt(2,3) /= 0 .or. kptrlatt(3,1) /= 0 .or. kptrlatt(3,2) /= 0 ) then
   write (message,'(7a)')ch10,' mknesting : WARNING-',ch10,                         &
&   ' kptrlatt should be diagonal in order to calculate the nesting factor,',ch10,&
&   ' skipping the nesting factor calculation ',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   return
 end if

 if (prtnest /= 1 .and. prtnest /= 2) then
   MSG_BUG(' prtnest should be 1 or 2')
 end if

 write(message,'(a,9(i0,1x))')' mknesting : kptrlatt = ',kptrlatt
 call wrtout(std_out,message,'COLL')

 nkptfull = kptrlatt(1,1)*kptrlatt(2,2)*kptrlatt(3,3)
 ABI_ALLOCATE(nestordered,(nkptfull))
 nestordered(:)=zero
 ABI_ALLOCATE(kptfull,(3,nkptfull))
 ikpt = 0
 do ik3 = 0, kptrlatt(3,3)-1
   do ik2 = 0, kptrlatt(2,2)-1
     do ik1 = 0, kptrlatt(1,1)-1
       ikpt = ikpt+1
       kptfull(:,ikpt) = (/dble(ik1)/dble(kptrlatt(1,1)), dble(ik2)/dble(kptrlatt(2,2)),dble(ik3)/dble(kptrlatt(3,3))/)
     end do
   end do
 end do

!NOTE: input weights are not normalised, the normalisation factor in introduced in bfactor
!new version now puts kptfull in correct order before bfactor, so no need to re-order...
 if (present(symrec)) then
   ABI_CHECK(present(nsym), "error - provide nsym and symrec arguments together")
   call mkkptrank (kpt,nkpt,kptrank_t, nsym, symrec)
 else
   call mkkptrank (kpt,nkpt,kptrank_t)
 end if

 call bfactor(nkptfull,kptfull,nkptfull,kptfull,kptrank_t,nkpt,weight,nband,nestordered)

!================================================================================================
!use linear interpolation to plot the bfactor along the given q-path
!1) order the kpoints of the grid putting them in increasing x, then y, then z (FORTRAN convention)
!2) make table from input kpts to ordered kpts
!3) perform interpolation
!================================================================================================

 call outnesting(base_name,gmet,gprimd,kptrlatt,nestordered,nkptfull,nqpath,prtnest,qpath_vertices)
 ABI_DEALLOCATE(nestordered)


!
!now do the same, but for the nesting factor over the phonon qpoints only
!
 ABI_ALLOCATE(nestfactor,(nqptfull))
 nestfactor(:)=zero
 call bfactor(nkptfull,kptfull,nqptfull,qptfull,kptrank_t,nkpt,weight,nband,nestfactor)

 call destroy_kptrank (kptrank_t)
 ABI_DEALLOCATE(kptfull)

!rank is used to order kpoints
 ABI_ALLOCATE(kptrank,(nqptfull))
 kptrank(:) = 0

 do ikpt=1,nqptfull
   call wrap2_zero_one(qptfull(1,ikpt),tmpkpt(1),res)
   call wrap2_zero_one(qptfull(2,ikpt),tmpkpt(2),res)
   call wrap2_zero_one(qptfull(3,ikpt),tmpkpt(3),res)
   kptrank(ikpt) = 100000000.0_dp*(tmpkpt(3)+one) + &
&   100000.0_dp*(tmpkpt(2)+one) + &
&   100.0_dp*(tmpkpt(1)+one)
 end do

 ABI_ALLOCATE(ktable,(nqptfull))
 ktable(:)=0

 kindex=nqptfull
 do ikpt=1,nqptfull
!  FIXME: this whole operation could be replaced by a sort(kptrank)
!  and then using the sorted indirect indexing ignoring the -1 rank values
   maxrank=maxval(kptrank)
   findmax2: do jkpt=1,nqptfull
     if(kptrank(jkpt)==maxrank) then
       kptrank(jkpt)=-kptrank(jkpt)
       ktable(jkpt)=kindex
       kindex=kindex-1
       exit findmax2
     end if
   end do findmax2
 end do !ikpt
 ABI_DEALLOCATE(kptrank)

!fill the datagrid for the nesting factor using the Fortran convention and the conventional unit cell
!NOTE: the Fortran convention is a must if we want to plot the data
!in the BXSF format, useful for the linear interpolation since we use interpol3d.F90

 ABI_ALLOCATE(nestordered,(nqptfull))
 nestordered(:)=zero
 do jkpt=1,nqptfull
   ikpt = ktable(jkpt)
   nestordered(ikpt)=nestfactor(jkpt)
 end do
 ABI_DEALLOCATE(nestfactor)
 ABI_DEALLOCATE(ktable)

 tmpname = trim(base_name)//"kplusq"
 call outnesting(tmpname,gmet,gprimd,qptrlatt,nestordered,nqptfull,nqpath,prtnest,qpath_vertices)

 ABI_DEALLOCATE(nestordered)

end subroutine mknesting
!!***

