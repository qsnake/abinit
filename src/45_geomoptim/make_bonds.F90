!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_bonds
!! NAME
!! make_bonds
!!
!! FUNCTION
!!  (to be completed)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JJ)
!!
!! INPUTS
!!  (to be completed)
!!
!! OUTPUT
!!  (to be completed)
!!
!! PARENTS
!!      make_prim_internals
!!
!! CHILDREN
!!      atmdata
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine make_bonds(deloc,natom,ntypat,icenter,rprimd,typat,xcart,znucl)

 use m_profiling

 use defs_basis
 use defs_mover
 use m_delocint

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_bonds'
 use interfaces_32_util
 use interfaces_45_geomoptim, except_this_one => make_bonds
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icenter
 integer,intent(in) :: natom,ntypat
 type(ab_delocint),intent(inout) :: deloc
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: znucl(:) ! znucl(ntypat) or
                                 ! znucl(npsp) ?
 real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)

!Local variables ------------------------------
!scalars
 integer :: iatom,ibond,irshift,itypat,jatom
 real(dp) :: amu,bl,bondfudge,rcov1,rcov2
 character(len=2) :: symbol
!arrays
 integer,allocatable :: bonds_tmp(:,:,:)
 real(dp) :: rcov(ntypat),rpt(3)

!************************************************************************

!DEBUG
!write(std_out,*) 'make_bonds: enter'
!ENDDEBUG
 do itypat=1,ntypat
   call atmdata(amu,rcov(itypat),symbol,znucl(itypat))
 end do

!DEBUG
!write(std_out,*) ' rcov =', rcov
!write(std_out,*) ' nrshift =', deloc%nrshift
!write(std_out,*) ' xcart =', xcart
!write(std_out,*) ' natom =',natom
!ENDDEBUG

!tentative first allocation: < 12 bonds per atom.
 ABI_ALLOCATE(bonds_tmp,(2,2,12*natom))

 bondfudge = 1.1_dp

 deloc%nbond = 0

 do iatom=1,natom
   rcov1 = rcov(typat(iatom))
   do jatom=iatom+1,natom
     rcov2 = rcov(typat(jatom))
     do irshift=1,deloc%nrshift
       rpt(:) = xcart(:,jatom) &
&       + deloc%rshift(1,irshift)*rprimd(:,1) &
&       + deloc%rshift(2,irshift)*rprimd(:,2) &
&       + deloc%rshift(3,irshift)*rprimd(:,3)
       bl =  bond_length(xcart(:,iatom),rpt)

!      DEBUG
!      write(std_out,*) ' bl, bondfudge*(rcov1+rcov2) = ',bl, bondfudge*(rcov1+rcov2)
!      ENDDEBUG

       if (bondfudge*(rcov1+rcov2) - bl > tol6) then
         deloc%nbond = deloc%nbond+1
         if (deloc%nbond > 12*natom) then
           write(std_out,*) 'make_bonds : error too many bonds !'
           stop
         end if
         bonds_tmp(1,1,deloc%nbond) = iatom
         bonds_tmp(2,1,deloc%nbond) = icenter
         bonds_tmp(1,2,deloc%nbond) = jatom
         bonds_tmp(2,2,deloc%nbond) = irshift

!        DEBUG
!        write(std_out,*) ' ibond bonds = ', deloc%nbond, bonds_tmp(:,:,deloc%nbond),xcart(:,iatom),rpt
!        ENDDEBUG

       end if
     end do
!    end jatom do
   end do
 end do
!end iatom do

 if (associated(deloc%bonds)) nullify(deloc%bonds)
 ABI_ALLOCATE(deloc%bonds,(2,2,deloc%nbond))
 do ibond=1,deloc%nbond
   deloc%bonds(:,:,ibond) = bonds_tmp(:,:,ibond)
 end do
 ABI_DEALLOCATE(bonds_tmp)

!DEBUG
!do ibond=1,deloc%nbond
!write(std_out,*) 'bond ', ibond, deloc%bonds(:,:,ibond)
!end do
!ENDDEBUG

end subroutine make_bonds
!!***
