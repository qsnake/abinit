!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawnabla_init
!! NAME
!! pawnabla_init
!!
!! FUNCTION
!! Evaluate onsite contributions of the nabla operator in cartesian coordinates.
!! Store values in the pawtab% data structure.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (SM,VR,FJ,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  indlmn(6,lmnmax,ntypat)=Array giving l,m,n,lm,ln,s for i=lmn
!!  lmnmax=Max number of (l,m,n) numbers over all types of atom
!!  mpsang=1+maximum angular momentum
!!  ntypat=Number of types of atoms in cell
!!  Pawrad(ntypat)<Pawrad_type>=PAW radial mesh and related data:
!!    %mesh_size=Dimension of radial mesh
!!    %rad(mesh_size)=The coordinates of all the points of the radial mesh
!!  Pawtab(ntypat) <type(pawtab_type>=PAW tabulated starting data:
!!    %lmn_size=Number of (l,m,n) elements for the PAW basis
!!
!! OUTPUT
!!  See side effects
!!
!! SIDE EFFECTS
!!  Pawtab(ntypat) <type(pawtab_type>=PAW tabulated starting data:
!!    %has_nabla=set to 1 in matrix elements are calculated and stored
!!    %nabla_ij(3,lmn_size,lmn_size)= <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j>
!!
!! NOTES
!!  MG extracted this piece of code from optics_paw.F90 in order to have something more 
!!  reusable! Note however the storage mode of nabla_ij differs from optics_paw 
!!  (here Cartesian coordinates run faster). Besides nabla_ij contains the matrix 
!!  elements of \nabla instead of the elements of the momentum operator p.
!!
!! PARENTS
!!      bethe_salpeter,screening
!!
!! CHILDREN
!!      deducer0,int_ang,nderiv_gen,simp_gen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawnabla_init(mpsang,lmnmax,ntypat,indlmn,pawrad,pawtab)

 use m_profiling
    
 use defs_basis
 use defs_datatypes
 use m_errors

 use m_radmesh,          only : simp_gen, nderiv_gen, deducer0

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawnabla_init'
 use interfaces_66_paw, except_this_one => pawnabla_init
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmnmax,mpsang,ntypat
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat)
 type(pawtab_type),intent(inout) :: pawtab(ntypat)
 type(pawrad_type),intent(in) :: pawrad(ntypat)

!Local variables-------------------------------
!scalars
 integer :: nln,il,ilm,ilmn,iln,itypat
 integer :: jl,jlm,jlmn,jln,lmn_size,mesh_size 
 real(dp) :: intg
 character(len=500) :: msg      
!arrays
 integer,allocatable :: indlmn_(:,:)
 real(dp) :: ang_phipphj(mpsang**2,mpsang**2,8)
 real(dp),allocatable :: dphi(:),dtphi(:),ff(:),int1(:,:),int2(:,:),rad(:)
 
! *************************************************************************

 DBG_ENTER("COLL")

 if (mpsang>4)then
   write(msg,'(3a)')&
&   ' Not designed for angular momentum greater than 3 ',ch10,&
&   ' Modification in the table defined in ang_int.F90 is required. '
   MSG_ERROR(msg)
 end if
!
!Integration of the angular part: all angular integrals have been computed 
!outside Abinit and tabulated for each (l,m) value up to l=2
 call int_ang(ang_phipphj,mpsang)

 do itypat=1,ntypat
!  
!  COMPUTE nabla_ij := <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j> for this type
   mesh_size=pawrad(itypat)%mesh_size
   lmn_size=pawtab(itypat)%lmn_size
   nln=pawtab(itypat)%basis_size

!  $ if (associated(pawtab(itypat)%nabla_ij)) deallocate(pawtab(itypat)%nabla_ij)
   ABI_ALLOCATE(pawtab(itypat)%nabla_ij,(3,lmn_size,lmn_size))
   pawtab(itypat)%has_nabla=1

   ABI_ALLOCATE(ff,(mesh_size))
   ABI_ALLOCATE(rad,(mesh_size))
   ABI_ALLOCATE(int2,(lmn_size,lmn_size))
   ABI_ALLOCATE(int1,(lmn_size,lmn_size))
   ABI_ALLOCATE(dphi,(mesh_size))
   ABI_ALLOCATE(dtphi,(mesh_size))
   ABI_ALLOCATE(indlmn_,(6,lmnmax))

   indlmn_(:,:)=indlmn(:,:,itypat)
   rad(1:mesh_size)=pawrad(itypat)%rad(1:mesh_size)
!  
!  int1=\int phi phj/r dr - \int tphi tphj /r dr
   do jln=1,nln
     do iln=1,nln
       ff(2:mesh_size)= ( &
&       pawtab(itypat)%phi (2:mesh_size,iln)*pawtab(itypat)%phi (2:mesh_size,jln) &
&       -pawtab(itypat)%tphi(2:mesh_size,iln)*pawtab(itypat)%tphi(2:mesh_size,jln) ) /rad(2:mesh_size)
       call deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int1(iln,jln)=intg
     end do
   end do
!  
!  int2=\int phi/r d/dr(phj/r) r^2dr - \int tphi/r d/dr(tphj/r)r^2 dr
   do jln=1,nln
     ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,jln)
     call nderiv_gen(dphi,ff,1,pawrad(itypat))
     ff(1:mesh_size)=pawtab(itypat)%tphi(1:mesh_size,jln)
     call nderiv_gen(dtphi,ff,1,pawrad(itypat))

     do iln=1,nln
       ff(2:mesh_size)= &
&       pawtab(itypat)%phi (2:mesh_size,iln)*dphi (2:mesh_size) &
&       -pawtab(itypat)%phi (2:mesh_size,iln)*pawtab(itypat)%phi (2:mesh_size,jln)/rad(2:mesh_size) &
&       -( pawtab(itypat)%tphi(2:mesh_size,iln)*dtphi(2:mesh_size) &
&       -pawtab(itypat)%tphi(2:mesh_size,iln)*pawtab(itypat)%tphi(2:mesh_size,jln)/rad(2:mesh_size) )
       call deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int2(iln,jln)=intg
     end do
   end do
!  
!  1-c Integration of the radial part, Note unpacked loop
   do jlmn=1,lmn_size
     jlm=indlmn_(4,jlmn)
     jl =indlmn_(5,jlmn)
     do ilmn=1,lmn_size
       ilm=indlmn_(4,ilmn)
       il =indlmn_(5,ilmn)

       pawtab(itypat)%nabla_ij(1,ilmn,jlmn)= &
&       int2(il,jl)* ang_phipphj(ilm,jlm,1) &
&       +int1(il,jl)*(ang_phipphj(ilm,jlm,2)+ang_phipphj(ilm,jlm,3))

       pawtab(itypat)%nabla_ij(2,ilmn,jlmn)= &
&       int2(il,jl)* ang_phipphj(ilm,jlm,4) &
&       +int1(il,jl)*(ang_phipphj(ilm,jlm,5)+ang_phipphj(ilm,jlm,6))

       pawtab(itypat)%nabla_ij(3,ilmn,jlmn)= &
&       int2(il,jl)* ang_phipphj(ilm,jlm,7) &
&       +int1(il,jl)* ang_phipphj(ilm,jlm,8)

     end do !ilmn
   end do !jlmn

   pawtab(itypat)%has_nabla=2
   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(rad)
   ABI_DEALLOCATE(int2)
   ABI_DEALLOCATE(int1)
   ABI_DEALLOCATE(dphi)
   ABI_DEALLOCATE(dtphi)
   ABI_DEALLOCATE(indlmn_)

 end do !itypat

 DBG_EXIT("COLL")

end subroutine pawnabla_init
!!***

