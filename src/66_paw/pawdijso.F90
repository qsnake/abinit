!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawdijso
!! NAME
!! pawdijso
!!
!! FUNCTION
!! Compute the spin-orbit contribution to the PAW
!! pseudopotential strength Dij
!! (for one atom only)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  iatom=index of current atom (note: this is the index on current proc, not the absolute index)
!!  itypat=type index of current atom
!!  natom=number of atoms in cell (note: this is the number on current proc, not the absolute number)
!!  ntypat=number of types of atoms in unit cell.
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!
!! OUTPUT
!!  paw_ij(iatom)%dijso(cplex_dij*lmn2_size,ndij)= spin-orbit Dij terms
!!  cplex_dij=2 must be 2
!!        dijso(:,:,1) contains Dij_SO^up-up
!!        dijso(:,:,2) contains Dij_SO^dn-dn
!!        dijso(:,:,3) contains Dij_SO^up-dn
!!        dijso(:,:,4) contains Dij_SO^dn-up
!!
!! PARENTS
!!      pawdenpot,pawdij
!!
!! CHILDREN
!!      deducer0,nderiv_gen,simp_gen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawdijso(iatom,itypat,natom,ntypat,paw_an,paw_ij,pawang,pawrad,pawtab,pawxcdev,spnorbscl)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

 use m_radmesh,          only : simp_gen, nderiv_gen, deducer0

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdijso'
!End of the abilint section

 implicit none
!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: iatom,itypat,natom,ntypat,pawxcdev
 real(dp), intent(in) :: spnorbscl
 type(pawang_type),intent(in) :: pawang
!arrays
 type(paw_an_type),intent(in) :: paw_an(natom)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: cplex_dij,idij,ij_size,ilm,ipts,ispden,jlm,klm,klmn,klmn1,kln,lmn2_size,mesh_size,nsploop
 real(dp), parameter :: HalfFineStruct2=half/InvFineStruct**2
 real(dp) :: fact
 character(len=500) :: msg
!arrays
 integer,allocatable :: indklmn(:,:)
 real(dp),allocatable :: dijso_rad(:),dv1dr(:),ff(:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (paw_ij(iatom)%ndij/=4) then
   msg='invalid ndij size for Dij with spin-orbit coupling !'
   MSG_BUG(msg)
 end if
 if (paw_ij(iatom)%cplex_dij/=2) then
   msg='invalid cplex size for Dij with spin-orbit coupling !'
   MSG_BUG(msg)
 end if
 if(paw_an(iatom)%has_vhartree==0) then
   msg='vh1 must be allocated !'
   MSG_BUG(msg)
 end if
 if(paw_an(iatom)%has_vxc==0) then
   msg='vxc1 must be allocated !'
   MSG_BUG(msg)
 end if
 if(paw_ij(iatom)%has_dijso==0) then
   msg='dijso must be allocated !'
   MSG_BUG(msg)
 end if
 if (pawxcdev==0.and.pawang%angl_size==0) then
   msg='  pawxcdev=0 and pawang%angl_size=0 !'
   MSG_BUG(msg)
 end if
 if (pawang%use_ls_ylm==0) then
   msg='  pawang%use_ls_ylm=0 !'
   MSG_BUG(msg)
 end if

!------------------------------------------------------------------------
!----------- Allocations and initializations
!------------------------------------------------------------------------

 mesh_size=pawrad(itypat)%mesh_size
 lmn2_size=paw_ij(iatom)%lmn2_size
 ij_size=pawtab(itypat)%ij_size
 cplex_dij=paw_ij(iatom)%cplex_dij
 nsploop=4

 ABI_ALLOCATE(indklmn,(6,lmn2_size))
 indklmn(:,:)=pawtab(itypat)%indklmn(:,:)

!Eventually compute <Phi_i|1/r.dV/dr|Phi_j>*alpha2/2*Y_00 (for spin-orbit)
 ABI_ALLOCATE(dv1dr,(mesh_size))
 ABI_ALLOCATE(dijso_rad,(ij_size))
 ABI_ALLOCATE(ff,(mesh_size))
 fact=one/sqrt(four_pi) ! Y_00
 if (pawxcdev/=0) then
   if (paw_an(iatom)%nspden==1) then
     ff(1:mesh_size)=paw_an(iatom)%vxc1(1:mesh_size,1,1)
   else
     ff(1:mesh_size)=half*(paw_an(iatom)%vxc1(1:mesh_size,1,1) &
&     +paw_an(iatom)%vxc1(1:mesh_size,1,2))
   end if
 else
   ff(1:mesh_size)=zero
   if (paw_an(iatom)%nspden==1) then
     do ipts=1,pawang%angl_size
       ff(1:mesh_size)=ff(1:mesh_size) &
&       +paw_an(iatom)%vxc1(1:mesh_size,ipts,1) &
&       *pawang%angwgth(ipts)
     end do
   else
     do ipts=1,pawang%angl_size
       ff(1:mesh_size)=ff(1:mesh_size) &
&       +half*(paw_an(iatom)%vxc1(1:mesh_size,ipts,1) &
&       +paw_an(iatom)%vxc1(1:mesh_size,ipts,2)) &
&       *pawang%angwgth(ipts)
     end do
   end if
   ff(1:mesh_size)=sqrt(four_pi)*ff(1:mesh_size)
 end if
 ff(1:mesh_size)=fact*(ff(1:mesh_size)+paw_an(iatom)%vh1(1:mesh_size,1,1))
 call nderiv_gen(dv1dr,ff,1,pawrad(itypat))
 dv1dr(2:mesh_size)=HalfFineStruct2*(one/(one-ff(2:mesh_size)/InvFineStruct**2)) &
& *dv1dr(2:mesh_size)/pawrad(itypat)%rad(2:mesh_size)
 call deducer0(dv1dr,mesh_size,pawrad(itypat))
 do kln=1,ij_size
   ff(1:mesh_size)= dv1dr(1:mesh_size)*pawtab(itypat)%phiphj(1:mesh_size,kln)
   call simp_gen(dijso_rad(kln),ff,pawrad(itypat))
 end do
 ABI_DEALLOCATE(dv1dr)
 ABI_DEALLOCATE(ff)
 dijso_rad(:)=spnorbscl*dijso_rad(:)

!------------------------------------------------------------------------
!----- Loop over density components
!------------------------------------------------------------------------
 do idij=1,nsploop

!  ------------------------------------------------------------------------
!  ----- Computation of Dij_so
!  ------------------------------------------------------------------------
   klmn1=1
   paw_ij(iatom)%dijso(:,idij)=zero
   if (mod(idij,2)==1) then
     ispden=(1+idij)/2
     do klmn=1,lmn2_size
       if (indklmn(3,klmn)==0) then   ! il==jl
         klm=indklmn(1,klmn);kln=indklmn(2,klmn)
         ilm=indklmn(5,klmn);jlm=indklmn(6,klmn)
         fact=dijso_rad(kln);if (ilm>jlm) fact=-fact
         paw_ij(iatom)%dijso(klmn1  ,idij)=fact*pawang%ls_ylm(1,klm,ispden)
         paw_ij(iatom)%dijso(klmn1+1,idij)=fact*pawang%ls_ylm(2,klm,ispden)
       end if
       klmn1=klmn1+cplex_dij
     end do
   else if (idij==2) then
     do klmn=1,lmn2_size
       if (indklmn(3,klmn)==0) then   ! il==jl
         paw_ij(iatom)%dijso(klmn1:klmn1+1,2)=-paw_ij(iatom)%dijso(klmn1:klmn1+1,1)
       end if
       klmn1=klmn1+cplex_dij
     end do
   else if (idij==4) then
     do klmn=1,lmn2_size
       if (indklmn(3,klmn)==0) then   ! il==jl
         paw_ij(iatom)%dijso(klmn1  ,4)=-paw_ij(iatom)%dijso(klmn1  ,3)
         paw_ij(iatom)%dijso(klmn1+1,4)= paw_ij(iatom)%dijso(klmn1+1,3)
       end if
       klmn1=klmn1+cplex_dij
     end do
   end if

!  ----- End loop over idij
 end do
 ABI_DEALLOCATE(indklmn)
 ABI_DEALLOCATE(dijso_rad)
 paw_ij(iatom)%has_dijso=2

 DBG_EXIT("COLL")

end subroutine pawdijso
!!***
