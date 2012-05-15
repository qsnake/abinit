!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawdijhartree
!! NAME
!! pawdijhartree
!!
!! FUNCTION
!! Compute the Hartree contribution to the PAW
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
!!  cplex: if 1, on-site quantities are REAL, if 2, COMPLEX (response function only)
!!  iatom=index of current atom (note: this is the index on current proc, not the absolute index)
!!  natom=number of atoms in cell (note: this is the number on current proc, not the absolute number)
!!  ntypat=number of types of atoms in unit cell.
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  paw_ij(iatom)%dijhartree(cplex*lmn2_size,ndij)= Hartree Dij terms
!!        dijhartree(:,:,1) contains Dij_H^up-up
!!        dijhartree(:,:,2) contains Dij_H^dn-dn
!!        dijhartree(:,:,3) contains Dij_H^up-dn
!!        dijhartree(:,:,4) contains Dij_H^dn-up
!!
!! PARENTS
!!      pawdenpot,pawnstd2e
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawdijhartree(cplex,iatom,natom,ntypat,paw_ij,pawrhoij,pawtab)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdijhartree'
!End of the abilint section

 implicit none
!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,iatom,natom,ntypat
!arrays
 type(paw_ij_type),intent(inout) :: paw_ij(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: irhoij,ispden,itypat,jrhoij,kklmn,kklmn1,klmn,klmn1,lmn2_size,nspdiag
 character(len=500) :: msg
!arrays
 real(dp) :: ro(cplex)

! *************************************************************************

 DBG_ENTER("COLL")

 if(paw_ij(iatom)%has_dijhartree==0) then
   msg='  dijhartree must be allocated !'
   MSG_BUG(msg)
 end if
 if (pawrhoij(1)%cplex<cplex) then
   msg='  pawrhoij()%cplex must be >=cplex  !'
   MSG_BUG(msg)
 end if

!------------------------------------------------------------------------
!----------- Allocations and initializations
!------------------------------------------------------------------------

 itypat=pawrhoij(iatom)%itypat
 lmn2_size=paw_ij(iatom)%lmn2_size
 nspdiag=1;if (paw_ij(iatom)%nspden==2) nspdiag=2

 paw_ij(iatom)%dijhartree=zero

!Real on-site quantities (ground-state calculation)
 if (cplex==1) then
   do ispden=1,nspdiag
     jrhoij=1
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
       klmn=pawrhoij(iatom)%rhoijselect(irhoij)
       ro(1)=pawrhoij(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
       paw_ij(iatom)%dijhartree(klmn)=paw_ij(iatom)%dijhartree(klmn)&
&       +ro(1)*pawtab(itypat)%eijkl(klmn,klmn)
       do klmn1=1,klmn-1
         paw_ij(iatom)%dijhartree(klmn1)=paw_ij(iatom)%dijhartree(klmn1)&
&         +ro(1)*pawtab(itypat)%eijkl(klmn1,klmn)
       end do
       do klmn1=klmn+1,lmn2_size
         paw_ij(iatom)%dijhartree(klmn1)=paw_ij(iatom)%dijhartree(klmn1)&
&         +ro(1)*pawtab(itypat)%eijkl(klmn,klmn1)
       end do
       jrhoij=jrhoij+pawrhoij(iatom)%cplex
     end do
   end do

!  Complex on-site quantities (response function calculation)
 else
   do ispden=1,nspdiag
     jrhoij=1
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
       klmn=pawrhoij(iatom)%rhoijselect(irhoij);kklmn=2*klmn-1
       ro(1:2)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+1,ispden)*pawtab(itypat)%dltij(klmn)
       paw_ij(iatom)%dijhartree(kklmn:kklmn+1)=paw_ij(iatom)%dijhartree(kklmn:kklmn+1)&
&       +ro(1:2)*pawtab(itypat)%eijkl(klmn,klmn)
       do klmn1=1,klmn-1
         kklmn1=2*klmn1-1
         paw_ij(iatom)%dijhartree(kklmn1:kklmn1+1)=paw_ij(iatom)%dijhartree(kklmn1:kklmn1+1)&
&         +ro(1:2)*pawtab(itypat)%eijkl(klmn1,klmn)
       end do
       do klmn1=klmn+1,lmn2_size
         kklmn1=2*klmn1-1
         paw_ij(iatom)%dijhartree(kklmn1:kklmn1+1)=paw_ij(iatom)%dijhartree(kklmn1:kklmn1+1)&
&         +ro(1:2)*pawtab(itypat)%eijkl(klmn,klmn1)
       end do
       jrhoij=jrhoij+pawrhoij(iatom)%cplex
     end do
   end do
 end if

 paw_ij(iatom)%has_dijhartree=2

 DBG_EXIT("COLL")

end subroutine pawdijhartree
!!***
