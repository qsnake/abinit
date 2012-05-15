!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawxpot
!! NAME
!! pawxpot
!!
!! FUNCTION
!! Local exact-exchange potential for PAW calculations
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (BA,FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors.
!!
!! INPUTS
!!  pawprtvol=control print volume and debugging output for PAW
!!  paw_ij <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawtab <type(pawtab_type)>=paw tabulated starting data:
!!  pawrhoij <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!
!! OUTPUT
!!  paw_ij%vpawx(pawtab%lexexch*2+1,pawtab%lexexch*2+1)=local exact-exchange potential
!!
!! PARENTS
!!      pawdenpot
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawxpot(pawprtvol,pawtab,paw_ij,pawrhoij)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawxpot'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: pawprtvol
 type(paw_ij_type),intent(inout) :: paw_ij
 type(pawrhoij_type),intent(in) :: pawrhoij
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: irhoij,irhoij1,ispden,jrhoij,jrhoij1,klmn,klmn1,lexexch,ll,m11,m21,m31,m41,n1
 integer :: n2,n3,n4,nk,nn1,nn2
 real(dp) :: tot
 character(len=500) :: message
!arrays
 integer :: indn(3,3)
 real(dp) :: factnk(6)

! *****************************************************

 DBG_ENTER("COLL")

 lexexch=pawtab%lexexch
 if (pawtab%nproju==1) nk=1
 if (pawtab%nproju==2) nk=6
 factnk(1)=one;factnk(2)=one;factnk(3)=one
 factnk(4)=two;factnk(5)=two;factnk(6)=two
 indn(1,1)=1;indn(1,2)=4;indn(1,3)=5
 indn(2,1)=4;indn(2,2)=2;indn(2,3)=6
 indn(3,1)=5;indn(3,2)=6;indn(3,3)=3

!======================================================
!Compute local exact exchange Potential on the basis of projectors.
!-----------------------------------------------------

 paw_ij%vpawx=zero
 do ispden=1,paw_ij%nspden
   jrhoij=1
   do irhoij=1,pawrhoij%nrhoijsel
     klmn=pawrhoij%rhoijselect(irhoij)
     if(pawtab%indklmn(3,klmn)==0.and.pawtab%indklmn(4,klmn)==2*lexexch) then
       m11=pawtab%klmntomn(1,klmn);m21=pawtab%klmntomn(2,klmn)
       n1=pawtab%klmntomn(3,klmn);n2=pawtab%klmntomn(4,klmn)
       nn1=(n1*n2)/2+1
       jrhoij1=1
       do irhoij1=1,pawrhoij%nrhoijsel
         klmn1=pawrhoij%rhoijselect(irhoij1)
         if(pawtab%indklmn(3,klmn1)==0.and.pawtab%indklmn(4,klmn1)==2*lexexch) then
           m31=pawtab%klmntomn(1,klmn1);m41=pawtab%klmntomn(2,klmn1)
           n3=pawtab%klmntomn(3,klmn1);n4=pawtab%klmntomn(4,klmn1)
           nn2=(n3*n4)/2+1
           do ll=1,lexexch+1
             paw_ij%vpawx(1,klmn,ispden)=paw_ij%vpawx(1,klmn,ispden)&
&             -pawtab%vex(m11,m31,m41,m21,ll)*pawtab%dltij(klmn1)*pawtab%fk(indn(nn1,nn2),ll)*&
&             pawrhoij%rhoijp(jrhoij1,ispden)
           end do

         end if
         jrhoij1=jrhoij1+pawrhoij%cplex
       end do !irhoij1
     end if
     jrhoij=jrhoij+pawrhoij%cplex
   end do !irhoij
 end do !ispden

 paw_ij%has_exexch_pot=2

!test
 if (abs(pawprtvol)>=2) then
   tot=zero
   do ispden=1,pawrhoij%nspden
     jrhoij=1
     do irhoij=1,pawrhoij%nrhoijsel
       klmn=pawrhoij%rhoijselect(irhoij)
       tot=tot+paw_ij%vpawx(1,klmn,ispden)*pawrhoij%rhoijp(jrhoij,ispden)*pawtab%dltij(klmn)
       jrhoij=jrhoij+pawrhoij%cplex
     end do
   end do
   write(message, '(a,es22.15)' )" Vpawx: tot=",tot*half
   call wrtout(std_out,  message,'COLL')
 end if

 DBG_EXIT("COLL")

 end subroutine pawxpot
!!***
