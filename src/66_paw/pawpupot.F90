!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawpupot
!! NAME
!! pawpupot
!!
!! FUNCTION
!! Compute LDA+U potential
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (BA,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors.
!!
!! INPUTS
!!  ispden=current spin component
!!  pawprtvol=control print volume and debugging output for PAW
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!    %noccmmp(2*lpawu+1,2*lpawu+1,nspden)=density matrix in the augm. region
!!    %nocctot(nspden)=number of electrons in the correlated subspace
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data:
!!     %lpawu=l used for lda+u
!!     %vee(2*lpawu+1*4)=screened coulomb matrix
!!
!! OUTPUT
!!  vpawu(cplex_dij,lpawu*2+1,lpawu*2+1)=lda+u potential (see eg PRB 52, 5467 (1995))
!!
!! SIDE EFFECTS
!!  VUKS=Sum over spins and atoms of nocc_mmp*Pot_(LDA+U)
!!
!! PARENTS
!!      ldau_self,m_paw_commutator,pawdij
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawpupot(ispden,paw_ij,pawprtvol,pawtab,vpawu,VUKS)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawpupot'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ispden,pawprtvol
 real(dp),intent(inout) :: VUKS
 type(paw_ij_type),intent(in) :: paw_ij
 type(pawtab_type),intent(in) :: pawtab
!arrays
 real(dp),intent(inout) :: vpawu(paw_ij%cplex_dij,pawtab%lpawu*2+1,pawtab%lpawu*2+1)

!Local variables ---------------------------------------
!scalars
!Option for interaction energy in case of non-collinear magnetism:
!           1: E_int=-J/4.N.(N-2)                   (better)
!           2: E_int=-J/2.(Nup.(Nup-1)+Ndn.(Ndn-1)) (Nup and Ndn are ill-defined)
 integer,parameter :: option_interaction=1

 integer :: iplex,jspden,lpawu,m1,m11,m2,m21,m3,m31,m4,m41,cplex_dij
 real(dp) :: mnorm,mx,my,mz,n_sig,n_msig,n_tot,VUKStemp,n_sigs,n_msigs
 real(dp),allocatable :: n34_msig(:),n34_sig(:),n43_sig(:)
 character(len=500) :: message
!array
 real(dp),parameter :: factcg(3:4)=(/one,-one/)

! *****************************************************

 DBG_ENTER("COLL")

 lpawu=pawtab%lpawu
 vpawu=zero
 cplex_dij=paw_ij%cplex_dij
 ABI_ALLOCATE(n34_sig,(cplex_dij))
 ABI_ALLOCATE(n34_msig,(cplex_dij))
 ABI_ALLOCATE(n43_sig,(cplex_dij))

!======================================================
!Compute LDA+U Potential on the basis of projectors.
!cf PRB 52 5467 (1995)
!-----------------------------------------------------

 if (ispden<=2) then   ! cases ndij=4, ispden=1,2 or ndij<4
!  jspden=min(paw_ij%nspden,2)-ispden+1
   jspden=min(paw_ij%ndij,2)-ispden+1   ! (ispden,ndij)=(1,4)=>jspden=2

   if (paw_ij%nspden<=2) then
     n_sig =paw_ij%nocctot(ispden)
     n_msig=paw_ij%nocctot(jspden)
     n_tot =n_sig+n_msig
   else
     n_tot=paw_ij%nocctot(1)
     mx=paw_ij%nocctot(2)
     my=paw_ij%nocctot(3)
     mz=paw_ij%nocctot(4)
     mnorm=sqrt(mx*mx+my*my+mz*mz)
     if (ispden==1) then
!      n_sig =half*(n_tot+mnorm)
!      n_msig=half*(n_tot-mnorm)
       n_sig =half*(n_tot+sign(mnorm,mz))
       n_msig=half*(n_tot-sign(mnorm,mz))
     else 
!      n_sig =half*(n_tot-mnorm)
!      n_msig=half*(n_tot+mnorm)
       n_sig =half*(n_tot-sign(mnorm,mz))
       n_msig=half*(n_tot+sign(mnorm,mz))
     end if
   end if

   n_sigs =n_sig/(float(2*lpawu+1))
   n_msigs =n_msig/(float(2*lpawu+1))
   do m1=-lpawu,lpawu
     m11=m1+lpawu+1
     do m2=-lpawu,lpawu
       m21=m2+lpawu+1
       do m3=-lpawu,lpawu
         m31=m3+lpawu+1
         do m4=-lpawu,lpawu
           m41=m4+lpawu+1
           n34_sig(:) =paw_ij%noccmmp(:,m31,m41,ispden) ! spin sigma
           n34_msig(:)=paw_ij%noccmmp(:,m31,m41,jspden) ! opposite spin (-sigma)
           if(m31==m41.and.pawtab%usepawu==3) then
             n34_sig(1)= n34_sig(1) - n_sigs
             n34_msig(1)= n34_msig(1) - n_msigs
           end if
           do iplex=1,cplex_dij
             vpawu(iplex,m11,m21)=vpawu(iplex,m11,m21) &
&             +n34_msig(iplex)*pawtab%vee(m11,m31,m21,m41) &
&             +n34_sig(iplex)*(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
           end do
           if(abs(pawprtvol)>=3.and.m11==1.and.m21==1) then
             write(message,'(a,i4,i4,2e20.10)') "m31,m41,vu=",m31,m41, vpawu(:,m11,m21)
             call wrtout(std_out,message,'COLL')
             write(message,'(a,4e20.10)') "vee",pawtab%vee(m11,m31,m21,m41),pawtab%vee(m11,m31,m41,m21)
             call wrtout(std_out,message,'COLL')
             write(message,'(a,4e20.10)') "n34_msig,n34_sig",n34_msig(1),n34_sig(1)
             call wrtout(std_out,message,'COLL')
           end if
         end do
       end do
       if(abs(pawprtvol)>=3) then
         if(m11/=m21) then
           write(message,'(a,i4,i4,2e20.10)') "vu=",m11,m21,vpawu(:,m11,m21)
           call wrtout(std_out,message,'COLL')
           write(message,'(a,2e20.10)') "vupred=",-pawtab%upawu*paw_ij%noccmmp(:,m21,m11,ispden)
           call wrtout(std_out,message,'COLL')
         end if
       end if
     end do ! m2

!    Full localized limit
     if(pawtab%usepawu==1) then ! not activated if usepawu=10 !!
       vpawu(1,m11,m11)=vpawu(1,m11,m11)-pawtab%upawu*(n_tot-half)
       if (paw_ij%ndij/=4.or.option_interaction==2) then
         vpawu(1,m11,m11)=vpawu(1,m11,m11)+pawtab%jpawu*(n_sig-half)
       else
         vpawu(1,m11,m11)=vpawu(1,m11,m11)+half*pawtab%jpawu*(n_tot-one)
       end if

!      Around mean field
     else if(pawtab%usepawu==2) then
       vpawu(1,m11,m11)=vpawu(1,m11,m11)-n_msig*pawtab%upawu &
&       -n_sig*(pawtab%upawu-pawtab%jpawu) &
&       *(dble(2*lpawu)/dble(2*lpawu+1))
     end if

     if (abs(pawprtvol)>=3) then
       write(message,'(a,i4,i4,2x,e20.10)') "vudiag= ",m11,m11,vpawu(1,m11,m11)
       call wrtout(std_out,  message,'COLL')
       write(message,'(a,2e20.10)') "vudiagpred= ",pawtab%upawu*(half-paw_ij%noccmmp(:,m11,m11,ispden))
       call wrtout(std_out,  message,'COLL')
     end if
   end do ! m1

 end if ! ispden<=2

!Non-collinear magnetism: add non-diagonal term; see (Eq 6) in PRB 72, 024458 (2005)
!BA Here, we compute the transpose --- with respect to spin indices --- of
!BA equation (6) of this reference, because of differences in notations,
!BA namely Eband=\sum rhoij^{alpha,beta}*Dij(beta,alpha) contrary to PRB 72, 024458 (2005)
 if (ispden>=3) then
   do m1=-lpawu,lpawu
     m11=m1+lpawu+1
     do m2=-lpawu,lpawu
       m21=m2+lpawu+1
       do m3=-lpawu,lpawu
         m31=m3+lpawu+1
         do m4=-lpawu,lpawu
           m41=m4+lpawu+1
!          n43_sig(:) =paw_ij%noccmmp(:,m41,m31,ispden)
!          vpawu(1,m11,m21)=vpawu(1,m11,m21)-n43_sig(1)*pawtab%vee(m11,m31,m41,m21)
!          vpawu(2,m11,m21)=vpawu(2,m11,m21)+n43_sig(2)*pawtab%vee(m11,m31,m41,m21)
           n34_sig(:) =paw_ij%noccmmp(:,m31,m41,ispden)
           vpawu(1,m11,m21)=vpawu(1,m11,m21)-n34_sig(1)*pawtab%vee(m11,m31,m41,m21)
           vpawu(2,m11,m21)=vpawu(2,m11,m21)-n34_sig(2)*pawtab%vee(m11,m31,m41,m21)
         end do
       end do
     end do
   end do
 end if
 if(abs(pawprtvol)>=3) then
   write(std_out,*) "vpawu, ispden",ispden
   do m11=1,2*lpawu+1
     write(message,'(12(1x,9(1x,"(",f10.7,",",f10.7,")")))') (vpawu(1:cplex_dij,m11,m21),m21=1,2*lpawu+1)
     call wrtout(std_out,message,'COLL')
   end do
 end if

!Printing for test
 if (abs(pawprtvol)>=3) then
   if (ispden==1) VUKS=zero
   VUKStemp=zero
   do m1=-lpawu,lpawu
     m11=m1+lpawu+1
     do m2=-lpawu,lpawu
       m21=m2+lpawu+1
       VUKStemp=VUKStemp+vpawu(1,m11,m21)*paw_ij%noccmmp(1,m11,m21,ispden)
       if (cplex_dij == 2) then
         VUKStemp=VUKStemp-vpawu(2,m11,m21)*paw_ij%noccmmp(2,m11,m21,ispden)

       end if
       write(message,'(a,2e20.10,2e20.10)') "m1,m2,vpawu,noccmmp= ", &
&       vpawu(:,m11,m21),paw_ij%noccmmp(:,m11,m21,ispden)
       call wrtout(std_out,  message,'COLL')
     end do
   end do
   VUKS=VUKS+VUKStemp
   write(message,*) "pawpupot: VUKStemp= ",ispden,VUKStemp
   call wrtout(std_out,  message,'COLL')
   if (ispden==paw_ij%nspden) then
     write(message,*) "pawpupot: VUKS= ",ispden,VUKS
     call wrtout(std_out,  message,'COLL')
   end if
 end if
 ABI_DEALLOCATE(n34_sig)
 ABI_DEALLOCATE(n34_msig)
 ABI_DEALLOCATE(n43_sig)

 DBG_EXIT("COLL")

 end subroutine pawpupot
!!***
