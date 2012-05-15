!{\src2tex{textfont=tt}}
!!****f* ABINIT/symmetrize_afm_chi0
!! NAME
!! symmetrize_afm_chi0
!!
!! FUNCTION
!!  Reconstruct the (down, down) component of the irreducible polarizability
!!  starting from the (up,up) element in case of systems with AFM symmetries 
!!  (i.e nspden==2 and nsppol=1). Return the trace (up,up)+(down,down) of the 
!!  matrix as required by GW calculations.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Cryst<Crystal_structure>= Information on symmetries and unit cell.
!!  Gsph<Gvectors_type>= The G-sphere used to descrive chi0.
!!  npwe=Number of G-vectors in chi0.
!!  nomega=number of frequencies.
!!  Ltg_q<Little_group>=Structure with useful table describing the little group of the q-point.
!!
!! SIDE EFFECTS
!! chi0(npwe,npwe,nomega)= In input the up-up component, in output the trace of chi0. 
!!   The value of matrix elements that should be zero due to AFM symmetry properties are 
!!   forced to be zero (see NOTES below).
!! [chi0_lwing(npwe,nomega,3)] = Lower wings, symmetrized in output.
!! [chi0_uwing(npwe,nomega,3)] = Upper wings, symmetrized in output.
!! [chi0_head(3,3,nomega)    ] = Head of chi0, symmetrized  in output.
!!
!! NOTES
!! In the case of magnetic group Shubnikov type III:
!!   For each set of paired FM-AFM symmetries, the down-down component of
!!   a generic response function in reciprocal space can be obtained according to:
!!
!!    chi^{down,down}_{G1,G2}(q) = chi^{up up}_{G1,G2}(q) e^{iS(G1-G2).(tnonsFM - tnonsAFM)}
!!
!!   where S is the rotational part common to the FM-AFM pair, tnonsFM and tnonsAFM 
!!   are the fractional translations associated to the ferromagnetic and antiferromagnetic symmetry, respectively.
!!   Note that, if for a given G1-G2 pair, the phase e^{iS(G1-G2).(tnonsFM - tnonsAFM) depends 
!!   on the FM-AFM symmetry pair, then the corresponding matrix element of chi0 must be zero.
!!   Actually this is manually enforced in the code because this property might not be 
!!   perfectly satisfied due to round-off errors.
!!
!! In the case of magnetic group Shubnikov type III:
!!   Only the AFM symmetries that preserve the external q-point (with or without time-reversal)
!!   are used to get the (down, down) component using the fact that:
!!
!!    chi^{down,down}_{G1,G2}(Sq) = chi^{up up}_{S^{-1}G1,S^{-1}G2}(q) e^{i(G2-G1).tnons_S }
!!
!!   Actually we perform an average over subset of the little group of q with AFM character in 
!!   order to reduce as much as possible errors due to round off errors. In brief we evaluate:
!!
!!    1/N_{Ltq} \sum_{S\in Ltg AFM} chi^{up up}_{S^{-1}G1,S^{-1}G2}(q) e^{i(G2-G1).tnons_S }
!!   
!!   where N_{Ltg} is the number of AFM operation in the little group (time reversal included)
!!
!! TODO 
!!  It is possible to symmetrize chi0 without any the extra allocation for afm_mat.
!!  More CPU demanding but safer in case of a large chi0 matrix. One might loop over G1 and G2 shells ...
!!
!! PARENTS
!!      cchi0,cchi0q0,cchi0q0_intraband,check_completeness
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine symmetrize_afm_chi0(Cryst,Gsph,Ltg_q,npwe,nomega,chi0,chi0_head,chi0_lwing,chi0_uwing)

 use m_profiling

 use defs_basis
 use m_errors

 use m_gwdefs,   only : czero_gw, cone_gw
 use m_crystal,  only : crystal_structure
 use m_gsphere,  only : gvectors_type
 use m_bz_mesh,  only : little_group

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symmetrize_afm_chi0'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwe,nomega
 type(Gvectors_type),intent(in) :: Gsph
 type(Crystal_structure),intent(in) :: Cryst
 type(Little_group),intent(in) :: Ltg_q
!arrays
 complex(gwpc),intent(inout) :: chi0(npwe,npwe,nomega)
 complex(dpc),optional,intent(inout) :: chi0_lwing(npwe,nomega,3)
 complex(dpc),optional,intent(inout) :: chi0_uwing(npwe,nomega,3)
 complex(dpc),optional,intent(inout) :: chi0_head(3,3,nomega)

!Local variables ------------------------------
!scalars
 integer :: io,ig1,ig2,isymf,isyma,isym,ipair,k0g,kg,npairs,nonzero
 integer :: iSmg1,iSmg2,itim,shubnikov,ntest
 complex(gwpc) :: phase,phase_old,sumchi,ctmp
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: rotfm(3,3),rotafm(3,3),pairs2sym(2,Cryst%nsym/2)
 real(dp) :: tfm(3),tafm(3)
 complex(gwpc),allocatable :: afm_mat(:),chi0_afm(:,:)

!************************************************************************

 ABI_CHECK(ANY(Cryst%symafm==-1),'Not magnetic space group')
 !
 ! ==== Find shubnikov type ====
 ! TODO This info should be stored in Cryst%
 shubnikov=4; npairs=0

 do isymf=1,Cryst%nsym
   if (Cryst%symafm(isymf)==-1) CYCLE
   rotfm = Cryst%symrec(:,:,isymf)
   tfm = Cryst%tnons(:,isymf)
   found = .FALSE.

   do isyma=1,Cryst%nsym
     if (Cryst%symafm(isyma)==1) CYCLE
     rotafm = Cryst%symrec(:,:,isyma)

     if (ALL(rotfm==rotafm)) then
       found=.TRUE.
       tafm = Cryst%tnons(:,isyma)
       npairs=npairs+1
       ABI_CHECK(npairs<=Cryst%nsym/2,'Wrong AFM group')
       pairs2sym(1,npairs)=isymf
       pairs2sym(2,npairs)=isyma
     end if
   end do !isyma

   if (.not.found) then 
    shubnikov=3; EXIT !isymf
   end if
 end do !isymf

 select case (shubnikov)

 case (4)
   call wrtout(std_out,' Found Magnetic group Shubnikov type IV','COLL')
   ABI_CHECK(npairs==Cryst%nsym/2,'Wrong AFM space group') 

   ABI_ALLOCATE(afm_mat,(npwe*(npwe+1)/2))

   do ig2=1,npwe
     k0g=ig2*(ig2-1)/2 
     do ig1=1,ig2
       kg=k0g+ig1  
       nonzero=1

       do ipair=1,Cryst%nsym/2
         isymf = pairs2sym(1,ipair)
         isyma = pairs2sym(2,ipair)
         phase = ( Gsph%phmSGt(ig1,isymf)*CONJG(Gsph%phmSGt(ig1,isyma)) ) * &
&                ( Gsph%phmSGt(ig2,isymf)*CONJG(Gsph%phmSGt(ig2,isyma)) )
         if (ipair>1 .and. (ABS(phase_old-phase) > tol6)) then 
           nonzero=0; EXIT
         end if
         phase_old=phase
       end do !ipair

       afm_mat(kg)=nonzero*(cone_gw + phase)
     end do !ig1
   end do !ig2
   ! 
   ! =======================================================================
   ! ==== Symmetrize chi0 constructing chi0^{\up\up} + chi0{\down\down} ====
   ! =======================================================================
   !
   !  head^{\down\down} = head^{\up\up}
   if (PRESENT(chi0_head)) chi0_head = two * chi0_head 
   !
   ! w^{\down\down}_{0 G'} =  e^{-iSG'.(tFM-tAFM)} w^{\up\up}_{0 G'}.
   ! w^{\down\down}_{G 0 } =  e^{+iSG .(tFM-tAFM)} w^{\up\up}_{G 0 }. 
   if (PRESENT(chi0_uwing)) then
     do io=1,nomega
       do ig2=1,npwe
         k0g=ig2*(ig2-1)/2 
         kg=k0g+1  
         chi0_uwing(ig2,io,:)=afm_mat(kg)*chi0_uwing(ig2,io,:)
       end do
     end do
   end if
   !                                                               
   if (PRESENT(chi0_lwing)) then
     do io=1,nomega
       do ig1=1,npwe
         k0g=ig1*(ig1-1)/2 
         kg=k0g+1  
         chi0_lwing(ig1,io,:)=CONJG(afm_mat(kg))*chi0_lwing(ig1,io,:)
       end do
     end do
   end if

   do io=1,nomega
     !
     do ig1=1,npwe ! Take care of diagonal.
       chi0(ig1,ig1,io)=two*chi0(ig1,ig1,io)
     end do
     !
     ! Upper and lower triangle are treated differently:
     ! We took advantage of the fact the afm_mat is hermitian to reduce memory.
     do ig2=2,npwe
       k0g=ig2*(ig2-1)/2 
       do ig1=1,ig2-1
         kg=k0g+ig1  
         chi0(ig1,ig2,io)=afm_mat(kg)*chi0(ig1,ig2,io)
       end do
     end do
     !
     do ig1=2,npwe
       k0g=ig1*(ig1-1)/2 
       do ig2=1,ig1-1
         kg=k0g+ig2 
         chi0(ig1,ig2,io)=CONJG(afm_mat(kg))*chi0(ig1,ig2,io)
       end do
     end do
     !
   end do !io 

   ABI_DEALLOCATE(afm_mat)

 case (3)
   call wrtout(std_out,' Found Magnetic group Shubnikov type III',"COLL")
   MSG_ERROR('Shubnikov type III not implemented')

   ntest=0
   do itim=1,ltg_q%timrev 
     do isym=1,ltg_q%nsym_sg
       ! use only afm sym preserving q with and without time-reversal
       if ( cryst%symafm(isym)==-1 .and. ltg_q%preserve(itim,isym)==1 ) ntest=ntest+1
     end do
   end do

   if (ntest==0) write(std_out,*)' WARNING: no symmetry can be used!!! '
   !RETURN
   ABI_ALLOCATE(chi0_afm,(npwe,npwe))

   do io=1,nomega

     do ig2=1,npwe
       do ig1=1,npwe
         sumchi=czero_gw

         do itim=1,ltg_q%timrev 
           do isym=1,ltg_q%nsym_sg
             ! use only afm sym preserving q with and without time-reversal
             if ( cryst%symafm(isym)==-1 .and. ltg_q%preserve(itim,isym)==1 ) then
               phase =  Gsph%phmGt(ig1,isym)*CONJG(Gsph%phmGt(ig2,isym)) 
               iSmg1=Gsph%rottbm1(ig1,itim,isym)
               iSmg2=Gsph%rottbm1(ig2,itim,isym)
               ctmp=chi0(iSmg1,iSmg2,io)*phase !; if (itim==2) ctmp=CONJG(ctmp) !check this
               sumchi=sumchi+ctmp !chi0(iSmg1,iSmg2,io)*phase
             end if
           end do ! isym
         end do !itim

         chi0_afm(ig1,ig2)=sumchi/Ltg_q%nsym_ltg  !has to be changed in case of time-reversal
       end do !ig1
     end do !ig2

     ! We want chi_{up,up} +chi_{dwn,dwn}.
     chi0(:,:,io)=chi0(:,:,io)+chi0_afm(:,:)
   end do !iomega

   ABI_DEALLOCATE(chi0_afm)

 case default
   write(msg,'(a,i4)')'Wrong value for shubnikov= ',shubnikov
   MSG_BUG(msg)
 end select

end subroutine symmetrize_afm_chi0
!!***
