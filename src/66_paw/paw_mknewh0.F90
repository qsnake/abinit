!!****f* ABINIT/paw_mknewh0
!! NAME
!! paw_mknewh0
!!
!! FUNCTION
!! Calculates the new bare PAW Hamiltonian in the case of quasi-particle self-consistent GW calculations.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nspden=number of spin-density components
!!  nfftf=(effective) number of FFT grid points (for this proc) for the "fine" grid
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated
!!  pawprtvol=control print volume and debugging output for PAW
!!  Cryst<Crystal_structure>=Info on unit cell and its symmetries
!!  Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  Pawtab(ntypat*usepaw)<type(pawtab_type)>=paw tabulated starting data
!!  Paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  Pawang<type(pawang_type)>=paw angular mesh and related data
!!  Pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  vxc(nfftf,nspden)=exchange-correlation potential
!!  vxc_val(nfftf,nspden)=valence only exchange-correlation potential
!!  vtrial(nfftf,nspden)=potential (Hartree+XC+loc)
!!
!! SIDE EFFECTS
!!  Paw_ij(natom*usepaw)<Paw_ij_type)>=paw arrays given on (i,j) channels
!!     At output: new value for Paw_ij()%dij
!!
!! PARENTS
!!      calc_vhxc_me
!!
!! CHILDREN
!!      pawgylm,symdij,symdij_all,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine paw_mknewh0(nsppol,nspden,nfftf,pawspnorb,pawprtvol,Cryst,Psps,&
&          Pawtab,Paw_an,Paw_ij,Pawang,Pawfgrtab,vxc,vxc_val,vtrial)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling
 use m_errors

 use m_crystal,only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_mknewh0'
 use interfaces_14_hidewrite
 use interfaces_66_paw, except_this_one => paw_mknewh0
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsppol,nspden,nfftf,pawprtvol,pawspnorb
!arrays
 real(dp),intent(in) :: vxc(nfftf,nspden),vxc_val(nfftf,nspden),vtrial(nfftf,nspden)
 type(Crystal_structure),intent(in) :: Cryst
 type(Pseudopotential_type),intent(in) :: Psps
 type(Pawang_type),intent(in) :: Pawang
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat)
 type(Paw_an_type),intent(in) :: Paw_an(Cryst%natom)
 type(Paw_ij_type),intent(inout) :: Paw_ij(Cryst%natom)
 type(Pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: ipert0=0
 integer :: iat,idij,cplex,ndij,option_dij
 integer :: itypat,lmn_size,j0lmn,jlmn,ilmn,klmn,klmn1,klm
 integer :: lmin,lmax,mm,isel,lm_size,lmn2_size,cplex_dij
 integer :: ils,ilslm,ic,lm0
 integer :: nsploop,is2fft
 real(dp) :: gylm,qijl
 logical :: ltest
 character(len=500) :: msg
!arrays
 integer,pointer :: indklmn_(:,:)
 real(dp) :: rdum(1)
 real(dp),allocatable :: prod_hloc(:,:),prodhxc_core(:,:)
 real(dp),allocatable :: dijhl_hat(:,:),dijhmxc_val(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call wrtout(std_out,'Assembling PAW strengths for the bare Hamiltonian','COLL')

!=== Test if required pointers in paw_ij are associated ===
 ltest = (associated(Paw_ij(1)%dijxc).and.associated(Paw_ij(1)%dijxc_val) )
!&        .and.(Paw_ij(1)%has_dijxc==2.and.paw_ij(1)%has_dijxc_val==2))
 ABI_CHECK(ltest,'dijxc or dijxc_val not calculated')

 ltest=(associated(Paw_ij(1)%dijhat)) !.and.Paw_ij(1)%has_dijhat==2)
 ABI_CHECK(ltest,'dijhat not calculated')

 ltest=(associated(Paw_ij(1)%dijhartree)) !.and.Paw_ij(1)%has_dijhartree==2)
 ABI_CHECK(ltest,'dijhartree not calculated')

 if (ANY(Pawtab(:)%usepawu>0)) then
   do iat=1,Cryst%natom
     itypat=Cryst%typat(iat)
     if (Pawtab(itypat)%usepawu>0) then
       ltest=(associated(Paw_ij(iat)%dijU) ) !.and.Paw_ij(iat)%has_dijU==2)
       write(msg,'(a,i3,a)')" For atom no. ",iat," %dijU(iat) has not been calculated."
       ABI_CHECK(ltest,msg)
     end if
   end do
 end if

 if (pawspnorb>0) then
   do iat=1,Cryst%natom
     ltest=(associated(Paw_ij(iat)%dijso) ) !.and.Paw_ij(iat)%has_dijso==2)
     write(msg,'(a,i3,a)')" For atom no. ",iat," %dijso(iat) has not been calculated."
     ABI_CHECK(ltest,msg)
   end do
 end if

!== Construct the new PAW H0 Hamiltonian ===
 do iat=1,Cryst%natom

   itypat    = Cryst%typat(iat)
   lmn_size  = Pawtab(itypat)%lmn_size
   lmn2_size = Pawtab(itypat)%lmn2_size
   lm_size   = Paw_an(iat)%lm_size
   cplex     = Paw_ij(iat)%cplex
   cplex_dij = Paw_ij(iat)%cplex_dij
   ndij      = Paw_ij(iat)%ndij

   ABI_CHECK(cplex==1,'cplex/=1 not implemented')
   ABI_CHECK(cplex_dij==1,'cplex_dij/=1 not implemented')
!  
!  Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
   if (Pawfgrtab(iat)%gylm_allocated==0) then
     if (associated(Pawfgrtab(iat)%gylm))  then
       ABI_DEALLOCATE(Pawfgrtab(iat)%gylm)
     end if
     ABI_ALLOCATE(Pawfgrtab(iat)%gylm,(Pawfgrtab(iat)%nfgd,lm_size))
     Pawfgrtab(iat)%gylm_allocated=2

     call pawgylm(Pawfgrtab(iat)%gylm,rdum,rdum,lm_size,&
&     Pawfgrtab(iat)%nfgd,1,0,0,Pawtab(itypat),&
&     Pawfgrtab(iat)%rfgd,Pawfgrtab(iat)%rfgd_allocated)
   end if

!  === Calculate LM contribution to dijhmxc_val for this atom ===
!  * Dijxc contains also the Hat term on the FFT mesh while Dijxc_val does not
!  contain neither the hat term nor the LM sum of onsite terms (they should cancel each other)
!  FIXME change paw_dij,  otherwise I miss tnc in vxc
!  * prodhxc_core is used to assemble $\int g_l Ylm (vtrial - vxc_val[tn+nhat] dr$ on the FFT mesh ===
!  * The following quantities do not depend on ij
   ABI_ALLOCATE(prod_hloc   ,(lm_size,ndij))
   ABI_ALLOCATE(prodhxc_core,(lm_size,ndij))
   prod_hloc   =zero
   prodhxc_core=zero
   do idij=1,ndij
     do ilslm=1,lm_size
       do ic=1,Pawfgrtab(iat)%nfgd
         is2fft=Pawfgrtab(iat)%ifftsph(ic)
         gylm=Pawfgrtab(iat)%gylm(ic,ilslm)
         prod_hloc (ilslm,idij)=prod_hloc (ilslm,idij) + (vtrial(is2fft,idij)-vxc(is2fft,idij))*gylm
!        prodhxc_core(ilslm,idij)=prodhxc_core(ilslm,idij) + (vxc_val(is2fft,idij))*gylm
         prodhxc_core(ilslm,idij)=prodhxc_core(ilslm,idij) + (vtrial(is2fft,idij)-vxc_val(is2fft,idij))*gylm
       end do
     end do
   end do !idij

!  === Assembly the "Hat" contribution for this atom ====
   ABI_ALLOCATE(dijhl_hat  ,(cplex_dij*lmn2_size,ndij))
   ABI_ALLOCATE(dijhmxc_val,(cplex_dij*lmn2_size,ndij))
   dijhl_hat  =zero
   dijhmxc_val=zero
   indklmn_ => Pawtab(itypat)%indklmn(1:6,1:lmn2_size)

   do idij=1,ndij
     do klmn=1,lmn2_size
       klm =indklmn_(1,klmn)
       lmin=indklmn_(3,klmn)
       lmax=indklmn_(4,klmn)

!      === $\sum_lm q_ij^l prod* for each idij$ ===
       do ils=lmin,lmax,2
         lm0=ils**2+ils+1
         do mm=-ils,ils
           ilslm=lm0+mm
           isel=Pawang%gntselect(lm0+mm,klm)
           if (isel>0) then
             qijl=Pawtab(itypat)%qijl(ilslm,klmn)
             dijhl_hat  (klmn,idij)=dijhl_hat  (klmn,idij) +  prod_hloc (ilslm,idij)*qijl
             dijhmxc_val(klmn,idij)=dijhmxc_val(klmn,idij) +prodhxc_core(ilslm,idij)*qijl
           end if
         end do
       end do
     end do
   end do

   ABI_DEALLOCATE(prod_hloc)
   ABI_DEALLOCATE(prodhxc_core)

!  * Normalization factor due to integration on the FFT mesh
   dijhl_hat  = dijhl_hat  *Cryst%ucvol/DBLE(nfftf)
   dijhmxc_val= dijhmxc_val*Cryst%ucvol/DBLE(nfftf)

!  === Now assembly the bare Hamiltonian ===
!  * Loop over density components overwriting %dij
   nsploop=nsppol; if (Paw_ij(iat)%ndij==4) nsploop=4

   do idij=1,nsploop
     klmn1=1

     do jlmn=1,lmn_size
       j0lmn=jlmn*(jlmn-1)/2
       do ilmn=1,jlmn
         klmn=j0lmn+ilmn

!        The following gives back the input dij.
!        since dijxc contains the hat term done on the FFT mesh
         if (.FALSE.) then
           Paw_ij(iat)%dij(klmn,idij) =        &
&           Pawtab(itypat)%dij0    (klmn)      &
&           +Paw_ij(iat)%dijhartree(klmn)      &
&           +Paw_ij(iat)%dijxc     (klmn,idij) &
&           +dijhl_hat   (klmn,idij)

         else
!          === Make nonlocal part of h0 removing the valence contribution ===
!          Remeber that XC contains already the Hat contribution
           Paw_ij(iat)%dij(klmn,idij) =        &
&           Pawtab(itypat)%dij0      (klmn)    &
&           +Paw_ij(iat)%dijhartree(klmn)      &
&           +Paw_ij(iat)%dijxc     (klmn,idij) &  ! 2 lines to get the d1-dt1 XC core contribution + XC hat (core+val)
&          -Paw_ij(iat)%dijxc_val (klmn,idij) &  ! I suppose that the "hat" term on the FFT mesh in included in both.
&          +dijhmxc_val(klmn,idij)               ! Local + Hartree - XC val contribution to the "hat" term.

!          Add the U contribution to the 
!          if (.FALSE. .and. Pawtab(itypat)%usepawu>0) then
           if (.TRUE. .and. Pawtab(itypat)%usepawu>0) then
             Paw_ij(iat)%dij(klmn,idij) = Paw_ij(iat)%dij(klmn,idij) + Paw_ij(iat)%dijU(klmn,idij)
           end if
         end if
!        TODO dijso, dijU, vpawx?
!        Just to be consistent, update some values.
!        $Paw_ij(iat)%dijhat(klmn,idij)=Paw_ij(iat)%dijhat(klmn,idij)-dijhmxc_val(klmn,idij)

       end do !ilmn
     end do !jlmn
   end do !idij

!  this is to be consistent?
!  deallocate(Paw_ij(iat)%dijvxc_val)
   ABI_DEALLOCATE(dijhl_hat)
   ABI_DEALLOCATE(dijhmxc_val)
 end do !iat

!=== Symmetrize total Dij ===
 option_dij=0 ! For total Dij.
#if 0
 call symdij(Cryst%gprimd,Psps%indlmn,Cryst%indsym,ipert0,Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,option_dij,&
& Paw_ij,Pawang,pawprtvol,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)
#else
 call symdij_all(Cryst%gprimd,Psps%indlmn,Cryst%indsym,ipert0,Psps%lmnmax,Cryst%natom,Cryst%nsym,Cryst%ntypat,&
& Paw_ij,Pawang,pawprtvol,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)
#endif

 DBG_EXIT("COLL")

end subroutine paw_mknewh0
!!***
