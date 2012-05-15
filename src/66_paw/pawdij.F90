!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawdij
!! NAME
!! pawdij
!!
!! FUNCTION
!! Compute the pseudopotential strengths Dij of the PAW non local operator.
!! Also compute several contributions to Dij.
!! Can compute first-order strenghts Dij for RF calculations.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex=(RF calculations only) - 1 if RF 1st-order quantities are REAL, 2 if COMPLEX
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   ! atvshift(16,nsppol,natom)=potential energy shift for specific lm channel, spin and atom
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  enunit=choice for units of output Dij
!!  fatvshift=factor that multiplies atvshift
!!  ipert=index of perturbation (used only for RF calculation ; set ipert<=0 for GS calculations.
!!  natom=number of atoms on current process, size of PAW arrays
!!  natom_tot=total number of atoms in cell
!!  nfft=total number of FFt grid
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in unit cell.
!!  paral_kgb=Flag related to the kpoint-band-fft parallelism
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  qphon(3)=wavevector of the phonon
!!  typat(natom)=type (integer) for each atom
!!  ucvol=unit cell volume
!!  vtrial(cplex*nfft,nspden)=GS potential
!!  vxc(cplex*nfft,nspden)=XC potential (Hartree) on the fine FFT mesh
!!  xred(3,natom)= reduced atomic coordinates
!!
!! OUTPUT
!!  paw_ij(iatom)%dij(cplex_dij*lmn2_size,ndij)= total Dij terms (GS calculation, ipert=0)
!!                                               total 1st-order Dij terms (RF ccalc., ipert>0)
!!  May be complex if cplex_dij=2
!!        dij(:,:,1) contains Dij^up-up
!!        dij(:,:,2) contains Dij^dn-dn
!!        dij(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!        dij(:,:,4) contains Dij^dn-up (only if nspinor=2)
!!  May also compute paw_ij(iatom)%dijU, paw_ij(iatom)%dijso,
!!                   paw_ij(iatom)%dijxc, paw_ij(iatom)%dijxc_val,
!!                   paw_ij(iatom)%dijhartree
!!
!! NOTES
!!  Response function calculations:
!!    In order to compute first-order Dij, paw_an (resp. paw_ij) datastructures
!!    must contain first-order quantities, namely paw_an1 (resp. paw_ij1).
!!
!! PARENTS
!!      bethe_salpeter,respfn,scfcv,scfcv3,screening,sigma
!!
!! CHILDREN
!!      leave_new,pawdijso,pawexpiqr,pawgylm,pawpupot,print_ij,simp_gen,timab
!!      wrtout,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawdij(cplex,dtset,enunit,fatvshift,gprimd,ipert,mpi_enreg,natom,natom_tot,nfft,ngfft,&
&                 nspden,ntypat,paral_kgb,paw_an,paw_ij,pawang,pawfgrtab,pawprtvol,pawrad,&
&                 pawspnorb,pawtab,pawxcdev,qphon,typat,ucvol,vtrial,vxc,xred,&
&                 electronpositron) ! optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors

 use m_radmesh,          only : simp_gen
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdij'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_66_paw, except_this_one => pawdij
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,enunit,ipert,natom,natom_tot,nfft,nspden,ntypat,paral_kgb,pawprtvol,pawspnorb,pawxcdev
 real(dp),intent(in) :: fatvshift,ucvol
 type(electronpositron_type),pointer,optional :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: gprimd(3,3),qphon(3),vtrial(cplex*nfft,nspden)
 real(dp),intent(in) :: vxc(cplex*nfft,nspden),xred(3,natom_tot)
 type(paw_an_type),intent(in) :: paw_an(natom)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)
 type(dataset_type),intent(in) :: dtset

!Local variables ---------------------------------------
!scalars
 integer :: bufdim,cplex_dij
 integer :: iatom,iatom_tot,ic,icinf,icount,icsup,idij,idijeff,idijend,ier,ij_size
 integer :: ilm,iln,ils,ils1,ilslm,ilslm1,im1,im2,in1,in2,ipositron,ipts,ir,ir1
 integer :: isel,ispden,itypat,j0lm,j0ln,jc,jlm,jln,klm,klm1,klmn,klmn1
 integer :: kln,l_size,lexexch,ll,lm0,lm_size,lmax,lmin
 integer :: lmn2_size,lmn_size,lpawu,mesh_size,mm
 integer :: need_dijhat,need_dijso,need_dijU,need_dijxc_val,need_dijxc
 integer :: nfftot,nfgd,npts,nsploop,nsppol,old_paral_level,optgr0,optgr1,shift,spaceComm,usexcnhat
 integer :: use_dijhat,use_dijso,use_dijU,use_dijxc,use_dijxc_val
 logical :: has_phase,qne0
 real(dp) :: vi,vr,VUKS,vxcij22,vxcij22_i,vxcij,vxcij_i,vxcijhat,vxcijhat_i,tmp
 character(len=3) :: pertstrg
 character(len=500) :: msg
!arrays
 integer,allocatable :: idum(:),indklmn(:,:)
 logical,allocatable :: lmselect(:)
 real(dp) :: rdum(1),tsec(2)
 real(dp),allocatable :: coeffpawu(:)
 real(dp),allocatable :: buffer(:)
 real(dp),allocatable :: dij0(:),dijexxc(:),dijhat(:),dijhat_tmp(:),dijpawu(:),dijsym(:),dijsymU(:,:)
 real(dp),allocatable :: dijxc(:),dijxc_hat(:),dijxchat_tmp(:),dijxc_val(:)
 real(dp),allocatable :: ff(:),gg(:),prod(:),prodxchat(:),vpawu(:,:,:)
 real(dp),allocatable :: vxcij1(:),vxcij2(:),vxcijtot(:),vxcij1_val(:),vxcval_ij(:)
 real(dp),allocatable :: yylmr(:,:)
 character(len=7),parameter :: dspin(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(561,1,tsec)

 if (paw_ij(1)%ndij==4.and.paw_ij(1)%cplex_dij/=2) then
   msg='  invalid cplex size for Dij (4 Dij components) !'
   MSG_BUG(msg)
 end if
 if (nspden==4.and.cplex==2) then
   msg='  nspden=4 probably not compatible with cplex=2 !'
   MSG_BUG(msg)
 end if
 if (paw_ij(1)%cplex/=paw_an(1)%cplex) then
   msg='  paw_ij()%cplex and paw_an()%cplex must be equal !'
   MSG_BUG(msg)
 end if
 if (ipert<=0.and.paw_ij(1)%cplex/=1) then
   msg='  cplex must be 1 for GS calculations !'
   MSG_BUG(msg)
 end if
 if (paw_ij(1)%cplex_dij<cplex) then
   msg='  cplex_dij must be >= cplex !'
   MSG_BUG(msg)
 end if
 if (paw_ij(1)%cplex/=cplex) then
   msg='  paw_ij()%cplex must be equal to cplex !'
   MSG_BUG(msg)
 end if
 if(paw_ij(1)%has_dij==0) then
   msg='  dij must be allocated !'
   MSG_BUG(msg)
 end if
 if(paw_ij(1)%has_dijhartree==0.and. .not.(ipert==natom_tot+1.or.ipert==natom_tot+5)) then
   msg='  dijhartree must be allocated !'
   MSG_BUG(msg)
 end if
 if ((paw_ij(1)%cplex==2).and.&
& ((paw_an(1)%has_vxcval>0.and.paw_ij(1)%has_dijxc_val==1).or.&
& (paw_an(1)%has_vxc>0.and.paw_ij(1)%has_dijxc==1))) then
   msg =' Computation of dijxc/dijxcval not compatible with cplex=2 !'
   MSG_BUG(msg)
 end if
 ipositron=0
 if (present(electronpositron)) then
   ipositron=electronpositron_calctype(electronpositron)
   if (ipositron==1.and.pawtab(1)%has_kij/=2) then
     msg=' kij must be in memory for electronpositron%calctype=1 !'
     MSG_BUG(msg)
   end if
 end if
 if (mpi_enreg%nproc_atom>1) then
   if (natom/=mpi_enreg%natom) then
     MSG_BUG("natom not equal to mpi_enreg%natom !")
   end if
 end if
 qne0=(qphon(1)**2+qphon(2)**2+qphon(3)**2>=1.d-15)
 if (pawfgrtab(1)%rfgd_allocated==0.and.ipert>0.and.ipert<=natom_tot.and.qne0) then
   MSG_BUG('  pawfgrtab()%rfgd array must be allocated  !')
 end if

!----- Various initializations
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 nsppol=paw_ij(1)%nsppol
 nsploop=nsppol;if (paw_ij(1)%ndij==4) nsploop=4
 usexcnhat=maxval(pawtab(1:ntypat)%usexcnhat)
 pertstrg=" ";if (ipert>0) pertstrg="(1)"
 VUKS=zero;npts=0

 if (abs(pawprtvol)>=1) then
   write(msg, '(4a)')ch10,' ==== In pawdij: several values of Dij',&
&   trim(pertstrg),' (Hartree) ============'
   call wrtout(std_out,msg,'COLL')
 end if

!----- Preliminary computation of Ylm,Ylpmp (only if angular mesh)
 if (ipert/=natom_tot+1.and.ipert/=natom_tot+5) then
   if (pawxcdev==0) then
     npts=pawang%angl_size
     ABI_ALLOCATE(yylmr,(pawang%l_max**2*(pawang%l_max**2+1)/2,npts))
     do ipts=1,npts
       do jlm=1,pawang%l_max**2
         j0lm=jlm*(jlm-1)/2
         do ilm=1,jlm
           klm=j0lm+ilm
           yylmr(klm,ipts)=pawang%ylmr(ilm,ipts)*pawang%ylmr(jlm,ipts)
         end do
       end do
     end do
   end if
 end if

!These statements have to be removed when temporary parallelism disappear
 use_dijxc=0;use_dijxc_val=0;use_dijso=0;use_dijU=0;use_dijhat=0

!------------------------------------------------------------------------
!----- Big loop over atoms
!------------------------------------------------------------------------

 do iatom=1,natom
   iatom_tot=iatom;if (mpi_enreg%nproc_atom>1) iatom_tot=mpi_enreg%atom_indx(iatom)

!  -----------------------------------------------------------------------
!  ----------- Allocations and initializations
!  -----------------------------------------------------------------------

   itypat=typat(iatom_tot)
   l_size=pawtab(itypat)%l_size
   mesh_size=pawrad(itypat)%mesh_size
   lmn_size=paw_ij(iatom)%lmn_size
   lmn2_size=paw_ij(iatom)%lmn2_size
   lm_size=paw_an(iatom)%lm_size
   ij_size=pawtab(itypat)%ij_size
   cplex_dij=paw_ij(iatom)%cplex_dij
   paw_ij(iatom)%dij(:,:)=zero
   nfgd=pawfgrtab(iatom)%nfgd

   need_dijxc=0;need_dijxc_val=0
   need_dijso=0;need_dijU=0
   need_dijhat=0
   has_phase=.false.

!  Nothing to do for some perturbations (RF case)
   if (ipert/=natom_tot+1.and.ipert/=natom_tot+5) then

!    if has_dijxc==2, dijxc has already been computed...
     if (paw_an(iatom)%has_vxc>0.and.paw_ij(iatom)%has_dijxc==1) then
       need_dijxc=1;paw_ij(iatom)%dijxc(:,:)=zero
       ABI_ALLOCATE(dijxc_hat,(cplex_dij*lmn2_size))
     end if
     if (paw_ij(iatom)%has_dijxc==2.or.need_dijxc==1) use_dijxc=1  !To be removed one day
     ABI_ALLOCATE(dijxc,(cplex_dij*lmn2_size))

!    if has_dijxc_val==2, dijxc_val has already been computed...
     if (paw_an(iatom)%has_vxcval>0.and.paw_ij(iatom)%has_dijxc_val==1) then
       need_dijxc_val=1;paw_ij(iatom)%dijxc_val(:,:)=zero
       ABI_ALLOCATE(dijxc_val,(cplex_dij*lmn2_size))
     end if
     if (paw_ij(iatom)%has_dijxc_val==2.or.need_dijxc_val==1) use_dijxc_val=1  !To be removed one day

!    if has_dijso==2, dijso has already been computed...
     if (paw_ij(iatom)%has_dijso==1.and.ipert<=0.and.ipositron/=1) then
       need_dijso=1;paw_ij(iatom)%dijso(:,:)=zero
     end if
     if (paw_ij(iatom)%has_dijso==2.or.need_dijso==1) use_dijso=1  !To be removed one day
     if (pawspnorb>0.and.paw_ij(iatom)%has_dijso<=1.and.ipert<=0.and.ipositron/=1) then
       if (paw_ij(iatom)%has_dijso==0)  then
         ABI_ALLOCATE(paw_ij(iatom)%dijso,(cplex_dij*lmn2_size,paw_ij(iatom)%ndij))
       end if
       if (paw_an(iatom)%has_vhartree>0) then
         call pawdijso(iatom,itypat,natom,ntypat,paw_an,paw_ij,pawang,pawrad,pawtab,pawxcdev,dtset%spnorbscl)
       else
         paw_ij(iatom)%dijso(:,:)=zero
       end if
     end if

!    if has_dijU==2, dijU has already been computed...
     if (paw_ij(iatom)%has_dijU==1.and.ipert<=0.and.ipositron/=1) then
       need_dijU=1;paw_ij(iatom)%dijU(:,:)=zero
     end if
     if (paw_ij(iatom)%has_dijU==2.or.need_dijU==1) use_dijU=1  !To be removed one day
     if (pawtab(itypat)%usepawu>0.and.ipert<=0.and.ipositron/=1)  then
       ABI_ALLOCATE(dijpawu,(cplex_dij*lmn2_size))
     end if
     if (pawtab(itypat)%useexexch>0.and.ipert<=0.and.ipositron/=1)  then
       ABI_ALLOCATE(dijexxc,(cplex_dij*lmn2_size))
     end if
     if (pawtab(itypat)%usepawu>0.and.paw_ij(iatom)%ndij==4.and.ipert<=0)  then
       ABI_ALLOCATE(dijsymU,(cplex_dij*lmn2_size,4))
     end if

!    if has_dijhat==2, dijhat has already been computed...
     if (paw_ij(iatom)%has_dijhat==1) then
       need_dijhat=1;paw_ij(iatom)%dijhat(:,:)=zero
     end if
     if (paw_ij(iatom)%has_dijhat==2.or.need_dijhat==1) use_dijhat=1  !To be removed one day
     ABI_ALLOCATE(dijhat,(cplex_dij*lmn2_size))

     if (paral_kgb==1.and.mpi_enreg%nproc_atom==1.and.mod(iatom-1,mpi_enreg%nproc_band)/=mpi_enreg%me_band) then
       ABI_DEALLOCATE(dijhat)
       ABI_DEALLOCATE(dijxc)
       if (need_dijxc==1)  then
         ABI_DEALLOCATE(dijxc_hat)
       end if
       if (need_dijxc_val==1)  then
         ABI_DEALLOCATE(dijxc_val)
       end if
       if (pawspnorb>0.and.paw_ij(iatom)%has_dijso==0.and.ipert<=0)  then
         ABI_DEALLOCATE(paw_ij(iatom)%dijso)
       end if
       if (pawtab(itypat)%usepawu>0.and.ipert<=0.and.ipositron/=1)  then
         ABI_DEALLOCATE(dijpawu)
       end if
       if (pawtab(itypat)%usepawu>0.and.paw_ij(iatom)%ndij==4.and.ipert<=0)  then
         ABI_DEALLOCATE(dijsymU)
       end if
       if (pawtab(itypat)%useexexch>0.and.ipert<=0.and.ipositron/=1)  then
         ABI_DEALLOCATE(dijexxc)
       end if
       cycle
     end if

     ABI_ALLOCATE(ff,(mesh_size))
     if (paw_ij(1)%cplex==2)  then
       ABI_ALLOCATE(gg,(mesh_size))
     end if
     ABI_ALLOCATE(indklmn,(6,lmn2_size))
     indklmn(:,:)=pawtab(itypat)%indklmn(:,:)

!    Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
     if ((pawfgrtab(iatom)%gylm_allocated==0).or.&
&     ((ipert==iatom_tot).and.(pawfgrtab(iatom)%gylmgr_allocated==0))) then
       optgr0=0;optgr1=0
       if (pawfgrtab(iatom)%gylm_allocated==0) then
         if (associated(pawfgrtab(iatom)%gylm))  then
           ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
         end if
         ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(nfgd,lm_size))
         pawfgrtab(iatom)%gylm_allocated=2;optgr0=1
       end if
       if ((ipert==iatom_tot).and.(pawfgrtab(iatom)%gylmgr_allocated==0)) then
         if (associated(pawfgrtab(iatom)%gylmgr))  then
           ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
         end if
         ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(3,nfgd,lm_size))
         pawfgrtab(iatom)%gylmgr_allocated=2;optgr1=1
       end if
       call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,rdum,&
&       lm_size,nfgd,optgr0,optgr1,0,pawtab(itypat),&
&       pawfgrtab(iatom)%rfgd,pawfgrtab(iatom)%rfgd_allocated)
     end if

!    Eventually compute exp(-i.q.r) factors for the current atom (if not already done)
     if ((ipert==iatom_tot).and.qne0.and.(pawfgrtab(iatom)%expiqr_allocated==0)) then
       if (associated(pawfgrtab(iatom)%expiqr))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%expiqr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%expiqr,(2,nfgd))
       pawfgrtab(iatom)%expiqr_allocated=2
       call pawexpiqr(gprimd,pawfgrtab(iatom),qphon,xred(:,iatom_tot))
     end if
     has_phase=(qne0.and.ipert>0.and.pawfgrtab(iatom)%expiqr_allocated/=0)

   end if ! ipert/=natom_tot+1.and.ipert/=natom_tot+5

!  ------------------------------------------------------------------------
!  ----- Loop over density components
!  ------------------------------------------------------------------------
   do idij=1,nsploop

!    Print title
     if (abs(pawprtvol)>=1) then
       if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
         if (nspden==2.and.nsppol==1) then
           write(msg, '(2a,i3,3a)') ch10,&
&           ' >>>>>>>>>> Atom ',iatom,':',ch10,&
&           ' (antiferromagnetism case: only one spin component)'
         else if (paw_ij(iatom)%ndij==1) then
           write(msg, '(2a,i3,a)') ch10,&
&           ' >>>>>>>>>> Atom ',iatom,':'
         else
           write(msg, '(2a,i3,3a)') ch10,&
&           ' >>>>>>>>>> Atom ',iatom,' (component ',trim(dspin(idij+2*(nsploop/4))),'):'
         end if
         call wrtout(std_out,msg,'COLL')
       end if
     end if

!    Nothing to do for some perturbations (RF case)
     if (ipert/=natom_tot+1.and.ipert/=natom_tot+5) then

!      ------------------------------------------------------------------------
!      ----------- Load atomic Dij0 into Dij
!      ------------------------------------------------------------------------
!      No contribution to 1st-order Dij

       if (idij<=2.and.ipert<=0) then

         ABI_ALLOCATE(dij0,(lmn2_size))
         if (ipositron/=1) then
           dij0(:)=pawtab(itypat)%dij0(:)
         else
           dij0(:)=two*pawtab(itypat)%kij(:)-pawtab(itypat)%dij0(:)
         end if

         klmn1=1
         do klmn=1,lmn2_size
           paw_ij(iatom)%dij(klmn1,idij)=dij0(klmn)
           klmn1=klmn1+cplex_dij
         end do

         if (abs(pawprtvol)>=1) then
           if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
             write(msg, '(a)') '   ************ Dij atomic (Dij0) ***********'
             call wrtout(std_out,msg,'COLL')
             call print_ij(dij0,lmn2_size,1,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
           end if
         end if

         ABI_DEALLOCATE(dij0)

       end if

!      ------------------------------------------------------------------------
!      ----------- Add Dij_Hartree to Dij
!      ------------------------------------------------------------------------

       if (idij<=2) then

         if (cplex==1) then
           klmn1=1
           do klmn=1,lmn2_size
             paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+paw_ij(iatom)%dijhartree(klmn)
             klmn1=klmn1+cplex_dij
           end do
         else
           paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij)+paw_ij(iatom)%dijhartree(:)
         end if

         if (abs(pawprtvol)>=1) then
           if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
             write(msg, '(3a)')'   ************** Dij Hartree',&
&             trim(pertstrg),' ***************'
             call wrtout(std_out,msg,'COLL')
             if (ipert==0) then
               call print_ij(paw_ij(iatom)%dijhartree,lmn2_size,cplex,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
             else
               call print_ij(paw_ij(iatom)%dijhartree,lmn2_size,cplex,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,&
&               opt_sym=1)
             end if
           end if
         end if

       end if

!      ------------------------------------------------------------------------
!      ----------- Add Dij_xc to Dij
!      ------------------------------------------------------------------------

!      cplex is for RF, cplex_dij is for non-collinear (nspinor==2)
       if (idij<=nsppol.or.(nspden==4.and.idij<=3).or.cplex==2) then
         ABI_ALLOCATE(vxcijtot,(cplex*lmn2_size))
         if (need_dijxc_val==1)  then
           ABI_ALLOCATE(vxcval_ij,(cplex*lmn2_size))
         end if

         idijend=idij+idij/3;if (cplex==2) idijend=idij
         do ispden=idij,idijend

           vxcijtot=zero
           if (need_dijxc_val==1) vxcval_ij=zero

!          ================================================
!          ===== First formalism: use (l,m) moments for vxc
!          ================================================
           if (pawxcdev/=0) then

             ABI_ALLOCATE(vxcij1,(cplex*ij_size))
             if (need_dijxc_val==1)  then
               ABI_ALLOCATE(vxcij1_val,(cplex*ij_size))
             end if
             ABI_ALLOCATE(lmselect,(lm_size))
             lmselect(:)=paw_an(iatom)%lmselect(:)
             if (ipositron/=0) lmselect(1:lm_size)=(lmselect(1:lm_size).or.electronpositron%lmselect_ep(1:lm_size,iatom))

             do klm=1,lm_size
!              Summing over klm moments.
               if (lmselect(klm)) then

!                ===== Vxc_ij_1 (tmp) =====
                 vxcij1=zero
                 if (cplex==1) then
                   do kln=1,ij_size
                     ff(1:mesh_size)= paw_an(iatom)%vxc1(1:mesh_size,klm,ispden)&
&                     *pawtab(itypat)%phiphj(1:mesh_size,kln)&
&                     - paw_an(iatom)%vxct1(1:mesh_size,klm,ispden)&
&                     *pawtab(itypat)%tphitphj(1:mesh_size,kln)
                     call simp_gen(vxcij1(kln),ff,pawrad(itypat))
                   end do
                 else
                   do kln=1,ij_size
                     do ir=1,mesh_size
                       ir1=2*ir
                       ff(ir)= paw_an(iatom)%vxc1(ir1-1,klm,ispden)&
&                       *pawtab(itypat)%phiphj(ir,kln)&
&                       - paw_an(iatom)%vxct1(ir1-1,klm,ispden)&
&                       *pawtab(itypat)%tphitphj(ir,kln)
                       gg(ir)= paw_an(iatom)%vxc1(ir1,klm,ispden)&
&                       *pawtab(itypat)%phiphj(ir,kln)&
&                       - paw_an(iatom)%vxct1(ir1,klm,ispden)&
&                       *pawtab(itypat)%tphitphj(ir,kln)
                     end do
                     call simp_gen(vxcij1(2*kln-1),ff,pawrad(itypat))
                     call simp_gen(vxcij1(2*kln  ),gg,pawrad(itypat))
                   end do
                 end if

!                ==== If required calculate valence-only onsite matrix elements ====
                 if (need_dijxc_val==1) then ! compatible only with cplex=1
                   vxcij1_val(:)=zero
                   do kln=1,ij_size
                     ff(1:mesh_size)= &
&                     paw_an(iatom)%vxc1_val (1:mesh_size,klm,ispden)*pawtab(itypat)%phiphj  (1:mesh_size,kln)&
&                     -paw_an(iatom)%vxct1_val(1:mesh_size,klm,ispden)*pawtab(itypat)%tphitphj(1:mesh_size,kln)
                     call simp_gen(vxcij1_val(kln),ff,Pawrad(itypat))
                   end do
                 end if

!                ===== Vxc_ij_2 (tmp) =====
                 vxcij22=zero;vxcij22_i=zero
                 if (usexcnhat/=0) then
                   ll=1+int(sqrt(dble(klm)-0.1))
                   if (cplex==1) then
                     ff(1:mesh_size)=paw_an(iatom)%vxct1(1:mesh_size,klm,ispden)&
&                     *pawtab(itypat)%shapefunc(1:mesh_size,ll)&
&                     *pawrad(itypat)%rad(1:mesh_size)**2
                     call simp_gen(vxcij22,ff,pawrad(itypat))
                   else
                     do ir=1,mesh_size
                       ir1=2*ir
                       tmp=pawtab(itypat)%shapefunc(ir,ll)*pawrad(itypat)%rad(ir)**2
                       ff(ir)=paw_an(iatom)%vxct1(ir1-1,klm,ispden)*tmp
                       gg(ir)=paw_an(iatom)%vxct1(ir1  ,klm,ispden)*tmp
                     end do
                     call simp_gen(vxcij22  ,ff,pawrad(itypat))
                     call simp_gen(vxcij22_i,gg,pawrad(itypat))
                   end if
                 end if

!                ===== Accumulate over klm moments Vxc_ij_1 and Vxc_ij_2 =====
!                ===== into total Vxc_ij                                 =====
                 if (cplex==1) then
                   do klmn=1,lmn2_size
                     klm1=indklmn(1,klmn);kln=indklmn(2,klmn)
                     vxcij=zero
                     isel=pawang%gntselect(klm,klm1)
                     if (isel>0) vxcij=vxcij1(kln)*pawang%realgnt(isel)
                     vxcijhat=pawtab(itypat)%qijl(klm,klmn)*vxcij22
!                    Accumulate into total Vxc_ij
                     vxcijtot(klmn)=vxcijtot(klmn)+vxcij-vxcijhat
!                    Store valence-only matrix elements
                     if (need_dijxc_val==1) then
                       if (isel>0) vxcval_ij(klmn)=vxcval_ij(klmn) + vxcij1_val(kln)*pawang%realgnt(isel)
                     end if
                   end do ! Loop klmn
                 else
                   klmn1=1
                   do klmn=1,lmn2_size
                     klm1=indklmn(1,klmn);kln=indklmn(2,klmn)
                     vxcij=zero;vxcij_i=zero
                     isel=pawang%gntselect(klm,klm1)
                     if (isel>0) then
                       vxcij  =vxcij1(2*kln-1)*pawang%realgnt(isel)
                       vxcij_i=vxcij1(2*kln  )*pawang%realgnt(isel)
                     end if
                     vxcijhat  =pawtab(itypat)%qijl(klm,klmn)*vxcij22
                     vxcijhat_i=pawtab(itypat)%qijl(klm,klmn)*vxcij22_i
!                    Accumulate into total Vxc_ij
                     vxcijtot(klmn1  )=vxcijtot(klmn1   )+vxcij  -vxcijhat
                     vxcijtot(klmn1+1)=vxcijtot(klmn1+1 )+vxcij_i-vxcijhat_i
                     klmn1=klmn1+cplex
                   end do ! Loop klmn
                 end if

               end if
             end do  ! Loop klm

             ABI_DEALLOCATE(lmselect)
             ABI_DEALLOCATE(vxcij1)
             if (need_dijxc_val==1)  then
               ABI_DEALLOCATE(vxcij1_val)
             end if

!            ================================================
!            ===== Second formalism: use vxc on r,theta,phi
!            ================================================
           else

             ABI_ALLOCATE(vxcij1,(cplex*ij_size))
             ABI_ALLOCATE(vxcij2,(cplex*l_size))
             if (need_dijxc_val==1)  then
               ABI_ALLOCATE(vxcij1_val,(cplex*ij_size))
             end if

!            ===== Loop on angular mesh =====
             do ipts=1,npts

!              ===== Vxc_ij_1 (tmp) =====
               vxcij1=zero
               if (cplex==1) then
                 do kln=1,ij_size
                   ff(1:mesh_size)= paw_an(iatom)%vxc1(1:mesh_size,ipts,ispden)&
&                   *pawtab(itypat)%phiphj(1:mesh_size,kln)&
&                   - paw_an(iatom)%vxct1(1:mesh_size,ipts,ispden)&
&                   *pawtab(itypat)%tphitphj(1:mesh_size,kln)
                   call simp_gen(vxcij1(kln),ff,pawrad(itypat))
                 end do
               else
                 do kln=1,ij_size
                   do ir=1,mesh_size
                     ir1=2*ir
                     ff(ir)= paw_an(iatom)%vxc1(ir1-1,ipts,ispden)&
&                     *pawtab(itypat)%phiphj(ir,kln)&
&                     - paw_an(iatom)%vxct1(ir1-1,ipts,ispden)&
&                     *pawtab(itypat)%tphitphj(ir,kln)
                     gg(ir)= paw_an(iatom)%vxc1(ir1,ipts,ispden)&
&                     *pawtab(itypat)%phiphj(ir,kln)&
&                     - paw_an(iatom)%vxct1(ir1,ipts,ispden)&
&                     *pawtab(itypat)%tphitphj(ir,kln)
                   end do
                   call simp_gen(vxcij1(2*kln-1),ff,pawrad(itypat))
                   call simp_gen(vxcij1(2*kln  ),gg,pawrad(itypat))
                 end do
               end if

!              ==== If required calculate valence-only matrix elements ====
               if (need_dijxc_val==1) then ! compatible only with cplex=1
                 vxcij1_val(:)=zero
                 do kln=1,ij_size
                   ff(1:mesh_size)= &
&                   paw_an(iatom)%vxc1_val (1:mesh_size,ipts,ispden)*pawtab(itypat)%phiphj  (1:mesh_size,kln)&
&                   -paw_an(iatom)%vxct1_val(1:mesh_size,ipts,ispden)*pawtab(itypat)%tphitphj(1:mesh_size,kln)
                   call simp_gen(vxcij1_val(kln),ff,Pawrad(itypat))
                 end do
               end if

!              ===== Vxc_ij_2 (tmp) =====
               vxcij2=zero
               if (usexcnhat/=0) then
                 if (cplex==1) then
                   do ils=1,l_size
                     ff(1:mesh_size)=paw_an(iatom)%vxct1(1:mesh_size,ipts,ispden)&
&                     *pawtab(itypat)%shapefunc(1:mesh_size,ils)&
&                     *pawrad(itypat)%rad(1:mesh_size)**2
                     call simp_gen(vxcij2(ils),ff,pawrad(itypat))
                   end do
                 else
                   do ils=1,l_size
                     do ir=1,mesh_size
                       ir1=2*ir
                       tmp=pawtab(itypat)%shapefunc(ir,ils)*pawrad(itypat)%rad(ir)**2
                       ff(ir)=paw_an(iatom)%vxct1(ir1-1,ipts,ispden)*tmp
                       gg(ir)=paw_an(iatom)%vxct1(ir1  ,ipts,ispden)*tmp
                     end do
                     call simp_gen(vxcij2(2*ils-1),ff,pawrad(itypat))
                     call simp_gen(vxcij2(2*ils  ),gg,pawrad(itypat))
                   end do
                 end if
               end if

!              ===== Integrate Vxc_ij_1 and Vxc_ij_2 over the angular mesh =====
!              ===== and accummulate in total Vxc_ij                       =====
               if (cplex==1) then
                 do klmn=1,lmn2_size
                   klm=indklmn(1,klmn);kln=indklmn(2,klmn)
                   lmin=indklmn(3,klmn);lmax=indklmn(4,klmn)
                   vxcij=vxcij1(kln)*pawang%angwgth(ipts)*yylmr(klm,ipts)*four_pi
                   vxcijhat=zero
                   if (usexcnhat/=0) then
                     do ils=lmin,lmax,2
                       lm0=ils**2+ils+1
                       do mm=-ils,ils
                         ilslm=lm0+mm;isel=pawang%gntselect(ilslm,klm)
                         if (isel>0) vxcijhat=vxcijhat+vxcij2(ils+1)*pawang%angwgth(ipts)&
&                         *pawtab(itypat)%qijl(ilslm,klmn)*four_pi&
&                         *pawang%ylmr(ilslm,ipts)
                       end do
                     end do
                   end if
!                  Accumulate into total Vxc_ij
                   vxcijtot(klmn)=vxcijtot(klmn)+vxcij-vxcijhat
!                  Store valence-only matrix elements
                   if (need_dijxc_val==1) then
                     tmp=vxcij1_val(kln)*pawang%angwgth(ipts)*yylmr(klm,ipts)*four_pi
                     vxcval_ij(klmn)=vxcval_ij(klmn)+tmp
                   end if
                 end do ! Loop klmn
               else
                 klmn1=1
                 do klmn=1,lmn2_size
                   klm=indklmn(1,klmn);kln=indklmn(2,klmn)
                   lmin=indklmn(3,klmn);lmax=indklmn(4,klmn)
                   vxcij  =vxcij1(2*kln-1)*pawang%angwgth(ipts)*yylmr(klm,ipts)*four_pi
                   vxcij_i=vxcij1(2*kln  )*pawang%angwgth(ipts)*yylmr(klm,ipts)*four_pi
                   vxcijhat=zero;vxcijhat_i=zero
                   if (usexcnhat/=0) then
                     do ils=lmin,lmax,2
                       lm0=ils**2+ils+1
                       do mm=-ils,ils
                         ilslm=lm0+mm;isel=pawang%gntselect(ilslm,klm)
                         if (isel>0) then
                           ils1=2*(ils+1)
                           tmp=pawang%angwgth(ipts)*pawang%ylmr(ilslm,ipts)&
&                           *pawtab(itypat)%qijl(ilslm,klmn)*four_pi
                           vxcijhat  =vxcijhat  +vxcij2(ils1-1)*tmp
                           vxcijhat_i=vxcijhat_i+vxcij2(ils1  )*tmp
                         end if
                       end do
                     end do
                   end if
!                  Accumulate into total Vxc_ij
                   vxcijtot(klmn1  )=vxcijtot(klmn1   )+vxcij  -vxcijhat
                   vxcijtot(klmn1+1)=vxcijtot(klmn1+1 )+vxcij_i-vxcijhat_i
                   klmn1=klmn1+cplex
                 end do ! Loop klmn
               end if
             end do  ! Loop ipts

             ABI_DEALLOCATE(vxcij1)
             ABI_DEALLOCATE(vxcij2)
             if (need_dijxc_val==1)  then
               ABI_DEALLOCATE(vxcij1_val)
             end if

           end if  ! choice XC

           if (ispden<3) then
             dijxc(1:cplex*lmn2_size)=vxcijtot(1:cplex*lmn2_size)
             if (need_dijxc_val==1) dijxc_val(1:cplex*lmn2_size)=vxcval_ij(1:cplex*lmn2_size)
           else
             if (cplex==cplex_dij) then
               dijxc(1:cplex*lmn2_size)=vxcijtot(1:cplex*lmn2_size)
               if (need_dijxc_val==1) dijxc_val(1:cplex*lmn2_size)=vxcval_ij(1:cplex*lmn2_size)
             else ! Note that cplex_dij>=cplex
               klmn1=max(1,ispden-2)
               do klmn=1,lmn2_size
                 dijxc(klmn1)=vxcijtot(klmn)
                 klmn1=klmn1+cplex_dij
               end do
             end if
             if (need_dijxc_val==1) then ! compatible only with cplex=1
               klmn1=max(1,ispden-2)
               do klmn=1,lmn2_size
                 dijxc_val(klmn1)=vxcval_ij(klmn)
                 klmn1=klmn1+cplex_dij
               end do
             end if
           end if

         end do ! ispden
         ABI_DEALLOCATE(vxcijtot)
         if (need_dijxc_val==1)  then
           ABI_DEALLOCATE(vxcval_ij)
         end if

       else if (nspden==4.and.idij==4) then
         klmn1=2
         do klmn=1,lmn2_size
           dijxc(klmn1)=-dijxc(klmn1)
           klmn1=klmn1+cplex_dij
         end do
         if (need_dijxc_val==1) then
           klmn1=2
           do klmn=1,lmn2_size
             dijxc_val(klmn1)=-dijxc_val(klmn1)
             klmn1=klmn1+cplex_dij
           end do
         end if
       end if
       if ((idij<=nsppol.or.idij==2).and.cplex==1) then
         klmn1=1
         do klmn=1,lmn2_size
           paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+dijxc(klmn)
           klmn1=klmn1+cplex_dij
         end do
         if (need_dijxc==1) then
           klmn1=1
           do klmn=1,lmn2_size
             paw_ij(iatom)%dijxc(klmn1,idij)=dijxc(klmn)
             klmn1=klmn1+cplex_dij
           end do
         end if
         if (need_dijxc_val==1) then
           klmn1=1
           do klmn=1,lmn2_size
             paw_ij(iatom)%dijxc_val(klmn1,idij)=dijxc_val(klmn)
             klmn1=klmn1+cplex_dij
           end do
         end if
       else if (nspden==4.or.cplex==2) then  ! cplex=cplex_dij
         paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij)+dijxc(:)
         if (need_dijxc    ==1) paw_ij(iatom)%dijxc    (:,idij)=dijxc(:)
         if (need_dijxc_val==1) paw_ij(iatom)%dijxc_val(:,idij)=dijxc_val(:)
       end if

       if ((abs(pawprtvol)>=1).and.(idij<=2.or.nspden==4)) then
         if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
           write(msg, '(3a)')'   ****************** Dij_xc',&
&           trim(pertstrg),' ****************'
           call wrtout(std_out,msg,'COLL')
           if ((idij<=nsppol.or.idij==2).and.cplex==1) then
             call print_ij(dijxc(1:lmn2_size),lmn2_size,1,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
           else
             call print_ij(dijxc,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,opt_sym=1)
           end if
         end if
       end if

!      ------------------------------------------------------------------------
!      ----------- Add Dij_hat to Dij
!      ------------------------------------------------------------------------

       if (idij<=nsppol.or.(nspden==4.and.idij<=3).or.cplex==2) then

         ABI_ALLOCATE(dijhat_tmp,(cplex*lmn2_size))
         if (need_dijxc==1)  then
           ABI_ALLOCATE(dijxchat_tmp,(cplex*lmn2_size))
         end if

         idijend=idij+idij/3;if (cplex==2) idijend=idij
         do ispden=idij,idijend

           dijhat_tmp=zero
!          Compute Int[V^(alpha,beta)(r).g_l(r).Y_lm(r)]
!          Remember: if nspden=4, V is stored as : V^11, V^22, V^12, i.V^21
           ABI_ALLOCATE(prod,(cplex*lm_size))
           prod=zero
           if (usexcnhat/=0) then
             if (has_phase) then
               if (cplex==1) then
                 do ilslm=1,lm_size
                   do ic=1,nfgd
                     vr=vtrial(pawfgrtab(iatom)%ifftsph(ic),ispden)
                     prod(ilslm)=prod(ilslm)+vr*pawfgrtab(iatom)%gylm(ic,ilslm)&
&                     *pawfgrtab(iatom)%expiqr(1,ic)
                   end do
                 end do
               else
                 ilslm1=1
                 do ilslm=1,lm_size
                   do ic=1,nfgd
                     jc=2*pawfgrtab(iatom)%ifftsph(ic)
                     vr=vtrial(jc-1,ispden);vi=vtrial(jc,ispden)
                     prod(ilslm1  )=prod(ilslm1  )+pawfgrtab(iatom)%gylm(ic,ilslm)&
&                     *(vr*pawfgrtab(iatom)%expiqr(1,ic)-vi*pawfgrtab(iatom)%expiqr(2,ic))
                     prod(ilslm1+1)=prod(ilslm1+1)+pawfgrtab(iatom)%gylm(ic,ilslm)&
&                     *(vr*pawfgrtab(iatom)%expiqr(2,ic)+vi*pawfgrtab(iatom)%expiqr(1,ic))
                   end do
                   ilslm1=ilslm1+cplex
                 end do
               end if
             else ! no phase
               if (cplex==1) then
                 do ilslm=1,lm_size
                   do ic=1,nfgd
                     vr=vtrial(pawfgrtab(iatom)%ifftsph(ic),ispden)
                     prod(ilslm)=prod(ilslm)+vr*pawfgrtab(iatom)%gylm(ic,ilslm)
                   end do
                 end do
               else
                 ilslm1=1
                 do ilslm=1,lm_size
                   do ic=1,nfgd
                     jc=2*pawfgrtab(iatom)%ifftsph(ic)
                     vr=vtrial(jc-1,ispden);vi=vtrial(jc,ispden)
                     prod(ilslm1  )=prod(ilslm1  )+vr*pawfgrtab(iatom)%gylm(ic,ilslm)
                     prod(ilslm1+1)=prod(ilslm1+1)+vi*pawfgrtab(iatom)%gylm(ic,ilslm)
                   end do
                   ilslm1=ilslm1+cplex
                 end do
               end if
             end if
           else ! usexcnhat=0
             if (has_phase) then
               if (cplex==1) then
                 do ilslm=1,lm_size
                   do ic=1,nfgd
                     jc=pawfgrtab(iatom)%ifftsph(ic)
                     vr=vtrial(jc,ispden)-vxc(jc,ispden)
                     prod(ilslm)=prod(ilslm)+vr*pawfgrtab(iatom)%gylm(ic,ilslm)&
&                     *pawfgrtab(iatom)%expiqr(1,ic)
                   end do
                 end do
               else
                 ilslm1=1
                 do ilslm=1,lm_size
                   do ic=1,nfgd
                     jc=2*pawfgrtab(iatom)%ifftsph(ic)
                     vr=vtrial(jc-1,ispden)-vxc(jc-1,ispden)
                     vi=vtrial(jc  ,ispden)-vxc(jc  ,ispden)
                     prod(ilslm1  )=prod(ilslm1  )+pawfgrtab(iatom)%gylm(ic,ilslm)&
&                     *(vr*pawfgrtab(iatom)%expiqr(1,ic)-vi*pawfgrtab(iatom)%expiqr(2,ic))
                     prod(ilslm1+1)=prod(ilslm1+1)+pawfgrtab(iatom)%gylm(ic,ilslm)&
&                     *(vr*pawfgrtab(iatom)%expiqr(2,ic)+vi*pawfgrtab(iatom)%expiqr(1,ic))
                   end do
                   ilslm1=ilslm1+cplex
                 end do
               end if
             else ! no phase
               if (cplex==1) then
                 do ilslm=1,lm_size
                   do ic=1,nfgd
                     jc=pawfgrtab(iatom)%ifftsph(ic)
                     vr=vtrial(jc,ispden)-vxc(jc,ispden)
                     prod(ilslm)=prod(ilslm)+vr*pawfgrtab(iatom)%gylm(ic,ilslm)
                   end do
                 end do
               else
                 ilslm1=1
                 do ilslm=1,lm_size
                   do ic=1,nfgd
                     jc=2*pawfgrtab(iatom)%ifftsph(ic)-1
                     vr=vtrial(jc  ,ispden)-vxc(jc  ,ispden)
                     vi=vtrial(jc+1,ispden)-vxc(jc+1,ispden)
                     prod(ilslm1  )=prod(ilslm1  )+vr*pawfgrtab(iatom)%gylm(ic,ilslm)
                     prod(ilslm1+1)=prod(ilslm1+1)+vi*pawfgrtab(iatom)%gylm(ic,ilslm)
                   end do
                   ilslm1=ilslm1+cplex
                 end do
               end if
             end if
           end if ! usexcnhat
           if(mpi_enreg%paral_fft==1)then
             old_paral_level= mpi_enreg%paral_level
             mpi_enreg%paral_level=3
             call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
             call xsum_mpi(prod,spaceComm,ier)
             mpi_enreg%paral_level=old_paral_level
           end if

           if (need_dijxc==1) then
             dijxchat_tmp(:)=zero
!            === Evaluate prodxchat i.e $\sum_{lm} \int g_l Ylm v_xc[tn+nhat+tnc]dr$ on the FFT mesh ===
!            * It does not depend on ij
             ABI_ALLOCATE(prodxchat,(lm_size))
             prodxchat(:)=zero
             do ilslm=1,lm_size
               do ic=1,nfgd
                 prodxchat(ilslm)=prodxchat(ilslm)&
&                 +vxc(pawfgrtab(iatom)%ifftsph(ic),ispden)*pawfgrtab(iatom)%gylm(ic,ilslm)
               end do
             end do
             if (mpi_enreg%paral_fft==1)then
               old_paral_level= mpi_enreg%paral_level
               mpi_enreg%paral_level=3
               call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
               call xsum_mpi(prodxchat,spaceComm,ier)
               mpi_enreg%paral_level=old_paral_level
             end if
           end if

!          Compute Sum_(i,j)_LM { q_ij^L Int[V^(alpha,beta)(r).g_l(r).Y_lm(r)] }
           if (cplex==1) then
             do klmn=1,lmn2_size
               klm=indklmn(1,klmn)
               lmin=indklmn(3,klmn);lmax=indklmn(4,klmn)
               do ils=lmin,lmax,2
                 lm0=ils**2+ils+1
                 do mm=-ils,ils
                   ilslm=lm0+mm;isel=pawang%gntselect(lm0+mm,klm)
                   if (isel>0) dijhat_tmp(klmn)=dijhat_tmp(klmn)+prod(ilslm)*pawtab(itypat)%qijl(ilslm,klmn)
                   if (isel>0.and.need_dijxc==1) then
                     dijxchat_tmp(klmn)=dijxchat_tmp(klmn)+prodxchat(ilslm)*pawtab(itypat)%qijl(ilslm,klmn)
                   end if
                 end do
               end do
             end do
           else
             do klmn=1,lmn2_size
               klm=indklmn(1,klmn);klmn1=2*klmn-1
               lmin=indklmn(3,klmn);lmax=indklmn(4,klmn)
               do ils=lmin,lmax,2
                 lm0=ils**2+ils+1
                 do mm=-ils,ils
                   ilslm=lm0+mm;ilslm1=2*ilslm;isel=pawang%gntselect(lm0+mm,klm)
                   if (isel>0) dijhat_tmp(klmn1:klmn1+1)=dijhat_tmp(klmn1:klmn1+1) &
&                   +prod(ilslm1-1:ilslm1)*pawtab(itypat)%qijl(ilslm,klmn)
                 end do
               end do
             end do
           end if
           ABI_DEALLOCATE(prod)
           if (need_dijxc==1)  then
             ABI_DEALLOCATE(prodxchat)
           end if

           if (cplex==1) then
             if (ispden<3) then
               dijhat(1:cplex*lmn2_size)=dijhat_tmp(1:cplex*lmn2_size)*ucvol/dble(nfftot)
               if (need_dijxc==1) dijxc_hat(1:cplex*lmn2_size)=dijxchat_tmp(1:cplex*lmn2_size)*ucvol/dble(nfftot)
             else
               klmn1=max(1,ispden-2)
               do klmn=1,lmn2_size
                 dijhat(klmn1)=dijhat_tmp(klmn)*ucvol/dble(nfftot)
                 klmn1=klmn1+cplex_dij
               end do
               if (need_dijxc==1) then
                 klmn1=max(1,ispden-2)
                 do klmn=1,lmn2_size
                   dijxc_hat(klmn1)=dijxchat_tmp(klmn)*ucvol/dble(nfftot)
                   klmn1=klmn1+cplex_dij
                 end do
               end if
             end if
           else !cplex=2
             if (ispden<=3) then
               dijhat(1:cplex*lmn2_size)=dijhat_tmp(1:cplex*lmn2_size)*ucvol/dble(nfftot)
               if (need_dijxc==1) dijxc_hat(1:cplex*lmn2_size)=dijxchat_tmp(1:cplex*lmn2_size)*ucvol/dble(nfftot)
             else
!              Remember V(4) contains i.V^21
               klmn1=1
               do klmn=1,lmn2_size
                 dijhat(klmn1  )= dijhat_tmp(klmn+1)*ucvol/dble(nfftot)
                 dijhat(klmn1+1)=-dijhat_tmp(klmn  )*ucvol/dble(nfftot)
                 klmn1=klmn1+cplex_dij
               end do
               if (need_dijxc==1) then
                 do klmn=1,lmn2_size
                   dijxc_hat(klmn1  )= dijxchat_tmp(klmn+1)*ucvol/dble(nfftot)
                   dijxc_hat(klmn1+1)=-dijxchat_tmp(klmn  )*ucvol/dble(nfftot)
                   klmn1=klmn1+cplex_dij
                 end do
               end if
             end if
           end if

         end do !ispden

         ABI_DEALLOCATE(dijhat_tmp)
         if (need_dijxc==1)  then
           ABI_DEALLOCATE(dijxchat_tmp)
         end if

       else if (nspden==4.and.idij==4) then  ! Note that cplex=1 here
         klmn1=2
         do klmn=1,lmn2_size
           dijhat(klmn1)=-dijhat(klmn1)
           klmn1=klmn1+cplex_dij
         end do
         if (need_dijxc==1) then
           klmn1=2
           do klmn=1,lmn2_size
             dijxc_hat(klmn1)=-dijxc_hat(klmn1)
             klmn1=klmn1+cplex_dij
           end do
         end if
       end if

       if ((idij<=nsppol.or.idij==2).and.cplex==1) then
         klmn1=1
         do klmn=1,lmn2_size
           paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+dijhat(klmn)
           klmn1=klmn1+cplex_dij
         end do
         if (need_dijhat>0) then
           klmn1=1
           do klmn=1,lmn2_size
             paw_ij(iatom)%dijhat(klmn1,idij)=dijhat(klmn)
             klmn1=klmn1+cplex_dij
           end do
         end if
         if (need_dijxc==1) then
           klmn1=1
           do klmn=1,lmn2_size
             paw_ij(iatom)%dijxc(klmn1,idij)=paw_ij(iatom)%dijxc(klmn1,idij)+dijxc_hat(klmn)
             klmn1=klmn1+cplex_dij
           end do
         end if
       else if (nspden==4.or.cplex==2) then
         paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij)+dijhat(:)
         if (need_dijhat>0) paw_ij(iatom)%dijhat(:,idij)=dijhat(:)
         if (need_dijxc==1) paw_ij(iatom)%dijxc(:,idij)=paw_ij(iatom)%dijxc(:,idij)+dijxc_hat(:)
       end if

       if ((abs(pawprtvol)>=1).and.(idij<=2.or.nspden==4)) then
         if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
           write(msg, '(3a)')'   ****************** Dij_hat',&
&           trim(pertstrg),' ***************'
           call wrtout(std_out,msg,'COLL')
           if ((idij<=nsppol.or.idij==2).and.cplex==1) then
             call print_ij(dijhat(1:lmn2_size),lmn2_size,1,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
           else
             call print_ij(dijhat,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,opt_sym=1)
           end if
         end if
       end if

!      ------------------------------------------------------------------------
!      ----------- Add RF frozen Dij to Dij
!      ------------------------------------------------------------------------
!      ----------- RF only

       if (ipert>0.and.paw_ij(iatom)%has_dijfr==2) then

         paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij)+paw_ij(iatom)%dijfr(:,idij)

         if (need_dijhat>0) paw_ij(iatom)%dijhat(:,idij)=paw_ij(iatom)%dijhat(:,idij)+paw_ij(iatom)%dijfr(:,idij)

         if (abs(pawprtvol)>=1) then
           if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
             write(msg, '(3a)')'   ****************** Dij frozen',&
&             trim(pertstrg),' ***************'
             call wrtout(std_out,msg,'COLL')
             call print_ij(paw_ij(iatom)%dijfr,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,&
&             opt_sym=1)
           end if
         end if

       end if

!      ------------------------------------------------------------------------
!      ----------- Add Dij spin-orbit to Dij
!      ------------------------------------------------------------------------
!      No contribution to 1st-order Dij
       if (pawspnorb>0.and.ipert<=0.and.ipositron/=1) then

         paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij)+paw_ij(iatom)%dijso(:,idij)

!        Printing
         if (abs(pawprtvol)>=1) then
           if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
             write(msg, '(a)')'   ************** Dij SpinOrbit ************'
             call wrtout(std_out,msg,'COLL')
             call print_ij(paw_ij(iatom)%dijso(:,idij), &
&             lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,opt_sym=3)
           end if
         end if

       end if

!      ------------------------------------------------------------------------
!      ----------- Add Dij_{LDA+U} to Dij
!      ------------------------------------------------------------------------
!      Dijpawu^{\sigma}_{mi,ni,mj,nj}=
!      \sum_{m,m'} [vpawu^{\sigma}_{m,m'}*phiphjint_{ni,nj}^{m,m'}]=
!      [vpawu^{\sigma}_{mi,mj}*phiphjint_{ni,nj}]
!      ------------------------------------------------------------------------
!      No contribution to 1st-order Dij

       if (pawtab(itypat)%usepawu>0.and.ipert<=0.and.ipositron/=1) then

         lpawu=pawtab(itypat)%lpawu
         if (idij<=nsppol.or.(paw_ij(iatom)%ndij==4.and.idij<=3).or.cplex==2) then
           idijend=idij+idij/3
           do idijeff=idij,idijend ! if ndij==4, idijeff is used to compute updn and dnup contributions

             if(abs(pawprtvol)>=3) then
               write(msg,*) "pawdij, LDA+U calculation for iatom,ispden,itypat= ",iatom_tot,ispden,itypat
               call wrtout(std_out,msg,'COLL')
             end if

             ABI_ALLOCATE(vpawu,(paw_ij(iatom)%cplex_dij,lpawu*2+1,lpawu*2+1))
             if(pawtab(itypat)%usepawu<10) then ! if dmft, do not apply U in LDA+U
               call pawpupot(idijeff,paw_ij(iatom),pawprtvol,pawtab(itypat),vpawu,VUKS) ! idij=1,max(nsppol,nspinor**2)
             else
               vpawu=zero
               VUKS=zero
             end if
             dijpawu=zero
             do klmn=1,lmn2_size
               if(cplex_dij==1) then
                 klmn1=klmn
               else
                 klmn1=cplex_dij*klmn-1  ! klmn1=cplex_dij*klmn-cplex_dij/2
               end if
               im1=pawtab(itypat)%klmntomn(1,klmn)
               im2=pawtab(itypat)%klmntomn(2,klmn)
               in1=pawtab(itypat)%klmntomn(3,klmn)
               in2=pawtab(itypat)%klmntomn(4,klmn)
               lmin=pawtab(itypat)%indklmn(3,klmn)
               lmax=pawtab(itypat)%indklmn(4,klmn)
               if(lmin==0.and.lmax==2*lpawu) then
                 icount=in1+(in2*(in2-1))/2
                 if(pawtab(itypat)%ij_proj<icount)  then
                   write(msg, '(4a)' ) ch10,&
&                   ' pawdij : BUG -',ch10,&
&                   '  PAW+U: Problem in the loop for calculating dijpawu'
                   call wrtout(std_out,msg,'COLL')
                   call leave_new('COLL')
                 end if
                 ABI_ALLOCATE(coeffpawu,(cplex_dij))
!                coeffpawu(:)=vpawu(:,im1,im2) ! use real and imaginary part
                 coeffpawu(:)=vpawu(:,im2,im1) ! because of transposition in setnoccmmp (for the cplex_dij==2)
                 if(dtset%natvshift/=0.and.idij<3.and.im1==im2) then
                   coeffpawu(1)=coeffpawu(1)+fatvshift*dtset%atvshift(im1,idij,iatom_tot)
                 end if
                 if(cplex_dij==1) then   !cplex_dij=nspinor=1
                   dijpawu(klmn1)=pawtab(itypat)%phiphjint(icount)*coeffpawu(1) ! *dtset%userra
                 elseif (cplex_dij==2) then   !cplex_dij=nspinor=2
                   dijpawu(klmn1)=pawtab(itypat)%phiphjint(icount)*coeffpawu(1)
                   dijpawu(klmn1+1)=pawtab(itypat)%phiphjint(icount)*coeffpawu(2) !  spinor==2
                 end if
!                write(std_out,*) "cplex_dij",cplex_dij,im1,im2
!                write(std_out,*) "vpawu",vpawu(:,im2,im1)
!                write(std_out,'(a,2i4,f8.5,f8.5)') "dijpawu1",idij,klmn1,dijpawu(klmn1),dijpawu(klmn1+1)
                 ABI_DEALLOCATE(coeffpawu)
               end if
             end do ! klmn
             ABI_DEALLOCATE(vpawu)
!            dijsymU useful for printing
             if (ipert<=0.and.paw_ij(iatom)%ndij==4) then
               do klmn=1,cplex_dij*lmn2_size,2
                 dijsymU(klmn,idijeff)=dijpawu(klmn)
                 dijsymU(klmn+1,idijeff)=dijpawu(klmn+1)
               end do
             end if
           end do ! idijeff

!          else if (nspden==4.and.idij==4) then
!          klmn1=2
!          do klmn=1,lmn2_size
!          dijpawu(klmn1)=-dijpawu(klmn1) !  dij(dnup)=dij(updn)^*: change sign of imaginary part, real part stays identical
!          klmn1=klmn1+cplex_dij
!          end do
         end if


!        if ((idij<=nsppol.or.idij==2).and.cplex==1) then
!        if (paw_ij(iatom)%ndij<4) then
         if (cplex_dij==1)then
           klmn1=1
           do klmn=1,lmn2_size
             paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+dijpawu(klmn)
             klmn1=klmn1+cplex_dij
           end do
           if (need_dijU==1) then
             klmn1=1
             do klmn=1,lmn2_size
               paw_ij(iatom)%dijU(klmn1,idij)=dijpawu(klmn)
               klmn1=klmn1+cplex_dij
             end do
           end if
         else if (nspden==4.or.cplex==2.or.paw_ij(iatom)%ndij==4.or.cplex_dij==2) then
           if(idij<=2)  then
             paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij)+dijpawu(:)
             if (need_dijU==1) paw_ij(iatom)%dijU(:,idij)=dijpawu(:)
           else if(idij>2)  then
             paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij)+dijsymU(:,idij)
             if (need_dijU==1) paw_ij(iatom)%dijU(:,idij)=dijsymU(:,idij)
           end if
         end if

         do klmn=1,lmn2_size
           klmn1=cplex_dij*klmn-1  ! klmn1=cplex_dij*klmn-cplex_dij/2
         end do
         if ((abs(pawprtvol)>=1).and.(idij<=2.or.nspden==4)) then
           if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
             write(msg, '(3a)')'   ***************** Dij_LDA+U',&
&             trim(pertstrg),' **************'
             call wrtout(std_out,msg,'COLL')
             if ((idij<=nsppol.or.idij==2).and.cplex==1.and.cplex_dij==1) then
               call print_ij(dijpawu(1:lmn2_size),lmn2_size,1,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
             else
               if(idij<=2) then
!                call print_ij(dijpawu,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,asym_ij=dijsymU)
                 call print_ij(dijpawu,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,opt_sym=2)
               else
                 call print_ij(dijsymU(:,idij),lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,&
&                 pawprtvol,idum,-1.d0,1,opt_sym=2,asym_ij=dijsymU(:,7-idij))
               end if
             end if
           end if
         end if

       end if

!      ------------------------------------------------------------------------
!      ----------- Add Dij_{local exact-exchange} to Dij
!      ------------------------------------------------------------------------
!      No contribution to 1st-order Dij

       if (pawtab(itypat)%useexexch>0.and.ipert<=0.and.ipositron/=1) then

         if(nspden==4)  then
           write(msg, '(4a)' ) ch10,&
&           '  pawdenpot : ERROR -',ch10,&
&           '  Local exact-exch. not implemented for nspden=4 !'
           call wrtout(ab_out,msg,'COLL')
           call wrtout(std_out,  msg,'COLL')
           call leave_new('COLL')
         end if

         if (idij<=nsppol) then

           if(pawprtvol>=3) then
             write(msg,*) "pawdij, local exact-exchange calculation for iatom,idij",iatom_tot,idij
             call wrtout(std_out,  msg,'COLL')
           end if

           lexexch=pawtab(itypat)%lexexch

           ABI_ALLOCATE(vxcij1,(ij_size))
           ABI_ALLOCATE(vxcijtot,(lmn2_size))
           vxcijtot=zero
           do klm=1,lm_size
             if (paw_an(iatom)%lmselect(klm)) then
!              ===== Vxc_ij_1 (tmp) =====
               vxcij1=zero
               do jln=pawtab(itypat)%lnproju(1),pawtab(itypat)%lnproju(pawtab(itypat)%nproju)
                 j0ln=jln*(jln-1)/2
                 do iln=pawtab(itypat)%lnproju(1),jln
                   kln=j0ln+iln
                   ff(1:mesh_size)=paw_an(iatom)%vxc_ex(1:mesh_size,klm,idij) &
&                   *pawtab(itypat)%phiphj(1:mesh_size,kln)
                   call simp_gen(vxcij1(kln),ff,pawrad(itypat))
                 end do
               end do
!              ===== Contribution to total Vxc_ij =====
               do klmn=1,lmn2_size
                 lmin=pawtab(itypat)%indklmn(3,klmn)
                 lmax=pawtab(itypat)%indklmn(4,klmn)
                 if(lmin==0.and.lmax==2*lexexch) then
                   klm1=indklmn(1,klmn);kln=indklmn(2,klmn)
                   isel=pawang%gntselect(klm,klm1)
                   if (isel>0) vxcijtot(klmn)=vxcijtot(klmn)+vxcij1(kln)*pawang%realgnt(isel)
                 end if
               end do ! Loop klmn
             end if
           end do  ! Loop klm
           ABI_DEALLOCATE(vxcij1)

           if (abs(pawprtvol)>=4) then
             if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
               write(msg, '(a)')'   ** INFO ***** Dij_local_exact-exchange (vxcij) **********'
               call wrtout(std_out,msg,'COLL')
               call print_ij(vxcijtot,lmn2_size,1,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
             end if
           end if

           dijexxc=zero
           do klmn=1,lmn2_size
             in1=pawtab(itypat)%klmntomn(3,klmn)
             in2=pawtab(itypat)%klmntomn(4,klmn)
             lmin=pawtab(itypat)%indklmn(3,klmn)
             lmax=pawtab(itypat)%indklmn(4,klmn)
             if(lmin==0.and.lmax==2*lexexch) then
               icount=in1+(in2*(in2-1))/2
               if(pawtab(itypat)%ij_proj<icount)  then
                 write(msg, '(4a)' ) ch10,&
&                 ' pawdij : BUG -',ch10,&
&                 '  PAW exact-exchange: Problem in the loop for calculating dijexxc'
                 call wrtout(std_out,msg,'COLL')
                 call leave_new('COLL')
               end if
               dijexxc(klmn)=(paw_ij(iatom)%vpawx(1,klmn,idij)-vxcijtot(klmn))*pawtab(itypat)%exchmix
             end if
           end do
           ABI_DEALLOCATE(vxcijtot)

         end if

         if ((idij<=nsppol.or.idij==2).and.cplex==1) then
           klmn1=1
           do klmn=1,lmn2_size
             paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+dijexxc(klmn)
             klmn1=klmn1+cplex_dij
           end do
         else if (nspden==4.or.cplex==2) then
           paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij)+dijexxc(:)
         end if

         if ((abs(pawprtvol)>=1).and.(idij<=2.or.nspden==4)) then
           if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
             write(msg, '(3a)')'   ************ Dij_Exact exchange',trim(pertstrg),' **********'
             call wrtout(std_out,msg,'COLL')
             if ((idij<=nsppol.or.idij==2).and.cplex==1) then
               call print_ij(dijexxc(1:lmn2_size),lmn2_size,1,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
             else
               call print_ij(dijexxc,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,opt_sym=1)
             end if
           end if
         end if

       end if

!      ------------------------------------------------------------------------
!      ----------- Final printing
!      ------------------------------------------------------------------------

     end if ! ipert/=natom_tot+1.and.ipert/=natom_tot+5

     if (idij>=3.and.pawspnorb>0.and.ipert<=0.and.ipositron/=1) then
       ABI_ALLOCATE(dijsym,(cplex_dij*lmn2_size))
       dijsym(:)=paw_ij(iatom)%dij(:,3)
       if (idij==3) then
         do klmn=1,cplex_dij*lmn2_size,2
           dijsym(klmn)=dijsym(klmn)-two*paw_ij(iatom)%dijso(klmn,3)
!          todo:  should do a similar thing for LDA+U (just for printing),
!          current printing is not correct for the upper triangular part for the LDA+U part in total dij.
         end do
       end if
     end if

     if (abs(pawprtvol)>=1) then
       if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
         write(msg, '(3a)' )'   **********    TOTAL Dij',&
&         trim(pertstrg),' in Ha   **********'
         call wrtout(std_out,msg,'COLL')
         if (idij<=2.or.pawspnorb==0) then
           if (ipert==0) then
             call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,&
&             1,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*idij),1)
           else
             call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,&
&             1,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*idij),1,opt_sym=1)
           end if
         else
           if (ipert==0) then
             call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,&
&             1,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*idij),1,asym_ij=dijsym)
           else
             call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,&
&             1,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*idij),1)
           end if
         end if
         if (enunit>0) then
           write(msg, '(3a)' )'   **********    TOTAL Dij',&
&           trim(pertstrg),' in eV   **********'
           call wrtout(std_out,msg,'COLL')
           if (idij<=2.or.pawspnorb==0) then
             if (ipert==0) then
               call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,&
&               1,-1,idum,0,pawprtvol,idum,-1.d0,2)
             else
               call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,&
&               1,-1,idum,0,pawprtvol,idum,-1.d0,2,opt_sym=1)
             end if
           else
             if (ipert==0) then
               call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,&
&               1,-1,idum,0,pawprtvol,idum,-1.d0,2,asym_ij=dijsym)
             else
               call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,&
&               1,-1,idum,0,pawprtvol,idum,-1.d0,2,opt_sym=1)
             end if
           end if
         end if
       end if
     end if

     if (abs(pawprtvol)<1) then
       if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
         if (nspden==2.and.nsppol==1) then
           write(msg, '(4a,i6,3a)') ch10,&
&           ' ******  TOTAL Dij', pertstrg, ' in Ha (atom ', iatom_tot,') *****',ch10,&
&           ' (antiferromagnetism case: only one spin component)'
         else if (paw_ij(iatom)%ndij==1) then
           write(msg, '(4a,i6,a)') ch10,&
&           ' ******  TOTAL Dij',pertstrg,' in Ha (atom ',iatom_tot,') *****'
         else
           write(msg, '(4a,i6,3a)') ch10,&
&           ' ******  TOTAL Dij',pertstrg,' in Ha (atom ',iatom_tot,&
&           ', component ',trim(dspin(idij+2*(nsploop/4))),') *****'
         end if
         call wrtout(std_out,msg,'COLL')
         if (idij<=2.or.pawspnorb==0) then
           if (ipert==0) then
             call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,1,-1,&
&             idum,0,pawprtvol,idum,50.d0*dble(3-2*idij),1)
           else
             call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,1,-1,&
&             idum,0,pawprtvol,idum,50.d0*dble(3-2*idij),1,opt_sym=1)
           end if
         else
           if (ipert==0) then
             call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,1,-1,&
&             idum,0,pawprtvol,idum,50.d0*dble(3-2*idij),1,asym_ij=dijsym)
           else
             call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,1,-1,&
&             idum,0,pawprtvol,idum,50.d0*dble(3-2*idij),1,opt_sym=1)
           end if
         end if
       end if
     end if

     if (idij>=3.and.pawspnorb>0.and.ipert<=0)  then
       ABI_DEALLOCATE(dijsym)
     end if

!    ----- End loops over iatom and idij
   end do

   paw_ij(iatom)%has_dij=2 ! Flag: total Dij is computed

   if (ipert/=natom_tot+1.and.ipert/=natom_tot+5) then
     ABI_DEALLOCATE(indklmn)
     ABI_DEALLOCATE(ff)
     ABI_DEALLOCATE(dijhat)
     ABI_DEALLOCATE(dijxc)
     if (paw_ij(1)%cplex==2)  then
       ABI_DEALLOCATE(gg)
     end if
     if (need_dijxc==1)  then
       ABI_DEALLOCATE(dijxc_hat)
     end if
     if (need_dijxc_val==1)  then
       ABI_DEALLOCATE(dijxc_val)
     end if
     if (pawspnorb>0.and.paw_ij(iatom)%has_dijso==0.and.ipert<=0.and.ipositron/=1)  then
       ABI_DEALLOCATE(paw_ij(iatom)%dijso)
     end if
     if (pawtab(itypat)%usepawu>0.and.ipert<=0.and.ipositron/=1)  then
       ABI_DEALLOCATE(dijpawu)
     end if
     if (pawtab(itypat)%useexexch>0.and.ipert<=0.and.ipositron/=1)  then
       ABI_DEALLOCATE(dijexxc)
     end if
     if (pawtab(itypat)%usepawu>0.and.paw_ij(iatom)%ndij==4.and.ipert<=0)  then
       ABI_DEALLOCATE(dijsymU)
     end if
     if (pawfgrtab(iatom)%gylm_allocated==2) then
       ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
       ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(0,0))
       pawfgrtab(iatom)%gylm_allocated=0
     end if
     if (pawfgrtab(iatom)%gylmgr_allocated==2) then
       ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
       ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(0,0,0))
       pawfgrtab(iatom)%gylmgr_allocated=0
     end if
     if (pawfgrtab(iatom)%expiqr_allocated==2) then
       ABI_DEALLOCATE(pawfgrtab(iatom)%expiqr)
       ABI_ALLOCATE(pawfgrtab(iatom)%expiqr,(0,0))
       pawfgrtab(iatom)%expiqr_allocated=0
     end if
   end if

 end do

 if (pawxcdev==0.and.ipert/=natom_tot+1.and.ipert/=natom_tot+5)  then
   ABI_DEALLOCATE(yylmr)
 end if

 call timab(561,2,tsec)

!----- Reduction of dij in case of parallelism
!! TO BE REMOVED
 if (paral_kgb==1.and.mpi_enreg%nproc_band>1.and.mpi_enreg%nproc_atom==1) then
   bufdim=sum(paw_ij(:)%cplex_dij*paw_ij(:)%lmn2_size)*nsploop &
&   *(1+use_dijxc+use_dijxc_val+use_dijso+use_dijU+use_dijhat)
   ABI_ALLOCATE(buffer,(bufdim))
   buffer(:)=zero
   ic=1
   do iatom=1,natom
     shift=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
     do idij=1,nsploop
       buffer(ic:ic+shift-1)=paw_ij(iatom)%dij(:,idij)
       if (use_dijxc==1) buffer(ic+shift:ic+(1+use_dijxc)*shift-1)=paw_ij(iatom)%dijxc(:,idij)
       icinf=ic+(1+use_dijxc)*shift
       icsup=ic+(1+use_dijxc+use_dijxc_val)*shift-1
       if (use_dijxc_val==1) buffer(icinf:icsup)=paw_ij(iatom)%dijxc_val(:,idij)
       icinf=icsup+1
       icsup=ic+(1+use_dijxc+use_dijxc_val+use_dijso)*shift-1
       if (use_dijso==1) buffer(icinf:icsup)=paw_ij(iatom)%dijso(:,idij)
       icinf=icsup+1
       icsup=ic+(1+use_dijxc+use_dijxc_val+use_dijso+use_dijU)*shift-1
       if (use_dijU==1)  buffer(icinf:icsup)=paw_ij(iatom)%dijU(:,idij)
       icinf=icsup+1
       icsup=ic+(1+use_dijxc+use_dijxc_val+use_dijso+use_dijU+use_dijhat)*shift-1
       if (use_dijhat==1) buffer(icinf:icsup)=paw_ij(iatom)%dijhat(:,idij)
       ic=ic+(1+use_dijxc+use_dijxc_val+use_dijso+use_dijU+use_dijhat)*shift
     end do
   end do
   call xsum_mpi(buffer,mpi_enreg%comm_band,ier)
   ic=1
   do iatom=1,natom
     shift=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
     do idij=1,nsploop
       paw_ij(iatom)%dij(:,idij)=buffer(ic:ic+shift-1)
       icsup = ic+(1+use_dijxc)*shift-1
       if (use_dijxc==1)paw_ij(iatom)%dijxc(:,idij)=buffer(ic:icsup)
       icinf=icsup+1
       icsup=ic+(1+use_dijxc+use_dijxc_val)*shift-1
       if (use_dijxc_val==1) paw_ij(iatom)%dijxc_val(:,idij)=buffer(icinf:icsup)
       icinf=icsup+1
       icsup=ic+(1+use_dijxc+use_dijxc_val+use_dijso)*shift-1
       if (use_dijso==1) paw_ij(iatom)%dijso(:,idij)=buffer(icinf:icsup)
       icinf=icsup+1
       icsup=ic+(1+use_dijxc+use_dijxc_val+use_dijso+use_dijU)*shift-1
       if (use_dijU==1) paw_ij(iatom)%dijU(:,idij)=buffer(icinf:icsup)
       icinf=icsup+1
       icsup=ic+(1+use_dijxc+use_dijxc_val+use_dijso+use_dijU+use_dijhat)*shift-1
       if (use_dijhat==1) paw_ij(iatom)%dijhat(:,idij)=buffer(icinf:icsup)
       ic=ic+(1+use_dijxc+use_dijxc_val+use_dijso+use_dijU+use_dijhat)*shift
     end do
   end do
   ABI_DEALLOCATE(buffer)
 end if

 DBG_EXIT("COLL")

end subroutine pawdij
!!***
