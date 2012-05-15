!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawmknhat
!! NAME
!! pawmknhat
!!
!! FUNCTION
!! PAW only:
!! Compute compensation charge density (and derivatives) on the fine FFT grid
!! Can also compute first-order compensation charge density (RF calculations)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!  ider= 0: nhat(r) is computed
!!        1: cartesian derivatives of nhat(r) are computed
!!        2: nhat(r) and derivatives are computed
!!        Note: ider>0 not compatible with ipert>0
!!  idir=direction of atomic displacement (in case of atomic displ. perturb.)
!!  ipert=index of perturbation; must be 0 for ground-state calculations
!!  izero=if 1, unbalanced components of nhat(g) have to be set to zero
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms on current process, size of PAW arrays
!!  natom_tot=total number of atoms in cell
!!  nfft=number of point on the rectangular fft grid
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhatgrdim= -PAW only- 0 if pawgrnhat array is not used ; 1 otherwise
!!  ntypat=number of types of atoms in unit cell.
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                         (1st-order occupancies if ipert>0)
!!  pawrhoij0(natom) <type(pawrhoij_type)>= GS paw rhoij occupancies and related data (used only if ipert>0)
!!                                          set equat to pawrhoij for GS calculations
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  qphon(3)=wavevector of the phonon (RF only)
!!  rprimd(3,3)=dimensional primitive translations for real space
!!  ucvol=volume of the unit cell
!!  xred(3,natom)= reduced atomic coordinates
!!
!! OUTPUT
!!  === if ider=0 or 2
!!    compch_fft=compensation charge inside spheres computed over fine fft grid
!!    pawnhat(nfft,ispden)=nhat on fine rectangular grid
!!  === if ider=1 or 2
!!    pawgrnhat(nfft,ispden,3)=derivatives of nhat on fine rectangular grid (and derivatives)
!!
!! PARENTS
!!      afterscfloop,bethe_salpeter,energy,nres2vres,odamix,paw_qpscgw,pawmkrho
!!      respfn,scfcv,scfcv3,screening,setup_positron,sigma
!!
!! CHILDREN
!!      fourdp,mean_fftr,pawexpiqr,pawfrnhat,pawgylm,timab,xsum_mpi,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,mpi_enreg,&
&          natom,natom_tot,nfft,ngfft,nhatgrdim,nspden,ntypat,paral_kgb,pawang,pawfgrtab,&
&          pawgrnhat,pawnhat,pawrhoij,pawrhoij0,pawtab,qphon,rprimd,ucvol,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawmknhat'
 use interfaces_18_timing
 use interfaces_53_ffts
 use interfaces_53_spacepar
 use interfaces_66_paw, except_this_one => pawmknhat
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,ider,idir,ipert,izero,natom,natom_tot,nfft,nhatgrdim,nspden,ntypat,paral_kgb
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: compch_fft
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),qphon(3),rprimd(3,3),xred(3,natom_tot)
 real(dp),intent(out) :: pawgrnhat(cplex*nfft,nspden,3*nhatgrdim)
 real(dp),intent(out) :: pawnhat(cplex*nfft,nspden)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom),pawrhoij0(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: dplex,iatom,iatom_tot,ic,ierr,ils,ilslm,irhoij,ispden,itypat,jc,jrhoij,kc,klm,klmn
 integer :: lmax,lmin,lm_size,mm,nfftot,nfgd,optfr,optgr0,optgr1
 logical :: compute_grad,compute_nhat,need_frozen,qeq0
!arrays
 real(dp) :: rdum(1),ro(cplex),ro_ql(cplex),tmp_compch_fft(nspden),tsec(2),vdum(1,1)
 real(dp),allocatable :: work(:,:)
 type(pawrad_type),allocatable :: pawrad_dum(:)
 type(paw_ij_type),allocatable :: paw_ij_dum(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 if(ider>0.and.nhatgrdim==0) then
   MSG_BUG('  Gradients of nhat required but not allocated !')
 end if
 if(ider>0.and.ipert/=0) then
   MSG_BUG('  Gradients of nhat not compatible with RF (ipert>0) !')
 end if
 if(ider>0.and.cplex==2) then
   MSG_BUG('  Gradients of nhat not compatible with cplex=2 !')
 end if
 if(nspden>1.and.nspden/=pawrhoij(1)%nspden) then
   MSG_BUG('  Wrong values for nspden and pawrhoij%nspden !')
 end if
 if(nspden>1.and.nspden/=pawfgrtab(1)%nspden) then
   MSG_BUG('  Wrong values for nspden and pawfgrtab%nspden !')
 end if
 if(pawrhoij(1)%cplex<cplex) then
   MSG_BUG('  Must have pawrhoij()%cplex >= cplex !')
 end if
 if (mpi_enreg%nproc_atom>1) then
   if (natom/=mpi_enreg%natom) then
     MSG_BUG("natom not equal to mpi_enreg%natom !")
   end if
 end if
 qeq0=(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15)
 if (ipert>0.and.ipert<=natom_tot.and.pawfgrtab(1)%rfgd_allocated==0.and.(.not.qeq0)) then
   MSG_BUG('  pawfgrtab()%rfgd array must be allocated  !')
 end if

!nhat1 does not have to be computed for ddk or magn. field perturbation
 if (ipert==natom_tot+1.or.ipert==natom_tot+5) return

!Initialisations
 compute_nhat=(ider==0.or.ider==2)
 compute_grad=(ider==1.or.ider==2)
 if ((.not.compute_nhat).and.(.not.compute_grad)) return
 dplex=cplex-1

 if (compute_nhat) pawnhat=zero
 if (compute_grad) pawgrnhat=zero

!------------------------------------------------------------------------
!----- Loop over atoms
!------------------------------------------------------------------------
 do iatom=1,natom
   iatom_tot=iatom;if (mpi_enreg%nproc_atom>1) iatom_tot=mpi_enreg%atom_indx(iatom)
   if (paral_kgb==1.and.mpi_enreg%nproc_atom==1.and.mod(iatom-1,mpi_enreg%nproc_band)/=mpi_enreg%me_band) cycle
   itypat=pawrhoij(iatom)%itypat
   lm_size=pawfgrtab(iatom)%l_size**2
   need_frozen=(ipert==iatom_tot)
   nfgd=pawfgrtab(iatom)%nfgd

!  Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
   if (((compute_nhat).and.(pawfgrtab(iatom)%gylm_allocated==0)).or.&
&   ((compute_grad).and.(pawfgrtab(iatom)%gylmgr_allocated==0))) then
     optgr0=0;optgr1=0
     if ((compute_nhat).and.(pawfgrtab(iatom)%gylm_allocated==0)) then
       if (associated(pawfgrtab(iatom)%gylm))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(nfgd,pawfgrtab(iatom)%l_size**2))
       pawfgrtab(iatom)%gylm_allocated=2;optgr0=1
     end if
     if ((compute_grad).and.(pawfgrtab(iatom)%gylmgr_allocated==0)) then
       if (associated(pawfgrtab(iatom)%gylmgr))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(3,nfgd,pawfgrtab(iatom)%l_size**2))
       pawfgrtab(iatom)%gylmgr_allocated=2;optgr1=1
     end if
     call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,rdum,&
&     lm_size,nfgd,optgr0,optgr1,0,pawtab(itypat),&
&     pawfgrtab(iatom)%rfgd,pawfgrtab(iatom)%rfgd_allocated)
   end if

!  Eventually compute exp(-i.q.r) factors for the current atom (if not already done)
   if (ipert>0.and.ipert<=natom_tot.and.pawfgrtab(1)%rfgd_allocated==0.and.(.not.qeq0)) then
     if (associated(pawfgrtab(iatom)%expiqr))  then
       ABI_DEALLOCATE(pawfgrtab(iatom)%expiqr)
     end if
     ABI_ALLOCATE(pawfgrtab(iatom)%expiqr,(2,nfgd))
     pawfgrtab(iatom)%expiqr_allocated=2
     call pawexpiqr(gprimd,pawfgrtab(iatom),qphon,xred(:,iatom_tot))
   end if

!  Eventually compute frozen part of nhat for the current atom (if not already done)
   if ((need_frozen).and.(pawfgrtab(iatom)%nhatfr_allocated==0)) then
     if (associated(pawfgrtab(iatom)%nhatfr))  then
       ABI_DEALLOCATE(pawfgrtab(iatom)%nhatfr)
     end if
     ABI_ALLOCATE(pawfgrtab(iatom)%nhatfr,(nfgd,pawfgrtab(iatom)%nspden))
     optfr=0;pawfgrtab(iatom)%nhatfr_allocated=2
     ABI_ALLOCATE(paw_ij_dum,(natom))
     ABI_ALLOCATE(pawrad_dum,(ntypat))
     call pawfrnhat(cplex,gprimd,idir,ipert,mpi_enreg,natom,natom_tot,nfft,ngfft,nspden,ntypat,&
&     optfr,paw_ij_dum,pawang,pawfgrtab,pawrad_dum,pawrhoij0,pawtab,qphon,rprimd,&
&     ucvol,vdum,vdum,xred)
     ABI_DEALLOCATE(paw_ij_dum)
     ABI_DEALLOCATE(pawrad_dum)
   end if

!  ------------------------------------------------------------------------
!  ----- Loop over density components
!  ------------------------------------------------------------------------

   do ispden=1,nspden

!    ------------------------------------------------------------------------
!    ----- Loop over ij channels (basis components)
!    ------------------------------------------------------------------------
     jrhoij=1
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
       klmn=pawrhoij(iatom)%rhoijselect(irhoij)
       klm =pawtab(itypat)%indklmn(1,klmn)
       lmin=pawtab(itypat)%indklmn(3,klmn)
       lmax=pawtab(itypat)%indklmn(4,klmn)

!      Retrieve rhoij
       if (nspden/=2) then
         ro(1:cplex)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
       else
         if (ispden==1) then
           ro(1:cplex)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,1)&
&           +pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,2)
         else if (ispden==2) then
           ro(1:cplex)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,1)
         end if
       end if
       ro(1:cplex)=pawtab(itypat)%dltij(klmn)*ro(1:cplex)

       if (compute_nhat) then
         if (cplex==1) then
           do ils=lmin,lmax,2
             do mm=-ils,ils
               ilslm=ils*ils+ils+mm+1
               if (pawang%gntselect(ilslm,klm)>0) then
                 ro_ql(1)=ro(1)*pawtab(itypat)%qijl(ilslm,klmn)
                 do ic=1,nfgd
                   jc=pawfgrtab(iatom)%ifftsph(ic)
                   pawnhat(jc,ispden)=pawnhat(jc,ispden)+ro_ql(1)*pawfgrtab(iatom)%gylm(ic,ilslm)
                 end do
               end if
             end do
           end do
         else
           do ils=lmin,lmax,2
             do mm=-ils,ils
               ilslm=ils*ils+ils+mm+1
               if (pawang%gntselect(ilslm,klm)>0) then
                 ro_ql(1:2)=ro(1:2)*pawtab(itypat)%qijl(ilslm,klmn)
                 do ic=1,nfgd
                   jc=2*pawfgrtab(iatom)%ifftsph(ic)-1
                   pawnhat(jc:jc+1,ispden)=pawnhat(jc:jc+1,ispden)+ro_ql(1:2)*pawfgrtab(iatom)%gylm(ic,ilslm)
                 end do
               end if
             end do
           end do
         end if
       end if

       if (compute_grad) then
!        Not compatible with cplex=2
         do ils=lmin,lmax,2
           do mm=-ils,ils
             ilslm=ils*ils+ils+mm+1
             if (pawang%gntselect(ilslm,klm)>0) then
               ro_ql(1)=ro(1)*pawtab(itypat)%qijl(ilslm,klmn)
               do ic=1,nfgd
                 jc=pawfgrtab(iatom)%ifftsph(ic)
                 pawgrnhat(jc,ispden,1)=pawgrnhat(jc,ispden,1)+ro_ql(1)*pawfgrtab(iatom)%gylmgr(1,ic,ilslm)
                 pawgrnhat(jc,ispden,2)=pawgrnhat(jc,ispden,2)+ro_ql(1)*pawfgrtab(iatom)%gylmgr(2,ic,ilslm)
                 pawgrnhat(jc,ispden,3)=pawgrnhat(jc,ispden,3)+ro_ql(1)*pawfgrtab(iatom)%gylmgr(3,ic,ilslm)
               end do
             end if
           end do
         end do
       end if

!      ------------------------------------------------------------------------
!      ----- End loop over ij channels
!      ------------------------------------------------------------------------
       jrhoij=jrhoij+pawrhoij(iatom)%cplex
     end do

!    If RF calculation, add eventually frozen part of 1st-order compensation density
     if (need_frozen) then
       if (cplex==1) then
         do ic=1,nfgd
           jc=pawfgrtab(iatom)%ifftsph(ic)
           pawnhat(jc,ispden)=pawnhat(jc,ispden)+pawfgrtab(iatom)%nhatfr(ic,ispden)
         end do
       else
         do ic=1,nfgd
           jc=2*pawfgrtab(iatom)%ifftsph(ic)-1
           pawnhat(jc,ispden)=pawnhat(jc,ispden)+pawfgrtab(iatom)%nhatfr(ic,ispden)
         end do
       end if
     end if

!    If needed, multiply eventually by exp(-i.q.r) phase
     if (ipert>0.and.ipert<=natom_tot.and.pawfgrtab(iatom)%expiqr_allocated/=0) then
       if (cplex==1) then
         do ic=1,nfgd
           jc=pawfgrtab(iatom)%ifftsph(ic)
           pawnhat(jc,ispden)=pawnhat(jc,ispden)*pawfgrtab(iatom)%expiqr(1,ic)
         end do
       else
         do ic=1,nfgd
           jc=2*ic-1;kc=2*pawfgrtab(iatom)%ifftsph(ic)-1
           ro(1:2)=pawnhat(kc:kc+1,ispden)
           pawnhat(kc  ,ispden)=ro(1)*pawfgrtab(iatom)%expiqr(1,ic)+ro(2)*pawfgrtab(iatom)%expiqr(2,ic)
           pawnhat(kc+1,ispden)=ro(2)*pawfgrtab(iatom)%expiqr(1,ic)-ro(1)*pawfgrtab(iatom)%expiqr(2,ic)
         end do
       end if
     end if

!    ------------------------------------------------------------------------
!    ----- End loop over density components
!    ------------------------------------------------------------------------
   end do

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
   if (pawfgrtab(iatom)%nhatfr_allocated==2) then
     ABI_DEALLOCATE(pawfgrtab(iatom)%nhatfr)
     ABI_ALLOCATE(pawfgrtab(iatom)%nhatfr,(0,0))
     pawfgrtab(iatom)%nhatfr_allocated=0
   end if
   if (pawfgrtab(iatom)%expiqr_allocated==2) then
     ABI_DEALLOCATE(pawfgrtab(iatom)%expiqr)
     ABI_ALLOCATE(pawfgrtab(iatom)%expiqr,(0,0))
     pawfgrtab(iatom)%expiqr_allocated=0
   end if

!  ------------------------------------------------------------------------
!  ----- End loop over atoms
!  ------------------------------------------------------------------------
 end do

!----- Reduction in case of parallelism
!! TO BE REMOVED
 if (paral_kgb==1.and.mpi_enreg%nproc_band>1.and.mpi_enreg%nproc_atom==1)then
   call timab(48,1,tsec)
   if (compute_nhat) then
     call xsum_mpi(pawnhat,mpi_enreg%comm_band,ierr)
   end if
   if (compute_grad) then
     call xsum_mpi(pawgrnhat,mpi_enreg%comm_band,ierr)
   end if
   call timab(48,2,tsec)
 end if

 if (mpi_enreg%nproc_atom>1)then
   call timab(48,1,tsec)
   if (compute_nhat) then
     call xsum_mpi(pawnhat,mpi_enreg%comm_atom,ierr)
   end if
   if (compute_grad) then
     call xsum_mpi(pawgrnhat,mpi_enreg%comm_atom,ierr)
   end if
   call timab(48,2,tsec)
 end if

!----- Avoid unbalanced g-components numerical errors
 if (izero==1.and.compute_nhat) then
   ABI_ALLOCATE(work,(2,nfft))
   do ispden=1,min(2,nspden)
     call fourdp(cplex,work,pawnhat(:,ispden),-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     call zerosym(work,2,mpi_enreg,ngfft(1),ngfft(2),ngfft(3))
     call fourdp(cplex,work,pawnhat(:,ispden),+1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   end do
   ABI_DEALLOCATE(work)
 end if

!----- Computation of compensation charge over FFT grid
 if (compute_nhat.and.ipert==0) then
   nfftot=ngfft(1)*ngfft(2)*ngfft(3)
   call mean_fftr(pawnhat,tmp_compch_fft,mpi_enreg,nfft,nfftot,1)
   compch_fft = tmp_compch_fft(1)
   compch_fft=compch_fft*ucvol
 end if

 DBG_EXIT("COLL")

end subroutine pawmknhat
!!***
