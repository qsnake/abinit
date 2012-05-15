!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawfrnhat
!!
!! NAME
!! pawfrnhat
!!
!! FUNCTION
!! PAW: Compute frozen part of charge compensation density nhat and
!!      frozen part of psp strength Dij due to 1st-order compensation density
!!      nhatfr(r)=Sum_ij,lm[rhoij_ij.q_ij^l.(g_l(r).Y_lm(r))^(1)]
!!      Dijfr    =Int_R^3{vtrial*Sum_LM[Q_ij_q^LM^(1)] + Vloc^(1)*Sum_LM[Q_ij_q^LM]}
!!      Depend on q wave vector but not on first-order wave-function.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL; if 2, COMPLEX
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!  idir=direction of atomic displacement (in case of phonons perturb.)
!!  ipert=nindex of perturbation
!!  natom=number of atoms on current process, size of PAW arrays
!!  natom_tot=total number of atoms in cell
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms
!!  optfr=0: computes only frozen part of compensation density
!!        1: computes only frozen part of Dij
!!        2: computes both
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= Ground-State paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  qphon(3)=wavevector of the phonon
!!  rprimd(3,3)=dimensional primitive translations for real space
!!  ucvol=unit cell volume (bohr^3)
!!  vpsp1(cplex*nfft)= first-order change of local potential
!!  vtrial(nfft,nspden)= total GS potential
!!  xred(3,natom)= reduced atomic coordinates
!!
!! OUTPUT
!!  === If optfr= 0 or 2
!!    pawfgrtab(iatom)%nhatfr(nfgd,nspden)
!!                    frozen part of charge compensation density (inside PAW spheres)
!!                    =Sum_ij,lm[rhoij_ij.q_ij^l.(g_l(r).Y_lm(r))^(1)]
!!  === If optfr= 0 or 1
!!    paw_ij1(iatom)%dijfr(cplex_dij*lmn2_size,nspden)=
!!                    frozen contribution to psp strength Dij
!!                    =Int_R^3{vtrial*Sum_LM[Q_ij_q^LM^(1)] + Vloc^(1)*Sum_LM[Q_ij_q^LM]}
!!
!! PARENTS
!!      dyfnl3,nstpaw3,pawmknhat,scfcv3
!!
!! CHILDREN
!!      deducer0,pawexpiqr,pawgylm,simp_gen,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawfrnhat(cplex,gprimd,idir,ipert,mpi_enreg,natom,natom_tot,nfft,ngfft,nspden,ntypat,&
&          optfr,paw_ij1,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,qphon,rprimd,&
&          ucvol,vpsp1,vtrial,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors

 use m_radmesh,          only : simp_gen, deducer0

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawfrnhat'
 use interfaces_51_manage_mpi
 use interfaces_66_paw, except_this_one => pawfrnhat
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,natom,natom_tot,nfft,nspden,ntypat,optfr
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),qphon(3),rprimd(3,3)
 real(dp),intent(in) :: vpsp1(cplex*nfft*(optfr+1)/2),vtrial(nfft,nspden*(optfr+1)/2)
 real(dp),intent(in) :: xred(3,natom_tot)
 type(paw_ij_type),intent(inout) :: paw_ij1(natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: dplex,iatom,iatom_tot,ic,ier,ils,ilslm,irhoij,isel,ispden,itypat,jc,jrhoij
 integer :: klm,klmn,klmn1,kln,lm_size,lmn2_size,lm0,lmax,lmin,mesh_size,mm,mu
 integer :: nfftot,nfgd,old_paral_level,optgr0,optgr1,spaceComm
 logical :: has_phase,need_dijfr_1,need_dijfr_2,need_dijfr_3,need_nhatfr,qne0
 real(dp) :: c1,fact,intg,rg1,ro
!arrays
 integer,parameter :: m_index(3)=(/1,-1,0/)
 real(dp) :: contrib(2),rdum(1)
 real(dp),allocatable :: ff(:),intv(:,:),intvloc(:,:),intv_tmp(:,:)
 real(dp),allocatable :: nhatfr_tmp(:,:),rg(:),vloc(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 qne0=(qphon(1)**2+qphon(2)**2+qphon(3)**2>=1.d-15)
 if (ipert<=natom_tot) then
   if (pawfgrtab(1)%rfgd_allocated==0.and.paw_ij1(1)%has_dijfr==1.and.qne0) then
     MSG_BUG('  pawfgrtab()%rfgd array must be allocated  !')
   end if
 end if
 if (optfr==1.or.optfr==2) then
   if (paw_ij1(1)%cplex_dij<cplex) then
     MSG_BUG('  paw_ij1()%cplex_dij must be >=cplex !')
   end if
 end if
 if (optfr==1.or.optfr==2) then
   if (paw_ij1(1)%cplex/=cplex) then
     MSG_BUG('  paw_ij1()%cplex and cplex must be equal !')
   end if
   if ((ipert<=natom_tot.or.ipert==natom_tot+2.or.ipert==natom_tot+5).and. &
&   paw_ij1(1)%has_dijfr==0) then
     MSG_BUG('  pawdij1()%dijfr must be allocated !')
   end if
 end if
 if (mpi_enreg%nproc_atom>1) then
   if (natom/=mpi_enreg%natom) then
     MSG_BUG("natom not equal to mpi_enreg%natom !")
   end if
 end if

!Nothing to be done for DDK or strain
 if (ipert==natom_tot+1.or.ipert==natom_tot+6 .or. &
& ipert==natom_tot+3.or.ipert==natom_tot+4) return

 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 fact=ucvol/dble(nfftot)
 dplex=cplex-1

!Loops over  atoms
 do iatom=1,natom
   iatom_tot=iatom;if (mpi_enreg%nproc_atom>1) iatom_tot=mpi_enreg%atom_indx(iatom)

!  Eventually allocate frozen nhat points
   if ((optfr==0.or.optfr==2).and.(ipert==iatom_tot)) then
     if (associated(pawfgrtab(iatom)%nhatfr))  then
       ABI_DEALLOCATE(pawfgrtab(iatom)%nhatfr)
     end if
     ABI_ALLOCATE(pawfgrtab(iatom)%nhatfr,(pawfgrtab(iatom)%nfgd,nspden))
     pawfgrtab(iatom)%nhatfr_allocated=1
   end if

!  Select which atom to treat
!  need_nhatfr=((ipert==iatom_tot).and.(pawfgrtab(iatom)%nhatfr_allocated==1))
   need_nhatfr=(ipert==iatom_tot);if (need_nhatfr) need_nhatfr=(pawfgrtab(iatom)%nhatfr_allocated==1)
   need_dijfr_1=(ipert==iatom_tot.and.(paw_ij1(iatom)%has_dijfr==1))
   need_dijfr_2=((ipert<=natom_tot).and.(paw_ij1(iatom)%has_dijfr==1))
   need_dijfr_3=((ipert==natom_tot+2.or.ipert==natom_tot+5).and.(paw_ij1(iatom)%has_dijfr==1))
   if ((.not.need_nhatfr).and.(.not.need_dijfr_1).and.(.not.need_dijfr_2).and.(.not.need_dijfr_3)) cycle

!  Some atom-dependent quantities
   itypat=pawrhoij(iatom)%itypat
   lm_size=pawtab(itypat)%l_size**2
   lmn2_size=pawtab(itypat)%lmn2_size
   mesh_size=pawrad(itypat)%mesh_size

!  Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
   nfgd=0
   if (need_nhatfr.or.need_dijfr_1.or.need_dijfr_2) then
     nfgd=pawfgrtab(iatom)%nfgd
     if (((need_dijfr_2).and.(pawfgrtab(iatom)%gylm_allocated==0)).or.&
&     ((need_dijfr_1.or.need_nhatfr).and.(pawfgrtab(iatom)%gylmgr_allocated==0))) then
       optgr0=0;optgr1=0
       if ((need_dijfr_2).and.(pawfgrtab(iatom)%gylm_allocated==0)) then
         if (associated(pawfgrtab(iatom)%gylm))  then
           ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
         end if
         ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(nfgd,lm_size))
         pawfgrtab(iatom)%gylm_allocated=2;optgr0=1
       end if
       if ((need_dijfr_1.or.need_nhatfr).and.(pawfgrtab(iatom)%gylmgr_allocated==0)) then
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
   end if

!  Eventually compute exp(-i.q.r) factors for the current atom (if not already done)
   has_phase=.false.
   if (need_dijfr_2) then
     if (qne0.and.(pawfgrtab(iatom)%expiqr_allocated==0)) then
       if (associated(pawfgrtab(iatom)%expiqr))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%expiqr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%expiqr,(2,nfgd))
       pawfgrtab(iatom)%expiqr_allocated=2
       call pawexpiqr(gprimd,pawfgrtab(iatom),qphon,xred(:,iatom_tot))
     end if
     has_phase=(pawfgrtab(iatom)%expiqr_allocated/=0)
   end if

!  Loop over spin components
   do ispden=1,nspden

!    Computation of frozen part of 1st-order compensation density
!    =Sum_ij,lm[rhoij_ij.q_ij^l.(g_l(r).Y_lm(r))^(1)]
!    ------------------------------------------------------------
     if (need_nhatfr) then
       jrhoij=1
       ABI_ALLOCATE(nhatfr_tmp,(3,nfgd))
       nhatfr_tmp=zero
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(irhoij)
         klm =pawtab(itypat)%indklmn(1,klmn)
         lmin=pawtab(itypat)%indklmn(3,klmn)
         lmax=pawtab(itypat)%indklmn(4,klmn)
         if (nspden/=2) then
           ro=pawrhoij(iatom)%rhoijp(jrhoij,ispden)
         else
           if (ispden==1) then
             ro=pawrhoij(iatom)%rhoijp(jrhoij,1)+pawrhoij(iatom)%rhoijp(jrhoij,2)
           else if (ispden==2) then
             ro=pawrhoij(iatom)%rhoijp(jrhoij,1)
           end if
         end if
         ro=pawtab(itypat)%dltij(klmn)*ro
         do ils=lmin,lmax,2
           lm0=ils**2+ils+1
           do mm=-ils,ils
             ilslm=lm0+mm;isel=pawang%gntselect(lm0+mm,klm)
             if (isel>0) then
               do ic=1,nfgd
                 do mu=1,3
                   contrib(1)=-ro*pawtab(itypat)%qijl(ilslm,klmn)*pawfgrtab(iatom)%gylmgr(mu,ic,ilslm)
                   nhatfr_tmp(mu,ic)=nhatfr_tmp(mu,ic)+contrib(1)
                 end do
               end do
             end if
           end do
         end do
         jrhoij=jrhoij+pawrhoij(iatom)%cplex
       end do

!      Convert from cartesian to reduced coordinates
       do ic=1,nfgd
         pawfgrtab(iatom)%nhatfr(ic,ispden)= &
&         rprimd(1,idir)*nhatfr_tmp(1,ic) &
&         +rprimd(2,idir)*nhatfr_tmp(2,ic) &
&         +rprimd(3,idir)*nhatfr_tmp(3,ic)
       end do
       ABI_DEALLOCATE(nhatfr_tmp)

     end if

!    Computation of frozen part of 1st-order psps strength Dij
!    ------------------------------------------------------------------

!    ============ Phonons ====================================
     if (ipert<=natom_tot) then

       if (need_dijfr_1.or.need_dijfr_2) then

         ABI_ALLOCATE(intv,(cplex,lm_size))
         intv=zero

!        First part: Int_R^3{vtrial*Sum_LM[Q_ij_q^LM^(1)]}
         if (need_dijfr_1) then

!          ----- Retrieve potential Vtrial (subtle if nspden=4 ;-)
           if (nspden/=4) then
             ABI_ALLOCATE(vloc,(1,nfgd))
             do ic=1,nfgd
               vloc(1,ic)=vtrial(pawfgrtab(iatom)%ifftsph(ic),ispden)
             end do
           else
             ABI_ALLOCATE(vloc,(2,nfgd))
             if (ispden<=2) then
               do ic=1,nfgd
                 jc=pawfgrtab(iatom)%ifftsph(ic)
                 vloc(1,ic)=vtrial(jc,ispden)
                 vloc(2,ic)=zero
               end do
             else if (ispden==3) then
               do ic=1,nfgd
                 jc=pawfgrtab(iatom)%ifftsph(ic)
                 vloc(1,ic)=vtrial(jc,3)
                 vloc(2,ic)=vtrial(jc,4)
               end do
             else ! ispden=4
               vloc(2,1:nfgd)=-vloc(2,1:nfgd)
             end if
           end if

!          ----- Compute Integral [ Vtrial(r).(g_l(r).Y_lm(r))^(1) dr ]
           ABI_ALLOCATE(intv_tmp,(cplex,3))
           do ilslm=1,lm_size
             intv_tmp=zero
             if (nspden/=4) then
               do ic=1,nfgd
                 do mu=1,3
!                  Minus sign because dg(r-R)/dR = -dg(r-R)/dr
                   contrib(1)=-vloc(1,ic)*pawfgrtab(iatom)%gylmgr(mu,ic,ilslm)
                   intv_tmp(1,mu)=intv_tmp(1,mu)+contrib(1)
                 end do
               end do
             else ! nspden=4
               do ic=1,nfgd
                 do mu=1,3
!                  Minus sign because dg(r-R)/dR = -dg(r-R)/dr
                   contrib(1:2)=-vloc(1:2,ic)*pawfgrtab(iatom)%gylmgr(mu,ic,ilslm)
                   intv_tmp(1:2,mu)=intv_tmp(1:2,mu)+contrib(1:2)
                 end do
               end do
             end if
!            Convert from cartesian to reduced coordinates
             intv(1:cplex,ilslm)=intv(1:cplex,ilslm) &
&             +(rprimd(1,idir)*intv_tmp(1:cplex,1) &
&             +rprimd(2,idir)*intv_tmp(1:cplex,2) &
&             +rprimd(3,idir)*intv_tmp(1:cplex,3))
           end do
           ABI_DEALLOCATE(vloc)
           ABI_DEALLOCATE(intv_tmp)
         end if ! need_dijfr_1

!        2nd part: Int_R^3{Vloc^(1)*Sum_LM[Q_ij_q^LM]}
         if (need_dijfr_2) then

           if (ispden==1) then

!            ----- Retrieve potential Vloc^(1)
             ABI_ALLOCATE(vloc,(cplex,nfgd))
             do ic=1,nfgd
               jc=cplex*pawfgrtab(iatom)%ifftsph(ic)-dplex
               vloc(1:cplex,ic)=vpsp1(jc:jc+dplex)
             end do

!            ----- Compute Integral [ Vloc^(1)(r).g_l(r).Y_lm(r) ]
             ABI_ALLOCATE(intvloc,(cplex,lm_size))
             intvloc=zero
             if (has_phase) then
               if (cplex==1) then
                 do ilslm=1,lm_size
                   do ic=1,nfgd
                     contrib(1)=vloc(1,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
                     intvloc(1,ilslm)=intvloc(1,ilslm)+contrib(1)*pawfgrtab(iatom)%expiqr(1,ic)
                   end do
                 end do
               else
                 do ilslm=1,lm_size
                   do ic=1,nfgd
                     contrib(1:cplex)=vloc(1:cplex,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
                     intvloc(1,ilslm)=intvloc(1,ilslm)+contrib(1)*pawfgrtab(iatom)%expiqr(1,ic) &
&                     -contrib(2)*pawfgrtab(iatom)%expiqr(2,ic)
                     intvloc(2,ilslm)=intvloc(2,ilslm)+contrib(1)*pawfgrtab(iatom)%expiqr(2,ic) &
&                     +contrib(2)*pawfgrtab(iatom)%expiqr(1,ic)
                   end do
                 end do
               end if
             else ! no phase
               do ilslm=1,lm_size
                 do ic=1,nfgd
                   contrib(1:cplex)=vloc(1:cplex,ic)*pawfgrtab(iatom)%gylm(ic,ilslm)
                   intvloc(1:cplex,ilslm)=intvloc(1:cplex,ilslm)+contrib(1:cplex)
                 end do
               end do
             end if
             ABI_DEALLOCATE(vloc)
           end if ! ispden=1

           if (ispden<=min(nspden,2)) then
             intv(1:cplex,1:lm_size)=intv(1:cplex,1:lm_size)+intvloc(1:cplex,1:lm_size)
             if (ispden==min(nspden,2))  then
               ABI_DEALLOCATE(intvloc)
             end if
           end if
         end if ! need_dijfr_2

!        Apply ucvol/nfft factor on integral
         intv(:,:)=fact*intv(:,:)

!        --- Reduction in case of parallelization ---
         if(mpi_enreg%paral_compil_fft==1)then
           old_paral_level= mpi_enreg%paral_level
           mpi_enreg%paral_level=3
           call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
           call xsum_mpi(intv,spaceComm,ier)
           mpi_enreg%paral_level=old_paral_level
         end if

         paw_ij1(iatom)%dijfr(:,ispden)=zero

!        ---- Loop over (i,j) components
         klmn1=1
         do klmn=1,lmn2_size
           klm =pawtab(itypat)%indklmn(1,klmn)
           lmin=pawtab(itypat)%indklmn(3,klmn)
           lmax=pawtab(itypat)%indklmn(4,klmn)
           do ils=lmin,lmax,2
             lm0=ils**2+ils+1
             do mm=-ils,ils
               ilslm=lm0+mm;isel=pawang%gntselect(lm0+mm,klm)
               if (isel>0) then
                 paw_ij1(iatom)%dijfr(klmn1:klmn1+dplex,ispden)= &
&                 paw_ij1(iatom)%dijfr(klmn1:klmn1+dplex,ispden) &
&                 +pawtab(itypat)%qijl(ilslm,klmn)*intv(1:cplex,ilslm)
               end if
             end do
           end do
           klmn1=klmn1+paw_ij1(iatom)%cplex_dij
         end do
         ABI_DEALLOCATE(intv)

!        Dijfr is marked as computed
         paw_ij1(iatom)%has_dijfr=2
       end if

!      ============ Electric field perturbation =======================
     else if (ipert==natom_tot+2) then

       if (need_dijfr_3) then

!        The following factor arises in expanding the angular dependence of the dipole
!        vector in terms of real spherical harmonics. The real spherical harmonics are as
!        in the routine initylmr.F90; see http://www.unioviedo.es/qcg/art/Theochem419-19-ov-BF97-rotation-matrices.pdf
         c1 = sqrt(four_pi/three)

         if (ispden==1) then

           ABI_ALLOCATE(ff,(mesh_size))
           ABI_ALLOCATE(rg,(3))

!          loop over basis state pairs for this atom
           klmn1=1
           do klmn = 1, paw_ij1(iatom)%lmn2_size
             klm =pawtab(itypat)%indklmn(1,klmn)
             kln =pawtab(itypat)%indklmn(2,klmn)
             lmin=pawtab(itypat)%indklmn(3,klmn)
             lmax=pawtab(itypat)%indklmn(4,klmn)

!            Select only l=1, because the dipole is a vector operator
             if (lmin==1) then
               lm0=3  ! (l^2+l+1) for l=1

!              Computation of <phi_i|r|phi_j>- <tphi_i|r|tphi_j>
!              the dipole vector has radial dependence r
               ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln)&
&               -pawtab(itypat)%tphitphj(2:mesh_size,kln))&
&               *pawrad(itypat)%rad(2:mesh_size)
               call deducer0(ff,mesh_size,pawrad(itypat))
               call simp_gen(intg,ff,pawrad(itypat))

!              Compute <S_li_mi|r-R|S_lj_mj>: use a real Gaunt expression (with selection rule)
               rg(1:3)=zero
               do ic=1,3
                 isel=pawang%gntselect(lm0+m_index(ic),klm)
                 if (isel>0) rg(ic)=pawang%realgnt(isel)
               end do

!              Translate from cartesian to reduced coordinates (in idir direction)
               rg1=gprimd(1,idir)*rg(1)+gprimd(2,idir)*rg(2)+gprimd(3,idir)*rg(3)

!              Build sqrt(4pi/3).<S_li_mi|r-R|S_lj_mj>.(<phi_i|r-R|phi_j>- <tphi_i|r-R|tphi_j>
               paw_ij1(iatom)%dijfr(klmn1,ispden)=c1*rg1*intg
               if (cplex==2) paw_ij1(iatom)%dijfr(klmn1+1,ispden)=zero

             else
               paw_ij1(iatom)%dijfr(klmn1,ispden)=zero
             end if ! end gaunt constraint

             klmn1=klmn1+paw_ij1(iatom)%cplex_dij
           end do ! end loop over lmn2_size pairs of basis states
           ABI_DEALLOCATE(ff)
           ABI_DEALLOCATE(rg)

!          Dijfr is spin-independent for electric field case
         else if (ispden==2) then
           paw_ij1(iatom)%dijfr(:,ispden)=paw_ij1(iatom)%dijfr(:,1)
         else
           paw_ij1(iatom)%dijfr(:,ispden)=zero
         end if

!        Dijfr is marked as computed
         paw_ij1(iatom)%has_dijfr=2
       end if

!      ============ Magnetic field perturbation =======================
     else if (ipert==natom_tot+5) then

       if (need_dijfr_3) then

!        Dijfr is marked as computed
         paw_ij1(iatom)%has_dijfr=2
       end if

     end if ! ipert

!    End loop over spin components
   end do ! ispden

!  Eventually free temporary space for g_l(r).Y_lm(r) gradients and exp(-i.q.r)
   if (need_nhatfr.or.need_dijfr_1.or.need_dijfr_2) then
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

!  End loop on atoms
 end do

 DBG_EXIT("COLL")

end subroutine pawfrnhat
!!***
