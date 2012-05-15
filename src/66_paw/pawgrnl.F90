!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawgrnl
!!
!! NAME
!! pawgrnl
!!
!! FUNCTION
!! PAW: Add to GRadients of total energy due to non-local term of Hamiltonian
!!      the contribution due to Dij derivatives
!! In particular, compute contribution to forces, stresses, dyn. matrix
!! Remember: Vnl=Sum_ij[|p_i>Dij<p_j|]
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dimnhat=second dimension of array nhat (0 or # of spin components)
!!  dimvtrial=second dimension of array vtrial (1 or # of spin components)
!!  dyfr_cplex=1 if dyfrnl is real, 2 if it is complex
!!  gsqcut_eff=Fourier cutoff on G^2 for "large sphere" of radius double that of the basis sphere
!!  mgfft=maximum size of 1D FFTs
!!  natom=number of atoms on current process, size of PAW arrays
!!  natom_tot=total number of atoms in cell
!!  nattyp(ntypat)=array describing how many atoms of each type in cell
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(nfft,dimnhat)=compensation charge density on rectangular grid in real space
!!  nspden=number of spin-density components
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms
!!  optgr= 1 if gradients with respect to atomic position(s) have to be computed
!!  optgr2= 1 if 2nd gradients with respect to atomic position(s) have to be computed
!!  optstr= 1 if gradients with respect to strain(s) have to be computed
!!  paral_kgb=flag for (kpt,FFT,bands) parallelism
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase (structure factor) information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qphon(3)=wavevector of the phonon
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  typat(natom_tot)=types of atoms
!!  vtrial(nfft,dimvtrial)= total potential
!!  xred(3,natom_tot)=reduced dimensionless atomic coordinates
!!
!! SIDE EFFECTS
!!  At input, this terms contain contribution from non-local projectors derivatives
!!  At output, they are updated with the contribution of Dij derivatives
!!  ==== if optgr=1 ====
!!   grnl(3*natom) =gradients of NL energy wrt atomic coordinates
!!  ==== if optstr=1 ====
!!   nlstr(6) =gradients of NL energy wrt strains
!!  ==== if optgr2=1 ====
!!   dyfrnl(dyfr_cplex,3,3,natom,natom) =2nd gradients of NL energy wrt atomic coordinates
!!
!! NOTES
!!   For the future parallelization over atoms:
!!     pawrhoij and pawfgrtab should be gathered...
!!
!! PARENTS
!!      dyfnl3,etotfor,forstr
!!
!! CHILDREN
!!      atm2fft3,metric,pawexpiqr,pawgylm,stresssym,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawgrnl(atindx1,dimnhat,dimvtrial,dyfrnl,dyfr_cplex,grnl,gsqcut,mgfft,mpi_enreg,natom,natom_tot,&
&                  nattyp,nfft,ngfft,nhat,nlstr,nspden,nsym,ntypat,optgr,optgr2,optstr,paral_kgb,&
&                  pawang,pawfgrtab,pawrhoij,pawtab,ph1d,psps,qphon,rprimd,symrec,typat,vtrial,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawgrnl'
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_65_psp
 use interfaces_66_paw, except_this_one => pawgrnl
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimnhat,dimvtrial,dyfr_cplex,mgfft,natom,natom_tot,nfft,nspden,nsym,ntypat
 integer,intent(in) :: optgr,optgr2,optstr,paral_kgb
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx1(natom_tot),nattyp(ntypat),ngfft(18)
 integer,intent(in) :: symrec(3,3,nsym),typat(natom_tot)
 real(dp),intent(in) :: nhat(nfft,dimnhat),ph1d(2,3*(2*mgfft+1)*natom_tot),qphon(3)
 real(dp),intent(in) :: rprimd(3,3),vtrial(nfft,dimvtrial),xred(3,natom_tot)
 real(dp),intent(inout) :: dyfrnl(dyfr_cplex,3,3,natom_tot,natom_tot*optgr2),grnl(3*natom_tot*optgr)
 real(dp),intent(inout) :: nlstr(6*optstr)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: bufind,bufsiz,cplex,iatm,iatom,iatom_tot,iatshft,ic,idiag,idir,ier,ilm,irhoij,isel,ishift_gr
 integer :: ishift_gr2,ishift_str,ispden,ispvtr,itypat,jatom,jatm,jc,jrhoij,jtypat,klm,klmn,klmn1
 integer :: ll,lm_size,lm_sizej,lmax,lmin,lmn2_size,mu,mua,mub,mushift,nfftot,nfgd
 integer :: ngrad,ngrad_nondiag,ngradp,ngradp_nondiag,ngrhat
 integer :: nsploop,old_paral_level,opt1,opt2,opt3,optv,optn,optn2,qne0,spaceComm
 logical :: has_phase
 real(dp) :: dlt_tmp,fact_ucvol,grhat_x,hatstr_diag,ro,ro_d,ucvol
 character(len=500) :: msg
!arrays
 integer,parameter :: alpha(9)=(/1,2,3,3,3,2,2,1,1/),beta(9)=(/1,2,3,2,1,1,3,3,2/)
 integer,parameter :: mu9(9)=(/1,2,3,4,5,6,4,5,6/)
 integer,allocatable :: atindx(:)
 real(dp) :: gmet(3,3),gprimd(3,3),hatstr(6),rdum(1),rmet(3,3),tmp(6)
 real(dp) :: work1(dyfr_cplex,3,3),work2(dyfr_cplex,3,3)
 real(dp),allocatable :: buf(:,:),dum_atmrho1(:,:),dum_gauss(:),dyfr(:,:,:,:,:),grhat_tmp(:,:)
 real(dp),allocatable :: prod(:,:),prodp(:,:),vloc(:),vpsp1(:,:)
 type(coeff2_type),allocatable :: prod_nondiag(:),prodp_nondiag(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 if (optgr2==1.and.pawrhoij(1)%ngrhoij==0) then
   msg='Inconsistency between variables optgr2 and ngrhoij !'
   MSG_BUG(msg)
 end if
 qne0=0;if (qphon(1)**2+qphon(2)**2+qphon(3)**2>=1.d-15) qne0=1
 if (optgr2==1.and.pawfgrtab(1)%rfgd_allocated==0.and.qne0==1) then
   MSG_BUG('  pawfgrtab()%rfgd array must be allocated  !')
 end if
 if (mpi_enreg%nproc_atom>1) then
   if (natom/=mpi_enreg%natom) then
     msg='natom not equal to mpi_enreg%natom !'
     MSG_BUG(msg)
   end if
 end if

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Initializations and allocations
 ngrhat=0;ngrad=0;ngradp=0;ngrad_nondiag=0;ngradp_nondiag=0
 ishift_gr=0;ishift_gr2=0;ishift_str=0
 cplex=1;if (qne0==1) cplex=2
 if (optgr==1) then
   ngrad=ngrad+3
   ngrhat=ngrhat+3
   ishift_gr2=ishift_gr2+3
 end if
 if (optgr2==1) then
   mu=min(dyfr_cplex,cplex)
   ngrad =ngrad +9
   ngradp=ngradp+3
   ngrad_nondiag =ngrad_nondiag +9*mu
   ngradp_nondiag=ngradp_nondiag+3*mu
   ngrhat=ngrhat+9*mu
 end if
 if (optstr==1) then
   hatstr=zero
   ngrad=ngrad+6
   ngrhat=ngrhat+6
   ishift_gr=ishift_gr+6
   ishift_gr2=ishift_gr2+6
 end if

 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 fact_ucvol=ucvol/dble(nfftot)
 nsploop=nspden;if (dimvtrial<nspden) nsploop=2
 if (optgr2/=1) then
   ABI_ALLOCATE(grhat_tmp,(ngrhat,1))
 else
   ABI_ALLOCATE(dyfr,(dyfr_cplex,3,3,natom_tot,natom_tot))
   dyfr=zero
   ABI_ALLOCATE(grhat_tmp,(ngrhat,natom_tot))
   ABI_ALLOCATE(prod_nondiag,(natom_tot))
   ABI_ALLOCATE(prodp_nondiag,(natom_tot))
   ABI_ALLOCATE(atindx,(natom_tot))
   ABI_ALLOCATE(vpsp1,(cplex*nfft,3))
   atindx(:)=0
   do iatom=1,natom_tot
     iatm=0
     do while (atindx(iatom)==0.and.iatm<natom_tot)
       iatm=iatm+1;if (atindx1(iatm)==iatom) atindx(iatom)=iatm
     end do
   end do
 end if

!The computation of dynamical matrix requires the knowledge of
!g_l(r).Y_lm(r) and derivatives for all atoms
 if (optgr2==1) then
   do jatom=1,natom_tot ! NOTE: Not compatible with parallelization over atoms
     opt1=0;opt2=0;opt3=0
     lm_sizej=pawfgrtab(jatom)%l_size**2
     if (pawfgrtab(jatom)%gylm_allocated==0) then
       if (associated(pawfgrtab(jatom)%gylm))  then
         ABI_DEALLOCATE(pawfgrtab(jatom)%gylm)
       end if
       ABI_ALLOCATE(pawfgrtab(jatom)%gylm,(pawfgrtab(jatom)%nfgd,lm_sizej))
       pawfgrtab(jatom)%gylm_allocated=2;opt1=1
     end if
     if (pawfgrtab(jatom)%gylmgr_allocated==0) then
       if (associated(pawfgrtab(jatom)%gylmgr))  then
         ABI_DEALLOCATE(pawfgrtab(jatom)%gylmgr)
       end if
       ABI_ALLOCATE(pawfgrtab(jatom)%gylmgr,(3,pawfgrtab(jatom)%nfgd,lm_sizej))
       pawfgrtab(jatom)%gylmgr_allocated=2;opt2=1
     end if
     call pawgylm(pawfgrtab(jatom)%gylm,pawfgrtab(jatom)%gylmgr,&
&     pawfgrtab(jatom)%gylmgr2,lm_sizej,pawfgrtab(jatom)%nfgd,&
&     opt1,opt2,opt3,pawtab(typat(jatom)),pawfgrtab(jatom)%rfgd,&
&     pawfgrtab(jatom)%rfgd_allocated)
   end do
 end if

!Loops over types and atoms
 iatshft=0
 do itypat=1,ntypat

   lmn2_size=pawtab(itypat)%lmn2_size
   do iatm=iatshft+1,iatshft+nattyp(itypat)
     iatom=atindx1(iatm)
     iatom_tot=iatom;if (mpi_enreg%nproc_atom>1) iatom_tot=mpi_enreg%atom_indx(iatom)
     idiag=1;if (optgr2==1) idiag=iatm
     lm_size=pawfgrtab(iatom)%l_size**2
     nfgd=pawfgrtab(iatom)%nfgd

     ABI_ALLOCATE(vloc,(nfgd))
     if (ngrad>0)  then
       ABI_ALLOCATE(prod,(ngrad,lm_size))
     end if
     if (ngradp>0)  then
       ABI_ALLOCATE(prodp,(ngradp,lm_size))
     end if
     if (optgr2==1) then
       do jatm=1,natom_tot ! NOTE: Not compatible with parallelization over atoms
         lm_sizej=pawfgrtab(atindx1(jatm))%l_size**2
         ABI_ALLOCATE(prod_nondiag (jatm)%value,(ngrad_nondiag ,lm_sizej))
         ABI_ALLOCATE(prodp_nondiag(jatm)%value,(ngradp_nondiag,lm_sizej))
       end do
     end if

     grhat_tmp=zero

!    Eventually compute g_l(r).Y_lm(r) derivatives for the current atom (if not already done)
     if ((optgr==1.or.optstr==1).and.(optgr2/=1)) then
       if (pawfgrtab(iatom)%gylmgr_allocated==0) then
         if (associated(pawfgrtab(iatom)%gylmgr))  then
           ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
         end if
         ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(3,pawfgrtab(iatom)%nfgd,lm_size))
         pawfgrtab(iatom)%gylmgr_allocated=2
         call pawgylm(rdum,pawfgrtab(iatom)%gylmgr,rdum,&
&         lm_size,pawfgrtab(iatom)%nfgd,0,1,0,pawtab(itypat),&
&         pawfgrtab(iatom)%rfgd,pawfgrtab(iatom)%rfgd_allocated)
       end if
     end if
     if (optgr2==1) then
       opt1=0;opt2=0;opt3=0
       if (pawfgrtab(iatom)%gylmgr_allocated==0) then
         if (associated(pawfgrtab(iatom)%gylmgr))  then
           ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
         end if
         ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(3,pawfgrtab(iatom)%nfgd,lm_size))
         pawfgrtab(iatom)%gylmgr_allocated=2;opt2=1
       end if
       if (pawfgrtab(iatom)%gylmgr2_allocated==0) then
         if (associated(pawfgrtab(iatom)%gylmgr2))  then
           ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr2)
         end if
         ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr2,(6,pawfgrtab(iatom)%nfgd,lm_size))
         pawfgrtab(iatom)%gylmgr2_allocated=2;opt3=1
       end if
       call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,&
&       pawfgrtab(iatom)%gylmgr2,lm_size,pawfgrtab(iatom)%nfgd,&
&       opt1,opt2,opt3,pawtab(itypat),pawfgrtab(iatom)%rfgd,&
&       pawfgrtab(iatom)%rfgd_allocated)
     end if

!    Eventually compute exp(-i.q.r) factors for the current atom (if not already done)
     if (optgr2==1.and.qne0==1.and.(pawfgrtab(iatom)%expiqr_allocated==0)) then
       if (associated(pawfgrtab(iatom)%expiqr))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%expiqr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%expiqr,(2,nfgd))
       pawfgrtab(iatom)%expiqr_allocated=2
       call pawexpiqr(gprimd,pawfgrtab(iatom),qphon,xred(:,iatom_tot))
     end if
     has_phase=(optgr2==1.and.pawfgrtab(iatom)%expiqr_allocated/=0)

!    Eventually compute 1st-order potential
     if (optgr2==1) then
       optv=1;optn=0;optn2=1;idir=0
       call atm2fft3(atindx,dum_atmrho1,vpsp1,cplex,dum_gauss,gmet,gsqcut,idir,iatom,&
&       mgfft,mpi_enreg,psps%mqgrid_vl,natom_tot,3,nfft,ngfft,ntypat,optn,optn2,optv,&
&       paral_kgb,pawtab,ph1d,psps%qgrid_vl,qphon,typat,ucvol,psps%usepaw,psps%vlspl,xred)
       if (cplex==1) then
         do ic=1,nfft
           tmp(1:3)=vpsp1(ic,1:3)
           do mu=1,3
             vpsp1(ic,mu)=-(gprimd(mu,1)*tmp(1)+gprimd(mu,2)*tmp(2)+gprimd(mu,3)*tmp(3))
           end do
         end do
       else ! cplex=2
         do ic=1,nfft
           jc=2*ic;tmp(1:3)=vpsp1(jc-1,1:3);tmp(4:6)=vpsp1(jc,1:3)
           do mu=1,3
             vpsp1(jc-1,mu)=-(gprimd(mu,1)*tmp(1)+gprimd(mu,2)*tmp(2)+gprimd(mu,3)*tmp(3))
             vpsp1(jc  ,mu)=-(gprimd(mu,1)*tmp(4)+gprimd(mu,2)*tmp(5)+gprimd(mu,3)*tmp(6))
           end do
         end do
       end if

     end if

!    Loop over spin components
     do ispden=1,nsploop

!      ----- Retrieve potential (subtle if nspden=4 ;-)
       if (nspden/=4) then
         ispvtr=min(dimvtrial,ispden)
         do ic=1,nfgd
           vloc(ic)=vtrial(pawfgrtab(iatom)%ifftsph(ic),ispvtr)
         end do
       else
         if (ispden==1) then
           ispvtr=min(dimvtrial,2)
           do ic=1,nfgd
             jc=pawfgrtab(iatom)%ifftsph(ic)
             vloc(ic)=half*(vtrial(jc,1)+vtrial(jc,ispvtr))
           end do
         else if (ispden==4) then
           ispvtr=min(dimvtrial,2)
           do ic=1,nfgd
             jc=pawfgrtab(iatom)%ifftsph(ic)
             vloc(ic)=half*(vtrial(jc,1)-vtrial(jc,ispvtr))
           end do
         else if (ispden==2) then
           ispvtr=min(dimvtrial,3)
           do ic=1,nfgd
             jc=pawfgrtab(iatom)%ifftsph(ic)
             vloc(ic)=vtrial(jc,ispvtr)
           end do
         else ! ispden=3
           ispvtr=min(dimvtrial,4)
           do ic=1,nfgd
             jc=pawfgrtab(iatom)%ifftsph(ic)
             vloc(ic)=-vtrial(jc,ispvtr)
           end do
         end if
       end if

!      ----- Compute projected scalars (integrals of vloc and Q_ij^hat)
!      ----- and/or their derivatives

       if (ngrad>0) prod=zero
       if (ngradp>0) prodp=zero

!      ==== Contribution to forces ====
       if (optgr==1) then
         do ilm=1,lm_size
           do ic=1,pawfgrtab(iatom)%nfgd
             do mu=1,3
               prod(mu+ishift_gr,ilm)=prod(mu+ishift_gr,ilm)-&
&               vloc(ic)*pawfgrtab(iatom)%gylmgr(mu,ic,ilm)
             end do
           end do
         end do
       end if
!      ==== Contribution to stresses ====
       if (optstr==1) then
         do ilm=1,lm_size
           do ic=1,pawfgrtab(iatom)%nfgd
             do mu=1,6
               mua=alpha(mu);mub=beta(mu)
               prod(mu+ishift_str,ilm)=prod(mu+ishift_str,ilm) &
&               +half*vloc(ic) &
&               *(pawfgrtab(iatom)%gylmgr(mua,ic,ilm)*pawfgrtab(iatom)%rfgd(mub,ic)&
&               +pawfgrtab(iatom)%gylmgr(mub,ic,ilm)*pawfgrtab(iatom)%rfgd(mua,ic))
             end do
           end do
         end do
       end if
!      ==== Contribution to frozen wf part of dyn. matrix ====
       if (optgr2==1) then
!        Diagonal contribution
         do ilm=1,lm_size
           do ic=1,pawfgrtab(iatom)%nfgd
             do mu=1,9
               prod(ishift_gr2+mu,ilm)=prod(ishift_gr2+mu,ilm) &
&               +half*vloc(ic)*pawfgrtab(iatom)%gylmgr2(mu9(mu),ic,ilm)
             end do
             do mu=1,3
               prodp(mu,ilm)=prodp(mu,ilm) &
&               -vloc(ic)*pawfgrtab(iatom)%gylmgr(mu,ic,ilm)
             end do
           end do
         end do
!        Off-diagonal contribution
         do jatm=1,natom_tot ! NOTE: Not compatible with parallelization over atoms
           jatom=atindx1(jatm)
           lm_sizej=pawfgrtab(jatom)%l_size**2
           prod_nondiag (jatm)%value=zero
           prodp_nondiag(jatm)%value=zero
           if (has_phase.or.cplex==2) then
             if (dyfr_cplex==1.or.cplex==1) then
               do ilm=1,lm_sizej
                 do ic=1,pawfgrtab(jatom)%nfgd
                   jc=2*pawfgrtab(jatom)%ifftsph(ic)
                   tmp(1:3)=vpsp1(jc-1,1:3)*pawfgrtab(jatom)%expiqr(1,ic) &
&                   -vpsp1(jc  ,1:3)*pawfgrtab(jatom)%expiqr(2,ic)
                   do mu=1,9
                     mua=alpha(mu);mub=beta(mu)
                     prod_nondiag(jatm)%value(mu,ilm)=prod_nondiag(jatm)%value(mu,ilm) &
&                     +tmp(mua)*pawfgrtab(jatom)%gylmgr(mub,ic,ilm)
                   end do
                   do mu=1,3
                     prodp_nondiag(jatm)%value(mu,ilm)=prodp_nondiag(jatm)%value(mu,ilm) &
&                     -tmp(mu)*pawfgrtab(jatom)%gylm(ic,ilm)
                   end do
                 end do
               end do
             else
               do ilm=1,lm_sizej
                 do ic=1,pawfgrtab(jatom)%nfgd
                   jc=2*pawfgrtab(jatom)%ifftsph(ic)
                   tmp(1:3)=vpsp1(jc-1,1:3)*pawfgrtab(jatom)%expiqr(1,ic) &
&                   -vpsp1(jc  ,1:3)*pawfgrtab(jatom)%expiqr(2,ic)
                   tmp(4:6)=vpsp1(jc-1,1:3)*pawfgrtab(jatom)%expiqr(2,ic) &
&                   +vpsp1(jc  ,1:3)*pawfgrtab(jatom)%expiqr(1,ic)
                   do mu=1,9
                     mua=alpha(mu);mub=beta(mu)
                     prod_nondiag(jatm)%value(mu,ilm)=prod_nondiag(jatm)%value(mu,ilm) &
&                     +tmp(mua)*pawfgrtab(jatom)%gylmgr(mub,ic,ilm)
                     prod_nondiag(jatm)%value(9+mu,ilm)=prod_nondiag(jatm)%value(9+mu,ilm) &
&                     +tmp(3+mua)*pawfgrtab(jatom)%gylmgr(mub,ic,ilm)
                   end do
                   do mu=1,3
                     prodp_nondiag(jatm)%value(mu,ilm)=prodp_nondiag(jatm)%value(mu,ilm) &
&                     -tmp(mu)*pawfgrtab(jatom)%gylm(ic,ilm)
                     prodp_nondiag(jatm)%value(3+mu,ilm)=prodp_nondiag(jatm)%value(3+mu,ilm) &
&                     -tmp(3+mu)*pawfgrtab(jatom)%gylm(ic,ilm)
                   end do
                 end do
               end do
             end if
           else ! no phase
             do ilm=1,lm_sizej
               do ic=1,pawfgrtab(jatom)%nfgd
                 jc=pawfgrtab(jatom)%ifftsph(ic)
                 do mu=1,9
                   mua=alpha(mu);mub=beta(mu)
                   prod_nondiag(jatm)%value(mu,ilm)=prod_nondiag(jatm)%value(mu,ilm) &
&                   +vpsp1(jc,mua)*pawfgrtab(jatom)%gylmgr(mub,ic,ilm)
                 end do
                 do mu=1,3
                   prodp_nondiag(jatm)%value(mu,ilm)=prodp_nondiag(jatm)%value(mu,ilm) &
&                   -vpsp1(jc,mu)*pawfgrtab(jatom)%gylm(ic,ilm)
                 end do
               end do
             end do
           end if
         end do
       end if

!      --- Reduction in case of parallelization ---
       if(mpi_enreg%paral_compil_fft==1)then
         old_paral_level= mpi_enreg%paral_level
         mpi_enreg%paral_level=3
         call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
         if (ngrad>0) then
           call xsum_mpi(prod,spaceComm,ier)
         end if
         if (ngradp>0) then
           call xsum_mpi(prodp,spaceComm,ier)
         end if
         if (optgr2==1) then
           bufsiz=0;bufind=0
           do jatm=1,natom_tot
             bufsiz=bufsiz+pawfgrtab(atindx1(jatm))%l_size**2
           end do
           ABI_ALLOCATE(buf,(ngrad_nondiag+ngradp_nondiag,bufsiz))
           do jatm=1,natom_tot ! NOTE: Not compatible with parallelization over atoms
             lm_sizej=pawfgrtab(atindx1(jatm))%l_size**2
             if (ngrad_nondiag> 0) buf(1:ngrad_nondiag,bufind+1:bufind+lm_sizej)= &
&             prod_nondiag (jatm)%value(:,:)
             if (ngradp_nondiag>0) buf(ngrad_nondiag+1:ngrad_nondiag+ngradp_nondiag, &
&             bufind+1:bufind+lm_sizej)=prodp_nondiag(jatm)%value(:,:)
             bufind=bufind+lm_sizej*(ngrad_nondiag+ngradp_nondiag)
           end do
           call xsum_mpi(buf,spaceComm,ier)
           bufind=0
           do jatm=1,natom_tot ! NOTE: Not compatible with parallelization over atoms
             lm_sizej=pawfgrtab(atindx1(jatm))%l_size**2
             if (ngrad> 0) prod_nondiag (jatm)%value(:,:)= &
&             buf(1:ngrad_nondiag,bufind+1:bufind+lm_sizej)
             if (ngradp>0) prodp_nondiag(jatm)%value(:,:)= &
&             buf(ngrad_nondiag+1:ngrad_nondiag+ngradp_nondiag,bufind+1:bufind+lm_sizej)
             bufind=bufind+lm_sizej*(ngrad_nondiag+ngradp_nondiag)
           end do
           ABI_DEALLOCATE(buf)
         end if
         mpi_enreg%paral_level=old_paral_level
       end if

!      ---- Compute all gradients
       jrhoij=1
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(irhoij)
         klm =pawtab(itypat)%indklmn(1,klmn)
         lmin=pawtab(itypat)%indklmn(3,klmn)
         lmax=pawtab(itypat)%indklmn(4,klmn)
         ro  =pawrhoij(iatom)%rhoijp(jrhoij,ispden)
         ro_d=ro*pawtab(itypat)%dltij(klmn)
         do ll=lmin,lmax,2
           do ilm=ll**2+1,(ll+1)**2
             isel=pawang%gntselect(ilm,klm)
             if (isel>0) then
               grhat_x=ro_d*pawtab(itypat)%qijl(ilm,klmn)
               do mu=1,ngrad
                 grhat_tmp(mu,idiag)=grhat_tmp(mu,idiag)+grhat_x*prod(mu,ilm)
               end do
             end if
           end do
         end do
         jrhoij=jrhoij+pawrhoij(iatom)%cplex
       end do ! irhoij

!      ---- Add additional terms for second gradients
       if (optgr2==1) then
!        Diagonal term including rhoij derivative
         klmn1=1
         do klmn=1,lmn2_size
           klm =pawtab(itypat)%indklmn(1,klmn)
           lmin=pawtab(itypat)%indklmn(3,klmn)
           lmax=pawtab(itypat)%indklmn(4,klmn)
           dlt_tmp=pawtab(itypat)%dltij(klmn)
           do ll=lmin,lmax,2
             do ilm=ll**2+1,(ll+1)**2
               isel=pawang%gntselect(ilm,klm)
               if (isel>0) then
                 ro_d=dlt_tmp*pawtab(itypat)%qijl(ilm,klmn)
                 do mu=1,9
                   mua=alpha(mu);mub=beta(mu)
                   grhat_tmp(ishift_gr2+mu,idiag)=grhat_tmp(ishift_gr2+mu,idiag)&
&                   +ro_d*pawrhoij(iatom)%grhoij(mua,klmn1,ispden)*prodp(mub,ilm)
                 end do
               end if
             end do
           end do
           klmn1=klmn1+pawrhoij(iatom)%cplex
         end do ! klmn
         do jatm=1,natom_tot ! NOTE: Not compatible with parallelization over atoms
           jatom=atindx1(jatm);jtypat=typat(jatom)
!          Off-diagonal term including rhoij
           if (dyfr_cplex==1.or.cplex==1) then
             jrhoij=1
             do irhoij=1,pawrhoij(jatom)%nrhoijsel
               klmn=pawrhoij(jatom)%rhoijselect(irhoij)
               klm =pawtab(jtypat)%indklmn(1,klmn)
               lmin=pawtab(jtypat)%indklmn(3,klmn)
               lmax=pawtab(jtypat)%indklmn(4,klmn)
               ro  =pawrhoij(jatom)%rhoijp(jrhoij,ispden)
               ro_d=ro*pawtab(jtypat)%dltij(klmn)
               do ll=lmin,lmax,2
                 do ilm=ll**2+1,(ll+1)**2
                   isel=pawang%gntselect(ilm,klm)
                   if (isel>0) then
                     grhat_x=ro_d*pawtab(jtypat)%qijl(ilm,klmn)
                     do mu=1,9
                       grhat_tmp(ishift_gr2+mu,jatm)=grhat_tmp(ishift_gr2+mu,jatm) &
&                       +grhat_x*prod_nondiag(jatm)%value(mu,ilm)
                     end do
                   end if
                 end do
               end do
               jrhoij=jrhoij+pawrhoij(jatom)%cplex
             end do ! irhoij
           else
             jrhoij=1;mushift=ishift_gr2+9
             do irhoij=1,pawrhoij(jatom)%nrhoijsel
               klmn=pawrhoij(jatom)%rhoijselect(irhoij)
               klm =pawtab(jtypat)%indklmn(1,klmn)
               lmin=pawtab(jtypat)%indklmn(3,klmn)
               lmax=pawtab(jtypat)%indklmn(4,klmn)
               ro  =pawrhoij(jatom)%rhoijp(jrhoij,ispden)
               ro_d=ro*pawtab(jtypat)%dltij(klmn)
               do ll=lmin,lmax,2
                 do ilm=ll**2+1,(ll+1)**2
                   isel=pawang%gntselect(ilm,klm)
                   if (isel>0) then
                     grhat_x=ro_d*pawtab(jtypat)%qijl(ilm,klmn)
                     do mu=1,9
                       grhat_tmp(ishift_gr2+mu,jatm)=grhat_tmp(ishift_gr2+mu,jatm) &
&                       +grhat_x*prod_nondiag(jatm)%value(mu,ilm)
                       grhat_tmp(mushift+mu,jatm)=grhat_tmp(mushift+mu,jatm) &
&                       +grhat_x*prod_nondiag(jatm)%value(9+mu,ilm)
                     end do
                   end if
                 end do
               end do
               jrhoij=jrhoij+pawrhoij(jatom)%cplex
             end do ! irhoij
           end if
!          Off-diagonal term including rhoij derivative
           if (dyfr_cplex==1.or.cplex==1) then
             klmn1=1
             do klmn=1,pawrhoij(jatom)%lmn2_size
               klm =pawtab(jtypat)%indklmn(1,klmn)
               lmin=pawtab(jtypat)%indklmn(3,klmn)
               lmax=pawtab(jtypat)%indklmn(4,klmn)
               dlt_tmp=pawtab(jtypat)%dltij(klmn)
               do ll=lmin,lmax,2
                 do ilm=ll**2+1,(ll+1)**2
                   isel=pawang%gntselect(ilm,klm)
                   if (isel>0) then
                     ro_d=dlt_tmp*pawtab(jtypat)%qijl(ilm,klmn)
                     do mu=1,9
                       mua=alpha(mu);mub=beta(mu)
                       grhat_tmp(ishift_gr2+mu,jatm)=grhat_tmp(ishift_gr2+mu,jatm)&
&                       +ro_d*pawrhoij(jatom)%grhoij(mub,klmn1,ispden)*prodp_nondiag(jatm)%value(mua,ilm)
                     end do
                   end if
                 end do
               end do
               klmn1=klmn1+pawrhoij(jatom)%cplex
             end do ! klmn
           else ! ngradp_nondiag>=6
             klmn1=1;mushift=ishift_gr2+9
             do klmn=1,pawrhoij(jatom)%lmn2_size
               klm =pawtab(jtypat)%indklmn(1,klmn)
               lmin=pawtab(jtypat)%indklmn(3,klmn)
               lmax=pawtab(jtypat)%indklmn(4,klmn)
               dlt_tmp=pawtab(jtypat)%dltij(klmn)
               do ll=lmin,lmax,2
                 do ilm=ll**2+1,(ll+1)**2
                   isel=pawang%gntselect(ilm,klm)
                   if (isel>0) then
                     ro_d=dlt_tmp*pawtab(jtypat)%qijl(ilm,klmn)
                     do mu=1,9
                       mua=alpha(mu);mub=beta(mu)
                       grhat_tmp(ishift_gr2+mu,jatm)=grhat_tmp(ishift_gr2+mu,jatm)&
&                       +ro_d*pawrhoij(jatom)%grhoij(mub,klmn1,ispden)*prodp_nondiag(jatm)%value(mua,ilm)
                       grhat_tmp(mushift+mu,jatm)=grhat_tmp(mushift+mu,jatm)&
&                       +ro_d*pawrhoij(jatom)%grhoij(mub,klmn1,ispden)*prodp_nondiag(jatm)%value(3+mua,ilm)
                     end do
                   end if
                 end do
               end do
               klmn1=klmn1+pawrhoij(jatom)%cplex
             end do ! klmn
           end if ! prodp_nondiag
         end do ! jatm
       end if ! optgr2==1

     end do ! ispden

!    Eventually free temporary space for g_l(r).Y_lm(r) factors
     if (pawfgrtab(iatom)%gylmgr_allocated==2) then
       ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
       ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(0,0,0))
       pawfgrtab(iatom)%gylmgr_allocated=0
     end if
     if (optgr2==1) then
       if (pawfgrtab(iatom)%gylmgr2_allocated==2) then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr2)
         ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr2,(0,0,0))
         pawfgrtab(iatom)%gylmgr2_allocated=0
       end if
     end if

!    ==== Forces ====
!    Convert from cartesian to reduced coordinates
     if (optgr==1) then
       mushift=3*(iatm-1)
       tmp(1:3)=grhat_tmp(ishift_gr+1:ishift_gr+3,idiag)
       do mu=1,3
         grnl(mu+mushift)=grnl(mu+mushift)&
&         +fact_ucvol*(rprimd(1,mu)*tmp(1)+rprimd(2,mu)*tmp(2)+rprimd(3,mu)*tmp(3))
       end do
     end if
!    ==== Stresses ====
     if (optstr==1) then
       hatstr(1:6)=hatstr(1:6)+ grhat_tmp(ishift_str+1:ishift_str+6,idiag)
     end if
!    ==== Frozen wf part of dyn. matrix ====
     if (optgr2==1) then
       do jatm=1,natom_tot ! NOTE: Not compatible with parallelization over atoms
         do mu=1,9
           mua=alpha(mu);mub=beta(mu)
           dyfr(1,mub,mua,jatm,iatm)=grhat_tmp(ishift_gr2+mu,jatm)
         end do
         if (dyfr_cplex==2.and.cplex==2) then
           mushift=ishift_gr2+9
           do mu=1,9
             mua=alpha(mu);mub=beta(mu)
             dyfr(2,mub,mua,jatm,iatm)=grhat_tmp(mushift+mu,jatm)
           end do
         end if
       end do
     end if

!    End loops on types and atoms
     ABI_DEALLOCATE(vloc)
     if (ngrad>0)  then
       ABI_DEALLOCATE(prod)
     end if
     if (ngradp>0)  then
       ABI_DEALLOCATE(prodp)
     end if
     if (optgr2==1) then
       do jatm=1,natom_tot
         ABI_DEALLOCATE(prod_nondiag(jatm)%value)
         ABI_DEALLOCATE(prodp_nondiag(jatm)%value)
       end do
     end if
   end do
   iatshft=iatshft+nattyp(itypat)
 end do

!Deallocate additional memory
 ABI_DEALLOCATE(grhat_tmp)
 if (optgr2==1) then
   ABI_DEALLOCATE(atindx)
   ABI_DEALLOCATE(vpsp1)
   ABI_DEALLOCATE(prod_nondiag)
   ABI_DEALLOCATE(prodp_nondiag)
   do jatom=1,natom_tot ! NOTE: Not compatible with parallelization over atoms
     if (pawfgrtab(jatom)%gylm_allocated==2) then
       ABI_DEALLOCATE(pawfgrtab(jatom)%gylm)
       ABI_ALLOCATE(pawfgrtab(jatom)%gylm,(0,0))
       pawfgrtab(jatom)%gylm_allocated=0
     end if
     if (pawfgrtab(jatom)%gylmgr_allocated==2) then
       ABI_DEALLOCATE(pawfgrtab(jatom)%gylmgr)
       ABI_ALLOCATE(pawfgrtab(jatom)%gylmgr,(0,0,0))
       pawfgrtab(jatom)%gylmgr_allocated=0
     end if
     if (pawfgrtab(jatom)%expiqr_allocated==2) then
       ABI_DEALLOCATE(pawfgrtab(jatom)%expiqr)
       ABI_ALLOCATE(pawfgrtab(jatom)%expiqr,(0,0))
       pawfgrtab(jatom)%expiqr_allocated=0
     end if
   end do
 end if

!===== Convert stresses (add diag and off-diag contributions) =====
 if (optstr==1) then
!  Has to compute int[nhat*vtrial]
   hatstr_diag=zero
   if (nspden==1.or.dimvtrial==1) then
     do ic=1,nfft
       hatstr_diag=hatstr_diag+vtrial(ic,1)*nhat(ic,1)
     end do
   else if (nspden==2) then
     do ic=1,nfft
       hatstr_diag=hatstr_diag+vtrial(ic,1)*nhat(ic,2)+vtrial(ic,2)*(nhat(ic,1)-nhat(ic,2))
     end do
   else if (nspden==4) then
     do ic=1,nfft
       hatstr_diag=hatstr_diag+half*(vtrial(ic,1)*(nhat(ic,1)+nhat(ic,4)) &
&       +vtrial(ic,2)*(nhat(ic,1)-nhat(ic,4))) &
&       +vtrial(ic,3)*nhat(ic,2)+vtrial(ic,4)*nhat(ic,3)
     end do
   end if
   if(mpi_enreg%paral_compil_fft==1)then
     old_paral_level= mpi_enreg%paral_level
     mpi_enreg%paral_level=3
     call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
     call xsum_mpi(hatstr_diag,spaceComm,ier)
     mpi_enreg%paral_level=old_paral_level
   end if
!  Convert hat contribution
   hatstr(1:3)=(hatstr(1:3)+hatstr_diag)/dble(nfftot)
   hatstr(4:6)= hatstr(4:6)/dble(nfftot)
!  Add to already computed NL contrib
   nlstr(1:6)=nlstr(1:6)+hatstr(1:6)
!  Apply symmetries
   call stresssym(gprimd,nsym,nlstr,symrec)
 end if

!===== Convert dynamical matrix (from cartesian to reduced coordinates) =====
 if (optgr2==1) then
   do iatm=1,natom_tot
     do jatm=1,natom_tot
       do mua=1,3
         do mub=1,3
           work1(1,mua,mub)=dyfr(1,mub,mua,jatm,iatm)+dyfr(1,mua,mub,iatm,jatm)
         end do
       end do
       if (dyfr_cplex==2) then
         do mua=1,3
           do mub=1,3
             work1(2,mua,mub)=dyfr(2,mub,mua,jatm,iatm)-dyfr(2,mua,mub,iatm,jatm)
           end do
         end do
       end if
       do mu=1,3
         work2(:,:,mu)=rprimd(1,mu)*work1(:,:,1)+rprimd(2,mu)*work1(:,:,2)+rprimd(3,mu)*work1(:,:,3)
       end do
       do mub=1,3
         do mua=1,3
           dyfrnl(:,mua,mub,jatm,iatm)=dyfrnl(:,mua,mub,jatm,iatm) &   ! Already contains NL projectors contribution
&          +fact_ucvol*(rprimd(1,mua)*work2(:,1,mub) &
&           +rprimd(2,mua)*work2(:,2,mub) &
&           +rprimd(3,mua)*work2(:,3,mub))
         end do
       end do
     end do
   end do
   ABI_DEALLOCATE(dyfr)
 end if

 DBG_ENTER("COLL")

end subroutine pawgrnl
!!***
