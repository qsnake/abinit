!{\src2tex{textfont=tt}}
!!****f* ABINIT/dyfnl3
!! NAME
!! dyfnl3
!!
!! FUNCTION
!! Compute the frozen-wavefunction non-local contribution to the
!! dynamical matrix.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GM, AR, MB, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>=Fourier coefficients of WF
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)=<p_lmn|C> coefficients for WF |C> (and 1st derivatives)
!!  dimcprj(natom*usepaw)=array of dimensions of array cprj (ordered by atom-type)
!!  dyfr_cplex=1 if dyfrnl is real, 2 if it is complex
!!  dyfr_nondiag=1 if dyfrnl is non diagonal with respect to atoms; 0 otherwise
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gsqcut_eff=Fourier cutoff on G^2 for "large sphere" of radius double that of the basis sphere
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!  kptns(3,nkpt)=coordinates of k points in terms of reciprocal space
!!   primitive translations
!!  kptopt=option for the generation of k points
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mgfftf=maximum size of 1D FFTs for the fine FFT grid (PAW)
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimension for number of planewaves
!!  natom=number of atoms in unit cell
!!  nattyp(ntypat)=array describing how many atoms of each type in cell
!!  nband(nkpt*nsppol)=number of bands being considered per k point
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngfftf(18)= -PAW ONLY- contain all needed information about 3D FFT for the fine grid
!!              (ngfftf=ngfft for norm-conserving potential runs)
!!  nkpt=number of k points
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for polarized
!!  nsym=number of symmetry elements in space group (at least 1)
!!  ntypat=integer specification of atom type (1, 2, ...)
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2) at each
!!    k point
!!  paral_kgb=flag for (kpt,FFT,bands) parallelism
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawbec= flag for the computation of Born Effective Charge within PAW ; set to 1 if yes
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase information related to structure factor
!!  ph1df(2,3*(2*mgfftf+1)*natom)=phase information related to structure factor on the fine FFT grid (PAW)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qphon(3)=wavevector of the phonon
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries in reciprocal space (dimensionless)
!!  typat(natom)=type integer for each atom in cell
!!  unkg=unit number for (k+G) sphere data
!!  unpaw=unit number for cprj PAW data (if used)
!!  unylm=unit number for disk file containing Ylm''s if mkmem==0
!!  usecprj=1 if cprj coefficients are already in memory (PAW only)
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vtrial(nfftf,nspden)=total potential (Hartree+XC+loc)
!!  vxc(nfftf,nspden)=XC potential
!!  wfftgs=struct info for disk file containing GS wavefunctions if mkmem==0
!!  wtk(nkpt)=k point weights
!!  xred(3,natom)=reduced coordinates of atoms (dimensionless)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!
!! OUTPUT
!!  becfrnl(3,natom,3*pawbec)=NL frozen contribution to Born Effective Charges (PAW only)
!!                            (3,natom) = derivative wr to the displ. of one atom in one direction
!!                            (3)       = derivative wr to electric field in one direction
!!  dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)=
!!         non-symmetrized non-local contribution to the dynamical matrix
!!         If NCPP, it depends on one atom
!!         If PAW,  it depends on two atoms
!!
!! SIDE EFFECTS
!!  ===== if psps%usepaw==1
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!                          pawfgrtab(:)%gylmgr2 are deallocated here
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!    (gradients of rhoij for each atom with respect to atomic positions are computed here)
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      cprj_alloc,cprj_copy,cprj_diskinit_r,cprj_free,cprj_get,destroy_paw_ij
!!      hdr_skip,init_paw_ij,leave_test,metric,mkffnl,mkkpg,nonlop
!!      nullify_paw_ij,pawaccrhoij,pawfrnhat,pawgrnl,ph1d3d,rdnpw,rwwf,symrhoij
!!      timab,xcomm_world,xdefineoff,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dyfnl3(atindx1,becfrnl,cg,cprj,dimcprj,dyfrnl,dyfr_cplex,dyfr_nondiag,eigen,gsqcut,indsym,&
&          istwfk,kg,kptns,kptopt,mband,mgfft,mgfftf,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nband,&
&          nfftf,ngfft,ngfftf,nkpt,nloalg,npwarr,nspden,nspinor,nsppol,nsym,ntypat,occ,&
&          paral_kgb,paw_ij,pawang,pawbec,pawprtvol,pawfgrtab,pawrad,pawrhoij,pawtab,ph1d,ph1df,psps,&
&          qphon,rprimd,symafm,symrec,typat,unkg,unpaw,unylm,usecprj,usexcnhat,wfftgs,vtrial,&
&          vxc,wtk,xred,ylm)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_paw_toolbox
 use m_xmpi
 use m_errors
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dyfnl3'
 use interfaces_18_timing
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_65_nonlocal
 use interfaces_66_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dyfr_cplex,dyfr_nondiag,kptopt,mband,mgfft,mgfftf,mkmem,mpsang,mpw,natom,nfftf
 integer,intent(in) :: nkpt,nspden,nspinor,nsppol,nsym,ntypat,paral_kgb,pawbec,pawprtvol,unkg
 integer,intent(in) :: unpaw,unylm,usecprj,usexcnhat
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wfftgs
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom*psps%usepaw)
 integer,intent(in) :: indsym(4,nsym,natom),istwfk(nkpt),kg(3,mpw*mkmem)
 integer,intent(in) :: nattyp(ntypat),nband(nkpt*nsppol),ngfft(18),ngfftf(18)
 integer,intent(in) :: nloalg(5),npwarr(nkpt),symafm(nsym),symrec(3,3,nsym)
 integer,intent(in) :: typat(ntypat)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),kptns(3,nkpt),occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),ph1df(2,3*(2*mgfftf+1)*natom)
 real(dp),intent(in) :: qphon(3),rprimd(3,3),vxc(nfftf,nspden),wtk(nkpt)
 real(dp),intent(in) :: xred(3,natom),ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),target,intent(in) :: vtrial(nfftf,nspden)
 real(dp),intent(out) :: becfrnl(3,natom,3*pawbec)
 real(dp),intent(out) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
 type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)
 type(paw_ij_type),intent(in) :: paw_ij(natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(inout) :: pawrhoij(natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,bufdim,choice,choice_bec,cplex,cplx,cpopt,cpopt_bec,dimdij,dimekb1,dimekb2
 integer :: dimffnl,dimnhat,dimvtmp,formeig,ia,iatom,iband,ibg,ibsp,icg,ider,idir,ierr
 integer :: ii,ikg,ikpt,ilm,iorder_cprj,ipw,isp,ispden,isppol
 integer :: istwf_k,itypat,jj,klmn,master,matblk,mcg_disk,me,mu,n1,n2,n3
 integer :: nband_k,nkpg,nnlout,nnlout_bec,npw_k,nsp,nsploop
 integer :: optgr,optgr2,option,option_rhoij,optstr,paw_opt,paw_opt_bec
 integer :: signs,spaceworld,tim_nonlop,tim_rwwf
 logical :: usetimerev
 real(dp) :: arg,eig_k,occ_k,ucvol,wtk_k
!arrays
 integer,allocatable :: dimlmn(:),kg_dum(:,:),kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),grhoij(3),kpoint(3),nonlop_dum(1,1)
 real(dp) :: rmet(3,3),tsec(2),ylmgr_dum(1)
 real(dp),allocatable :: becfrnl_tmp(:,:,:),becij(:,:,:,:)
 real(dp),allocatable :: buffer1(:),buffer2(:),cg_disk(:,:),cwavef(:,:)
 real(dp),allocatable :: dummy(:),dyfrnlk(:,:),eig_dum(:),ekb(:,:,:),enlout(:),enlout_bec(:)
 real(dp),allocatable :: ffnl(:,:,:,:),kpg_k(:,:),nhat_dum(:,:),occ_dum(:)
 real(dp),allocatable :: ph3d(:,:,:),phkxred(:,:),sij(:,:),ylm_k(:,:)
 real(dp),pointer :: vtmp(:,:)
 type(cprj_type),allocatable :: cprj_disk(:,:),cwaveprj(:,:)
 type(pawfgrtab_type),allocatable :: pawfgrtab_tmp(:)
 type(paw_ij_type),allocatable :: paw_ij_tmp(:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(159,1,tsec)

!Default for sequential use
 master=0
 call xme_init(mpi_enreg,me)
!Init spaceworld
 call xcomm_world(mpi_enreg,spaceworld)




!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Prepare temporary files if mkmem==0
!Wavefunction file
 if (mkmem==0) then
!  Read header
   call hdr_skip(wfftgs,ierr)
!  Define offsets, in case of MPI I/O
   formeig=0
   call xdefineOff(formeig,wfftgs,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)
   mcg_disk=mpw*nspinor*mband
   ABI_ALLOCATE(cg_disk,(2,mcg_disk))
 end if
!PAW file
 iorder_cprj=0
 if (usecprj==1) then
   call cprj_diskinit_r(atindx1,natom,iorder_cprj,mkmem,natom,3,dimcprj,nspinor,unpaw)
 end if

 dyfrnl(:,:,:,:,:)=zero
 bdtot_index=0
 ibg=0;icg=0
 nsploop=nsppol;if (nspden==4) nsploop=4

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(cwavef,(2,mpw*nspinor))
 ABI_ALLOCATE(dyfrnlk,(6,natom))
 ABI_ALLOCATE(phkxred,(2,natom))

!Common data for "nonlop" routine
 signs=1 ; idir=0 ; tim_nonlop=6

!Common data for "nonlop" routine
 signs=1 ; idir=0  ; tim_nonlop=4
 choice=4 ; nnlout=max(1,6*natom)
 ABI_ALLOCATE(enlout,(nnlout))
 if (psps%usepaw==0) then
   paw_opt=0 ; cpopt=-1
 else
   paw_opt=2 ; cpopt=1+2*usecprj
 end if
 if (pawbec==1) then
   choice_bec=2 ; nnlout_bec=max(1,3*natom) ; paw_opt_bec=1 ; cpopt_bec=4
   ABI_ALLOCATE(enlout_bec,(nnlout_bec))
 else
   choice_bec=0 ; nnlout_bec=0 ; paw_opt_bec=0 ; cpopt_bec=0
 end if

!Non-local factors:
!Norm-conserving: kleimann-Bylander energies
!PAW: Dij coefficients and overlap coefficients
 if (psps%usepaw==0) then   !NCPP
   dimekb1=psps%dimekb;dimekb2=ntypat
   ABI_ALLOCATE(ekb,(psps%dimekb,ntypat,nspinor**2))
   ekb(:,:,1)=psps%ekb(:,:)
   if (nspinor==2) then
     ekb(:,:,2)=psps%ekb(:,:)
     ekb(:,:,3:4)=zero
   end if
 else                       !PAW
!  Dij coefficients
   dimekb1=psps%dimekb*paw_ij(1)%cplex_dij;dimekb2=natom
   ABI_ALLOCATE(ekb,(dimekb1,dimekb2,nspinor**2))
   ABI_ALLOCATE(sij,(dimekb1,ntypat))
   do itypat=1,ntypat
     if (paw_ij(1)%cplex_dij==1) then
       sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
     else
       do klmn=1,pawtab(itypat)%lmn2_size
         sij(2*klmn-1,itypat)=pawtab(itypat)%sij(klmn)
         sij(2*klmn  ,itypat)=zero
       end do
     end if
   end do
!  If PAW and Born Eff., Charges has to compute some additional data:
!  For each atom and for electric field direction k:
!  becij(k)=<Phi_i|r_k-R_k|Phi_j>-<tPhi_i|r_k-R_k|tPhi_j> + sij.R_k
   if (pawbec==1) then
     ABI_ALLOCATE(becij,(psps%dimekb,dimekb2,nspinor**2,3))
     becij=zero
     ABI_ALLOCATE(paw_ij_tmp,(natom))
     ABI_ALLOCATE(pawfgrtab_tmp,(natom))
     call nullify_paw_ij(paw_ij_tmp)
     call init_paw_ij(paw_ij_tmp,1,1,1,1,1,0,natom,ntypat,typat,pawtab,has_dijfr=1)
     do ii=1,3 ! Loop over direction of electric field
       call pawfrnhat(1,gprimd,ii,natom+2,mpi_enreg,natom,natom,nfftf,ngfftf,nspden,&
&       ntypat,1,paw_ij_tmp,pawang,pawfgrtab_tmp,pawrad,pawrhoij,pawtab,&
&       (/zero,zero,zero/),rprimd,ucvol,vtrial,vtrial,xred) ! vtrial not used here
       do iatom=1,natom
         itypat=typat(iatom);dimdij=pawtab(itypat)%lmn2_size
!        Add contribution from Phi_i/Phi_j and tPhi_i/tPhi_j to becij
         do klmn=1,dimdij
           becij(klmn,iatom,1,ii)=paw_ij_tmp(iatom)%dijfr(klmn,1)
         end do
!        Add contribution from sij to becij
!        xred are atomic positions in reduced cooridinates of real space
!        need them in reduced coordinates of reciprocal space
         arg=gmet(ii,1)*xred(1,iatom)+gmet(ii,2)*xred(2,iatom)+gmet(ii,3)*xred(3,iatom)
         if (paw_ij(1)%cplex_dij==1) then
           becij(1:dimdij,iatom,1,ii)=becij(1:dimdij,iatom,1,ii)+arg*sij(1:dimdij,itypat)
         else
           do klmn=1,dimdij
             becij(klmn,iatom,1,ii)=becij(dimdij,iatom,1,ii)+arg*sij(2*klmn-1,itypat)
           end do
         end if
       end do
     end do
     if (nspinor==2) becij(:,:,2,:)=becij(:,:,1,:)
     call destroy_paw_ij(paw_ij_tmp)
     ABI_DEALLOCATE(paw_ij_tmp)
     ABI_DEALLOCATE(pawfgrtab_tmp)
   end if
!  Projected WF (cprj)
   ABI_ALLOCATE(cwaveprj,(natom,nspinor))
   call cprj_alloc(cwaveprj,3,dimcprj)
   do iatom=1,natom
     pawrhoij(iatom)%ngrhoij=3
     ABI_ALLOCATE(pawrhoij(iatom)%grhoij,(3,pawrhoij(iatom)%cplex*pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
     pawrhoij(iatom)%grhoij=zero
   end do
   option_rhoij=3
   usetimerev=(kptopt>0.and.kptopt<3)
 end if

!LOOP OVER SPINS
 do isppol=1,nsppol


!  Rewind kpgsph data file if needed:
   if (mkmem==0) rewind(unkg)
   if (mkmem==0.and.psps%useylm==1) rewind unylm

!  PAW: retrieve Dij coefficients
   if (psps%usepaw==1) then
     do ispden=1,nspinor**2
       isp=isppol;if (nspinor==2) isp=ispden
       do iatom=1,natom
         dimdij=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
         do klmn=1,dimdij
           ekb(klmn,iatom,ispden)=paw_ij(iatom)%dij(klmn,isp)
         end do
         if(dimdij+1<=dimekb1) ekb(dimdij+1:dimekb1,iatom,ispden)=zero
       end do
     end do
   end if

   ikg=0

!  Loop over k points
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     istwf_k=istwfk(ikpt)
     npw_k=npwarr(ikpt)
     wtk_k=wtk(ikpt)

     if(mpi_enreg%paral_compil_kpt==1)then
!      Skip this k-point if not the proper processor
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))/=0) then
         bdtot_index=bdtot_index+nband_k
         cycle
       end if
     end if

     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     kpoint(:)=kptns(:,ikpt)

     kg_k(:,:) = 0
     if (mkmem==0) then

       nsp=nspinor
       call rdnpw(ikpt,isppol,nband_k,npw_k,nsp,0,unkg)
!      Read k+g data
       read (unkg) ((kg_k(ii,ipw),ii=1,3),ipw=1,npw_k)
!      Eventually read spherical harmonics
       if (psps%useylm==1) then
         read(unylm)
         read(unylm) ((ylm_k(ipw,ilm),ipw=1,npw_k),ilm=1,mpsang*mpsang)
       end if
!      Read the wavefunction block for ikpt,isppol
       tim_rwwf=14
       ABI_ALLOCATE(eig_dum,(mband))
       ABI_ALLOCATE(kg_dum,(3,0))
       ABI_ALLOCATE(occ_dum,(mband))
       call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,mband,mcg_disk,mpi_enreg,&
&       nband_k,nband_k,npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wfftgs)
       ABI_DEALLOCATE(eig_dum)
       ABI_DEALLOCATE(kg_dum)
       ABI_DEALLOCATE(occ_dum)

     else

!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(ikg,kg,kg_k,npw_k)
       do ipw=1,npw_k
         kg_k(1,ipw)=kg(1,ipw+ikg)
         kg_k(2,ipw)=kg(2,ipw+ikg)
         kg_k(3,ipw)=kg(3,ipw+ikg)
       end do
!      $OMP END PARALLEL DO
       if (psps%useylm==1) then
!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(ikg,npw_k,ylm,ylm_k)
         do ilm=1,mpsang*mpsang
           do ipw=1,npw_k
             ylm_k(ipw,ilm)=ylm(ipw+ikg,ilm)
           end do
         end do
!        $OMP END PARALLEL DO
       end if

!      End if for choice governed by mkmem
     end if

     cplex=2;if (istwf_k>1) cplex=1

!    Extract PAW cprj quantities according to mkmem
     if (mkmem==0.and.usecprj==1) then
       ABI_ALLOCATE(cprj_disk,(natom,nspinor*nband_k))
       call cprj_alloc(cprj_disk,3,dimcprj)
       call cprj_get(atindx1,cprj_disk,cprj,natom,1,ibg,ikpt,iorder_cprj,isppol,&
&       mband,mkmem,mpi_enreg,natom,nband_k,nband_k,nspinor,nsppol,unpaw)
     end if

!    Compute nonlocal psp energy

!    Compute (k+G) vectors (only if useylm=1)
     nkpg=9*nloalg(5)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

!    Compute nonlocal form factors ffnl at all (k+G):
     ider=0;idir=0;dimffnl=1
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&     gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,&
&     ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)

     dyfrnlk(:,:)=zero
     if (pawbec==1) becfrnl(:,:,:)=zero

!    Compute phkxred and eventually ph3d.
     do iatom=1,natom
       ia=atindx1(iatom)
       arg=two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
       phkxred(1,iatom)=cos(arg)
       phkxred(2,iatom)=sin(arg)
     end do
     if(nloalg(1)<=0)then
!      Only the allocation, not the precomputation.
       matblk=nloalg(4)
       ABI_ALLOCATE(ph3d,(2,npw_k,matblk))
     else
!      Here, allocation as well as precomputation
       matblk=natom
       ABI_ALLOCATE(ph3d,(2,npw_k,matblk))
       call ph1d3d(1,natom,kg_k,matblk,natom,npw_k,n1,n2,n3,&
&       phkxred,ph1d,ph3d)
     end if

     do iband=1,nband_k

       if(mpi_enreg%paral_compil_kpt==1)then
         if(mpi_enreg%proc_distrb(ikpt, iband,isppol) /= me) then
           cycle
         end if
       end if

       occ_k=occ(iband+bdtot_index)

       if(mkmem/=0)then
         cwavef(:,1:npw_k*nspinor)=&
&         cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
       else
         cwavef(:,1:npw_k*nspinor)=&
&         cg_disk(:,1+(iband-1)*npw_k*nspinor:iband*npw_k*nspinor)
       end if

       if (psps%usepaw==1.and.usecprj==1) then
         if (mkmem/=0) then
           ibsp=(iband-1)*nspinor+ibg
           call cprj_copy(cprj     (:,ibsp+1:ibsp+nspinor),cwaveprj)
         else
           ibsp=(iband-1)*nspinor
           call cprj_copy(cprj_disk(:,ibsp+1:ibsp+nspinor),cwaveprj)
         end if
       end if

!      Compute non-local contributions from n,k
       if (psps%usepaw==1) eig_k=eigen(iband+bdtot_index)
       call nonlop(atindx1,choice,cpopt,cwaveprj,dimekb1,dimekb2,dimffnl,dimffnl,ekb,&
&       enlout,ffnl,ffnl,gmet,gprimd,idir,psps%indlmn,istwf_k,kg_k,kg_k,kpg_k,kpg_k,kpoint,&
&       kpoint,eig_k,psps%lmnmax,matblk,mgfft,mpi_enreg,mpsang,psps%mpssoang,natom,nattyp,&
&       ngfft,nkpg,nkpg,nloalg,nnlout,npw_k,npw_k,nspinor,nspinor,ntypat,0,paw_opt,phkxred,&
&       phkxred,ph1d,ph3d,ph3d,signs,sij,nonlop_dum,tim_nonlop,ucvol,&
&       psps%useylm,cwavef,cwavef)

!      Accumulate non-local contributions from n,k
       dyfrnlk(:,:)=dyfrnlk(:,:)+occ_k*reshape(enlout(:),(/6,natom/))

!      PAW: accumulate gradients of rhoij
       if (psps%usepaw==1) then
         call pawaccrhoij(atindx1,cplex,cwaveprj,cwaveprj,0,isppol,natom,&
&         nspinor,occ_k,option_rhoij,pawrhoij,usetimerev,wtk_k)
       end if

!      PAW: Compute frozen contribution to Born Effective Charges
       if (pawbec==1) then
         do ii=1,3 ! Loop over elect. field directions
           call nonlop(atindx1,choice_bec,cpopt_bec,cwaveprj,psps%dimekb,dimekb2,dimffnl,dimffnl,&
&           becij(:,:,:,ii),enlout_bec,ffnl,ffnl,gmet,gprimd,idir,psps%indlmn,istwf_k,kg_k,kg_k,&
&           kpg_k,kpg_k,kpoint,kpoint,eig_k,psps%lmnmax,matblk,mgfft,mpi_enreg,mpsang,psps%mpssoang,&
&           natom,nattyp,ngfft,nkpg,nkpg,nloalg,nnlout_bec,npw_k,npw_k,nspinor,nspinor,ntypat,0,paw_opt_bec,phkxred,&
&           phkxred,ph1d,ph3d,ph3d,signs,sij,nonlop_dum,tim_nonlop,ucvol,&
&           psps%useylm,cwavef,cwavef)
           becfrnl(:,:,ii)=becfrnl(:,:,ii)+occ_k*reshape(enlout_bec(:),(/3,natom/))
         end do
       end if

!      End of loop on bands
     end do

     do iatom=1,natom
       ia=iatom;if (dyfr_nondiag==0) ia=1
       dyfrnl(1,1,1,iatom,ia)=dyfrnl(1,1,1,iatom,ia)+wtk_k*dyfrnlk(1,iatom)/ucvol
       dyfrnl(1,2,2,iatom,ia)=dyfrnl(1,2,2,iatom,ia)+wtk_k*dyfrnlk(2,iatom)/ucvol
       dyfrnl(1,3,3,iatom,ia)=dyfrnl(1,3,3,iatom,ia)+wtk_k*dyfrnlk(3,iatom)/ucvol
       dyfrnl(1,2,3,iatom,ia)=dyfrnl(1,2,3,iatom,ia)+wtk_k*dyfrnlk(4,iatom)/ucvol
       dyfrnl(1,1,3,iatom,ia)=dyfrnl(1,1,3,iatom,ia)+wtk_k*dyfrnlk(5,iatom)/ucvol
       dyfrnl(1,1,2,iatom,ia)=dyfrnl(1,1,2,iatom,ia)+wtk_k*dyfrnlk(6,iatom)/ucvol
     end do

!    Incremente indexes
     bdtot_index=bdtot_index+nband_k
     if (mkmem/=0) then
       ibg=ibg+nband_k*nspinor
       icg=icg+npw_k*nspinor*nband_k
       ikg=ikg+npw_k
     else if (usecprj==1) then
       call cprj_free(cprj_disk)
       ABI_DEALLOCATE(cprj_disk)
     end if

     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(ylm_k)

!    End loops on isppol and ikpt
   end do
 end do

 ABI_DEALLOCATE(phkxred)
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(dyfrnlk)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(enlout)
 ABI_DEALLOCATE(ekb)
 if (pawbec==1)  then
   ABI_DEALLOCATE(enlout_bec)
   ABI_DEALLOCATE(becij)
 end if
 if (psps%usepaw==1) then
   call cprj_free(cwaveprj)
   ABI_DEALLOCATE(cwaveprj)
   ABI_DEALLOCATE(sij)
 end if
 if(mkmem==0)  then
   ABI_DEALLOCATE(cg_disk)
 end if

 do iatom=1,natom
   ia=iatom;if (dyfr_nondiag==0) ia=1
   dyfrnl(1,3,2,iatom,ia)=dyfrnl(1,2,3,iatom,ia)
   dyfrnl(1,3,1,iatom,ia)=dyfrnl(1,1,3,iatom,ia)
   dyfrnl(1,2,1,iatom,ia)=dyfrnl(1,1,2,iatom,ia)
 end do

!Parallel case: accumulate (n,k) contributions
 if( mpi_enreg%paral_compil_kpt==1) then
   call timab(48,1,tsec)
   call leave_test()
!  Accumulate dyfrnl
   call xsum_mpi(dyfrnl(1,:,:,:,:),spaceworld,ierr)
!  Accumulate becfrnl
   if (pawbec==1) then
     call xsum_mpi(becfrnl(:,:,:),spaceworld,ierr)
   end if
   if (psps%usepaw==1) then
!    PAW: accumulate gradients of rhoij
     ABI_ALLOCATE(dimlmn,(natom))
     dimlmn(1:natom)=pawrhoij(1:natom)%cplex*pawrhoij(1:natom)%lmn2_size
     bufdim=3*sum(dimlmn)*nsploop
     ABI_ALLOCATE(buffer1,(bufdim))
     ABI_ALLOCATE(buffer2,(bufdim))
     ii=0
     do iatom=1,natom
       do isppol=1,nsploop
         do mu=1,3
           buffer1(ii+1:ii+dimlmn(iatom))=pawrhoij(iatom)%grhoij(mu,1:dimlmn(iatom),isppol)
           ii=ii+dimlmn(iatom)
         end do
       end do
     end do
     call xsum_mpi(buffer1,buffer2,bufdim,spaceworld,ierr)
     ii=0
     do iatom=1,natom
       do isppol=1,nsploop
         do mu=1,3
           pawrhoij(iatom)%grhoij(mu,1:dimlmn(iatom),isppol)=buffer2(ii+1:ii+dimlmn(iatom))
           ii=ii+dimlmn(iatom)
         end do
       end do
     end do
     ABI_DEALLOCATE(buffer1)
     ABI_DEALLOCATE(buffer2)
     ABI_DEALLOCATE(dimlmn)
   end if
   call timab(48,2,tsec)
 end if

 if (psps%usepaw==1) then
!  PAW: symmetrize rhoij gradients and transfer to cartesian (reciprocal space) coord.
   choice=2;option=0  ! This symetrization is necessary in the antiferromagnetic case...
   call symrhoij(choice,gprimd,psps%indlmn,indsym,0,psps%lmnmax,natom,natom,nsym,ntypat,option,&
&   pawang,pawprtvol,pawrhoij,rprimd,symafm,symrec,typat)
   do iatom=1,natom
     cplx=pawrhoij(iatom)%cplex
     do isppol=1,nsploop
       do klmn=1,pawrhoij(iatom)%lmn2_size
         do ii=1,cplx
           grhoij(1:3)=pawrhoij(iatom)%grhoij(1:3,cplx*(klmn-1)+ii,isppol)
           do mu=1,3
             pawrhoij(iatom)%grhoij(mu,cplx*(klmn-1)+ii,isppol)=gprimd(mu,1)*grhoij(1) &
&             +gprimd(mu,2)*grhoij(2)+gprimd(mu,3)*grhoij(3)
           end do
         end do
       end do
     end do
   end do

!  PAW: Add gradients due to Dij derivatives to dynamical matrix
   if (usexcnhat==0) then
     dimvtmp=1
     ABI_ALLOCATE(vtmp,(nfftf,dimvtmp))
     vtmp(:,1)=vtrial(:,1)-vxc(:,1)
   else
     dimvtmp=nspden
     vtmp => vtrial
   end if
   dimnhat=0
   optgr=0
   optgr2=1
   optstr=0
   ABI_ALLOCATE(nhat_dum,(1,0))
   call pawgrnl(atindx1,dimnhat,dimvtmp,dyfrnl,dyfr_cplex,dummy,gsqcut,mgfftf,mpi_enreg,natom,natom,&
&   nattyp,nfftf,ngfftf,nhat_dum,dummy,nspden,nsym,ntypat,optgr,optgr2,optstr,paral_kgb,&
&   pawang,pawfgrtab,pawrhoij,pawtab,ph1df,psps,qphon,rprimd,symrec,typat,vtmp,xred)
   ABI_DEALLOCATE(nhat_dum)
   if (usexcnhat==0)  then
     ABI_DEALLOCATE(vtmp)
   end if

 end if

!Born Effective Charges and PAW:
!1-Re-order atoms -- 2-Add contribution from rhoij
 if (pawbec==1) then
   ABI_ALLOCATE(becfrnl_tmp,(3,natom,3))
   becfrnl_tmp=becfrnl
   do ia=1,natom         ! Atom (sorted by type)
     iatom=atindx1(ia)   ! Atom (not sorted)
     arg=zero            ! Computation of Sum [Rhoij.Sij]
     do isp=1,nsppol
       do mu=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(mu)
         arg=arg+pawrhoij(iatom)%rhoijp(mu,isp)*pawtab(itypat)%sij(klmn)
       end do
     end do
     do ii=1,3           ! Direction of electric field
       do jj=1,3         ! Direction of atom
         becfrnl(jj,iatom,ii)=becfrnl_tmp(jj,ia,ii)+arg*gmet(ii,jj)
       end do
     end do
   end do
   ABI_DEALLOCATE(becfrnl_tmp)
!  TEMPORARY - for testing purpose
   becfrnl=zero
 end if

 call timab(159,2,tsec)

 DBG_EXIT("COLL")

end subroutine dyfnl3
!!***
