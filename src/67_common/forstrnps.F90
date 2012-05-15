!{\src2tex{textfont=tt}}
!!****f* ABINIT/forstrnps
!! NAME
!! forstrnps
!!
!! FUNCTION
!! Compute nonlocal pseudopotential energy contribution to forces and/or stress tensor
!! as well as kinetic energy contribution to stress tensor.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, AF, AR, MB, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cg(2,mcg)=wavefunctions (may be read from disk file)
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!       and electric field calculations (see initberry.f).
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass=effective mass for electrons (1. in common case)
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced coordinates (integers) of G vecs in basis
!!  kpt(3,nkpt)=k points in reduced coordinates
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=maximum number of k points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw= maximum number of plane waves
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nband(nkpt)=number of bands at each k point
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points in Brillouin zone
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt)=number of planewaves in basis and boundary at each k
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of elements in symmetry group
!!  ntypat=number of types of atoms
!!  occ(mband*nkpt*nsppol)=occupation numbers for each band over all k points
!!  optfor=1 if computation of forces is required
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  stress_needed=1 if computation of stress tensor is required
!!  symrec(3,3,nsym)=symmetries in reciprocal space (dimensionless)
!!  typat(natom)=type integer for each atom in cell
!!  unkg=unit number for (k+G) data (if used)
!!  unylm=unit number for Ylm(k) data (if used)
!!  use_gpu_cuda= 0 or 1 to know if we use cuda for nonlop call
!!  wffnow=unit number of disk file for wf if used
!!  wtk(nkpt)=weight associated with each k point
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  if (optfor==1)
!!   grnl(3*natom)=stores grads of nonlocal energy wrt atomic coordinates
!!  if (stress_needed==1)
!!   kinstr(6)=kinetic energy part of stress tensor (hartree/bohr^3)
!!   Store 6 unique components of symmetric 3x3 tensor in the order
!!   11, 22, 33, 32, 31, 21
!!   npsstr(6)=nonlocal pseudopotential energy part of stress tensor
!!    (hartree/bohr^3)
!!
!! PARENTS
!!      forstr
!!
!! CHILDREN
!!      hdr_skip,leave_test,meanvalue_g,metric,mkffnl,mkkpg,nonlop,ph1d3d
!!      prep_bandfft_tabs,prep_nonlop,rdnpw,rwwf,stresssym,timab,wrtout
!!      xcomm_init,xdefineoff,xmaster_init,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine forstrnps(atindx1,cg,dtefield,ecut,ecutsm,effmass,eigen,electronpositron,&
&  grnl,istwfk,kg,kinstr,npsstr,kpt,mband,mcg,mgfft,mkmem,mpi_enreg,mpsang,&
&  mpw,natom,nattyp,nband,ngfft,nkpt,nloalg,npwarr,nspinor,nsppol,nsym,&
&  ntypat,occ,optfor,paw_ij,pawtab,ph1d,psps,rprimd,&
&  stress_needed,symrec,typat,unkg,unylm,use_gpu_cuda,wffnow,wtk,xred,ylm,ylmgr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_efield
 use m_xmpi
 use m_wffile

 use m_electronpositron, only : electronpositron_type,electronpositron_calctype

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'forstrnps'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_spacepar
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_65_nonlocal
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mgfft,mkmem,mpsang,mpw,natom,nkpt
 integer,intent(in) :: nsppol,nspinor,nsym,ntypat,optfor,stress_needed
 integer,intent(in) :: unkg,unylm,use_gpu_cuda
 real(dp),intent(in) :: ecut,ecutsm,effmass
 type(efield_type), intent(in) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow
!arrays
 integer,intent(in) :: atindx1(natom),istwfk(nkpt)
 integer,intent(in) :: kg(3,mpw*mkmem),nattyp(ntypat),nband(nkpt*nsppol)
 integer,intent(in) :: ngfft(18),nloalg(5),npwarr(nkpt)
 integer,intent(in) :: symrec(3,3,nsym),typat(natom)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),kpt(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rprimd(3,3),wtk(nkpt),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang*psps%useylm)
 real(dp),intent(out) :: grnl(3*natom),kinstr(6),npsstr(6)
 type(paw_ij_type),intent(in) :: paw_ij(natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: bdtot_index,choice,cpopt,dimdij,dimekb1,dimekb2,dimffnl
 integer :: formeig,ia,iatom,iband,icg,ider,idir,ierr,ii,ikg,ikpt,ilm,ilmn
 integer :: ipositron,ipw,ishift,isp,ispden,isppol,istwf_k,itypat
 integer :: master,matblk,mcg_disk,me_distrb,my_nspinor,n1,n2,n3,nband_k,nkpg
 integer :: nnlout,npw_k,nspinor_,paw_opt,signs,spaceComm,tim_nonlop,tim_nonlop_prep
 integer :: tim_rwwf
 real(dp) :: ar,arg,renorm_factor,dfsm,ecutsm_inv,eig_k,fact_kin,fsm,htpisq,kgc1
 real(dp) :: kgc2,kgc3,kin,ucvol,xx
 character(len=500) :: message
!arrays
 integer,allocatable :: kg_dum(:,:),kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3)
 real(dp) :: nonlop_dum(1,1),rmet(3,3),tsec(2)
 real(dp),allocatable :: cg_disk(:,:),cwavef(:,:),eig_dum(:),ekb(:,:,:)
 real(dp),allocatable :: enlout(:),ffnl(:,:,:,:),kpg_k(:,:),kstr1(:),kstr2(:)
 real(dp),allocatable :: kstr3(:),kstr4(:),kstr5(:),kstr6(:),occ_dum(:)
 real(dp),allocatable :: ph3d(:,:,:),phkxred(:,:),sij(:,:),ylm_k(:,:)
 real(dp),allocatable :: ylmgr_k(:,:,:)
 type(cprj_type) :: cprj_dum(1,1)
!BandFFT parallelization
 integer :: blocksize,iblock,iblocksize,ibs,nblockbd
 integer :: gbound_dum(1,1)
 real(dp) :: kinpw_dum(1,1)
 real(dp),allocatable :: lambda_loc(:),occblock(:),weight(:)

!*************************************************************************

!DEBUG
!write(std_out,*)' forstrnps : enter '
!if(.true.)stop
!ENDDEBUG

 call timab(920,1,tsec)
 call timab(921,1,tsec)

!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
!Init me
 call xme_init(mpi_enreg,me_distrb)

!Init master
 call xmaster_init(mpi_enreg,master)

!Some constants
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 ipositron=abs(electronpositron_calctype(electronpositron))
 my_nspinor=max(1,nspinor/mpi_enreg%nproc_spin)
!Smearing of plane wave kinetic energy
 ecutsm_inv=zero;if( ecutsm>1.0d-20) ecutsm_inv=1/ecutsm
!htpisq is (1/2) (2 Pi) **2:
 htpisq=0.5_dp*(two_pi)**2

!Arrays initializations
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(phkxred,(2,natom))
 if (optfor==1) grnl(:)=zero
 if (stress_needed==1) then
   kinstr(:)=zero;npsstr(:)=zero
 end if

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Determine whether wf are being read from disk file or not
 if (mkmem==0) then
!  Skip wavefunction file header
   call hdr_skip(wffnow,ierr)
   mcg_disk=mpw*my_nspinor*mband
   ABI_ALLOCATE(cg_disk,(2,mcg_disk))
!  Define offsets, in case of MPI I/O
   formeig=0
   call xdefineOff(formeig,wffnow,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)
 end if

!Common data for "nonlop" routine
 signs=1 ; idir=0  ; ishift=0 ; tim_nonlop=4 ; tim_nonlop_prep=12
 choice=2*optfor;if (stress_needed==1) choice=10*choice+3
 if (optfor==1.and.stress_needed==1)  ishift=6
 nnlout=max(1,6*stress_needed+3*natom*optfor)
 if (psps%usepaw==0) then
   paw_opt=0 ; cpopt=-1
 else
   paw_opt=2 ; cpopt=-1
 end if

!Non-local factors:
!Norm-conserving: kleimann-Bylander energies
!PAW: Dij coefficients and overlap coefficients
 if (psps%usepaw==0) then
   dimekb1=psps%dimekb;dimekb2=ntypat
   ABI_ALLOCATE(ekb,(psps%dimekb,ntypat,nspinor**2))
   ekb(:,:,1)=psps%ekb(:,:)
   if (nspinor==2) then
     ekb(:,:,2)=psps%ekb(:,:)
     ekb(:,:,3:4)=zero
   end if
   if (ipositron==1) ekb(:,:,:)=-ekb(:,:,:)
 else
   dimekb1=psps%dimekb*paw_ij(1)%cplex_dij;dimekb2=natom
   ABI_ALLOCATE(ekb,(dimekb1,dimekb2,nspinor**2))
   ABI_ALLOCATE(sij,(dimekb1,ntypat))
   do itypat=1,ntypat
     if (paw_ij(1)%cplex_dij==1) then
       sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
     else
       do ilmn=1,pawtab(itypat)%lmn2_size
         sij(2*ilmn-1,itypat)=pawtab(itypat)%sij(ilmn)
         sij(2*ilmn  ,itypat)=zero
       end do
     end if
   end do
 end if

!DEBUG ! Do not remove this line : needed for the gfortran compiler ?!
!write(std_out,*)' forstrnps : usepaw=',psps%usepaw
!ENDDEBUG

 call timab(921,2,tsec)

!LOOP OVER SPINS
 bdtot_index=0;icg=0
 do isppol=1,nsppol

   call timab(927,1,tsec)

   if (nsppol==2) then
     write(message, '(a,i3)' )' ****  In forstrnps for isppol=',isppol
     call wrtout(std_out,message,'COLL')
   end if

!  Rewind temporary disk files
   if (mkmem==0) rewind unkg
   if (mkmem==0.and.psps%useylm==1) rewind unylm

!  PAW: retrieve Dij coefficients for this spin component
   if (psps%usepaw==1) then
     do ispden=1,nspinor**2
       do iatom=1,natom
         isp=isppol;if (nspinor==2) isp=ispden
         dimdij=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
         do ilmn=1,dimdij
           ekb(ilmn,iatom,ispden)=paw_ij(iatom)%dij(ilmn,isp)
         end do
         if(dimdij+1<=dimekb1) ekb(dimdij+1:dimekb1,iatom,ispden)=zero
       end do
     end do
     if (dtefield%berryopt==4) then ! finite field with paw has on-site dipole force
       do ispden=1,nspinor**2
         do iatom=1,natom
           itypat = typat(iatom)
           isp=isppol;if (nspinor==2) isp=ispden
           dimdij=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
           do idir = 1, 3
             if (abs(dtefield%efield_dot(idir)) < tol12) cycle
             do ilmn=1,dimdij
               ekb(ilmn,iatom,ispden)=ekb(ilmn,iatom,ispden)-&
&               dtefield%efield_dot(idir)*dtefield%rij(ilmn,itypat,idir)/two_pi
             end do
           end do ! end loop over directions of e field
           if(dimdij+1<=dimekb1) ekb(dimdij+1:dimekb1,iatom,ispden)=zero
         end do
       end do
     end if ! end on-site dipole term
   end if ! end PAW Dij retrieval

   call timab(927,2,tsec)

!  Loop over k points
   ikg=0
   do ikpt=1,nkpt

     nband_k=nband(ikpt+(isppol-1)*nkpt)
     istwf_k=istwfk(ikpt)
     npw_k=npwarr(ikpt)

     if(mpi_enreg%paral_compil_kpt==1)then
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me_distrb))/=0) then
         bdtot_index=bdtot_index+nband_k
         cycle
       end if
     end if ! mpi_enreg%paral_compil_kpt

     call timab(922,1,tsec)

!    Define several block values for band parallelization
     if (mpi_enreg%mode_para/='b') then
       nblockbd=nband_k
       blocksize=1
     else
       nblockbd=nband_k/mpi_enreg%nproc_fft
       if (nband_k/=nblockbd*mpi_enreg%nproc_fft) nblockbd=nblockbd+1
       if(mpi_enreg%mode_para=='b') then
         nblockbd=nband_k/(mpi_enreg%nproc_band*mpi_enreg%bandpp)
       end if
       blocksize=nband_k/nblockbd
     end if
     ABI_ALLOCATE(enlout,(nnlout*blocksize))
     ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor*blocksize))

     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     if (stress_needed==1) then
       if (psps%useylm==1)  then
         ABI_ALLOCATE(ylmgr_k,(npw_k,3,mpsang*mpsang*psps%useylm))
       end if
       ABI_ALLOCATE(kstr1,(npw_k))
       ABI_ALLOCATE(kstr2,(npw_k))
       ABI_ALLOCATE(kstr3,(npw_k))
       ABI_ALLOCATE(kstr4,(npw_k))
       ABI_ALLOCATE(kstr5,(npw_k))
       ABI_ALLOCATE(kstr6,(npw_k))
     end if

     kpoint(:)=kpt(:,ikpt)

     kg_k(:,:) = 0
     if (mkmem==0) then

       nspinor_=nspinor
       call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor_,0,unkg)
!      Skip sphere data centered at k in unkg, then read k+g data
       read (unkg) ((kg_k(ii,ipw),ii=1,3),ipw=1,npw_k)

!      Read the wavefunction block for ikpt,isppol
       tim_rwwf=7
       ABI_ALLOCATE(eig_dum,(mband))
       ABI_ALLOCATE(kg_dum,(3,0))
       ABI_ALLOCATE(occ_dum,(mband))
       call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,mband,mcg_disk,mpi_enreg,nband_k,nband_k,&
&       npw_k,my_nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
       ABI_DEALLOCATE(eig_dum)
       ABI_DEALLOCATE(kg_dum)
       ABI_DEALLOCATE(occ_dum)

!      Eventually read spherical harmonics
       if (psps%useylm==1) then
         read(unylm)
         if (stress_needed==1) then
           read(unylm) ((ylm_k(ipw,ilm),ipw=1,npw_k),ilm=1,mpsang*mpsang),&
&           (((ylmgr_k(ipw,ii,ilm),ipw=1,npw_k),ii=1,3),ilm=1,mpsang*mpsang)
         else
           read(unylm) ((ylm_k(ipw,ilm),ipw=1,npw_k),ilm=1,mpsang*mpsang)
         end if
       end if

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
!        $OMP PARALLEL DO PRIVATE(ilm,ipw) &
!        $OMP&SHARED(ikg,mpsang,npw_k,ylm,ylm_k)
         do ilm=1,mpsang*mpsang
           do ipw=1,npw_k
             ylm_k(ipw,ilm)=ylm(ipw+ikg,ilm)
           end do
         end do
!        $OMP END PARALLEL DO
         if (stress_needed==1) then
!          $OMP PARALLEL DO PRIVATE(ilm,ipw) &
!          $OMP&SHARED(ikg,mpsang,npw_k,ylmgr,ylmgr_k)
           do ilm=1,mpsang*mpsang
             do ii=1,3
               do ipw=1,npw_k
                 ylmgr_k(ipw,ii,ilm)=ylmgr(ipw+ikg,ii,ilm)
               end do
             end do
           end do
!          $OMP END PARALLEL DO
         end if
       end if

!      End if for choice governed by mkmem
     end if

!    Prepare kinetic contribution to stress tensor (Warning : the symmetry
!    has not been broken, like in mkkin.f or kpg3.f . It should
!    be, in order to be coherent).
!    $OMP PARALLEL DO PRIVATE(fact_kin,ipw,kgc1,kgc2,kgc3,kin,xx) &
!    $OMP&SHARED(ecut,ecutsm,ecutsm_inv) &
!    $OMP&SHARED(gprimd,htpisq,kg_k,kpoint,kstr1,kstr2,kstr3,kstr4,kstr5,kstr6,npw_k)
     if (stress_needed==1) then
       do ipw=1,npw_k
!        Compute Cartesian coordinates of (k+G)
         kgc1=gprimd(1,1)*(kpoint(1)+kg_k(1,ipw))+&
&         gprimd(1,2)*(kpoint(2)+kg_k(2,ipw))+&
&         gprimd(1,3)*(kpoint(3)+kg_k(3,ipw))
         kgc2=gprimd(2,1)*(kpoint(1)+kg_k(1,ipw))+&
&         gprimd(2,2)*(kpoint(2)+kg_k(2,ipw))+&
&         gprimd(2,3)*(kpoint(3)+kg_k(3,ipw))
         kgc3=gprimd(3,1)*(kpoint(1)+kg_k(1,ipw))+&
&         gprimd(3,2)*(kpoint(2)+kg_k(2,ipw))+&
&         gprimd(3,3)*(kpoint(3)+kg_k(3,ipw))
         kin=htpisq* ( kgc1**2 + kgc2**2 + kgc3**2 )
         fact_kin=1.0_dp
         if(kin>ecut-ecutsm)then
           if(kin>ecut)then
             fact_kin=0.0_dp
           else
!            See the routine mkkin.f, for the smearing procedure
             xx=(ecut-kin)*ecutsm_inv
!            This kinetic cutoff smoothing function and its xx derivatives
!            were produced with Mathematica and the fortran code has been
!            numerically checked against Mathematica.
             fsm=1.0_dp/(xx**2*(3+xx*(1+xx*(-6+3*xx))))
             dfsm=-3.0_dp*(-1+xx)**2*xx*(2+5*xx)*fsm**2
!            d2fsm=6.0_dp*xx**2*(9+xx*(8+xx*(-52+xx*(-3+xx*(137+xx*&
!            &                         (-144+45*xx))))))*fsm**3
             fact_kin=fsm+kin*(-ecutsm_inv)*dfsm
           end if
         end if
         kstr1(ipw)=fact_kin*kgc1*kgc1
         kstr2(ipw)=fact_kin*kgc2*kgc2
         kstr3(ipw)=fact_kin*kgc3*kgc3
         kstr4(ipw)=fact_kin*kgc3*kgc2
         kstr5(ipw)=fact_kin*kgc3*kgc1
         kstr6(ipw)=fact_kin*kgc2*kgc1
       end do ! ipw
!      $OMP END PARALLEL DO
     end if

!    Compute (k+G) vectors (only if useylm=1)
     nkpg=3*nloalg(5)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

!    Compute nonlocal form factors ffnl at all (k+G)
!    (ider=1 computes gradients needed for stress tensor)
     ider=0;idir=0;dimffnl=1
     if (stress_needed==1) then
       ider=1;dimffnl=2+2*psps%useylm
     end if
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&     gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&     npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&     psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

!    Allocate the arrays phkxred and ph3d, compute phkxred and eventually ph3d
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

     call timab(922,2,tsec)

!    Loop over bands; accumulate forces and/or stresses
!    Note that in sequential mode iblock=iband, nblockbd=nband_k and blocksize=1
     ABI_ALLOCATE(occblock,(blocksize))
     ABI_ALLOCATE(weight,(blocksize))
     do iblock=1,nblockbd

       iband=(iblock-1)*blocksize+1
       if(mpi_enreg%paral_compil_kpt==1)then
!        Skip this band if not the proper processor
         if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/= me_distrb) cycle
       end if

!      Select occupied bands
       occblock(:)=occ(1+(iblock-1)*blocksize+bdtot_index:iblock*blocksize+bdtot_index)
       if( abs(maxval(occblock))>=tol8 ) then
         call timab(923,1,tsec)
         weight(:)=wtk(ikpt)*occblock(:)

!        Load contribution from n,k
         if(mkmem/=0)cwavef(:,1:npw_k*my_nspinor*blocksize)=&
&         cg(:,1+(iblock-1)*npw_k*my_nspinor*blocksize+icg:iblock*npw_k*my_nspinor*blocksize+icg)
         if(mkmem==0)cwavef(:,1:npw_k*my_nspinor)=&
&         cg_disk(:,1+(iblock-1)*npw_k*my_nspinor:iblock*npw_k*my_nspinor)

         call timab(923,2,tsec)
         call timab(926,1,tsec)

         if (mpi_enreg%mode_para/='b') then
           if (psps%usepaw==1) eig_k=eigen(iblock+bdtot_index)
           call nonlop(atindx1,choice,cpopt,cprj_dum,dimekb1,dimekb2,dimffnl,dimffnl,ekb,&
&           enlout,ffnl,ffnl,gmet,gprimd,idir,psps%indlmn,istwf_k,kg_k,kg_k,kpg_k,kpg_k,kpoint,&
&           kpoint,eig_k,psps%lmnmax,matblk,mgfft,mpi_enreg,mpsang,psps%mpssoang,natom,&
&           nattyp,ngfft,nkpg,nkpg,nloalg,nnlout,npw_k,npw_k,my_nspinor,nspinor,ntypat,0,paw_opt,&
&           phkxred,phkxred,ph1d,ph3d,ph3d,signs,sij,nonlop_dum,tim_nonlop,&
&           ucvol,psps%useylm,cwavef,cwavef,use_gpu_cuda=use_gpu_cuda)
         else
!          Transpose the ffnl, kinpw, kpg and ph3d arrays.
           call prep_bandfft_tabs(dimffnl,ffnl,gbound_dum,ikpt,kinpw_dum,kpoint,&
&           psps%lmnmax,matblk,mgfft,mkmem,mpi_enreg,nkpg,npw_k,ntypat,0,ph3d)
           ABI_ALLOCATE(lambda_loc,(blocksize))
           if (psps%usepaw==1)lambda_loc(1:blocksize)=&
&           eigen(1+(iblock-1)*blocksize+bdtot_index:iblock*blocksize+bdtot_index)
           call prep_nonlop(atindx1,choice,cpopt,cprj_dum,dimekb1,dimekb2,dimffnl,ekb,&
&           enlout,gmet,gprimd,idir,ikpt,psps%indlmn,istwf_k,&
&           kpoint,lambda_loc,psps%lmnmax,matblk,blocksize,mgfft,mpi_enreg,mpsang,psps%mpssoang,natom,&
&           nattyp,ngfft,nkpg,nloalg,nnlout,npw_k,my_nspinor,nspinor,ntypat,paw_opt,&
&           phkxred,ph1d,signs,sij,nonlop_dum,tim_nonlop_prep,&
&           ucvol,psps%useylm,cwavef,cwavef,use_gpu_cuda=use_gpu_cuda)
           ABI_DEALLOCATE(lambda_loc)
         end if

         call timab(926,2,tsec)

!        Accumulate non-local contributions from n,k
         if (optfor==1) then
           do iblocksize=1,blocksize
             ibs=nnlout*(iblocksize-1)
             grnl(1:3*natom)=grnl(1:3*natom)+weight(iblocksize)*enlout(ibs+1+ishift:ibs+3*natom+ishift)
           end do
         end if
         if (stress_needed==1) then
           do iblocksize=1,blocksize
             ibs=nnlout*(iblocksize-1)
             npsstr(1:6)    =npsstr(1:6)    +weight(iblocksize)*enlout(ibs+1       :ibs+6)
           end do
         end if

!        Accumulate stress tensor kinetic contributions
         if (stress_needed==1) then
           call timab(924,1,tsec)
           do iblocksize=1,blocksize
             call meanvalue_g(ar,kstr1,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),0)
             kinstr(1)=kinstr(1)+weight(iblocksize)*ar
             call meanvalue_g(ar,kstr2,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),0)
             kinstr(2)=kinstr(2)+weight(iblocksize)*ar
             call meanvalue_g(ar,kstr3,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),0)
             kinstr(3)=kinstr(3)+weight(iblocksize)*ar
             call meanvalue_g(ar,kstr4,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),0)
             kinstr(4)=kinstr(4)+weight(iblocksize)*ar
             call meanvalue_g(ar,kstr5,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),0)
             kinstr(5)=kinstr(5)+weight(iblocksize)*ar
             call meanvalue_g(ar,kstr6,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),0)
             kinstr(6)=kinstr(6)+weight(iblocksize)*ar
           end do
           call timab(924,2,tsec)
         end if

!        End of loop on block of bands
       end if
     end do

     ABI_DEALLOCATE(occblock)
     ABI_DEALLOCATE(weight)
     ABI_DEALLOCATE(enlout)
     ABI_DEALLOCATE(cwavef)

!    Incremente indexes
     bdtot_index=bdtot_index+nband_k
     if (mkmem/=0) then
       icg=icg+npw_k*my_nspinor*nband_k
       ikg=ikg+npw_k
     end if

     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(ylm_k)
     if (stress_needed==1) then
       ABI_DEALLOCATE(kstr1)
       ABI_DEALLOCATE(kstr2)
       ABI_DEALLOCATE(kstr3)
       ABI_DEALLOCATE(kstr4)
       ABI_DEALLOCATE(kstr5)
       ABI_DEALLOCATE(kstr6)
       if (psps%useylm==1)  then
         ABI_DEALLOCATE(ylmgr_k)
       end if
     end if

!    End k point loop
   end do
!  End loop over spins
 end do

 if(mkmem==0) then
   ABI_DEALLOCATE(cg_disk)
 end if

!Parallel case: accumulate (n,k) contributions
 if( mpi_enreg%paral_compil_kpt==1) then
   write(message, '(a)' ) 'forstrnps: loop on k-points and spins done in parallel'
   call wrtout(std_out,message,'COLL')
!  Forces
   if (optfor==1) then
     call timab(65,1,tsec)
     call leave_test()

     if ((mpi_enreg%paral_compil_kpt==1) .and. &
&     (mpi_enreg%paral_compil_fft==1)) then
       call xsum_mpi(grnl,mpi_enreg%comm_kpt,ierr)
     else
       call xsum_mpi(grnl,spaceComm,ierr)
     end if


     call timab(65,2,tsec)
   end if
!  Stresses
   if (stress_needed==1) then
     call timab(65,1,tsec)
     call leave_test()

!    PATCH forstrnps // KPT & FFT spacecomm --> comm_kpt
     if ((mpi_enreg%paral_compil_kpt==1) .and. &
&     (mpi_enreg%paral_compil_fft==1)) then
       call xsum_mpi(kinstr,mpi_enreg%comm_kpt,ierr)
       call xsum_mpi(npsstr,mpi_enreg%comm_kpt,ierr)
     else
       call xsum_mpi(kinstr,spaceComm,ierr)
       call xsum_mpi(npsstr,spaceComm,ierr)
     end if
     call timab(65,2,tsec)
   end if
 end if

 call timab(925,1,tsec)

!Deallocate temporary space
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(phkxred)
 ABI_DEALLOCATE(ekb)
 if (psps%usepaw==1)  then
   ABI_DEALLOCATE(sij)
 end if

!Do final normalisation of nl contribution to forces
 if (optfor==1) grnl(:)=grnl(:)/ucvol

!Do final normalizations and symmetrizations of stress tensor contributions
 if (stress_needed==1) then
   renorm_factor=-(two_pi**2)/effmass/ucvol
   kinstr(:)=kinstr(:)*renorm_factor
   if (nsym>1) then
     call stresssym(gprimd,nsym,kinstr,symrec)
     call stresssym(gprimd,nsym,npsstr,symrec)
   end if
 end if

!DEBUG
!write(std_out,*)' forstrnps : exit '
!stop
!ENDDEBUG

 call timab(925,2,tsec)
 call timab(920,2,tsec)

end subroutine forstrnps

!!***
