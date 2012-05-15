!{\src2tex{textfont=tt}}
!!****f* ABINIT/cgwf
!! NAME
!! cgwf
!!
!! FUNCTION
!! Update all wavefunction |C>, non self-consistently.
!! also compute the corresponding H|C> and Vnl|C> (and S|C> if paw).
!! Uses a conjugate-gradient algorithm.
!! In case of paw, resolves a generalized eigenproblem using an
!!  overlap matrix (not used for norm conserving psps).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  berryopt == 4: electric field is on;
!!              5: magnetic field is on;
!!              all other values, no field is present
!!  chkexit= if non-zero, check whether the user wishes to exit
!!  cpus = CPU time limit
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  filnam_ds1=name of input file (used for exit checking)
!!  filstat=name of the status file
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array cg
!!  ikpt=number of the k-point
!!  inonsc=index of non self-consistent loop
!!  isppol=spin polarization currently treated
!!  kg_k(3,npw)=coordinates of planewaves in basis sphere.
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mband =maximum number of bands
!!  mcg=second dimension of the cg array
!!  mcgq=second dimension of the cgq array
!!  mgfft=maximum size of 1D FFTs
!!  mgsc=second dimension of the gsc array
!!  mkgq = second dimension of pwnsfacq
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw
!!  natom=number of atoms in cell.
!!  nband=number of bands.
!!  nbdblock=number of bands in a block
!!  nkpt=number of k points
!!  nline=number of line minimizations per band.
!!  nloalg(5) data concerning nonlop application
!!  npw=number of planewaves in basis sphere at given k.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=number of spin polarizations
!!  ntypat=number of types of atoms in cell.
!!  nvloc=final dimension of vlocal (usually 1, but 4 for non-collinear
!!  n4,n5,n6 used for dimensionning of vlocal
!!  ortalg=governs the choice of the algorithm for orthogonalisation.
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  pwnsfacq(2,mkgq) = phase factors for the nearest neighbours of the
!!                     current k-point (electric field, MPI //)
!!  tolwfr=tolerance on largest wf residual
!!  use_subovl=1 if the overlap matrix is not identity in WFs subspace
!!  vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!!  vxctaulocal(n4,n5,n6,nvloc,4)= local potential corresponding to the derivative of XC energy with respect to
!!    kinetic energy density, in real space, on the augmented fft grid. (optional argument)
!!    This array contains also the gradient of vxctaulocal (gvxctaulocal) in vxctaulocal(:,:,:,:,2:4).
!!  wfoptalg=govern the choice of algorithm for wf optimisation
!!   (0, 1, 10 and 11 : in the present routine, usual CG algorithm ;
!!   (2 and 3 : use shifted square Hamiltonian)
!!  zshift(nband)=in case wfoptalg is 2 or 3, shift of the Hamiltonian
!!
!! OUTPUT
!!  dphase_k(3) = change in Zak phase for the current k-point in case berryopt = 4 (electric field)
!!  resid(nband)=wf residual for new states=|(H-e)|C>|^2 (hartree^2)
!!  subham(nband*(nband+1))=Hamiltonian expressed in sthe WFs subspace
!!  subovl(nband*(nband+1)*use_subovl)=overlap matrix expressed in sthe WFs subspace
!!  subvnl(nband*(nband+1)*(1-gs_hamk%usepaw))=non-local Hamiltonian expressed in sthe WFs subspace
!!
!! SIDE EFFECTS
!!  cg(2,mcg)
!!    at input =wavefunction <G|C band,k> coefficients for ALL bands
!!    at output same as input except that
!!      the current band, with number 'band' has been updated
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  quit= if 1, proceeds to smooth ending of the job.
!!  if(gs_hamk%usepaw==1)
!!   gsc(2,mgsc)=<G|S|C band,k> coefficients for ALL bands
!!               where S is the overlap matrix (used only for paw)
!!
!! NOTES
!!  cg should not be filtered and normalized : it should already
!!   be OK at input !
!!  Not sure that that the generalized eigenproblem (when gs_hamk%usepaw=1)
!!   is compatible with wfoptalg=2 or 3 (use of shifted square
!!   Hamiltonian) - to be verified
!!
!! PARENTS
!!      vtowfk
!!
!! CHILDREN
!!      bestwfs,chkexi,cprj_alloc,cprj_copy,cprj_free,cprj_get
!!      cprj_mpi_allgather,cprj_put,dotprod_g,etheta,getcprj,getghc,leave_new
!!      linemin,make_grad_berry,mksubham,precon,projbd,smatrix,smatrix_k_paw
!!      sqnorm_g,status,sym_cprj_kn,timab,wrtout,xallgather_mpi,xcast_mpi
!!      xcomm_init,xme_init,xrecv_mpi,xsend_mpi,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine cgwf(berryopt,cg,cgq,chkexit,cpus,dimffnl,dphase_k,dtefield,&
&                ffnl,filnam_ds1,filstat,&
&                gsc,gs_hamk,icg,igsc,ikpt,inonsc,&
&                isppol,kg_k,kinpw,lmnmax,matblk,mband,&
&                mcg,mcgq,mgfft,mgsc,mkgq,mpi_enreg,mpsang,&
&                mpssoang,mpw,natom,nband,nbdblock,nkpt,nline,nloalg,npw,npwarr,&
&                nspinor,nsppol,ntypat,nvloc,n4,n5,n6,ortalg,&
&                paral_kgb,ph3d,prtvol,pwind,pwind_alloc,pwnsfac,&
&                pwnsfacq,quit,resid,subham,subovl,subvnl,tolwfr,&
&                use_subovl,vlocal,wfoptalg,zshift,&
&                vxctaulocal) ! optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_xmpi
 use m_efield

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgwf'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_53_spacepar
 use interfaces_59_io_mpi
 use interfaces_65_nonlocal
 use interfaces_66_paw
 use interfaces_66_wfs
 use interfaces_67_common, except_this_one => cgwf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: berryopt,chkexit,dimffnl,icg,igsc,ikpt,inonsc,isppol,lmnmax,matblk
 integer,intent(in) :: mband,mcg,mcgq,mgfft,mgsc,mkgq,mpsang,mpssoang,mpw,n4
 integer,intent(in) :: n5,n6,natom,nband,nbdblock,nkpt,nline,npw,nspinor,nsppol,ntypat
 integer,intent(in) :: nvloc,ortalg,paral_kgb,prtvol,pwind_alloc,use_subovl,wfoptalg
 integer,intent(inout) :: quit
 real(dp),intent(in) :: cpus,tolwfr
 character(len=fnlen),intent(in) :: filnam_ds1,filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(efield_type),intent(inout) :: dtefield
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 integer,intent(in) :: kg_k(3,npw),nloalg(5),npwarr(nkpt),pwind(pwind_alloc,2,3)
 real(dp),intent(in) :: cgq(2,mcgq),ffnl(npw,dimffnl,lmnmax,ntypat),kinpw(npw)
 real(dp),intent(in) :: pwnsfac(2,pwind_alloc),pwnsfacq(2,mkgq),zshift(nband)
 real(dp),intent(inout) :: cg(2,mcg),gsc(2,mgsc)
 real(dp),intent(inout) :: ph3d(2,npw,matblk),vlocal(n4,n5,n6,nvloc)
 real(dp),intent(out) :: dphase_k(3)
 real(dp),intent(out) :: subham(nband*(nband+1)),subovl(nband*(nband+1)*use_subovl)
 real(dp),intent(out) :: subvnl(nband*(nband+1)*(1-gs_hamk%usepaw))
 real(dp),intent(out) :: resid(nband)
 real(dp), intent(inout), optional :: vxctaulocal(n4,n5,n6,nvloc,4)

!Local variables-------------------------------
 integer,parameter :: level=113,tim_getghc=1
 integer,save :: nskip=0
 integer :: choice,counter,cpopt,ddkflag,dimenlc1,dimenlr1,dimenl2,iat,iatom,itypat
 integer :: iband,ibandmin,ibandmax
 integer :: ibdblock,iblock,icg1,icg_shift,icp1,icp2,idir,idum1,ierr,iexit,ifor,igs,igsc_shift,ii,ikgf
 integer :: ikpt2,ikpt2f,ikptf,iline,iproc,iprint,ipw,ispinor,istwf_k,isubh,isubo,itrs
!integer :: jband ! used below in debugging
 integer :: job,mcg_q,me_distrb,ncpgr,nkpg,nbdblock_eff,nblock,nproc_distrb,npw_k2
 integer :: openexit,optekin,outofplace,paw_opt,signs,shiftbd,sij_opt,spaceComm_distrb
 integer :: tag,tim_projbd,use_vnl,useoverlap,wfopta10
 real(dp) :: chc,costh,deltae,deold,dhc,dhd,diff,dotgg,dotgp,doti,dotr
 real(dp) :: dphase_aux2,e0,e0_old,e1,e1_old,eval,gamma
 real(dp) :: lam0,lamold,root,sinth,sintn,swap,tan2th,theta,thetam
 real(dp) :: xnorm
 logical :: gen_eigenpb, finite_field
 character(len=500) :: message
 integer :: hel(2,3)
 integer,allocatable :: dimlmn(:),dimlmn_srt(:),ikptf_recv(:),pwind_k(:),sflag_k(:)
 real(dp) :: bcut(2,3),dphase_aux1(3),dtm_k(2),phase_end(3)
 real(dp) :: phase_init(3),tsec(2)
 real(dp),allocatable :: cg1_k(:,:),cgq_k(:,:),conjgr(:,:),cwavef(:,:)
 real(dp),allocatable :: detovc(:,:,:),detovd(:,:,:),direc(:,:),direc_tmp(:,:),g_dummy(:,:,:),g_dummy_block(:,:,:,:)
 real(dp),allocatable :: gcc(:,:),gcc_block(:,:,:)
 real(dp),allocatable :: gh_direc(:,:),gh_direcws(:,:),ghc(:,:),ghc_all(:,:),ghc_block(:,:,:),ghcws(:,:)
 real(dp),allocatable :: grad_berry(:,:),grad_total(:,:),gs_direc(:,:)
 real(dp),allocatable :: gscc(:,:),gscc_block(:,:,:),gsc_dummy(:,:)
 real(dp),allocatable :: gvnlc(:,:),gvnlc_block(:,:,:),gvnl_direc(:,:),gvnl_dummy(:,:),kpg(:,:)
 real(dp),allocatable :: pcon(:),pwnsfac_k(:,:),scprod(:,:),scwavef(:,:)
 real(dp),allocatable :: smat_inv(:,:,:),smat_k(:,:,:),smat_k_paw(:,:,:),swork(:,:),vresid(:,:),work(:,:)
 real(dp),allocatable :: sub_tmp(:)
 type(cprj_type) :: cprj_dum(1,1)
 type(cprj_type),allocatable :: cprj_k(:,:),cprj_kb(:,:)
 type(cprj_type),allocatable :: cprj_direc(:,:),cprj_band_srt(:,:),cprj_gat(:,:)
 type(cprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Touching chkexit, cpus,filnam_ds to avoid warning for abirules. This is dirty...
 if(chkexit<0)then
   message='   chkexit should be positive !'
   MSG_BUG(message)
 end if
 if(cpus<0 .and. filnam_ds1=='a')then
   message='   cpus should be positive !'
   MSG_BUG(message)
 end if

!======================================================================
!========= LOCAL VARIABLES DEFINITIONS AND ALLOCATIONS ================
!======================================================================

!Starting the routine
 call timab(22,1,tsec)
 if(prtvol<0) then
   call status(0,filstat,iexit,level,'enter cgwf    ')
 end if

!Structured debugging if prtvol==-level
 if(prtvol==-level)then
   write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' cgwf : enter '
   call wrtout(std_out,message,'PERS')
 end if

!MPI data
 call xcomm_init(mpi_enreg,spaceComm_distrb)
 call xme_init(mpi_enreg,me_distrb)
 nproc_distrb=xcomm_size(spaceComm_distrb)

!if PAW, one has to solve a generalized eigenproblem
!(H|Psi>=Lambda.S|Psi>)
!else,   one has to solve a classical eigenproblem
!(H|Psi>=Lambda.|Psi>)
 gen_eigenpb=(gs_hamk%usepaw==1)
 useoverlap=0;if (gen_eigenpb) useoverlap=1

!if PAW, no need to compute Vnl contributions
 use_vnl=0;if (gs_hamk%usepaw==0) use_vnl=1
 if (gen_eigenpb.and.(use_vnl==1)) stop "Error in cgwf: contact Abinit group"

!Initializations and allocations
 tim_projbd=1
 isubh=1;isubo=1
 nblock=(nband-1)/nbdblock+1
 istwf_k=gs_hamk%istwf_k
 iprint=0;if(prtvol==-level)iprint=1
 wfopta10=mod(wfoptalg,10)
 optekin=0;if (wfoptalg>=10) optekin=1
 cpopt=-1
 ABI_ALLOCATE(pcon,(npw))
 ABI_ALLOCATE(ghc,(2,npw*nspinor))
 ABI_ALLOCATE(gvnlc,(2,npw*nspinor))
 ABI_ALLOCATE(conjgr,(2,npw*nspinor))
 ABI_ALLOCATE(cwavef,(2,npw*nspinor))
 ABI_ALLOCATE(direc,(2,npw*nspinor))
 ABI_ALLOCATE(scprod,(2,nband))
 if (gen_eigenpb)  then
   ABI_ALLOCATE(scwavef,(2,npw*nspinor))
   ABI_ALLOCATE(direc_tmp,(2,npw*nspinor))
 end if
 if (gen_eigenpb.and.(inonsc==1))  then
   ABI_ALLOCATE(ghc_all,(2,nband*npw*nspinor))
 end if
 if (wfopta10==2.or.wfopta10==3)  then
   ABI_ALLOCATE(work,(2,npw*nspinor))
 end if
 if (gen_eigenpb.and.(wfopta10==2.or.wfopta10==3))  then
   ABI_ALLOCATE(swork,(2,npw*nspinor))
 end if
 if(wfopta10==2 .or. wfopta10==3)then
   ABI_ALLOCATE(ghcws,(2,npw*nspinor))
   ABI_ALLOCATE(gh_direcws,(2,npw*nspinor))
   ABI_ALLOCATE(gvnl_dummy,(2,npw*nspinor))
 end if
 outofplace=0;if(wfopta10==1)outofplace=1
 if(outofplace==1) then
   ABI_ALLOCATE(gcc,(2,npw*nspinor))
   if(gen_eigenpb) ABI_ALLOCATE(gscc,(2,npw*nspinor))
 end if

!if "generalized eigenproblem", not sure of wfoptalg=2,3 algorithms
 if ((gen_eigenpb).and.(wfopta10==2.or.wfopta10==3)) then
   write(message, '(a,a,a,a,a)' ) &
&   '  Conjugate gradient algorithm not tested with',ch10,&
&   '  wfoptalg=2 or 3 and usepaw==1 !',ch10,&
&   '  Program will continue at your own risk...'
   MSG_WARNING(message)
 end if

!Define blocks of vectors
 if(wfopta10==1)then
   ABI_ALLOCATE(gcc_block,(2,npw*nspinor,nbdblock))
   ABI_ALLOCATE(ghc_block,(2,npw*nspinor,nbdblock))
   ABI_ALLOCATE(gvnlc_block,(2,npw*nspinor,nbdblock*use_vnl))
   ABI_ALLOCATE(gscc_block,(2,npw*nspinor,nbdblock*useoverlap))
 end if

!Electric field: definition of local variables:
!detovc(1:2,ifor,idir) determinant of the overlap matrix
!S_{nm}(k,k+dk)=<u_{n,k}|u_{m,k+dk}>, with the states at
!k as bras (neighbor is specified by ifor and idir)
!detovd                same as detovc but with <u_{n,k}| replaced by
!<D| (search direction) in the band-th line
!grad_berry(1:2,ipw)   Berry phase term contribution to the gradient vector
!hel(ifor,idir)        helicity of the ellipse associated w/ (ifor,idir)
!bcut(ifor,idir)       branch cut of the ellipse associated w/ (ifor,idir)
!theta_min             optimal angle theta in line_minimization when electric
!field is on
!grad_total(1:2,ipw)   total gradient (zero field term plus Berry phase term)

 finite_field = ( (berryopt == 4) .or. (berryopt == 5) )
 ncpgr = 0 ! do not think the cprj's here need gradients (no force computation in cgwf)

 if (finite_field) then

!  ji : These could be a couple of new input variables (but it is OK to define them here)
   ikptf = dtefield%i2fbz(ikpt)
   ikgf = dtefield%fkgindex(ikptf)  ! this is the shift for pwind
   mcg_q = mpw*mband*nspinor
   ABI_ALLOCATE(detovc,(2,2,3))
   ABI_ALLOCATE(detovd,(2,2,3))
   ABI_ALLOCATE(grad_berry,(2,npw*nspinor))
   ABI_ALLOCATE(cg1_k,(2,mpw))
   ABI_ALLOCATE(cgq_k,(2,mcg_q))
   ABI_ALLOCATE(grad_total,(2,npw*nspinor))
   ABI_ALLOCATE(sflag_k,(dtefield%nband_occ))
   ABI_ALLOCATE(pwind_k,(mpw))
   ABI_ALLOCATE(pwnsfac_k,(4,mpw))
   ABI_ALLOCATE(smat_k,(2,dtefield%nband_occ,dtefield%nband_occ))
   ABI_ALLOCATE(smat_inv,(2,dtefield%nband_occ,dtefield%nband_occ))
!  now the special features if using PAW
   if (gs_hamk%usepaw /= 0) then
     ABI_ALLOCATE(smat_k_paw,(2,gs_hamk%usepaw*dtefield%nband_occ,gs_hamk%usepaw*dtefield%nband_occ))
     smat_k_paw(:,:,:) = zero
!    the following are arguments to nonlop used to apply the on-site dipole to direc vector
     choice = 1 ! only apply projectors
     paw_opt = 1 ! only apply Dij
     signs = 2 ! apply nonlop to vector in k space
!    following two items are the nonlocal potential strength dij due to the on-site dipoles
     dimenlc1 = 2*lmnmax*(lmnmax+1)/2
     dimenlr1 = lmnmax*(lmnmax+1)/2
     dimenl2 = natom
!    here are some items defined by nonlop and getcprj calls but ignored in them
     nkpg = 0
     ABI_ALLOCATE(kpg,(npw,nkpg))
!    cprj structures for berryopt 4 case
     ABI_ALLOCATE(dimlmn,(natom))
     do iatom = 1, natom
       itypat = gs_hamk%typat(iatom)
       dimlmn(iatom)=dtefield%lmn_size(itypat)
     end do
     ABI_ALLOCATE(dimlmn_srt,(natom))
     iatom = 0
     do itypat = 1, ntypat
       do iat = 1, gs_hamk%nattyp(itypat)
         iatom = iatom + 1
         dimlmn_srt(iatom) = dtefield%lmn_size(itypat)
       end do
     end do
     ABI_ALLOCATE(ikptf_recv,(nproc_distrb))
     ABI_ALLOCATE(cprj_k,(natom,dtefield%nband_occ))
     ABI_ALLOCATE(cprj_kb,(natom,dtefield%nband_occ))
     ABI_ALLOCATE(cprj_direc,(natom,dtefield%nband_occ))
     ABI_ALLOCATE(cprj_band_srt,(natom,1))
     ABI_ALLOCATE(cprj_gat,(natom,dtefield%nband_occ*nproc_distrb))
     call cprj_alloc(cprj_k,ncpgr,dimlmn)
     call cprj_alloc(cprj_kb,ncpgr,dimlmn)
     call cprj_alloc(cprj_direc,ncpgr,dimlmn)
     call cprj_alloc(cprj_band_srt,ncpgr,dimlmn_srt)
     call cprj_alloc(cprj_gat,ncpgr,dimlmn)
     if (nkpt /= dtefield%fnkpt) then
       ABI_ALLOCATE(cprj_fkn,(natom,dtefield%nband_occ))
       ABI_ALLOCATE(cprj_ikn,(natom,dtefield%nband_occ))
       call cprj_alloc(cprj_fkn,ncpgr,dimlmn)
       call cprj_alloc(cprj_ikn,ncpgr,dimlmn)
     end if
   end if
 end if

!======================================================================
!If generalized eigenproblem: has to know <g|S|c> for all
!bands (for orthogonalization purpose); take benefit of this
!calculation to compute <g|H|c> at the same time, and cprj_k if
!berryopt==4
!======================================================================
 if (gen_eigenpb.and.(inonsc==1)) then
   do iblock=1,nblock
     ibandmin=1+(iblock-1)*nbdblock
     ibandmax=min(iblock*nbdblock,nband)
     do iband=ibandmin,ibandmax
       ibdblock=iband-(iblock-1)*nbdblock
       icg_shift=npw*nspinor*(iband-1)+icg
       igsc_shift=npw*nspinor*(iband-1)+igsc

       if ((wfopta10==1).and.(mpi_enreg%paralbd>=1)) then
         if ((mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me_distrb).and.&
         (mpi_enreg%me_group/=0)) cycle
       end if

       if( wfopta10/=1 .or.&
&       mpi_enreg%paralbd == 0 .or.&
&       ((wfopta10==1).and. (mpi_enreg%paralbd >= 1))) then
!        
!        the following test is already enforced with the previous cycle clause!
!        .and. (mpi_enreg%proc_distrb(ikpt,iband,isppol)==me_distrb))) then

!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(cg,cwavef,icg_shift,npw,nspinor)
         do ipw=1,npw*nspinor
           cwavef(1:2,ipw)=cg(1:2,ipw+icg_shift)
         end do
!        $OMP END PARALLEL DO
!        Compute <g|H|c>
         sij_opt=1
         if(prtvol<0) then
           call status(0,filstat,iexit,level,'call getghc   ')
         end if
         if (finite_field .and. gs_hamk%usepaw == 1) then
           if(present(vxctaulocal))then
             call getghc(0,cwavef,cprj_band_srt,dimffnl,ffnl,filstat,ghc,scwavef,gs_hamk,gvnlc,kg_k,&
&             kinpw,eval,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,1,&
&             npw,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal,vxctaulocal=vxctaulocal)
           else
             call getghc(0,cwavef,cprj_band_srt,dimffnl,ffnl,filstat,ghc,scwavef,gs_hamk,gvnlc,kg_k,&
&             kinpw,eval,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,1,&
&             npw,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
           end if
           call cprj_put(gs_hamk%atindx,cprj_band_srt,cprj_k,natom,iband,0,ikpt,1,1,&
&           dtefield%nband_occ,1,mpi_enreg,natom,1,dtefield%nband_occ,&
&           dimlmn,1,1,spaceComm_distrb,0)
         else
           if(present(vxctaulocal))then
             call getghc(cpopt,cwavef,cprj_dum,dimffnl,ffnl,filstat,ghc,scwavef,gs_hamk,gvnlc,kg_k,&
&             kinpw,eval,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,1,&
&             npw,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal,vxctaulocal=vxctaulocal)
           else
             call getghc(cpopt,cwavef,cprj_dum,dimffnl,ffnl,filstat,ghc,scwavef,gs_hamk,gvnlc,kg_k,&
&             kinpw,eval,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,1,&
&             npw,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
           end if
         end if
         if(prtvol<0) then
           call status(0,filstat,iexit,level,'after getghc  ')
         end if
       end if

       if(wfopta10==1)then
         if (mpi_enreg%paralbd>=1) then
           tag=nbdblock*(mpi_enreg%proc_distrb(ikpt,iband,isppol))+ibdblock
           if (mpi_enreg%me_group==0) then
             if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me_distrb) then
               ABI_ALLOCATE(g_dummy,(2,npw*nspinor,2))
               call timab(48,1,tsec)
               call xrecv_mpi(g_dummy,xmpi_any_source,tag,mpi_enreg%kpt_comm(mpi_enreg%num_group),ierr)
               call timab(48,2,tsec)
               ghc=g_dummy(:,:,1)
               scwavef=g_dummy(:,:,2)
               ABI_DEALLOCATE(g_dummy)
             end if
             do ipw=1,npw*nspinor
               ghc_block(1:2,ipw,ibdblock)=ghc(1:2,ipw)
               gscc_block(1:2,ipw,ibdblock)=scwavef(1:2,ipw)
             end do
           else
             ABI_ALLOCATE(g_dummy,(2,npw*nspinor,2))
             g_dummy(:,:,1)=ghc
             g_dummy(:,:,2)=scwavef
             call timab(48,1,tsec)
             call xsend_mpi(g_dummy,0,tag,mpi_enreg%kpt_comm(mpi_enreg%num_group),ierr)
             call timab(48,2,tsec)
             ABI_DEALLOCATE(g_dummy)
           end if
         else
           do ipw=1,npw*nspinor
             ghc_block(1:2,ipw,ibdblock)=ghc(1:2,ipw)
             gscc_block(1:2,ipw,ibdblock)=scwavef(1:2,ipw)
           end do
         end if
       else
!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(gsc,scwavef,icg,icg_shift,igsc_shift,npw,nspinor)
         do ipw=1,npw*nspinor
           ghc_all(1:2,ipw+icg_shift-icg)=ghc(1:2,ipw)
           gsc(1:2,ipw+igsc_shift)=scwavef(1:2,ipw)
         end do
!        $OMP END PARALLEL DO
       end if
     end do
     if(wfopta10==1)then
       if (mpi_enreg%paralbd >= 1) then
         ABI_ALLOCATE(g_dummy_block,(2,npw*nspinor,nbdblock,2))
         g_dummy_block(:,:,:,1)=ghc_block;g_dummy_block(:,:,:,2)=gscc_block
         call timab(48,1,tsec)
         call xcast_mpi(g_dummy_block,0,mpi_enreg%kpt_comm(mpi_enreg%num_group),ierr)
         call timab(48,2,tsec)
         ghc_block=g_dummy_block(:,:,:,1);gscc_block=g_dummy_block(:,:,:,2)
         ABI_DEALLOCATE(g_dummy_block)
       end if
       do iband=1+(iblock-1)*nbdblock,min(iblock*nbdblock,nband)
         ibdblock=iband-(iblock-1)*nbdblock
         do ipw=1,npw*nspinor
           ghc_all(1:2,ipw+icg_shift-icg)=ghc_block(1:2,ipw,ibdblock)
           gsc(1:2,ipw+igsc_shift)=gscc_block(1:2,ipw,ibdblock)
         end do
       end do
     end if
   end do
 end if

!======================================================================
!====================== LOOP OVER BANDS ===============================
!======================================================================

!Loop over blocks of bands.
!In the standard band-sequential algorithm, nblock=nband.
 do iblock=1,nblock
   counter=100*iblock*nbdblock+inonsc

!  Because in this loop, the CPU time matters, the writing
!  in the STATUS file is usually inhibited
   if(prtvol>=10) then
     call status(counter,filstat,iexit,level,'loop iband    ')
   end if

!  Not too often, check whether the run must be stopped.
!  If so, iexit will be non-zero.
!  Note that when the number of bands becomes large, the check
!  must be done more often, because treating one band takes also longer ...
   if (mpi_enreg%paral_compil==0) then
     if(iblock==1 .or. (nband>=16 .and. mod(iblock,8)==1) &
&     .or. (nband>=32 .and. mod(iblock,4)==1) &
&     .or. (nband>=64 .and. mod(iblock,2)==1) &
&     .or. (nband>=128)                        )then
       openexit=1 ; if(chkexit<=1)openexit=0
       call chkexi(cpus,filnam_ds1,iexit,6,mpi_enreg,openexit)
       if(iexit/=0)quit=1
     end if
   end if

!  Loop over bands in a block
!  This loop can be MPI-parallelized, over processors attached to the same k point
   ibandmin=1+(iblock-1)*nbdblock
   ibandmax=min(iblock*nbdblock,nband)

!  =====================================================================================
!  Big iband loop

   do iband=ibandmin,ibandmax
     ibdblock=iband-(iblock-1)*nbdblock
     counter=100*iband+inonsc
     icg_shift=npw*nspinor*(iband-1)+icg
     igsc_shift=npw*nspinor*(iband-1)+igsc

!    In MPI-parallelisation, one should determine if the present
!    band is treated by the present processor ...
     if(wfopta10==1 .and. mpi_enreg%paralbd >= 1) then
       if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/= me_distrb) then
         resid(iband)=zero
         if (mpi_enreg%me_group/=0) cycle
       end if
     end if

     if( (wfopta10/=1) .or.&
&     ((wfopta10==1) .and. (mpi_enreg%paralbd == 0)) .or. &
&     ((wfopta10==1) .and. (mpi_enreg%paralbd >= 1)) ) then
!      again, the following clause is necessarily true following the previous nested if()
!      .and. (mpi_enreg%proc_distrb(ikpt,iband,isppol)== me_distrb))) then


!      ======================================================================
!      ========== INITIALISATION OF MINIMIZATION ITERATIONS =================
!      ======================================================================


!      Tell us what is going on:
       if(prtvol>=10)then
         write(message, '(a,i6,2x,a,i3,a)' ) &
&         ' --- cgwf is called for band',iband,'for',nline,' lines'
         call  wrtout(std_out,message,'PERS')
       end if

       dotgp=one
       if (berryopt == 4) then
         detovc(:,:,:) = zero ; detovd(:,:,:) = zero
         phase_init(:) = zero
         dphase_aux1(:) = zero
         phase_end(:) = zero
         bcut(:,:) = zero
         hel(:,:) = zero
       end if


!      Extraction of the vector that is iteratively updated
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(cg,cwavef,icg_shift,npw,nspinor)
       do ipw=1,npw*nspinor
         cwavef(1,ipw)=cg(1,ipw+icg_shift)
         cwavef(2,ipw)=cg(2,ipw+icg_shift)
       end do
!      $OMP END PARALLEL DO


!      If generalized eigenproblem: extraction of the overlap information
       if (gen_eigenpb) then
!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(gsc,scwavef,igsc_shift,npw,nspinor)
         do ipw=1,npw*nspinor
           scwavef(1,ipw)=gsc(1,ipw+igsc_shift)
           scwavef(2,ipw)=gsc(2,ipw+igsc_shift)
         end do
!        $OMP END PARALLEL DO
       end if


!      Normalize incoming wf (and S.wf, if generalized eigenproblem):
!      WARNING : It might be interesting to skip the following operation.
!      The associated routines should be reexamined to see whether cwavef
!      is not already normalized.
       if (gen_eigenpb) then
         call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw*nspinor,2,cwavef,scwavef)
         dotr=sqrt(dotr**2+doti**2);xnorm=1._dp/sqrt(dotr)
!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(scwavef,cwavef,npw,nspinor,xnorm)
         do ipw=1,npw*nspinor
           cwavef(1,ipw)=cwavef(1,ipw)*xnorm
           cwavef(2,ipw)=cwavef(2,ipw)*xnorm
           scwavef(1,ipw)=scwavef(1,ipw)*xnorm
           scwavef(2,ipw)=scwavef(2,ipw)*xnorm
         end do
!        $OMP END PARALLEL DO
       else
         call sqnorm_g(dotr,istwf_k,mpi_enreg,npw*nspinor,cwavef)
         xnorm=1._dp/sqrt(abs(dotr))
!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(scwavef,cwavef,npw,nspinor,xnorm)
         do ipw=1,npw*nspinor
           cwavef(1,ipw)=cwavef(1,ipw)*xnorm
           cwavef(2,ipw)=cwavef(2,ipw)*xnorm
         end do
!        $OMP END PARALLEL DO
       end if
       if (prtvol==-level)then
         write(message,'(a,f14.6)') 'cgwf : xnorm=',xnorm
         call wrtout(std_out,message,'PERS')
       end if

!      Compute (or extract) <g|H|c>
       if (gen_eigenpb.and.(inonsc==1)) then
!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(ghc,ghc_all,icg,icg_shift,npw,nspinor,xnorm)
         do ipw=1,npw*nspinor
           ghc(1,ipw)=xnorm*ghc_all(1,ipw+icg_shift-icg)
           ghc(2,ipw)=xnorm*ghc_all(2,ipw+icg_shift-icg)
         end do
!        $OMP END PARALLEL DO
       else
         sij_opt=0
         if(prtvol<0)then
           call status(0,filstat,iexit,level,'call getghc   ')
         end if
         if(present(vxctaulocal))then
           call getghc(cpopt,cwavef,cprj_dum,dimffnl,ffnl,filstat,ghc,gsc_dummy,gs_hamk,gvnlc,kg_k,&
&           kinpw,eval,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,1,&
&           npw,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal,vxctaulocal=vxctaulocal)
         else
           call getghc(cpopt,cwavef,cprj_dum,dimffnl,ffnl,filstat,ghc,gsc_dummy,gs_hamk,gvnlc,kg_k,&
&           kinpw,eval,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,1,&
&           npw,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
         end if
         if(prtvol<0)then
           call status(0,filstat,iexit,level,'after getghc  ')
         end if
       end if


!      In case of outofplace=1, must save the (nearly) original vector
       if(outofplace==1)then
!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(cwavef,gcc,npw,nspinor)
         do ipw=1,npw*nspinor
!          gcc is used as a temporary space for the old state
           gcc(1,ipw)=cwavef(1,ipw)
           gcc(2,ipw)=cwavef(2,ipw)
         end do
!        $OMP END PARALLEL DO
         if (gen_eigenpb.and.(inonsc>1)) then
!          $OMP PARALLEL DO PRIVATE(ipw) &
!          $OMP&SHARED(scwavef,gscc,npw,nspinor)
           do ipw=1,npw*nspinor
!            gscc is used as a temporary space for the old S|c>
             gscc(1,ipw)=scwavef(1,ipw)
             gscc(2,ipw)=scwavef(2,ipw)
           end do
!          $OMP END PARALLEL DO
         end if
       end if


!      Minimisation of the residual: compute <G|(H-zshift)^2|C iband,k>
       if(wfopta10==2 .or. wfopta10==3) then
         ghcws(:,:)=ghc(:,:)
         if (gen_eigenpb) then
           sij_opt=1
           work(:,:)=ghc(:,:)-zshift(iband)*scwavef(:,:)
         else
           sij_opt=0
           work(:,:)=ghc(:,:)-zshift(iband)*cwavef(:,:)
         end if
         if(present(vxctaulocal))then
           call getghc(cpopt,work,cprj_dum,dimffnl,ffnl,filstat,ghc,swork,gs_hamk,gvnl_dummy,kg_k,&
&           kinpw,eval,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,1,npw,nspinor,ntypat,&
&           nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal,vxctaulocal=vxctaulocal)
         else
           call getghc(cpopt,work,cprj_dum,dimffnl,ffnl,filstat,ghc,swork,gs_hamk,gvnl_dummy,kg_k,&
&           kinpw,eval,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,1,npw,nspinor,ntypat,&
&           nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
         end if
         if (gen_eigenpb) then
           ghc(:,:)=ghc(:,:)-zshift(iband)*swork(:,:)
         else
           ghc(:,:)=ghc(:,:)-zshift(iband)*work(:,:)
         end if
       end if

!      ======================================================================
!      ====== BEGIN LOOP FOR A GIVEN BAND: MINIMIZATION ITERATIONS ==========
!      ======================================================================

       if(nline/=0)then
         do iline=1,nline
           if(prtvol<0) then
             call status(iline,filstat,iexit,level,'iline         ')
           end if

!          ======================================================================
!          ================= COMPUTE THE RESIDUAL ===============================
!          ======================================================================

!          Compute lambda = <C|H|C> or <C|(H-zshift)**2|C>
           call dotprod_g(chc,doti,istwf_k,mpi_enreg,npw*nspinor,1,cwavef,ghc)
           lam0=chc

!          Check that lam0 is decreasing on succeeding lines:
           if (berryopt /= 4) then
             if (iline==1) then
               lamold=lam0
             else
               if (lam0>lamold+1.d-12) then
                 write(message, '(a,i8,a,1p,e14.6,a1,3x,a,1p,e14.6,a1)')&
&                 '  New trial energy at line',iline,' = ',lam0,ch10,&
&                 '  is higher than former:',lamold,ch10
                 MSG_WARNING(message)
               end if
               lamold=lam0
             end if
           end if

!          Compute residual vector:
!          Note that vresid is precomputed to garantee cancellation of errors
!          and allow residuals to reach values as small as 1.0d-24 or better.
           ABI_ALLOCATE(vresid,(2,npw*nspinor))
           if (wfopta10<=1) then
             eval=chc
             if (gen_eigenpb) then
!              $OMP PARALLEL DO PRIVATE(ipw) &
!              $OMP&SHARED(scwavef,chc,ghc,npw,nspinor,vresid)
               do ipw=1,npw*nspinor
                 vresid(1,ipw)=ghc(1,ipw)-chc*scwavef(1,ipw)
                 vresid(2,ipw)=ghc(2,ipw)-chc*scwavef(2,ipw)
               end do
!              $OMP END PARALLEL DO
             else
!              $OMP PARALLEL DO PRIVATE(ipw) &
!              $OMP&SHARED(cwavef,chc,ghc,npw,nspinor,vresid)
               do ipw=1,npw*nspinor
                 vresid(1,ipw)=ghc(1,ipw)-chc*cwavef(1,ipw)
                 vresid(2,ipw)=ghc(2,ipw)-chc*cwavef(2,ipw)
               end do
!              $OMP END PARALLEL DO
             end if
           else
             call dotprod_g(eval,doti,istwf_k,mpi_enreg,npw*nspinor,1,cwavef,ghcws)
             if (gen_eigenpb) then
!              $OMP PARALLEL DO PRIVATE(ipw) &
!              $OMP&SHARED(scwavef,eval,ghcws,npw,nspinor,vresid)
               do ipw=1,npw*nspinor
                 vresid(1,ipw)=ghcws(1,ipw)-eval*scwavef(1,ipw)
                 vresid(2,ipw)=ghcws(2,ipw)-eval*scwavef(2,ipw)
               end do
!              $OMP END PARALLEL DO
             else
!              $OMP PARALLEL DO PRIVATE(ipw) &
!              $OMP&SHARED(cwavef,eval,ghcws,npw,nspinor,vresid)
               do ipw=1,npw*nspinor
                 vresid(1,ipw)=ghcws(1,ipw)-eval*cwavef(1,ipw)
                 vresid(2,ipw)=ghcws(2,ipw)-eval*cwavef(2,ipw)
               end do
!              $OMP END PARALLEL DO
             end if
           end if

!          Compute residual (squared) norm
           call sqnorm_g(resid(iband),istwf_k,mpi_enreg,npw*nspinor,vresid)
           if(prtvol==-level)then
             write(message,'(a,a,i3,2f14.6)') ch10,&
&             'cgwf : iline,eval,resid=',iline,eval,resid(iband)
             call wrtout(std_out,message,'PERS')
           end if

!          ======================================================================
!          ============== CHECK FOR CONVERGENCE CRITERIA ========================
!          ======================================================================

!          If residual sufficiently small stop line minimizations
           if (resid(iband)<tolwfr) then
             if(prtvol>=10)then
               write(message, '(a,i4,a,i2,a,es12.4)' ) &
&               ' cgwf: band',iband,' converged after ',iline,&
&               ' line minimizations : resid =',resid(iband)
               call wrtout(std_out,message,'PERS')
             end if
!            Number of two-way 3D ffts skipped
             nskip=nskip+(nline-iline+1)
!            Exit from the loop on iline
             ABI_DEALLOCATE(vresid)
             exit
           end if

!          If user require exiting the job, stop line minimisations
           if (quit==1) then
             write(message, '(a,i4)' ) &
&             ' cgwf : user require exiting => skip update of band ',iband
             call wrtout(std_out,message,'PERS')
!            Number of two-way 3D ffts skipped
             nskip=nskip+(nline-iline+1)
!            Exit from the loop on iline
             ABI_DEALLOCATE(vresid)
             exit
           end if


!          ======================================================================
!          =========== COMPUTE THE STEEPEST DESCENT DIRECTION ===================
!          ======================================================================

!          Compute the steepest descent direction
           if (gen_eigenpb) then
!            Store <G|H-lambda.S|C> in direc
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(direc,npw,nspinor,vresid)
             direc(1:2,1:npw*nspinor)=vresid(1:2,1:npw*nspinor)
!            do ipw=1,npw*nspinor
!            direc(1,ipw)=vresid(1,ipw)
!            direc(2,ipw)=vresid(2,ipw)
!            end do
!            $OMP END PARALLEL DO
           else
!            Store <G|H|C> in direc
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(direc,ghc,npw,nspinor)
             direc(1:2,1:npw*nspinor)=ghc(1:2,1:npw*nspinor)
!            do ipw=1,npw*nspinor
!            direc(1,ipw)=ghc(1,ipw)
!            direc(2,ipw)=ghc(2,ipw)
!            end do
!            $OMP END PARALLEL DO
           end if
           ABI_DEALLOCATE(vresid)
!          
!          Electric field: compute the gradient of the Berry phase part of the energy functional.
!          See PRL 89, 117602 (2002), grad_berry(:,:) is the second term of Eq. (4)
           if (berryopt == 4) then

             call make_grad_berry(cg,cgq,cprj_k,detovc,dimffnl,dimlmn,dimlmn_srt,direc,&
&             dtefield,ffnl,grad_berry,gs_hamk,iband,icg,ikpt,isppol,kg_k,lmnmax,matblk,&
&             mband,mcg,mcgq,mgfft,mkgq,mpi_enreg,mpsang,mpssoang,mpw,natom,nkpt,nloalg,&
&             npw,npwarr,nspinor,nsppol,ntypat,pwind,ph3d,pwind_alloc,pwnsfac,pwnsfacq)

!            Add grad_berry to direc and store original gradient
             direc(:,:) = direc(:,:) + grad_berry(:,:)
             grad_total(:,:) = direc(:,:)
!            DEBUG: check that grad_berry is orthogonal to the occupied manifold at k
!            do jband = 1, dtefield%nband_occ
!            dotr = zero  ;  doti = zero
!            do ipw = 1, npw*nspinor
!            if(.not.gen_eigenpb) then
!            dotr = dotr + cg(1,icg + (jband-1)*npw*nspinor + ipw)*grad_berry(1,ipw) + &
!            &                 cg(2,icg + (jband-1)*npw*nspinor + ipw)*grad_berry(2,ipw)
!            doti = doti + cg(1,icg + (jband-1)*npw*nspinor + ipw)*grad_berry(2,ipw) - &
!            &                 cg(2,icg + (jband-1)*npw*nspinor + ipw)*grad_berry(1,ipw)
!            end if
!            end do
!            if ((abs(dotr) > tol12).or.(abs(doti) > tol12)) then
!            write(std_out,'(a)')'cgwf-berry : ERROR (orthogonality)'
!            write(std_out,'(3(2x,i3),2(5x,e16.9))')ikpt,iband,jband,dotr,doti
!            stop
!            end if
!            end do
!            ENDDEBUG
           end if   ! berryopt == 4

!          ======================================================================
!          =========== PROJECT THE STEEPEST DESCENT DIRECTION ===================
!          ========= OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS ================
!          ======================================================================

!          The following projection over the subspace orthogonal to occupied bands
!          is optional. It is a bit more accurate, but doubles the number of N^3 ops.
!          It is done only if ortalg>=0.

!          Project the steepest descent direction:
!          direc(2,npw)=<G|H|Cnk> - \sum_{(i<=n)} <G|H|Cik> , normalized.

           if(ortalg>=0)then
             if(prtvol<0) then
               call status(iline,filstat,iexit,level,'projbd(1)     ')
             end if
             if (gen_eigenpb) then
               call projbd(cg,direc,iband,icg,igsc,istwf_k,mcg,mpi_enreg,mgsc,nband,npw,nspinor,&
&               ortalg,iprint,gsc,scprod,0,tim_projbd,useoverlap)
             else
               call projbd(cg,direc,-1   ,icg,igsc,istwf_k,mcg,mpi_enreg,mgsc,nband,npw,nspinor,&
&               ortalg,iprint,gsc,scprod,0,tim_projbd,useoverlap)
             end if
           else
!            For negative ortalg must still project current band out
!            of conjugate vector (unneeded if gen_eigenpb)
             if (.not.gen_eigenpb) then
               call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw*nspinor,3,cwavef,direc)
               if(istwf_k==1)then
!                $OMP PARALLEL DO PRIVATE(ipw) &
!                $OMP&SHARED(cwavef,direc,doti,dotr,npw,nspinor)
                 do ipw=1,npw*nspinor
                   direc(1,ipw)=direc(1,ipw)-(dotr*cwavef(1,ipw)-doti*cwavef(2,ipw))
                   direc(2,ipw)=direc(2,ipw)-(dotr*cwavef(2,ipw)+doti*cwavef(1,ipw))
                 end do
!                $OMP END PARALLEL DO
               else
!                $OMP PARALLEL DO PRIVATE(ipw) &
!                $OMP&SHARED(cwavef,direc,dotr,npw,nspinor)
                 do ipw=1,npw*nspinor
                   direc(1,ipw)=direc(1,ipw)-dotr*cwavef(1,ipw)
                   direc(2,ipw)=direc(2,ipw)-dotr*cwavef(2,ipw)
                 end do
!                $OMP END PARALLEL DO
               end if
             end if
           end if

!          For a generalized eigenpb, store the steepest descent direction
           if (gen_eigenpb) direc_tmp=direc

!          ======================================================================
!          ======== PRECONDITION THE STEEPEST DESCENT DIRECTION =================
!          ======================================================================

!          If wfoptalg>=10, the precondition matrix is kept constant
!          during iteration ; otherwise it is recomputed

           if (wfoptalg<10.or.iline==1) then

             if(prtvol<0) then
               call status(iline,filstat,iexit,level,'call precon   ')
             end if
             call precon(cwavef,zero,istwf_k,kinpw,mpi_enreg,npw,nspinor,optekin,pcon,direc)

!            Minimisation of the residual: must precondition twice
!            (might make only one call, with modified precon routine - might also make a shift !!!)
             if(wfopta10==2 .or. wfopta10==3)then
               call precon(cwavef,zero,istwf_k,kinpw,mpi_enreg,npw,nspinor,optekin,pcon,direc)
               if(iline==1)then
!                $OMP PARALLEL DO PRIVATE(ipw) &
!                $OMP&SHARED(npw,pcon)
                 do ipw=1,npw
                   pcon(ipw)=pcon(ipw)**2
                   pcon(ipw)=pcon(ipw)**2
                 end do
!                $OMP END PARALLEL DO
               end if
             end if
           else

             do ispinor=1,nspinor
               igs=(ispinor-1)*npw
!              $OMP PARALLEL DO PRIVATE(ipw) &
!              $OMP&SHARED(igs,npw,direc,pcon)
               do ipw=1+igs,npw+igs
                 direc(1,ipw)=direc(1,ipw)*pcon(ipw-igs)
                 direc(2,ipw)=direc(2,ipw)*pcon(ipw-igs)
               end do
!              $OMP END PARALLEL DO
             end do
           end if

!          ======================================================================
!          ======= PROJECT THE PRECOND. STEEPEST DESCENT DIRECTION ==============
!          ========= OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS ================
!          ======================================================================
!          Projecting again out all bands (not normalized).
           if(prtvol<0) then
             call status(0,filstat,iexit,level,'prjbd(2)      ')
           end if
           call projbd(cg,direc,-1,icg,igsc,istwf_k,mcg,mpi_enreg,mgsc,nband,npw,nspinor,&
&           ortalg,iprint,gsc,scprod,0,tim_projbd,useoverlap)

!          ======================================================================
!          ================= COMPUTE THE CONJUGATE-GRADIENT =====================
!          ======================================================================

           if (berryopt == 4) then
             call dotprod_g(dotgg,doti,istwf_k,mpi_enreg,npw*nspinor,1,direc,grad_total)
!            DEBUG (electric field)
!            check that the dotproduct is real
!            if (abs(doti) > tol8) then
!            write(std_out,*) ' cgwf-berry: ERROR'
!            write(std_out,*) ' doti = ',doti
!            stop
!            end if
!            ENDDEBUG
           else
             if (gen_eigenpb) then
               call dotprod_g(dotgg,doti,istwf_k,mpi_enreg,npw*nspinor,1,direc,direc_tmp)
             else
               call dotprod_g(dotgg,doti,istwf_k,mpi_enreg,npw*nspinor,1,direc,ghc)
             end if
           end if

!          MJV: added 5 Feb 2012 - causes divide by 0 on next iteration of iline
           if (abs(dotgg) < TINY(0.0_dp)*1.e50_dp) dotgg = TINY(0.0_dp)*1.e50_dp

!          At first iteration, gamma is set to zero
           if (iline==1) then
             gamma=zero
             dotgp=dotgg
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(conjgr,direc,npw,nspinor)
             do ipw=1,npw*nspinor
               conjgr(1,ipw)=direc(1,ipw)
               conjgr(2,ipw)=direc(2,ipw)
             end do
!            $OMP END PARALLEL DO

           else

             gamma=dotgg/dotgp
             dotgp=dotgg
             if(prtvol==-level)then
               write(message,'(a,2es16.6)') 'cgwf : dotgg,gamma =',dotgg,gamma
               call wrtout(std_out,message,'PERS')
             end if

!            Note: another way to compute gamma: Polak, Ribiere
!            no real improvment ; to be more carrefully tested
!            call dotprod_g(dotgg,doti,istwf_k,mpi_enreg,npw*nspinor,1,direc,direc_tmp)
!            !direcp must be set to zero at the beginning
!            direcp=direc-direcp
!            call dotprod_g(dotgmg,doti,istwf_k,mpi_enreg,npw*nspinor,1,direcp,direc_tmp)
!            direcp=direc;gamma=dotgmg/dotgp;dotgp=dotgmg

!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(conjgr,direc,gamma,npw,nspinor)
             do ipw=1,npw*nspinor
               conjgr(1,ipw)=direc(1,ipw)+gamma*conjgr(1,ipw)
               conjgr(2,ipw)=direc(2,ipw)+gamma*conjgr(2,ipw)
             end do
!            $OMP END PARALLEL DO
           end if

!          ======================================================================
!          ============ PROJECTION OF THE CONJUGATED GRADIENT ===================
!          ======================================================================

           if (gen_eigenpb) then
             call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw*nspinor,3,scwavef,conjgr)
           else
             call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw*nspinor,3,cwavef,conjgr)
           end if

!          Project the conjugated gradient onto the current band
           if(istwf_k==1)then
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(conjgr,cwavef,direc,doti,dotr,npw,nspinor)
             do ipw=1,npw*nspinor
               direc(1,ipw)=conjgr(1,ipw)-(dotr*cwavef(1,ipw)-doti*cwavef(2,ipw))
               direc(2,ipw)=conjgr(2,ipw)-(dotr*cwavef(2,ipw)+doti*cwavef(1,ipw))
             end do
!            $OMP END PARALLEL DO
           else
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(conjgr,cwavef,direc,dotr,npw,nspinor)
             do ipw=1,npw*nspinor
               direc(1,ipw)=conjgr(1,ipw)-dotr*cwavef(1,ipw)
               direc(2,ipw)=conjgr(2,ipw)-dotr*cwavef(2,ipw)
             end do
!            $OMP END PARALLEL DO
           end if

!          In case of generalized eigenproblem, normalization of direction vector
!          cannot be done here (because S|D> is not known here).
           if (.not.gen_eigenpb) then
             call sqnorm_g(dotr,istwf_k,mpi_enreg,npw*nspinor,direc)
             xnorm=1._dp/sqrt(abs(dotr))
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(direc,npw,nspinor,xnorm)
             do ipw=1,npw*nspinor
               direc(1,ipw)=direc(1,ipw)*xnorm
               direc(2,ipw)=direc(2,ipw)*xnorm
             end do
!            $OMP END PARALLEL DO
             xnorm=one
           end if

!          ======================================================================
!          ===== COMPUTE CONTRIBUTIONS TO 1ST AND 2ND DERIVATIVES OF ENERGY =====
!          ======================================================================

!          Compute gh_direc = <G|H|D> and eventually gs_direc = <G|S|D>
           ABI_ALLOCATE(gh_direc,(2,npw*nspinor))
           ABI_ALLOCATE(gvnl_direc,(2,npw*nspinor))
           if (gen_eigenpb)  then
             ABI_ALLOCATE(gs_direc,(2,npw*nspinor))
           end if
           if(prtvol<0) then
             call status(iline,filstat,iexit,level,'call getghc   ')
           end if
           sij_opt=0;if (gen_eigenpb) sij_opt=1
           if(present(vxctaulocal))then
             call getghc(cpopt,direc,cprj_dum,dimffnl,ffnl,filstat,gh_direc,gs_direc,gs_hamk,gvnl_direc,&
&             kg_k,kinpw,eval,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,&
&             1,npw,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal,vxctaulocal=vxctaulocal)
           else
             call getghc(cpopt,direc,cprj_dum,dimffnl,ffnl,filstat,gh_direc,gs_direc,gs_hamk,gvnl_direc,&
&             kg_k,kinpw,eval,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,&
&             1,npw,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
           end if
           if(prtvol<0) then
             call status(iline,filstat,iexit,level,'after getghc  ')
           end if

           if(wfopta10==2 .or. wfopta10==3)then
!            Minimisation of the residual, so compute <G|(H-zshift)^2|D>
             gh_direcws(:,:)=gh_direc(:,:)
             if (gen_eigenpb) then
               sij_opt=1
               work(:,:)=gh_direc(:,:)-zshift(iband)*gs_direc(:,:)
             else
               sij_opt=0
               work(:,:)=gh_direc(:,:)-zshift(iband)*direc(:,:)
             end if
             if(prtvol<0) then
               call status(iline,filstat,iexit,level,'call getghc   ')
             end if
             if(present(vxctaulocal))then
               call getghc(cpopt,work,cprj_dum,dimffnl,ffnl,filstat,gh_direc,swork,gs_hamk,gvnl_dummy,&
&               kg_k,kinpw,eval,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,&
&               1,npw,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,0,tim_getghc,0,vlocal,vxctaulocal=vxctaulocal)
             else
               call getghc(cpopt,work,cprj_dum,dimffnl,ffnl,filstat,gh_direc,swork,gs_hamk,gvnl_dummy,&
&               kg_k,kinpw,eval,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,&
&               1,npw,nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d,prtvol,0,tim_getghc,0,vlocal)
             end if
             if(prtvol<0) then
               call status(iline,filstat,iexit,level,'after getghc  ')
             end if
             if (gen_eigenpb) then
               gh_direc(:,:)=gh_direc(:,:)-zshift(iband)*swork(:,:)
             else
               gh_direc(:,:)=gh_direc(:,:)-zshift(iband)*work(:,:)
             end if
           end if

!          In case of generalized eigenproblem, compute now the norm of the conjugated gradient
           if (gen_eigenpb) then
             call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw*nspinor,1,direc,gs_direc)
             xnorm=1._dp/sqrt(abs(dotr))
           end if

!          Compute dhc = Re{<D|H|C>}
           call dotprod_g(dhc,doti,istwf_k,mpi_enreg,npw*nspinor,1,direc,ghc)
           dhc=dhc*xnorm

!          Compute <D|H|D> or <D|(H-zshift)^2|D>
           call dotprod_g(dhd,doti,istwf_k,mpi_enreg,npw*nspinor,1,direc,gh_direc)
           dhd=dhd*xnorm**2

           if(prtvol==-level)then
             write(message,'(a,3f14.6)') 'cgwf : chc,dhc,dhd=',chc,dhc,dhd
             call wrtout(std_out,message,'PERS')
           end if

!          ======================================================================
!          ======= COMPUTE MIXING FACTORS - CHECK FOR CONVERGENCE ===============
!          ======================================================================

           if (berryopt /=4) then

!            Compute tan(2 theta),sin(theta) and cos(theta)
             tan2th=2.0_dp*dhc/(chc-dhd)

             if (abs(tan2th)<1.d-05) then

               costh=1.0_dp-0.125_dp*tan2th**2
               sinth=0.5_dp*tan2th*(1.0_dp-0.375_dp*tan2th**2)

!              Check that result is above machine precision
               if (abs(sinth)<epsilon(0._dp)) then
                 write(message, '(a,es16.4)' ) ' cgwf: converged with tan2th=',tan2th
                 call wrtout(std_out,message,'PERS')
!                Number of one-way 3D ffts skipped
                 nskip=nskip+2*(nline-iline)
!                Exit from the loop on iline
                 ABI_DEALLOCATE(gh_direc)
                 ABI_DEALLOCATE(gvnl_direc)
                 if (gen_eigenpb) ABI_DEALLOCATE(gs_direc)
                 exit
               end if

             else
               root=sqrt(1.0_dp+tan2th**2)
               costh=sqrt(0.5_dp+0.5_dp/root)
               sinth=sign(sqrt(0.5_dp-0.5_dp/root),tan2th)
             end if

!            Check for lower of two possible roots (same sign as curvature at theta where slope is zero)
             diff=(chc-dhd)
!            Swap c and d if value of diff is positive
             if (diff>zero) then
               swap=costh
               costh=-sinth
               sinth=swap
               if(prtvol<0 .or. prtvol>=10)then
                 write(message,*)'   Note: swap roots, iline,diff=',iline,diff
                 call wrtout(std_out,message,'PERS')
               end if
             end if

           else

!            In case the eletric field is on, the line minimization has to be done numerically

!            Compute determinant of the overlap matrix where in the band-th line
!            the wavefunction is replaced by the search direction
             job = 10 ; shiftbd = 0
             do idir = 1, 3
!              do not do this for efield_dot(idir)=0
               if (abs(dtefield%efield_dot(idir)) < tol12) cycle
               do ifor = 1, 2
                 ikpt2f = dtefield%ikpt_dk(ikptf,ifor,idir)
                 if (dtefield%indkk_f2ibz(ikpt2f,6) == 1) then
                   itrs = 10
                 else
                   itrs = 0
                 end if
                 ikpt2 = dtefield%indkk_f2ibz(ikpt2f,1)
                 npw_k2 = npwarr(ikpt2)
                 pwind_k(1:npw) = pwind(ikgf+1:ikgf+npw,ifor,idir)
                 pwnsfac_k(1:2,1:npw) = pwnsfac(1:2,ikgf+1:ikgf+npw)
                 sflag_k(:) = dtefield%sflag(:,ikpt+(isppol-1)*nkpt,ifor,idir)
                 smat_k(:,:,:) = dtefield%smat(:,:,:,ikpt+(isppol-1)*nkpt,ifor,idir)
                 if (mpi_enreg%paral_compil_kpt == 1) then
                   icg1 = dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt)
                   cgq_k(:,1:dtefield%nband_occ*nspinor*npw_k2) = &
&                   cgq(:,icg1+1:icg1+dtefield%nband_occ*nspinor*npw_k2)
                   idum1 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt)
                   pwnsfac_k(3:4,1:npw_k2) = pwnsfacq(1:2,idum1+1:idum1+npw_k2)
                 else
                   icg1 = dtefield%cgindex(ikpt2,isppol)
                   cgq_k(:,1:dtefield%nband_occ*nspinor*npw_k2) = &
&                   cg(:,icg1+1:icg1+dtefield%nband_occ*nspinor*npw_k2)
                   idum1=dtefield%fkgindex(ikpt2f)
                   pwnsfac_k(3:4,1:npw_k2) = pwnsfac(1:2,idum1+1:idum1+npw_k2)
                 end if
                 icg1 = 0 ; ddkflag = 0
                 if (gen_eigenpb) then
!                  $OMP PARALLEL DO PRIVATE(ipw) &
!                  $OMP&SHARED(direc,npw,nspinor,xnorm)
                   do ipw=1,npw*nspinor
                     direc_tmp(1,ipw)=direc(1,ipw)*xnorm
                     direc_tmp(2,ipw)=direc(2,ipw)*xnorm
                   end do
!                  $OMP END PARALLEL DO
!                  
!                  need cprj corresponding to direc_tmp in order to make smat_k_paw properly
!                  
                   call getcprj(1,0,direc_tmp,cprj_band_srt,gs_hamk%dimekb1,gs_hamk%dimekb2,&
&                   dimffnl,gs_hamk%ekb,ffnl,0,gs_hamk%indlmn,istwf_k,kg_k,kpg,gs_hamk%kpoint,&
&                   lmnmax,matblk,mgfft,mpi_enreg,natom,gs_hamk%nattyp,gs_hamk%ngfft,0,&
&                   nloalg,npw,nspinor,ntypat,gs_hamk%phkxred,gs_hamk%ph1d,ph3d,gs_hamk%ucvol,&
&                   gs_hamk%usepaw,gs_hamk%useylm)
!                  icp1=dtefield%nband_occ*(ikptf-1)
                   icp2=dtefield%nband_occ*(ikpt2-1)
                   call cprj_copy(cprj_k,cprj_direc)
                   call cprj_put(gs_hamk%atindx,cprj_band_srt,cprj_direc,natom,iband,0,ikpt,1,1,&
&                   dtefield%nband_occ,1,mpi_enreg,natom,1,dtefield%nband_occ,&
&                   dimlmn,1,1,spaceComm_distrb,0)
                   call cprj_get(gs_hamk%atindx1,cprj_kb,dtefield%cprj,natom,1,icp2,ikpt,0,isppol,&
&                   dtefield%nband_occ,dtefield%fnkpt,mpi_enreg,natom,&
&                   dtefield%nband_occ,dtefield%nband_occ,nspinor,nsppol,0)
                   if (ikpt2 /= ikpt2f) then ! construct cprj_kb by symmetry
                     call cprj_copy(cprj_kb,cprj_ikn)
                     call sym_cprj_kn(cprj_fkn,cprj_ikn,dtefield%atom_indsym,dimlmn,-1,gs_hamk%indlmn,&
&                     dtefield%indkk_f2ibz(ikpt2f,2),dtefield%indkk_f2ibz(ikpt2f,6),&
&                     dtefield%fkptns(:,dtefield%i2fbz(ikpt2)),&
&                     dtefield%lmax,dtefield%lmnmax,natom,dtefield%nband_occ,nspinor,&
&                     dtefield%nsym,ntypat,gs_hamk%typat,dtefield%zarot)
                     call cprj_copy(cprj_fkn,cprj_kb)
                   end if
                   call smatrix_k_paw(cprj_direc,cprj_kb,dtefield,idir,ifor,natom,smat_k_paw,gs_hamk%typat)
                   call smatrix(direc_tmp,cgq_k,cg1_k,ddkflag,dtm_k,icg1,icg1,&
&                   itrs,job,iband,npw,mcg_q,mpw,iband,&
&                   mpw,dtefield%nband_occ,&
&                   npw,npw_k2,nspinor,pwind_k,pwnsfac_k,sflag_k,&
&                   shiftbd,smat_inv,smat_k,smat_k_paw,gs_hamk%usepaw)
                 else
                   call smatrix(direc,cgq_k,cg1_k,ddkflag,dtm_k,icg1,icg1,&
&                   itrs,job,iband,npw,mcg_q,mpw,iband,&
&                   mpw,dtefield%nband_occ,&
&                   npw,npw_k2,nspinor,pwind_k,pwnsfac_k,sflag_k,&
&                   shiftbd,smat_inv,smat_k,smat_k_paw,gs_hamk%usepaw)
                 end if
                 detovd(:,ifor,idir) = dtm_k(:) ! Store the determinant of the overlap
!                matrix (required to compute theta_min)
!                DEBUG
!                write(std_out,*)'cgwf-berry: detovc and detovd'
!                write(std_out,*)detovc(:,ifor,idir)
!                write(std_out,*)detovd(:,ifor,idir)
!                write(std_out,*)'smat_k'
!                do jband = 1, 4
!                write(std_out,'(4(2x,e14.6))')smat_k(1,jband,:)
!                write(std_out,'(4(2x,e14.6))')smat_k(2,jband,:)
!                write(std_out,*)
!                end do
!                ENDDEBUG
               end do  ! ifor
             end do    ! idir
             call linemin(bcut,chc,costh,detovc,detovd,dhc,dhd,&
&             dphase_aux1,dtefield%efield_dot,iline,&
&             dtefield%fnkpt,dtefield%nstr,hel,phase_end,&
&             phase_init,dtefield%sdeg,sinth,thetam)
!            DEBUG
!            if (mpi_enreg%me == 1) then
!            write(std_out,*)'after linemin '
!            write(std_out,'(a,3(2x,f16.9))')'phase_init  = ',phase_init(:)
!            write(std_out,'(a,3(2x,f16.9))')'phase_end   = ',phase_end(:)
!            write(std_out,'(a,3(2x,f16.9))')'dphase_aux1 = ',dphase_aux1(:)
!            write(std_out,*) 'thetam',thetam
!            end if
!            ENDDEBUG
           end if  ! berryopt/=4

!          ======================================================================
!          =========== GENERATE NEW |wf>, H|wf>, Vnl|Wf>, S|Wf> ... =============
!          ======================================================================

           sintn=sinth*xnorm

!          $OMP PARALLEL DO PRIVATE(ipw) &
!          $OMP&SHARED(cg,costh,cwavef,direc,ghc,gh_direc,icg_shift,npw,nspinor,sintn)
           do ipw=1,npw*nspinor
             cwavef(1,ipw)=cwavef(1,ipw)*costh+direc(1,ipw)*sintn
             cwavef(2,ipw)=cwavef(2,ipw)*costh+direc(2,ipw)*sintn
             cg(1,ipw+icg_shift)=cwavef(1,ipw)
             cg(2,ipw+icg_shift)=cwavef(2,ipw)
             ghc(1,ipw)  =ghc(1,ipw)  *costh + gh_direc(1,ipw)*sintn
             ghc(2,ipw)  =ghc(2,ipw)  *costh + gh_direc(2,ipw)*sintn
           end do
!          $OMP END PARALLEL DO
           if (use_vnl==1) then
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(costh,gvnlc,gvnl_direc,icg,npw,nspinor,sintn)
             do ipw=1,npw*nspinor
               gvnlc(1,ipw)=gvnlc(1,ipw)*costh + gvnl_direc(1,ipw)*sintn
               gvnlc(2,ipw)=gvnlc(2,ipw)*costh + gvnl_direc(2,ipw)*sintn
             end do
!            $OMP END PARALLEL DO
           end if
           if (gen_eigenpb) then
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(costh,scwavef,gsc,gs_direc,igsc_shift,npw,nspinor,sintn)
             do ipw=1,npw*nspinor
               scwavef(1,ipw)=scwavef(1,ipw)*costh+gs_direc(1,ipw)*sintn
               scwavef(2,ipw)=scwavef(2,ipw)*costh+gs_direc(2,ipw)*sintn
               gsc(1,ipw+igsc_shift)=scwavef(1,ipw)
               gsc(2,ipw+igsc_shift)=scwavef(2,ipw)
             end do
             if (berryopt == 4) then  ! must update cprj for the new wavefunction
               call getcprj(1,0,cwavef,cprj_band_srt,gs_hamk%dimekb1,gs_hamk%dimekb2,&
&               dimffnl,gs_hamk%ekb,ffnl,0,gs_hamk%indlmn,istwf_k,kg_k,kpg,gs_hamk%kpoint,&
&               lmnmax,matblk,mgfft,mpi_enreg,natom,gs_hamk%nattyp,gs_hamk%ngfft,0,&
&               nloalg,npw,nspinor,ntypat,gs_hamk%phkxred,gs_hamk%ph1d,ph3d,gs_hamk%ucvol,&
&               gs_hamk%usepaw,gs_hamk%useylm)
               call cprj_put(gs_hamk%atindx,cprj_band_srt,cprj_k,natom,iband,0,ikpt,1,1,&
&               dtefield%nband_occ,1,mpi_enreg,natom,1,dtefield%nband_occ,&
&               dimlmn,1,1,spaceComm_distrb,0)
             end if
!            $OMP END PARALLEL DO
           end if
           if(wfopta10==2 .or. wfopta10==3)then
!            Need to keep track of ghcws, in order to avoid recomputing it
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(costh,ghcws,gh_direcws,icg,npw,nspinor,sintn)
             do ipw=1,npw*nspinor
               ghcws(1,ipw)=ghcws(1,ipw)*costh + gh_direcws(1,ipw)*sintn
               ghcws(2,ipw)=ghcws(2,ipw)*costh + gh_direcws(2,ipw)*sintn
             end do
!            $OMP END PARALLEL DO
           end if

           ABI_DEALLOCATE(gh_direc)
           ABI_DEALLOCATE(gvnl_direc)
           if (gen_eigenpb) ABI_DEALLOCATE(gs_direc)

!          ======================================================================
!          =========== CHECK CONVERGENCE AGAINST TRIAL ENERGY ===================
!          ======================================================================

!          Compute delta(E)
           if (berryopt /= 4) then

             deltae=chc*(costh**2-1._dp)+dhd*sinth**2+2._dp*costh*sinth*dhc

           else

!            Compute deltae
             call etheta(bcut,chc,detovc,detovd,dhc,dhd,dtefield%efield_dot,e0,e1,&
&             hel,dtefield%fnkpt,dtefield%nstr,dtefield%sdeg,thetam)
             theta = zero
             call etheta(bcut,chc,detovc,detovd,dhc,dhd,&
&             dtefield%efield_dot,e0_old,e1_old,&
&             hel,dtefield%fnkpt,dtefield%nstr,dtefield%sdeg,theta)
             deltae = e0 - e0_old
!            DEBUG
!            write(std_out,*) 'e0, e0_old, deltae', e0, e0_old, deltae
!            ENDDEBUG
!            Check that e0 is decreasing on succeeding lines:
!            if (deltae > zero) then
             if (deltae > tol12) then ! exploring different checks for berryopt = 4
               write(message, '(3a,i8,a,1p,e14.6,a1,3x,a,1p,e14.6,a1)')&
&               '  (electric field)',ch10,&
&               '  New trial energy at line',iline,' = ',e0,ch10,&
&               '  is higher than former:',e0_old,ch10
               MSG_WARNING(message)
             end if

           end if         ! berryopt /= 4

!          Check convergence and eventually exit
           if (iline==1) then
             deold=deltae
           else if (abs(deltae)<0.005_dp*abs(deold) .and. iline/=nline .and. wfopta10<2)then
             if(prtvol>=10)then
               write(message, '(a,i4,1x,a,1p,e12.4,a,e12.4,a)' ) &
&               ' cgwf: line',iline,&
&               ' deltae=',deltae,' < 0.005*',deold,' =>skip lines'
               call wrtout(std_out,message,'PERS')
             end if
!            Number of one-way 3D ffts skipped
             nskip=nskip+2*(nline-iline)
!            Exit from the loop on iline
             exit
           end if

!          ======================================================================
!          ================ END LOOP FOR A GIVEN BAND ===========================
!          ======================================================================

!          Note that there are three "exit" instructions inside
         end do

!        Additionnal computations in case of electric field
         if (berryopt == 4) then
!          Bring present contribution to dphasek(idir) into [-pi,pi]
           do idir = 1, 3
             dphase_aux2 = mod(phase_end(idir) - phase_init(idir) + 100*two_pi,two_pi)
             if (dphase_aux2 > pi) dphase_aux2 = dphase_aux2 - two_pi
!            DEBUG
!            dphase_aux1(idir)=mod(dphase_aux1(idir)+100*two_pi,two_pi)
!            if(dphase_aux1(idir)>pi) dphase_aux1(idir)=dphase_aux1(idir)-two_pi
!            diff = dphase_aux2 - dphase_aux1(idir)
!            if (abs(diff) > tol10) then
!            write(std_out,*)'cgwf-berry: ERROR'
!            write(std_out,'(a,3(2x,i3),f16.9)')'ikpt,iband,idir,diff',ikpt,iband,idir,diff
!            stop
!            end if
!            write(100,*) idir, dphase_aux2
!            ENDDEBUG
             dphase_k(idir) = dphase_k(idir) + dphase_aux2
!            DEBUG
!            write(std_out,*) 'idir,phase_init,phase_end,dphase_k'
!            write(std_out,*) idir,phase_init(idir),phase_end(idir),dphase_k(idir)
!            ENDDEBUG
           end do
         end if   ! berryopt == 4

       else

!        nline==0 , needs to provide a residual
!        ===================================================
         resid(iband)=-one
!        End nline==0 case
       end if

!      ======================================================================
!      =============== END OF CURRENT BAND: CLEANING ========================
!      ======================================================================

!      It was checked that getghc is NOT needed here : equivalent results with the copy below.
       if(wfopta10==2 .or. wfopta10==3) ghc(:,:)=ghcws(:,:)

       if (berryopt == 4) dtefield%sflag(:,ikpt + (isppol-1)*nkpt,:,:) = 0

!      Special treatment for outofplace==1 :
       if(outofplace==1)then
!        Put (nearly) original cg back in place, and free gcc for
!        the new vector
!        $OMP PARALLEL DO PRIVATE(ipw) &
!        $OMP&SHARED(cg,cwavef,gcc,icg_shift,npw,nspinor)
         do ipw=1,npw*nspinor
           cg(1,ipw+icg_shift)=gcc(1,ipw)
           cg(2,ipw+icg_shift)=gcc(2,ipw)
           gcc(1,ipw)=cwavef(1,ipw)
           gcc(2,ipw)=cwavef(2,ipw)
         end do
!        $OMP END PARALLEL DO
!        Eventually do the same with gsc and gcc
         if (useoverlap==1) then
           if (inonsc>1) then
!            $OMP PARALLEL DO PRIVATE(ipw) &
!            $OMP&SHARED(gsc,gscc,igsc_shift,npw,nspinor)
             do ipw=1,npw*nspinor
               gsc(1,ipw+igsc_shift)=gscc(1,ipw)
               gsc(2,ipw+igsc_shift)=gscc(2,ipw)
             end do
!            $OMP END PARALLEL DO
           end if
!          $OMP PARALLEL DO PRIVATE(ipw) &
!          $OMP&SHARED(scwavef,gscc,icg,npw,nspinor)
           do ipw=1,npw*nspinor
             gscc(1,ipw)=scwavef(1,ipw)
             gscc(2,ipw)=scwavef(2,ipw)
           end do
!          $OMP END PARALLEL DO
         end if
       end if ! outofplace==1

!      At the end of the treatment of a set of bands, write the number
!      of one-way 3D ffts skipped
       if (mpi_enreg%paral_compil==0) then
         if(paral_kgb == 0) then
           if(iband==nband .and. prtvol/=0)then
             write(message, '(a,i8)' )&
&             ' cgwf : number of one-way 3D ffts skipped in cgwf until now =',nskip
             call wrtout(std_out,message,'PERS')
           end if
         end if
       end if

       if(prtvol<0) then
         call status(0,filstat,iexit,level,'after iline   ')
       end if

!      Structured debugging : if prtvol=-level, stop here.
       if(prtvol==-level)then
         write(message,'(a1,a,a1,a,i1,a)') ch10,&
&         ' cgwf : exit ',&
&         ch10,'  prtvol=-',level,', debugging mode => stop '
         call wrtout(std_out,message,'PERS')
         call leave_new('PERS')
       end if

!      ======================================================================
!      ============== END OF CURRENT BLOCK: CLEANING ========================
!      ======================================================================

     end if

!    If blocked algorithm, save (or transmit) the output of cgwf
!    Note: in MPI-parallel, one will have to send
!    gcc,ghc,gscc and gvnlc to the master processor of this k point.
!    resid is NOT transmitted to all procs
     if(wfopta10==1)then
       if (mpi_enreg%paralbd >= 1) then
!        020827 : iband could not be larger than 100
!        tag=100*(mpi_enreg%proc_distrb(ikpt,iband,isppol)+1) + iband
         tag=nbdblock*(mpi_enreg%proc_distrb(ikpt,iband,isppol)) + ibdblock
         if (mpi_enreg%me_group==0) then
           if (mpi_enreg%proc_distrb(ikpt,iband,isppol)== me_distrb) then
!            I am the master and I have the data
           else
!            I have to receive the data
             ABI_ALLOCATE(g_dummy,(2,npw*nspinor,(2+useoverlap+use_vnl)))
             call timab(48,1,tsec)
             call xrecv_mpi(g_dummy,xmpi_any_source,tag,mpi_enreg%kpt_comm(mpi_enreg%num_group),ierr)
             call timab(48,2,tsec)
!            One should directly copy to the *_block arrays
             gcc=g_dummy(:,:,1)
             ghc=g_dummy(:,:,2)
             if(useoverlap==1)gscc=g_dummy(:,:,3)
             if(use_vnl==1)gvnlc=g_dummy(:,:,3+useoverlap)
             ABI_DEALLOCATE(g_dummy)
           end if
           do ipw=1,npw*nspinor
             gcc_block(1,ipw,ibdblock)=gcc(1,ipw)
             gcc_block(2,ipw,ibdblock)=gcc(2,ipw)
             ghc_block(1,ipw,ibdblock)=ghc(1,ipw)
             ghc_block(2,ipw,ibdblock)=ghc(2,ipw)
           end do
           if (useoverlap==1) then
             do ipw=1,npw*nspinor
               gscc_block(1,ipw,ibdblock)=gscc(1,ipw)
               gscc_block(2,ipw,ibdblock)=gscc(2,ipw)
             end do
           end if
           if (use_vnl==1) then
             do ipw=1,npw*nspinor
               gvnlc_block(1,ipw,ibdblock)=gvnlc(1,ipw)
               gvnlc_block(2,ipw,ibdblock)=gvnlc(2,ipw)
             end do
           end if
         else
!          I am not the master and I have the data to send to the master
           ABI_ALLOCATE(g_dummy,(2,npw*nspinor,(2+useoverlap+use_vnl)))
           g_dummy(:,:,1)=gcc
           g_dummy(:,:,2)=ghc
           if(useoverlap==1)g_dummy(:,:,3)=gscc
           if(use_vnl==1)g_dummy(:,:,3+useoverlap)=gvnlc
           call timab(48,1,tsec)
           call xsend_mpi(g_dummy,0,tag,mpi_enreg%kpt_comm(mpi_enreg%num_group),ierr)
           call timab(48,2,tsec)
           ABI_DEALLOCATE(g_dummy)
         end if
       else !(mpi_enreg%paralbd == 0)
         do ipw=1,npw*nspinor
           gcc_block(1,ipw,ibdblock)=gcc(1,ipw)
           gcc_block(2,ipw,ibdblock)=gcc(2,ipw)
           ghc_block(1,ipw,ibdblock)=ghc(1,ipw)
           ghc_block(2,ipw,ibdblock)=ghc(2,ipw)
         end do
         if (useoverlap==1) then
           do ipw=1,npw*nspinor
             gscc_block(1,ipw,ibdblock)=gscc(1,ipw)
             gscc_block(2,ipw,ibdblock)=gscc(2,ipw)
           end do
         end if
         if (use_vnl==1) then
           do ipw=1,npw*nspinor
             gvnlc_block(1,ipw,ibdblock)=gvnlc(1,ipw)
             gvnlc_block(2,ipw,ibdblock)=gvnlc(2,ipw)
           end do
         end if
       end if

     end if !(end if wfoptalg==1)

   end do ! iband in a block

!  End big iband loop
!  ===========================================================================


!  If blocked algorithm, determine now the "best" wavefunctions
!  among the subspace. If MPI-parallel, only the master should do it,
!  then transmit gcc_block to all processors of this k point,
!  for storage in their cg
   if(wfopta10==1)then
     if ((mpi_enreg%paralbd == 0) .or. ((mpi_enreg%paralbd >= 1) .and. (mpi_enreg%me_group==0))) then
       nbdblock_eff=min(nbdblock,nband-(iblock-1)*nbdblock)
       call bestwfs(gcc_block,ghc_block,gscc_block,useoverlap,gvnlc_block,use_vnl,&
&       istwf_k,mpi_enreg,nbdblock,npw,nspinor,nbdblock_eff,nbdblock_eff,wfoptalg)
     end if
!    Proc master bcast data to others proc of the group
     if (mpi_enreg%paralbd >= 1) then
       ABI_ALLOCATE(g_dummy_block,(2,npw*nspinor,nbdblock,(2+useoverlap+use_vnl)))
       g_dummy_block(:,:,:,1)=gcc_block
       g_dummy_block(:,:,:,2)=ghc_block
       if(useoverlap==1)g_dummy_block(:,:,:,3)=gscc_block
       if(use_vnl==1)g_dummy_block(:,:,:,3+useoverlap)=gvnlc_block
       call timab(48,1,tsec)
       call xcast_mpi(g_dummy_block,0,mpi_enreg%kpt_comm(mpi_enreg%num_group),ierr)
       call timab(48,2,tsec)
       gcc_block=g_dummy_block(:,:,:,1)
       ghc_block=g_dummy_block(:,:,:,2)
       if(useoverlap==1)gscc_block=g_dummy_block(:,:,:,3)
       if(use_vnl==1)gvnlc_block=g_dummy_block(:,:,:,3+useoverlap)
       ABI_DEALLOCATE(g_dummy_block)
     end if

!    If MPI-parallel, all "blocked" arrays (gcc_block,gscc_block,ghc_block,gvnlc_block)
!    (for ibdblock=1 to nbdblock) have been transmitted
!    to all processors of this k point.
!    All processors must have all data for all bands of the block
     do iband=1+(iblock-1)*nbdblock,min(iblock*nbdblock,nband)
       ibdblock=iband-(iblock-1)*nbdblock
       icg_shift=npw*nspinor*(iband-1)+icg
       do ipw=1,npw*nspinor
         cg(1,ipw+icg_shift)=gcc_block(1,ipw,ibdblock)
         cg(2,ipw+icg_shift)=gcc_block(2,ipw,ibdblock)
       end do
       if (useoverlap==1) then
         igsc_shift=npw*nspinor*(iband-1)+igsc
         do ipw=1,npw*nspinor
           gsc(1,ipw+igsc_shift)=gscc_block(1,ipw,ibdblock)
           gsc(2,ipw+igsc_shift)=gscc_block(2,ipw,ibdblock)
         end do
       end if
     end do
   end if !(end if wfoptalg==1)

!  ======================================================================
!  ============= COMPUTE HAMILTONIAN IN WFs SUBSPACE ====================
!  ======================================================================

   call mksubham(cg,ghc,ghc_block,gsc,gvnlc,gvnlc_block,iblock,icg,igsc,ikpt,isppol,istwf_k,&
&   isubh,isubo,mcg,mgsc,mpi_enreg,nband,nbdblock,npw,&
&   nspinor,subham,subovl,subvnl,use_subovl,use_vnl,wfoptalg)

!  DEBUG
!  write(std_out,*)' cgwf : iblock,subham(1)=',iblock,subham(1)
!  ENDDDEBUG

   if(prtvol>=10)then
     call status(counter,filstat,iexit,level,'end loop iband')
   end if

!  End loop over block of bands
 end do ! iblock

 if (berryopt == 4 .and. gs_hamk%usepaw == 1) then ! store updated cprjs for this kpt
!  switch from ikptf to ikpt
   ikptf = ikpt
   call xallgather_mpi(ikptf,ikptf_recv,spaceComm_distrb,ierr)
   call cprj_mpi_allgather(cprj_k,cprj_gat,natom,dtefield%nband_occ,dimlmn,ncpgr,nproc_distrb,spaceComm_distrb,ierr)
   do iproc = 1, nproc_distrb
     icp2=dtefield%nband_occ*(iproc-1)
     call cprj_get(gs_hamk%atindx1,cprj_k,cprj_gat,natom,1,icp2,ikpt,0,1,&
&     dtefield%nband_occ,nproc_distrb,mpi_enreg,natom,&
&     dtefield%nband_occ,dtefield%nband_occ,1,1,0)
!    ikptf = ikptf_recv(iproc)
     icp1 = dtefield%nband_occ*(ikptf_recv(iproc)-1)
     call cprj_put(gs_hamk%atindx1,cprj_k,dtefield%cprj,natom,1,icp1,ikpt,0,isppol,&
&     dtefield%nband_occ,nkpt,mpi_enreg,natom,dtefield%nband_occ,&
&     dtefield%nband_occ,dimlmn,nspinor,nsppol,spaceComm_distrb,0)
   end do
 end if

!Parallelism
!If MPI-parallel, subham, subovl and subvnl will have to be transmitted to all processors of current k point
 if (mpi_enreg%paral_compil==1) then

!  Pack the needed data in sub_tmp
   ABI_ALLOCATE(sub_tmp,((1+use_vnl+use_subovl)*nband*(nband+1)))
   sub_tmp(1:nband*(nband+1))=subham
   if (use_vnl==1) sub_tmp(1+nband*(nband+1):2*nband*(nband+1))=subvnl
   if (use_subovl==1) sub_tmp(1+(1+use_vnl)*nband*(nband+1):(2+use_vnl)*nband*(nband+1))=subovl

   if ((wfopta10==1) .and. (mpi_enreg%paralbd >= 1)) then
     call timab(48,1,tsec)
     call xsum_mpi(sub_tmp,mpi_enreg%kpt_comm(mpi_enreg%num_group),ierr)
     call timab(48,2,tsec)
   end if

!  Unpack the data
   subham=sub_tmp(1:nband*(nband+1))
   if (use_vnl==1) subvnl=sub_tmp(1+nband*(nband+1):2*nband*(nband+1))
   if (use_subovl==1) subovl=sub_tmp(1+(1+use_vnl)*nband*(nband+1):(2+use_vnl)*nband*(nband+1))
   ABI_DEALLOCATE(sub_tmp)
 end if

!Debugging ouputs
 if(prtvol==-level)then
   isubh=1
   if (use_vnl==1) write(message,'(a)') ' cgwf : isubh  subham(isubh:isubh+1)  subvnl(isubh:isubh+1)'
   if (use_vnl==0) write(message,'(a)') ' cgwf : isubh  subham(isubh:isubh+1)'
   do iband=1,nband
     do ii=1,iband
       if (use_vnl==1) then
         write(message,'(i5,4es16.6)')isubh,subham(isubh:isubh+1),subvnl(isubh:isubh+1)
       else
         write(message,'(i5,2es16.6)')isubh,subham(isubh:isubh+1)
       end if
       call wrtout(std_out,message,'PERS')
       isubh=isubh+2
     end do
   end do
 end if

!======================================================================
!==================== FINAL DEALLOCATIONS =============================
!======================================================================

 ABI_DEALLOCATE(conjgr)
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(direc)
 ABI_DEALLOCATE(pcon)
 ABI_DEALLOCATE(scprod)
 ABI_DEALLOCATE(ghc)
 ABI_DEALLOCATE(gvnlc)
 if (gen_eigenpb)  then
   ABI_DEALLOCATE(scwavef)
   ABI_DEALLOCATE(direc_tmp)
 end if
 if (gen_eigenpb.and.(inonsc==1))  then
   ABI_DEALLOCATE(ghc_all)
 end if
 if (outofplace==1)  then
   ABI_DEALLOCATE(gcc)
 end if
 if (outofplace==1.and.gen_eigenpb)  then
   ABI_DEALLOCATE(gscc)
 end if
 if(wfopta10==1) then
   ABI_DEALLOCATE(gcc_block)
   ABI_DEALLOCATE(ghc_block)
   ABI_DEALLOCATE(gscc_block)
   ABI_DEALLOCATE(gvnlc_block)
 end if
 if(wfopta10==2.or.wfopta10==3)  then
   ABI_DEALLOCATE(ghcws)
   ABI_DEALLOCATE(gh_direcws)
   ABI_DEALLOCATE(gvnl_dummy)
 end if
 if(wfopta10==2.or.wfopta10==3)  then
   ABI_DEALLOCATE(work)
 end if
 if(gen_eigenpb.and.(wfopta10==2.or.wfopta10==3))  then
   ABI_DEALLOCATE(swork)
 end if
 if (finite_field) then
   ABI_DEALLOCATE(cg1_k)
   ABI_DEALLOCATE(cgq_k)
   ABI_DEALLOCATE(detovc)
   ABI_DEALLOCATE(detovd)
   ABI_DEALLOCATE(grad_berry)
   ABI_DEALLOCATE(sflag_k)
   ABI_DEALLOCATE(smat_inv)
   ABI_DEALLOCATE(smat_k)
   ABI_DEALLOCATE(pwind_k)
   ABI_DEALLOCATE(pwnsfac_k)
   ABI_DEALLOCATE(grad_total)
   if (gs_hamk%usepaw /= 0) then
     ABI_DEALLOCATE(smat_k_paw)
     ABI_DEALLOCATE(kpg)
     call cprj_free(cprj_k)
     call cprj_free(cprj_kb)
     call cprj_free(cprj_direc)
     call cprj_free(cprj_band_srt)
     call cprj_free(cprj_gat)
     if (nkpt /= dtefield%fnkpt) then
       call cprj_free(cprj_fkn)
       call cprj_free(cprj_ikn)
       ABI_DEALLOCATE(cprj_fkn)
       ABI_DEALLOCATE(cprj_ikn)
     end if
     ABI_DEALLOCATE(dimlmn)
     ABI_DEALLOCATE(cprj_k)
     ABI_DEALLOCATE(cprj_kb)
     ABI_DEALLOCATE(cprj_direc)
     ABI_DEALLOCATE(cprj_gat)
     ABI_DEALLOCATE(ikptf_recv)
     ABI_DEALLOCATE(cprj_band_srt)
   end if
 end if

 if(prtvol<0) then
   call status(0,filstat,iexit,level,'exit          ')
 end if
 call timab(22,2,tsec)

 DBG_EXIT("COLL")

end subroutine cgwf
!!***
