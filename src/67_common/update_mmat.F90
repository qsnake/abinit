!{\src2tex{textfont=tt}}
!!****f* ABINIT/update_mmat
!! NAME
!! update_mmat
!!
!! FUNCTION
!! This routine updates the M matrix for a Berrys phase magnetization calculation.
!!
!! COPYRIGHT
!!
!! INPUTS
!!  berryopt 5 is finite B field, -5 is just magnetization
!!  cg(2,mcg) wavefunction <G|C band,k> coefficients for ALL bands
!!  cgq(2,mcgq) wavefunctions at neighbouring k points
!!  dimffnl = 2nd dimension of ffnl (1 + number of derivatives)
!!  ffnl(npw_k,dimffnl,lmnmax,ntypat) = nonlocal form factors
!!  filstat=name of the status file
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  ikpt=number of the k-point
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kinpw(npw_k)=(modified) kinetic energy for each plane wave (Hartree)
!!  lmnmax = max number of (l,m,n) components over all atom types
!!  matblk = dimension of array ph3d
!!  mband =maximum number of bands
!!  mcg=second dimension of the cg array
!!  mcgq=second dimension of the cgq array
!!  mgfft = max size of 1D FFTs
!!  mkgq = second dimension of pwnsfacq
!!  mpi_enreg=information about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+max(spin*angular momentum) for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw
!!  natom=number of atoms in cell.
!!  nkpg = 3*optforces*dtset%nloalg(5)
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  npw_k = number of plane waves for this k point
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  ntypat=number of types of atoms in cell.
!!  nvloc=final dimension of vlocal (usually 1, but 4 for non-collinear
!!  n4,n5,n6 used for dimensionning of vlocal
!!  paral_kgb = flag defining parallelism in getghc
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph3d(2,npw_k,matblk) = 3D structure factors for each atom and planewave
!!  prtvol = flag controlling verbosity of output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  pwnsfacq(2,mkgq) = phase factors for the nearest neighbours of the
!!                     current k-point
!!  rmet(3,3)=real space metric (bohr**2)
!!  ucvol=unit cell volume in bohr**3.
!!  vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtefield <type(efield_type)> = variables related to Berry phase. The contribution to the bare
!!   magnetization at the given k point is stored in dtefield%mag_bare_k(idir,ikpt)
!!
!! NOTES
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      cprj_alloc,cprj_copy,cprj_free,cprj_get,getghc,mkffnl,mkkpg,smatrix
!!      smatrix_k_paw,sym_cprj_kn
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine update_mmat(berryopt,cg,cgq,dimffnl,dtefield,ffnl,filstat,&
&                   gmet,gprimd,gs_hamk,icg,ikpt,kg,kinpw,lmnmax,matblk,&
&                   mband,mcg,mcgq,mgfft,mkgq,mkmem,mpi_enreg,mpsang,mpssoang,mpw,&
&                   natom,nkpg,nkpt,npw_k,npwarr,nspinor,ntypat,&
&                   nvloc,n4,n5,n6,pawtab,psps,pwind,pwind_alloc,&
&                   paral_kgb,ph3d,prtvol,pwnsfac,pwnsfacq,rmet,ucvol,vlocal,ylm,ylmgr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_efield

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'update_mmat'
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_65_nonlocal
 use interfaces_66_paw
 use interfaces_66_wfs
 use interfaces_67_common, except_this_one => update_mmat
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: berryopt,icg,ikpt,dimffnl,lmnmax,matblk,mband,mcg,mcgq,mgfft
 integer,intent(in) :: mkgq,mkmem,mpsang,mpssoang,mpw,natom
 integer,intent(in) :: nkpg,nkpt,npw_k,nspinor,ntypat,nvloc,n4,n5,n6
 integer,intent(in) :: paral_kgb,prtvol,pwind_alloc
 real(dp),intent(in) :: ucvol
 character(len=fnlen),intent(in) :: filstat
 type(efield_type),intent(inout) :: dtefield
 type(gs_hamiltonian_type), intent(in) :: gs_hamk
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(in) :: psps

!arrays
 integer, intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt),pwind(pwind_alloc,2,3)
 real(dp),intent(in) :: cg(2,mcg),cgq(2,mcgq),gmet(3,3),gprimd(3,3),kinpw(npw_k)
 real(dp),intent(in) :: ffnl(npw_k,dimffnl,lmnmax,ntypat)
 real(dp),intent(inout) :: ph3d(2,npw_k,matblk)
 real(dp),intent(in) :: pwnsfac(2,pwind_alloc),pwnsfacq(2,mkgq),rmet(3,3)
 real(dp), intent(inout) :: vlocal(n4,n5,n6,nvloc) ! this variable is inout in getghc
 real(dp), intent(in) :: ylm(mpw*mkmem,mpsang*mpsang)
 real(dp), intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang)
 type(pawtab_type),intent(in)  :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: bdx,bsig,cpopt,dij_ind,ddkflag,iat,iatom,icg1,icp,icp2,idir,idum1
 integer :: ifor,ihcg,ikg,ikg2,ikgf,ikptf,il,ilm,ipw,jpw
 integer :: bband,bdir,bfor,kdir,kfor, blmn, klmn
 integer :: ikpt2,ikpt2f,istep,itrs,itypat,job,kband
 integer :: kdx,ksig,mcg1_k,mcg_q,nband_k
 integer :: ndat,nn,nnp,nnpp,nnppp,npw_k2,shiftbd,sij_opt,tim_getghc,tind,type_calc
 real(dp) :: lambda
 real(dp) :: sfac
 complex(dpc) :: IA, IA1, IA2, IA3, IB, IB1, IB2, IB3, IB4
 complex(dpc) :: IIA, IIA1, IIA2, IIA3
 complex(dpc) :: IIIA, IIIA1, IIIA2, IIIA3, IIIA4
 complex(dpc) :: bffnlfac
 complex(dpc) :: cpb, cpk, dij, onfac, onsite
 
!arrays
 integer,allocatable :: dimlmn(:),dimlmn_srt(:),kg_k(:,:),kg_k2(:,:),pwind_k(:),sflag_k(:)
 real(dp) :: dk(3),dtm_k(2),dotri(2),kpoint2(3)
 real(dp),allocatable :: bwave(:,:),kwave(:,:)
 real(dp),allocatable :: bffnl(:,:,:,:,:),bffnl_(:,:,:,:),cg1_k(:,:),cgq_k(:,:)
 real(dp),allocatable :: omat(:,:,:,:,:,:,:),ghc(:,:),gvnlc(:,:),gsc(:,:)
 real(dp),allocatable :: kpg_k2(:,:),tcg(:,:,:,:),pwnsfac_k(:,:)
 real(dp),allocatable :: smat_inv(:,:,:),smat_inv_all(:,:,:,:,:)
 real(dp),allocatable :: smat_k(:,:,:),smat_k_paw(:,:,:)
 real(dp),allocatable :: ylm_k2(:,:)
! the following is (-i)^L mod 4.
 complex(dpc),dimension(0:3) :: iml(0:3)=(/(1.0,0.0),(0.0,-1.0),(-1.0,0.0),(0.0,1.0)/)
 complex(dpc),allocatable :: dbkcpk(:,:)
 type(cprj_type),allocatable :: cprj_k(:,:),cprj_kb(:,:),tcprj(:,:,:,:)
 type(cprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)
 type(cprj_type),allocatable :: bcprj(:,:),kcprj(:,:)

! ************************************************************************
!write(std_out,*)' JWZ Debug enter update_mmat '

 dtefield%twh(:,:,:,ikpt,:) = zero
 nband_k = dtefield%nband_occ

 if (berryopt == 5) then
   ihcg = dtefield%cgindex(ikpt,1)
   dtefield%hcg(:,ihcg+1:ihcg+nband_k*npw_k,:) = zero
   ABI_ALLOCATE(bffnl,(MAXVAL(npwarr),dimffnl,lmnmax,ntypat,6))
 end if

!=============================================================================
!update smatrix for neighboring k points. Note that S^{-1} is not used in this
!implementation, unlike what is done in electric field case.
!
!the following code is largely copied from cgwf, where the same thing is done
!=============================================================================

 mcg1_k = mpw*nband_k

 ABI_ALLOCATE(tcg,(2,mcg1_k,3,2))
 ABI_ALLOCATE(cg1_k,(2,mcg1_k))
 tcg(:,:,:,:) = zero
 
 mcg_q = mpw*mband*nspinor
 ABI_ALLOCATE(cgq_k,(2,mcg_q))
 ABI_ALLOCATE(sflag_k,(nband_k))
 ABI_ALLOCATE(pwind_k,(mpw))
 ABI_ALLOCATE(pwnsfac_k,(4,mpw))

 ABI_ALLOCATE(cprj_k,(natom,nband_k))
 ABI_ALLOCATE(cprj_kb,(natom,nband_k))
 ABI_ALLOCATE(dimlmn_srt,(natom))
 iatom = 0
 do itypat = 1, ntypat
   do iat = 1, gs_hamk%nattyp(itypat)
     iatom = iatom + 1
     dimlmn_srt(iatom)=pawtab(itypat)%lmn_size
   end do
 end do
 ABI_ALLOCATE(dimlmn,(natom))
 do iatom = 1, natom
   itypat = gs_hamk%typat(iatom)
   dimlmn(iatom)=dtefield%lmn_size(itypat)
 end do
 call cprj_alloc(cprj_k,0,dimlmn)
 call cprj_alloc(cprj_kb,0,dimlmn)

 ABI_ALLOCATE(tcprj,(natom,nband_k,3,2))

 do idir = 1, 3
   do ifor = 1, 2
     call cprj_alloc(tcprj(:,:,idir,ifor),0,dimlmn)
   end do
 end do

 ABI_ALLOCATE(cprj_fkn,(natom,dtefield%nband_occ))
 ABI_ALLOCATE(cprj_ikn,(natom,dtefield%nband_occ))
 call cprj_alloc(cprj_fkn,0,dimlmn)
 call cprj_alloc(cprj_ikn,0,dimlmn)

 ikptf = dtefield%i2fbz(ikpt)
 ikgf = dtefield%fkgindex(ikptf)  ! this is the shift for pwind

 icp=nband_k*(ikpt-1)

 call cprj_get(gs_hamk%atindx1,cprj_k,dtefield%cprj,natom,1,icp,ikpt,0,1,&
& nband_k,nkpt,mpi_enreg,natom,nband_k,nband_k,1,1,0)

 ABI_ALLOCATE(smat_k,(2,nband_k,nband_k))
 ABI_ALLOCATE(smat_inv,(2,nband_k,nband_k))
 ABI_ALLOCATE(smat_inv_all,(2,nband_k,nband_k,3,2))
 ABI_ALLOCATE(smat_k_paw,(2,nband_k,nband_k))

 job = 20 ! update overlap matrix and transfer cgq to cg1_k without S^{-1}
 shiftbd = 1 ! cg contains all bands
 ddkflag = 0 ! do not multiply wavefunctions at k by S^{-1}

 do idir = 1, 3

   dk(:) = dtefield%dkvecs(:,idir)

   do ifor = 1, 2

     ikpt2f = dtefield%ikpt_dk(ikptf,ifor,idir)
     if (dtefield%indkk_f2ibz(ikpt2f,6) == 1) then
       itrs = 10
     else
       itrs = 0
     end if
     ikpt2 = dtefield%indkk_f2ibz(ikpt2f,1)
     npw_k2 = npwarr(ikpt2)
     pwind_k(1:npw_k) = pwind(ikgf+1:ikgf+npw_k,ifor,idir)
     pwnsfac_k(1:2,1:npw_k) = pwnsfac(1:2,ikgf+1:ikgf+npw_k)
     sflag_k(:) = dtefield%sflag(:,ikpt,ifor,idir)
     smat_k(:,:,:) = dtefield%smat(:,:,:,ikpt,ifor,idir)

     if (mpi_enreg%paral_compil_kpt == 1) then
       icg1 = dtefield%cgqindex(2,ifor+2*(idir-1),ikpt) ! nsppol implicitly = 1 
       cgq_k(:,1:nband_k*nspinor*npw_k2) = &
&       cgq(:,icg1+1:icg1+nband_k*nspinor*npw_k2)
       idum1 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt) ! nsppol implicitly = 1
       pwnsfac_k(3:4,1:npw_k2) = pwnsfacq(1:2,idum1+1:idum1+npw_k2)
     else
       icg1 = dtefield%cgindex(ikpt2,1) ! nsppol implicitly = 1
       cgq_k(:,1:nband_k*nspinor*npw_k2) = &
&       cg(:,icg1+1:icg1+nband_k*nspinor*npw_k2)
       idum1 = dtefield%fkgindex(ikpt2f)
       pwnsfac_k(3:4,1:npw_k2) = pwnsfac(1:2,idum1+1:idum1+npw_k2)
     end if

     icp2=nband_k*(ikpt2-1)
     call cprj_get(gs_hamk%atindx1,cprj_kb,dtefield%cprj,natom,1,icp2,ikpt,0,1,&
&     nband_k,nkpt,mpi_enreg,natom,nband_k,nband_k,1,1,0)

     if (ikpt2 /= ikpt2f) then ! construct cprj_kb by symmetry
       call cprj_copy(cprj_kb,cprj_ikn)
       call sym_cprj_kn(cprj_fkn,cprj_ikn,dtefield%atom_indsym,dimlmn,-1,gs_hamk%indlmn,&
&       dtefield%indkk_f2ibz(ikpt2f,2),dtefield%indkk_f2ibz(ikpt2f,6),&
&       dtefield%fkptns(:,dtefield%i2fbz(ikpt2)),&
&       dtefield%lmax,dtefield%lmnmax,natom,dtefield%nband_occ,nspinor,&
&       dtefield%nsym,ntypat,gs_hamk%typat,dtefield%zarot)
       call cprj_copy(cprj_fkn,cprj_kb)
     end if

     call smatrix_k_paw(cprj_k,cprj_kb,dtefield,idir,ifor,natom,smat_k_paw,gs_hamk%typat)

     icg1 = 0
     call smatrix(cg,cgq_k,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,&
&     job,nband_k,mcg,mcg_q,mcg1_k,1,mpw,nband_k,&
&     npw_k,npw_k2,nspinor,pwind_k,pwnsfac_k,sflag_k,&
&     shiftbd,smat_inv,smat_k,smat_k_paw,gs_hamk%usepaw)
     dtefield%sflag(:,ikpt,ifor,idir) = sflag_k(:)
     dtefield%smat(:,:,:,ikpt,ifor,idir) = smat_k(:,:,:)

!    now cg1_k contains |u}_k,k+b>, and these have also been shifted by pwind so 
!    that they are expanded in exactly the same G vectors as |u_k>

     tcg(:,:,idir,ifor) = cg1_k(:,:)
     call cprj_copy(cprj_kb,tcprj(:,:,idir,ifor))

!    berryopt 5 needs ffnl for the neighbor k points
     if (berryopt == 5) then
!      Compute (k+G) vectors
       ikg2 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt) ! this may or may not be correct in // case
       ABI_ALLOCATE(kpg_k2,(npw_k2,nkpg))
       ABI_ALLOCATE(kg_k2,(3,npw_k2))
       kpoint2(1:3) = dtefield%fkptns(1:3,ikpt2f)
       kg_k2(:,1:npw_k2) = kg(:,1+ikg2:npw_k2+ikg2)
       if (nkpg>0) call mkkpg(kg_k2,kpg_k2,kpoint2,nkpg,npw_k2)
       ABI_ALLOCATE(ylm_k2,(npw_k2,mpsang*mpsang))
       do ilm=1,mpsang*mpsang
         ylm_k2(1:npw_k2,ilm)=ylm(1+ikg2:npw_k2+ikg2,ilm) ! check carefully!!
       end do
       tind = 2*idir-ifor+1
!      ider=0;idir=0;
       ABI_ALLOCATE(bffnl_,(1:npw_k2,1:dimffnl,1:lmnmax,1:ntypat))
       call mkffnl(psps%dimekb,dimffnl,psps%ekb,&
&       bffnl_(1:npw_k2,1:dimffnl,1:lmnmax,1:ntypat),&
&       psps%ffspl,gmet,gprimd,0,0,psps%indlmn,&
&       kg_k2,kpg_k2,kpoint2,psps%lmnmax,&
&       psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&       npw_k2,ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&       psps%usepaw,psps%useylm,ylm_k2,ylmgr)

!      bffnl_ was generated on k+G for k2, need to convert to expansion on G sphere for k
!      TODO : should pwnsfac_k be involved here? Maybe when symmetry is used (kptopt /= 3)?
       do ipw = 1, npw_k
         jpw = pwind_k(ipw)
         if(jpw > 0) bffnl(ipw,1:dimffnl,1:lmnmax,1:ntypat,tind) = &
&         bffnl_(jpw,1:dimffnl,1:lmnmax,1:ntypat)
       end do

       ABI_DEALLOCATE(kg_k2)
       ABI_DEALLOCATE(kpg_k2)
       ABI_DEALLOCATE(ylm_k2)
       ABI_DEALLOCATE(bffnl_)
     end if ! end bffnl computation

   end do ! end loop over ifor

 end do ! end loop over idir

!compute emat(2,kband,ikpt) = <u_k|H_k|u_k>

 sij_opt = 0 ! compute <G|H|u> only, not gsc
 cpopt = 2 ! cpopt in memory
 lambda = 0.0 ! no shift to H used here
 ndat = 1 ! number of FFTs to do in parallel
 tim_getghc=8 
 
 type_calc = 0 ! use entire Hamiltonian
!
 ABI_ALLOCATE(bwave,(2,npw_k))
 ABI_ALLOCATE(kwave,(2,npw_k))
 ABI_ALLOCATE(ghc,(2,npw_k))
 ABI_ALLOCATE(gvnlc,(2,npw_k))
 ABI_ALLOCATE(gsc,(2,npw_k*ndat*(sij_opt+1)/2))

 ghc(:,:) = zero
 gvnlc(:,:) = zero

 ABI_ALLOCATE(kcprj,(natom,nspinor*((cpopt+5)/5)))
 call cprj_alloc(kcprj,0,dimlmn_srt)
 ABI_ALLOCATE(bcprj,(natom,nspinor*((cpopt+5)/5)))
 call cprj_alloc(bcprj,0,dimlmn_srt)

 ABI_ALLOCATE(kg_k,(3,mpw))
 ikg = dtefield%kgindex(ikpt)
 kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)

 do kband = 1, nband_k

   kwave(:,:) = cg(:,icg+(kband-1)*npw_k+1:icg+kband*npw_k)
   bwave(:,:) = cg(:,icg+(kband-1)*npw_k+1:icg+kband*npw_k)

!  copy cprj for kband into kcprj structure, change order from input atom order
!  to atom type sort order as needed by getghc
   call cprj_get(gs_hamk%atindx,kcprj,cprj_k,natom,kband,0,ikpt,1,1,&
&   nband_k,1,mpi_enreg,natom,1,nband_k,1,1,0)

   call  getghc(cpopt,kwave,kcprj,dimffnl,ffnl,filstat,ghc,gsc,&
&   gs_hamk,gvnlc,kg_k,kinpw,lambda,lmnmax,&
&   matblk,mgfft,mpi_enreg,mpsang,mpssoang,&
&   natom,ndat,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,&
&   paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,type_calc,vlocal)

   dtefield%emat(1,kband,ikpt) = dot_product(bwave(1,:),ghc(1,:)) + dot_product(bwave(2,:),ghc(2,:)) 
   dtefield%emat(2,kband,ikpt) = dot_product(bwave(1,:),ghc(2,:)) - dot_product(bwave(2,:),ghc(1,:)) 

 end do

!compute omat(2,bband,kband,bdir,bfor,kdir,kfor) = <_bband,k+bdir*bfor|u_kband,k+kdir*kfor>
 ABI_ALLOCATE(omat,(2,nband_k,nband_k,3,2,3,2))
 omat(:,:,:,:,:,:,:) = zero
 do bdir = 1, 3
   do kdir = 1, 3
     if (bdir == kdir) cycle ! never need omat for bdir // kdir
     do bfor = 1, 2
       do kfor = 1, 2

         call smatrix_k_paw(tcprj(:,:,bdir,bfor),tcprj(:,:,kdir,kfor),dtefield,kdir,kfor,natom,&
&         smat_k_paw,gs_hamk%typat,bdir,bfor)

         do bband = 1, nband_k
           bwave(:,:) = tcg(:,(bband-1)*npw_k+1:bband*npw_k,bdir,bfor)

           do kband = 1, nband_k

             kwave(:,:) = tcg(:,(kband-1)*npw_k+1:kband*npw_k,kdir,kfor)

             dotri(1) = dot_product(bwave(1,:),kwave(1,:))+dot_product(bwave(2,:),kwave(2,:))
             dotri(2) = -dot_product(bwave(2,:),kwave(1,:))+dot_product(bwave(1,:),kwave(2,:))

             omat(1,bband,kband,bdir,bfor,kdir,kfor) = dotri(1)+smat_k_paw(1,bband,kband)
             omat(2,bband,kband,bdir,bfor,kdir,kfor) = dotri(2)+smat_k_paw(2,bband,kband)

           end do ! end loop over kband
         end do ! end loop over bband
       end do ! end loop over kfor
     end do ! end loop over bfor
   end do ! end loop over kdir
 end do ! end loop over bdir
!
!!compute chern_k(2,ikpt,idir)
!
 dtefield%chern_k(:,ikpt,:) = zero

 do idir = 1, 3
   bdir = mod(idir,3)+1
   kdir = mod(bdir,3)+1
   sfac = 0.25
   do istep = 1, 2
     do bsig = -1, 1, 2
       do ksig= -1, 1, 2
         bfor = (-bsig+3)/2; kfor = (-ksig+3)/2

         IA = cmplx(zero,zero)
         IB = cmplx(zero,zero)
         do nn = 1, nband_k
           do nnp = 1, nband_k
             IA1 = cmplx(dtefield%smat(1,nn,nnp,ikpt,bfor,bdir),&
&             dtefield%smat(2,nn,nnp,ikpt,bfor,bdir))
             IB1 = IA1
             do nnpp = 1, nband_k
               IA2 = cmplx(omat(1,nnp,nnpp,bdir,bfor,kdir,kfor),&
&               omat(2,nnp,nnpp,bdir,bfor,kdir,kfor))
               IA3 = cmplx(dtefield%smat(1,nn,nnpp,ikpt,kfor,kdir),&
&               dtefield%smat(2,nn,nnpp,ikpt,kfor,kdir))
               IA = IA + IA1*IA2*conjg(IA3)
               
               IB2 = cmplx(dtefield%smat(1,nnpp,nnp,ikpt,bfor,bdir),&
&               dtefield%smat(2,nnpp,nnp,ikpt,bfor,bdir))
               
               do nnppp = 1, nband_k
                 IB3 = cmplx(dtefield%smat(1,nnpp,nnppp,ikpt,kfor,kdir),&
&                 dtefield%smat(2,nnpp,nnppp,ikpt,kfor,kdir))
                 IB4 = cmplx(dtefield%smat(1,nn,nnppp,ikpt,kfor,kdir),&
&                 dtefield%smat(2,nn,nnppp,ikpt,kfor,kdir))
                 IB = IB + IB1*conjg(IB2)*IB3*conjg(IB4)
               end do ! end loop over nnppp
             end do ! end loop over nnpp
           end do ! end loop over nnp
         end do ! end loop over nn

         dtefield%chern_k(1,ikpt,idir) = &
&         dtefield%chern_k(1,ikpt,idir) + sfac*bsig*ksig*real(IA-IB)
         dtefield%chern_k(2,ikpt,idir) = &
&         dtefield%chern_k(2,ikpt,idir) + sfac*bsig*ksig*aimag(IA-IB)

       end do ! end loop over ksig
     end do ! end loop over bsig
     
     sfac = -sfac
     idum1 = bdir; bdir = kdir; kdir = idum1

   end do ! loop over istep

 end do ! loop over idir

 
!update dtefield%twh, giving Hamiltonian matrix element between 
!different k points 
!<u_n,k+sig_b b_b|H_k|u_m,k+sig_k b_k>

 dtefield%twh(:,:,:,ikpt,:) = zero
!twh(2,nband_occ,nband_occ,nkpt*nsppol,tind)
!!
 sij_opt = 0 ! compute <G|H|u> only, not gsc
 cpopt = 2 ! use cprj in memory
 lambda = 0.0 ! no shift to H used here
 ndat = 1 ! number of FFTs to do in parallel
 tim_getghc=8 ! as in cgwf.F90

 do kdir = 1, 3
   do ksig = -1, 1, 2
     kfor = (-ksig+3)/2
     do kband = 1, nband_k

!      !     copy cprj for kband into kcprj structure, change order from input atom order
!      !     to atom type sort order as needed by getghc
       call cprj_get(gs_hamk%atindx,kcprj,tcprj(:,:,kdir,kfor),natom,kband,0,ikpt,1,1,&
&       nband_k,1,mpi_enreg,natom,1,nband_k,1,1,0)
       kwave(:,:) = tcg(:,(kband-1)*npw_k+1:kband*npw_k,kdir,kfor)
       type_calc = 3 ! use kinetic + local part only
!      !   note that tilde cg wavefunctions tcg have been adjusted by smatrix.F90 such that they are
!      !   expanded on the same set of G vectors as |u_nk>, therefore use of kg_k is still correct
       call getghc(cpopt,kwave,kcprj,dimffnl,ffnl,filstat,ghc,gsc,&
&       gs_hamk,gvnlc,kg_k,kinpw,lambda,lmnmax,&
&       matblk,mgfft,mpi_enreg,mpsang,mpssoang,&
&       natom,ndat,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,&
&       paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,type_calc,vlocal)
       do bdir = 1, 3
         if (bdir == kdir) cycle ! never need bdir // kdir case
         do bsig = -1, 1, 2
           bfor = (-bsig+3)/2
           bdx = 2*bdir-bfor+1; kdx = 2*kdir-kfor+1
           dij_ind = dtefield%indhk(bdx,kdx)

!          contract D_bk cp_k
           ABI_ALLOCATE(dbkcpk,(natom,psps%lmnmax))
           dbkcpk(:,:) = czero
           do iatom = 1, natom
             itypat = gs_hamk%typat(iatom)           
             do blmn = 1, pawtab(itypat)%lmn_size
               do klmn = 1, pawtab(itypat)%lmn_size
                 cpk = cmplx(tcprj(iatom,kband,kdir,kfor)%cp(1,klmn),&
&                 tcprj(iatom,kband,kdir,kfor)%cp(2,klmn))
                 dij = cmplx(dtefield%twdij(1,blmn,klmn,iatom,dij_ind),&
&                 dtefield%twdij(2,blmn,klmn,iatom,dij_ind))
                 dbkcpk(iatom,blmn) = dbkcpk(iatom,blmn) + dij*cpk
               end do ! end loop over klmn
             end do ! end loop over blmn
           end do ! end loop over atoms

           if (berryopt == 5) then ! compute H_k2|u_k3> for use in cgwf

             do ipw = 1, npw_k
               onfac = czero
               do iatom = 1, natom
                 itypat = gs_hamk%typat(iatom)
                 do blmn = 1, pawtab(itypat)%lmn_size
                   il = gs_hamk%indlmn(1,blmn,itypat)
                   bffnlfac = four_pi*iml(mod(il,4))/sqrt(ucvol)
                   onfac = onfac + bffnlfac*bffnl(ipw,1,blmn,itypat,bdx)*dbkcpk(iatom,blmn)
                 end do ! end loop over blmn
               end do ! end loop over iatom

               dtefield%hcg(1,ihcg+(kband-1)*npw_k+ipw,dij_ind) = ghc(1,ipw)+real(onfac)
               dtefield%hcg(2,ihcg+(kband-1)*npw_k+ipw,dij_ind) = ghc(2,ipw)+aimag(onfac)

             end do  ! end loop over ipw

           end if ! end berryopt == 5

           do bband = 1, nband_k
             bwave(:,:) = tcg(:,(bband-1)*npw_k+1:bband*npw_k,bdir,bfor)

!            contract with cp_b to get full onsite term
             onsite = cmplx(0.d0,0.d0)
             do iatom = 1, natom
               itypat = gs_hamk%typat(iatom)
               do blmn = 1, pawtab(itypat)%lmn_size
                 cpb = cmplx(tcprj(iatom,bband,bdir,bfor)%cp(1,blmn),&
                 tcprj(iatom,bband,bdir,bfor)%cp(2,blmn))
                 onsite = onsite + conjg(cpb)*dbkcpk(iatom,blmn)
                 
               end do ! end loop over blmn
             end do ! end loop over atoms

!            !           ghc contains (K.E. + v_eff)|ket> , no nonlocal parts
             dotri(1) = dot_product(bwave(1,:),ghc(1,:))+dot_product(bwave(2,:),ghc(2,:))
             dotri(2) = -dot_product(bwave(2,:),ghc(1,:))+dot_product(bwave(1,:),ghc(2,:))

             dtefield%twh(1,bband,kband,ikpt,dij_ind) = dotri(1)+real(onsite)
             dtefield%twh(2,bband,kband,ikpt,dij_ind) = dotri(2)+aimag(onsite)
             
           end do ! end loop over bband
           ABI_DEALLOCATE(dbkcpk)

         end do ! end loop over bsig
       end do ! end loop over bdir
     end do ! end loop over kband
   end do ! end loop over ksig
 end do ! end loop over kdir

!
!
!!compute mag_k(2,ikpt,idir)
!
 dtefield%mag_k(:,ikpt,:) = zero

 do idir = 1, 3
   bdir = mod(idir,3)+1
   kdir = mod(bdir,3)+1
   sfac = 0.25
   do istep = 1, 2
     do bsig = -1, 1, 2
       do ksig= -1, 1, 2
         bfor = (-bsig+3)/2; kfor = (-ksig+3)/2
         bdx = 2*bdir-bfor+1; kdx = 2*kdir-kfor+1
         dij_ind = dtefield%indhk(bdx,kdx)

         IIA = cmplx(zero,zero)
         IIIA = cmplx(zero,zero)
         do nn = 1, nband_k
           IIIA1 = cmplx(dtefield%emat(1,nn,ikpt),dtefield%emat(2,nn,ikpt))
           do nnp = 1, nband_k
             IIA1 = cmplx(dtefield%smat(1,nn,nnp,ikpt,bfor,bdir),&
&             dtefield%smat(2,nn,nnp,ikpt,bfor,bdir))
             IIIA2 = cmplx(dtefield%smat(1,nn,nnp,ikpt,kfor,kdir),&
&             dtefield%smat(2,nn,nnp,ikpt,kfor,kdir))
             do nnpp = 1, nband_k
               IIA2 = cmplx(dtefield%twh(1,nnp,nnpp,ikpt,dij_ind),&
&               dtefield%twh(2,nnp,nnpp,ikpt,dij_ind))
               IIA3 = cmplx(dtefield%smat(1,nn,nnpp,ikpt,kfor,kdir),&
&               dtefield%smat(2,nn,nnpp,ikpt,kfor,kdir))
               IIIA3 = cmplx(omat(1,nnp,nnpp,kdir,kfor,bdir,bfor),&
&               omat(2,nnp,nnpp,kdir,kfor,bdir,bfor))
               IIIA4 = cmplx(dtefield%smat(1,nn,nnpp,ikpt,bfor,bdir),&
&               dtefield%smat(2,nn,nnpp,ikpt,bfor,bdir))
               IIA  = IIA   + IIA1*IIA2*conjg(IIA3)
               IIIA = IIIA3 - IIIA1*IIIA2*IIIA3*conjg(IIIA4)
             end do ! end loop on nnpp
           end do ! end loop on nnp
         end do ! end loop on nn

         dtefield%mag_k(1,ikpt,idir) = &
&         dtefield%mag_k(1,ikpt,idir) + sfac*bsig*ksig*real(IIA+IIIA)
         dtefield%mag_k(2,ikpt,idir) = &
&         dtefield%mag_k(2,ikpt,idir) + sfac*bsig*ksig*aimag(IIA+IIIA)

       end do ! end loop over ksig
     end do ! end loop over bsig
     
     sfac = -sfac
     idum1 = bdir; bdir = kdir; kdir = idum1

   end do ! loop over istep
 end do ! loop over idir

 ABI_DEALLOCATE(bwave)
 ABI_DEALLOCATE(kwave)
 ABI_DEALLOCATE(ghc)
 ABI_DEALLOCATE(gvnlc)
 ABI_DEALLOCATE(gsc)
 ABI_DEALLOCATE(omat)
 call cprj_free(kcprj)
 ABI_DEALLOCATE(kcprj)
 call cprj_free(bcprj)
 ABI_DEALLOCATE(bcprj)
!
!write(std_out,*)' JWZ Debug deallocating update_mmat '

 ABI_DEALLOCATE(tcg)
 ABI_DEALLOCATE(cg1_k)
 ABI_DEALLOCATE(cgq_k)
 ABI_DEALLOCATE(sflag_k)
 ABI_DEALLOCATE(pwind_k)
 ABI_DEALLOCATE(pwnsfac_k)
 call cprj_free(cprj_k)
 ABI_DEALLOCATE(cprj_k)
 call cprj_free(cprj_kb)
 ABI_DEALLOCATE(cprj_kb)
 do idir =1, 3
   do ifor = 1, 2
     call cprj_free(tcprj(:,:,idir,ifor))
   end do
 end do
 ABI_DEALLOCATE(tcprj)

 call cprj_free(cprj_fkn)
 ABI_DEALLOCATE(cprj_fkn)
 call cprj_free(cprj_ikn)
 ABI_DEALLOCATE(cprj_ikn)

 if (berryopt == 5)  then
   ABI_DEALLOCATE(bffnl)
 end if


 ABI_DEALLOCATE(smat_k)
 ABI_DEALLOCATE(smat_inv)
 ABI_DEALLOCATE(smat_inv_all)
 ABI_DEALLOCATE(smat_k_paw)
 ABI_DEALLOCATE(dimlmn)
 ABI_DEALLOCATE(dimlmn_srt)
 ABI_DEALLOCATE(kg_k)

!write(std_out,*)' JWZ Debug leaving update_mmat '

end subroutine update_mmat
!!***
