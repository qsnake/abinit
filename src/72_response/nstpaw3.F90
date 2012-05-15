!!{\src2tex{textfont=tt}}
!!****f* ABINIT/nstpaw3
!! NAME
!! nstpaw3
!!
!! FUNCTION
!! This routine compute the additional PAW contributions
!! to the non-stationary expression for the second derivative of the total energy,
!! for a whole row of mixed derivatives (including diagonal terms contributing
!! to non-stationnary 2nd-order total energy).
!! Compared with NC-pseudopotentials, these contributions include
!! the changes of the overlap between 0-order wave-functions.
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions at k
!!  cgq(2,mpw1*nspinor*mband*mkqmem*nsppol)=pw coefficients of GS wavefunctions at k+q.
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF wavefunctions at k,q.
!!  cplex=if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!                                                  projected with non-local projectors
!!  cprjq(natom,nspinor*mband*mkqmem*nsppol*usecprj)= wave functions at k+q
!!                                                    projected with non-local projectors
!!  dimcprj(natom*usepaw)=array of dimensions of arrays cprj, cprjq (ordered by atom-type)
!!  docckqde(mband*nkpt_rbz*nsppol)=derivative of occkq wrt the energy
!!  doccde_rbz(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigenq(mband*nkpt_rbz*nsppol)=GS eigenvalues at k+q (hartree)
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  eigen1(2*mband*mband*nkpt_rbz*nsppol)=1st-order eigenvalues at k,q (hartree)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  idir=direction of the perturbation
!!  indkpt1(nkpt_rbz)=non-symmetrized indices of the k-points
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  ipert=type of the perturbation
!!  irrzon1(nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data for RF symmetries
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points in the reduced BZ
!!  kxc(nfftf,nkxc)=exchange and correlation kernel
!!  mgfftf=maximum size of 1D FFTs for the "fine" grid (see NOTES in respfn.F90)
!!  mpert =maximum number of ipert
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  nfftf=(effective) number of FFT grid points (for this proc) for the "fine" grid
!!  ngfftf(1:18)=integer array with FFT box dimensions and other for the "fine" grid
!!  nhat1(cplex*nfftf,nspden*psps%usepaw)=1st-order compensation charge density (PAW)
!!  nkpt=number of k points in the full BZ
!!  nkpt_rbz=number of k points in the reduced BZ for this perturbation
!!  nkxc=second dimension of the kxc array
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym1=number of symmetry elements in space group consistent with i perturbation
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used otherwise, cplex*nfftf
!!  occkq(mband*nkpt_rbz*nsppol)=occupation number for each band at each k+q point of the reduced BZ
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band and k in the reduced BZ
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh for the GS
!!  paw_an1(natom) <type(paw_an_type)>=1st-order paw arrays given on angular mesh for the perturbation (j1)
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  paw_ij1(natom) <type(paw_ij_type)>=1st-order paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawang1 <type(pawang_type)>=pawang datastructure containing only the symmetries preserving the perturbation
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid for the GS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data for the GS
!!  pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  phnons1(2,nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic transl. phases, for RF symmetries
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  ph1df(2,3*(2*mgfftf+1)*natom)=one-dimensional structure factor information for the "fine" grid
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhor1(cplex*nfftf,nspden)=RF electron density in electrons/bohr**3.
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symaf1(nsym1)=anti(ferromagnetic) part of symmetry operations
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  symrl1(3,3,nsym1)=symmetry operations in real space in terms
!!  ucvol=unit cell volume in bohr**3.
!!  usecprj= 1 if cprj, cprjq arrays are stored in memory
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  useylmgr1= 1 if ylmgr1 array is allocated
!!  vhartr1(cplex*nfft)=1-order Hartree potential
!!  vpsp1(cplex*nfftf)=first-order derivative of the ionic potential
!!  vtrial(nfftf,nspden)=GS potential (Hartree).
!!  vtrial1(cplex*nfftf,nspden)= RF 1st-order potential (Hartree).
!!  vxc(nfftf,nspden)=XC GS potential
!!  wffnow=struct info for 1st-order WF disk file
!!  wfftgs,wfftkq=struct info for ground-state WF disk files
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point in the reduced BZ
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density, see n3xccc
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k+q point
!!  ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics at k+q
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert)=flags for each element of the 2DTE (=1 if computed)
!!  d2lo(2,3,mpert,3,mpert)=local contributions to the 2DTEs
!!  d2nl(2,3,mpert,3,mpert)=non-local contributions to the 2DTEs
!!  d2ovl(2,3,mpert,3,mpert)=overlap contributions to the 2DTEs
!!  eovl1=1st-order change of wave-functions overlap, part of 2nd-order energy
!!        PAW only - Eq(79) and Eq(80) of PRB 78, 035105 (2008)
!!
!! NOTES
!   We perform here the computation of
!!     delta_u^(j1)=-1/2 Sum_{j}[<u0_k+q_j|S^(j1)|u0_k_i>.|u0_k+q_j>]
!      see PRB 78, 035105 (2008), Eq. (42)
!!
!! PARENTS
!!      scfcv3
!!
!! CHILDREN
!!      accrho3,appdig,atm2fft3,clsopn,cprj_alloc,cprj_copy,cprj_diskinit_r
!!      cprj_free,cprj_get,destroy_hamiltonian,destroy_paw_ij,dotprod_g
!!      dotprod_vn,fftpac,getdc1,getgh1c,hdr_skip,init_hamiltonian,init_paw_ij
!!      initylmg,kpg3,leave_test,mkffnl,mkkpg,mkvxc3,nullify_paw_ij,occeig
!!      pawfrnhat,pawmkrho,pawnstd2e,ph1d3d,projbd,rdnpw,rhoij_alloc,rhoij_free
!!      sphereboundary,sygra3,symrhg,timab,wffclose,wffkg,wffopen
!!      wffreaddatarec,wffreadnpwrec,wffreadskipk,wffreadskiprec,wrtout
!!      xcomm_world,xdefineoff,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine nstpaw3(blkflg,cg,cgq,cg1,cplex,cprj,cprjq,dimcprj,docckqde,doccde_rbz,dtfil,dtset,d2lo,d2nl,d2ovl,&
&                  eigenq,eigen0,eigen1,eovl1,gmet,gprimd,gsqcut,idir,indkpt1,indsy1,ipert,irrzon1,istwfk_rbz,&
&                  kg,kg1,kpt_rbz,kxc,mgfftf,mpert,mpi_enreg,mpw,mpw1,nband_rbz,ncpgr,nfftf,ngfftf,nhat1,&
&                  nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,nspinor,nsppol,nsym1,n3xccc,occkq,occ_rbz,&
&                  paw_an,paw_an1,paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoij,&
&                  pawrhoij1,pawtab,phnons1,ph1d,ph1df,psps,rhor1,rmet,rprimd,symaf1,symrc1,symrl1,&
&                  ucvol,usecprj,usexcnhat,useylmgr1,vhartr1,vpsp1,vtrial,vtrial1,vxc,&
&                  wffnow,wfftgs,wfftkq,wtk_rbz,xccc3d1,xred,ylm,ylm1,ylmgr1)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_paw_toolbox
 use m_wffile

 use m_hamiltonian, only : init_hamiltonian, destroy_hamiltonian

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nstpaw3'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_53_spacepar
 use interfaces_56_recipspace
 use interfaces_56_xc
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_62_occeig
 use interfaces_65_nonlocal
 use interfaces_65_psp
 use interfaces_66_paw
 use interfaces_66_wfs
 use interfaces_67_common
 use interfaces_72_response, except_this_one => nstpaw3
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mgfftf,mpert,mpw,mpw1,ncpgr,nfftf,nkpt,nkpt_rbz,nkxc
 integer,intent(in) :: nspden,nspinor,nsppol,nsym1,n3xccc,usecprj,usexcnhat,useylmgr1
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(out) :: eovl1
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang,pawang1
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow,wfftgs,wfftkq
!arrays
 integer,intent(in) :: dimcprj(dtset%natom*psps%usepaw),nband_rbz(nkpt_rbz*nsppol)
 integer,intent(in) :: indkpt1(nkpt_rbz),indsy1(4,nsym1,dtset%natom)
 integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))
 integer,intent(in) :: istwfk_rbz(nkpt_rbz),kg(3,mpw*dtset%mkmem),kg1(3,mpw1*dtset%mk1mem)
 integer,intent(in) :: ngfftf(18),npwarr(nkpt_rbz),npwar1(nkpt_rbz)
 integer,intent(in) :: symaf1(nsym1),symrc1(3,3,nsym1),symrl1(3,3,nsym1)
 integer,intent(out) :: blkflg(3,mpert,3,mpert)
 real(dp),intent(in) :: cg(2,mpw*nspinor*dtset%mband*dtset%mkmem*nsppol)
 real(dp),intent(in) :: cgq(2,mpw1*nspinor*dtset%mband*dtset%mkqmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*dtset%mband*dtset%mk1mem*nsppol)
 real(dp),intent(in) :: docckqde(dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: doccde_rbz(dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigenq(dtset%mband*nkpt_rbz*nsppol),eigen0(dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),kpt_rbz(3,nkpt_rbz)
 real(dp),intent(in) :: kxc(nfftf,nkxc),nhat1(cplex*nfftf,nspden*psps%usepaw)
 real(dp),intent(in) :: occkq(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: occ_rbz(dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
 real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*dtset%natom)
 real(dp),intent(in) :: rhor1(cplex*nfftf,nspden),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: vhartr1(cplex*nfftf),vtrial1(cplex*nfftf,nspden),vxc(nfftf,nspden)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xred(3,dtset%natom)
 real(dp),intent(in) :: ylm(mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*dtset%mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr1(mpw1*dtset%mk1mem,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 real(dp),target,intent(in) :: vpsp1(cplex*nfftf),vtrial(nfftf,nspden),xccc3d1(cplex*n3xccc)
 real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert),d2nl(2,3,mpert,3,mpert),d2ovl(2,3,mpert,3,mpert)
 type(cprj_type),intent(in) :: cprj(dtset%natom,nspinor*dtset%mband*dtset%mkmem*nsppol*usecprj)
 type(cprj_type),intent(in) :: cprjq(dtset%natom,nspinor*dtset%mband*dtset%mkqmem*nsppol*usecprj)
 type(paw_an_type),intent(in) :: paw_an(dtset%natom*psps%usepaw)
 type(paw_an_type),intent(inout) :: paw_an1(dtset%natom*psps%usepaw)
 type(paw_ij_type),intent(in) :: paw_ij(dtset%natom*psps%usepaw)
 type(paw_ij_type),intent(inout) :: paw_ij1(dtset%natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij1(dtset%natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourwf=18,tim_getgh1c=3,tim_projbd=3
 integer :: bd2tot_index,bdtot_index,berryopt,bufsz,counter,cplx
 integer :: ddkcase,dimdij,dimdij1,dime1kb,dimffnl,dimffnl1
 integer :: dimffnl1_idir1,dimphkxred,dimylmgr1,formeig,ia,iatom,iband,ibg,ibgq,ibg1
 integer :: icg,icgq,icg1,ider,idir0,idir1,ierr,ii,ikg,ikg1,ikpt,ikpt_,ikpt_me,ilmn
 integer :: iorder_cprj,iorder_cprj1,ipert1,isp,ispden,isppol,istwf_k,itypat,jband,jj,kdir1
 integer :: kpert1,master,matblk,mcgq,mcprjq,mdir1,me,mpert1,nband_,nband_k,nband_kocc
 integer :: need_ylmgr1,nfftot,nkpg,nkpg1,nkpt_me,npw_,npw_k,npw1_k,nskip,nsp,nspinor_
 integer :: nvxc1,nzlmopt_ipert,nzlmopt_ipert1,optfr,optlocal,optn,optn2,optnl
 integer :: option,optv,opt_gvnl1,sij_opt,spaceworld,t_iostat,usee1kb,usevnl,v_size,wfcorr
 real(dp) :: arg,doti,dotr,dot1i,dot1r,dot2i,dot2r,invocc,lambda,wtk_k
 logical :: has_dcwf,has_dcwf2,has_drho,has_ddk_file,has_ipert_dijh,has_ipert_vxc
 logical :: is_metal,is_metal_or_qne0,need_ddk_file,need_wfk,need_wfq,need_wf1,t_test
 character(len=500) :: msg
 character(len=fnlen) :: fiwfddk(3)
 type(gs_hamiltonian_type) :: gs_hamkq
!arrays
 integer :: ddkfil(3),ikpt_fbz(3),ikpt_fbz_previous(3),nband_tmp(1)
 integer :: npwar1_tmp(1),pspso_typ(1),skipddk(3)
 integer,allocatable :: gbound(:,:),indlmn_typ(:,:,:),jpert1(:),jdir1(:)
 integer,allocatable :: kg1_k(:,:),kg_k(:,:),nlmn(:)
 real(dp) :: dum1(1,1),dum2(1,1),dum3(1,1),dum4(1,1,1),epawnst(2),kpoint(3)
 real(dp) :: sumelfd(2),summgfd(2),tsec(2),ylmgr_dum(1,3,1)
 real(dp),allocatable :: buffer(:),cgq_disk(:,:),ch1c(:,:,:,:),cs1c(:,:,:,:)
 real(dp),allocatable :: cwave0(:,:),cwavef(:,:),dcwavef(:,:)
 real(dp),allocatable :: doccde_k(:),doccde_kq(:)
 real(dp),allocatable :: dnhat1(:,:),drhoaug1(:,:,:,:)
 real(dp),allocatable :: drhor1(:,:),drho1wfg(:,:),drho1wfr(:,:,:)
 real(dp),allocatable :: d2nl_elfd(:,:),d2nl_mgfd(:,:),dkinpw(:)
 real(dp),allocatable :: d2nl_k(:,:),d2ovl_drho(:,:,:,:,:),d2ovl_k(:,:),eig_k(:),eig_kq(:),eig1_k(:)
 real(dp),allocatable :: ekb_typ(:,:,:),e1kbfr(:,:,:),ffnlk(:,:,:,:)
 real(dp),allocatable,target :: ffnl1(:,:,:,:)
 real(dp),allocatable :: gh1(:,:),gs1(:,:),gvnl1(:,:),kinpw1(:),kpg_k(:,:),kpg1_k(:,:)
 real(dp),allocatable :: occ_k(:),occ_kq(:),phkxred(:,:),ph3d(:,:,:),rhotmp(:,:),rocceig(:,:),sij_typ(:,:)
 real(dp),allocatable :: ylm_k(:,:),ylm1_k(:,:),ylmgr1_k(:,:,:),vtmp1(:,:),vxc10(:,:),work(:,:,:)
 real(dp),pointer :: ffnlkq(:,:,:),ffnl1_idir1(:,:,:,:),vpsp1_idir1(:),vtmp(:,:),xccc3d1_idir1(:)
 type(cprj_type),allocatable :: cprjq_disk(:,:),dcwaveprj(:,:)
 type(cprj_type),allocatable,target :: cwaveprj0(:,:)
 type(cprj_type),pointer :: cwaveprj0_idir1(:,:)
 type(paw_ij_type),allocatable :: paw_ij10(:,:),paw_ij1_tmp(:)
 type(pawrhoij_type),allocatable :: pawdrhoij1(:,:)
 type(wffile_type) :: wffddk(3)

! *********************************************************************

 DBG_ENTER("COLL")

!Keep track of total time spent in nstpaw3
 call timab(566,1,tsec)

!Only valid for PAW
 if (psps%usepaw==0) then
   msg='  should be called with usepaw=1 !'
   MSG_BUG(msg)
 end if

!Not valid for strain perturbation
 if (ipert==dtset%natom+3.or.ipert==dtset%natom+4) then
   msg='  not yet valid for strain perturbation !'
   MSG_BUG(msg)
 end if

!Not valid for magnetic field perturbation
 if (ipert==dtset%natom+5) then
   msg='  not yet valid for magnetic field perturbation !'
   MSG_BUG(msg)
 end if

!Not valid for PrintBandByBand
 if (dtset%prtbbb/=0) then
   msg='  not yet valid for prtbbb/=0 !'
   MSG_BUG(msg)
 end if

!Test on FFT grid sizes
 if (pawfgr%nfft/=nfftf) then
   msg='  wrong values for nfft, nfftf !'
   MSG_BUG(msg)
 end if

!Init data for parallelism
 master=0
 call xcomm_world(mpi_enreg,spaceworld)
 call xme_init(mpi_enreg,me)
 nkpt_me=nkpt_rbz
 if(mpi_enreg%paral_compil_kpt==1)then
   nkpt_me=0
   do isppol=1,nsppol
     do ikpt=1,nkpt_rbz
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       if (minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))==0) nkpt_me=nkpt_me+1
     end do
   end do
 end if

!Check for d/dk files (needed to compute electric field perturbations)
 ddkfil(:)=0
 do idir1=1,3
   ddkcase=idir1+dtset%natom*3
   call appdig(ddkcase,dtfil%fnamewffddk,fiwfddk(idir1))
!  Check that ddk file exists
   inquire(file=fiwfddk(idir1),iostat=t_iostat,exist=t_test)
   if (t_iostat/=0) then
     write(unit=msg,fmt='(5a,i8)' ) '  Check for existence of file ',trim(fiwfddk(idir1)),',',ch10,&
&     '  but INQUIRE statement returns error code',t_iostat
     MSG_ERROR(msg)
   else if (t_test) then
     ddkfil(idir1)=20+idir1 ! Note the use of unit numbers 21, 22 and 23
   end if
 end do
 has_ddk_file=(any(ddkfil(:)>0))

!Define the set of perturbations (j1)=(ipert1,idir1)
!The first perturbation must be (j1)=(j2)=(ipert,idir)
!because we need to compute <g|H^(j2)-Eps.S^(j2)|u0> first.
 if (ipert/=dtset%natom+1) then
   mpert1=0
   ABI_ALLOCATE(jpert1,(mpert))
   if ((ipert/=dtset%natom+2.and.ipert/=dtset%natom+5).or.has_ddk_file) then
     mpert1=mpert1+1;jpert1(mpert1)=ipert
   end if
   do ipert1=1,mpert
     if (ipert1/=ipert.and.&
&     (ipert1<=dtset%natom.or.&
&     ((ipert1==dtset%natom+2).and.has_ddk_file))) then
!      &      ((ipert1==dtset%natom+2.or.ipert1==dtset%natom+5).and.has_ddk_file))) then
       mpert1=mpert1+1;jpert1(mpert1)=ipert1
     end if
   end do
 else
   mpert1=1
   ABI_ALLOCATE(jpert1,(mpert1))
   jpert1(1)=dtset%natom+1
 end if
 mdir1=3
 ABI_ALLOCATE(jdir1,(mdir1))
 jdir1(1:3)= (/ (idir1,idir1=1,3) /)
 jdir1(1)=idir;jdir1(idir)=1

!Open ddk WF file(s)
 if (has_ddk_file) then
   do kdir1=1,mdir1
     idir1=jdir1(kdir1)
     if (ddkfil(idir1)/=0) then
       write(msg, '(a,a)') '-open ddk wf file :',fiwfddk(idir1)
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       call WffOpen(dtset%accesswff,spaceworld,fiwfddk(idir1),ierr,wffddk(idir1),master,me,ddkfil(idir1))
     end if
   end do
 end if

!Zero only portion of matrix to be computed here
 d2nl (:,:,1:dtset%natom+2,idir,ipert)=zero
 d2nl (:,:,  dtset%natom+5,idir,ipert)=zero
 d2ovl(:,:,1:dtset%natom+2,idir,ipert)=zero
 d2ovl(:,:,  dtset%natom+5,idir,ipert)=zero
 if (ipert==dtset%natom+3.or.ipert==dtset%natom+4) then
   d2nl (2,:,:,idir,ipert)=zero
   d2ovl(2,:,:,idir,ipert)=zero
 end if

!Update list of computed matrix elements
 do kpert1=1,mpert1
   ipert1=jpert1(kpert1)
   do kdir1=1,mdir1
     idir1=jdir1(kdir1)
     if ((ipert1<=dtset%natom).or.&
&     (ipert1==dtset%natom+1.and.((ddkfil(idir1)/=0).or.(dtset%rfdir(idir1)/=0.and.idir1<=idir))).or.&
&     ((ipert1==dtset%natom+2.or.ipert1==dtset%natom+5).and.ddkfil(idir1)/=0)) then
       blkflg(idir1,ipert1,idir,ipert)=1
     end if
   end do
 end do
 blkflg(:,dtset%natom+5,idir,ipert)=0 ! ipert=natom+5 to be activated later

!Initialize most of the (1st-order) Hamiltonian
!1) Allocate all arrays and initialize quantities that do not depend on k and spin.
!2) Perform the setup needed for the non-local factors:
 call init_hamiltonian(gs_hamkq,psps,paw_ij,pawtab,nspinor,nspden,dtset%natom,dtset%ntypat,dtset%typat,xred,&
& dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,ph1d=ph1d,use_gpu_cuda=dtset%use_gpu_cuda)

!Allocate on-site Hxc potentials for the ipert perturbation
 has_ipert_dijh=(paw_ij1(1)%has_dijhartree>0)
 has_ipert_vxc =(paw_an1(1)%has_vxc>0)
 if (ipert/=dtset%natom+1) then
   if (.not.(has_ipert_dijh.and.has_ipert_vxc)) then
     do iatom=1,dtset%natom
       itypat=dtset%typat(iatom)
       if (paw_ij1(iatom)%has_dijhartree==0) then
         ABI_ALLOCATE(paw_ij1(iatom)%dijhartree,(cplex*pawtab(itypat)%lmn2_size))
         paw_ij1(iatom)%has_dijhartree=1
       end if
       if (paw_an1(iatom)%has_vxc==0) then
         if (dtset%pawxcdev==0) v_size=paw_an1(iatom)%angl_size
         if (dtset%pawxcdev/=0) v_size=paw_an1(iatom)%lm_size
         ABI_ALLOCATE(paw_an1(iatom)%vxc1 ,(cplex*pawtab(itypat)%mesh_size,v_size,paw_an1(iatom)%nspden))
         ABI_ALLOCATE(paw_an1(iatom)%vxct1,(cplex*pawtab(itypat)%mesh_size,v_size,paw_an1(iatom)%nspden))
         paw_an1(iatom)%has_vxc=1
       end if
     end do
   end if
 end if

!Variables common to all perturbations
 arg=maxval(occ_rbz)-minval(occ_rbz)
 is_metal=((dtset%occopt>=3.and.dtset%occopt<=8).or.(abs(arg)>tol8))
 is_metal_or_qne0=((is_metal).or.(dtset%qptn(1)**2+dtset%qptn(2)**2+dtset%qptn(3)**2>=tol14))
 ABI_ALLOCATE(ch1c,(2,dtset%mband,dtset%mband,nkpt_me))
 ch1c(:,:,:,:)=zero
 if (is_metal_or_qne0) then
   ABI_ALLOCATE(cs1c,(2,dtset%mband,dtset%mband,nkpt_me))
   cs1c(:,:,:,:)=zero
 end if
 ABI_ALLOCATE(d2ovl_drho,(2,3,mpert,3,mpert))
 d2ovl_drho=zero
 if (dtset%pawnzlm==0) then
   nzlmopt_ipert=0;nzlmopt_ipert1=0
 else
   nzlmopt_ipert=1;if (dtset%nstep<2) nzlmopt_ipert=-1
   nzlmopt_ipert1=-1
 end if

!LOOP OVER PERTURBATION TYPES (j1)
 do kpert1=1,mpert1
   ipert1=jpert1(kpert1)

!  Flag for use of DDK file
   need_ddk_file=(has_ddk_file.and.(ipert1==dtset%natom+1.or.ipert1==dtset%natom+2.or.ipert1==dtset%natom+5))

!  We want to compute delta_u^(j1))=-1/2 Sum_{j}[<u0_k+q_j|S^(j1)|u0_k_i>.|u0_k+q_j>]
!  see PRB 78, 035105 (2008), Eq. (42)
   has_dcwf =(ipert1/=dtset%natom+2.and.ipert1/=dtset%natom+5)
   has_dcwf2=(ipert /=dtset%natom+2.and.ipert /=dtset%natom+5)
   has_drho =(has_dcwf.and.ipert1/=dtset%natom+1)

!  Select which WF are needed
   need_wfk=(.true.)
   need_wfq=(has_dcwf)
   need_wf1=(.true.)

!  Allocate arrays depending on the perturbation (j1): (factors for NL Hamiltonian and overlap)
   usee1kb=0;dime1kb=0
   if (ipert1/=dtset%natom+1) then
     usee1kb=1;dime1kb=psps%dimekb*paw_ij1(1)%cplex_dij
   end if
   ABI_ALLOCATE(ekb_typ,(gs_hamkq%dimekb1,1,nspinor**2))
   ABI_ALLOCATE(sij_typ,(gs_hamkq%dimekb1,1))
   ABI_ALLOCATE(indlmn_typ,(6,psps%lmnmax,1))
   ABI_ALLOCATE(e1kbfr,(dime1kb,gs_hamkq%dimekb2,nspinor**2))
   if (dime1kb>0) e1kbfr=zero

!  The following contributions are needed only for non-DDK perturbation:
!  - Frozen part of 1st-order Dij
!  - Contribution from local potential to dynamical matrix (due to Vxc^(j1)(tild_nc)+VH^(j1)(tild_nZc))
   if (ipert/=dtset%natom+1.and.ipert1/=dtset%natom+1) then

!    Allocations
     nfftot=ngfftf(1)*ngfftf(2)*ngfftf(3)
     nvxc1=0;if (ipert1<=dtset%natom) nvxc1=nspden
     ABI_ALLOCATE(vxc10,(cplex*nfftf,nvxc1))
     ABI_ALLOCATE(paw_ij10,(dtset%natom,mdir1))

!    LOOP OVER PERTURBATION DIRECTIONS
     do kdir1=1,mdir1
       idir1=jdir1(kdir1)

!      Get first-order local potential and first-order pseudo core density
       if (ipert==ipert1.and.idir==idir1) then
         vpsp1_idir1 => vpsp1
         xccc3d1_idir1 => xccc3d1
       else
         ABI_ALLOCATE(vpsp1_idir1,(cplex*nfftf))
         ABI_ALLOCATE(xccc3d1_idir1,(cplex*n3xccc))
         optv=1;optn=n3xccc/nfftf;optn2=1
         call atm2fft3(gs_hamkq%atindx,xccc3d1_idir1,vpsp1_idir1,cplex,dum1,gmet,gsqcut,idir1,ipert1,&
&         mgfftf,mpi_enreg,psps%mqgrid_vl,dtset%natom,1,nfftf,ngfftf,dtset%ntypat,&
&         optn,optn2,optv,dtset%paral_kgb,pawtab,ph1df,psps%qgrid_vl,dtset%qptn,&
&         dtset%typat,ucvol,psps%usepaw,psps%vlspl,xred)
       end if

!      Compute 1st-order non-local factors (Dij^(j1)_fr)
       call nullify_paw_ij(paw_ij10(:,idir1))
       call init_paw_ij(paw_ij10(:,idir1),cplex,paw_ij1(1)%cplex_dij,nspinor,dtset%nsppol,dtset%nspden,&
&       0,dtset%natom,dtset%ntypat,dtset%typat,pawtab,has_dijfr=1)
       if (ipert/=ipert1.or.idir/=idir1) then
         optfr=1
         if (usexcnhat==0) then
           ABI_ALLOCATE(vtmp,(nfftf,nspden))
           vtmp=vtrial-vxc
         else
           vtmp => vtrial
         end if
         call pawfrnhat(cplex,gprimd,idir1,ipert1,mpi_enreg,dtset%natom,dtset%natom,nfftf,ngfftf,&
&         nspden,dtset%ntypat,optfr,paw_ij10(:,idir1),pawang,pawfgrtab,pawrad,pawrhoij,pawtab,&
&         dtset%qptn,rprimd,ucvol,vpsp1_idir1,vtmp,xred)
         if (usexcnhat==0)  then
           ABI_DEALLOCATE(vtmp)
         end if
       else
         do iatom=1,dtset%natom
           paw_ij10(iatom,idir1)%has_dijfr=paw_ij1(iatom)%has_dijfr
           if (paw_ij1(iatom)%has_dijfr==2) paw_ij10(iatom,idir1)%dijfr=paw_ij1(iatom)%dijfr
         end do
       end if

!      Get first-order exchange-correlation potential (core-correction contribution only)
       if (nvxc1>0) then
         if (n3xccc/=0)then
           call mkvxc3(cplex,kxc,mpi_enreg,nfftf,ngfftf,nkxc,nspden,n3xccc,&
&           0,dtset%paral_kgb,dtset%qptn,dum1,rprimd,vxc10,xccc3d1_idir1)
         else
           vxc10=zero
         end if
       end if

!      Get local contribution to dynamical matrix
       if (ipert1<=dtset%natom) then
         if (usexcnhat/=0) then
!          vxc1 is integrated with the total 1st-order density (rhor1 including nhat1)
!          vpsp1 is integrated with the 1st-order pseudo density (rhor1 without nhat1)
           ABI_ALLOCATE(rhotmp,(cplex*nfftf,1))
           rhotmp(:,1)=rhor1(:,1)-nhat1(:,1)
           call dotprod_vn(cplex,rhor1,dot1r,dot1i,mpi_enreg,nfftf,nfftot,nspden,2,vxc10,ucvol)
           call dotprod_vn(cplex,rhotmp,dot2r,dot2i,mpi_enreg,nfftf,nfftot,1,2,vpsp1_idir1,ucvol)
         else
!          vxc1 is integrated with the 1st-order pseudo density (rhor1 without nhat1)
!          vpsp1 is integrated with the 1st-order pseudo density (rhor1 without nhat1)
           ABI_ALLOCATE(rhotmp,(cplex*nfftf,nspden))
           rhotmp(:,:)=rhor1(:,:)-nhat1(:,:)
           call dotprod_vn(cplex,rhotmp,dot1r,dot1i,mpi_enreg,nfftf,nfftot,nspden,2,vxc10,ucvol)
           call dotprod_vn(cplex,rhotmp,dot2r,dot2i,mpi_enreg,nfftf,nfftot,1,2,vpsp1_idir1,ucvol)
         end if
         ABI_DEALLOCATE(rhotmp)
!        Note: factor 2 (from d2E/dj1dj2=2E^(j1j2)) eliminated by factor 1/2
         dotr=dot1r+dot2r;doti=dot1i+dot2i
!        In case ipert = natom+2, these lines compute the local part
!        of the Born effective charges from phonon and electric
!        field type perturbations, see eq. 43 of X. Gonze and C. Lee, PRB 55, 10355 (1997)
!        The minus sign is due to the fact that the effective charges
!        are minus the second derivatives of the energy
!        This has to be checked for ipert=natom+5 (magnetic field pert.)
         if (ipert==dtset%natom+2) then
           d2lo(1,idir1,ipert1,idir,ipert)=-dotr
           d2lo(2,idir1,ipert1,idir,ipert)=-doti
         else if (ipert/=dtset%natom+1) then
           d2lo(1,idir1,ipert1,idir,ipert)=dotr
           d2lo(2,idir1,ipert1,idir,ipert)=doti
         end if
       end if ! ipert1<=natom

       if (ipert/=ipert1.or.idir/=idir1)  then
         ABI_DEALLOCATE(vpsp1_idir1)
         ABI_DEALLOCATE(xccc3d1_idir1)
       end if

!      End loop on directions
     end do

!    Free memory
     ABI_DEALLOCATE(vxc10)

   else ! ddk perturbation
     d2lo(1:2,1:mdir1,ipert1,idir,ipert)=zero
   end if

!  Prepare GS k wf file for reading if mkmem==0
   if (dtset%mkmem==0.and.need_wfk) then
     formeig=0
     call clsopn(wfftgs)
     call hdr_skip(wfftgs,ierr)
     call WffKg(wfftgs,1)
     call xdefineOff(formeig,wfftgs,mpi_enreg,nband_rbz,npwarr,nspinor,nsppol,nkpt_rbz)
   end if

!  Prepare GS k+q wf file for reading if mkqmem==0
   if (dtset%mkqmem==0.and.need_wfq) then
     formeig=0
     call clsopn(wfftkq)
     call hdr_skip(wfftkq,ierr)
     call WffKg(wfftkq,1)
     call xdefineOff(formeig,wfftkq,mpi_enreg,nband_rbz,npwar1,nspinor,nsppol,nkpt_rbz)
   end if

!  Prepare RF wf files for reading and writing if mk1mem==0
   if (dtset%mk1mem==0.and.need_wf1) then
     formeig=1
     call clsopn(wffnow)
     call hdr_skip(wffnow,ierr)
     call WffKg(wffnow,1)
     call xdefineOff(formeig,wffnow,mpi_enreg,nband_rbz,npwar1,nspinor,nsppol,nkpt_rbz)
   end if

!  Prepare RF PAW files for reading and writing if mkmem, mkqmem or mk1mem==0
   iorder_cprj=0;iorder_cprj1=0
   if (need_wfk) then
     call cprj_diskinit_r(gs_hamkq%atindx1,dtset%natom,iorder_cprj,dtset%mkmem,&
&     dtset%natom,ncpgr,dimcprj,nspinor,dtfil%unpaw)
   end if
   if (need_wfq) then
     call cprj_diskinit_r(gs_hamkq%atindx1,dtset%natom,iorder_cprj,dtset%mkqmem,&
&     dtset%natom,0,dimcprj,nspinor,dtfil%unpawq)
   end if
   if (need_wf1) then
     call cprj_diskinit_r(gs_hamkq%atindx1,dtset%natom,iorder_cprj1,dtset%mk1mem,&
&     dtset%natom,0,dimcprj,nspinor,dtfil%unpaw1)
   end if

!  Prepare DDK files for reading
   skipddk(:)=0;ikpt_fbz(1:3)=0;ikpt_fbz_previous(1:3)=0
   if (need_ddk_file) then
     do kdir1=1,mdir1
       idir1=jdir1(kdir1)
       if (ddkfil(idir1)/=0) then
         call clsopn(wffddk(idir1))
         call hdr_skip(wffddk(idir1),ierr)
       end if
     end do
   end if

!  Allocate arrays used to accumulate density change due to overlap
   if (has_drho) then
     ABI_ALLOCATE(drhoaug1,(cplex*dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),mdir1))
     ABI_ALLOCATE(drho1wfr,(cplex*dtset%nfft,dtset%nspden,mdir1))
     ABI_ALLOCATE(pawdrhoij1,(dtset%natom,mdir1))
     do kdir1=1,mdir1
       idir1=jdir1(kdir1)
       drho1wfr(:,:,idir1)=zero
       ABI_ALLOCATE(nlmn,(dtset%ntypat))
       do itypat=1,dtset%ntypat
         nlmn(itypat)=pawtab(itypat)%lmn_size
       end do
       call rhoij_alloc(pawrhoij1(1)%cplex,nlmn,pawrhoij1(1)%nspden,pawrhoij1(1)%nspinor,&
&       pawrhoij1(1)%nsppol,pawdrhoij1(:,idir1),dtset%typat,mpi_enreg=mpi_enreg,use_rhoij_=1)
       ABI_DEALLOCATE(nlmn)
     end do
   end if

!  Initialize shifts for global arrays
   bdtot_index=0
   bd2tot_index=0
   ibg=0;icg=0
   ibg1=0;icg1=0
   ibgq=0;icgq=0
   ikpt_me=0

!  LOOP OVER SPINS
   do isppol=1,nsppol

!    Rewind (k+G) data if needed
     ikg=0;ikg1=0;ikpt_fbz(1:3)=0
     if (dtset%mkmem==0) rewind(dtfil%unkg)
     if (dtset%mkmem==0.and.psps%useylm==1) rewind(dtfil%unylm)
     if (dtset%mk1mem==0) rewind(dtfil%unkg1)
     if (dtset%mk1mem==0.and.psps%useylm==1) rewind(dtfil%unylm1)

!    DDK files: skip the remaining isppol=1 records
     if (isppol==2.and.need_ddk_file) then
       do kdir1=1,mdir1
         idir1=jdir1(kdir1)
         if ((ddkfil(idir1)/=0).and.(skipddk(idir1)<nkpt)) then
           do ikpt=1,(nkpt-skipddk(idir1))
             call WffReadNpwRec(ierr,ikpt,isppol,nband_,npw_,nspinor_,wffddk(idir1))
             call WffReadSkipRec(ierr,1+2*nband_,wffddk(idir1))
           end do
         end if
       end do
     end if

!    Retrieve factors for 0-order and 1st-order NL Hamiltonian (PAW Dij coefficients) for this spin
!    1-Retrieve factors for 0-order NL Hamiltonian
     do ispden=1,nspinor**2
       isp=isppol;if (nspinor==2) isp=ispden
       do iatom=1,dtset%natom
         dimdij=paw_ij(iatom)%cplex*paw_ij(iatom)%lmn2_size
         do ilmn=1,dimdij
           gs_hamkq%ekb(ilmn,iatom,ispden)=paw_ij(iatom)%dij(ilmn,isp)
         end do
         if(dimdij+1<=gs_hamkq%dimekb1) gs_hamkq%ekb(dimdij+1:gs_hamkq%dimekb1,iatom,ispden)=zero
       end do
     end do
!    2-Retrieve factors for 1st-order NL Hamiltonian
     if (ipert1<=dtset%natom) then
       itypat=gs_hamkq%typat(ipert1)
       dimdij=paw_ij(ipert1)%cplex_dij*paw_ij(ipert1)%lmn2_size
       do ispden=1,nspinor**2
         isp=isppol;if (nspinor==2) isp=ispden
         ekb_typ(1:dimdij,1,ispden)=paw_ij(ipert1)%dij(1:dimdij,isp)
         if (paw_ij(ipert1)%cplex_dij==1) then
           sij_typ(1:dimdij,1)=gs_hamkq%sij(1:dimdij,itypat)
         else
           do ilmn=1,dimdij/2
             sij_typ(2*ilmn-1,1)=gs_hamkq%sij(ilmn,itypat)
             sij_typ(2*ilmn  ,1)=zero
           end do
         end if
         if (dimdij<gs_hamkq%dimekb1) then
           ekb_typ(dimdij+1:gs_hamkq%dimekb1,1,ispden)=zero
           sij_typ(dimdij+1:gs_hamkq%dimekb1,1)=zero
         end if
       end do
       indlmn_typ(:,:,1)=psps%indlmn(:,:,itypat)
       pspso_typ(1)=psps%pspso(itypat)
     else
       ekb_typ=zero;sij_typ=zero;pspso_typ=1
     end if

!    Initialize accumulation of density
     if (has_drho) drhoaug1(:,:,:,:)=zero

!    LOOP OVER K-POINTS
     do ikpt=1,nkpt_rbz

!      Load dimensions for this k-point
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       istwf_k=istwfk_rbz(ikpt)
       npw_k=npwarr(ikpt)
       npw1_k=npwar1(ikpt)

!      Skip loop if this k-point is not to be treated by this proc
       if(mpi_enreg%paral_compil_kpt==1)then
         if (minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))/=0) then
           bdtot_index=bdtot_index+nband_k
           bd2tot_index=bd2tot_index+2*nband_k**2
           cycle ! Skip the rest of the k-point loop
         end if
       end if

!      Allocate/initialize local arrays and scalars for this k-point
       ABI_ALLOCATE(d2nl_k,(2,3))
       ABI_ALLOCATE(d2ovl_k,(2,3))
       ABI_ALLOCATE(eig_k,(nband_k))
       ABI_ALLOCATE(eig_kq,(nband_k))
       ABI_ALLOCATE(eig1_k,(2*nband_k**2))
       ABI_ALLOCATE(occ_k,(nband_k))
       d2nl_k(:,:)=zero;d2ovl_k(:,:)=zero
       eig_k (:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
       eig_kq(:)=eigenq(1+bdtot_index:nband_k+bdtot_index)
       eig1_k(:)=eigen1(1+bd2tot_index:2*nband_k**2+bd2tot_index)
       occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)
       nband_kocc=count(abs(occ_k(:))>tol8)
       kpoint(:)=kpt_rbz(:,ikpt)
       wtk_k=wtk_rbz(ikpt)
       need_ylmgr1=0;dimylmgr1=0
       nkpg=0;nkpg1=0;dimphkxred=0
       ikpt_me=ikpt_me+1
       if (ipert1<=dtset%natom) then
         dimphkxred=dtset%natom
         ABI_ALLOCATE(phkxred,(2,dimphkxred))
         do ia=1,dtset%natom
           iatom=gs_hamkq%atindx(ia)
           arg=two_pi*(kpoint(1)*gs_hamkq%xred(1,ia)+kpoint(2)*gs_hamkq%xred(2,ia) &
&           +kpoint(3)*gs_hamkq%xred(3,ia))
           phkxred(1,iatom)=cos(arg) ; phkxred(2,iatom)=sin(arg)
         end do
       else
         ABI_ALLOCATE(phkxred,(2,dimphkxred))
       end if
       if (is_metal) then
!        For each pair of active bands (m,n), generates the ratios
!        rocceig(m,n)=(occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n))
         ABI_ALLOCATE(doccde_k,(nband_k))
         ABI_ALLOCATE(doccde_kq,(nband_k))
         ABI_ALLOCATE(occ_kq,(nband_k))
         ABI_ALLOCATE(rocceig,(nband_k,nband_k))
         doccde_k(:)=doccde_rbz(1+bdtot_index:nband_k+bdtot_index)
         doccde_kq(:)=docckqde(1+bdtot_index:nband_k+bdtot_index)
         occ_kq(:)=occkq(1+bdtot_index:nband_k+bdtot_index)
         call occeig(doccde_k,doccde_kq,eig_k,eig_kq,nband_k,dtset%occopt,occ_k,occ_kq,rocceig)
       end if

!      Continue to initialize the Hamiltonian at k+q
       gs_hamkq%istwf_k=istwf_k
       gs_hamkq%npw    =npw1_k
       gs_hamkq%kpoint(:)=kpoint(:)+dtset%qptn(:)
       do ia=1,dtset%natom
         iatom=gs_hamkq%atindx(ia)
         arg=two_pi*(gs_hamkq%kpoint(1)*xred(1,ia)+gs_hamkq%kpoint(2)*xred(2,ia)&
&         +gs_hamkq%kpoint(3)*xred(3,ia))
         gs_hamkq%phkxred(1,iatom)=cos(arg);gs_hamkq%phkxred(2,iatom)=sin(arg)
       end do

!      Take care of the npw and kg records in WF and DDK files
       if (dtset%mkmem==0.and.need_wfk) then
         call WffReadNpwRec(ierr,ikpt,isppol,nband_,npw_,nspinor_,wfftgs)
         call WffReadSkipRec(ierr,2,wfftgs)
       end if
       if (dtset%mkqmem==0.and.need_wfq) then
         call WffReadNpwRec(ierr,ikpt,isppol,nband_,npw_,nspinor_,wfftkq)
         call WffReadSkipRec(ierr,2,wfftkq)
       end if
       if (dtset%mk1mem==0.and.need_wf1) then
         call WffReadNpwRec(ierr,ikpt,isppol,nband_,npw_,nspinor_,wffnow)
         call WffReadSkipRec(ierr,1,wffnow)
       end if
       if (need_ddk_file) then
         do kdir1=1,mdir1
           idir1=jdir1(kdir1)
           if (ddkfil(idir1)/=0)then
!            Skip records in DDK file
             ikpt_fbz_previous(idir1)=ikpt_fbz(idir1)
             ikpt_fbz(idir1)=indkpt1(ikpt)
!            Number of k points to skip in the full set of k pointsp
             nskip=ikpt_fbz(idir1)-ikpt_fbz_previous(idir1)-1
             skipddk(idir1)=skipddk(idir1)+nskip+1
             if (nskip/=0) then
               do ikpt_=1+ikpt_fbz_previous(idir1),ikpt_fbz(idir1)-1
                 call WffReadSkipK(1,0,ikpt_,isppol,mpi_enreg,wffddk(idir1))
               end do
             end if
!            Begin to read current record (k+G)
             call WffReadNpwRec(ierr,ikpt,isppol,nband_,npw_,nspinor_,wffddk(idir1))
             if (npw_/=npw_k) then
               write(unit=msg,fmt='(a,i3,a,i5,a,i3,a,a,i5,a,a,i5)')&
&               ' For isppol = ',isppol,', ikpt = ',ikpt,' and idir = ',idir,ch10,&
&               ' the number of plane waves in the ddk file is equal to', npw_,ch10,&
&               ' while it should be ',npw_k
               MSG_ERROR(msg)
             end if
             call WffReadSkipRec(ierr,1,wffddk(idir1))
           end if
         end do
       end if

!      Eventually load WF at k+q (needed for the computation of delta_u^(j1))
       mcgq=0;mcprjq=0
       if (has_dcwf) then
         if (dtset%mkqmem==0) then
           mcgq=npw1_k*nspinor*nband_k
           ABI_ALLOCATE(cgq_disk,(2,mcgq))
           do iband=1,nband_k
             call WffReadDataRec(cgq_disk(:,(iband-1)*npw1_k*nspinor+1:iband*npw1_k*nspinor),&
&             ierr,2,npw1_k*nspinor,wfftkq)
           end do
           mcprjq=nspinor*nband_k*usecprj
           ABI_ALLOCATE(cprjq_disk,(dtset%natom,mcprjq))
           if (mcprjq>0) then
             call cprj_alloc(cprjq_disk,0,dimcprj)
             call cprj_get(gs_hamkq%atindx1,cprjq_disk,cprjq,dtset%natom,1,ibgq,ikpt,iorder_cprj,isppol,&
&             dtset%mband,dtset%mkqmem,mpi_enreg,dtset%natom,nband_k,nband_k,nspinor,nsppol,dtfil%unpawq)
           end if
         else
           mcgq=mpw1*nspinor*dtset%mband*dtset%mkqmem*nsppol
           mcprjq=nspinor*dtset%mband*dtset%mkqmem*nsppol*usecprj
         end if
       end if

!      Allocate arrays used for NL form factors
       ABI_ALLOCATE(kg_k,(3,npw_k))
       ABI_ALLOCATE(kg1_k,(3,npw1_k))
       ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
       ABI_ALLOCATE(ylm1_k,(npw1_k,psps%mpsang*psps%mpsang*psps%useylm))
       if (need_ddk_file.or.ipert1==dtset%natom+1) need_ylmgr1=1
       dimylmgr1=max(useylmgr1,need_ylmgr1)
       ABI_ALLOCATE(ylmgr1_k,(npw1_k,3,psps%mpsang*psps%mpsang*psps%useylm*dimylmgr1))
       ABI_ALLOCATE(gbound,(2*dtset%mgfft+8,2))

!      Read plane-wave vectors and related data at k
       if (dtset%mkmem==0) then
         nspinor_=nspinor
         call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor_,0,dtfil%unkg)
         read(dtfil%unkg) ((kg_k(ii,jj),ii=1,3),jj=1,npw_k)
         if (psps%useylm==1) then
           read(dtfil%unylm)
           read(dtfil%unylm) ((ylm_k(ii,jj),ii=1,npw_k),jj=1,psps%mpsang*psps%mpsang)
         end if
       else
         kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
         if (psps%useylm==1) then
           do jj=1,psps%mpsang*psps%mpsang
             ylm_k(1:npw_k,jj)=ylm(1+ikg:npw_k+ikg,jj)
           end do
         end if
       end if
       call sphereboundary(gbound,istwf_k,kg_k,dtset%mgfft,npw_k)

!      Read plane-wave vectors and related data at k+q
       if (dtset%mk1mem==0) then
         nspinor_=nspinor
         call rdnpw(ikpt,isppol,nband_k,npw1_k,nspinor_,0,dtfil%unkg1)
         read(dtfil%unkg1) ((kg1_k(ii,jj),ii=1,3),jj=1,npw1_k)
         if (psps%useylm==1) then
           read(dtfil%unylm1)
           if (need_ylmgr1==1.and.useylmgr1/=0) then
             read(dtfil%unylm1) ((ylm1_k(ii,jj),ii=1,npw1_k),jj=1,psps%mpsang*psps%mpsang),&
&             (((ylmgr1_k(ii,ia,jj),ii=1,npw1_k),ia=1,3),jj=1,psps%mpsang*psps%mpsang)
           else
             read(dtfil%unylm1) ((ylm1_k(ii,jj),ii=1,npw1_k),jj=1,psps%mpsang*psps%mpsang)
           end if
         end if
       else
         kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
         if (psps%useylm==1) then
           do jj=1,psps%mpsang*psps%mpsang
             ylm1_k(1:npw1_k,jj)=ylm1(1+ikg1:npw1_k+ikg1,jj)
           end do
           if (need_ylmgr1==1.and.useylmgr1/=0) then
             do jj=1,psps%mpsang*psps%mpsang
               do ia=1,3
                 ylmgr1_k(1:npw1_k,ia,jj)=ylmgr1(1+ikg1:npw1_k+ikg1,ia,jj)
               end do
             end do
           end if
         end if
       end if
       call sphereboundary(gs_hamkq%gbound,istwf_k,kg1_k,dtset%mgfft,npw1_k)

!      If Ylm gradients at k+q are needed and not in memory, compute them
       if (need_ylmgr1==1.and.useylmgr1==0) then
         option=-1;npwar1_tmp(1)=npw1_k;nband_tmp(1)=nband_k
         call initylmg(gprimd,kg1_k,gs_hamkq%kpoint,1,mpi_enreg,psps%mpsang,&
&         npw1_k,nband_tmp,1,npwar1_tmp,nsppol,option,rprimd,&
&         tmp_unit,tmp_unit,ylm1_k,ylmgr1_k)
       end if

!      Compute (k+G) vectors
       nkpg=0;if(ipert1<=dtset%natom) nkpg=3*dtset%nloalg(5)
       ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
       if (nkpg>0) then
         call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
       end if

!      Compute (k+q+G) vectors
       nkpg=0;if(ipert1<=dtset%natom) nkpg=3*dtset%nloalg(5)
       ABI_ALLOCATE(kpg1_k,(npw1_k,nkpg1))
       if (nkpg1>0) then
         call mkkpg(kg1_k,kpg1_k,gs_hamkq%kpoint,nkpg1,npw1_k)
       end if

!      Compute nonlocal form factors ffnl at (k+G), for all atoms
       dimffnl=0;if (ipert1/=dtset%natom+1) dimffnl=1
       ABI_ALLOCATE(ffnlk,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
       ider=0;idir0=0
       if (ipert1<=dtset%natom) then
         call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnlk,psps%ffspl,gs_hamkq%gmet,gs_hamkq%gprimd,&
&         ider,idir0,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,psps%lnmax,psps%mpsang,&
&         psps%mqgrid_ff,nkpg,npw_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,&
&         psps%useylm,ylm_k,ylmgr_dum)
       end if

!      Compute nonlocal form factors ffnl1 at (k+q+G), for all atoms
       dimffnl1=1+3*need_ylmgr1;ider=need_ylmgr1
       idir0=0;if (ider>0) idir0=4  ! This will be -7 for strain perturbation
       ABI_ALLOCATE(ffnl1,(npw1_k,dimffnl1,psps%lmnmax,psps%ntypat))
       call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1,psps%ffspl,gs_hamkq%gmet,gs_hamkq%gprimd,&
&       ider,idir0,psps%indlmn,kg1_k,kpg1_k,gs_hamkq%kpoint,psps%lmnmax,psps%lnmax,psps%mpsang,&
&       psps%mqgrid_ff,nkpg1,npw1_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,&
&       psps%useylm,ylm1_k,ylmgr1_k)

!      Extract non-local form factor for H^(j1)
       if (ipert1<=dtset%natom) then
         ffnlkq => ffnl1(:,:,:,dtset%typat(ipert1))
       else
         ffnlkq => ffnl1(:,:,:,1)  !  unused in that case
       end if
       if (ipert1<=dtset%natom) then
         ffnl1_idir1 => ffnl1(:,:,:,:)
         dimffnl1_idir1=dimffnl1
       else
         dimffnl1_idir1=1+need_ylmgr1
         ABI_ALLOCATE(ffnl1_idir1,(npw1_k,dimffnl1_idir1,psps%lmnmax,psps%ntypat))
         do itypat=1,psps%ntypat
           do ilmn=1,psps%lmnmax
             ffnl1_idir1(1:npw1_k,1,ilmn,itypat)=ffnl1(1:npw1_k,1,ilmn,itypat)
           end do
         end do
       end if
       if (ipert1==dtset%natom+2) then
         do itypat=1,psps%ntypat
           do ilmn=1,psps%lmnmax
             do ii=1,min(dimffnl,dimffnl1)
               ffnlk(1:npw1_k,ii,ilmn,itypat)=ffnl1(1:npw1_k,ii,ilmn,itypat)
             end do
           end do
         end do
       else if (ipert1==dtset%natom+5) then
         ffnlk=zero  ! to be activated later
       end if

!      Allocate (and eventually compute) ph3d at k+q+G
       if (gs_hamkq%nloalg(1)<=0) then
         matblk=gs_hamkq%nloalg(4)
         ABI_ALLOCATE(ph3d,(2,npw1_k,matblk))
       else
         matblk=dtset%natom
         ABI_ALLOCATE(ph3d,(2,npw1_k,matblk))
         call ph1d3d(1,dtset%natom,kg1_k,matblk,dtset%natom,npw1_k,&
&         dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),&
&         gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d)
       end if

!      Allocate memory space for one band
       if (need_wfk)  then
         ABI_ALLOCATE(cwave0,(2,npw_k*nspinor))
       end if
       if (need_wf1)  then
         ABI_ALLOCATE(cwavef,(2,npw1_k*nspinor))
       end if
       ABI_ALLOCATE(gh1,(2,npw1_k*nspinor))
       nullify(cwaveprj0_idir1)
       if (usecprj==1) then
         ABI_ALLOCATE(cwaveprj0,(dtset%natom,nspinor))
         call cprj_alloc(cwaveprj0,ncpgr,dimcprj)
         if (ncpgr>1) then
           ABI_ALLOCATE(cwaveprj0_idir1,(dtset%natom,nspinor))
           call cprj_alloc(cwaveprj0_idir1,1,dimcprj)
         end if
       end if
       if (has_dcwf) then
         ABI_ALLOCATE(gs1,(2,npw1_k*nspinor))
       else
         ABI_ALLOCATE(gs1,(0,0))
       end if

!      Allocate unused kinetic contributions
       ABI_ALLOCATE(dkinpw,(npw_k))
       ABI_ALLOCATE(kinpw1,(npw1_k))
       kinpw1=zero;dkinpw=zero

!      LOOP OVER BANDS
       do iband=1,nband_k

!        Skip band if not to be treated by this proc
         if (mpi_enreg%paral_compil_kpt==1) then
           if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me) then
             if (dtset%mkmem==0.and.need_wfk) then
               call WffReadSkipRec(ierr,1,wfftgs)
             end if
             if (dtset%mk1mem==0.and.need_wf1) then
               call WffReadSkipRec(ierr,2,wffnow)
             end if
             if (need_ddk_file) then
               do kdir1=1,mdir1
                 idir1=jdir1(kdir1)
                 if (ddkfil(idir1)/=0) then
                   call WffReadSkipRec(ierr,2,wffddk(idir1))
                 end if
               end do
             end if
             cycle
           end if
         end if

!        Extract GS wavefunctions
         if (need_wfk) then
           if (dtset%mkmem/=0) then
             cwave0(:,:)=cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
           else
             call WffReadDataRec(cwave0,ierr,2,npw_k*nspinor,wfftgs)
           end if
           if (usecprj==1) then
             call cprj_get(gs_hamkq%atindx1,cwaveprj0,cprj,dtset%natom,iband,ibg,ikpt,iorder_cprj,&
&             isppol,dtset%mband,dtset%mkmem,mpi_enreg,dtset%natom,1,nband_k,nspinor,nsppol,dtfil%unpaw)
           end if
         end if

!        Extract 1st-order wavefunctions
         if (need_wf1) then
           if (dtset%mk1mem/=0) then
             cwavef(:,:)=cg1(:,1+(iband-1)*npw1_k*nspinor+icg1:iband*npw1_k*nspinor+icg1)
           else
             call WffReadSkipRec(ierr,1,wffnow)
             call WffReadDataRec(cwavef,ierr,2,npw1_k*nspinor,wffnow)
           end if
         end if

!        LOOP OVER PERTURBATION DIRECTIONS
         do kdir1=1,mdir1
           idir1=jdir1(kdir1)

!          Not able to compute if ipert1=(Elect. field) and no ddk WF file
           if ((ipert1==dtset%natom+2.or.ipert1==dtset%natom+5).and.ddkfil(idir1)==0) cycle

!          Extract ground state projected WF and derivatives in idir1 direction
           if ((ipert1<=dtset%natom+1.or.ipert1==dtset%natom+2.or.ipert1==dtset%natom+5) &
&           .and.usecprj==1.and.ncpgr>1) then
             call cprj_copy(cwaveprj0,cwaveprj0_idir1,icpgr=idir1)
           else
             cwaveprj0_idir1 => cwaveprj0
           end if

!          Extract 1st-order NL form factors for this idir1
           if (dimffnl1_idir1>=2.and.ipert1>dtset%natom) then
             do itypat=1,psps%ntypat
               do ilmn=1,psps%lmnmax
                 ffnl1_idir1(1:npw1_k,2,ilmn,itypat)=ffnl1(1:npw1_k,1+idir1,ilmn,itypat)
               end do
             end do
           end if

!          Eventually compute 1st-order kinetic operator
           if (ipert1==dtset%natom+1) then
             call kpg3(dkinpw,dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,idir1,kg_k,kpoint,npw_k)
           end if

!          Extract frozen part of 1st-order Dij
           if (ipert/=dtset%natom+1.and.ipert1/=dtset%natom+1) then
             do ispden=1,nspinor**2
               isp=isppol;if (nspinor==2) isp=ispden
               do iatom=1,dtset%natom
                 dimdij1=paw_ij10(iatom,idir1)%cplex_dij*paw_ij10(iatom,idir1)%lmn2_size
                 e1kbfr(1:dimdij1,iatom,ispden)=paw_ij10(iatom,idir1)%dijfr(1:dimdij1,isp)
               end do
             end do
           end if

!          Read DDK wave function (if ipert1=electric field)
           if (ipert1==dtset%natom+2.or.ipert1==dtset%natom+5) then
             usevnl=1
             ABI_ALLOCATE(gvnl1,(2,npw1_k*nspinor*usevnl))
             if (need_ddk_file) then
               if (ddkfil(idir1)/=0) then
                 call WffReadSkipRec(ierr,1,wffddk(idir1))
                 call WffReadDataRec(gvnl1,ierr,2,npw1_k*nspinor,wffddk(idir1))
               else
                 gvnl1=zero
               end if
               if (ipert1==dtset%natom+2) then
                 do ii=1,npw1_k*nspinor ! Multiply ddk by +i (to be consistent with getgh1c)
                   arg=gvnl1(1,ii)
                   gvnl1(1,ii)=-gvnl1(2,ii)
                   gvnl1(2,ii)=arg
                 end do
               end if
             else
               gvnl1=zero
             end if
           else
             usevnl=0
             ABI_ALLOCATE(gvnl1,(2,npw1_k*nspinor*usevnl))
           end if

!          Get |H^(j2)-Eps_k_i.S^(j2)|u0_k_i> (VHxc-dependent part not taken into account) and S^(j2)|u0>
           lambda=eig_k(iband);berryopt=1;cplx=1;optlocal=0
           optnl=0;if (ipert1/=dtset%natom+1.or.idir==idir1) optnl=1
           opt_gvnl1=0;if (ipert1==dtset%natom+2) opt_gvnl1=2
           sij_opt=-1;if (has_dcwf) sij_opt=1
           call getgh1c(berryopt,cplx,cwave0,cwaveprj0_idir1,&
&           dimcprj,gs_hamkq%dimekb1,dime1kb,dimffnl,dimffnl,dimffnl1_idir1,dimphkxred,&
&           dkinpw,ekb_typ,e1kbfr,e1kbfr,ffnlk,ffnlkq,ffnl1_idir1,dtfil%filstat,&
&           gs_hamkq%gbound,gh1,dum1,gs1,gs_hamkq,gvnl1,idir1,indlmn_typ,&
&           ipert1,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpoint,lambda,psps%lmnmax,matblk,dtset%mgfft,&
&           mpi_enreg,psps%mpsang,psps%mpssoang,dtset%natom,nkpg,nkpg1,npw_k,npw1_k,nspinor,&
&           gs_hamkq%ntypat,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,optlocal,optnl,opt_gvnl1,&
&           dtset%paral_kgb,ph3d,phkxred,dtset%prtvol,sij_opt,sij_typ,tim_getgh1c,usecprj,&
&           usee1kb,usevnl,dum2,dum3)
           if (sij_opt==1.and.optnl==1) gh1=gh1-lambda*gs1
           ABI_DEALLOCATE(gvnl1)

!          If needed, compute here <delta_u^(j1)_k_i|H^(j2)-Eps_k_i.S^(j2)|u0_k_i>
!          with delta_u^(j1)=-1/2 Sum_{j}[<u0_k+q_j|S^(j1)|u0_k_i>.|u0_k+q_j>]
!          (see PRB 78, 035105 (2008), Eq. (42))
!          This can be rewritten as:
!          -1/2.<u0_k_i|S^(j1)| Sum_{j}[<u0_k+q_j|H^(j2)-Eps_k_i.S^(j2)|u0_k_i>.|u0_k+q_j>
!          The sum over j can be computed with a single call to projbd routine
!          At first call (when j1=j2), ch1c=<u0_k+q_j|H^(j2)-Eps_k_i.S^(j2)|u0_k_i> is stored
!          For the next calls, it is re-used.
           if (has_dcwf) then
!            note: gvnl1 used as temporary space
             ABI_ALLOCATE(gvnl1,(2,npw1_k*nspinor))
             if (ipert==ipert1.and.idir==idir1) then
               option=0;gvnl1=gh1
             else
               option=1;gvnl1=zero
             end if
!            Compute -Sum_{j}[<u0_k+q_j|H^(j2)-Eps_k_i.S^(j2)|u0_k_i>.|u0_k+q_j>
             if (dtset%mkqmem==0) then
               call projbd(cgq_disk,gvnl1,-1,icgq,0,istwf_k,mcgq,mpi_enreg,0,nband_k,npw1_k,nspinor,&
&               dtset%ortalg,0,dum1,ch1c(:,1:nband_k,iband,ikpt_me),option,tim_projbd,0)
             else
               call projbd(cgq,gvnl1,-1,icgq,0,istwf_k,mcgq,mpi_enreg,0,nband_k,npw1_k,nspinor,&
&               dtset%ortalg,0,dum1,ch1c(:,1:nband_k,iband,ikpt_me),option,tim_projbd,0)
             end if
             if (ipert==ipert1.and.idir==idir1) gvnl1=gvnl1-gh1
             if (abs(occ_k(iband))>tol8) then
!              Compute: -<u0_k_i|S^(j1)| Sum_{j}[<u0_k+q_j|H^(j2)-Eps_k_i.S^(j2)|u0_k_i>.|u0_k+q_j>
               call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1_k*nspinor,2,gs1,gvnl1)
!              Add contribution to DDB
!              Note: factor 2 (from d2E/dj1dj2=2E^(j1j2)) eliminated by factor 1/2
!              (-1) factor already present
               d2ovl_k(1,idir1)=d2ovl_k(1,idir1)+wtk_k*occ_k(iband)*dotr
               d2ovl_k(2,idir1)=d2ovl_k(2,idir1)+wtk_k*occ_k(iband)*doti
             end if
             ABI_DEALLOCATE(gvnl1)
           end if

!          If needed, compute here <delta_u^(j1)_k_i|H-Eps_k_i.S|u^(j2)_k_i>
!          This is equal to <delta_u^(j1)_k_i|H-Eps_k_i.S|delta_u^(j2)_k_i>  (I)
!          +<delta_u^(j1)_k_i|H-Eps_k_i.S|u^paral^(j2)_k_i>  (II)
!          (u^paral^(j2)_k_i is the part of u^(j2)_k_i parallel to active space : metals)
!          (I) can be rewritten as:
!          Sum_j{ 1/4.<u0_k_i|S^(j1)|u0_k+q_j>.<u0_k+q_j|S^(j2)|u0_k_i>.(Eps_k+q_j-Eps_k_i) }
!          (II) can be rewritten as:
!          Sum_j{1/2.(occ_kq_j-occ_k_i).Eps1_k,q_ij.<u0_k_i|S^(j1)|u0_k+q_j> }
!          where Eps1_k,q_ij=<u0_k+q_j|H^(j2)-1/2(Eps_k+q_j-Eps_k_i)S^(j2)|u0_k_i>
!          At first call (when j1=j2), cs1c=<u0_k_i|S^(j1)|u0_k+q_j> is stored
!          For the next calls, it is re-used.
           if (has_dcwf.and.has_dcwf2.and.is_metal_or_qne0) then
             ABI_ALLOCATE(gvnl1,(2,npw1_k*nspinor))
             dotr=zero;doti=zero
             invocc=two/occ_k(iband)
             do jband=1,nband_k
               if ((ipert==ipert1.and.idir==idir1).or.(abs(occ_k(iband))>tol8)) then
                 gvnl1(:,1:npw1_k*nspinor)=cgq(:,1+npw1_k*nspinor*(jband-1)+icgq:npw1_k*nspinor*jband+icgq)
                 call dotprod_g(dot1r,dot1i,istwf_k,mpi_enreg,npw1_k*nspinor,2,gs1,gvnl1)
                 if (ipert==ipert1.and.idir==idir1) then
                   cs1c(1,jband,iband,ikpt_me)=dot1r
                   cs1c(2,jband,iband,ikpt_me)=dot1i
                 end if
               end if
               if (abs(occ_k(iband))>tol8) then
                 arg=eig_kq(jband)-eig_k(iband)
                 dot2r=cs1c(1,jband,iband,ikpt_me)
                 dot2i=cs1c(2,jband,iband,ikpt_me)
                 dotr=dotr+(dot1r*dot2r+dot1i*dot2i)*arg
                 doti=doti+(dot1i*dot2r-dot1r*dot2i)*arg
               end if
             end do
             if (is_metal.and.abs(occ_k(iband))>tol8) then
               do jband=1,nband_k
                 if (abs(rocceig(jband,iband))>tol8) then
                   ii=2*jband-1+(iband-1)*2*nband_k
                   arg=invocc*rocceig(jband,iband)*(eig_k(iband)-eig_kq(jband))
                   dot1r=eig1_k(ii);dot1i=eig1_k(ii+1)
                   dot2r=cs1c(1,jband,iband,ikpt_me)
                   dot2i=cs1c(2,jband,iband,ikpt_me)
                   dotr=dotr+arg*(dot1r*dot2r-dot1i*dot2i)
                   doti=doti+arg*(dot1i*dot2r+dot1r*dot2i)
                 end if
               end do
             end if
             dotr=quarter*dotr;doti=quarter*doti
!            Note: factor 2 (from d2E/dj1dj2=2E^(j1j2))
             d2ovl_k(1,idir1)=d2ovl_k(1,idir1)+wtk_k*occ_k(iband)*two*dotr
             d2ovl_k(2,idir1)=d2ovl_k(2,idir1)+wtk_k*occ_k(iband)*two*doti
             ABI_DEALLOCATE(gvnl1)
           end if

!          Build the matrix element <u0_k_i|H^(j1)-Eps_k_i.S^(j1)|u^(j2)_k,q_i>
!          and add contribution to DDB
           if (ipert1/=dtset%natom+1) then
             if (abs(occ_k(iband))>tol8) then
               call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1_k*nspinor,2,gh1,cwavef)
!              Case ipert1=natom+2 (electric field):
!              gh1 contains H^(j1)|u0_k_i> (VHxc constant) which corresponds
!              to i.d/dk in Eq. (38) of Gonze, PRB 55, 10355 (1997).
!              * if ipert==natom+2, we apply directly Eq. (38)
!              * if ipert/=natom+2, Born effective charges are minus D2E
               if ((ipert1==dtset%natom+2.and.ipert <=dtset%natom).or. &
&               (ipert ==dtset%natom+2.and.ipert1<=dtset%natom)) then
                 dotr=-dotr;doti=-doti
               end if
               d2nl_k(1,idir1)=d2nl_k(1,idir1)+wtk_k*occ_k(iband)*two*dotr
               d2nl_k(2,idir1)=d2nl_k(2,idir1)+wtk_k*occ_k(iband)*two*doti
             end if

!            Or compute localisation tensor (ddk)
!            See M. Veithen thesis Eq(2.5)
!            MT jan-2010: this is probably not correctly implemented for PAW !!!
!            missing terms due to S^(1) and S^(2)
           else
!            note: gh1 used as temporary space (to store idir ddk WF)
             if (idir==idir1) then
               gh1=cwavef
               if (need_ddk_file.and.ddkfil(idir1)/=0) then
                 call WffReadSkipRec(ierr,2,wffddk(idir1))
               end if
             else
               if (need_ddk_file.and.ddkfil(idir1)/=0) then
                 call WffReadSkipRec(ierr,1,wffddk(idir1))
                 call WffReadDataRec(gh1,ierr,2,npw1_k*nspinor,wffddk(idir1))
               else
                 gh1=zero
               end if
             end if
             if (abs(occ_k(iband))>tol8) then
               call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1_k*nspinor,2,gh1,cwavef)
               call dotprod_g(dot1r,dot1i,istwf_k,mpi_enreg,npw1_k*nspinor,2,cwave0,gh1)
               call dotprod_g(dot2r,dot2i,istwf_k,mpi_enreg,npw1_k*nspinor,2,cwavef,cwave0)
               dotr=dotr-(dot1r*dot2r-dot1i*dot2i)
               doti=doti-(dot1r*dot2i+dot1i*dot2r)
               d2nl_k(1,idir1)=d2nl_k(1,idir1)+wtk_k*occ_k(iband)*dotr/(nband_kocc*two)
               d2nl_k(2,idir1)=d2nl_k(2,idir1)+wtk_k*occ_k(iband)*doti/(nband_kocc*two)
             end if
           end if

!          Accumulate here 1st-order density change due to overlap operator changes (if any)
           if (has_drho) then
             if (abs(occ_k(iband))>tol8) then
!              Compute here delta_u^(j1)=-1/2 Sum_{j}[<u0_k+q_j|S^(j1)|u0_k_i>.|u0_k+q_j>]
!              (see PRB 78, 035105 (2008), Eq. (42))
               ABI_ALLOCATE(dcwavef,(2,npw1_k*nspinor))
               ABI_ALLOCATE(dcwaveprj,(dtset%natom,nspinor))
               call cprj_alloc(dcwaveprj,0,dimcprj)
               if (dtset%mkqmem==0) then
                 call getdc1(cgq_disk,cprjq_disk,dcwavef,dcwaveprj,ibgq,icgq,istwf_k,mcgq,&
&                 mcprjq,mpi_enreg,dtset%natom,nband_k,npw1_k,nspinor,1,dtset%ortalg,gs1)
               else
                 call getdc1(cgq     ,cprjq     ,dcwavef,dcwaveprj,ibgq,icgq,istwf_k,mcgq,&
&                 mcprjq,mpi_enreg,dtset%natom,nband_k,npw1_k,nspinor,1,dtset%ortalg,gs1)
               end if

!              Accumulate 1st-order density due to delta_u^(j1)
               counter=500*iband;option=1;wfcorr=0
               call accrho3(counter,cplex,cwave0,dcwavef,dcwavef,cwaveprj0_idir1,dcwaveprj,&
&               dimcprj,dimffnl,ffnlk,dimphkxred,lambda,dtfil%filstat,gbound,&
&               gs_hamkq,iband,idir1,ipert1,isppol,kg_k,kg1_k,kpg_k,kpoint,&
&               dtset%kptopt,psps%lmnmax,matblk,dtset%mgfft,mpi_enreg,dtset%natom,&
&               nband_k,1,nkpg,npw_k,npw1_k,nspinor,dtset%ntypat,dtset%ngfft(4),dtset%ngfft(5),&
&               dtset%ngfft(6),occ_k,option,dtset%paral_kgb,pawdrhoij1(:,idir1),ph3d,phkxred,&
&               dtset%prtvol,drhoaug1(:,:,:,idir1),tim_fourwf,usecprj,dum4,wfcorr,wtk_k)
               call cprj_free(dcwaveprj)
               ABI_DEALLOCATE(dcwavef)
               ABI_DEALLOCATE(dcwaveprj)
             end if
           end if

!          End of loops
         end do   ! idir1
       end do     ! iband

!      Accumulate contribution of this k-point
       d2nl (:,:,ipert1,idir,ipert)=d2nl (:,:,ipert1,idir,ipert)+d2nl_k (:,:)
       d2ovl(:,:,ipert1,idir,ipert)=d2ovl(:,:,ipert1,idir,ipert)+d2ovl_k(:,:)

!      Deallocations of arrays used for this k-point
       ABI_DEALLOCATE(gh1)
       ABI_DEALLOCATE(gs1)
       if (need_wfk)  then
         ABI_DEALLOCATE(cwave0)
       end if
       if (need_wf1)  then
         ABI_DEALLOCATE(cwavef)
       end if
       ABI_DEALLOCATE(kg_k)
       ABI_DEALLOCATE(kg1_k)
       ABI_DEALLOCATE(ylm_k)
       ABI_DEALLOCATE(ylm1_k)
       ABI_DEALLOCATE(ylmgr1_k)
       ABI_DEALLOCATE(kpg_k)
       ABI_DEALLOCATE(kpg1_k)
       ABI_DEALLOCATE(gbound)
       ABI_DEALLOCATE(d2nl_k)
       ABI_DEALLOCATE(d2ovl_k)
       ABI_DEALLOCATE(eig_k)
       ABI_DEALLOCATE(eig_kq)
       ABI_DEALLOCATE(eig1_k)
       ABI_DEALLOCATE(occ_k)
       if (is_metal)  then
         ABI_DEALLOCATE(doccde_k)
         ABI_DEALLOCATE(doccde_kq)
         ABI_DEALLOCATE(occ_kq)
         ABI_DEALLOCATE(rocceig)
       end if
       ABI_DEALLOCATE(dkinpw)
       ABI_DEALLOCATE(kinpw1)
       ABI_DEALLOCATE(phkxred)
       ABI_DEALLOCATE(ph3d)
       if (has_dcwf.and.dtset%mkqmem==0)  then
         ABI_DEALLOCATE(cgq_disk)
       end if
       if (ipert1>dtset%natom)  then
         ABI_DEALLOCATE(ffnl1_idir1)
       end if
       ABI_DEALLOCATE(ffnlk)
       ABI_DEALLOCATE(ffnl1)
       nullify(ffnlkq,ffnl1_idir1)
       if (usecprj==1) then
         call cprj_free(cwaveprj0)
         ABI_DEALLOCATE(cwaveprj0)
         if (ncpgr>1) then
           call cprj_free(cwaveprj0_idir1)
           ABI_DEALLOCATE(cwaveprj0_idir1)
         end if
       end if
       nullify(cwaveprj0_idir1)
!      Shift arrays
       bdtot_index=bdtot_index+nband_k
       bd2tot_index=bd2tot_index+2*nband_k**2
       if (dtset%mkmem/=0) then
         ibg=ibg+nspinor*nband_k
         icg=icg+npw_k*nspinor*nband_k
         ikg=ikg+npw_k
       end if
       if (dtset%mkqmem/=0) then
         ibgq=ibgq+nspinor*nband_k
         icgq=icgq+npw1_k*nspinor*nband_k
       end if
       if (dtset%mk1mem/=0) then
         ibg1=ibg1+nspinor*nband_k
         icg1=icg1+npw1_k*nspinor*nband_k
         ikg1=ikg1+npw1_k
       end if

!      End loop over KPTS
     end do ! End loop over K-POINTS

!    Transfer 1st-order density change due to overlap; also take into account the spin.
     if(has_drho) then
       do kdir1=1,mdir1
         idir1=jdir1(kdir1)
         call fftpac(isppol,nspden,cplex*dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),&
&         cplex*dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&         dtset%ngfft,drho1wfr(:,:,idir1),drhoaug1(:,:,:,idir1),1)
       end do
     end if

   end do   ! En loop over SPINS

!  Free memory used for this type of perturbation
   ABI_DEALLOCATE(ekb_typ)
   ABI_DEALLOCATE(sij_typ)
   ABI_DEALLOCATE(indlmn_typ)
   ABI_DEALLOCATE(e1kbfr)
   if (has_drho)  then
     ABI_DEALLOCATE(drhoaug1)
   end if
   if (ipert/=dtset%natom+1.and.ipert1/=dtset%natom+1) then
     do kdir1=1,mdir1
       idir1=jdir1(kdir1)
       call destroy_paw_ij(paw_ij10(:,idir1))
     end do
     ABI_DEALLOCATE(paw_ij10)
   end if

!  In case of parallelism, sum 1st-order density and occupation matrix over processors
   if (has_drho.and.mpi_enreg%paral_compil_kpt==1) then
     ABI_ALLOCATE(nlmn,(dtset%natom))
     nlmn(1:dtset%natom)=pawdrhoij1(1:dtset%natom,jdir1(1))%cplex &
&     *pawdrhoij1(1:dtset%natom,jdir1(1))%lmn2_size
     nsp=pawdrhoij1(1,jdir1(1))%nsppol;if (pawdrhoij1(1,jdir1(1))%nspden==4) nsp=4
     bufsz=(cplex*dtset%nfft*nspden+sum(nlmn)*nsp)*mdir1
     ABI_ALLOCATE(buffer,(bufsz))
     ii=cplex*dtset%nfft*nspden*mdir1
     buffer(1:ii)=reshape(drho1wfr,(/ii/))
     do kdir1=1,mdir1
       idir1=jdir1(kdir1)
       do iatom=1,dtset%natom
         do isppol=1,nsp
           buffer(ii+1:ii+nlmn(iatom))=pawdrhoij1(iatom,idir1)%rhoij_(:,isppol)
           ii=ii+nlmn(iatom)
         end do
       end do
     end do
     call timab(48,1,tsec)
     call xsum_mpi(buffer,bufsz,spaceworld,ierr)
     call timab(48,2,tsec)
     ii=cplex*dtset%nfft*nspden*mdir1
     drho1wfr(:,:,:)=reshape(buffer(1:ii),(/cplex*dtset%nfft,nspden,mdir1/))
     do kdir1=1,mdir1
       idir1=jdir1(kdir1)
       do iatom=1,dtset%natom
         do isppol=1,nsp
           pawdrhoij1(iatom,idir1)%rhoij_(:,isppol)=buffer(ii+1:ii+nlmn(iatom))
           ii=ii+nlmn(iatom)
         end do
       end do
     end do
     ABI_DEALLOCATE(buffer)
     ABI_DEALLOCATE(nlmn)
   end if

!  Compute second part of overlap contribution (due to VHxc^(j2)(tild_n+hat_n))
   if (has_drho) then

     ABI_ALLOCATE(drhor1,(cplex*nfftf,nspden))
     ABI_ALLOCATE(dnhat1,(cplex*nfftf,nspden))

!    LOOP OVER PERTURBATION DIRECTIONS
     do kdir1=1,mdir1
       idir1=jdir1(kdir1)

!      Build and symmetrize 1st-order density change due to change of overlap
       ABI_ALLOCATE(drho1wfg,(2,dtset%nfft))
       call symrhg(cplex,gprimd,irrzon1,mpi_enreg,dtset%nfft,dtset%nfft,dtset%ngfft,&
&       nspden,nsppol,nsym1,dtset%paral_kgb,phnons1,drho1wfg,drho1wfr(:,:,idir1),&
&       rprimd,symaf1,symrl1)
       if (dtset%pawstgylm/=0) then
         optfr=0
         ABI_ALLOCATE(paw_ij1_tmp,(dtset%natom))
         call nullify_paw_ij(paw_ij1_tmp)
         call pawfrnhat(cplex,gprimd,idir1,ipert1,mpi_enreg,dtset%natom,dtset%natom,nfftf,ngfftf,&
&         nspden,dtset%ntypat,optfr,paw_ij1_tmp,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,&
&         dtset%qptn,rprimd,ucvol,vpsp1,vtrial,xred)  ! vpsp1 is unused
         ABI_DEALLOCATE(paw_ij1_tmp)
       end if
       call pawmkrho(arg,cplex,gprimd,idir1,psps%indlmn,indsy1,ipert1,psps%lmnmax,&
&       mpi_enreg,dtset%natom,nspden,nsym1,dtset%ntypat,dtset%paral_kgb,pawang,pawfgr,&
&       pawfgrtab,-10001,pawdrhoij1(:,idir1),pawtab,dtset%qptn,drho1wfg,&
&       drho1wfr(:,:,idir1),drhor1,rprimd,symaf1,symrc1,dtset%typat,&
&       ucvol,xred,pawang_sym=pawang1,pawnhat=dnhat1,pawrhoij0=pawrhoij)
       ABI_DEALLOCATE(drho1wfg)

!      Compute plane-wave contribution to overlap contribution
!      This is subtle as it is a mix of Eq(79) and Eq(80) of PRB 78, 035105 (2008)
!      Details:
!      The VH(tild_nZc)^(1) term of Eq(79) is:
!      <VH(tild_nZc)^(j2)|delta_tild_rho^(j1)>            = <vpsp1|drhor1-dnhat1>
!      The first term of Eq(80) is:
!      <VHxc^(j2)|delta_tild_rho^(j1)+delta_hat_rho^(j1)> = <vtrial1-vpsp1|drhor1>
!      The addition of these two terms gives:
!      <vtrial1|drhor1>-<vpsp1|dnhat1>
!      And this is more subtle when usexcnhat=0
       call dotprod_vn(cplex,drhor1,dot1r,dot1i,mpi_enreg,nfftf,nfftot,nspden,2,vtrial1,ucvol)
       if (usexcnhat/=0) then
         call dotprod_vn(cplex,dnhat1,dot2r,dot2i,mpi_enreg,nfftf,nfftot,1   ,2,vpsp1  ,ucvol)
       else
         ABI_ALLOCATE(vtmp1,(cplex*nfftf,nspden))
         do ispden=1,nspden
           vtmp1(:,ispden)=vtrial1(:,ispden)-vhartr1(:)
         end do
         call dotprod_vn(cplex,dnhat1,dot2r,dot2i,mpi_enreg,nfftf,nfftot,nspden,2,vtmp1,ucvol)
         ABI_DEALLOCATE(vtmp1)
       end if
       dotr=dot1r-dot2r;doti=dot1i-dot2i

!      Compute on-site contributions to overlap contribution
!      (two last terms of Eq(80) of PRB 78, 035105 (2008))
!      (note: Dij^(j2) and Vxc^(j2) are computed for ipert at first call)
       call pawnstd2e(epawnst,ipert,ipert1,mpi_enreg,dtset%natom,dtset%natom,dtset%ntypat,&
&       nzlmopt_ipert,nzlmopt_ipert1,paw_an,paw_an1,paw_ij1,pawang,dtset%pawprtvol,&
&       pawrad,pawrhoij1,pawdrhoij1(:,idir1),pawtab,dtset%pawxcdev,dtset%xclevel)

!      Accumulate in 2nd-order matrix:
!      Note: factor 2 (from d2E/dj1dj2=2E^(j1j2)) eliminated by factor 1/2
!      has to take the complex conjugate because we want here Int[VHxc^(j1)^*.delta_rho^(j2)]
       dotr=dotr+epawnst(1);doti=-(doti+epawnst(2))
       if (ipert==dtset%natom+2) then
         d2ovl_drho(1,idir1,ipert1,idir,ipert)=-dotr
         d2ovl_drho(2,idir1,ipert1,idir,ipert)=-doti
       else
         d2ovl_drho(1,idir1,ipert1,idir,ipert)=dotr
         d2ovl_drho(2,idir1,ipert1,idir,ipert)=doti
       end if

     end do ! End loop over perturbation directions

!    Free no more needed memory
     ABI_DEALLOCATE(drhor1)
     ABI_DEALLOCATE(dnhat1)
     ABI_DEALLOCATE(drho1wfr)
     do iatom=1,dtset%natom
       if (pawfgrtab(iatom)%nhatfr_allocated>0)  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%nhatfr)
       end if
       pawfgrtab(iatom)%nhatfr_allocated=0
     end do
     do kdir1=1,mdir1
       idir1=jdir1(kdir1)
       call rhoij_free(pawdrhoij1(:,idir1))
     end do
     ABI_DEALLOCATE(pawdrhoij1)
   end if ! has_drho

!  End loop over perturbations (j1)
 end do

!Final deallocations
 ABI_DEALLOCATE(ch1c)
 if (is_metal_or_qne0) ABI_DEALLOCATE(cs1c)
 call destroy_hamiltonian(gs_hamkq)
 if (ipert/=dtset%natom+1) then
   if (.not.has_ipert_dijh) then
     do iatom=1,dtset%natom
       if (associated(paw_ij1(iatom)%dijhartree))  then
         ABI_DEALLOCATE(paw_ij1(iatom)%dijhartree)
       end if
       paw_ij1(iatom)%has_dijhartree=0
     end do
   end if
   if (.not.has_ipert_vxc) then
     do iatom=1,dtset%natom
       if (associated(paw_an1(iatom)%vxc1))   then
         ABI_DEALLOCATE(paw_an1(iatom)%vxc1)
       end if
       if (associated(paw_an1(iatom)%vxct1))  then
         ABI_DEALLOCATE(paw_an1(iatom)%vxct1)
       end if
       paw_an1(iatom)%has_vxc=0
     end do
   end if
 end if

!In case of parallelism, sum over processors
 if (mpi_enreg%paral_compil_kpt==1)then
   call timab(161,1,tsec)
   call leave_test()
   call timab(161,2,tsec)
   ABI_ALLOCATE(buffer,(12*mpert))
   buffer(        1:        6*mpert)=reshape(d2nl (:,:,:,idir,ipert),(/6*mpert/))
   buffer(6*mpert+1:6*mpert+6*mpert)=reshape(d2ovl(:,:,:,idir,ipert),(/6*mpert/))
   call timab(48,1,tsec)
   call xsum_mpi(buffer,12*mpert,spaceworld,ierr)
   call timab(48,2,tsec)
   d2nl (:,:,:,idir,ipert)=reshape(buffer(        1:        6*mpert),(/2,3,mpert/))
   d2ovl(:,:,:,idir,ipert)=reshape(buffer(6*mpert+1:6*mpert+6*mpert),(/2,3,mpert/))
   ABI_DEALLOCATE(buffer)
 end if

!Build complete d2ovl matrix
 d2ovl=d2ovl+d2ovl_drho

!Close the ddk WF files
 if (has_ddk_file) then
   do kdir1=1,mdir1
     idir1=jdir1(kdir1)
     if (ddkfil(idir1)/=0)then
       call WffClose(wffddk(idir1),ierr)
     end if
   end do
 end if
 ABI_DEALLOCATE(jpert1)
 ABI_DEALLOCATE(jdir1)

!Symmetrize the contributions, as was needed for the forces in a GS calculation
 ABI_ALLOCATE(work,(2,3,dtset%natom))
 do ipert1=1,dtset%natom
   do idir1=1,3
     work(:,idir1,ipert1)=d2nl(:,idir1,ipert1,idir,ipert)
   end do
 end do
 call sygra3(dtset%natom,d2nl(:,:,:,idir,ipert),work,indsy1,ipert,nsym1,dtset%qptn,symrc1)
 do ipert1=1,dtset%natom
   do idir1=1,3
     work(:,idir1,ipert1)=d2ovl(:,idir1,ipert1,idir,ipert)
   end do
 end do
 call sygra3(dtset%natom,d2ovl(:,:,:,idir,ipert),work,indsy1,ipert,nsym1,dtset%qptn,symrc1)
 ABI_DEALLOCATE(work)

!Must also symmetrize the electric/magnetic field perturbation response !
!Note: d2ovl is not symetrized because it is zero for electric/magnetic field perturbation
 if (has_ddk_file) then
   ABI_ALLOCATE(d2nl_elfd,(2,3))
   ABI_ALLOCATE(d2nl_mgfd,(2,3))
!  There should not be any imaginary part, but stay general (for debugging)
   d2nl_elfd (:,:)=d2nl(:,:,dtset%natom+2,idir,ipert)
   d2nl_mgfd (:,:)=d2nl(:,:,dtset%natom+5,idir,ipert)
   do ii=1,3
     sumelfd(:)=zero
     summgfd(:)=zero
     do ia=1,nsym1
       do jj=1,3
         if(symrl1(ii,jj,ia)/=0)then
           if(ddkfil(jj)==0)then
             blkflg(ii,dtset%natom+2,idir,ipert)=0
             blkflg(ii,dtset%natom+5,idir,ipert)=0
           end if
         end if
       end do
       sumelfd(:)=sumelfd(:)+dble(symrl1(ii,1,ia))*d2nl_elfd(:,1) &
&       +dble(symrl1(ii,2,ia))*d2nl_elfd(:,2) &
&       +dble(symrl1(ii,3,ia))*d2nl_elfd(:,3)
       summgfd(:)=summgfd(:)+dble(symrl1(ii,1,ia))*d2nl_mgfd(:,1) &
&       +dble(symrl1(ii,2,ia))*d2nl_mgfd(:,2) &
&       +dble(symrl1(ii,3,ia))*d2nl_mgfd(:,3)
     end do
     d2nl(:,ii,dtset%natom+2,idir,ipert)=sumelfd(:)/dble(nsym1)
     d2nl(:,ii,dtset%natom+5,idir,ipert)=summgfd(:)/dble(nsym1)
   end do
   ABI_DEALLOCATE(d2nl_elfd)
   ABI_DEALLOCATE(d2nl_mgfd)
 end if

!Store the diagonal part of the matrix in the 2nd-order energy non-stationnary expression
 eovl1=d2ovl(1,idir,ipert,idir,ipert)

 call timab(566,2,tsec)

 DBG_EXIT("COLL")

end subroutine nstpaw3
!!***
