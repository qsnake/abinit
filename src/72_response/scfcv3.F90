!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfcv3
!! NAME
!! scfcv3
!!
!! FUNCTION
!! Conducts set of passes or overall iterations of preconditioned
!! conjugate gradient algorithm to converge wavefunctions to
!! optimum and optionally to compute mixed derivatives of energy.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG, DRH, MB, XW, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=pw coefficients of GS wavefunctions at k.
!!  cgq(2,mpw1*nspinor*mband*mkqmem*nsppol)=pw coefficients of GS wavefunctions at k+q.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk>
!!  cprjq(natom,nspinor*mband*mkqmem*nsppol*usecprj)= wave functions at k+q
!!              projected with non-local projectors: cprjq=<p_i|Cnk+q>
!!  cpus= cpu time limit in seconds
!!  dimcprj(natom*usepaw)=array of dimensions of arrays cprj, cprjq (ordered by atom-type)
!!  doccde_rbz(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy
!!  docckqde(mband*nkpt_rbz*nsppol)=derivative of occkq wrt the energy
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eew=2nd derivative of Ewald energy (hartree)
!!  efrhar=Contribution from frozen-wavefunction, hartree energy,
!!           to the second-derivative of total energy.
!!  efrkin=Contribution from frozen-wavefunction, kinetic energy,
!!           to the second-derivative of total energy.
!!  efrloc=Contribution from frozen-wavefunction, local potential,
!!           to the second-derivative of total energy.
!!  efrnl=Contribution from frozen-wavefunction, non-local potential,
!!           to the second-derivative of total energy.
!!  efrx1=Contribution from frozen-wavefunction, xc core correction(1),
!!           to the second-derivative of total energy.
!!  efrx2=Contribution from frozen-wavefunction, xc core correction(2),
!!           to the second-derivative of total energy.
!!  eigenq(mband*nkpt_rbz*nsppol)=GS eigenvalues at k+q (hartree)
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  eii=2nd derivative of pseudopotential core energy (hartree)
!!  fermie=fermi energy (Hartree)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  idir=direction of the current perturbation
!!  indkpt1(nkpt_rbz)=non-symmetrized indices of the k-points
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  ipert=type of the perturbation
!!  irrzon1(nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data for RF symmetries
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates at k
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points.
!!  kxc(nfftf,nkxc)=exchange and correlation kernel (see rhohxc.f)
!!  mgfftf=maximum size of 1D FFTs for the "fine" grid (see NOTES in respfn.F90)
!!  mkmem =number of k points which can fit in memory (GS data); 0 if use disk
!!  mkqmem =number of k+q points which can fit in memory (GS data); 0 if use disk
!!  mk1mem =number of k points which can fit in memory (RF data); 0 if use disk
!!  mpert=maximum number of ipert
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw for wfs at k.
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  nattyp(ntypat)= # atoms of each type.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point, for each polarization
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  nfftf=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!!  ngfftf(1:18)=integer array with FFT box dimensions and other for the "fine" grid (see NOTES in respfn.F90)
!!  nkpt=number of k points in the full BZ
!!  nkpt_rbz=number of k points in the reduced BZ for this perturbation
!!  nkxc=second dimension of the kxc array.
!!  mpi_enreg=informations about MPI parallelization
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsym1=number of symmetry elements in space group consistent with perturbation
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used otherwise, nfftf
!!  occkq(mband*nkpt_rbz*nsppol)=occupation number for each band (often 2)
!!   at each k+q point of the reduced Brillouin zone.
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band (often 2)
!!   at each k point of the reduced Brillouin zone.
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh for the GS
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawang1 <type(pawang_type)>=pawang datastructure containing only the symmetries preserving the perturbation
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid for the GS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data for the GS
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pertcase=fuill index of the perturbation
!!  phnons1(2,nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic transl. phases, for RF symmetries
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  ph1df(2,3*(2*mgfftf+1)*natom)=one-dimensional structure factor information for the "fine" grid
!!  prtbbb=if 1, band-by-band decomposition (also dim of d2bbb)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qphon(3)=reduced coordinates for the phonon wavelength
!!  rhog(2,nfftf)=array for Fourier transform of GS electron density
!!  rhor(nfftf,nspden)=array for GS electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symaf1(nsym1)=anti(ferromagnetic) part of symmetry operations
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  symrl1(3,3,nsym1)=symmetry operations in real space in terms
!!   of primitive translations
!!  usecprj= 1 if cprj, cprjq arrays are stored in memory
!!  useylmgr = 1 if ylmgr  array is allocated
!!  useylmgr1= 1 if ylmgr1 array is allocated
!!  wffddk=struct info for ddk file
!!  wffnew=struct info for 1WF at exit, if mk1mem=0
!!  wffnow=struct info for 1WF at start, if mk1mem=0
!!  wfftgs=struct info for GS WF at start, if mkmem=0
!!  wfftkq=struct info for k+q GS WF at start, if mkqmem=0
!!  vpsp1(cplex*nfftf)=first-order derivative of the ionic potential
!!  vtrial(nfftf,nspden)=GS potential (Hartree).
!!  vxc(nfftf,nspden)=Exchange-Correlation GS potential (Hartree)
!!  wtk_rbz(nkpt_rbz)=weight for each k point in the reduced Brillouin zone
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density, see n3xccc
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k+q point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm*useylmgr)= gradients of real spherical harmonics at k
!!  ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*useylm*useylmgr1)= gradients of real spherical harmonics at k+q
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert)=flags for each element of the 2DTE (=1 if computed)
!!  cg1_active(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF
!!    wavefunctions at k,q. They are orthogonalized to the active.
!!  d2bbb(2,3,3,mpert,mband,mband*prtbbb)=band by band decomposition of some
!!       second order derivatives
!!  d2lo(2,mpert,3,mpert)=local contributions to the 2DTEs
!!  d2nl(2,mpert,3,mpert)=non-local contributions to the 2DTEs
!!  d2ovl(2,mpert,3,mpert*usepaw)=1st-order change of WF overlap contributions to the 2DTEs
!!  eberry=energy associated with Berry phase
!!  edocc=correction to 2nd-order total energy coming from changes of occupation
!!  eeig0=0th-order eigenenergies part of 2nd-order total energy
!!  ehart01=inhomogeneous 1st-order Hartree part of 2nd-order total energy
!!    for strain perturbation only (zero otherwise, and not used)
!!  ehart1=1st-order Hartree part of 2nd-order total energy
!!  eigen1(2*mband*mband*nkpt_rbz*nsppol)=array for holding eigenvalues (hartree)
!!  ek0=0th-order kinetic energy part of 2nd-order total energy.
!!  ek1=1st-order kinetic energy part of 2nd-order total energy.
!!  eloc0=0th-order local (psp+vxc+Hart) part of 2nd-order total energy
!!  elpsp1=1st-order local pseudopot. part of 2nd-order total energy.
!!  enl0=0th-order nonlocal pseudopot. part of 2nd-order total energy.
!!  enl1=1st-order nonlocal pseudopot. part of 2nd-order total energy.
!!  eovl1=1st-order change of wave-functions overlap, part of 2nd-order energy
!!        PAW only - Eq(79) and Eq(80) of PRB 78, 035105 (2008)
!!  epaw1=1st-order PAW on-site part of 2nd-order total energy.
!!  etotal=total energy (sum of 7 contributions) (hartree)
!!  exc1=1st-order exchange-correlation part of 2nd-order total energy.
!!  gh1c_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)= set of <G|H^{(1)}|nK>
!!  gh0c1_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)= set of <G|H^{(0)}|\Psi^{(1)}>
!!      The wavefunction is orthogonal to the active space (for metals). It is not
!!      coherent with cg1.
!!  resid(mband*nkpt_rbz*nsppol)=residuals for each band over all k points
!!   of the reduced Brillouin zone, and spins
!!  residm=maximum value from resid array (except for nbdbuf highest bands)
!!
!! SIDE EFFECTS
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=updated wavefunctions (ortho. to occ. states);
!!    if mk1mem<=nkpt_rbz, these are kept in a disk file, see wffnow and wffnew.
!!  initialized= if 0 the initialization of the RF run is not yet finished
!!  mpi_enreg=informations about MPI parallelization
!!  rhog1(2,nfftf)=array for Fourier transform of RF electron density
!!  rhor1(cplex*nfftf,nspden)=array for RF electron density in electrons/bohr**3.
!!  === if psps%usepaw==1
!!    pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!
!! TODO
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      ab6_mixing_deallocate,ab6_mixing_new,ab6_mixing_use_disk_cache,appdig
!!      bec3,cprj_alloc,cprj_free,destroy_paw_an,destroy_paw_ij,die3,ebp3,edie3
!!      etot3,fourdp,getcut,init_paw_an,init_paw_ij,initberry3,ioarr,leave_new
!!      leave_test,metric,newfermie1,newvtr3,nselt3,nstdy3,nstpaw3
!!      nullify_paw_an,nullify_paw_ij,pawdenpot,pawdij,pawfrnhat,pawmknhat
!!      qmatrix,rhofermi3,rhotov3,scprqt,status,symdij,timab,vtorho3,wffclose
!!      wrtout,xcomm_world,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine scfcv3(atindx,atindx1,blkflg,cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cpus,&
&  dimcprj,dim_eig2rf,doccde_rbz,docckqde,dtfil,dtset,&
&  d2bbb,d2lo,d2nl,d2ovl,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&
&  ehart01,ehart1,eigenq,eigen0,eigen1,eii,ek0,ek1,eloc0,elpsp1,&
&  enl0,enl1,eovl1,epaw1,etotal,exc1,fermie,gh0c1_set,gh1c_set,hdr,idir,indkpt1,&
&  indsy1,initialized,ipert,irrzon1,istwfk_rbz,&
&  kg,kg1,kpt_rbz,kxc,mgfftf,mkmem,mkqmem,mk1mem,&
&  mpert,mpi_enreg,mpsang,mpw,mpw1,nattyp,nband_rbz,ncpgr,&
&  nfftf,ngfftf,nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,&
&  nsym1,n3xccc,occkq,occ_rbz,&
&  paw_an,paw_ij,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoij,pawrhoij1,pawtab,&
&  pertcase,phnons1,ph1d,ph1df,&
&  prtbbb,psps,qphon,resid,residm,rhog,rhog1,&
&  rhor,rhor1,rprimd,symaf1,symrc1,symrl1,&
&  usecprj,useylmgr,useylmgr1,wffddk,wffnew,wffnow,wfftgs,wfftkq,vpsp1,vtrial,vxc,&
&  wtk_rbz,wvl,xccc3d1,xred,ylm,ylm1,ylmgr,ylmgr1)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_errors
 use m_ab6_mixing
 use m_paw_toolbox
 use m_wffile
 use m_efield

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scfcv3'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_66_paw
 use interfaces_67_common
 use interfaces_72_response, except_this_one => scfcv3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!no_abirules
!Needed for integer arrays
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
!---
 integer,intent(in) :: cplex,dim_eig2rf,idir,ipert,mgfftf,mk1mem,mkmem,mkqmem
 integer,intent(in) :: mpert,mpsang,mpw,mpw1,n3xccc,ncpgr,nfftf
 integer,intent(in) :: nkpt,nkpt_rbz,nkxc,nspden
 integer,intent(in) :: nsym1,pertcase,prtbbb,usecprj,useylmgr,useylmgr1
 integer,intent(inout) :: initialized
! nfft**(1-1/nsym1) is 1 if nsym1==1, and nfft otherwise
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom),dimcprj(dtset%natom*psps%usepaw),ngfftf(18)
 integer,intent(out) :: blkflg(3,mpert,3,mpert)
 integer,intent(in) :: indkpt1(nkpt_rbz),indsy1(4,nsym1,dtset%natom)
 integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 integer,intent(in) :: istwfk_rbz(nkpt_rbz)
 integer,intent(in) :: kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem),nattyp(psps%ntypat)
 integer,intent(in) :: nband_rbz(nkpt_rbz*dtset%nsppol)
 integer,intent(in) :: npwar1(nkpt_rbz),npwarr(nkpt_rbz)
 integer,intent(in) :: symaf1(nsym1),symrc1(3,3,nsym1),symrl1(3,3,nsym1)
 real(dp),intent(in) :: cpus,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2
 real(dp),intent(out) :: eberry,edocc,eeig0,ehart01,ehart1,ek0,ek1,eloc0,elpsp1,enl0
 real(dp),intent(out) :: enl1,eovl1,epaw1,etotal,exc1,residm
 real(dp),intent(in) :: eii
 real(dp),intent(inout) :: fermie
 real(dp),intent(in) :: qphon(3)
! nfft**(1-1/nsym1) is 1 if nsym1==1, and nfft otherwise
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*mkmem*dtset%nsppol)
 real(dp),intent(inout) :: cg1(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol)
 real(dp),intent(out) :: cg1_active(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*dim_eig2rf)
 real(dp),intent(out) :: gh1c_set(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*dim_eig2rf)
 real(dp),intent(out) :: gh0c1_set(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*dim_eig2rf)
 real(dp),intent(in) :: cgq(2,mpw1*dtset%nspinor*dtset%mband*mkqmem*dtset%nsppol)
 real(dp),intent(out) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb)
 real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert),d2nl(2,3,mpert,3,mpert)
 real(dp),intent(out) :: d2ovl(2,3,mpert,3,mpert*psps%usepaw)
 real(dp),intent(in) :: doccde_rbz(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: docckqde(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: eigen0(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(out) :: eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: eigenq(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz),kxc(nfftf,nkxc)
 real(dp),intent(in) :: occ_rbz(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: occkq(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom),ph1df(2,3*(2*mgfftf+1)*dtset%natom)
 real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp),intent(out) :: resid(dtset%mband*nkpt_rbz*nspden)
 real(dp),intent(in) :: rhog(2,nfftf),rhor(nfftf,nspden),rprimd(3,3)
 real(dp),intent(inout) :: rhog1(2,nfftf),rhor1(cplex*nfftf,nspden),xred(3,dtset%natom)
 real(dp),target,intent(inout) :: vtrial(nfftf,nspden)
 real(dp),intent(in) :: vpsp1(cplex*nfftf),vxc(nfftf,nspden)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xccc3d1(cplex*n3xccc)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang*psps%useylm*useylmgr)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*psps%useylm*useylmgr1)
 type(cprj_type),intent(in) :: cprj(dtset%natom,dtset%nspinor*dtset%mband*mkmem*dtset%nsppol*usecprj)
 type(cprj_type),intent(in) :: cprjq(dtset%natom,dtset%nspinor*dtset%mband*mkqmem*dtset%nsppol*usecprj)
 type(datafiles_type),intent(in) :: dtfil
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang,pawang1
 type(pawfgr_type),intent(in) :: pawfgr
 type(paw_an_type),intent(inout) :: paw_an(dtset%natom*psps%usepaw)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij1(dtset%natom*psps%usepaw)
 type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(MPI_type),intent(inout) :: mpi_enreg
 type(wffile_type),intent(inout) :: wffddk,wfftgs,wfftkq
 type(wffile_type),intent(inout) :: wffnow,wffnew

!Local variables-------------------------------
!scalars
 integer,parameter :: level=12,response=1
 integer :: accessfil,afford,choice,cplex_dij,dbl_nnsclo,fformr,fformv
 integer :: has_dijfr,iatom,ider,ierr,iexit,ii,errid,denpot
 integer :: iprcel,ir,iscf10_mod,iscf_mod,ispden,ispmix
 integer :: istep,itypat,izero,lmn2_size,me,mgfftdiel,mvdum
 integer :: nfftdiel,nfftmix,nfftot,npawmix,npwdiel,nstep,nzlmopt
 integer :: optene,optfr,option,optres,prtfor,quit,quit_sum,qzero,rdwr
 integer :: rdwrpaw,spaceComm,usexcnhat,v_size
 real(dp) :: boxcut,deltae,diffor,dum,ecut,ecutf,elast
 real(dp) :: epawdc1_dum,evar,fe1fixed,fermie1,gsqcut,maxfor,res2
 real(dp) :: ucvol,vxcavg
 logical :: ex
 character(len=500) :: msg
 character(len=fnlen) :: fi1o
 type(ab6_mixing_object) :: mix
! TODO : this field is not used - why?
 type(efield_type) :: dtefield
!arrays
 integer :: ngfftmix(18)
 real(dp) :: dielar(7),favg(3),gmet(3,3),gprimd(3,3),k0(3)
 real(dp) :: rmet(3,3),tollist(12),tsec(2)
 real(dp),allocatable :: dielinv(:,:,:,:,:)
 real(dp),allocatable :: fcart(:,:),nhat1(:,:),nhatgr_dum(:,:,:)
 real(dp),allocatable :: nvresid1(:,:),rhorfermi(:,:),susmat(:,:,:,:,:)
 real(dp),allocatable :: vhartr1(:)
 real(dp),pointer :: vtmp(:,:),vxc1(:,:)
 real(dp),allocatable,target :: vtrial1(:,:)
 type(cprj_type),allocatable :: cprj1(:,:)
 type(paw_an_type),allocatable :: paw_an1(:)
 type(paw_ij_type),allocatable :: paw_ij1(:)
!no_abirules
 integer,allocatable :: pwindall(:,:,:)
 real(dp),allocatable ::qmat(:,:,:,:,:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 call timab(120,1,tsec)
 call timab(154,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if dtset%prtvol==-level
 if(dtset%prtvol==-level)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' scfcv3: enter '
   call wrtout(std_out,msg,'COLL')
 end if

!Init me
 call xme_init(mpi_enreg,me)

!
!If dtset%accesswff == 2 set all array outputs to netcdf format
!
 accessfil = 0
 if (dtset%accesswff == IO_MODE_NETCDF) then
   accessfil = 1
 end if
 if (dtset%accesswff == IO_MODE_ETSF) then
   accessfil = 3
 end if

!Save some variables from dataset definition
 ecut=dtset%ecut
 ecutf=ecut;if (psps%usepaw==1.and.pawfgr%usefinegrid==1) ecutf=dtset%pawecutdg
 iprcel=dtset%iprcel
 tollist(1)=dtset%tolmxf;tollist(2)=dtset%tolwfr
 tollist(3)=dtset%toldff;tollist(4)=dtset%toldfe
 tollist(6)=dtset%tolvrs;tollist(7)=dtset%tolrff
 nstep=dtset%nstep
 iscf_mod=dtset%iscf
 iscf10_mod=mod(iscf_mod,10)
 qzero=0;if(qphon(1)**2+qphon(2)**2+qphon(3)**2 < tol14) qzero=1

!The value of iscf must be modified if ddk perturbation, (maybe eventually natom+5 too) see loper3.f
 if (ipert==dtset%natom+1) iscf_mod=-3

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Some variables need to be initialized/nullify at start
 quit=0 ; dbl_nnsclo=0 ; elast=zero
!This might be taken away later
 edocc=zero ; eeig0=zero ; ehart01=zero ; ehart1=zero ; ek0=zero ; ek1=zero
 eloc0=zero ; elpsp1=zero ; enl0=zero ; enl1=zero ; eovl1=zero; exc1=zero
 deltae=zero ; fermie1=zero ; epaw1=zero
 optres=merge(0,1,abs(iscf_mod)<10)
 usexcnhat=0

!Examine tolerance criteria, and eventually  print a line to the output
!file (with choice=1, the only non-dummy arguments of scprqt are
!nstep, tollist and iscf - still, diffor,res2,prtfor,fcart are here initialized to 0)
 choice=1 ; prtfor=0 ; diffor=zero ; res2=zero
 ABI_ALLOCATE(fcart,(3,dtset%natom))
 call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
& etotal,favg,fcart,fermie,dtfil%fnametmp_eig,dtfil%filnam_ds(1),&
& 1,iscf_mod,istep,kpt_rbz,maxfor,&
& mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
& nstep,occ_rbz,0,prtfor,0,&
& quit,res2,resid,residm,response,&
& tollist,psps%usepaw,vxcavg,wtk_rbz,xred)

!Various allocations (potentials, gradients, ...)
 call status(0,dtfil%filstat,iexit,level,'allocate      ')
 ABI_ALLOCATE(vhartr1,(cplex*nfftf))
 ABI_ALLOCATE(vtrial1,(cplex*nfftf,nspden))

!Allocations/initializations for PAW only
 if(psps%usepaw==1) then
   usexcnhat=maxval(pawtab(:)%usexcnhat)
!  1st-order compensation density
   ABI_ALLOCATE(nhat1,(cplex*nfftf,dtset%nspden))
   if (nstep==0) nhat1=zero
!  Projections of 1-st order WF on nl projectors
   ABI_ALLOCATE(cprj1,(dtset%natom,dtset%nspinor*dtset%mband*dtset%mk1mem*dtset%nsppol*usecprj))
   if (usecprj==1.and.dtset%mk1mem/=0) then
     call cprj_alloc(cprj1,0,dimcprj)
   end if
!  1st-order arrays/variables related to the PAW spheres
   ABI_ALLOCATE(paw_an1,(dtset%natom))
   ABI_ALLOCATE(paw_ij1,(dtset%natom))
   call nullify_paw_an(paw_an1)
   call nullify_paw_ij(paw_ij1)
   cplex_dij=max(cplex,dtset%nspinor)
   has_dijfr=0;if (ipert<=dtset%natom.or.ipert==dtset%natom+2.or.ipert==dtset%natom+5) has_dijfr=1
   call init_paw_an(dtset%natom,dtset%ntypat,0,dtset%nspden,cplex,dtset%pawxcdev,dtset%typat,&
&   pawang,pawtab,paw_an1)
   call init_paw_ij(paw_ij1,cplex,cplex_dij,dtset%nspinor,dtset%nsppol,dtset%nspden,0,dtset%natom,&
&   dtset%ntypat,dtset%typat,pawtab,has_dij=1,has_dijfr=has_dijfr)
 end if ! PAW

!Several parameters and arrays for the SCF mixing:
!These arrays are needed only in the self-consistent case
 if (iscf_mod>0.or.iscf_mod==-3) then
   ABI_ALLOCATE(nvresid1,(cplex*nfftf,dtset%nspden))
   if (nstep==0) nvresid1=zero
 end if
 if(nstep>0 .and. iscf_mod>0) then
   dielar(1)=dtset%diecut;dielar(2)=dtset%dielng
   dielar(3)=dtset%diemac;dielar(4)=dtset%diemix
   dielar(5)=dtset%diegap;dielar(6)=dtset%dielam
   dielar(7)=dtset%diemix;if (dtset%iscf>=10) dielar(7)=dtset%diemixmag
!  Additional allocation for mixing within PAW
   npawmix=0
   if(psps%usepaw==1) then
     do iatom=1,dtset%natom
       itypat=dtset%typat(iatom)
       lmn2_size=pawtab(itypat)%lmn2_size
       pawrhoij1(iatom)%use_rhoijres=1
       ABI_ALLOCATE(pawrhoij1(iatom)%rhoijres,(pawrhoij1(iatom)%cplex*lmn2_size,pawrhoij1(iatom)%nspden))
       do ispden=1,pawrhoij1(iatom)%nspden
         pawrhoij1(iatom)%rhoijres(:,ispden)=zero
       end do
       ABI_ALLOCATE(pawrhoij1(iatom)%kpawmix,(pawtab(itypat)%lmnmix_sz))
       pawrhoij1(iatom)%lmnmix_sz=pawtab(itypat)%lmnmix_sz
       pawrhoij1(iatom)%kpawmix=pawtab(itypat)%kmix
       npawmix=npawmix+pawrhoij1(iatom)%nspden*pawtab(itypat)%lmnmix_sz*pawrhoij1(iatom)%cplex
     end do
   end if
   denpot = AB6_MIXING_POTENTIAL
   if (dtset%iscf > 10) denpot = AB6_MIXING_DENSITY
   if (psps%usepaw==1.and.dtset%pawmixdg==0) then
     ispmix=AB6_MIXING_FOURRIER_SPACE;nfftmix=dtset%nfft;ngfftmix(:)=dtset%ngfft(:)
   else
     ispmix=AB6_MIXING_REAL_SPACE;nfftmix=nfftf;ngfftmix(:)=ngfftf(:)
   end if
   if (iscf10_mod == 5 .or. iscf10_mod == 6) then
     call ab6_mixing_new(mix, iscf10_mod, denpot, cplex, &
&     nfftf, dtset%nspden, npawmix, errid, msg, dtset%npulayit)
   else
     call ab6_mixing_new(mix, iscf10_mod, denpot, max(cplex, ispmix), &
&     nfftmix, dtset%nspden, npawmix, errid, msg, dtset%npulayit)
   end if
   if (errid /= AB6_NO_ERROR) then
     call wrtout(std_out, msg, 'COLL')
     call leave_new('COLL')
   end if
   if (dtset%mffmem == 0) then
     call ab6_mixing_use_disk_cache(mix, dtfil%fnametmp_fft)
   end if
 end if ! iscf, nstep

!Here, allocate arrays for computation of susceptibility and dielectric matrix or for TDDFT
 if( (nstep>0 .and. iscf_mod>0) .or. iscf_mod==-1 ) then
!  Here, for TDDFT, artificially set iprcel . Also set a variable to reduce
!  the memory needs.
   afford=1
   if(iscf_mod==-1) then
     iprcel=21
     afford=0
   end if
   npwdiel=1
   mgfftdiel=1
   nfftdiel=1
!  Now, performs allocation
!  CAUTION : the dimensions are still those of GS, except for phnonsdiel
   ABI_ALLOCATE(dielinv,(2,npwdiel*afford,nspden,npwdiel,nspden))
   ABI_ALLOCATE(susmat,(2,npwdiel*afford,nspden,npwdiel,nspden))
 end if

!Initialize Berry-phase related stuffs
 if (dtset%berryopt == 4) then
   ABI_ALLOCATE(pwindall,(max(mpw,mpw1)*mkmem,8,3))
   call initberry3(dtefield,dtfil,dtset,gmet,kg,kg1,dtset%mband,mkmem,mpi_enreg,&
&   mpw,mpw1,nkpt,npwarr,npwar1,dtset%nsppol,occ_rbz,pwindall,rprimd)
!  calculate inverse of the overlap matrix
   ABI_ALLOCATE(qmat,(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3))
   call qmatrix(cg,dtefield,qmat,mpw,mpw1,mkmem,dtset%mband,npwarr,nkpt,dtset%nspinor,dtset%nsppol,pwindall)
 end if

!Compute large sphere cut-off gsqcut
 if (psps%usepaw==1) then
   write(msg,'(2a)') ch10,' FFT (fine) grid used for densities/potentials:'
   call wrtout(std_out,msg,'COLL')
 end if
 k0(:)=zero
 call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,std_out,k0,ngfftf)

 call timab(154,2,tsec)

!######################################################################
!PERFORM ELECTRONIC ITERATIONS
!######################################################################

!Offer option of computing 2nd-order total energy with existing
!wavefunctions when nstep<=0, else do nstep iterations
!Note that for non-self-consistent calculations, this loop will be exited
!after the first call to vtorho3

!Pass through the first routines even when nstep==0
 write(std_out,*) 'scfcv3, nstep=', max(1,nstep)
 do istep=1,max(1,nstep)

!  ######################################################################
!  The following steps are done once
!  ----------------------------------------------------------------------
   if (istep==1)then

!    PAW only: compute frozen part of 1st-order compensation density
!    and frozen part of psp strengths Dij
!    ----------------------------------------------------------------------
     if (psps%usepaw==1) then
       optfr=1;if (iscf_mod>=0.or.usexcnhat==0) optfr=optfr+dtset%pawstgylm
       if (usexcnhat==0) then
         ABI_ALLOCATE(vtmp,(nfftf,nspden))
         vtmp=vtrial-vxc
       else
         vtmp => vtrial
       end if
       call pawfrnhat(cplex,gprimd,idir,ipert,mpi_enreg,dtset%natom,dtset%natom,nfftf,ngfftf,nspden,&
&       psps%ntypat,optfr,paw_ij1,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,qphon,&
&       rprimd,ucvol,vpsp1,vtmp,xred)
       if (usexcnhat==0)  then
         ABI_DEALLOCATE(vtmp)
       end if
       nullify(vtmp)
     end if

!    PAW only: we sometimes have to compute 1st-order compensation density
!    and eventually add it to density from 1st-order WFs
!    ----------------------------------------------------------------------
     if (psps%usepaw==1.and.(ipert/=dtset%natom+1.and.ipert/=dtset%natom+5).and. &
&     ((usexcnhat==0) &
&     .or.(dtfil%ireadwf/=0.and.dtset%get1den==0.and.dtset%ird1den==0.and.&
&     iscf_mod>0.and.initialized==0))) then
       call timab(564,1,tsec)
       ider=0;izero=0
       call pawmknhat(dum,cplex,ider,idir,ipert,izero,gprimd,mpi_enreg,dtset%natom,dtset%natom,&
&       nfftf,ngfftf,ider,nspden,psps%ntypat,dtset%paral_kgb,pawang,pawfgrtab,&
&       nhatgr_dum,nhat1,pawrhoij1,pawrhoij,pawtab,qphon,rprimd,ucvol,xred)
       if (dtfil%ireadwf/=0.and.dtset%get1den==0.and.dtset%ird1den==0.and.initialized==0) then
         rhor1(:,:)=rhor1(:,:)+nhat1(:,:)
         call fourdp(cplex,rhog1,rhor1(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
       end if
       call timab(564,2,tsec)
     end if

!    Set initial guess for 1st-order potential
!    ----------------------------------------------------------------------
     call status(istep,dtfil%filstat,iexit,level,'get vtrial1   ')
     optene=-1;option=1
     call rhotov3(cplex,ehart01,ehart1,elpsp1,exc1,gmet,gprimd,gsqcut,idir,ipert,&
&     kxc,mpi_enreg,dtset%natom,nfftf,ngfftf,nhat1,nkxc,nspden,n3xccc,&
&     optene,option,dtset%paral_kgb,dtset%qptn,rhog,rhog1,rhor,rhor1,&
&     rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,nvresid1,res2,vtrial1,xccc3d1)

!    For Q=0 and metallic occupation, initialize quantities needed to
!    compute the first-order Fermi energy
!    ----------------------------------------------------------------------
     if(qzero==1 .and. (dtset%occopt>=3 .and. dtset%occopt <=8) .and.&
&     (ipert<=dtset%natom .or. ipert==dtset%natom+3 .or. ipert==dtset%natom+4) .and. dtset%frzfermi==0) then
       ABI_ALLOCATE(rhorfermi,(cplex*nfftf,nspden))
       if (psps%usepaw==1)then ! This is temporary !!!!
         rhorfermi=zero;fe1fixed=zero;eigen1=zero
       else
         call rhofermi3(atindx,atindx1,cg,cgq,cplex,&
&         doccde_rbz,docckqde,dtfil,dtset,edocc,eeig0,eigenq,eigen0,eigen1,&
&         fe1fixed,gmet,gprimd,idir,&
&         ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,dtset%mband,mgfftf,&
&         mkmem,mkqmem,mk1mem,mpi_enreg,mpsang,mpw,mpw1,&
&         dtset%natom,nattyp,nband_rbz,nfftf,nkpt_rbz,npwarr,npwar1,nspden,&
&         dtset%nsppol,nsym1,dtset%ntypat,occkq,occ_rbz,phnons1,&
&         ph1d,dtset%prtvol,psps,rhorfermi,rmet,rprimd,symaf1,symrl1,ucvol,&
&         wfftgs,wfftkq,wtk_rbz,xred,ylm,ylm1,ylmgr1)
       end if !First-order fermi energy setup
     end if
!    End the condition of istep==1
   end if

!  ######################################################################
!  The following steps are done at every iteration
!  ----------------------------------------------------------------------

   if (psps%usepaw==1)then
!    Computation of "on-site" 2nd-order energy, first-order potentials, first-order densities
     nzlmopt=0;if (istep==2.and.dtset%pawnzlm>0) nzlmopt=-1
     if (istep>2) nzlmopt=dtset%pawnzlm
     option=0;if (dtset%iscf>0.and.dtset%iscf<10.and.nstep>0) option=1
     if (ipert/=dtset%natom+1.and.ipert/=dtset%natom+5) then
       do iatom=1,dtset%natom
         itypat=dtset%typat(iatom)
         v_size=paw_an1(iatom)%lm_size;if (dtset%pawxcdev==0) v_size=paw_an1(iatom)%angl_size
         paw_ij1(iatom)%has_dij=1;paw_ij1(iatom)%has_dijhartree=1
         paw_an1(iatom)%has_vxc=1
         ABI_ALLOCATE(paw_ij1(iatom)%dijhartree,(cplex*pawtab(itypat)%lmn2_size))
         paw_ij1(iatom)%dijhartree(:)=zero
         ABI_ALLOCATE(paw_an1(iatom)%vxc1 ,(cplex*pawtab(itypat)%mesh_size,v_size,paw_an1(iatom)%nspden))
         paw_an1(iatom)%vxc1(:,:,:)=zero
         ABI_ALLOCATE(paw_an1(iatom)%vxct1,(cplex*pawtab(itypat)%mesh_size,v_size,paw_an1(iatom)%nspden))
         paw_an1(iatom)%vxct1(:,:,:)=zero
       end do
     end if
     call status(istep,dtfil%filstat,iexit,level,'call pawdenpot')
     call pawdenpot(dum,epaw1,epawdc1_dum,ipert,dtset%ixc,mpi_enreg,dtset%natom,dtset%natom,dtset%nspden,&
&     psps%ntypat,nzlmopt,option,dtset%paral_kgb,paw_an1,paw_an,paw_ij1,pawang,dtset%pawprtvol,&
&     pawrad,pawrhoij1,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%spnorbscl,&
&     dtset%xclevel,dtset%xc_denpos,psps%znuclpsp)
!    First-order Dij computation
     call status(istep,dtfil%filstat,iexit,level,'call pawdij   ')
!    vpsp1 contribution to Dij already stored in forzen part of Dij
     if (ipert/=dtset%natom+1.and.ipert/=dtset%natom+5) then
       if (paw_ij1(1)%has_dijfr==2) then
         do ispden=1,min(dtset%nspden,2)
           vtrial1(:,ispden)=vtrial1(:,ispden)-vpsp1(:)
         end do
       end if
       if (usexcnhat/=0) then
         vxc1 => vtrial1   ! Unused
       else
         if (paw_ij1(1)%has_dijfr==2) then
           ABI_ALLOCATE(vxc1,(cplex*nfftf,dtset%nspden))
           do ispden=1,min(dtset%nspden,2)
             vxc1(:,ispden)=vtrial1(:,ispden)-vhartr1(:)
           end do
         else
           do ispden=1,min(dtset%nspden,2)
             vxc1(:,ispden)=vtrial1(:,ispden)-vhartr1(:)-vpsp1(:)
           end do
         end if
       end if
     else
       vxc1 => vtrial1   ! Unused
     end if
     call pawdij(cplex,dtset,dtset%enunit,one,gprimd,ipert,mpi_enreg,dtset%natom,dtset%natom,nfftf,ngfftf,&
&     dtset%nspden,psps%ntypat,dtset%paral_kgb,paw_an1,paw_ij1,pawang,pawfgrtab,dtset%pawprtvol,pawrad,&
&     dtset%pawspnorb,pawtab,dtset%pawxcdev,qphon,dtset%typat,ucvol,vtrial1,vxc1,xred)

     if (ipert/=dtset%natom+1.and.ipert/=dtset%natom+5) then
       if (paw_ij1(1)%has_dijfr==2) then
         do ispden=1,min(dtset%nspden,2)
           vtrial1(:,ispden)=vtrial1(:,ispden)+vpsp1(:)
         end do
       end if
       if (usexcnhat==0)  then
         ABI_DEALLOCATE(vxc1)
       end if
       nullify(vxc1)
       do iatom=1,dtset%natom
         ABI_DEALLOCATE(paw_ij1(iatom)%dijhartree)
         ABI_DEALLOCATE(paw_an1(iatom)%vxc1)
         ABI_DEALLOCATE(paw_an1(iatom)%vxct1)
         paw_an1(iatom)%has_vxc=0;paw_ij1(iatom)%has_dijhartree=0
       end do
     end if
     call status(istep,dtfil%filstat,iexit,level,'call symdij   ')
     call symdij(gprimd,psps%indlmn,indsy1,ipert,psps%lmnmax,dtset%natom,nsym1,psps%ntypat,0,&
&     paw_ij1,pawang1,dtset%pawprtvol,rprimd,symaf1,symrc1,dtset%typat)
     if (ipert/=dtset%natom+1.and.ipert/=dtset%natom+5) then
       if (paw_ij(1)%has_dijhat>0) &
&       call symdij(gprimd,psps%indlmn,indsy1,ipert,psps%lmnmax,dtset%natom,nsym1,psps%ntypat,1,&
&       paw_ij1,pawang1,dtset%pawprtvol,rprimd,symaf1,symrc1,dtset%typat)
     end if
   end if ! end usepaw section

!  No need to continue and call vtorho3, when nstep==0
   if(nstep==0)exit

!  ######################################################################
!  The following steps are done only when nstep>0
!  ----------------------------------------------------------------------
   call status(istep,dtfil%filstat,iexit,level,'loop istep    ')

   if(iscf_mod>0)then
     write(msg, '(a,a,i4)' )ch10,' ITER STEP NUMBER  ',istep
     call wrtout(std_out,msg,'COLL')
   end if

!  For Q=0 and metallic occupation, calculate the first-order Fermi energy
   if(qzero==1 .and. (dtset%occopt>=3 .and. dtset%occopt <=8) .and.&
&   (ipert<=dtset%natom .or. ipert==dtset%natom+3 .or. ipert==dtset%natom+4) .and. dtset%frzfermi==0) then
     nfftot=ngfftf(1)*ngfftf(2)*ngfftf(3)
!    MT Fev. 27 2009: this routine has to be checked for PAW !
     if (psps%usepaw==1) then ! This is temporary !!!
       fermie1=zero
     else
       call newfermie1(cplex,fermie1,fe1fixed,istep,&
&       mpi_enreg,nfftf,nfftot,nspden,dtset%occopt,&
&       dtset%prtvol,rhorfermi,ucvol,vtrial1)
     end if
   end if

!  DEBUG
!  write(std_out,*)' scfcv3 : before vtorho3, vtrial1(1,1)=',vtrial1(1,1)
!  ENDDEBUG

!  ######################################################################
!  Compute the 1st-order density rho1 from the 1st-order trial potential
!  ----------------------------------------------------------------------
   call status(istep,dtfil%filstat,iexit,level,'call vtorho3  ')
   call vtorho3(atindx,atindx1,cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cprj1,cpus,dbl_nnsclo,&
&   dimcprj,dim_eig2rf,doccde_rbz,docckqde,dtefield,&
&   dtfil,dtset,edocc,eeig0,eigenq,eigen0,eigen1,ek0,ek1,eloc0,&
&   enl0,enl1,fermie1,gh0c1_set,gh1c_set,gmet,gprimd,hdr,idir,indsy1,&
&   ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,dtset%mband,&
&   mkmem,mkqmem,mk1mem,mpi_enreg,mpsang,mpw,mpw1,&
&   dtset%natom,nband_rbz,ncpgr,nfftf,nhat1,nkpt_rbz,npwarr,npwar1,res2,nspden,&
&   dtset%nsppol,nsym1,dtset%ntypat,nvresid1,occkq,occ_rbz,optres,&
&   paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrhoij,pawrhoij1,pawtab,&
&   phnons1,ph1d,dtset%prtvol,psps,pwindall,qmat,resid,residm,rhog1,rhor1,rmet,&
&   rprimd,symaf1,symrc1,symrl1,ucvol,usecprj,useylmgr1,&
&   wffddk,wffnew,wffnow,wfftgs,wfftkq,vtrial,vtrial1,wtk_rbz,xred,ylm,ylm1,ylmgr1)
   call status(istep,dtfil%filstat,iexit,level,'after vtorho3 ')

   if (dtset%berryopt == 4) then
!    calculate \Omega E \cdot P term
     if (ipert<=dtset%natom) then
!      phonon perturbation
       call  ebp3(cg,cg1,dtefield,eberry,dtset%mband,mkmem,&
&       mpw,mpw1,nkpt,npwarr,npwar1,dtset%nsppol,dtset%nspinor,pwindall,qmat)
     else if (ipert==dtset%natom+2) then
!      electric field perturbation
       call  edie3(cg,cg1,dtefield,eberry,idir,dtset%mband,mkmem,&
&       mpw,mpw1,nkpt,npwarr,npwar1,dtset%nsppol,dtset%nspinor,pwindall,qmat,rprimd)
     end if
   end if

!  ######################################################################
!  Skip out of step loop if non-SCF (completed)
!  ----------------------------------------------------------------------

!  Indeed, nstep loops have been done inside vtorho3
   if (iscf_mod<=0 .and. iscf_mod/=-3) exit

!  ######################################################################
!  In case of density mixing , compute the total 2nd-order energy,
!  check the exit criterion,
!  then mix the 1st-order density
!  ----------------------------------------------------------------------

   if (iscf_mod>=10) then
     optene = 1 ! use double counting scheme
     call etot3(dtset%berryopt,deltae,eberry,edocc,eeig0,eew,efrhar,efrkin,&
&     efrloc,efrnl,efrx1,efrx2,ehart1,ek0,ek1,eii,elast,eloc0,elpsp1,&
&     enl0,enl1,epaw1,etotal,evar,exc1,ipert,dtset%natom,optene)

     call timab(152,1,tsec)
     choice=2
     call status(istep,dtfil%filstat,iexit,level,'print info    ')
     call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
&     etotal,favg,fcart,fermie,dtfil%fnametmp_eig,dtfil%filnam_ds(1),&
&     1,iscf_mod,istep,kpt_rbz,maxfor,&
&     mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
&     nstep,occ_rbz,0,prtfor,0,&
&     quit,res2,resid,residm,response,&
&     tollist,psps%usepaw,vxcavg,wtk_rbz,xred)
     call timab(152,2,tsec)
     if (istep==nstep) quit=1
!    If criteria in scprqt say to quit, then exit the loop over istep
     quit_sum=quit
     if (mpi_enreg%paral_compil_kpt==1) then
       call xcomm_world(mpi_enreg,spaceComm)
       call xsum_mpi(quit_sum,spaceComm,ierr)
     end if
     if (quit_sum>0) exit
     call status(istep,dtfil%filstat,iexit,level,'call newrho   ')
!    INSERT HERE CALL TO NEWRHO3 : to be implemented
     if (psps%usepaw==1) stop " newrho3 not implemented: use potential mixing !"
     initialized=1
   end if

!  ######################################################################
!  Compute the new 1st-order potential from the 1st-order density
!  ----------------------------------------------------------------------

   call status(istep,dtfil%filstat,iexit,level,'call rhotov3   ')
   optene=2*optres; !if(psps%usepaw==1) optene=2
   call rhotov3(cplex,ehart01,ehart1,elpsp1,exc1,gmet,gprimd,gsqcut,idir,ipert,&
&   kxc,mpi_enreg,dtset%natom,nfftf,ngfftf,nhat1,nkxc,nspden,n3xccc,&
&   optene,optres,dtset%paral_kgb,dtset%qptn,rhog,rhog1,rhor,rhor1,&
&   rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,nvresid1,res2,vtrial1,xccc3d1)

!  ######################################################################
!  In case of potential mixing , compute the total 2nd-order energy,
!  check the exit criterion, then mix the 1st-order potential
!  ----------------------------------------------------------------------

   if (iscf_mod<10) then

!    PAW: has to compute here the "on-site" 2nd-order energy
     if (psps%usepaw==1) then
       nzlmopt=0;if (istep==1.and.dtset%pawnzlm>0) nzlmopt=-1
       if (istep>1) nzlmopt=dtset%pawnzlm
       option=2
       if (ipert/=dtset%natom+1.and.ipert/=dtset%natom+5) then
         do iatom=1,dtset%natom
           ABI_ALLOCATE(paw_ij1(iatom)%dijhartree,(cplex*pawtab(dtset%typat(iatom))%lmn2_size))
           paw_ij1(iatom)%has_dijhartree=1
         end do
       end if
       call pawdenpot(dum,epaw1,epawdc1_dum,ipert,dtset%ixc,mpi_enreg,dtset%natom,dtset%natom,dtset%nspden,&
&       psps%ntypat,nzlmopt,option,dtset%paral_kgb,paw_an1,paw_an,paw_ij1,pawang,dtset%pawprtvol,&
&       pawrad,pawrhoij1,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%spnorbscl,&
&       dtset%xclevel,dtset%xc_denpos,psps%znuclpsp)
       if (ipert/=dtset%natom+1.and.ipert/=dtset%natom+5) then
         do iatom=1,dtset%natom
           ABI_DEALLOCATE(paw_ij1(iatom)%dijhartree)
           paw_ij1(iatom)%has_dijhartree=0
         end do
       end if
     end if

     optene = 0 ! use direct scheme
     call etot3(dtset%berryopt,deltae,eberry,edocc,eeig0,eew,efrhar,efrkin,&
&     efrloc,efrnl,efrx1,efrx2,ehart1,ek0,ek1,eii,elast,eloc0,elpsp1,&
&     enl0,enl1,epaw1,etotal,evar,exc1,ipert,dtset%natom,optene)

     call timab(152,1,tsec)
     choice=2
     call status(istep,dtfil%filstat,iexit,level,'print info    ')
     call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
&     etotal,favg,fcart,fermie,dtfil%fnametmp_eig,dtfil%filnam_ds(1),&
&     1,iscf_mod,istep,kpt_rbz,maxfor,&
&     mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
&     nstep,occ_rbz,0,prtfor,0,&
&     quit,res2,resid,residm,response,&
&     tollist,psps%usepaw,vxcavg,wtk_rbz,xred)
     call timab(152,2,tsec)
!    If criteria in scprqt say to quit, then exit the loop over istep
     quit_sum=quit
     if (mpi_enreg%paral_compil_kpt==1) then
       call xcomm_world(mpi_enreg,spaceComm)
       call xsum_mpi(quit_sum,spaceComm,ierr)
     end if
     if (quit_sum>0) exit
     if(iscf_mod/=-3)then
!      Note that nvresid1 and vtrial1 are called vresid and vtrial inside this routine
       call newvtr3(cplex,dbl_nnsclo,dielar,dtset,etotal,pawfgr%fintocoa,&
&       initialized,iscf_mod,ispmix,istep,mix,pawfgr%coatofin,&
&       mpi_enreg,nfftf,nfftmix,ngfftf,ngfftmix,npawmix,pawrhoij1,&
&       qphon,rhor1,rprimd,psps%usepaw,nvresid1,vtrial1)
       initialized=1
     end if
   end if

!  ######################################################################
!  END MINIMIZATION ITERATIONS
!  ######################################################################

!  Note that there are different "exit" instructions within the loop
 end do ! istep

 if (iscf_mod>0.or.iscf_mod==-3)  then
   ABI_DEALLOCATE(nvresid1)
 end if

!DEBUG
!write(std_out,*)' scfcv3 : after istep loop, continue'
!stop
!ENDDEBUG

!if(nstep==0) then (not implemented)
!end if


!######################################################################
!Additional steps after SC iterations
!----------------------------------------------------------------------

 call timab(160,1,tsec)
 call status(0,dtfil%filstat,iexit,level,'endloop istep ')

!Delete eventual _FFT file
 if(dtset%mffmem==0)then
   inquire (file=dtfil%fnametmp_fft,exist=ex)
   if(ex)then
     open(unit=tmp_unit,file=dtfil%fnametmp_fft,form='unformatted',status='old')
     close(unit=tmp_unit,status='DELETE')
   end if
 end if

!Eventually close the dot file, before calling nstdy3
 if(ipert==dtset%natom+2 .and. sum( (dtset%qptn(1:3)) **2 ) <= 1.0d-7 .and. (dtset%berryopt .ne. 4) )then
   call WffClose(wffddk,ierr)
 end if

!Deallocate the no more needed arrays
 if (iscf_mod>0) then
   call ab6_mixing_deallocate(mix)
 end if
 if( (nstep>0 .and. iscf_mod>0) .or. iscf_mod==-1 ) then
   ABI_DEALLOCATE(dielinv)
   ABI_DEALLOCATE(susmat)
 end if
 if(allocated(rhorfermi))  then
   ABI_DEALLOCATE(rhorfermi)
 end if
 if(psps%usepaw==1) then
   if (mk1mem/=0.and.usecprj==1) then
     call cprj_free(cprj1)
   end if
   ABI_DEALLOCATE(cprj1)
   do iatom=1,dtset%natom
     if (pawfgrtab(iatom)%nhatfr_allocated>0)  then
       ABI_DEALLOCATE(pawfgrtab(iatom)%nhatfr)
     end if
     pawfgrtab(iatom)%nhatfr_allocated=0
   end do
   if (nstep>0.and.iscf_mod>0) then
     do iatom=1,dtset%natom
       pawrhoij1(iatom)%lmnmix_sz=0
       pawrhoij1(iatom)%use_rhoijres=0
       ABI_DEALLOCATE(pawrhoij1(iatom)%kpawmix)
       ABI_DEALLOCATE(pawrhoij1(iatom)%rhoijres)
     end do
   end if
 end if ! PAW

 call timab(160,2,tsec)
 call timab(150,1,tsec)

 if(ipert==dtset%natom+3 .or. ipert==dtset%natom+4) then
   call status(0,dtfil%filstat,iexit,level,'enter nselt3  ')
   call nselt3(atindx,atindx1,blkflg,cg,cg1,cplex,&
&   d2bbb,d2lo,d2nl,ecut,dtset%ecutsm,dtset%effmass,&
&   gmet,gprimd,gsqcut,idir,&
&   ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,dtset%mband,mgfftf,&
&   mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,mpw1,&
&   dtset%natom,nattyp,nband_rbz,nfftf,ngfftf,&
&   nkpt_rbz,nkxc,dtset%nloalg,&
&   npwarr,npwar1,nspden,dtset%nspinor,dtset%nsppol,&
&   nsym1,dtset%ntypat,occ_rbz,&
&   dtset%paral_kgb,ph1d,dtset%prtbbb,psps,dtset%qptn,rhog,&
&   rhor,rhor1,rmet,rprimd,symrc1,dtset%typat,ucvol,&
&   dtfil%unkg,dtfil%unkg1,&
&   wffnow,wfftgs,dtfil%unylm,dtfil%unylm1,&
&   wtk_rbz,xred,ylm,ylm1,ylmgr,ylmgr1)
 end if
 if(ipert<=dtset%natom+5)then
   if (psps%usepaw==1) then
     call status(0,dtfil%filstat,iexit,level,'enter pawnst3 ')
     call nstpaw3(blkflg,cg,cgq,cg1,cplex,cprj,cprjq,dimcprj,docckqde,doccde_rbz,dtfil,dtset,&
&     d2lo,d2nl,d2ovl,eigenq,eigen0,eigen1,&
&     eovl1,gmet,gprimd,gsqcut,idir,indkpt1,indsy1,ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,&
&     kxc,mgfftf,mpert,mpi_enreg,mpw,mpw1,nband_rbz,ncpgr,nfftf,ngfftf,nhat1,&
&     nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,dtset%nspinor,dtset%nsppol,nsym1,n3xccc,occkq,occ_rbz,&
&     paw_an,paw_an1,paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoij,pawrhoij1,&
&     pawtab,phnons1,ph1d,ph1df,psps,rhor1,rmet,rprimd,symaf1,symrc1,symrl1,ucvol,usecprj,&
&     usexcnhat,useylmgr1,vhartr1,vpsp1,vtrial,vtrial1,vxc,wffnow,wfftgs,wfftkq,wtk_rbz,&
&     xccc3d1,xred,ylm,ylm1,ylmgr1)
   else
     call status(0,dtfil%filstat,iexit,level,'enter nstdy3  ')
     call nstdy3(atindx,blkflg,cg,cg1,cplex,dtfil,dtset,d2bbb,d2lo,d2nl,eigen0,eigen1,gmet,&
&     gsqcut,idir,indkpt1,indsy1,ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mpert,mpi_enreg,&
&     mpw,mpw1,nattyp,nband_rbz,nfftf,ngfftf,nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,&
&     dtset%nsppol,nsym1,occ_rbz,ph1d,psps,rhor1,rmet,rprimd,symrc1,ucvol,&
&     wffnow,wfftgs,wtk_rbz,xred,ylm,ylm1)
   end if
 end if

 call timab(150,2,tsec)
 call timab(160,1,tsec)

!calculate Born effective charge and store it in d2lo
 if (dtset%berryopt == 4 .and. ipert <= dtset%natom) then
   call bec3(cg,cg1,dtefield,dtset%natom,d2lo,idir,ipert,dtset%mband,mkmem,&
&   mpw,mpw1,mpert,nkpt,npwarr,npwar1,dtset%nsppol,dtset%nspinor,pwindall,qmat,rprimd)
   blkflg(:,dtset%natom+2,:,1:dtset%natom)=1
 end if

!calculate dielectric tensor and store it in d2lo
 if (dtset%berryopt == 4 .and. ipert == dtset%natom+2 ) then
   call die3(cg,cg1,dtefield,d2lo,idir,ipert,dtset%mband,mkmem,&
&   mpw,mpw1,mpert,nkpt,npwarr,npwar1,dtset%nsppol,dtset%nspinor,pwindall,qmat,rprimd)
   blkflg(:,dtset%natom+2,:,dtset%natom+2)=1
 end if

!If SCF convergence was not reached (for nstep>0),
!print a warning to the output file (non-dummy arguments: nstep,
!residm, diffor - infos from tollist have been saved inside )
 call status(0,dtfil%filstat,iexit,level,'call scprqt(en')
 choice=3
 call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
& etotal,favg,fcart,fermie,dtfil%fnametmp_eig,dtfil%filnam_ds(1),&
& 1,iscf_mod,istep,kpt_rbz,maxfor,&
& mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
& nstep,occ_rbz,0,prtfor,0,&
& quit,res2,resid,residm,response,&
& tollist,psps%usepaw,vxcavg,wtk_rbz,xred)

!Optionally provide output of charge density and/or potential in real space,
!as well as analysis of geometrical factors (bond lengths and bond angles).
!Warnings :
!- core charge is excluded from the charge density;
!- the potential is the INPUT vtrial.
 if(me==0) then
   if (dtset%prtden>0) then
     call status(0,dtfil%filstat,iexit,level,'call ioarr-den')
     rdwr=2 ; fformr=52 ; rdwrpaw=0
     call appdig(pertcase,dtfil%fnameabo_den,fi1o)
     call ioarr(accessfil,rhor1, dtset, etotal,fformr,fi1o,hdr, mpi_enreg, &
&     cplex*nfftf,pawrhoij1,rdwr,rdwrpaw,wvl)
   end if
   if (dtset%prtpot>0) then
     call status(0,dtfil%filstat,iexit,level,'call ioarr-pot')
     rdwr=2 ; fformv=102 ; rdwrpaw=0
     call appdig(pertcase,dtfil%fnameabo_pot,fi1o)
     call ioarr(accessfil,vtrial1, dtset, fermie,fformv,fi1o,hdr, mpi_enreg, &
&     cplex*nfftf,pawrhoij1,rdwr,rdwrpaw,wvl)
   end if
 end if

!All procs waiting here...
 if(mpi_enreg%paral_compil_kpt==1 .or. mpi_enreg%paral_compil_fft==1)then
   call timab(61,1,tsec)
   call leave_test()
   call timab(61,2,tsec)
 end if

!Debugging : print the different parts of rhor1
!MPIWF Warning : this should not be parallelized over space, leave this debugging feature as such.
 if(dtset%prtvol==-level)then
   write(msg,'(a)') '   ir       rhor1(ir)    '
   call wrtout(std_out,msg,'COLL')
   do ir=1,nfftf
     if(ir<=11 .or. mod(ir,301)==0 )then
       write(msg,'(i5,a,es13.6)')ir,' ',rhor1(ir,1)
       call wrtout(std_out,msg,'COLL')
       if(nspden==2)then
         write(msg,'(a,es13.6)')'      ',rhor1(ir,2)
         call wrtout(std_out,msg,'COLL')
       end if
     end if
   end do
 end if

!Deallocate arrays
 ABI_DEALLOCATE(fcart)
 ABI_DEALLOCATE(vtrial1)
 ABI_DEALLOCATE(vhartr1)
 if (dtset%berryopt == 4) then
   ABI_DEALLOCATE(pwindall)
   ABI_DEALLOCATE(qmat)
   ABI_DEALLOCATE(dtefield%ikpt_dk)
   ABI_DEALLOCATE(dtefield%cgindex)
   ABI_DEALLOCATE(dtefield%idxkstr)
   ABI_DEALLOCATE(dtefield%kgindex)
   if(associated(dtefield%fkgindex))  then
     ABI_DEALLOCATE(dtefield%fkgindex)
   end if
   if(associated(mpi_enreg%kpt_loc2ibz_sp))  then
     ABI_DEALLOCATE(mpi_enreg%kpt_loc2ibz_sp)
   end if
 end if
 if(psps%usepaw==1) then
   ABI_DEALLOCATE(nhat1)
   call destroy_paw_an(paw_an1)
   ABI_DEALLOCATE(paw_an1)
   call destroy_paw_ij(paw_ij1)
   ABI_DEALLOCATE(paw_ij1)
 end if

!Structured debugging : if dtset%prtvol=-level, stop here.
 if(dtset%prtvol==-level)then
   write(msg,'(a1,a,a1,a,i2,a)') ch10,' scfcv3: exit ',&
&   ch10,'  prtvol=-',level,', debugging mode => stop '
   call wrtout(std_out,msg,'COLL')
   call leave_new('COLL')
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(160,2,tsec)
 call timab(120,2,tsec)

 DBG_EXIT("COLL")

end subroutine scfcv3
!!***
