!{\src2tex{textfont=tt}}
!!****f* ABINIT/vtorho
!! NAME
!! vtorho
!!
!! FUNCTION
!! This routine compute the new density from a fixed potential (vtrial)
!! but might also simply compute eigenvectors and eigenvalues.
!! The main part of it is a wf update over all k points.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MF, AR, MM, MT, FJ, MB, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  afford=used to dimension susmat
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cpus= cpu time limit in seconds
!!  dbl_nnsclo=if 1, will double the value of dtset%nnsclo
!!  dielop= if positive, the dielectric matrix must be computed.
!!  dielstrt=number of the step at which the dielectric preconditioning begins.
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem =number of k points which can fit in memory; set to 0 if use disk
!!   | mpw=maximum dimensioned size of npw
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | nkpt=number of k points.
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!   | typat= array of types of the natoms
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  etotal=total energy (Ha) - only needed for tddft
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for the dielectric matrix
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!   (3x3 tensor) and grads wrt atomic coordinates (3*natom)
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  irrzondiel(nfftdiel**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data for diel matrix
!!                                     nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
!!  istep=index of the number of steps in the routine scfcv
!!  istep_mix=index of the number of steps for the SCF mixing (can be <istep)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kxc(nfftf,nkxc)=exchange-correlation kernel, needed only if nkxc/=0 .
!!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of the dielectric matrix
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  nfftdiel=number of fft grid points for the computation of the diel matrix
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!                see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  npwdiel=size of the susmat array.
!!  ntypat=number of types of atoms in unit cell.
!!  optforces=option for the computation of forces (0: no force;1: forces)
!!  optres=0: the new value of the density is computed in place of the input value
!!         1: only the density residual is computed ; the input density is kept
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!                                    nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
!!  phnonsdiel(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases,
!!   for diel matr
!!                                     nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  ph1ddiel(2,3*(2*mgfftdiel+1)*natom*usepaw)=one-dimensional structure factor information
!!                                             for the dielectric matrix
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive vectors
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  ucvol=unit cell volume in bohr**3.
!!  wffnew,wffnow=unit numbers for wf disk files.
!!  vtrial(nfftf,nspden)=INPUT potential Vtrial(r).
!!  vxctau=(only for meta-GGA): derivative of XC energy density with respect to
!!    kinetic energy density (depsxcdtau). The arrays vxctau(nfft,nspden,4) contains also
!!    the gradient of vxctau (gvxctau) in vxctau(:,:,2:4)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!  ylmdiel(npwdiel,lmax_diel**2)= real spherical harmonics for each G and k point
!!                                 for the dielectric matrix
!!
!! OUTPUT
!!  compch_fft=-PAW only- compensation charge inside spheres computed over fine fft grid
!!  dphase(3) : dphase(idir) = accumulated change in the string-averaged
!!     Zak phase along the idir-th direction caused by the update of all
!!     the occupied Bloch states at all the k-points (only if finite electric field)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points.
!!  residm=maximum value from resid array (except for nbdbuf highest bands)
!!  susmat(2,npwdiel*afford,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!  === if optforces>0 ===
!!    grnl(3*natom)=stores grads of nonlocal energy wrt length scales
!!  ==== if optres==1
!!    nres2=square of the norm of the residual
!!    nvresid(nfftf,nspden)=density residual
!!  ==== if psps%usepaw==1
!!    nhat(nfftf,nspden*psps%usepaw)=compensation charge density on rectangular grid in real space
!!
!! SIDE EFFECTS
!!  cg(2,mpw*dtset%nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!   At output contains updated wavefunctions coefficients;
!!    if nkpt>1, these are kept in a disk file.
!!  energies <type(energies_type)>=storage for energies computed here :
!!   | e_eigenvalues=Sum of the eigenvalues - Band energy (Hartree)
!!   | e_kinetic=kinetic energy part of total energy
!!   | e_nonlocalpsp=nonlocal pseudopotential part of total energy
!!   | e_fermie=fermi energy (Hartree)
!!  occ(mband*nkpt*nsppol)=occupation number for each band for each k.
!!      (input if insulator - occopt<3 - ; output if metallic)
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  rhog(2,nfftf)=Fourier transform of total electron density
!!  rhor(nfftf,nspden)=total electron density (el/bohr**3)
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  wvl <type(wvl_data)>=wavelets structures in case of wavelets basis.
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      calcdensph,clsopn,cprj_alloc,cprj_diskinit_w,cprj_free,cprj_get
!!      cprj_mpi_allgather,cprj_put,ctocprj,datafordmft,destroy_crystal
!!      destroy_dmft,destroy_hamiltonian,destroy_oper,dmft_solve,fftpac
!!      gpu_finalize_ffnl_ph3d,gpu_finalize_ham_data,gpu_update_ffnl_ph3d
!!      gpu_update_ham_data,hdr_io,hdr_io_netcdf,hdr_skip,hdr_update
!!      ini_wf_netcdf,init_crystal,init_dmft,init_hamiltonian,init_oper
!!      leave_new,leave_test,mag_loc_k,magcart,mkffnl,mkkin,mkkpg,mkrho,newocc
!!      pawmkrho,pawmkrhoij,ph1d3d,prep_bandfft_tabs,print_dmft,prteigrs
!!      prtrhomxmn,rdnpw,rhoij_init_unpacked,rwwf,sphereboundary,sqnorm_v
!!      status,suscep_stat,symrhg,tddft,testsusmat,timab,transgrid,update_mmat
!!      vtowfk,wffkg,wrtout,wvl_nl_gradient,wvl_vtorho,xallgather_mpi
!!      xbarrier_mpi,xcomm_init,xdefineoff,xmax_mpi,xme_init,xrecv_mpi
!!      xredxcart,xsend_mpi,xsum_mpi
!!
!! NOTES
!!  Be careful to the meaning of nfft (size of FFT grids):
!!   - In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!   - In case of PAW calculations:
!!     Two FFT grids are used; one with nfft points (coarse grid) for
!!     the computation of wave functions ; one with nfftf points
!!     (fine grid) for the computation of total density.
!!
!!  The total electronic density (rhor,rhog) is divided into two terms:
!!   - The density related to WFs =Sum[Psi**2]
!!   - The compensation density (nhat) - only in PAW
!!
!!  The parallelisation needed for the electric field should be
!!  made an independent subroutine, so that this routine could be put
!!  back in the 95_drive directory.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine vtorho(afford,atindx,atindx1,cg,compch_fft,cpus,dbl_nnsclo,&
&           dielop,dielstrt,dphase,dtefield,dtfil,dtset,&
&           eigen,electronpositron,energies,etotal,gbound_diel,&
&           gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
&           istep,istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
&           mpsang,natom,nattyp,nfftf,nfftdiel,ngfftdiel,nhat,nkxc,&
&           npwarr,npwdiel,nres2,ntypat,nvresid,occ,optforces,&
&           optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
&           pwind,pwind_alloc,pwnsfac,resid,residm,rhog,rhor,&
&           rmet,rprimd,susmat,symrec,taug,taur,&
&           ucvol,wffnew,wffnow,vtrial,wvl,xred,ylm,ylmgr,ylmdiel,&
&           vxctau) ! optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_wffile

 use m_energies,           only : energies_type
 use m_hamiltonian,        only : init_hamiltonian, destroy_hamiltonian, finalize_hamiltonian
 use m_electronpositron,   only : electronpositron_type,electronpositron_calctype
 use m_paw_dmft,           only : paw_dmft_type,init_dmft,destroy_dmft,print_dmft
 use m_crystal,            only : init_crystal, destroy_crystal, crystal_structure
 use m_oper,               only : oper_type,init_oper,destroy_oper
 use m_efield

#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vtorho'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_53_spacepar
 use interfaces_56_recipspace
 use interfaces_59_io_mpi
 use interfaces_61_ionetcdf
 use interfaces_62_iowfdenpot
 use interfaces_62_occeig
 use interfaces_62_wvl_wfs
 use interfaces_65_nonlocal
 use interfaces_66_paw
 use interfaces_66_wfs
 use interfaces_67_common
 use interfaces_68_dmft
 use interfaces_77_suscep
 use interfaces_79_seqpar_mpi, except_this_one => vtorho
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments -------------------------------
 integer, intent(in) :: afford,dbl_nnsclo,dielop,dielstrt,istep,istep_mix,lmax_diel,mcg,mgfftdiel
 integer, intent(in) :: mpsang,natom,nfftf,nfftdiel,nkxc,npwdiel
 integer, intent(in) :: ntypat,optforces,optres,pwind_alloc
 real(dp), intent(in) :: cpus,etotal,gsqcut,ucvol
 real(dp), intent(out) :: compch_fft,nres2,residm
 type(MPI_type), intent(inout) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in) :: dtset
 type(efield_type), intent(inout) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(energies_type), intent(inout) :: energies
 type(hdr_type), intent(inout) :: hdr
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 type(pawang_type), intent(in) :: pawang
 type(pawfgr_type), intent(in) :: pawfgr
 type(pseudopotential_type), intent(in) :: psps
 type(wffile_type), intent(inout) :: wffnew,wffnow
 type(wvl_data), intent(inout) :: wvl
 integer, intent(in) :: atindx(natom),atindx1(natom),gbound_diel(2*mgfftdiel+8,2)
 integer, intent(in) :: indsym(4,dtset%nsym,natom)
 integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem),kg_diel(3,npwdiel),nattyp(ntypat),ngfftdiel(18),npwarr(dtset%nkpt)
 integer, intent(in) :: pwind(pwind_alloc,2,3),symrec(3,3,dtset%nsym)
 real(dp), intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*dtset%mgfft+1)*natom)
 real(dp), intent(in) :: ph1ddiel(2,(3*(2*mgfftdiel+1)*natom)*psps%usepaw)
 real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp), intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc),rmet(3,3),rprimd(3,3)
 real(dp), intent(inout) :: vtrial(nfftf,dtset%nspden)
 real(dp), intent(inout) :: xred(3,natom)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,mpsang*mpsang*psps%useylm)
 real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
 real(dp), intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol),dphase(3),grnl(3*natom)
 real(dp), intent(out) :: nhat(nfftf,dtset%nspden*psps%usepaw)
 real(dp), intent(out) :: nvresid(nfftf,dtset%nspden),resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(out) :: susmat(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden)
 real(dp), intent(inout) :: cg(2,mcg)
 real(dp), intent(inout) :: kxc(nfftf,nkxc),occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden)
 real(dp), intent(inout) :: taug(2,nfftf*dtset%usekden),taur(nfftf,dtset%nspden*dtset%usekden)
 real(dp), intent(inout),optional :: vxctau(nfftf,dtset%nspden*dtset%usekden,4)
 type(paw_ij_type),intent(in) :: paw_ij(natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type),intent(in)  :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
 integer,parameter :: level=111
 integer,save :: nwarning=0
 integer :: bantot,bdtot_index,count,count1,counter
 integer :: cplex,dest,dimdij,dimffnl,enunit
! integer :: dimenl1 ! used in magnetization
 integer :: fform,formeig,i1,i2,i3,ia,iatom,iband,iband1,ibdkpt
 integer :: ibg,icg,icg1,icg2,icp1,icp2,ider,idir
 integer :: ierr,iexit,ifft,ifor,ifor1,ii,ikg,ikg1,ikg2,ikpt
 integer :: ikptf,ikpt1f,ikpt1i
 integer :: ikpt_loc,ikpt1,ikpt_this_proc,ikxc,ilm,ilmn,index1,imagn,iorder_cprj,ipert
 integer :: iproc,ir,iscf,isp,ispden,isppol,istwf_k,itypat,jkpt
 integer :: jsppol,matblk,mband_cprj,mbdkpsp,mb2dkpsp
 integer :: mcgq,mcg_disk,mcprj,mcprj_tmp,me_distrb,mkgq,muig
 integer :: mwarning,my_nspinor,n1,n2,n3,n4,n5,n6,nband_eff
!integer :: jj,kk
 integer :: nband_k,nbuf,ndatarecv,neglect_pawhat,nfftot,nkpg,nkpt1,nnn,nnsclo_now
 integer :: nproc_distrb,npw_k,npw_k1,nsp,nvloc,option,prtvol
 integer :: rdwr,spaceComm_distrb,tag,tim_mkrho,tim_rwwf,usecprj,usetimerev
 integer :: my_source, his_source, jkptf, jkpt1f, jkpt1i
 logical :: berryflag,computesusmat,fixed_occ
 logical :: locc_test,remove_inv,with_vxctau
 real(dp) :: arg,dmft_ldaocc,dummy
!real(dp) :: sum
 real(dp) :: edmft,ebandlda,ebanddmft,ebandldatot,ekindmft,ekindmft2,ekinlda
 real(dp) :: emax,min_occ,vxcavg_dum
 character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamk
 type(wffile_type) :: wfftmp
 integer,allocatable :: dimcprj(:),dimcprj_u(:),ikpt_recv(:),kg_dum(:,:),kg_k(:,:)
 integer,allocatable :: flag_send(:,:), flag_receive(:)
 real(dp) :: dielar(7),dphase_k(3),kpoint(3),qpt(3),rhodum(1),tsec(2),ylmgr_dum(1)
 real(dp),allocatable :: EigMin(:,:),buffer(:,:),buffer1(:),buffer2(:),cgq(:,:)
 real(dp),allocatable :: cg_disk(:,:),cgrkxc(:,:),cgrvtrial(:,:),doccde(:)
 real(dp),allocatable :: dphasek(:,:),eig_dum(:),eig_k(:),ek_k(:),ek_k_nd(:,:),eknk(:),eknk_nd(:,:,:,:)
 real(dp),allocatable :: enl_k(:),enlnk(:),ffnl(:,:,:,:),grnl_k(:,:), xcart(:,:)
 real(dp),allocatable :: grnlnk(:,:),kinpw(:),kpg_k(:,:),occ_dum(:),occ_k(:),ph3d(:,:,:)
 real(dp),allocatable :: pwnsfacq(:,:),resid_k(:),rhoaug(:,:,:,:),rhowfg(:,:),rhowfr(:,:)
 real(dp),allocatable :: vlocal(:,:,:,:),vlocal_tmp(:,:,:),vxctaulocal(:,:,:,:,:),ylm_k(:,:),zshift(:)
 type(cprj_type),allocatable :: cprj(:,:),cprj_tmp(:,:),cprj_berry(:,:),cprj_gat(:,:)
 type(oper_type) :: lda_occup
 type(crystal_structure) :: cryst_struc

! *********************************************************************

!DEBUG
!write(std_out,*)' vtorho : enter'
!sum=0.0
!do jj=1,dtset%nspden
!do kk=1,nfftf
!sum=sum+abs(rhor(kk,jj))
!end do
!end do
!write(std_out,*) 'SUM(rhor)',sum
!write(std_out,*)' vtorho : enter cg2',cg(2,1:10)
!ENDDEBUG

!Keep track of total time spent in vtorho
 call timab(980,1,tsec)
 call timab(981,1,tsec)

!Check that usekden is not 0 if want to use vxctau
 with_vxctau = .false.
 if (present(vxctau) .and. dtset%usekden /= 0) with_vxctau = .true.

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
   write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' vtorho : enter '
   call wrtout(std_out,message,'COLL')
 end if

!Init MPI for kpt communicator (.i.e. communicator for that k-point)
 call xcomm_init(mpi_enreg,spaceComm_distrb)
 call xme_init(mpi_enreg,me_distrb)
 nproc_distrb=xcomm_size(spaceComm_distrb)
 if (mpi_enreg%me_img/=0) nwarning=mwarning+1

!Test size of FFT grids (1 grid in norm-conserving, 2 grids in PAW)
 if ((psps%usepaw==1.and.pawfgr%nfft/=nfftf).or.(psps%usepaw==0.and.dtset%nfft/=nfftf)) then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' vtorho :  BUG -',ch10,&
&   '  wrong values for nfft, nfftf !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Test optforces (to prevent memory overflow)
 if (optforces/=0.and.optforces/=1) then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' vtorho :  BUG -',ch10,&
&   '  wrong value for optforces !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Debugging : print vtrial and rhor
!MPIWF Warning : this should not be parallelized over space, leave this debugging feature as such.
 if(prtvol==-level)then
   if (psps%usepaw==0) then
     n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
   else
     n1=pawfgr%ngfft(1) ; n2=pawfgr%ngfft(2) ; n3=pawfgr%ngfft(3)
   end if
   write(message,'(a)') '   ir              vtrial(ir)     rhor(ir) '
   call wrtout(std_out,message,'COLL')
   do ir=1,nfftf
!    if(ir<=11 .or. mod(ir,301)==0 )then
     i3=(ir-1)/n1/n2
     i2=(ir-1-i3*n1*n2)/n1
     i1=ir-1-i3*n1*n2-i2*n1
     write(message,'(i5,3i3,a,2es13.6)')ir,i1,i2,i3,' ',vtrial(ir,1),rhor(ir,1)
     call wrtout(std_out,message,'COLL')
     if(dtset%nspden>=2)then
       write(message,'(a,2es13.6)')'               ',vtrial(ir,2),rhor(ir,2)
       call wrtout(std_out,message,'COLL')
     end if
!    end if
   end do
 end if

!WVL - Branching with a separate vtorho procedure
!in wavelet. Should be merge in vtorho later on.
 if (dtset%usewvl == 1) then
   call wvl_vtorho(dtset, energies, irrzon, istep, mpi_enreg, &
&   phnons, residm, rhor, rprimd, vtrial, wvl, xred)
   if (optforces == 1) then
     ABI_ALLOCATE(xcart,(3, dtset%natom))
     call xredxcart(dtset%natom, 1, rprimd, xcart, xred)
     call wvl_nl_gradient(grnl, mpi_enreg, dtset%natom, rprimd, wvl, xcart)
     ABI_DEALLOCATE(xcart)
   end if
   return
 end if
!WVL - Following is done in plane waves.

 iscf=dtset%iscf
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 fixed_occ=(dtset%occopt<3.or.electronpositron_calctype(electronpositron)==1)

 energies%e_eigenvalues = zero
 energies%e_kinetic     = zero
 energies%e_nonlocalpsp = zero
 grnl(:)=zero
 resid(:) = zero ! JWZ 13 May 2010. resid and eigen need to be fully zeroed each time before use
 eigen(:) = zero
 bdtot_index=0
 ibg=0;icg=0
 mbdkpsp=dtset%mband*dtset%nkpt*dtset%nsppol
 if(paw_dmft%use_dmft==1) mb2dkpsp=dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol

 ABI_ALLOCATE(eknk,(mbdkpsp))
 ABI_ALLOCATE(kg_k,(3,dtset%mpw))
 ABI_ALLOCATE(eknk_nd,(dtset%nsppol,dtset%nkpt,dtset%mband,dtset%mband*paw_dmft%use_dmft))
 ABI_ALLOCATE(EigMin,(2,dtset%mband))
 ABI_ALLOCATE(grnlnk,(3*natom,mbdkpsp*optforces))
 if (psps%usepaw==0)  then
   ABI_ALLOCATE(enlnk,(mbdkpsp))
 end if

 if(paw_dmft%use_dmft==1) eknk_nd=zero
 eknk(:)=zero;if (optforces>0) grnlnk(:,:)=zero
 if (psps%usepaw==0) enlnk(:)=zero

!Initialize rhor if needed; store old rhor
 if(iscf>0 .or. iscf==-3) then
   if (optres==1) nvresid=rhor
   if (psps%usepaw==0) then
     rhor=zero
   else
     ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
     ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
     rhowfr(:,:)=zero
   end if
 end if

!Set max number of non-self-consistent loops nnsclo_now for use in vtowfk
 if(iscf<=0)then
   nnsclo_now=dtset%nstep
 else if(iscf>0)then
   if(dtset%nnsclo>0) then
     nnsclo_now=dtset%nnsclo
   else if(dtset%nnsclo<=0)then
     nnsclo_now=1
     if(istep<=2)nnsclo_now=2
   end if
   if(dbl_nnsclo==1)then
!    DEBUG
!    write(std_out,*)' vtorho : use doubled nnsclo '
!    ENDDEBUG
     nnsclo_now=nnsclo_now*2
   end if
 end if
 if(dtset%wfoptalg==2)nnsclo_now=40  ! UNDER DEVELOPMENT

 write(message, '(a,i3,a,i3,i2,i3)' ) ' vtorho : nnsclo_now=',nnsclo_now,&
& ', note that nnsclo,dbl_nnsclo,istep=',dtset%nnsclo,dbl_nnsclo,istep
 call wrtout(std_out,message,'COLL')

 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)
 nvloc=1;if(dtset%nspden==4)nvloc=4
 ABI_ALLOCATE(rhoaug,(n4,n5,n6,nvloc))
 ABI_ALLOCATE(vlocal,(n4,n5,n6,nvloc))
 if(with_vxctau)ABI_ALLOCATE(vxctaulocal,(n4,n5,n6,nvloc,4))
 rhoaug(:,:,:,:)=zero

!Prepare wf files for reading if dtset%mkmem==0
 if (dtset%mkmem==0) then

!  Close files, and then reopen them
!  (this is supposedly helpful for use of networked workstations
!  and also sets up for later addition of a checkpoint facility
!  for restarting crashed jobs)
!  clsopn automatically checks to see whether file is scratch
!  file and if so, does not close and open it.
   call clsopn(wffnow)

!  Read wffnow header
   call hdr_skip(wffnow,ierr)

!  Define offsets, in case of MPI I/O
   formeig=0
   call xdefineOff(formeig,wffnow,mpi_enreg,dtset%nband,npwarr,dtset%nspinor,dtset%nsppol,dtset%nkpt)

   call clsopn(wffnew)

!  Update the content of the header (evolving variables)
   bantot=hdr%bantot ; dummy=1.0d20
!  WARNING! fermie is used before set.
   call hdr_update(bantot,dummy,energies%e_fermie,hdr,natom,&
&   residm,rprimd,occ,mpi_enreg,pawrhoij,psps%usepaw,xred)

!  Write the content of hdr to the new wf file
   rdwr=2 ; fform=2
   if (wffnew%accesswff /= IO_MODE_NETCDF) then
     call hdr_io(fform,hdr,rdwr,wffnew)
#if defined HAVE_TRIO_NETCDF
   else if (wffnew%accesswff == IO_MODE_NETCDF) then
     call hdr_io_netcdf(fform,hdr,rdwr,wffnew)

     call ini_wf_netcdf(dtset%mpw,wffnew%unwff,0)
#endif
   end if


!  Define offsets, in case of MPI I/O
   formeig=0
   call WffKg(wffnew,1)
   call xdefineOff(formeig,wffnew,mpi_enreg,dtset%nband,npwarr,dtset%nspinor,dtset%nsppol,dtset%nkpt)

   mcg_disk=dtset%mpw*my_nspinor*dtset%mband
   ABI_ALLOCATE(cg_disk,(2,mcg_disk))

 end if

!============================================
!==== Initialize most of the Hamiltonian ====
!============================================
!1) Allocate all arrays and initialize quantities that do not depend on k and spin.
!2) Perform the setup needed for the non-local factors:
!* Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
!* PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_hamk,psps,paw_ij,pawtab,dtset%nspinor,dtset%nspden,natom,ntypat,dtset%typat,xred,&
& dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,&
& ph1d=ph1d,electronpositron=electronpositron,use_gpu_cuda=dtset%use_gpu_cuda)

 usecprj=0; mband_cprj=0; mcprj=0
 if (psps%usepaw==1) then
   ABI_ALLOCATE(dimcprj,(natom))
   ia=0
   do itypat=1,ntypat
     dimcprj(ia+1:ia+nattyp(itypat))=pawtab(itypat)%lmn_size
     ia=ia+nattyp(itypat)
   end do
   usecprj=1
   mband_cprj=dtset%mband;if (dtset%paral_kgb/=0) mband_cprj=mband_cprj/mpi_enreg%nproc_band
   mcprj=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol*usecprj
   ABI_ALLOCATE(cprj,(natom,mcprj))
   iorder_cprj=0
   if (usecprj==1) then
     if (dtset%mkmem/=0) then
       call cprj_alloc(cprj,0,dimcprj)
     else
       call cprj_diskinit_w(atindx,natom,iorder_cprj,dtset%mkmem,natom,0,dimcprj,my_nspinor,dtfil%unpaw)
     end if
   end if
 end if

!Electric and magnetic fields: set flag to turn on various behaviors
 berryflag = .FALSE.
 if (dtset%berryopt == 4 .or. (abs(dtset%berryopt) == 5)) berryflag = .TRUE.

!allocate cprj in berryopt +/- 5 case
 if (abs(dtset%berryopt) == 5) then
   ABI_ALLOCATE(cprj_berry,(natom,dtset%mband))
   ABI_ALLOCATE(dimcprj_u,(natom))
   do iatom = 1, natom
     itypat = dtset%typat(iatom)
     dimcprj_u(iatom) = pawtab(itypat)%lmn_size
   end do
   call cprj_alloc(cprj_berry,0,dimcprj_u)
   ABI_ALLOCATE(ikpt_recv,(nproc_distrb))
   ABI_ALLOCATE(cprj_gat,(natom,nproc_distrb*dtefield%nband_occ))
   call cprj_alloc(cprj_gat,0,dimcprj_u)
 end if

!Electric field: allocate dphasek
 if (dtset%berryopt==4) then
   ABI_ALLOCATE(dphasek,(3,dtset%nkpt*dtset%nsppol))
   dphasek(:,:) = zero
   nkpt1 = dtefield%mkmem_max
 else
   nkpt1 = dtset%nkpt
 end if

 ikpt_loc = 0

 call timab(981,2,tsec)

!LOOP OVER SPINS
 do isppol=1,dtset%nsppol

   call timab(982,1,tsec)

   if (dtset%nsppol==2) then
     write(message,*)' ****  In vtorho for isppol=',isppol
     call wrtout(std_out,message,'COLL')
   end if

!  original line
!  if ((mpi_enreg%paral_compil_kpt == 0).or.(.not.berryflag)) ikpt_loc = 0
   ikpt_loc = 0

!  DEBUG feature
!  write(std_out,'(a,i4,L8,i4)')' debug--paral_compil_kpt, berryflag, ikpt_loc : ',mpi_enreg%paral_compil_kpt,&
!  &        berryflag,ikpt_loc
!  END DEBUG feature

!  Rewind kpgsph data file if needed:
   if (dtset%mkmem==0) rewind dtfil%unkg
   if (dtset%mkmem==0.and.psps%useylm==1) rewind dtfil%unylm
   ikg=0

!  Set up local potential vlocal with proper dimensioning, from vtrial
!  Also take into account the spin.
   if(dtset%nspden/=4)then
     if (psps%usepaw==0.or.pawfgr%usefinegrid==0) then
       call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vtrial,vlocal,2)
       if(with_vxctau) then
         do ii=1,4
           call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vxctau(:,:,ii),vxctaulocal(:,:,:,:,ii),2)
         end do
       end if
     else
       ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
       call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
       call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal,2)
       ABI_DEALLOCATE(cgrvtrial)
     end if
   else
     ABI_ALLOCATE(vlocal_tmp,(n4,n5,n6))
     if (psps%usepaw==0.or.pawfgr%usefinegrid==0) then
       do ispden=1,dtset%nspden
         call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vtrial,vlocal_tmp,2)
         vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
       end do
     else
       ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
       call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
       do ispden=1,dtset%nspden
         call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal_tmp,2)
         vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
       end do
       ABI_DEALLOCATE(cgrvtrial)
     end if
     ABI_DEALLOCATE(vlocal_tmp)
   end if
   rhoaug(:,:,:,:)=zero

!  PAW: retrieve Dij coefficients for this spin component
   if (psps%usepaw==1) then
     do ispden=1,dtset%nspinor**2
       isp=isppol;if (dtset%nspinor==2) isp=ispden
       do iatom=1,dtset%natom
         dimdij=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
         do ilmn=1,dimdij
           gs_hamk%ekb(ilmn,iatom,ispden)=paw_ij(iatom)%dij(ilmn,isp)
         end do
         if(dimdij+1<=gs_hamk%dimekb1) gs_hamk%ekb(dimdij+1:gs_hamk%dimekb1,iatom,ispden)=zero
       end do
     end do
   end if

!  Update enl and sij on GPU
#if defined HAVE_GPU_CUDA
   if(dtset%use_gpu_cuda==1) then
     call gpu_update_ham_data(gs_hamk%ekb,gs_hamk%sij,gprimd,gs_hamk%dimekb1,&
&     gs_hamk%dimekb2,my_nspinor,ntypat,psps%usepaw)
   end if
#endif

   call timab(982,2,tsec)

!  BIG FAT k POINT LOOP
!  MVeithen: I had to modify the structure of this loop in order to implement
!  MPI // of the electric field

!  note that the loop here differs from the similar one in berryphase_new.F90.
!  here, ikpt_loc numbers the kpts treated by the current processor.
!  in berryphase_new.F90, ikpt_loc ALSO includes info about value of isppol.

   ikpt = 0

   do while (ikpt_loc < nkpt1)

     call timab(997,1,tsec)

     if (.not.berryflag) then
       ikpt_loc = ikpt_loc + 1
       ikpt = ikpt_loc
     else
       if (ikpt_loc < dtset%mkmem) ikpt = ikpt + 1
       if ((ikpt > dtset%nkpt).and.(ikpt_loc < dtset%mkmem)) exit
     end if

     dphase_k(:) = zero
     counter=100*ikpt+isppol
     call status(counter,dtfil%filstat,iexit,level,'loop ikpt     ')
     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     istwf_k=dtset%istwfk(ikpt)
     npw_k=npwarr(ikpt)
     mcgq = 1 ; mkgq = 1

     if(mpi_enreg%paral_compil_kpt==1)then

       if (.not.berryflag) then

         if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me_distrb))/=0) then
           eigen(1+bdtot_index : nband_k+bdtot_index) = zero
           resid(1+bdtot_index : nband_k+bdtot_index) = zero
           bdtot_index=bdtot_index+nband_k
!          Skip the rest of the k-point loop
           cycle
         end if

       else      !

         if ((minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) - &
&         me_distrb)) /= 0).and.(ikpt_loc <= dtset%mkmem)) then
           eigen(1+bdtot_index : nband_k+bdtot_index) = zero
           resid(1+bdtot_index : nband_k+bdtot_index) = zero
           bdtot_index = bdtot_index + nband_k
           cycle
         end if

         mcgq = dtset%mpw*my_nspinor*nband_k*dtefield%nneigh(ikpt)
         mkgq = 6*dtset%mpw
         ikg = dtefield%kgindex(ikpt)

       end if     ! berryflag

     end if ! parallel kpt

     if (berryflag) ikpt_loc = ikpt_loc + 1
     ABI_ALLOCATE(cgq,(2,mcgq))
     ABI_ALLOCATE(pwnsfacq,(2,mkgq))

     call timab(997,2,tsec)

!    In case of MPI // of a finite field calculation
!    build the cgq array that stores the wavefunctions for the
!    neighbours of ikpt, and the pwnsfacq array that stores the
!    corresponding phase factors (in case of tnons)

     if (berryflag .and. mpi_enreg%paral_compil_kpt==1 ) then

       call timab(983,1,tsec)

       ABI_ALLOCATE(flag_send,(0:nproc_distrb-1,dtefield%fnkpt))
       ABI_ALLOCATE(flag_receive,(dtset%nkpt))
       flag_send(:,:) = 0
       flag_receive(:) = 0

       ikptf = dtefield%i2fbz(ikpt)

       do idir = 1, 3

!        skip idir values for which efield_dot(idir) = 0
         if (abs(dtefield%efield_dot(idir)) < tol12 .and. dtset%berryopt == 4) cycle

!        however, for the magnetization we need the following computations even when B = 0, because
!        cgq will be needed to construct density operator at neighboring k points

         do ifor = 1, 2

           dtefield%sflag(:,ikpt + dtset%nkpt*(isppol - 1),ifor,idir) = 0

           ikpt1f = dtefield%ikpt_dk(ikptf,ifor,idir)
           ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)
           npw_k1 = npwarr(ikpt1i)
           count = npw_k1*my_nspinor*nband_k
           my_source = mpi_enreg%proc_distrb(ikpt1i,1,isppol)

           do dest = 0, nproc_distrb-1

             if ((dest==me_distrb).and.(ikpt_loc <= dtset%mkmem)) then
!              I am dest and have something to do

               if ( my_source == me_distrb ) then
!                I am destination and source

                 ikg1 = dtefield%fkgindex(ikpt1f)
                 ikg2 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
                 pwnsfacq(:,ikg2 + 1:ikg2 + npw_k1) = pwnsfac(:,ikg1 + 1:ikg1 + npw_k1)

                 icg1 = dtefield%cgindex(ikpt1i,isppol)
                 icg2 = dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
                 cgq(:,icg2 + 1:icg2 + count) = cg(:,icg1 + 1:icg1 + count)

               else
!                I am the destination but not the source -> receive

!                receive pwnsfacq
                 ABI_ALLOCATE(buffer,(2,npw_k1))
                 tag = ikpt1f + (isppol - 1)*dtefield%fnkpt
                 call xrecv_mpi(buffer,my_source,tag,spaceComm_distrb,ierr)
                 ikg1 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
                 pwnsfacq(:,ikg1+1:ikg1+npw_k1) = buffer(:,1:npw_k1)
                 ABI_DEALLOCATE(buffer)

!                receive cgq if necessary
                 if(flag_receive(ikpt1i) == 0) then

                   ABI_ALLOCATE(buffer,(2,count))
                   tag = ikpt1i + (isppol - 1)*dtset%nkpt
                   call xrecv_mpi(buffer,my_source,tag,spaceComm_distrb,ierr)
                   icg1 = dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
                   cgq(:,icg1+1:icg1+count) = buffer(:,1:count)
                   ABI_DEALLOCATE(buffer)

                   flag_receive(ikpt1i) = 1
                 end if

               end if

             else if (ikpt_loc <= mpi_enreg%mkmem(dest)) then  ! dest != me and the dest has a k-point to treat
!              else if (dest/=me_distrb) then  

!              jkpt is the kpt which is being treated by dest (in ibz)
!              jsppol is his isppol

               jkpt = mpi_enreg%kpt_loc2ibz_sp(dest, ikpt_loc,1)
               jsppol = mpi_enreg%kpt_loc2ibz_sp(dest, ikpt_loc,2)

               if(jkpt > 0 .and. jsppol > 0) then

                 jkptf = dtefield%i2fbz(jkpt)
                 jkpt1f = dtefield%ikpt_dk(jkptf,ifor,idir)
                 jkpt1i = dtefield%indkk_f2ibz(jkpt1f,1)
                 his_source = mpi_enreg%proc_distrb(jkpt1i,1,jsppol)

                 if (his_source == me_distrb) then
!                  send

!                  pwnsfacq
                   ikg1 = dtefield%fkgindex(jkpt1f)
                   tag = jkpt1f + (jsppol - 1)*dtefield%fnkpt
                   count1 = npwarr(jkpt1i)
                   ABI_ALLOCATE(buffer,(2,count1))
                   buffer(:,1:count1)  = pwnsfac(:,ikg1+1:ikg1+count1)
                   call xsend_mpi(buffer,dest,tag,spaceComm_distrb,ierr)
                   ABI_DEALLOCATE(buffer)

!                  send cgq if necessary
                   if(flag_send(dest, jkpt1i)==0) then

                     icg1 = dtefield%cgindex(jkpt1i,jsppol)
                     tag = jkpt1i + (jsppol - 1)*dtset%nkpt
                     count1 = npwarr(jkpt1i)*nband_k*my_nspinor
                     ABI_ALLOCATE(buffer,(2,count1))
                     buffer(:,1:count1)  = cg(:,icg1+1:icg1+count1)
                     call xsend_mpi(buffer,dest,tag,spaceComm_distrb,ierr)
                     ABI_DEALLOCATE(buffer)

                     flag_send(dest, jkpt1i)=1
                   end if

                 end if ! end check that his_source == me

               end if ! end check on jkpt > 0 and jsppol > 0

             end if

           end do !dest

         end do !ifor

       end do !idir

       call timab(983,2,tsec)

       ABI_DEALLOCATE(flag_send)
       ABI_DEALLOCATE(flag_receive)
       if (ikpt_loc > dtset%mkmem) then
         ABI_DEALLOCATE(cgq)
         ABI_DEALLOCATE(pwnsfacq)
         cycle
       end if

     end if !berryopt and parallel over k pts

     if (abs(dtset%berryopt) == 5) then
       do idir = 1, 3
         do ifor = 1, 2
           dtefield%sflag(:,ikpt + dtset%nkpt*(isppol - 1),ifor,idir) = 0
         end do
       end do
     end if

     call timab(984,1,tsec)

!    Continue to initialize the Hamiltonian
     gs_hamk%istwf_k    =istwf_k
     gs_hamk%npw        =npw_k

     ABI_ALLOCATE(eig_k,(nband_k))
     ABI_ALLOCATE(ek_k,(nband_k))
     ABI_ALLOCATE(ek_k_nd,(nband_k,nband_k*paw_dmft%use_dmft))
     ABI_ALLOCATE(occ_k,(nband_k))
     ABI_ALLOCATE(resid_k,(nband_k))
     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     ABI_ALLOCATE(zshift,(nband_k))
     ABI_ALLOCATE(grnl_k,(3*natom,nband_k*optforces))
     if (psps%usepaw==0)  then
       ABI_ALLOCATE(enl_k,(nband_k))
     end if

     eig_k(:)=zero
     ek_k(:)=zero
     if(paw_dmft%use_dmft==1) ek_k_nd(:,:)=zero
     if (optforces>0) grnl_k(:,:)=zero
     if (psps%usepaw==0) enl_k(:)=zero
     kpoint(:)=dtset%kptns(:,ikpt)
     gs_hamk%kpoint(:)=dtset%kptns(:,ikpt)
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
     resid_k(:)=zero
     zshift(:)=dtset%eshift

     if (dtset%mkmem==0) then
!      Read (k+G) basis sphere data (same for each spin)
       nsp=dtset%nspinor
       call rdnpw(ikpt,isppol,nband_k,npw_k,nsp,0,dtfil%unkg)
!      Read k+g data
       read (dtfil%unkg) kg_k(1:3,1:npw_k)
       call sphereboundary(gs_hamk%gbound,istwf_k,kg_k,dtset%mgfft,npw_k)
!      Eventually read spherical harmonics
       if (psps%useylm==1) then
         read(dtfil%unylm)
         read(dtfil%unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,mpsang*mpsang)
       end if

!      Read the wavefunction block for ikpt,isppol
       call status(counter,dtfil%filstat,iexit,level,'read wfs      ')
       tim_rwwf=1
       ABI_ALLOCATE(eig_dum,(dtset%mband))
       ABI_ALLOCATE(kg_dum,(3,0))
       ABI_ALLOCATE(occ_dum,(dtset%mband))
       call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,dtset%mband,mcg_disk,mpi_enreg,nband_k,nband_k,&
&       npw_k,my_nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
       ABI_DEALLOCATE(eig_dum)
       ABI_DEALLOCATE(kg_dum)
       ABI_DEALLOCATE(occ_dum)

     else

       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       if (mpi_enreg%mode_para/='b'.or.istep<=1)&
&       call sphereboundary(gs_hamk%gbound,istwf_k,kg_k,dtset%mgfft,npw_k)

       if (psps%useylm==1) then
         do ilm=1,mpsang*mpsang
           ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
         end do
       end if

!      End if for choice governed by dtset%mkmem
     end if

!    Set up remaining of the Hamiltonian

!    Compute (1/2) (2 Pi)**2 (k+G)**2:
     call status(0,dtfil%filstat,iexit,level,'call mkkin    ')
     ABI_ALLOCATE(kinpw,(npw_k))
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg_k,kinpw,kpoint,npw_k)

!    Allocate the arrays phkxred and ph3d, compute phkxred
!    and eventually ph3d.
     do ia=1,natom
       iatom=atindx(ia)
       arg=two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
       gs_hamk%phkxred(1,iatom)=cos(arg)
       gs_hamk%phkxred(2,iatom)=sin(arg)
!      DEBUG
!      write(std_out,'(a,i4,2es16.6)' )&
!      &  'vtorho : iatom, phkxred',iatom,phkxred(1,iatom),phkxred(2,iatom)
!      ENDDEBUG
     end do
     if(dtset%nloalg(1)<=0)then
!      Here, only the allocation, not the precomputation.
       matblk=dtset%nloalg(4)
       ABI_ALLOCATE(ph3d,(2,npw_k,matblk))
     else
!      Here, allocation as well as precomputation
       matblk=natom
       ABI_ALLOCATE(ph3d,(2,npw_k,matblk))
       if (mpi_enreg%mode_para/='b'.or.istep<=1)&
&       call ph1d3d(1,natom,kg_k,matblk,natom,npw_k,n1,n2,n3,&
&       gs_hamk%phkxred,ph1d,ph3d)
     end if
     gs_hamk%matblk=matblk

!    Compute (k+G) vectors (only if useylm=1)
     nkpg=3*optforces*dtset%nloalg(5)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if ((mpi_enreg%mode_para/='b'.or.istep<=1).and.nkpg>0)&
&     call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

!    Compute nonlocal form factors ffnl at all (k+G):
     call status(0,dtfil%filstat,iexit,level,'call mkffnl   ')

     ider=0;idir=0;dimffnl=1
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
     if (mpi_enreg%mode_para/='b'.or.istep<=1)&
&     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&     gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&     npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&     psps%usepaw,psps%useylm,ylm_k,ylmgr)

!    Transpose the ffnl, kinpw, kpg and ph3d arrays.
     if (mpi_enreg%mode_para=='b'.and.istep<=1) then
       call prep_bandfft_tabs(dimffnl,ffnl,gs_hamk%gbound,ikpt,kinpw,&
&       kpoint,psps%lmnmax,matblk,dtset%mgfft,dtset%mkmem,mpi_enreg,nkpg,npw_k,ntypat,1,ph3d)
     end if

     call status(counter,dtfil%filstat,iexit,level,'call vtowfk   ')

     if(dtset%use_gpu_cuda==1) then
       if(mpi_enreg%mode_para=='b') then
         ikpt_this_proc = mpi_enreg%tab_kpt_distrib(ikpt)
         ndatarecv      = mpi_enreg%bandfft_kpt(ikpt_this_proc)%ndatarecv
#if defined HAVE_GPU_CUDA
         call gpu_update_ffnl_ph3d(mpi_enreg%bandfft_kpt(ikpt_this_proc)%ph3d_gather, &
&         mpi_enreg%bandfft_kpt(ikpt_this_proc)%ffnl_gather, &
&         ndatarecv,dimffnl,psps%lmnmax,ntypat,natom)
#endif
       else
#if defined HAVE_GPU_CUDA
         call gpu_update_ffnl_ph3d(ph3d,ffnl,npw_k,dimffnl,psps%lmnmax,ntypat,natom)
#endif
       end if
     end if

     call timab(984,2,tsec)

!    Compute the eigenvalues, wavefunction, residuals,
!    contributions to kinetic energy, nonlocal energy, forces,
!    and update of rhor to this k-point and this spin polarization.
     if(dtset%mkmem/=0)then
       if(with_vxctau)then
         call vtowfk(cg,cgq,cprj,cpus,dimcprj,dimffnl,dphase_k,dtefield,dtfil,&
&         dtset,eig_k,ek_k,ek_k_nd,enl_k,fixed_occ,ffnl,grnl_k,gs_hamk,&
&         ibg,icg,ikpt,iscf,isppol,kg_k,kinpw,kpg_k,&
&         psps%lmnmax,matblk,mband_cprj,mcg,mcgq,mcprj,dtset%mgfft,mkgq,&
&         mpi_enreg,psps%mpsang,psps%mpssoang,&
&         dtset%mpw,natom,nband_k,nkpg,dtset%nkpt,&
&         nnsclo_now,npw_k,npwarr,ntypat,nvloc,n4,n5,n6,&
&         occ_k,optforces,ph3d,prtvol,pwind,pwind_alloc,pwnsfac,pwnsfacq,resid_k,&
&         rhoaug,paw_dmft%use_dmft,usecprj,vlocal,dtset%wtk(ikpt),zshift,vxctaulocal=vxctaulocal)
       else
         call vtowfk(cg,cgq,cprj,cpus,dimcprj,dimffnl,dphase_k,dtefield,dtfil,&
&         dtset,eig_k,ek_k,ek_k_nd,enl_k,fixed_occ,ffnl,grnl_k,gs_hamk,&
&         ibg,icg,ikpt,iscf,isppol,kg_k,kinpw,kpg_k,&
&         psps%lmnmax,matblk,mband_cprj,mcg,mcgq,mcprj,dtset%mgfft,mkgq,&
&         mpi_enreg,psps%mpsang,psps%mpssoang,&
&         dtset%mpw,natom,nband_k,nkpg,dtset%nkpt,&
&         nnsclo_now,npw_k,npwarr,ntypat,nvloc,n4,n5,n6,&
&         occ_k,optforces,ph3d,prtvol,pwind,pwind_alloc,pwnsfac,pwnsfacq,resid_k,&
&         rhoaug,paw_dmft%use_dmft,usecprj,vlocal,dtset%wtk(ikpt),zshift)
       end if
     else if(dtset%mkmem==0)then
       if(with_vxctau)then
         call vtowfk(cg_disk,cgq,cprj,cpus,dimcprj,dimffnl,dphase_k,dtefield,dtfil,&
&         dtset,eig_k,ek_k,ek_k_nd,enl_k,fixed_occ,ffnl,grnl_k,gs_hamk,&
&         ibg,icg,ikpt,iscf,isppol,kg_k,kinpw,kpg_k,&
&         psps%lmnmax,matblk,mband_cprj,mcg_disk,mcprj,mcgq,dtset%mgfft,mkgq,&
&         mpi_enreg,psps%mpsang,psps%mpssoang,&
&         dtset%mpw,natom,nband_k,nkpg,dtset%nkpt,nnsclo_now,npw_k,npwarr,ntypat,nvloc,n4,n5,n6,&
&         occ_k,optforces,ph3d,prtvol,pwind,pwind_alloc,pwnsfac,pwnsfacq,resid_k,&
&         rhoaug,paw_dmft%use_dmft,usecprj,vlocal,dtset%wtk(ikpt),zshift,vxctaulocal=vxctaulocal)
       else
         call vtowfk(cg_disk,cgq,cprj,cpus,dimcprj,dimffnl,dphase_k,dtefield,dtfil,&
&         dtset,eig_k,ek_k,ek_k_nd,enl_k,fixed_occ,ffnl,grnl_k,gs_hamk,&
&         ibg,icg,ikpt,iscf,isppol,kg_k,kinpw,kpg_k,&
&         psps%lmnmax,matblk,mband_cprj,mcg_disk,mcgq,mcprj,dtset%mgfft,mkgq,&
&         mpi_enreg,psps%mpsang,psps%mpssoang,&
&         dtset%mpw,natom,nband_k,nkpg,dtset%nkpt,nnsclo_now,npw_k,npwarr,ntypat,nvloc,n4,n5,n6,&
&         occ_k,optforces,ph3d,prtvol,pwind,pwind_alloc,pwnsfac,pwnsfacq,resid_k,&
&         rhoaug,paw_dmft%use_dmft,usecprj,vlocal,dtset%wtk(ikpt),zshift)
       end if
     end if
     call status(counter,dtfil%filstat,iexit,level,'after vtowfk  ')

     call timab(985,1,tsec)

!    update magnetizations at this k point
!    do this here BEFORE ffnl, kpg_k, ph3d are deallocated
     if (abs(dtset%berryopt) == 5) then

!      get the cprj on this processor for this kpt, store in input order
       call cprj_get(atindx1,cprj_berry,cprj,natom,1,ibg,ikpt,1,1,dtset%mband,&
&       dtset%mkmem,mpi_enreg,natom,dtefield%nband_occ,dtefield%nband_occ,&
&       my_nspinor,dtset%nsppol,0)

!      gather on each processor the list of kpts all processors are currently handling
       call xallgather_mpi(ikpt,ikpt_recv,spaceComm_distrb,ierr)

!      gather the current cprjs to all processes
       call cprj_mpi_allgather(cprj_berry,cprj_gat,natom,dtefield%nband_occ,&
&       dimcprj_u,0,nproc_distrb,spaceComm_distrb,ierr)
       do iproc = 1, nproc_distrb
         icp2=dtefield%nband_occ*(iproc-1)
         call cprj_get(atindx1,cprj_berry,cprj_gat,natom,1,icp2,ikpt,0,1,&
&         dtefield%nband_occ,nproc_distrb,mpi_enreg,natom,&
&         dtefield%nband_occ,dtefield%nband_occ,1,1,0)
         icp1 = dtefield%nband_occ*(ikpt_recv(iproc)-1)
         call cprj_put(atindx1,cprj_berry,dtefield%cprj,natom,1,icp1,ikpt,0,1,&
&         dtefield%nband_occ,dtset%nkpt,mpi_enreg,natom,dtefield%nband_occ,&
&         dtefield%nband_occ,dimcprj_u,1,1,spaceComm_distrb,0)
       end do

       call mag_loc_k(atindx,atindx1,cg,cprj_berry,dimffnl,dtefield,ffnl,gmet,gprimd,icg,&
&       ikpt,psps%indlmn,istwf_k,kg_k,kpg_k,kpoint,psps%lmnmax,matblk,mcg,&
&       dtset%mgfft,mpi_enreg,psps%mpsang,psps%mpssoang,natom,nattyp,nband_k,&
&       dtset%ngfft,nkpg,dtset%nloalg,npw_k,&
&       my_nspinor,ntypat,pawtab,gs_hamk%phkxred,ph1d,ph3d,ucvol)

!      
!      dimenl1 = psps%lmnmax*(psps%lmnmax+1)/2
!      call mag_nonloc_k(atindx1,cg,dimenl1,dimffnl,dtefield,gs_hamk%ekb,ffnl,gmet,gprimd,icg,&
!      &             ikpt,psps%indlmn,istwf_k,kg_k,kpg_k,kpoint,psps%lmnmax,matblk,mcg,&
!      &             dtset%mgfft,mpi_enreg,psps%mpsang,psps%mpssoang,natom,nattyp,nband_k,&
!      &             dtset%ngfft,nkpg,dtset%nloalg,npw_k,&
!      &               my_nspinor,ntypat,gs_hamk%phkxred,ph1d,ph3d,ucvol)
!      
       call update_mmat(dtset%berryopt,cg,cgq,dimffnl,dtefield,ffnl,dtfil%filstat,&
&       gmet,gprimd,gs_hamk,icg,ikpt,kg,kinpw,psps%lmnmax,matblk,dtset%mband,&
&       mcg,mcgq,dtset%mgfft,mkgq,dtset%mkmem,mpi_enreg,&
&       psps%mpsang,psps%mpssoang,dtset%mpw,natom,nkpg,dtset%nkpt,npw_k,npwarr,&
&       dtset%nspinor,ntypat,nvloc,n4,n5,n6,pawtab,psps,&
&       pwind,pwind_alloc,dtset%paral_kgb,ph3d,dtset%prtvol,&
&       pwnsfac,pwnsfacq,rmet,ucvol,vlocal,ylm,ylmgr)
!      
       call magcart(dtefield,mpi_enreg,dtset%nkpt,rprimd,dtset%wtk)
!      
     end if

#if defined HAVE_GPU_CUDA
     if(dtset%use_gpu_cuda==1) then
       call gpu_finalize_ffnl_ph3d()
     end if
#endif
     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(kinpw)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(cgq)
     ABI_DEALLOCATE(pwnsfacq)

!    electric field
     if (dtset%berryopt == 4) then

       dphasek(:,ikpt + (isppol - 1)*dtset%nkpt) = dphase_k(:)

!      The overlap matrices for all first neighbours of ikpt
!      are no more up to date
       do idir = 1, 3
         do ifor = 1, 2
           ikpt1 = dtefield%ikpt_dk(dtefield%i2fbz(ikpt),ifor,idir)
           ikpt1 = dtefield%indkk_f2ibz(ikpt1,1)
           ifor1 = -1*ifor + 3   ! ifor = 1 -> ifor1 = 2 & ifor = 2 -> ifor1 = 1
           dtefield%sflag(:,ikpt1+(isppol-1)*dtset%nkpt,ifor1,idir) = 0
         end do
       end do

     end if  ! berryopt

!    Save eigenvalues (hartree), residuals (hartree**2)
     eigen(1+bdtot_index : nband_k+bdtot_index) = eig_k(:)
     eknk (1+bdtot_index : nband_k+bdtot_index) = ek_k (:)
     if(paw_dmft%use_dmft==1) eknk_nd(isppol,ikpt,:,:) = ek_k_nd(:,:)
     resid(1+bdtot_index : nband_k+bdtot_index) = resid_k(:)
     if (optforces>0) grnlnk(:,1+bdtot_index : nband_k+bdtot_index) = grnl_k(:,:)
     if (psps%usepaw==0) enlnk(1+bdtot_index : nband_k+bdtot_index) = enl_k(:)

     if(iscf>0 .or. iscf==-3)then
!      Accumulate sum over k points for band, nonlocal and kinetic energies,
!      also accumulate gradients of Enonlocal:
       do iband=1,nband_k
         if (abs(occ_k(iband))>tol8) then
           energies%e_kinetic = energies%e_kinetic + &
&           dtset%wtk(ikpt)*occ_k(iband)*ek_k(iband)
           energies%e_eigenvalues = energies%e_eigenvalues + &
&           dtset%wtk(ikpt)*occ_k(iband)*eig_k(iband)
           if (optforces>0) grnl(:)=grnl(:)+dtset%wtk(ikpt)/ucvol*occ_k(iband)*grnl_k(:,iband)
           if (psps%usepaw==0) then
             energies%e_nonlocalpsp = energies%e_nonlocalpsp + &
&             dtset%wtk(ikpt)*occ_k(iband)*enl_k(iband)
           end if
         end if
       end do
     end if

     call timab(985,2,tsec)

!    Write new wavefunctions to disk if needed
!    (header records were written earlier in loopcv)
     if (dtset%mkmem==0) then
       tim_rwwf=1
       call rwwf(cg_disk,eig_k,0,0,0,ikpt,isppol,kg_k,dtset%mband,mcg_disk,mpi_enreg,nband_k,nband_k,&
&       npw_k,my_nspinor,occ_k,2,1,tim_rwwf,wffnew)
     end if

     ABI_DEALLOCATE(eig_k)
     ABI_DEALLOCATE(ek_k)
     ABI_DEALLOCATE(ek_k_nd)
     ABI_DEALLOCATE(grnl_k)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(resid_k)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(zshift)
     if (psps%usepaw==0)  then
       ABI_DEALLOCATE(enl_k)
     end if

!    Keep track of total number of bands (all k points so far, even for
!    k points not treated by me)
     bdtot_index=bdtot_index+nband_k

!    Also shift array memory if dtset%mkmem/=0
     if (dtset%mkmem/=0) then
       ibg=ibg+my_nspinor*nband_k
       icg=icg+npw_k*my_nspinor*nband_k
       ikg=ikg+npw_k
     end if


!    End big k point loop
   end do

   call status(counter,dtfil%filstat,iexit,level,'after k loop  ')

   call timab(986,1,tsec)


   if (fixed_occ .and. mpi_enreg%mode_para=='b') then
     call xsum_mpi(rhoaug,mpi_enreg%commcart_3d,ierr) !Sum the contributions over bands/FFT/spinors
   end if

!  Transfer density on augmented fft grid to normal fft grid in real space
!  Also take into account the spin.
   if(iscf>0.or.iscf==-3)then
     if( mpi_enreg%paral_compil_kpt==0 .or.       &
&     mpi_enreg%paralbd <= 1        .or.       &
&     (mpi_enreg%paralbd >= 1 .and. mpi_enreg%me_group==0)) then
       if (psps%usepaw==0) then
         call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug(:,:,:,1),1)
         if(dtset%nspden==4)then
           do imagn=2,4
             call fftpac(imagn,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug(:,:,:,imagn),1)
           end do
         end if
       else
         call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhowfr,rhoaug(:,:,:,1),1)
         if(dtset%nspden==4)then
           do imagn=2,4
             call fftpac(imagn,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhowfr,rhoaug(:,:,:,imagn),1)
           end do
         end if
       end if
     end if
   end if

   call timab(986,2,tsec)

!  End loop over spins
 end do

 call status(counter,dtfil%filstat,iexit,level,'after spinloop')

 if(mpi_enreg%paral_compil_kpt==1)then
   call timab(987,1,tsec)
   call leave_test()
   write(message,*) 'vtorho: loop on k-points and spins done in parallel'
   call wrtout(std_out,message,'COLL')
   call timab(987,2,tsec)
 end if


 call timab(988,1,tsec)
!electric field: compute string-averaged change in Zak phase
!along each direction, store it in dphase(idir)

!ji: it is not convenient to do this anymore. Remove. Set dphase(idir)=0.0_dp.
!eventually, dphase(idir) will have to go...

 if (dtset%berryopt == 4)  dphase(:) = 0.0_dp


!In case of MPI // of a finite field calculation, send dphasek to all cpus
 if ((mpi_enreg%paral_compil_kpt == 1).and.(dtset%berryopt == 4)) then
   call xsum_mpi(dphasek,spaceComm_distrb,ierr)
 end if

 if (dtset%berryopt == 4)  then
   ABI_DEALLOCATE(dphasek)
 end if

 call destroy_hamiltonian(gs_hamk)

#if defined HAVE_GPU_CUDA
 if(dtset%use_gpu_cuda==1) then
   call gpu_finalize_ham_data()
 end if
#endif

 ABI_DEALLOCATE(rhoaug)
 ABI_DEALLOCATE(vlocal)
 if(dtset%mkmem==0)ABI_DEALLOCATE(cg_disk)

 call status(0,dtfil%filstat,iexit,level,'after loops   ')

 ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
 doccde(:)=zero !MF initialize

 call timab(988,2,tsec)

!Treat now varying occupation numbers, in the self-consistent case
 if((.not.fixed_occ) .and. (iscf>0.or.iscf==-3)) then

!  Parallel case
   if( mpi_enreg%paral_compil_kpt==1 .or. mpi_enreg%paral_compil_fft==1)then

     call timab(989,1,tsec)

!    If needed, exchange the values of eigen,resid,eknk,enlnk,grnlnk
     ABI_ALLOCATE(buffer1,((4+3*natom*optforces-psps%usepaw)*mbdkpsp))
     if(paw_dmft%use_dmft==1) then
       ABI_ALLOCATE(buffer2,(mb2dkpsp*paw_dmft%use_dmft))
     end if
!    Pack eigen,resid,eknk,enlnk,grnlnk in buffer1
     buffer1(1          :  mbdkpsp)=eigen(:)
     buffer1(1+  mbdkpsp:2*mbdkpsp)=resid(:)
     buffer1(1+2*mbdkpsp:3*mbdkpsp)=eknk(:)
     if(paw_dmft%use_dmft==1) then
       nnn=0
       do ikpt=1,dtset%nkpt
         do isppol=1,dtset%nsppol
           do iband=1,dtset%mband
             do iband1=1,dtset%mband
               nnn=nnn+1
               buffer2(nnn)=eknk_nd(isppol,ikpt,iband,iband1)
             end do
           end do
         end do
       end do
       if(nnn.ne.mb2dkpsp) stop
     end if
     index1=3*mbdkpsp
     if (psps%usepaw==0) then
       buffer1(index1+1:index1+mbdkpsp)=enlnk(:)
       index1=index1+mbdkpsp
     end if
     if (optforces>0) then
       buffer1(index1+1:index1+3*natom*mbdkpsp)=reshape(grnlnk, (/(3*natom)*mbdkpsp/) )
     end if
!    Build sum of everything
     call timab(48,1,tsec)

!    PATCH vtorho // KPT & FFT sum fft_master_comm --> comm_kpt
     if(mpi_enreg%mode_para/='b') then
       call xsum_mpi(buffer1,mpi_enreg%fft_master_comm,ierr)
       if(paw_dmft%use_dmft==1) call xsum_mpi(buffer2,mpi_enreg%fft_master_comm,ierr)
     else
       if ((mpi_enreg%paral_compil_kpt==1) .and. &
&       (mpi_enreg%paral_compil_fft==1)) then
         call xsum_mpi(buffer1,mpi_enreg%comm_kpt ,ierr)
       end if
     end if

     call timab(48,2,tsec)
!    Unpack eigen,resid,eknk,enlnk,grnlnk in buffer1
     eigen(:) =buffer1(1          :  mbdkpsp)
     resid(:) =buffer1(1+  mbdkpsp:2*mbdkpsp)
     eknk(:)  =buffer1(1+2*mbdkpsp:3*mbdkpsp)
     if(paw_dmft%use_dmft==1) then
       nnn=0
       do ikpt=1,dtset%nkpt
         do isppol=1,dtset%nsppol
           do iband=1,dtset%mband
             do iband1=1,dtset%mband
               nnn=nnn+1
               eknk_nd(isppol,ikpt,iband,iband1)=buffer2(nnn)
             end do
           end do
         end do
       end do
     end if
     index1=3*mbdkpsp
     if (psps%usepaw==0) then
       enlnk(:) =buffer1(index1+1:index1+mbdkpsp)
       index1=index1+mbdkpsp
     end if
     if (optforces>0) then
       grnlnk(:,:)=reshape(buffer1(index1+1:index1+3*natom*mbdkpsp),&
&       (/ 3*natom , mbdkpsp /) )
     end if
     if(allocated(buffer2))  then
       ABI_DEALLOCATE(buffer2)
     end if
     ABI_DEALLOCATE(buffer1)
     call timab(989,2,tsec)

   end if ! parallel

   call timab(990,1,tsec)

!  Compute the new occupation numbers from eigen
   call status(0,dtfil%filstat,iexit,level,'call newocc   ')
   call newocc(doccde,eigen,energies%entropy,energies%e_fermie,dtset%fixmom,&
&   dtset%mband,dtset%nband,dtset%nelect,dtset%nkpt,dtset%nspinor,&
&   dtset%nsppol,occ,dtset%occopt,prtvol,dtset%stmbias,dtset%tphysel,&
&   dtset%tsmear,dtset%wtk)

   call timab(990,2,tsec)

!  !=========  DMFT call begin ============================================
   dmft_ldaocc=0
   if(paw_dmft%use_dmft==1.and.psps%usepaw==1) then

     call timab(991,1,tsec)

!    ==  0 to a dmft calculation and do not use lda occupations
!    ==  1 to a lda calculation with the dmft loop
     if(dtset%dmftcheck==-1) dmft_ldaocc=1
!    write(std_out,*) "dmft_ldaocc",dmft_ldaocc

!    ==  initialise occnd
     paw_dmft%occnd=zero

     if(dmft_ldaocc==0) then
       if(dtset%occopt/=3) then
         write(message,'(a,a,a,a,a,a)')  ch10,&
&         ' occopt should be equal to 3 in dmft',ch10
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if

!      ==  initialise edmft
       if(paw_dmft%use_dmft>=1) edmft = zero

!      !=================================================================
!      ==  allocate paw_dmft%psichi and paw_dmft%eigen_lda
       call init_dmft(dtset,energies%e_fermie,dtfil%fnameabo_app,dtset%nspinor,paw_dmft,pawtab,psps,dtset%typat)
       call print_dmft(paw_dmft,dtset%pawprtvol)

!      ==  gather crystal structure date into data "cryst_struc"
       remove_inv=.false.
       if(dtset%nspden==4) remove_inv=.true.
       call init_crystal(cryst_struc,dtset%spgroup,natom,dtset%npsp,ntypat, &
&       dtset%nsym,rprimd,dtset%typat,xred,dtset%ziontypat,dtset%znucl,1,&
&       dtset%nspden==2.and.dtset%nsppol==1,remove_inv,hdr%title,&
&       dtset%symrel,dtset%tnons,dtset%symafm)

!      ==  compute psichi
       call xbarrier_mpi(spaceComm_distrb)
       call init_oper(paw_dmft,lda_occup)
       call datafordmft(cryst_struc,cprj,dimcprj,dtset,eigen,energies%e_fermie,&
&       lda_occup,dtset%mband,dtset%mkmem,mpi_enreg,&
&       dtset%nkpt,my_nspinor,dtset%nsppol,occ,&
&       paw_dmft,paw_ij,pawang,pawtab,psps,dtfil%unpaw)

!      ==  solve dmft loop
       call xbarrier_mpi(spaceComm_distrb)

       call dmft_solve(cryst_struc,istep,lda_occup,mpi_enreg,paw_dmft,pawang,pawtab,dtset%pawprtvol)
       edmft=paw_dmft%edmft
       energies%e_paw=energies%e_paw+edmft
       energies%e_pawdc=energies%e_pawdc+edmft
!      paw_dmft%occnd(:,:,:,:)=0.5_dp

!      call print_dmft(paw_dmft,dtset%pawprtvol)

       call xbarrier_mpi(spaceComm_distrb)
!      call leave_new('COLL')
       call destroy_dmft(paw_dmft)
       call xbarrier_mpi(spaceComm_distrb)

!      ==  destroy crystal_structure cryst_struc
       call destroy_crystal(cryst_struc)
       call destroy_oper(lda_occup)
     end if ! dmft_ldaocc

     call timab(991,2,tsec)

   end if ! usedmft
!  !=========  DMFT call end   ============================================

   call timab(992,1,tsec)

!  Compute eeig, ek,enl and grnl from the new occ, and the shared eknk,enlnk,grnlnk
   energies%e_eigenvalues = zero
   energies%e_kinetic     = zero
   energies%e_nonlocalpsp = zero
   if(paw_dmft%use_dmft>=1) then
     ebandlda               = zero
     ebanddmft              = zero
     ebandldatot            = zero
     ekindmft               = zero
     ekindmft2              = zero
     ekinlda                = zero
   end if
!  compute new energy terms due to non diagonal occupations and DMFT.
!  
   if (psps%usepaw==0) energies%e_nonlocalpsp = zero
   if (optforces>0) grnl(:)=zero
   bdtot_index=1
   do isppol=1,dtset%nsppol
     do ikpt=1,dtset%nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       do iband=1,nband_k

         locc_test = abs(occ(bdtot_index))>tol8
!        dmft
         if(paw_dmft%use_dmft>=1) then
           if(paw_dmft%band_in(iband)) then
             if( paw_dmft%use_dmft == 1 .and. dmft_ldaocc == 1 ) then ! test of the code
               paw_dmft%occnd(iband,iband,ikpt,isppol)=occ(bdtot_index)
             end if
             locc_test = abs(paw_dmft%occnd(iband,iband,ikpt,isppol))>tol8
           end if
         end if

         if (locc_test) then
!          dmft
           if(paw_dmft%use_dmft==1) then
             ebandldatot=ebandldatot+dtset%wtk(ikpt)*occ(bdtot_index)*eigen(bdtot_index)
             if(paw_dmft%band_in(iband)) then
               ebandlda=ebandlda+dtset%wtk(ikpt)*occ(bdtot_index)*eigen(bdtot_index)
               ekinlda=ekinlda+dtset%wtk(ikpt)*occ(bdtot_index)*eknk(bdtot_index)
               occ(bdtot_index)=paw_dmft%occnd(iband,iband,ikpt,isppol)
               ebanddmft=ebanddmft+dtset%wtk(ikpt)*occ(bdtot_index)*eigen(bdtot_index)
               ekindmft=ekindmft+dtset%wtk(ikpt)*occ(bdtot_index)*eknk(bdtot_index)
             end if
           end if

           energies%e_eigenvalues = energies%e_eigenvalues + &
&           dtset%wtk(ikpt)*occ(bdtot_index)*eigen(bdtot_index)
           energies%e_kinetic = energies%e_kinetic + &
&           dtset%wtk(ikpt)*occ(bdtot_index)*eknk(bdtot_index)
           if (optforces>0) grnl(:)=grnl(:)+dtset%wtk(ikpt)/ucvol*occ(bdtot_index)*grnlnk(:,bdtot_index)
           if (psps%usepaw==0) energies%e_nonlocalpsp = energies%e_nonlocalpsp + &
&           dtset%wtk(ikpt)*occ(bdtot_index)*enlnk(bdtot_index)
         end if
         bdtot_index=bdtot_index+1
         if(paw_dmft%use_dmft==1) then
           do iband1=1,nband_k
             if(paw_dmft%band_in(iband).and.paw_dmft%band_in(iband1)) then
!              write(std_out,*) "II+", isppol,ikpt,iband,iband1
               ekindmft2=ekindmft2+dtset%wtk(ikpt)*paw_dmft%occnd(iband,iband1,ikpt,isppol)*&
&               eknk_nd(isppol,ikpt,iband,iband1)
!              write(std_out,*) "II", occnd(iband,iband1,ikpt,isppol),eknk_nd(isppol,ikpt,iband,iband1)
             end if
           end do
         end if
       end do
     end do
   end do

   if(paw_dmft%use_dmft==1) then
     energies%e_kinetic = energies%e_kinetic -ekindmft+ekindmft2
     if(abs(dtset%pawprtvol)>=2) then
       write(message,'(4a,7(a,2x,e12.5,a),a)') &
&       "-----------------------------------------------",ch10,&
&       "--- Energy for DMFT and tests (in Ha)  ",ch10,&
&       "--- Ebandldatot    (Ha.) = ",ebandldatot,ch10,&
&       "--- Ebandlda       (Ha.) = ",ebandlda,ch10,&
&       "--- Ebanddmft      (Ha.) = ",ebanddmft,ch10,&
&       "--- Ekinlda        (Ha.) = ",ekinlda,ch10, &
&       "--- Ekindmftdiag   (Ha.) = ",ekindmft,ch10,&
&       "--- Ekindmftnondiag(Ha.) = ",ekindmft2,ch10,&
&       "--- Edmft=         (Ha.) = ",edmft,ch10,&
&       "-----------------------------------------------"
       call wrtout(std_out,message,'COLL')
     end if
   end if


   call status(0,dtfil%filstat,iexit,level,'call mkrho    ')
   tim_mkrho=2
   if (psps%usepaw==0) then
     call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,&
&     rhog  ,rhor  ,rprimd,tim_mkrho,ucvol,dtfil%unkg,wffnew,wvl%wfs,wvl%descr)
   else
     call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,&
&     rhowfg,rhowfr,rprimd,tim_mkrho,ucvol,dtfil%unkg,wffnew,wvl%wfs,wvl%descr)
   end if
   call timab(992,2,tsec)

!  Treat fixed occupation numbers or non-self-consistent case
 else

   if(mpi_enreg%paral_compil_kpt==1 .or. mpi_enreg%paral_compil_fft==1)then

     call timab(989,1,tsec)

     nbuf=2*mbdkpsp+dtset%nfft*dtset%nspden+3-psps%usepaw+3*natom*optforces
     if(iscf==-1 .or. iscf==-2)nbuf=2*mbdkpsp
     ABI_ALLOCATE(buffer1,(nbuf))
!    Pack eigen,resid,rho[wf]r,grnl,enl,ek
     buffer1(1:mbdkpsp)=eigen(:)
     buffer1(1+mbdkpsp:2*mbdkpsp)=resid(:)
     index1=2*mbdkpsp
     if(iscf/=-1 .and. iscf/=-2)then
       if (psps%usepaw==0) then
         buffer1(index1+1:index1+dtset%nfft*dtset%nspden)=reshape(rhor  ,&
&         (/dtset%nfft*dtset%nspden/))
       else
         buffer1(index1+1:index1+dtset%nfft*dtset%nspden)=reshape(rhowfr,&
&         (/dtset%nfft*dtset%nspden/))
       end if
       index1=index1+dtset%nfft*dtset%nspden
       buffer1(index1+1) = energies%e_kinetic
       buffer1(index1+2) = energies%e_eigenvalues
       if (psps%usepaw==0) buffer1(index1+3) = energies%e_nonlocalpsp
       index1=index1+3-psps%usepaw
       if (optforces>0) buffer1(index1+1:index1+3*natom)=grnl(1:3*natom)
     end if
!    Build sum of everything
     call timab(48,1,tsec)

!    PATCH vtorho // KPT & FFT sum fft_master_comm --> comm_kpt
     if(mpi_enreg%mode_para/='b') then
       call xsum_mpi(buffer1,nbuf,mpi_enreg%fft_master_comm,ierr)
     else
       if ((mpi_enreg%paral_compil_kpt==1) .and. &
&       (mpi_enreg%paral_compil_fft==1)) then
         call xsum_mpi(buffer1,nbuf,mpi_enreg%comm_kpt ,ierr)
       end if
     end if
     call timab(48,2,tsec)
!    Unpack the final result
     eigen(:)=buffer1(1:mbdkpsp)
     resid(:)=buffer1(1+mbdkpsp:2*mbdkpsp)
     index1=2*mbdkpsp
     if(iscf/=-1 .and. iscf/=-2)then
       if (psps%usepaw==0) then
         ii=1
         do ispden=1,dtset%nspden
           do ifft=1,dtset%nfft
             rhor(ifft,ispden)=buffer1(index1+ii)
             ii=ii+1
           end do
         end do
       else
         ii=1
         do ispden=1,dtset%nspden
           do ifft=1,dtset%nfft
             rhowfr(ifft,ispden)=buffer1(index1+ii)
             ii=ii+1
           end do
         end do
       end if
       index1=index1+dtset%nfft*dtset%nspden
       energies%e_kinetic = buffer1(index1+1)
       energies%e_eigenvalues = buffer1(index1+2)
       if (psps%usepaw==0) energies%e_nonlocalpsp = buffer1(index1+3)
       index1=index1+3-psps%usepaw
       if (optforces>0) grnl(1:3*natom)=buffer1(index1+1:index1+3*natom)
     end if
     ABI_DEALLOCATE(buffer1)
     call timab(989,2,tsec)

   end if ! parallel

   call timab(993,1,tsec)

!  Compute the highest occupied eigenenergy 
   if(iscf/=-1 .and. iscf/=-2)then
     energies%e_fermie = -huge(one)
     bdtot_index=1
     do isppol=1,dtset%nsppol
       do ikpt=1,dtset%nkpt
         nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
         do iband=1,nband_k
           if(abs(occ(bdtot_index))>tol8 .and. eigen(bdtot_index)>energies%e_fermie+tol10) then
             energies%e_fermie=eigen(bdtot_index)
           end if
           bdtot_index=bdtot_index+1
         end do
       end do
     end do
     if(mpi_enreg%mode_para/='b') then
       call xmax_mpi(energies%e_fermie,emax,spaceComm_distrb,ierr)
       energies%e_fermie=emax
     else
       if ((mpi_enreg%paral_compil_kpt==1) .and. &
&       (mpi_enreg%paral_compil_fft==1)) then
         call xmax_mpi(energies%e_fermie,emax,mpi_enreg%comm_kpt,ierr)
         energies%e_fermie=emax
       end if
     end if

   end if

   call timab(993,2,tsec)

!  If needed, compute rhog, and symmetrizes the density
   if (iscf > 0 .or. iscf==-3 ) then

!    energies%e_fermie=zero  ! Actually, should determine the maximum of the valence band XG20020802

     call timab(994,1,tsec)

     call status(0,dtfil%filstat,iexit,level,'compute rhog  ')
     nfftot=dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)
     if (psps%usepaw==0) then
       call symrhg(1,gprimd,irrzon,mpi_enreg,dtset%nfft,nfftot,dtset%ngfft,dtset%nspden,&
&       dtset%nsppol,dtset%nsym,dtset%paral_kgb,phnons,rhog  ,rhor  ,rprimd,dtset%symafm,dtset%symrel)
     else
       call symrhg(1,gprimd,irrzon,mpi_enreg,dtset%nfft,nfftot,dtset%ngfft,dtset%nspden,&
&       dtset%nsppol,dtset%nsym,dtset%paral_kgb,phnons,rhowfg,rhowfr,rprimd,dtset%symafm,dtset%symrel)
     end if
!    We now have both rho(r) and rho(G), symmetrized, and if dtset%nsppol=2
!    we also have the spin-up density, symmetrized, in rhor(:,2).

     call timab(994,2,tsec)

   end if

!  End of test on varying or fixed occupation numbers
 end if

 call timab(994,1,tsec)

!Compute the kinetic energy density
 if(dtset%usekden==1 .and. (iscf > 0 .or. iscf==-3 ) )then
   tim_mkrho=2
   call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,npwarr,occ,paw_dmft,phnons,&
&   taug  ,taur  ,rprimd,tim_mkrho,ucvol,dtfil%unkg,wffnew,wvl%wfs,wvl%descr,option=1)
 end if

 ABI_DEALLOCATE(eknk)
 ABI_DEALLOCATE(eknk_nd)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(grnlnk)
 if (psps%usepaw==0)  then
   ABI_DEALLOCATE(enlnk)
 end if

!In the self-consistent case, diagnose lack of unoccupied
!state (for each spin and k-point).
!Print a warning if the number of such messages already written
!does not exceed mwarning.
 mwarning=5
 if(nwarning<mwarning .and. iscf>0)then
   nwarning=nwarning+1
   bdtot_index=1
   do isppol=1,dtset%nsppol
     do ikpt=1,dtset%nkpt
       min_occ=two
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       do iband=1,nband_k
         if(occ(bdtot_index)<min_occ)min_occ=occ(bdtot_index)
         bdtot_index=bdtot_index+1
       end do
       if(min_occ>0.01_dp)then
         if(dtset%nsppol==1)then
           write(message, '(a,a,a,a,i4,a,a,a,f7.3,a,a,a,a,a,a,a)' )ch10,&
&           ' vtorho : WARNING -',ch10,&
&           '  For k-point number ',ikpt,',',ch10,&
&           '  The minimal occupation factor is',min_occ,'.',ch10,&
&           '  An adequate monitoring of convergence requires it to be',&
&           ' at most 0.01_dp.',ch10,&
&           '  Action : increase slightly the number of bands.',ch10
         else
           write(message, '(a,a,a,a,i4,a,a,a,i3,a,f7.3,a,a,a,a,a,a,a)' )ch10,&
&           ' vtorho : WARNING -',ch10,&
&           '  For k-point number ',ikpt,', and',ch10,&
&           '  for spin polarization',isppol,&
&           ', the minimal occupation factor is',min_occ,'.',ch10,&
&           '  An adequate monitoring of convergence requires it to be',&
&           ' at most 0.01_dp.',ch10,&
&           '  Action : increase slightly the number of bands.',ch10
         end if
         call wrtout(std_out,message,'COLL')
!        It is enough if one lack of adequate occupation is identified, so exit.
         exit
       end if
     end do
   end do
 end if
!In the non-self-consistent case, print eigenvalues and residuals
 if(iscf<=0)then
   option=2 ; enunit=1 ; vxcavg_dum=zero
   call prteigrs(eigen,enunit,energies%e_fermie,dtfil%fnameabo_app_eig,&
&   ab_out,iscf,dtset%kptns,dtset%kptopt,dtset%mband,dtset%nband,&
&   dtset%nkpt,nnsclo_now,dtset%nsppol,occ,dtset%occopt,option,&
&   dtset%prteig,prtvol,resid,dtset%tolwfr,vxcavg_dum,dtset%wtk)
 end if

!Find largest residual over bands, k points, and spins,
!except for nbdbuf highest bands
 ibdkpt=1
 residm=zero
 do isppol=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     nband_eff=max(1,nband_k-dtset%nbdbuf)
     residm=max(residm,maxval(resid(ibdkpt:ibdkpt+nband_eff-1)))
     ibdkpt=ibdkpt+nband_k
   end do
 end do

 if (iscf>0.or.iscf==-3) then

!  PAW: Build new rhoij quantities from new occ then symetrize them
!  Compute and add the compensation density to rhowfr to get the total density
   if (psps%usepaw==1) then
     call timab(555,1,tsec)
!    Build unpacked rhoij
     call status(istep,dtfil%filstat,iexit,level,'call pawmkrhoij')
     call rhoij_init_unpacked(pawrhoij)
     if (usecprj==1) then
       call pawmkrhoij(atindx1,cprj,dimcprj,dtset%istwfk,dtset%kptopt,dtset%mband,mband_cprj,&
&       mcprj,dtset%mkmem,mpi_enreg,natom,dtset%nband,dtset%nkpt,dtset%nspinor,dtset%nsppol,&
&       occ,dtset%paral_kgb,paw_dmft,dtset%pawprtvol,pawrhoij,dtfil%unpaw,dtset%wtk)
     else
       mcprj_tmp=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol*usecprj
       ABI_ALLOCATE(cprj_tmp,(natom,mcprj_tmp))
       call cprj_alloc(cprj_tmp,0,dimcprj)
       call ctocprj(atindx,cg,1,cprj_tmp,gmet,gprimd,0,0,0,dtset%istwfk,kg,dtset%kptns,&
&       dtset%mband,mcg,mcprj_tmp,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,&
&       dtset%natom,nattyp,dtset%nband,dtset%natom,dtset%ngfft,dtset%nkpt,dtset%nloalg,&
&       npwarr,dtset%nspinor,dtset%nsppol,ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,&
&       ucvol,dtfil%unpaw,dtfil%unkg,dtfil%unylm,0,wffnew,xred,ylm,ylmgr_dum)
       call pawmkrhoij(atindx1,cprj_tmp,dimcprj,dtset%istwfk,dtset%kptopt,dtset%mband,mband_cprj,&
&       mcprj_tmp,dtset%mkmem,mpi_enreg,natom,dtset%nband,dtset%nkpt,dtset%nspinor,dtset%nsppol,&
&       occ,dtset%paral_kgb,paw_dmft,dtset%pawprtvol,pawrhoij,dtfil%unpaw,dtset%wtk)
       call cprj_free(cprj_tmp)
       ABI_DEALLOCATE(cprj_tmp)
     end if
     call timab(555,2,tsec)
!    Build symetrized packed rhoij and compensated pseudo density
     cplex=1;ipert=0;idir=0;qpt(:)=zero
     call  pawmkrho(compch_fft,cplex,gprimd,idir,psps%indlmn,indsym,ipert,psps%lmnmax,mpi_enreg,&
&     natom,dtset%nspden,dtset%nsym,ntypat,dtset%paral_kgb,pawang,pawfgr,pawfgrtab,&
&     dtset%pawprtvol,pawrhoij,pawtab,qpt,rhowfg,rhowfr,rhor,rprimd,dtset%symafm,symrec,&
&     dtset%typat,ucvol,xred,pawnhat=nhat,rhog=rhog)
   end if

!  Find and print minimum and maximum total electron density and locations
!  Compute density residual (if required) and its squared norm
   if (iscf>0) then
     if (psps%usepaw==0) then
       write(message,'(a)') ' vtorho: echo density'
       call wrtout(std_out,message,'COLL')
       call prtrhomxmn(std_out,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%nspden,1,rhor,ucvol=ucvol)
     else
       write(message,'(a)') ' vtorho: echo density (fine grid) '
       call wrtout(std_out,message,'COLL')
       call prtrhomxmn(std_out,mpi_enreg,nfftf,pawfgr%ngfft,dtset%nspden,1,rhor,ucvol=ucvol)
     end if
     if (dtset%prtdensph==1)then
       if (psps%usepaw==0) then
         call calcdensph(gmet,mpi_enreg,natom,dtset%nfft,dtset%ngfft,dtset%nspden,&
&         ntypat,std_out,dtset%ratsph,rhor,rprimd,dtset%typat,ucvol,xred)
       else
         call calcdensph(gmet,mpi_enreg,natom,nfftf,pawfgr%ngfft,dtset%nspden,&
&         ntypat,std_out,dtset%ratsph,rhor,rprimd,dtset%typat,ucvol,xred)
       end if
     end if

     if (optres==1) then
       nvresid=rhor-nvresid
       call sqnorm_v(1,mpi_enreg,nfftf,nres2,dtset%nspden,optres,nvresid)
     end if
   end if

 end if ! iscf>0 or iscf=-3

 if(psps%usepaw==1.and.(iscf>0.or.iscf==-3))  then
   ABI_DEALLOCATE(rhowfr)
   ABI_DEALLOCATE(rhowfg)
 end if

 call timab(994,2,tsec)

 if(iscf==-1)then

   call timab(995,1,tsec)
   call status(0,dtfil%filstat,iexit,level,'call tddft    ')

!  Eventually compute the excited states within tddft
   if (psps%usepaw==1) then
!    In case of PAW calculation, have to transfer kxc from the fine to the coarse grid:
     ABI_ALLOCATE(cgrkxc,(dtset%nfft,nkxc))
     do ikxc=1,nkxc
       call transgrid(1,mpi_enreg,1,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrkxc(:,ikxc),kxc(:,ikxc))
     end do
     call tddft(cg,dtfil,dtset,eigen,etotal,gmet,gprimd,gsqcut,&
&     kg,cgrkxc,dtset%mband,mgfftdiel,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%nfft,&
&     ngfftdiel,dtset%nkpt,nkxc,npwarr,dtset%nspinor,dtset%nsppol,occ,ucvol,wffnew)
     ABI_DEALLOCATE(cgrkxc)
   else
     call tddft(cg,dtfil,dtset,eigen,etotal,gmet,gprimd,gsqcut,&
&     kg,kxc,dtset%mband,mgfftdiel,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%nfft,&
&     ngfftdiel,dtset%nkpt,nkxc,npwarr,dtset%nspinor,dtset%nsppol,occ,ucvol,wffnew)
   end if
   call timab(995,2,tsec)

 else

   call timab(996,1,tsec)

!  Eventually compute the susceptibility matrix and the
!  dielectric matrix when istep_mix is equal to 1 or dielstrt
!  if( (istep_mix==1        .and. dielop>=2) .or. &
!  &     (istep_mix==dielstrt .and. dielop>=1) .or. &
!  &       computesusmat       )then
   call testsusmat(computesusmat,dielop,dielstrt,dtset,istep_mix) !test if the matrix is to be computed
   if(computesusmat) then
     dielar(1)=dtset%diecut;dielar(2)=dtset%dielng
     dielar(3)=dtset%diemac;dielar(4)=dtset%diemix
     dielar(5)=dtset%diegap;dielar(6)=dtset%dielam
     dielar(7)=dtset%diemix;if (iscf>=10) dielar(7)=dtset%diemixmag
     usetimerev=1
     if (psps%usepaw==1.and.dtset%pawspnorb>0.and.dtset%kptopt/=1.and.dtset%kptopt/=2) usetimerev=0

     call status(0,dtfil%filstat,iexit,level,'call suscep_st')
     neglect_pawhat=1-dtset%pawsushat
     call suscep_stat(atindx1,cg,cprj,dielar,&
&     dimcprj,doccde,eigen,gbound_diel,gprimd,&
&     irrzondiel,dtset%istwfk,kg,kg_diel,lmax_diel,&
&     dtset%mband,mcg,mcprj,mgfftdiel,dtset%mkmem,mpi_enreg,dtset%mpw,natom,dtset%nband,&
&     neglect_pawhat,nfftdiel,ngfftdiel,&
&     dtset%nkpt,npwarr,npwdiel,dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%nsym,ntypat,&
&     occ,dtset%occopt,pawang,pawtab,phnonsdiel,ph1ddiel,rprimd,&
&     susmat,dtset%symafm,dtset%symrel,dtset%tnons,dtset%typat,&
&     ucvol,dtfil%unkg,dtfil%unpaw,usecprj,psps%usepaw,usetimerev,wffnew,dtset%wtk,ylmdiel)
!    GMR
!    Print the susceptibility matrix
!    do isp1=1,dtset%nspden
!    do isp2=1,dtset%nspden
!    write(std_out,'(5x,a,2i2)') 'Susceptibility matrix for spins=',isp1,isp2
!    write(std_out,'(9x,a,13x,a,10x,a,10x,a)') "g","g'","real","imag"
!    do ipw1=1,10
!    do ipw2=ipw1,10
!    write(std_out,'(2x,3i4,2x,3i4,2x,f12.8,2x,f12.8)') &
!    &      kg_diel(1:3,ipw1),kg_diel(1:3,ipw2),&
!    &      susmat(1,ipw1,isp1,ipw2,isp2),susmat(2,ipw1,isp1,ipw2,isp2)
!    end do
!    end do
!    end do
!    end do

   end if
   call timab(996,2,tsec)

 end if ! end condition on iscf

 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(EigMin)

!PAW: deallocate <p|c> (cprj)
 if (psps%usepaw==1) then
   if(dtset%mkmem/=0.and.usecprj==1) then
     call cprj_free(cprj)
   end if
   usecprj=0
   ABI_DEALLOCATE(dimcprj)
   ABI_DEALLOCATE(cprj)
 end if

!deallocations in berryopt +/- 5 case
 if (abs(dtset%berryopt)==5) then
   call cprj_free(cprj_berry)
   ABI_DEALLOCATE(cprj_berry)
   call cprj_free(cprj_gat)
   ABI_DEALLOCATE(cprj_gat)
   ABI_DEALLOCATE(ikpt_recv)
   ABI_DEALLOCATE(dimcprj_u)
 end if

!Rotate labels of disk files when wf i/o is used
 if (dtset%mkmem==0) then
   wfftmp=wffnow ; wffnow=wffnew ; wffnew=wfftmp
 end if

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
   write(message,'(a1,a,a1,a,i1,a)') ch10,' vtorho : exit ',&
&   ch10,'  prtvol=-',level,', debugging mode => stop '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(980,2,tsec)

!DEBUG
!sum=0.0
!do jj=1,dtset%nspden
!do kk=1,nfftf
!sum=sum+abs(rhor(kk,jj))
!end do
!end do
!write(std_out,*) 'SUM(rhor)',sum
!write(std_out,*)' vtorho : exit, residm=',residm
!ENDDEBUG

end subroutine vtorho
!!***
