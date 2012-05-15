!{\src2tex{textfont=tt}}
!!****f* ABINIT/respfn
!! NAME
!! respfn
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations of Response functions.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG, DRH, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  codvsn=code version
!!  cpui=initial cpu time
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum single fft dimension
!!   | mkmem=maximum number of k points which can fit in core memory
!!   | mpw=maximum number of planewaves in basis sphere (large number)
!!   | natom=number of atoms in unit cell
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=number of channels for spin-polarization (1 or 2)
!!   | nsym=number of symmetry elements in space group
!!  mkmems(3)=array containing the tree values of mkmem (see above) (k-GS, k+q-GS and RF)
!!  mpi_enreg=informations about MPI parallelization
!!  npwtot(nkpt)=number of planewaves in basis and boundary at each k point
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  etotal=total energy (sum of 7 or 8 contributions) (hartree)
!!
!! SIDE EFFECTS
!!  iexit=index of "exit" on first line of file (0 if not found)
!!  occ(mband*nkpt*nsppol)=occup number for each band (often 2) at each k point
!!    Occupations number may have been read from a previous dataset...
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!    Some dimensions in pawrad have been set in driver.f
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!    Some dimensions in pawtab have been set in driver.f
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!    Before entering the first time in respfn, a significant part of psps
!!    has been initialized: the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,
!!    mpsso,mgrid,ntypat,n1xccc,usepaw,useylm, and the arrays dimensioned to npsp
!!    All the remaining components of psps are to be initialized in the call
!!    to pspini.  The next time the code enters respfn, psps might be identical
!!    to the one of the previous dtset, in which case, no reinitialisation
!!    is scheduled in pspini.f .
!!  results_respfn <type(results_respfn_type)>=stores some results of respfn calls
!!
!! NOTES
!! USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      alloc_hamilt_gpu,atm2fft,beginprtscphon,bstruct_clean,bstruct_init
!!      chkexi,chkpawovlp,chkph3,clnmpi_atom,clnmpi_fft,d2sym3
!!      dealloc_hamilt_gpu,destroy_paw_an,destroy_paw_ij,distrb2,dyfnl3,dyfro3
!!      dyout3,dyxc13,eigen_meandege,elph2_fanddw,eltfrhar3,eltfrkin3,eltfrloc3
!!      eltfrnl3,eltfrxc3,endprtscphon,ewald3,ewald4,fourdp,gath3,getcut,getph
!!      hdr_clean,hdr_init,hdr_update,init_paw_an,init_paw_ij,init_pawfgr
!!      initmpi_atom,initmpi_fft,initrhoij,initylmg,inwffil,ioarr,ioddb8_out
!!      irred_perts,kpgio,leave_new,loper3,mkcore,mklocl,mkrho,newocc,nhatgrid
!!      nullify_paw_an,nullify_paw_ij,pawdenpot,pawdij,pawexpiqr,pawfgrtab_free
!!      pawfgrtab_init,pawinit,pawmknhat,pawpuxinit,phfrq3,prteigrs,prtph3
!!      prtscphon,psddb8,pspini,q0dy3_apply,q0dy3_calc,rhohxc,rhoij_copy
!!      rhoij_free,setsym,setsymrhoij,setup1,status,symdij,symkchk,symq3,symzat
!!      syper3,timab,transgrid,wffclose,wffdelete,wings3,wrtloctens,wrtout
!!      xbarrier_mpi,xcast_mpi,xcomm_world,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine respfn(codvsn,cpui,dtfil,dtset,etotal,iexit,&
&  mkmems,mpi_enreg,npwtot,&
&  occ,pawang,pawrad,pawtab,psps,results_respfn,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_paw_toolbox
 use m_wffile
 use m_errors
 use m_prtscphon
 use m_results_respfn
 use m_paw_dmft, only: paw_dmft_type
 use m_header,   only : hdr_init, hdr_clean
 use m_ebands,   only : bstruct_init, bstruct_clean

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'respfn'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
#if defined HAVE_GPU_CUDA
 use interfaces_52_manage_cuda
#endif
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_56_xc
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_62_occeig
 use interfaces_65_psp
 use interfaces_66_paw
 use interfaces_67_common
 use interfaces_72_response
 use interfaces_79_seqpar_mpi
 use interfaces_95_drive, except_this_one => respfn
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(inout) :: iexit
 real(dp),intent(in) :: cpui
 real(dp),intent(out) :: etotal
 character(len=6),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
 integer,intent(in) :: mkmems(3)
 integer,intent(inout) :: npwtot(dtset%nkpt)
 real(dp),intent(inout) :: xred(3,dtset%natom)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(results_respfn_type),intent(inout) :: results_respfn

!Local variables-------------------------------
!---- Local variables : integer scalars
!Define file format for different type of files. Presently,
!only one file format is supported for each type of files, but this might
!change soon ...
!1   for wavefunction file, old format (version prior to 2.0)
!2   for wavefunction file, new format (version 2.0 and after)    (fform)
!(51 or 52   for density rho(r)       (fformr)
!101 or 102 for potential V(r) file. (fformv)
 integer,parameter :: fform=2,fformv=102,formeig=0,level=10
 integer,parameter :: response=1,syuse=0
 integer,save :: nsym_old=-1
 integer,save :: paw_gencond(6)=(/-1,-1,-1,-1,-1,-1/)
 integer :: fformr=52
 integer :: nk3xc
 integer :: accessfil,analyt,ask_accurate,band_index,bantot,choice,cplex,cplex_dij
 integer :: dim_eig2nkq,dim_eigbrd,dyfr_cplex,dyfr_nondiag,fullinit
 integer :: gscase,has_dijso,iatom,iband,idir,ider,ierr,ifft,ii,ikpt,indx
 integer :: initialized,ipert,ipert2,ir,ireadwf0,iscf,iscf_eff,ispden,isppol
 integer :: itypat,izero,master,mcg,me,mgfftf,mk1mem,mkqmem,mpert,mpsang,mu
 integer :: natom,n3xccc,nband_k,nblok,nfftf,nfftot,nhatdim,nhatgrdim
 integer :: nkpt_eff,nkpt_max,nkxc,nspden_rhoij,ntypat,nzlmopt,openexit
 integer :: optcut,option,optgr0,optgr1,optgr2,optorth,optrad
 integer :: optatm,optdyfr,optgr,optn,optn2,optstr,optv
 integer :: outd2,pawbec,prtbbb,psp_gencond,qzero,rdwr,rdwrpaw
 integer :: rfasr,rfddk,rfelfd,rfphon,rfstrs,rfuser
 integer :: spaceworld,sumg0,tim_mkrho,timrev,usecprj,usexcnhat,v_size,vrsddb
 logical :: qeq0
 real(dp) :: boxcut,compch_fft,compch_sph,cpus,ecore,ecut_eff,ecutdg_eff,ecutf,eei,eew,ehart,eii,ek,enl,entropy,enxc
 real(dp) :: epaw,epawdc,etot,fermie,gsqcut,gsqcut_eff,gsqcutc_eff,qphnrm,residm
 real(dp) :: tolwfr
 real(dp) :: ucvol,vxcavg
 character(len=fnlen) :: dscrpt
 character(len=fnlen) :: phonon_freq_filename, phonon_vec_filename
 character(len=500) :: message
 type(bandstructure_type) :: bstruct
 type(hdr_type) :: hdr
 type(paw_dmft_type) :: paw_dmft
 type(pawfgr_type) :: pawfgr
 type(wffile_type) :: wffgs,wfftgs
 type(wvl_data) :: wvl
 integer :: ddkfil(3),ngfft(18),ngfftf(18),rfdir(3)
 integer,allocatable :: atindx(:),atindx1(:),blkflg(:,:,:,:),blkflgfrx1(:,:,:,:),blkflg1(:,:,:,:)
 integer,allocatable :: blkflg2(:,:,:,:),carflg(:,:,:,:),dimcprj(:),indsym(:,:,:)
 integer,allocatable :: irrzon(:,:,:),kg(:,:),l_size_atm(:),nattyp(:),npwarr(:)
 integer,allocatable :: pertsy(:,:),rfpert(:),symq(:,:,:),symrec(:,:,:)
 real(dp) :: dummy6(6),gmet(3,3),gprimd(3,3),qphon(3)
 real(dp) :: rmet(3,3),rprimd(3,3),strsxc(6),tsec(2)
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp),allocatable :: amass(:),becfrnl(:,:,:),cg(:,:),d2bbb(:,:,:,:,:,:),d2cart(:,:,:,:,:)
 real(dp),allocatable :: d2cart_bbb(:,:,:,:,:,:),d2eig0(:,:,:,:,:)
 real(dp),allocatable :: d2k0(:,:,:,:,:),d2lo(:,:,:,:,:),d2loc0(:,:,:,:,:)
 real(dp),allocatable :: d2matr(:,:,:,:,:),d2nfr(:,:,:,:,:),d2nl(:,:,:,:,:),d2ovl(:,:,:,:,:)
 real(dp),allocatable :: d2nl0(:,:,:,:,:),d2nl1(:,:,:,:,:),d2tmp(:,:,:,:,:),d2vn(:,:,:,:,:)
 real(dp),allocatable :: displ(:),doccde(:),dyew(:,:,:,:,:)
 real(dp),allocatable :: dum_gauss(:),dum_dyfrn(:),dum_dyfrv(:)
 real(dp),allocatable :: dum_grn(:),dum_grv(:),dum_rhog(:),dum_vg(:)
 real(dp),allocatable :: dyewq0(:,:,:),dyfrlo(:,:,:),dyfrlo_indx(:,:,:)
 real(dp),allocatable :: dyfrnl(:,:,:,:,:),dyfrwf(:,:,:,:,:),dyfrx1(:,:,:,:,:)
 real(dp),allocatable :: dyfrx2(:,:,:),eigen0(:),eigval(:),eigvec(:)
 real(dp),allocatable :: eig2nkq(:,:,:,:,:,:,:),eigbrd(:,:,:,:,:,:,:)
 real(dp),allocatable :: eigen_fan(:),eigen_ddw(:),eigen_fanddw(:)
 real(dp),allocatable :: eigen_fan_mean(:),eigen_ddw_mean(:)
 real(dp),allocatable :: eltcore(:,:),elteew(:,:),eltfrhar(:,:),eltfrkin(:,:)
 real(dp),allocatable :: eltfrloc(:,:),eltfrnl(:,:),eltfrxc(:,:),grtn_indx(:,:)
 real(dp),allocatable :: grxc(:,:),kxc(:,:),nhat(:,:),nhatgr(:,:,:)
 real(dp),allocatable :: ph1d(:,:),ph1df(:,:),phfrq(:),phnons(:,:,:)
 real(dp),allocatable :: rhog(:,:),rhor(:,:),rhowfg(:,:),rhowfr(:,:)
 real(dp),allocatable :: vhartr(:),vpsp(:),vtrial(:,:)
 real(dp),allocatable :: vxc(:,:),work(:),xccc3d(:),ylm(:,:),ylmgr(:,:,:)
 type(cprj_type),allocatable :: cprj_dum(:,:)
 type(paw_an_type),allocatable :: paw_an(:)
 type(paw_ij_type),allocatable :: paw_ij(:)
 type(pawfgrtab_type),allocatable,save :: pawfgrtab(:)
 type(pawrhoij_type),allocatable :: pawrhoij(:)

! local GIPAW variables
! real(dp),allocatable :: cs(:,:,:),gcart(:,:,:,:),jvec(:,:,:)

! ***********************************************************************

 DBG_ENTER("COLL")

 call timab(132,1,tsec)
 call timab(133,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter         ')
!
!Print the list of irreducible perturbations then exit.
 if (dtset%prtvol==-level) then
   call irred_perts(Dtset,ab_out)
   call leave_new('COLL')
 end if

 mpi_enreg%paralbd=1
 mpi_enreg%me_fft=0
 mpi_enreg%nproc_fft=1
 mpi_enreg%paral_fft=0
 mpi_enreg%paral_level=2
 mpi_enreg%flag_ind_kg_mpi_to_seq=0
 if(mpi_enreg%paral_compil==1) then
   ABI_ALLOCATE(mpi_enreg%proc_distrb,(dtset%nkpt,dtset%mband,dtset%nsppol))
   call distrb2(dtset%mband, dtset%nband, dtset%nkpt, dtset%nsppol, mpi_enreg)
 end if
 nkpt_max=50;if(mpi_enreg%paral_compil==1)nkpt_max=-1
 call initmpi_fft(dtset,mpi_enreg)
 call initmpi_atom(dtset,mpi_enreg)

!Define FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 call init_pawfgr(dtset,pawfgr,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfft,ngfftf)

!If dtset%accesswff == 2 set all array outputs to netcdf format
 accessfil = 0
 if (dtset%accesswff == IO_MODE_NETCDF) then
   accessfil = 1
 end if
 if (dtset%accesswff == IO_MODE_ETSF) then
   accessfil = 3
 end if

!Structured debugging if dtset%prtvol==-level
 if(dtset%prtvol==-level)then
   write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&   ' respfn : enter , debug mode '
   call wrtout(std_out,message,'COLL')
 end if

!Option input variables
 iscf=dtset%iscf

!Respfn input variables
 rfasr=dtset%rfasr   ; rfdir(1:3)=dtset%rfdir(1:3)
 rfddk=dtset%rfddk   ; rfelfd=dtset%rfelfd
 rfphon=dtset%rfphon
 rfstrs=dtset%rfstrs ; rfuser=dtset%rfuser

!mkmem variables (mkmem is already argument)
 mkqmem=mkmems(2) ; mk1mem=mkmems(3)

 ntypat=psps%ntypat
 natom=dtset%natom

 call status(0,dtfil%filstat,iexit,level,'call setup1   ')

 ecore=zero

!LIKELY TO BE TAKEN AWAY
 initialized=0
 ek=zero ; ehart=zero ; enxc=zero ; eei=zero ; enl=zero
 eii=zero ; eew=zero

!Set up for iterations
 ABI_ALLOCATE(amass,(natom))
 call setup1(dtset%acell_orig(1:3,1),amass,bantot,dtset,&
& ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,&
& natom,ngfftf,ngfft,dtset%nkpt,dtset%nsppol,&
& response,rmet,dtset%rprim_orig(1:3,1:3,1),rprimd,ucvol,psps%usepaw)

!Define the set of admitted perturbations
 mpert=natom+6

!Initialize the list of perturbations rfpert
 ABI_ALLOCATE(rfpert,(mpert))
 rfpert(:)=0
 if(rfphon==1)rfpert(dtset%rfatpol(1):dtset%rfatpol(2))=1

 if(rfddk==1)rfpert(natom+1)=1
 if(rfddk==2)rfpert(natom+6)=1

 if(rfelfd==1.or.rfelfd==2)rfpert(natom+1)=1
 if(rfelfd==1.or.rfelfd==3)rfpert(natom+2)=1

 if(rfstrs==1.or.rfstrs==3)rfpert(natom+3)=1
 if(rfstrs==2.or.rfstrs==3)rfpert(natom+4)=1

 if(rfuser==1.or.rfuser==3)rfpert(natom+5)=1
 if(rfuser==2.or.rfuser==3)rfpert(natom+6)=1

 qeq0=(dtset%qptn(1)**2+dtset%qptn(2)**2+dtset%qptn(3)**2<1.d-14)

!Init spaceworld
 call xcomm_world(mpi_enreg,spaceworld)

!Default for sequential use
 master=0
 me = xcomm_rank(spaceworld)

 call status(0,dtfil%filstat,iexit,level,'call kpgio(1) ')

!Set up the basis sphere of planewaves
 ABI_ALLOCATE(kg,(3,dtset%mpw*dtset%mkmem))
 ABI_ALLOCATE(npwarr,(dtset%nkpt))
 call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg,dtfil%fnametmp_kg,&
& dtset%kptns,dtset%mkmem,dtset%nband,dtset%nkpt,'PERS',mpi_enreg,dtset%mpw,npwarr,npwtot,&
& dtset%nsppol,dtfil%unkg)

!Set up the Ylm for each k point
 mpsang=psps%mpsang
 ABI_ALLOCATE(ylm,(dtset%mpw*dtset%mkmem,mpsang*mpsang*psps%useylm))
 if (rfstrs/=0) then
   ABI_ALLOCATE(ylmgr,(dtset%mpw*dtset%mkmem,9,mpsang*mpsang*psps%useylm))
 else
   ABI_ALLOCATE(ylmgr,(1,1,psps%useylm))
 end if
 if (psps%useylm==1) then
   if(dtset%mkmem==0) open(dtfil%unylm,file=dtfil%fnametmp_ylm,form='unformatted',status='unknown')
   call status(0,dtfil%filstat,iexit,level,'call initylmg ')
   option=0
   if (rfstrs/=0) option=2
   call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&   npwarr,dtset%nsppol,option,rprimd,dtfil%unkg,dtfil%unylm,ylm,ylmgr)
 end if

 call timab(133,2,tsec)
 call timab(134,1,tsec)

!Open and read pseudopotential files
 call status(0,dtfil%filstat,iexit,level,'call pspini(1)')
 call pspini(dtset,ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,level,&
& pawrad,pawtab,psps,rprimd)

 call timab(134,2,tsec)
 call timab(135,1,tsec)

!Initialize band structure datatype
 ABI_ALLOCATE(doccde,(bantot))
 ABI_ALLOCATE(eigen0,(bantot))
 doccde(:)=zero ; eigen0(:)=zero
 call bstruct_init(bantot,bstruct,dtset%nelect,doccde,eigen0,dtset%istwfk,dtset%kptns,&
& dtset%nband,dtset%nkpt,npwarr,dtset%nsppol,dtset%nspinor,dtset%tphysel,dtset%tsmear,&
& dtset%occopt,dtset%occ_orig,dtset%wtk)
 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(eigen0)

!Initialize PAW atomic occupancies
 if (psps%usepaw==1) then
   ABI_ALLOCATE(pawrhoij,(natom))
   nspden_rhoij=dtset%nspden;if (dtset%pawspnorb>0.and.dtset%nspinor==2) nspden_rhoij=4
   call initrhoij(dtset%pawcpxocc,psps%indlmn,dtset%lexexch,psps%lmnmax,dtset%lpawu,mpi_enreg,&
&   natom,mpi_enreg%natom,nspden_rhoij,dtset%nspinor,dtset%nsppol,dtset%ntypat,&
&   pawrhoij,pawtab,dtset%spinat,dtset%typat)
 end if

!Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,hdr,pawtab,gscase,psps,wvl%descr)

!Update header, with evolving variables, when available
!Here, rprimd, xred and occ are available
 etot=hdr%etot ; fermie=hdr%fermie ; residm=hdr%residm
 call hdr_update(bantot,etot,fermie,hdr,natom,&
& residm,rprimd,occ,mpi_enreg,pawrhoij,psps%usepaw,xred)

!Clean band structure datatype (should use it more in the future !)
 call bstruct_clean(bstruct)

 call status(0,dtfil%filstat,iexit,level,'call inwffil(1')

!Initialize wavefunction files and wavefunctions.
 ireadwf0=1

 mcg=dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
 ABI_ALLOCATE(cg,(2,mcg))
 ABI_ALLOCATE(eigen0,(dtset%mband*dtset%nkpt*dtset%nsppol))
 eigen0(:)=zero ; ask_accurate=1
 optorth=1;if (psps%usepaw==1.and.ireadwf0==1) optorth=0

 call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen0,dtset%exchn2n3d,&
& formeig,gmet,hdr,ireadwf0,dtset%istwfk,kg,dtset%kptns,&
& dtset%localrdwf,dtset%mband,mcg,dtset%mkmem,mpi_enreg,dtset%mpw,&
& dtset%nband,ngfft,dtset%nkpt,npwarr,dtset%nsppol,dtset%nsym,&
& occ,optorth,rprimd,dtset%symafm,dtset%symrel,dtset%tnons,&
& dtfil%unkg,wffgs,wfftgs,dtfil%unwffgs,dtfil%unwftgs,dtfil%fnamewffk,dtfil%fnametmp_wfgs,wvl)

 if (psps%usepaw==1.and.ireadwf0==1) then
   call rhoij_copy(hdr%pawrhoij,pawrhoij,mpi_enreg=mpi_enreg)
 end if

 call timab(135,2,tsec)
 call timab(136,1,tsec)

!Close wffgs, if it was ever opened (in inwffil)
 if (ireadwf0==1) then
   call WffClose(wffgs,ierr)
 end if

!Report on eigen0 values   ! Should use prteigrs.F90
 write(message, '(a,a)' )ch10,' respfn : eigen0 array'
 call wrtout(std_out,message,'COLL')
 nkpt_eff=dtset%nkpt
 if( (dtset%prtvol==0.or.dtset%prtvol==1.or.dtset%prtvol==2) .and. dtset%nkpt>nkpt_max ) nkpt_eff=nkpt_max
 band_index=0
 do isppol=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     if(ikpt<=nkpt_eff)then
       write(message, '(a,i2,a,i5)' )&
&       '  isppol=',isppol,', k point number',ikpt
       call wrtout(std_out,message,'COLL')
       do iband=1,nband_k,4
         write(message, '(a,4es16.6)')&
&         '  ',eigen0(iband+band_index:min(iband+3,nband_k)+band_index)
         call wrtout(std_out,message,'COLL')
       end do
     else if(ikpt==nkpt_eff+1)then
       write(message,'(a,a)' )&
&       '  respfn : prtvol=0, 1 or 2, stop printing eigen0.',ch10
       call wrtout(std_out,message,'COLL')
     end if
     band_index=band_index+nband_k
   end do
 end do

!Allocation for forces and atomic positions (should be taken away, also argument ... )
 ABI_ALLOCATE(grxc,(3,natom))

!Do symmetry stuff
 call status(0,dtfil%filstat,iexit,level,'call setsym   ')
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 ABI_ALLOCATE(irrzon,(nfftot**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(phnons,(2,nfftot**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(indsym,(4,dtset%nsym,natom))
 ABI_ALLOCATE(symrec,(3,3,dtset%nsym))
 iscf_eff=0
!If the density is to be computed by mkrho, need irrzon and phnons
 if(dtset%getden==0)iscf_eff=1
 call setsym(indsym,irrzon,iscf_eff,natom,&
& nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
& phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

!Symmetrize atomic coordinates over space group elements:
 call status(0,dtfil%filstat,iexit,level,'call symzat   ')
 call symzat(indsym,natom,dtset%nsym,dtset%symrel,dtset%tnons,xred)

!Examine the symmetries of the q wavevector
 ABI_ALLOCATE(symq,(4,2,dtset%nsym))
 timrev=1
 call status(0,dtfil%filstat,iexit,level,'call symq3    ')
 call symq3(dtset%nsym,dtset%qptn,symq,symrec,timrev)

!Generate an index table of atoms, in order for them to be used
!type after type.
 ABI_ALLOCATE(atindx,(natom))
 ABI_ALLOCATE(atindx1,(natom))
 ABI_ALLOCATE(nattyp,(ntypat))
 indx=1
 do itypat=1,ntypat
   nattyp(itypat)=0
   do iatom=1,natom
     if(dtset%typat(iatom)==itypat)then
       atindx(iatom)=indx
       atindx1(indx)=iatom
       indx=indx+1
       nattyp(itypat)=nattyp(itypat)+1
     end if
   end do
 end do

!Here allocation of GPU for fft calculations
#if defined HAVE_GPU_CUDA
 if (dtset%use_gpu_cuda==1) then
   call alloc_hamilt_gpu(atindx1,dtset,gprimd,mpi_enreg,nattyp,0,psps,dtset%use_gpu_cuda)
 end if
#endif

!Compute structure factor phases for current atomic pos:
 ABI_ALLOCATE(ph1d,(2,3*(2*dtset%mgfft+1)*natom))
 ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*natom))
 call status(0,dtfil%filstat,iexit,level,'call getph    ')
 call getph(atindx,natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)
 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
   call getph(atindx,natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
 else
   ph1df(:,:)=ph1d(:,:)
 end if

!Compute occupation numbers and fermi energy, in case
!occupation scheme is metallic.
 ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
 if( dtset%occopt>=3.and.dtset%occopt<=8 ) then
   call status(0,dtfil%filstat,iexit,level,'call newocc   ')
   call newocc(doccde,eigen0,entropy,fermie,dtset%fixmom,dtset%mband,dtset%nband,&
&   dtset%nelect,dtset%nkpt,dtset%nspinor,dtset%nsppol,occ,dtset%occopt,dtset%prtvol,dtset%stmbias,&
&   dtset%tphysel,dtset%tsmear,dtset%wtk)

!  Update fermie and occ
   etot=hdr%etot ; residm=hdr%residm
   call hdr_update(bantot,etot,fermie,hdr,natom,&
&   residm,rprimd,occ,mpi_enreg,pawrhoij,psps%usepaw,xred)

 else
!  doccde is irrelevant in this case
   doccde(:)=zero

 end if

!Recompute first large sphere cut-off gsqcut,
!without taking into account dilatmx
 if (psps%usepaw==1) then
   write(message,'(2a)') ch10,' FFT (fine) grid used in SCF cycle:'
   call wrtout(std_out,message,'COLL')
 end if
 ecutf=dtset%ecut;if (pawfgr%usefinegrid==1) ecutf=dtset%pawecutdg
 call status(0,dtfil%filstat,iexit,level,'call getcut   ')
 call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,std_out,k0,ngfftf)

!PAW: 1- Initialize values for several arrays depending only on atomic data
!2- Check overlap
!3- Identify FFT points in spheres and compute g_l(r).Y_lm(r) (and exp(-i.q.r) if needed)
!4- Allocate PAW specific arrays
!5- Compute perturbed local potential inside spheres
!6- Eventually open temporary storage files
 if(psps%usepaw==1) then
!  1-Initialize values for several arrays depending only on atomic data
   if (psp_gencond==1.or.&
&   paw_gencond(1)/=dtset%pawlcutd .or.paw_gencond(2)/=dtset%pawlmix  .or.&
&   paw_gencond(3)/=dtset%pawnphi  .or.paw_gencond(4)/=dtset%pawntheta.or.&
&   paw_gencond(5)/=dtset%pawspnorb.or.paw_gencond(6)/=dtset%pawxcdev) then
     call timab(553,1,tsec)
     call pawinit(0._dp,psps%indlmn,dtset%pawlcutd,dtset%pawlmix,psps%lmnmax,psps%mpsang,&
&     dtset%pawnphi,dtset%nsym,dtset%pawntheta,psps%ntypat,&
&     pawang,pawrad,dtset%pawspnorb,pawtab,dtset%pawxcdev)
     paw_gencond(1)=dtset%pawlcutd ; paw_gencond(2)=dtset%pawlmix
     paw_gencond(3)=dtset%pawnphi  ; paw_gencond(4)=dtset%pawntheta
     paw_gencond(5)=dtset%pawspnorb; paw_gencond(6)=dtset%pawxcdev
     call timab(553,2,tsec)
   else
     if (pawtab(1)%has_kij  ==1) pawtab(1:psps%ntypat)%has_kij  =2
     if (pawtab(1)%has_nabla==1) pawtab(1:psps%ntypat)%has_nabla=2
   end if
   if (psp_gencond==1.or.nsym_old/=dtset%nsym) then
     call setsymrhoij(gprimd,pawang%l_max-1,dtset%nsym,dtset%pawprtvol,&
&     rprimd,symrec,pawang%zarot)
     nsym_old=dtset%nsym
   end if
   psps%n1xccc=maxval(pawtab(1:psps%ntypat)%usetcore)
   pawtab(:)%usepawu=0
   pawtab(:)%useexexch=0
   pawtab(:)%exchmix=zero
!  if (dtset%usepawu>0.or.dtset%useexexch>0) then
   call pawpuxinit(dtset%dmatpuopt,dtset%exchmix,dtset%jpawu,dtset%lexexch,dtset%lpawu,&
&   psps%indlmn,psps%lmnmax,ntypat,pawang,dtset%pawprtvol,pawrad,pawtab,dtset%upawu,&
&   dtset%usedmft,dtset%useexexch,dtset%usepawu)
!  end if
   compch_fft=-1.d5;compch_sph=-1.d5
   usexcnhat=maxval(pawtab(:)%usexcnhat)
   usecprj=dtset%pawusecp
!  2-Check overlap
   call status(0,dtfil%filstat,iexit,level,'call chkpawovlp')
   call chkpawovlp(natom,psps%ntypat,dtset%pawovlp,pawtab,rmet,dtset%typat,xred)
!  3-Identify FFT points in spheres and compute g_l(r).Y_lm(r) and exp(-i.q.r)
   ABI_ALLOCATE(pawfgrtab,(natom))
   ABI_ALLOCATE(l_size_atm,(natom))
   do iatom=1,natom
     l_size_atm(iatom)=pawtab(dtset%typat(iatom))%lcut_size
   end do
   call pawfgrtab_init(pawfgrtab,1,l_size_atm,pawrhoij(1)%nspden)
   ABI_DEALLOCATE(l_size_atm)
   optcut=0;optgr0=dtset%pawstgylm;optgr1=0;optgr2=0;optrad=1-dtset%pawstgylm
   if (dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0) optgr1=dtset%pawstgylm
   if (rfphon==1) then
     if (optgr1==0) optgr1=dtset%pawstgylm
     if (optgr2==0) optgr2=dtset%pawstgylm
     if(optrad==0.and.(.not.qeq0)) optrad=1
   end if
   if (rfelfd==1.or.rfelfd==3) then
     if (optgr1==0) optgr1=dtset%pawstgylm
   end if
   call status(0,dtfil%filstat,iexit,level,'call nhatgrid ')
   call nhatgrid(atindx1,gmet,mpi_enreg,natom,natom,nattyp,ngfftf,psps%ntypat,&
&   optcut,optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,ucvol,xred)
   do iatom=1,natom
     ii=iatom
     if (mpi_enreg%nproc_atom>1) ii=mpi_enreg%atom_indx(iatom)
     call pawexpiqr(gprimd,pawfgrtab(iatom),dtset%qptn,xred(:,ii))
   end do
!  4-Allocate PAW specific arrays
   ABI_ALLOCATE(paw_an,(natom))
   ABI_ALLOCATE(paw_ij,(natom))
   call nullify_paw_an(paw_an)
   call nullify_paw_ij(paw_ij)
   cplex_dij=max(1+dtset%pawspnorb,dtset%nspinor)
   has_dijso=0;if (dtset%pawspnorb>0) has_dijso=1
   call init_paw_an(natom,dtset%ntypat,0,dtset%nspden,1,dtset%pawxcdev,dtset%typat,&
&   pawang,pawtab,paw_an)
   call init_paw_ij(paw_ij,1,cplex_dij,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%pawspnorb,natom,&
&   dtset%ntypat,dtset%typat,pawtab,has_dij=1,has_exexch_pot=1,has_pawu_occ=1)
   ABI_ALLOCATE(dimcprj,(natom))
   do iatom=1,natom
     dimcprj(atindx(iatom))=pawtab(dtset%typat(iatom))%lmn_size ! Be careful : ordered by atom-type
   end do
!  5-Eventually open temporary storage files
   if(dtset%mkmem==0) then
     open(dtfil%unpaw ,file=dtfil%fnametmp_paw,form='unformatted',status='unknown')
     rewind(unit=dtfil%unpaw)
   end if

 else ! PAW vs NCPP
   usexcnhat=0;usecprj=0
 end if

 ABI_ALLOCATE(rhog,(2,nfftf))
 ABI_ALLOCATE(rhor,(nfftf,dtset%nspden))

!Read ground-state charge density from diskfile in case getden /= 0
!or compute it from wfs that were
!read previously : rhor as well as rhog

 if (dtset%getden /= 0) then

   if (me == 0) then
     rdwr=1;rdwrpaw=psps%usepaw;if(ireadwf0/=0) rdwrpaw=0
     call status(0,dtfil%filstat,iexit,level,'call ioarr    ')
     call ioarr (accessfil,rhor, dtset, etotal,fformr,dtfil%fildensin,hdr, mpi_enreg, &
&     nfftf,pawrhoij,rdwr,rdwrpaw,wvl%descr)
     if (psps%usepaw==1.and.rdwrpaw/=0) then
       call hdr_update(bantot,etot,fermie,hdr,natom,&
&       residm,rprimd,occ,mpi_enreg,pawrhoij,psps%usepaw,xred)
     end if
   end if

   call xbarrier_mpi(spaceworld)
   call xcast_mpi(rhor,0,spaceworld,ierr)

!  Compute up+down rho(G) by fft
   call status(0,dtfil%filstat,iexit,level,'call fourdp   ')
   ABI_ALLOCATE(work,(nfftf))
   work(:)=rhor(:,1)
   call fourdp(1,rhog,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
   ABI_DEALLOCATE(work)

 else
   izero=0
!  Obtain the charge density from read wfs
!  Be careful: in PAW, compensation density has to be added !
   call status(0,dtfil%filstat,iexit,level,'call mkrho    ')
   tim_mkrho=4
   paw_dmft%use_dmft=0 ! respfn with dmft not implemented
   if (psps%usepaw==1) then
     ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
     ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
     call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&     mpi_enreg,npwarr,occ,paw_dmft,phnons,rhowfg,rhowfr,rprimd,tim_mkrho,ucvol,&
&     dtfil%unkg,wfftgs,wvl%wfs,wvl%descr)
     call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,rhog,rhowfr,rhor)
     ABI_DEALLOCATE(rhowfg)
     ABI_DEALLOCATE(rhowfr)
   else
     call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&     mpi_enreg,npwarr,occ,paw_dmft,phnons,rhog,rhor,rprimd,tim_mkrho,ucvol,&
&     dtfil%unkg,wfftgs,wvl%wfs,wvl%descr)
   end if

 end if    ! getden

!In PAW, compensation density has eventually to be added
 nhatgrdim=0;nhatdim=0
 if (psps%usepaw==1.and. &
& ((usexcnhat==0).or.(dtset%getden==0).or.dtset%xclevel==2)) then
   nhatdim=1
   ABI_ALLOCATE(nhat,(nfftf,dtset%nspden))
   call timab(558,1,tsec)
   nhatgrdim=0;if (dtset%xclevel==2.and.dtset%pawnhatxc>0) nhatgrdim=usexcnhat
   ider=2*nhatgrdim
   if (nhatgrdim>0)  then
     ABI_ALLOCATE(nhatgr,(nfftf,dtset%nspden,3))
   end if
   ider=0;izero=0;cplex=1;ipert=0;idir=0;qphon(:)=zero
   call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,mpi_enreg,natom,natom,&
&   nfftf,ngfftf,nhatgrdim,dtset%nspden,psps%ntypat,dtset%paral_kgb,pawang,pawfgrtab,&
&   nhatgr,nhat,pawrhoij,pawrhoij,pawtab,qphon,rprimd,ucvol,xred)
   if (dtset%getden==0) then
     rhor(:,:)=rhor(:,:)+nhat(:,:)
     call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
   end if
   call timab(558,2,tsec)
 end if

!The GS irrzon and phnons were only needed to symmetrize the GS density
 ABI_DEALLOCATE(irrzon)
 ABI_DEALLOCATE(phnons)

!Debugging : print the different parts of rhor
 if(dtset%prtvol==-level)then
   write(message,'(a)') '   ir     rhor(ir)     '
   call wrtout(std_out,message,'COLL')
   do ir=1,nfftf
     if(ir<=11 .or. mod(ir,301)==0 )then
       write(message,'(i5,a,es13.6)')ir,' ',rhor(ir,1)
       call wrtout(std_out,message,'COLL')
       if(dtset%nspden==2)then
         write(message,'(a,es13.6)')'      ',rhor(ir,2)
         call wrtout(std_out,message,'COLL')
       end if
     end if
   end do
 end if
 write(std_out,'(a)')' ' ! needed to make ibm6_xlf12 pass tests. No idea why this works. JWZ 5 Sept 2011
!Will compute now the total potential

!Compute local ionic pseudopotential vpsp
!and core electron density xccc3d:
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=nfftf
 ABI_ALLOCATE(xccc3d,(n3xccc))
 ABI_ALLOCATE(vpsp,(nfftf))
 if (psps%usepaw==1) then
!  PAW: compute Vloc and core charge together in reciprocal space
   call timab(562,1,tsec)
   optatm=1;optdyfr=0;optgr=0;optstr=0;optv=1;optn=n3xccc/nfftf;optn2=1
   call atm2fft(atindx1,xccc3d,vpsp,dum_dyfrn,dum_dyfrv,eei,dum_gauss,gmet,gprimd,&
&   dum_grn,dum_grv,gsqcut,&
&   mgfftf,mpi_enreg,psps%mqgrid_vl,natom,nattyp,nfftf,ngfftf,ntypat,&
&   optatm,optdyfr,optgr,optn,optn2,optstr,optv,dtset%paral_kgb,&
&   pawtab,ph1df,psps%qgrid_vl,dtset%qprtrb,dum_rhog,dummy6,dummy6,&
&   ucvol,psps%usepaw,dum_vg,dtset%vprtrb,psps%vlspl)
   call timab(562,2,tsec)
 else
!  Norm-cons.: compute Vloc in reciprocal space and core charge in real space
   option=1
   ABI_ALLOCATE(dyfrlo_indx,(3,3,natom))
   ABI_ALLOCATE(grtn_indx,(3,natom))
   call status(0,dtfil%filstat,iexit,level,'call mklocl   ')
   call mklocl(dtset,dyfrlo_indx,eei,gmet,gprimd,&
&   grtn_indx,gsqcut,dummy6,mgfftf,mpi_enreg,natom,nattyp,&
&   nfftf,ngfftf,dtset%nspden,ntypat,option,ph1df,psps,&
&   dtset%qprtrb,rhog,rhor,rprimd,ucvol,dtset%vprtrb,vpsp,wvl%descr,xred)
   ABI_DEALLOCATE(dyfrlo_indx)
   ABI_DEALLOCATE(grtn_indx)
   if (psps%n1xccc/=0) then
     call status(0,dtfil%filstat,iexit,level,'call mkcore   ')
     ABI_ALLOCATE(dyfrx2,(3,3,natom))
     call mkcore(dummy6,dyfrx2,grxc,mpi_enreg,natom,nfftf,dtset%nspden,ntypat,&
&     ngfftf(1),psps%n1xccc,ngfftf(2),ngfftf(3),option,rprimd,dtset%typat,ucvol,vxc,&
&     psps%xcccrc,psps%xccc1d,xccc3d,xred)
     ABI_DEALLOCATE(dyfrx2)
   end if
 end if

!Set up hartree and xc potential. Compute kxc here.
 option=2
 nkxc=2*min(dtset%nspden,2)-1
 if(dtset%xclevel==2)nkxc=23
 ABI_ALLOCATE(kxc,(nfftf,nkxc))
 ABI_ALLOCATE(vhartr,(nfftf))
 ABI_ALLOCATE(vxc,(nfftf,dtset%nspden))
 call status(0,dtfil%filstat,iexit,level,'call rhohxc   ')

!adjusted to call rhohxc
 nk3xc=1

 call rhohxc(dtset,enxc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,ngfftf,&
& nhat,nhatdim,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,option,rhog,rhor,&
& rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d)

!Compute local + Hxc potential, and subtract mean potential.
 ABI_ALLOCATE(vtrial,(nfftf,dtset%nspden))
 do ispden=1,min(dtset%nspden,2)
   do ifft=1,nfftf
     vtrial(ifft,ispden)=vhartr(ifft)+vxc(ifft,ispden)+vpsp(ifft)
   end do
 end do
 if (dtset%nspden==4) then
   do ispden=3,4
     do ifft=1,nfftf
       vtrial(ifft,ispden)=vxc(ifft,ispden)
     end do
   end do
 end if
 ABI_DEALLOCATE(vhartr)

 call status(0,dtfil%filstat,iexit,level,'end respfn(1) ')

 if(dtset%prtvol==-level)then
   write(message,'(a,a)') ch10,&
&   ' respfn : ground-state density and potential set up. '
   call wrtout(std_out,message,'COLL')
 end if

!PAW: compute Dij quantities (psp strengths)
 if (psps%usepaw==1)then
   option=1;nzlmopt=0;if (dtset%pawnzlm>0) nzlmopt=-1
   do iatom=1,natom
     itypat=dtset%typat(iatom)
     v_size=paw_an(iatom)%lm_size;if (dtset%pawxcdev==0) v_size=paw_an(iatom)%angl_size
     cplex=paw_ij(1)%cplex;idir=0;ipert=0
     paw_an(iatom)%nkxc1=0;if (rfphon/=0.or.rfelfd==1.or.rfelfd==3) paw_an(iatom)%nkxc1=nkxc
     paw_an(iatom)%has_vxc=1
     paw_an(iatom)%has_kxc=0;if (rfphon/=0.or.rfelfd==1.or.rfelfd==3) paw_an(iatom)%has_kxc=1
     paw_ij(iatom)%has_dijhartree=1
     ABI_ALLOCATE(paw_ij(iatom)%dijhartree,(pawtab(itypat)%lmn2_size))
     ABI_ALLOCATE(paw_an(iatom)%vxc1 ,(pawtab(itypat)%mesh_size,v_size,paw_an(iatom)%nspden))
     ABI_ALLOCATE(paw_an(iatom)%vxct1,(pawtab(itypat)%mesh_size,v_size,paw_an(iatom)%nspden))
     ABI_ALLOCATE(paw_an(iatom)%kxc1 ,(pawtab(itypat)%mesh_size,v_size,paw_an(iatom)%nkxc1))
     ABI_ALLOCATE(paw_an(iatom)%kxct1,(pawtab(itypat)%mesh_size,v_size,paw_an(iatom)%nkxc1))
     if (pawtab(itypat)%useexexch>0)  then
       ABI_ALLOCATE(paw_an(iatom)%vxc_ex,(pawtab(itypat)%mesh_size,v_size,paw_an(iatom)%nspden))
     end if
   end do
   call status(0,dtfil%filstat,iexit,level,'call pawdenpot')
   call pawdenpot(compch_sph,epaw,epawdc,ipert,dtset%ixc,mpi_enreg,natom,natom,dtset%nspden,ntypat,&
&   nzlmopt,option,dtset%paral_kgb,paw_an,paw_an,paw_ij,pawang,dtset%pawprtvol,&
&   pawrad,pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,&
&   dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,psps%znuclpsp)
   call status(0,dtfil%filstat,iexit,level,'call pawdij   ')
   call pawdij(cplex,dtset,dtset%enunit,one,gprimd,ipert,mpi_enreg,&
&   natom,natom,nfftf,ngfftf,dtset%nspden,ntypat,&
&   dtset%paral_kgb,paw_an,paw_ij,pawang,pawfgrtab,dtset%pawprtvol,pawrad,dtset%pawspnorb,pawtab,&
&   dtset%pawxcdev,k0,dtset%typat,ucvol,vtrial,vxc,xred)
   call symdij(gprimd,psps%indlmn,indsym,ipert,psps%lmnmax,natom,dtset%nsym,ntypat,0,paw_ij,pawang,&
&   dtset%pawprtvol,rprimd,dtset%symafm,symrec,dtset%typat)
   do iatom=1,natom
     ABI_DEALLOCATE(paw_ij(iatom)%dijhartree)
     ABI_DEALLOCATE(paw_an(iatom)%vxc1)
     ABI_DEALLOCATE(paw_an(iatom)%vxct1)
     paw_an(iatom)%has_vxc=0
     paw_ij(iatom)%has_dijhartree=0
     if (pawtab(itypat)%useexexch>0)  then
       ABI_DEALLOCATE(paw_an(iatom)%vxc_ex)
     end if
   end do
 end if

!-----2. Frozen-wavefunctions and Ewald(q=0) parts of 2DTE

 ABI_ALLOCATE(eltcore,(6,6))
 ABI_ALLOCATE(elteew,(6+3*natom,6))
 ABI_ALLOCATE(eltfrhar,(6,6))
 ABI_ALLOCATE(eltfrnl,(6+3*natom,6))
 ABI_ALLOCATE(eltfrloc,(6+3*natom,6))
 ABI_ALLOCATE(eltfrkin,(6,6))
 ABI_ALLOCATE(eltfrxc,(6+3*natom,6))
 eltcore(:,:)=zero
 elteew(:,:)=zero
 eltfrnl(:,:)=zero
 eltfrloc(:,:)=zero
 eltfrkin(:,:)=zero
 eltfrhar(:,:)=zero
 eltfrxc(:,:)=zero

 dyfr_nondiag=0;if (psps%usepaw==1.and.rfphon==1) dyfr_nondiag=1
 dyfr_cplex=1;if (psps%usepaw==1.and.rfphon==1.and.(.not.qeq0)) dyfr_cplex=2
 pawbec=0;if (psps%usepaw==1.and.rfphon==1.and.rfpert(natom+2)==1) pawbec=1
 ABI_ALLOCATE(dyew,(2,3,natom,3,natom))
 ABI_ALLOCATE(dyewq0,(3,3,natom))
 ABI_ALLOCATE(dyfrlo,(3,3,natom))
 ABI_ALLOCATE(dyfrx2,(3,3,natom))
 ABI_ALLOCATE(dyfrnl,(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag))
 ABI_ALLOCATE(dyfrwf,(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag))
 ABI_ALLOCATE(becfrnl,(3,natom,3*pawbec))
 dyew(:,:,:,:,:)=zero
 dyewq0(:,:,:)=zero
 dyfrnl(:,:,:,:,:)=zero
 dyfrwf(:,:,:,:,:)=zero
 dyfrlo(:,:,:)=zero
 dyfrx2(:,:,:)=zero
 if (pawbec==1) becfrnl(:,:,:)=zero

 if (rfphon==1) then

   call dyfnl3(atindx1,becfrnl,cg,cprj_dum,dimcprj,dyfrnl,dyfr_cplex,dyfr_nondiag,eigen0,gsqcut,indsym,&
&   dtset%istwfk,kg,dtset%kptns,dtset%kptopt,dtset%mband,dtset%mgfft,mgfftf,dtset%mkmem,mpi_enreg,mpsang,&
&   dtset%mpw,natom,nattyp,dtset%nband,nfftf,ngfft,ngfftf,dtset%nkpt,dtset%nloalg,npwarr,&
&   dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%nsym,ntypat,occ,dtset%paral_kgb,paw_ij,pawang,pawbec,&
&   dtset%pawprtvol,pawfgrtab,pawrad,pawrhoij,pawtab,ph1d,ph1df,psps,dtset%qptn(:),rprimd,dtset%symafm,symrec,&
&   dtset%typat,dtfil%unkg,dtfil%unpaw,dtfil%unylm,0,usexcnhat,wfftgs,vtrial,vxc,dtset%wtk,xred,ylm)

!  No more need of these local derivatives
   if (psps%usepaw==1) then
     do iatom=1,natom
       if (associated(pawfgrtab(iatom)%gylmgr2)) then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr2)
       end if
       pawfgrtab(iatom)%gylmgr2_allocated=0
     end do
   end if

!  dyfrnl has not yet been symmetrized, but will be in the next routine
   call status(0,dtfil%filstat,iexit,level,'call dyfro3    ')
   call dyfro3(atindx1,dyfrnl,dyfrlo,dyfrwf,dyfrx2,dyfr_cplex,dyfr_nondiag,&
&   gmet,gprimd,gsqcut,indsym,mgfftf,mpi_enreg,psps%mqgrid_vl,natom,nattyp,&
&   nfftf,ngfftf,dtset%nspden,dtset%nsym,ntypat,psps%n1xccc,n3xccc,dtset%paral_kgb,pawtab,ph1df,psps%qgrid_vl,&
&   dtset%qptn,rhog,rprimd,symq,symrec,dtset%typat,ucvol,psps%usepaw,psps%vlspl,&
&   vxc,psps%xcccrc,psps%xccc1d,xccc3d,xred)

!  The frozen-wavefunction part of the dynamical matrix is now:
!  dyfrnl:  non-local contribution
!  dyfrlo:  local contribution
!  dyfrx2:  2nd order xc core correction contribution
!  dyfrwf:  all      contributions
!  In case of PAW, it misses a term coming from the pertubed overlap operator

!  Compute Ewald (q=0) contribution
   qphon(:)=zero
   sumg0=0
   call status(0,dtfil%filstat,iexit,level,'call ewald3(1)')
   call ewald3(dyew,gmet,natom,qphon,rmet,sumg0,dtset%typat,ucvol,xred,psps%ziontypat)
   option=1
   call q0dy3_calc(natom,dyewq0,dyew,option)

!  End of the frozen-wavefunction and Ewald(q=0) parts of the dynamical matrix
 end if

!Section for the strain perturbation - frozen-wavefunction, Ewald, etc.
!parts of the elastic tensor

 if(rfstrs/=0) then

!  Verify that k-point set has full space-group symmetry; otherwise exit
   call status(0,dtfil%filstat,iexit,level,'call symkchk ')
   timrev=1
   call symkchk(dtset%kptns,dtset%nkpt,dtset%nsym,symrec,timrev)

!  Calculate the nonlocal part of the elastic tensor
   call status(0,dtfil%filstat,iexit,level,'call eltfrnl3 ')
   call eltfrnl3(atindx,atindx1,cg,eltfrnl,dtset%istwfk,&
&   kg,dtset%kptns,dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,mpsang,&
&   dtset%mpw,natom,nattyp,dtset%nband,dtset%nkpt,ngfft,dtset%nloalg,npwarr,&
&   dtset%nspinor,dtset%nsppol,ntypat,occ,ph1d,&
&   psps,rprimd,dtfil%unkg,wfftgs,dtfil%unylm,&
&   dtset%useylm,dtset%wtk,xred,ylm,ylmgr)
!  Calculate the kinetic part of the elastic tensor
   call status(0,dtfil%filstat,iexit,level,'call eltfrkin3')
   call eltfrkin3(cg,eltfrkin,dtset%ecut,dtset%ecutsm,dtset%effmass,&
&   dtset%istwfk,kg,dtset%kptns,dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,&
&   dtset%mpw,dtset%nband,dtset%nkpt,ngfft,npwarr,&
&   dtset%nspinor,dtset%nsppol,occ,rprimd,dtfil%unkg,wfftgs,dtset%wtk)

!  Calculate the hartree part of the elastic tensor
   call status(0,dtfil%filstat,iexit,level,'call eltfrhar3')
   call eltfrhar3(eltfrhar,rprimd,gsqcut,mpi_enreg,nfftf,ngfftf,rhog)

!  Calculate the xc part of the elastic tensor
   call status(0,dtfil%filstat,iexit,level,'call eltfrxc3 ')
   call eltfrxc3(eltfrxc,enxc,kxc,mpi_enreg,natom,&
&   nfftf,ngfftf,nkxc,dtset%nspden,ntypat,psps%n1xccc,n3xccc,dtset%paral_kgb,rhor,rprimd,&
&   dtset%typat,vxc,psps%xcccrc,psps%xccc1d,xccc3d,xred)

!  Calculate the local potential part of the elastic tensor
   call status(0,dtfil%filstat,iexit,level,'call eltfrloc3')
   call eltfrloc3(atindx,eltfrloc,gmet,gprimd,gsqcut,mgfftf,mpi_enreg,psps%mqgrid_vl,&
&   natom,nattyp,nfftf,ngfftf,ntypat,ph1df,psps%qgrid_vl,rhog,psps%vlspl)

!  Calculate the Ewald part of the elastic tensor
   call status(0,dtfil%filstat,iexit,level,'call ewald4')
   call ewald4(elteew,gmet,gprimd,natom,ntypat,rmet,rprimd,&
&   dtset%typat,ucvol,xred,psps%ziontypat)

!  Calculate the psp core energy part of elastic tensor (trivial)
   eltcore(1:3,1:3)=ecore/ucvol

 end if !rfstrs/=0
!End section for strain perturbation

 ABI_DEALLOCATE(vpsp)
 ABI_DEALLOCATE(xccc3d)

 if(dtset%prtvol==-level)then
   write(message,'(a,a)') ch10,&
&   ' respfn : frozen wavef. and Ewald(q=0) part of 2DTE done. '
   call wrtout(std_out,message,'COLL')
 end if

 call timab(136,2,tsec)

!-----3. Initialisation of 1st response, taking into account the q vector.

 call timab(137,1,tsec)

 write(message, '(a,a,a)' )ch10,&
& ' ==>  initialize data related to q vector <== ',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 qphon(:)=dtset%qptn(:)
 sumg0=1

!Treat the case of q=0 or q too close to 0
 qzero=0
 if(qeq0)then
   qphon(:)=zero
   write(message, '(a,a,a)' )&
&   ' respfn : the norm of the phonon wavelength (as input) was small (<1.d-7).',&
&   ch10,'  q has been set exactly to (0 0 0)'
   call wrtout(std_out,message,'COLL')
   sumg0=0
   qzero=1
 else
   if(rfelfd/=0 .or. rfstrs/=0 .or. rfddk /= 0)then
!    Temporarily, ...
     write(message, '(a,a,a,3es16.6,a,a,a,i2,a,i2,a,i2,a,i2,a,a,a)' )ch10,&
&     '  The treatment of non-zero wavevector q is restricted to phonons.',&
&     '  However, the input normalized qpt is',qphon(:),',',ch10,&
&     '  while rfelfd=',rfelfd,', rfddk=',rfddk,', and rfstrs=',rfstrs,'.',ch10,&
&     '  Action : change qpt, or rfelfd, or rfstrs in the input file.'
     MSG_ERROR(message)
   else if(rfasr.eq.2)then
     write(message, '(a,a)' )ch10,&
&     '  rfasr=2 not allowed with q/=0 => rfasr was reset to 0.'
     MSG_WARNING(message)
     rfasr=0
   end if
 end if

!Determine the symmetrical perturbations
 ABI_ALLOCATE(pertsy,(3,mpert))
 call status(0,dtfil%filstat,iexit,level,'call syper3   ')
 call syper3(indsym,mpert,natom,dtset%nsym,pertsy,rfdir,rfpert,symq,symrec,dtset%symrel)
 write(message, '(a)' ) &
& ' The list of irreducible perturbations for this q vector is:'
 call wrtout(ab_out,message,'COLL')
 ii=1
 do ipert=1,mpert
   do idir=1,3
     if(rfpert(ipert)==1.and.rfdir(idir)==1)then
       if( pertsy(idir,ipert)==1 )then
         write(message, '(i5,a,i2,a,i4)' )ii,')    idir=',idir,'    ipert=',ipert
         call wrtout(ab_out,message,'COLL')
         ii=ii+1
       end if
     end if
   end do
 end do

!test if the user left default rfdir 0 0 0
 if (ii==1) then
   write(message,'(5a)')ch10,' WARNING: no perturbations to be done at this q-point.', &
&   ch10, ' You may have forgotten to set the rfdir or rfatpol variables. Continuing normally.',ch10
   call wrtout(std_out,message,'COLL')
 end if

!Contribution to the dynamical matrix from ion-ion energy
 if(rfphon==1)then
   call status(0,dtfil%filstat,iexit,level,'call ewald3(2)')
   call ewald3(dyew,gmet,natom,qphon,rmet,sumg0,dtset%typat,ucvol,xred,psps%ziontypat)
   call q0dy3_apply(natom,dyewq0,dyew)
 end if

!1-order contribution of the xc core correction to the
!dynamical matrix
 ABI_ALLOCATE(dyfrx1,(2,3,natom,3,natom))
 dyfrx1(:,:,:,:,:)=zero
 if(rfphon==1.and.psps%n1xccc/=0)then
   ABI_ALLOCATE(blkflgfrx1,(3,natom,3,natom))
   call dyxc13(atindx,blkflgfrx1,dyfrx1,gmet,gsqcut,kxc,mgfftf,mpert,mpi_enreg,&
&   psps%mqgrid_vl,natom,nfftf,ngfftf,nkxc,dtset%nspden,&
&   ntypat,psps%n1xccc,dtset%paral_kgb,pawtab,ph1df,psps%qgrid_vl,qphon,&
&   rfdir,rfpert,rprimd,timrev,dtset%typat,ucvol,psps%usepaw,psps%xcccrc,psps%xccc1d,xred)
 end if

!Deallocate the arrays that were needed only for the frozen wavefunction part
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(ph1df)
 ABI_DEALLOCATE(cg)
 ABI_DEALLOCATE(kg)
 ABI_DEALLOCATE(npwarr)

!Close the unneeded temporary data files, if any
 if (dtset%mkmem==0) then
   close (unit=dtfil%unkg,status='delete')
   if (psps%useylm==1) close (unit=dtfil%unylm,status='delete')
   if (psps%usepaw==1) close (unit=dtfil%unpaw ,status='delete')
   call WffDelete(wfftgs,ierr)
 end if

 if(mpi_enreg%paral_compil==1) then
   ABI_DEALLOCATE(mpi_enreg%proc_distrb)
 end if

 ABI_ALLOCATE(blkflg,(3,mpert,3,mpert))
 ABI_ALLOCATE(d2eig0,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2k0,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2lo,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2loc0,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2nfr,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2nl,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2nl0,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2nl1,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2vn,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2ovl,(2,3,mpert,3,mpert*psps%usepaw))
 blkflg(:,:,:,:)=0
 d2eig0(:,:,:,:,:)=zero ; d2k0(:,:,:,:,:)=zero
 d2lo(:,:,:,:,:)=zero   ; d2loc0(:,:,:,:,:)=zero
 d2nfr(:,:,:,:,:)=zero  ; d2nl(:,:,:,:,:)=zero
 d2nl0(:,:,:,:,:)=zero  ; d2nl1(:,:,:,:,:)=zero
 d2vn(:,:,:,:,:)=zero
 if (psps%usepaw==1) d2ovl(:,:,:,:,:)=zero

 prtbbb=dtset%prtbbb
 ABI_ALLOCATE(d2bbb,(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb))
 ABI_ALLOCATE(d2cart_bbb,(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb))
 if(prtbbb==1)then
   d2cart_bbb(:,:,:,:,:,:)=zero
   d2bbb(:,:,:,:,:,:)=zero
 end if

 dim_eig2nkq = 0
 if(dtset%ieig2rf /= 0) dim_eig2nkq = 1
 ABI_ALLOCATE(eig2nkq,(2,dtset%mband*dtset%nsppol,dtset%nkpt,3,natom,3,natom*dim_eig2nkq))
 dim_eigbrd=0
 if(dtset%ieig2rf /= 0 .and. dtset%smdelta>0 ) dim_eigbrd = 1
 ABI_ALLOCATE(eigbrd,(2,dtset%mband*dtset%nsppol,dtset%nkpt,3,natom,3,natom*dim_eigbrd))

 call timab(137,2,tsec)


!Check whether exiting was required by the user.
!If found then do not start minimization steps
!At this first call to chkexi, initialize cpus
 cpus=dtset%cpus
 if(abs(cpus)>1.0d-5)cpus=cpus+cpui
 openexit=1 ; if(dtset%chkexit==0) openexit=0
 call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
 if (iexit==0) then

!  #######################################################################

   write(message,'(a,80a)')ch10,('=',mu=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   call status(0,dtfil%filstat,iexit,level,'call loper3   ')

   ddkfil(:)=0
!  Note that kg, cg, eigen0, mpw and npwarr are NOT passed to loper3 :
!  they will be reinitialized for each perturbation, with an eventually
!  reduced set of k point, thanks to the use of symmetry operations.
   call loper3(atindx,atindx1,blkflg,codvsn,cpus,dimcprj,dim_eigbrd,dim_eig2nkq,doccde,&
&   ddkfil,dtfil,dtset,dyew,dyfrlo,dyfrnl,dyfrx1,dyfrx2,dyfr_cplex,dyfr_nondiag,&
&   d2bbb,d2lo,d2nl,d2ovl,eigbrd,eig2nkq,&
&   eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,&
&   etotal,fermie,gsqcut_eff,iexit,indsym,kxc,&
&   dtset%mkmem,mkqmem,mk1mem,mpert,mpi_enreg,mpsang,nattyp,&
&   nfftf,dtset%nkpt,nkxc,dtset%nspden,dtset%nsym,occ,&
&   paw_an,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,&
&   pertsy,prtbbb,psps,rfpert,rhog,rhor,symq,symrec,timrev,&
&   usecprj,vtrial,vxc,vxcavg,xred)

!  #####################################################################

!  End of the check of hasty exit
 end if

 call timab(138,1,tsec)

 write(message, '(80a,a,a,a,a)' ) ('=',mu=1,80),ch10,ch10,&
& ' ---- first-order wavefunction calculations are completed ----',&
& ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 ABI_DEALLOCATE(vxc)
 ABI_DEALLOCATE(doccde)

!in paral_respfn-case the masters have to reconstruct some array's
 if(mpi_enreg%paral_compil_respfn==1) then
!  gather arrays from all cpu's
!  call xsum_master(blkflg,0,mpi_enreg%respfn_master_comm,ierr) ! Does not work on some machines
!  call xsum_master(d2lo,0,mpi_enreg%respfn_master_comm,ierr)
!  call xsum_master(d2nl,0,mpi_enreg%respfn_master_comm,ierr)
!  call xsum_master(vtrial,0,mpi_enreg%respfn_master_comm,ierr)
!  if (psps%usepaw==1) then
!  call xsum_master(d2ovl,0,mpi_enreg%respfn_master_comm,ierr)
!  endif
   call xsum_mpi(blkflg,mpi_enreg%respfn_master_comm,ierr)
   call xsum_mpi(d2lo,mpi_enreg%respfn_master_comm,ierr)
   call xsum_mpi(d2nl,mpi_enreg%respfn_master_comm,ierr)
   call xsum_mpi(vtrial,mpi_enreg%respfn_master_comm,ierr)
   if (psps%usepaw==1) then
     call xsum_mpi(d2ovl,mpi_enreg%respfn_master_comm,ierr)
   end if
 end if
!DEBUG
!write(std_out,*) "blkflg",mpi_enreg%paral_compil_respfn,":", blkflg
!write(std_out,*) "ddkfil",mpi_enreg%paral_compil_respfn,":", ddkfil
!write(std_out,*) "d2bbb-array",mpi_enreg%paral_compil_respfn,":", d2bbb
!write(std_out,*) "d2lo-array",mpi_enreg%paral_compil_respfn,":", d2lo
!write(std_out,*) "d2nl-array",mpi_enreg%paral_compil_respfn,":", d2nl
!if (psps%usepaw==1) write(std_out,*) "d2ovl-array",mpi_enreg%paral_compil_respfn,":", d2ovl
!write(std_out,*) "etotal",mpi_enreg%paral_compil_respfn,":", etotal
!write(std_out,*) "fermie",mpi_enreg%paral_compil_respfn,":", fermie
!write(std_out,*) "vtrial",mpi_enreg%paral_compil_respfn,":", vtrial
!write(std_out,*) "xred",mpi_enreg%paral_compil_respfn,":", xred
!ENDDEBUG

!Output of the localization tensor
 if ( rfpert(natom+1) /= 0 .and. (me == 0) .and. dtset%occopt<=2) then
   call wrtloctens(blkflg,d2bbb,d2nl,dtset%mband,mpert,natom,dtset%prtbbb,rprimd,psps%usepaw)
 end if

!The perturbation  natom+1 was only an auxiliary perturbation,
!needed to construct the electric field response, so its flag
!is now set to 0.
!rfpert(natom+1)=0

!Were 2DTE computed ?
 if(rfphon==0 .and. (rfddk/=0 .or. rfelfd==2) .and. rfstrs==0 .and. rfuser==0)then

   write(message,'(a,a)' )ch10,&
&   ' respfn : d/dk was computed, but no 2DTE, so no DDB output.'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  If 2DTE were computed, only one processor must output them and compute
!  frequencies.
 else if(me==0)then

   write(message,'(a,a)' )ch10,&
&   ' ==> Compute Derivative Database <== '
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  Open the formatted derivative database file, and write the
!  preliminary information
   call status(0,dtfil%filstat,iexit,level,'call ioddb8_ou')
   vrsddb=100401
   dscrpt=' Note : temporary (transfer) database '
!  tolwfr must be initialized here, but it is a dummy value
   tolwfr=1.0_dp
   call ioddb8_out (dscrpt,dtfil%fnameabo_ddb,natom,dtset%mband,&
&   dtset%nkpt,dtset%nsym,ntypat,dtfil%unddb,vrsddb,&
&   dtset%acell_orig(1:3,1),dtset%amu,dtset%dilatmx,dtset%ecut,dtset%ecutsm,&
&   dtset%intxc,iscf,dtset%ixc,dtset%kpt,dtset%kptnrm,&
&   natom,dtset%nband,ngfft,dtset%nkpt,dtset%nspden,dtset%nspinor,&
&   dtset%nsppol,dtset%nsym,ntypat,occ,dtset%occopt,dtset%pawecutdg,&
&   dtset%rprim_orig(1:3,1:3,1),dtset%sciss,dtset%spinat,dtset%symafm,dtset%symrel,&
&   dtset%tnons,tolwfr,dtset%tphysel,dtset%tsmear,&
&   dtset%typat,dtset%usepaw,dtset%wtk,xred,psps%ziontypat,dtset%znucl)

   nblok=1 ; fullinit=1 ; choice=2
   call psddb8 (choice,psps%dimekb,psps%ekb,fullinit,psps%indlmn,&
&   psps%lmnmax,nblok,ntypat,dtfil%unddb,pawtab,&
&   psps%pspso,psps%usepaw,psps%useylm,vrsddb)

!  In the RESPFN code, nstdy3 and stady3 were called here
   d2nfr(:,:,:,:,:)=d2lo(:,:,:,:,:)+d2nl(:,:,:,:,:)
   if (psps%usepaw==1) d2nfr(:,:,:,:,:)=d2nfr(:,:,:,:,:)+d2ovl(:,:,:,:,:)

!  In case of bbb decomposition
   if(prtbbb==1)then
     ABI_ALLOCATE(blkflg1,(3,mpert,3,mpert))
     ABI_ALLOCATE(blkflg2,(3,mpert,3,mpert))
     blkflg2(:,:,:,:) = blkflg(:,:,:,:)
     do ipert = 1, mpert
       do ipert2 = 1, mpert
         if ((ipert /= natom + 2).and.(ipert>natom).and.(ipert2/=natom+2)) then
           blkflg2(:,ipert2,:,ipert) = 0
         end if
       end do
     end do
     ABI_ALLOCATE(d2tmp,(2,3,mpert,3,mpert))
     do iband = 1,dtset%mband
       d2tmp(:,:,:,:,:)=zero
       blkflg1(:,:,:,:) = blkflg2(:,:,:,:)
       d2tmp(:,:,natom+2,:,:) = d2bbb(:,:,:,:,iband,iband)
       call d2sym3(blkflg1,d2tmp,indsym,mpert,natom,dtset%nsym,qphon,symq,&
&       symrec,dtset%symrel,timrev)
       d2bbb(:,:,:,:,iband,iband) = d2tmp(:,:,natom+2,:,:)
     end do
     ABI_DEALLOCATE(blkflg1)
     ABI_DEALLOCATE(blkflg2)
     ABI_DEALLOCATE(d2tmp)
   end if

!  Complete the d2nfr matrix by symmetrization of the existing elements
   call d2sym3(blkflg,d2nfr,indsym,mpert,natom,dtset%nsym,qphon,symq,symrec,dtset%symrel,timrev)

   if(rfphon==1.and.psps%n1xccc/=0)then
!    Complete the dyfrx1 matrix by symmetrization of the existing elements
     call d2sym3(blkflgfrx1,dyfrx1,indsym,natom,natom,dtset%nsym,qphon,symq,symrec,dtset%symrel,timrev)
   end if

!  Note that there is a bug in d2sym3 which will set some elements of
!  blkflg to 1 even when no corresponding symmetry-related element
!  has been computed.  This has the effect of producing spurious extra
!  output lines in the 2nd-order matrix listing in the .out file
!  and in the DDB file. The suprious matrix elements are all zero,
!  so this is primarily an annoyance.(DRH)


!  Add the frozen-wf (dyfrwf) part to the ewald part (dyew),
!  the part 1 of the frozen wf part of the xc core correction
!  (dyfrx1) and the non-frozen part (dynfr) to get the second-order
!  derivative matrix (d2matr), then
!  take account of the non-cartesian coordinates (d2cart).
   ABI_ALLOCATE(d2cart,(2,3,mpert,3,mpert))
   ABI_ALLOCATE(carflg,(3,mpert,3,mpert))
   ABI_ALLOCATE(d2matr,(2,3,mpert,3,mpert))
   outd2=1

   call status(0,dtfil%filstat,iexit,level,'call gath3    ')
   call gath3(becfrnl,dtset%berryopt,blkflg,carflg,&
&   dyew,dyfrwf,dyfrx1,dyfr_cplex,dyfr_nondiag,d2bbb,d2cart,d2cart_bbb,d2matr,d2nfr,&
&   eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,&
&   gprimd,dtset%mband,mpert,natom,ntypat,outd2,pawbec,dtset%prtbbb,&
&   rfasr,rfpert,rprimd,dtset%typat,ucvol,psps%ziontypat)

!  Output of the dynamical matrix
!  (Note : remember, previously, the processor me=0 has been selected)
   call status(0,dtfil%filstat,iexit,level,'call dyout3   ')
   call dyout3(becfrnl,dtset%berryopt,blkflg,carflg,dtfil%unddb,ddkfil,dyew,dyfrlo,&
&   dyfrnl,dyfrx1,dyfrx2,dyfr_cplex,dyfr_nondiag,d2cart,d2cart_bbb,d2eig0,&
&   d2k0,d2lo,d2loc0,d2matr,d2nl,d2nl0,d2nl1,d2ovl,d2vn,&
&   eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,&
&   ab_out,dtset%mband,mpert,natom,ntypat,&
&   outd2,pawbec,dtset%prtbbb,dtset%prtvol,qphon,qzero,dtset%typat,rfdir,rfpert,rfphon,&
&   rfstrs,psps%usepaw,psps%ziontypat)

   close(dtfil%unddb)

!  In case of phonons, diagonalize the dynamical matrix
   if(rfphon==1)then

!    First, suppress the 'wings' elements,
!    for which the diagonal element is not known
     call wings3(carflg,d2cart,mpert)

!    Check the analyticity of the dynamical matrix
     analyt=0
     if (rfpert(natom+2)==0 .or. rfpert(natom+2)==2 .or. sumg0==1 ) analyt=1

!    Diagonalize the analytic part
     ABI_ALLOCATE(displ,(2*3*natom*3*natom))
     ABI_ALLOCATE(eigval,(3*natom))
     ABI_ALLOCATE(eigvec,(2*3*natom*3*natom))
     ABI_ALLOCATE(phfrq,(3*natom))
     qphnrm=one
     call phfrq3(dtset%amu,displ,d2cart,eigval,eigvec,indsym,mpert,&
&     dtset%nsym,natom,dtset%nsym,ntypat,phfrq,qphnrm,qphon,&
&     dtset%rprimd_orig(1:3,1:3,1),0,dtset%symrel,dtset%typat,ucvol)

!    Print the phonon frequencies
     call prtph3(displ,0,dtset%enunit,-1,ab_out,natom,phfrq,qphnrm,qphon)

!    Check the completeness of the dynamical matrix and eventually send a warning
     call chkph3(carflg,0,mpert,natom)

!    Compute and print the T=0 Fan, and possibly DDW contributions to the eigenenergies.
     if(dtset%ieig2rf > 0) then
       write(message, '(80a,9a)' ) ('=',mu=1,80),ch10,ch10,&
&       ' ---- T=0 shift of eigenenergies due to electron-phonon interation at q ---- ',ch10,&
&       ' Warning : the total shift must be computed through anaddb,                  ',ch10,&
&       ' here, only the contribution of one q point is printed.                      ',ch10,&
&       ' Print first the electronic eigenvalues, then the q-dependent Fan shift of eigenvalues.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')

       if(qeq0)then
         write(message, '(a)' )&
&         ' Phonons at gamma, also compute the Diagonal Debye-Waller shift of eigenvalues.'
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,  message,'COLL')
       end if

       write(message, '(a)' ) ' '
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')

       call prteigrs(eigen0,dtset%enunit,fermie,dtfil%fnameabo_eig,ab_out,0,dtset%kptns,dtset%kptopt,&
&       dtset%mband,dtset%nband,dtset%nkpt,1,dtset%nsppol,occ,dtset%occopt,3,0,dtset%prtvol,&
&       eigen0,zero,zero,dtset%wtk)

       write(message, '(a)' ) ch10
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')

!      Compute and print Fan contribution
       ABI_ALLOCATE(eigen_fan,(dtset%mband*dtset%nkpt*dtset%nsppol))
       ABI_ALLOCATE(eigen_fan_mean,(dtset%mband*dtset%nkpt*dtset%nsppol))
       call elph2_fanddw(dim_eig2nkq,displ,eig2nkq,eigen_fan,gprimd,&
&       dtset%mband,natom,dtset%nkpt,dtset%nsppol,1,phfrq)
       call eigen_meandege(eigen0,eigen_fan,eigen_fan_mean,dtset%mband,dtset%nband,dtset%nkpt,dtset%nsppol,2)
       call prteigrs(eigen_fan_mean,dtset%enunit,fermie,dtfil%fnameabo_eig,ab_out,0,dtset%kptns,dtset%kptopt,&
&       dtset%mband,dtset%nband,dtset%nkpt,1,dtset%nsppol,occ,dtset%occopt,5,0,dtset%prtvol,&
&       eigen0,zero,zero,dtset%wtk)

       if(qeq0 .or. dtset%getgam_eig2nkq>0)then

         write(message, '(a)' ) ch10
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,  message,'COLL')

!        Compute and print Diagonal Debye-Waller contribution
         ABI_ALLOCATE(eigen_ddw,(dtset%mband*dtset%nkpt*dtset%nsppol))
         ABI_ALLOCATE(eigen_ddw_mean,(dtset%mband*dtset%nkpt*dtset%nsppol))
         if(qeq0)then
           call elph2_fanddw(dim_eig2nkq,displ,eig2nkq,eigen_ddw,gprimd,&
&           dtset%mband,natom,dtset%nkpt,dtset%nsppol,2,phfrq)
           if(results_respfn%gam_jdtset == -dtset%jdtset)then
             ABI_ALLOCATE(results_respfn%gam_eig2nkq,(2,dtset%mband*dtset%nsppol,dtset%nkpt,3,natom,3,natom*dim_eig2nkq))
             results_respfn%gam_eig2nkq(:,:,:,:,:,:,:)=eig2nkq(:,:,:,:,:,:,:)
             results_respfn%gam_jdtset=dtset%jdtset
           end if
         else if(dtset%getgam_eig2nkq>0)then
           if(results_respfn%gam_jdtset==dtset%getgam_eig2nkq)then
             call elph2_fanddw(dim_eig2nkq,displ,results_respfn%gam_eig2nkq,eigen_ddw,&
&             gprimd,dtset%mband,natom,dtset%nkpt,dtset%nsppol,2,phfrq)
           else
             write(message, '(4a,i5,2a,i5,2a)' )ch10,&
&             ' respfn : BUG -',ch10,&
&             '  results_respfn%gam_jdtset=',results_respfn%gam_jdtset,ch10,&
&             '  dtset%getgam_eig2nkq=',dtset%getgam_eig2nkq,ch10,&
&             '  So, it seems eig2nkq at gamma has not yet been computed, while it is needed now.'
             call wrtout(std_out,message,'COLL')
             call leave_new('COLL')
           end if
         end if
         call eigen_meandege(eigen0,eigen_ddw,eigen_ddw_mean,dtset%mband,dtset%nband,dtset%nkpt,dtset%nsppol,2)
         call prteigrs(eigen_ddw_mean,dtset%enunit,fermie,dtfil%fnameabo_eig,ab_out,0,dtset%kptns,dtset%kptopt,&
&         dtset%mband,dtset%nband,dtset%nkpt,1,dtset%nsppol,occ,dtset%occopt,6,0,dtset%prtvol,&
&         eigen0,zero,zero,dtset%wtk)

         write(message, '(a)' ) ch10
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,  message,'COLL')

!        Print sum of mean Fan and DDW
         ABI_ALLOCATE(eigen_fanddw,(dtset%mband*dtset%nkpt*dtset%nsppol))
         eigen_fanddw=eigen_fan_mean+eigen_ddw_mean
         call prteigrs(eigen_fanddw,dtset%enunit,fermie,dtfil%fnameabo_eig,ab_out,0,dtset%kptns,dtset%kptopt,&
&         dtset%mband,dtset%nband,dtset%nkpt,1,dtset%nsppol,occ,dtset%occopt,7,0,dtset%prtvol,&
&         eigen0,zero,zero,dtset%wtk)

         ABI_DEALLOCATE(eigen_ddw)
         ABI_DEALLOCATE(eigen_ddw_mean)
         ABI_DEALLOCATE(eigen_fanddw)

       end if

       ABI_DEALLOCATE(eigen_fan)
       ABI_DEALLOCATE(eigen_fan_mean)
     end if

!    for SCphon print out necessary files: _PHFRQ _PHVEC _PCINFO
     if (dtset%prepscphon == 1) then
!      close unit for phonon frequencies
       phonon_freq_filename = trim(dtfil%filnam_ds(4))//"_PHFRQ"
       phonon_vec_filename = trim(dtfil%filnam_ds(4))//"_PHVEC"
       call beginprtscphon(phonon_freq_filename, phonon_vec_filename,1)
       call prtscphon(eigvec, natom, phfrq, qphon)
       call endprtscphon()
     end if

!    In case of a non-analytical part,
!    get the phonon frequencies for three different directions
!    (in cartesian coordinates)
     if(analyt==0)then
       qphnrm=zero
       do idir=1,3
!        Need to know the corresponding dielectric constant
         if(carflg(idir,natom+2,idir,natom+2)==1)then
           qphon(:)=zero ; qphon(idir)=one
!          Get the phonon frequencies
           call phfrq3(dtset%amu,displ,d2cart,eigval,eigvec,indsym,mpert,&
&           dtset%nsym,natom,dtset%nsym,ntypat,phfrq,qphnrm,qphon,&
&           dtset%rprimd_orig(1:3,1:3,1),0,dtset%symrel,dtset%typat,ucvol)
!          Print the phonon frequencies
           call prtph3(displ,0,dtset%enunit,-1,ab_out,natom,phfrq,qphnrm,qphon)
!          Check the completeness of the dynamical matrix
!          and eventually send a warning
           call chkph3(carflg,idir,mpert,natom)
         end if
       end do
       if (idir < 4) then
         qphon(idir)=zero
       end if
     end if

     ABI_DEALLOCATE(displ)
     ABI_DEALLOCATE(eigval)
     ABI_DEALLOCATE(eigvec)
     ABI_DEALLOCATE(phfrq)

!    End condition on phonons
   end if

   ABI_DEALLOCATE(carflg)
   ABI_DEALLOCATE(d2cart)
   ABI_DEALLOCATE(d2matr)

!  End of the DDB output part
 end if

!Deallocate arrays
 ABI_DEALLOCATE(amass)
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(atindx1)
 ABI_DEALLOCATE(blkflg)
 ABI_DEALLOCATE(dyew)
 ABI_DEALLOCATE(dyewq0)
 ABI_DEALLOCATE(dyfrlo)
 ABI_DEALLOCATE(dyfrnl)
 ABI_DEALLOCATE(dyfrwf)
 ABI_DEALLOCATE(dyfrx1)
 ABI_DEALLOCATE(dyfrx2)
 ABI_DEALLOCATE(d2bbb)
 ABI_DEALLOCATE(d2cart_bbb)
 ABI_DEALLOCATE(d2eig0)
 ABI_DEALLOCATE(d2k0)
 ABI_DEALLOCATE(d2lo)
 ABI_DEALLOCATE(d2loc0)
 ABI_DEALLOCATE(d2nfr)
 ABI_DEALLOCATE(d2nl)
 ABI_DEALLOCATE(d2nl0)
 ABI_DEALLOCATE(d2nl1)
 ABI_DEALLOCATE(d2ovl)
 ABI_DEALLOCATE(d2vn)
 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(eig2nkq)
 ABI_DEALLOCATE(eigbrd)
 ABI_DEALLOCATE(eltcore)
 ABI_DEALLOCATE(elteew)
 ABI_DEALLOCATE(eltfrhar)
 ABI_DEALLOCATE(eltfrnl)
 ABI_DEALLOCATE(eltfrloc)
 ABI_DEALLOCATE(eltfrkin)
 ABI_DEALLOCATE(eltfrxc)
 ABI_DEALLOCATE(grxc)
 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(kxc)
 ABI_DEALLOCATE(nattyp)
 ABI_DEALLOCATE(pertsy)
 ABI_DEALLOCATE(rfpert)
 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(rhor)
 ABI_DEALLOCATE(symq)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(vtrial)
 ABI_DEALLOCATE(ylm)
 ABI_DEALLOCATE(ylmgr)
 ABI_DEALLOCATE(pawfgr%fintocoa)
 ABI_DEALLOCATE(pawfgr%coatofin)
 if (pawbec==1)  then
   ABI_DEALLOCATE(becfrnl)
 end if
 if (psps%usepaw==1) then
   if (nhatdim>0)  then
     ABI_DEALLOCATE(nhat)
   end if
   if (nhatgrdim>0)  then
     ABI_DEALLOCATE(nhatgr)
   end if
   call rhoij_free(pawrhoij)
   ABI_DEALLOCATE(pawrhoij)
   call pawfgrtab_free(pawfgrtab)
   call destroy_paw_an(paw_an)
   call destroy_paw_ij(paw_ij)
   ABI_DEALLOCATE(pawfgrtab)
   ABI_DEALLOCATE(paw_an)
   ABI_DEALLOCATE(paw_ij)
   ABI_DEALLOCATE(dimcprj)
 end if
 if(rfphon==1.and.psps%n1xccc/=0)then
   ABI_DEALLOCATE(blkflgfrx1)
 end if

!Clean the header
 call hdr_clean(hdr)

!Clean GPU data
#if defined HAVE_GPU_CUDA
 if (dtset%use_gpu_cuda==1) then
   call dealloc_hamilt_gpu(0,dtset%use_gpu_cuda)
 end if
#endif

!Clean MPI datas
 call clnmpi_fft(mpi_enreg)
 call clnmpi_atom(mpi_enreg)

 write(message, '(a,a)' ) ch10,' respfn : exiting '
 call wrtout(std_out,message,'COLL')

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(138,2,tsec)
 call timab(132,2,tsec)

 DBG_EXIT("COLL")

end subroutine respfn
!!***
