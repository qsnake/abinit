!{\src2tex{textfont=tt}}
!!****f* ABINIT/energy
!! NAME
!! energy
!!
!! FUNCTION
!! Compute electronic energy terms
!! energies%e_eigenvalues, ek and enl from arbitrary (orthonormal) provided wf,
!! ehart, enxc, and eei from provided density and potential,
!! energies%e_eigenvalues=Sum of the eigenvalues - Band energy (Hartree)
!! ek=kinetic energy, ehart=Hartree electron-electron energy,
!! enxc,enxcdc=exchange-correlation energies, eei=local pseudopotential energy,
!! enl=nonlocal pseudopotential energy
!! Also, compute new density from provided wfs, after the evaluation
!! of ehart, enxc, and eei.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, AR, MB, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>=Fourier coefficients of wavefunction
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem=maximum number of k points which can fit in core memory
!!   | mpw=maximum dimension for number of planewaves
!!   | natom=number of atoms in unit cell
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for polarized
!!   | nspinor=number of spinorial components
!!   | nsym=number of symmetry elements in space group (at least 1)
!!   | occopt=option for occupancies
!!   | tsmear=smearing energy or temperature (if metal)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gsqcut=G^2 cutoff from gsqcut=ecut/(2 Pi^2)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)=array describing how many atoms of each type in cell
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  ngfftf(18)= -PAW ONLY- contain all needed information about 3D FFT for the fine grid
!!              (ngfftf=ngfft for norm-conserving potential runs)
!!  nhatgr(nfft,nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  n3xccc=dimension of the xccc3d array (0 or nfftf).
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2) at each k point
!!  optene=option for the computation of total energy (direct scheme or double-counting scheme)
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr(natom) <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase information related to structure factor
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!   | ntypat=number of types of atoms in cell
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vpsp(nfftf)=local pseudopotential in real space (hartree)
!!  wffnow=structured array giving all information about wavefunction file
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  xred(3,natom)=reduced coordinates of atoms (dimensionless)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!
!! OUTPUT
!!  compch_fft=-PAW only- compensation charge inside spheres computed over fine fft grid
!!  etotal=total energy (hartree):
!!    - computed by direct scheme if optene=0 or 2
!!    - computed by double-counting scheme if optene=1 or 3
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points (hartree^2)
!!  strsxc(6)=exchange-correlation contribution to stress tensor
!!  vhartr(nfftf)=work space to hold Hartree potential in real space (hartree)
!!  vtrial(nfftf,nspden)=total local potential (hartree)
!!  vxc(nfftf,nspden)=work space to hold Vxc(r) in real space (hartree)
!!  vxctau(nfftf,dtset%nspden*dtset%usekden,4)=derivative of XC energy density with respect to
!!    kinetic energy density (metaGGA cases) (optional output)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  energies <type(energies_type)>=all part of total energy.
!!   | entropy(IN)=entropy due to the occupation number smearing (if metal)
!!   | e_ewald(IN)=Ewald energy (hartree)
!!   | e_corepsp(IN)=psp core-core energy
!!   | e_paw(IN)=PAW spherical part energy
!!   | e_pawdc(IN)=PAW spherical part double-counting energy
!!   | e_eigenvalues(OUT)=Sum of the eigenvalues - Band energy (Hartree)
!!   | e_hartree(OUT)=Hartree part of total energy (hartree units)
!!   | e_kinetic(OUT)=kinetic energy part of total energy.
!!   | e_nonlocalpsp(OUT)=nonlocal pseudopotential part of total energy.
!!   | e_xc(OUT)=exchange-correlation energy (hartree)
!!  ==== if optene==0, 2 or 3
!!   | e_localpsp(OUT)=local psp energy (hartree)
!!  ==== if optene==1, 2 or 3
!!   | e_xcdc(OUT)=exchange-correlation double-counting energy (hartree)
!!  rhog(2,nfftf)=work space for rho(G); save intact on return (? MT 08-12-2008: is that true now ?)
!!  rhor(nfftf,nspden)=work space for rho(r); save intact on return (? MT 08-12-2008: is that true now ?)
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  nspinor should not be modified in the call of rdnpw
!!  === if psps%usepaw==1 ===
!!    nhat(nfftf,nspden*usepaw)= compensation charge density
!!    pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!
!! NOTES
!!  Be careful to the meaning of nfft (size of FFT grids):
!!   - In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!   - In case of PAW calculations:
!!     Two FFT grids are used; one with nfft points (coarse grid) for
!!     the computation of wave functions ; one with nfftf points
!!     (fine grid) for the computation of total density.
!!
!!  There is a large amount of overhead in the way this routine do the computation of the energy !
!!  For example, the density has already been precomputed, so why to compute it again here ??
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      cprj_alloc,cprj_free,cprj_gather_spin,dotprod_vn,fftpac,fourdp,hdr_skip
!!      leave_new,leave_test,meanvalue_g,metric,mkffnl,mkkin,mkresi,mkrho
!!      nonlop,pawaccrhoij,pawmknhat,ph1d3d,psolver_rhohxc,rdnpw,rhohxc
!!      rhohxcpositron,rwwf,sphereboundary,symrhoij,timab,transgrid,wrtout
!!      xcomm_init,xdefineoff,xmaster_init,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine energy(atindx,atindx1,cg,compch_fft,dtfil,dtset,electronpositron,&
& energies,eigen,etotal,gsqcut,indsym,irrzon,kg,mcg,mpi_enreg,nattyp,nfftf,ngfftf,nhat,&
& nhatgr,nhatgrdim,npwarr,n3xccc,occ,optene,paw_ij,pawang,pawfgr,&
& pawfgrtab,pawrhoij,pawtab,phnons,ph1d,psps,resid,rhog,rhor,rprimd,strsxc,symrec,&
& taug,taur,usexcnhat,vhartr,vtrial,vpsp,vxc,wffnow,wfs,wvl,xccc3d,xred,ylm,&
& vxctau) ! optional arguments

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_energies, only : energies_type
 use m_xmpi
 use m_wffile

 use m_electronpositron, only : electronpositron_type,electronpositron_calctype
 use m_paw_dmft,         only: paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'energy'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_53_spacepar
 use interfaces_56_recipspace
 use interfaces_56_xc
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_62_poisson
 use interfaces_65_nonlocal
 use interfaces_66_paw
 use interfaces_67_common, except_this_one => energy
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,n3xccc,nfftf,nhatgrdim,optene,usexcnhat
 real(dp),intent(in) :: gsqcut
 real(dp),intent(out) :: compch_fft,etotal
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow
 type(wvl_wf_type),intent(inout) :: wfs
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
!no_abirules
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)
 integer :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)),kg(3,dtset%mpw*dtset%mkmem)
 integer, intent(in) :: nattyp(psps%ntypat),ngfftf(18),npwarr(dtset%nkpt),symrec(3,3,dtset%nsym)
 real(dp), intent(in) :: cg(2,mcg),eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol),ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
 real(dp), intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
 real(dp),intent(in) :: nhatgr(nfftf,dtset%nspden,3*nhatgrdim)
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden)
 real(dp), intent(inout) :: taug(2,nfftf*dtset%usekden),taur(nfftf,dtset%nspden*dtset%usekden)
 real(dp), intent(out) :: strsxc(6)
 real(dp), intent(in) :: rprimd(3,3),vpsp(nfftf),xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp), intent(out) :: vhartr(nfftf),vtrial(nfftf,dtset%nspden),vxc(nfftf,dtset%nspden)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(out),optional :: vxctau(nfftf,dtset%nspden*dtset%usekden,4)
 type(paw_ij_type), intent(in) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type), intent(in)  :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,bufdim,choice,cplex,cplex_rf,cpopt,dimdij,dimffnl,formeig,ia
 integer :: iatom,iband,icg,ider,idir,ierr,ifft,ig,ii,ikg,ikpt,ilm,ilmn
 integer :: indx,ipert,ipositron,iresid,isp,ispden,isppol,istwf_k,itypat,izero
 integer :: master,matblk,mcg_disk,me_distrb,muig,my_nspinor,n1,n2,n3,n4,n5,n6
 integer :: nband_k,nfftotf,nkpg,nkxc,nk3xc,nnlout,npw_k,nsp,nsp2,nvloc,option
 integer :: option_rhoij,paw_opt,signs,spaceComm,tim_mkrho,tim_nonlop,tim_rwwf
 logical :: usetimerev, with_vxctau
 real(dp) :: arg,doti,dum,e_xcdc_vxctau,eeigk,ekk,enlk,ucvol,vxcavg
 character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamk
 type(paw_dmft_type) :: paw_dmft
!arrays
 integer,allocatable :: dimlmn(:),kg_dum(:,:),kg_k(:,:)
 real(dp) :: enlout(1),gmet(3,3),gprimd(3,3),kpoint(3),nonlop_dum(1,1)
 real(dp) :: qpt(3),rhodum(1),rmet(3,3),tsec(2),ylmgr_dum(1)
 real(dp) :: vzeeman(4)
 real(dp),allocatable :: buffer(:),buffer2(:),cg_disk(:,:),cgrvtrial(:,:)
 real(dp),allocatable :: cwavef(:,:),eig_dum(:),eig_k(:),ffnl(:,:,:,:),kinpw(:)
 real(dp),allocatable :: kpg_dum(:,:),kxc(:,:),occ_dum(:),occ_k(:)
 real(dp),allocatable :: ph3d(:,:,:)
 real(dp),allocatable :: resid_k(:),rhowfg(:,:),rhowfr(:,:),vlocal(:,:,:,:)
 real(dp),allocatable :: vlocal_tmp(:,:,:),ylm_k(:,:)
 type(cprj_type),target,allocatable :: cwaveprj(:,:)
 type(cprj_type),pointer :: cwaveprj_gat(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' energy : enter '
!stop
!ENDDEBUG

!Check that usekden is not 0 if want to use vxctau
 with_vxctau = .false.
 if (present(vxctau) .and. dtset%usekden /= 0) with_vxctau = .true.

!Test size of FFT grids (1 grid in norm-conserving, 2 grids in PAW)
 nfftotf=ngfftf(1)*ngfftf(2)*ngfftf(3)
 if ((psps%usepaw==1.and.pawfgr%nfft/=nfftf).or.(psps%usepaw==0.and.dtset%nfft/=nfftf)) then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' energy :  BUG -',ch10,&
&   '  wrong values for nfft, nfftf !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 call timab(59,1,tsec)

!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_kpt)
!Init me
 call xme_init(mpi_enreg,me_distrb)
!Init master
 call xmaster_init(mpi_enreg,master)

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Compute Hxc potential from density
 option=1;nkxc=0
 ipositron=electronpositron_calctype(electronpositron)

 if (ipositron/=1) then

   if (dtset%icoulomb == 0) then
!    Use the periodic solver to compute Hxc.
     ABI_ALLOCATE(kxc,(1,nkxc))
!    to be adjusted for the call to rhohxc
     nk3xc=1
     if (ipositron==0) then
       if(with_vxctau)then
         call rhohxc(dtset,energies%e_xc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,&
&         ngfftf,nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,option,rhog,rhor,rprimd,strsxc,&
&         usexcnhat,vhartr,vxc,vxcavg,xccc3d,taug=taug,taur=taur,vxctau=vxctau)
       else
         call rhohxc(dtset,energies%e_xc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,&
&         ngfftf,nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,option,rhog,rhor,rprimd,strsxc,&
&         usexcnhat,vhartr,vxc,vxcavg,xccc3d,taug=taug,taur=taur)
       end if
     else
       if(with_vxctau)then
         call rhohxc(dtset,energies%e_xc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,&
&         ngfftf,nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,option,rhog,rhor,rprimd,strsxc,&
&         usexcnhat,vhartr,vxc,vxcavg,xccc3d,electronpositron=electronpositron,taug=taug,taur=taur,vxctau=vxctau)
       else
         call rhohxc(dtset,energies%e_xc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,&
&         ngfftf,nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,option,rhog,rhor,rprimd,strsxc,&
&         usexcnhat,vhartr,vxc,vxcavg,xccc3d,electronpositron=electronpositron,taug=taug,taur=taur)
       end if
     end if
     ABI_DEALLOCATE(kxc)
   else
!    Use the free boundary solver.
     call PSolver_rhohxc(dtset, energies%e_hartree, energies%e_xc, energies%e_vxc, &
&     mpi_enreg, rhor, rprimd, vhartr, vxc, vxcavg, wvl)
   end if
 else
   energies%e_xc=zero
   call rhohxcpositron(electronpositron,gprimd,kxc,mpi_enreg,nfftf,ngfftf,nhat,nkxc,dtset%nspden,n3xccc,&
&   dtset%paral_kgb,rhor,strsxc,ucvol,usexcnhat,psps%usepaw,vhartr,vxc,vxcavg,xccc3d,dtset%xc_denpos)
 end if
 if (ipositron/=0) then
   call dotprod_vn(1,rhor,electronpositron%e_hartree,doti,mpi_enreg,&
&   nfftf,nfftotf,1,1,electronpositron%vha_ep,ucvol)
   vhartr=vhartr+electronpositron%vha_ep
 end if

!Total local potential (for either spin channel) is
!Hartree + local psp + Vxc(spin), minus its mean
!(Note : this potential should agree with the input vtrial)

 vzeeman(:) = zero
 if(dtset%nspden==2)vzeeman(2) = dtset%zeemanfield(3) ! For collinear ispden=2 is rho_up only
 if(dtset%nspden==4)then
   do ispden=2,4
     vzeeman(ispden)=dtset%zeemanfield(ispden-1)
   end do !ispden
 end if


 do ispden=1,min(dtset%nspden,2)
   do ifft=1,nfftf
     vtrial(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)+vzeeman(ispden)
   end do
 end do
 if (dtset%nspden==4) then
   do ifft=1,nfftf
     vtrial(ifft,3:4)=vxc(ifft,3:4)+vzeeman(3:4)
   end do
 end if

!Compute Hartree energy - use up+down rhor
 if (ipositron/=1) then
   call dotprod_vn(1,rhor,energies%e_hartree ,doti,mpi_enreg,nfftf,nfftotf,1,1,vhartr,ucvol)
   if (ipositron==0) energies%e_hartree=half*energies%e_hartree
   if (ipositron==2) energies%e_hartree = half *(energies%e_hartree-electronpositron%e_hartree)
 else
   energies%e_hartree=zero
 end if

!Compute local psp energy - use up+down rhor
 if (optene/=1) then
   call dotprod_vn(1,rhor,energies%e_localpsp,doti,mpi_enreg,nfftf,nfftotf,1,1,vpsp,ucvol)
 end if

!Compute DC-xc energy - use up+down rhor
 if (optene>0) then
   if (ipositron/=1) then
     if (psps%usepaw==0.or.usexcnhat/=0) then
       call dotprod_vn(1,rhor,energies%e_xcdc,doti,mpi_enreg,nfftf,nfftotf,dtset%nspden,1,vxc,ucvol)
       if(with_vxctau)then
         call dotprod_vn(1,taur,e_xcdc_vxctau,doti,mpi_enreg,nfftf,nfftotf,dtset%nspden,1,vxctau(:,:,1),ucvol)
         energies%e_xcdc=energies%e_xcdc+e_xcdc_vxctau
       end if
     else
       ABI_ALLOCATE(rhowfr,(nfftf,dtset%nspden))
       rhowfr=rhor-nhat
       call dotprod_vn(1,rhowfr,energies%e_xcdc,doti,mpi_enreg,nfftf,nfftotf,dtset%nspden,1,vxc,ucvol)
       ABI_DEALLOCATE(rhowfr)
     end if
     if (ipositron==2) energies%e_xcdc=energies%e_xcdc-electronpositron%e_xcdc
   else
     energies%e_xcdc=zero
   end if
 end if

 if (dtset%mkmem==0) then

!  Read wavefunction file header
   call hdr_skip(wffnow,ierr)

!  Define offsets, in case of MPI I/O
   formeig=0
   call xdefineOff(formeig,wffnow,mpi_enreg,dtset%nband,npwarr,dtset%nspinor,dtset%nsppol,dtset%nkpt)

   mcg_disk=dtset%mpw*my_nspinor*dtset%mband
   ABI_ALLOCATE(cg_disk,(2,mcg_disk))

 end if

 energies%e_eigenvalues=zero
 energies%e_kinetic=zero
 energies%e_nonlocalpsp=zero
 bdtot_index=0
 icg=0


!DEBUG
!write(std_out,*)' energy : before loop over spins '
!stop
!ENDDEBUG

 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)
 nvloc=1
 if(dtset%nspden==4)nvloc=4
 ABI_ALLOCATE(vlocal,(n4,n5,n6,nvloc))
 ABI_ALLOCATE(kg_k,(3,dtset%mpw))
 ABI_ALLOCATE(cwavef,(2,dtset%mpw*my_nspinor))

!Allocate the arrays of the Hamiltonian whose dimensions do not depend on k
 ABI_ALLOCATE(gs_hamk%atindx,(dtset%natom))
 ABI_ALLOCATE(gs_hamk%atindx1,(dtset%natom))
 ABI_ALLOCATE(gs_hamk%gbound,(2*dtset%mgfft+8,2))
 ABI_ALLOCATE(gs_hamk%indlmn,(6,psps%lmnmax,psps%ntypat))
 ABI_ALLOCATE(gs_hamk%nattyp,(psps%ntypat))
 ABI_ALLOCATE(gs_hamk%phkxred,(2,dtset%natom))
 ABI_ALLOCATE(gs_hamk%ph1d,(2,3*(2*dtset%mgfft+1)*dtset%natom))
 ABI_ALLOCATE(gs_hamk%pspso,(psps%ntypat))
 ABI_ALLOCATE(gs_hamk%xred,(3,dtset%natom))

!Initialize most of the Hamiltonian
 gs_hamk%atindx(:)  =atindx(:)
 gs_hamk%atindx1(:) =atindx1(:)
 gs_hamk%gmet(:,:)  =gmet(:,:)
 gs_hamk%gprimd(:,:)=gprimd(:,:)
 gs_hamk%indlmn(:,:,:)=psps%indlmn(:,:,:)
 gs_hamk%lmnmax     =psps%lmnmax
 gs_hamk%mgfft      =dtset%mgfft
 gs_hamk%mpsang     =psps%mpsang
 gs_hamk%mpssoang   =psps%mpssoang
 gs_hamk%natom      =dtset%natom
 gs_hamk%nattyp(:)  =nattyp(:)
 gs_hamk%nfft       =dtset%nfft
 gs_hamk%ngfft(:)   =dtset%ngfft(:)
 gs_hamk%nloalg(:)  =dtset%nloalg(:)
 gs_hamk%nspinor    =dtset%nspinor
 gs_hamk%ntypat     =psps%ntypat
 gs_hamk%nvloc      =nvloc
 gs_hamk%n4         =n4
 gs_hamk%n5         =n5
 gs_hamk%n6         =n6
 gs_hamk%usepaw     =psps%usepaw
 gs_hamk%use_gpu_cuda=dtset%use_gpu_cuda
 gs_hamk%ph1d(:,:)  =ph1d(:,:)
 gs_hamk%pspso(:)   =psps%pspso(:)
 gs_hamk%ucvol      =ucvol
 gs_hamk%useylm     =psps%useylm
 gs_hamk%xred(:,:)  =xred(:,:)

!Non-local factors:
!Norm-conserving: kleimann-Bylander energies
!PAW: Dij coefficients and overlap coefficients
 if (psps%usepaw==0) then
   gs_hamk%dimekb1=psps%dimekb
   gs_hamk%dimekb2=dtset%ntypat
   ABI_ALLOCATE(gs_hamk%ekb,(psps%dimekb,dtset%ntypat,dtset%nspinor**2))
   ABI_ALLOCATE(gs_hamk%sij,(0,0))
   gs_hamk%ekb(:,:,1)=psps%ekb(:,:)
   if (dtset%nspinor==2) then
     gs_hamk%ekb(:,:,2)=psps%ekb(:,:)
     gs_hamk%ekb(:,:,3:4)=zero
   end if
   if (ipositron==1) gs_hamk%ekb(:,:,:)=-gs_hamk%ekb(:,:,:)
 else
   gs_hamk%dimekb1=psps%dimekb*paw_ij(1)%cplex_dij
   gs_hamk%dimekb2=dtset%natom
   ABI_ALLOCATE(gs_hamk%ekb,(gs_hamk%dimekb1,gs_hamk%dimekb2,dtset%nspinor**2))
   ABI_ALLOCATE(gs_hamk%sij,(gs_hamk%dimekb1,dtset%ntypat))
   ABI_ALLOCATE(dimlmn,(dtset%natom))
   ia=0
   do itypat=1,dtset%ntypat
     if (paw_ij(1)%cplex_dij==1) then
       gs_hamk%sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
     else
       do ilmn=1,pawtab(itypat)%lmn2_size
         gs_hamk%sij(2*ilmn-1,itypat)=pawtab(itypat)%sij(ilmn)
         gs_hamk%sij(2*ilmn  ,itypat)=zero
       end do
     end if
     dimlmn(ia+1:ia+nattyp(itypat))=pawtab(itypat)%lmn_size
     ia=ia+nattyp(itypat)
   end do
   ABI_ALLOCATE(cwaveprj,(dtset%natom,my_nspinor))
   call cprj_alloc(cwaveprj,0,dimlmn)
   if (mpi_enreg%paral_spin==1) then
     ABI_ALLOCATE(cwaveprj_gat,(dtset%natom,dtset%nspinor))
     call cprj_alloc(cwaveprj_gat,0,dimlmn)
   else
     cwaveprj_gat => cwaveprj
   end if
   ABI_DEALLOCATE(dimlmn)
   do iatom=1,dtset%natom
     ABI_ALLOCATE(pawrhoij(iatom)%rhoij_,(pawrhoij(iatom)%cplex*pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
     pawrhoij(iatom)%rhoij_(:,:)=zero
     pawrhoij(iatom)%use_rhoij_=1
   end do
   option_rhoij=1
   usetimerev=(dtset%kptopt>0.and.dtset%kptopt<3)
 end if

!LOOP OVER SPINS
 do isppol=1,dtset%nsppol

!  Rewind kpgsph data file if needed:
   if (dtset%mkmem==0) rewind dtfil%unkg
   if (dtset%mkmem==0.and.psps%useylm==1) rewind dtfil%unylm
   ikg=0

!  Set up local potential vlocal with proper dimensioning, from vtrial
!  Also take into account the spin.
   if(dtset%nspden/=4)then
     if (psps%usepaw==0) then
       call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vtrial,vlocal,2)
     else
       ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
       call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
       call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal,2)
       ABI_DEALLOCATE(cgrvtrial)
     end if
   else
     ABI_ALLOCATE(vlocal_tmp,(n4,n5,n6))
     if (psps%usepaw==0) then
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

!  PAW: retrieve Dij coefficients for this spin component
   if (psps%usepaw==1) then
     do ispden=1,dtset%nspinor**2
       isp=isppol
       if (dtset%nspinor==2) isp=ispden
       do iatom=1,dtset%natom
         dimdij=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
         do ilmn=1,dimdij
           gs_hamk%ekb(ilmn,iatom,ispden)=paw_ij(iatom)%dij(ilmn,isp)
         end do
         if(dimdij+1<=gs_hamk%dimekb1) gs_hamk%ekb(dimdij+1:gs_hamk%dimekb1,iatom,ispden)=zero
       end do
     end do
   end if

!  Loop over k points
   do ikpt=1,dtset%nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     istwf_k=dtset%istwfk(ikpt)
     npw_k=npwarr(ikpt)

     if(mpi_enreg%paral_compil_kpt==1)then
!      Skip this k-point if not the proper processor
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me_distrb))/=0) then
         resid(1+bdtot_index : nband_k+bdtot_index) = zero
         bdtot_index=bdtot_index+nband_k
         cycle
       end if
     end if

!    Continue to initialize the Hamiltonian
     gs_hamk%istwf_k    =istwf_k
     gs_hamk%npw        =npw_k

     ABI_ALLOCATE(eig_k,(nband_k))
     ABI_ALLOCATE(occ_k,(nband_k))
     ABI_ALLOCATE(resid_k,(nband_k))
     ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
     resid_k(:)=0.0_dp
     kpoint(:)=dtset%kptns(:,ikpt)
     gs_hamk%kpoint(:)  =kpoint(:)
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
     eig_k(:)=eigen(1+bdtot_index:nband_k+bdtot_index)
     if (minval(eig_k)>1.d100) eig_k=zero
     cplex=2
     if (istwf_k>1) cplex=1

     if (dtset%mkmem==0) then

!      Read sphere data centered at k in dtfil%unkg, then k+g data
       nsp=dtset%nspinor
       call rdnpw(ikpt,isppol,nband_k,npw_k,nsp,0,dtfil%unkg)
       read (dtfil%unkg) kg_k(1:3,1:npw_k)
       call sphereboundary(gs_hamk%gbound,istwf_k,kg_k,dtset%mgfft,npw_k)

!      Eventually read spherical harmonics
       if (psps%useylm==1) then
         read(dtfil%unylm)
         read(dtfil%unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,psps%mpsang*psps%mpsang)
       end if

!      Read the wavefunction block for ikpt,isppol
       tim_rwwf=3
       ABI_ALLOCATE(eig_dum,(dtset%mband))
       ABI_ALLOCATE(kg_dum,(3,0))
       ABI_ALLOCATE(occ_dum,(dtset%mband))
       call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,dtset%mband,mcg_disk,mpi_enreg,nband_k,&
&       nband_k,npw_k,my_nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
       ABI_DEALLOCATE(eig_dum)
       ABI_DEALLOCATE(kg_dum)
       ABI_DEALLOCATE(occ_dum)

     else

       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       call sphereboundary(gs_hamk%gbound,istwf_k,kg_k,dtset%mgfft,npw_k)
       if (psps%useylm==1) then
         do ilm=1,psps%mpsang*psps%mpsang
           ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
         end do
       end if

!      End if for choice governed by dtset%mkmem
     end if

!    Compute kinetic energy
     ABI_ALLOCATE(kinpw,(npw_k))
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg_k,kinpw,kpoint,npw_k)

     indx=1+icg
     do iband=1,nband_k
!      Compute kinetic energy of each band

       if(mpi_enreg%paral_compil_kpt==1)then
!        Skip this band if not the proper processor
         if (mpi_enreg%paralbd >1) then
           if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/= me_distrb) then
             indx=indx+npw_k*my_nspinor
             cycle
           end if
         end if
       end if

       if (abs(occ_k(iband))>tol8) then
         if(dtset%mkmem/=0)then
           do ig=1,npw_k*my_nspinor
             cwavef(1,ig)=cg(1,indx)
             cwavef(2,ig)=cg(2,indx)
             indx=indx+1
           end do
         else
           do ig=1,npw_k*my_nspinor
             cwavef(1,ig)=cg_disk(1,indx)
             cwavef(2,ig)=cg_disk(2,indx)
             indx=indx+1
           end do
         end if
         call meanvalue_g(ekk,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,cwavef,cwavef,0)
         energies%e_kinetic=energies%e_kinetic+dtset%wtk(ikpt)*occ_k(iband)*ekk
       else
         indx=indx+npw_k*my_nspinor
       end if
     end do

     enlk=zero
     eeigk=zero

!    Compute nonlocal form factors ffnl at all (k+G):
     ider=0;dimffnl=1;nkpg=0
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&     gmet,gprimd,ider,ider,psps%indlmn,kg_k,kpg_dum,kpoint,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&     npw_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&     psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)

!    Allocate the arrays phkxred and ph3d, compute phkxred
!    and eventually ph3d.
     do ia=1,dtset%natom
       iatom=atindx(ia)
       arg=two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
       gs_hamk%phkxred(1,iatom)=cos(arg)
       gs_hamk%phkxred(2,iatom)=sin(arg)
     end do
     if(dtset%nloalg(1)<=0)then
!      Only the allocation, not the precomputation.
       matblk=dtset%nloalg(4)
       ABI_ALLOCATE(ph3d,(2,npw_k,matblk))
     else
!      Here, allocation as well as precomputation
       matblk=dtset%natom
       ABI_ALLOCATE(ph3d,(2,npw_k,matblk))
       call ph1d3d(1,dtset%natom,kg_k,matblk,dtset%natom,npw_k,n1,n2,n3,&
&       gs_hamk%phkxred,ph1d,ph3d)
     end if
     gs_hamk%matblk=matblk

!    DEBUG
!    write(std_out,*)' energy : before nonlop '
!    stop
!    ENDDEBUG

!    Compute nonlocal psp energy (NCPP) or Rhoij (PAW)
     do iband=1,nband_k
       if(mpi_enreg%paral_compil_kpt==1)then
!        Skip this band if not the proper processor
         if (mpi_enreg%paralbd >1) then
           if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/= me_distrb) then
             cycle
           end if
         end if
       end if

!      Select occupied bands
       if (abs(occ_k(iband))>tol8) then

         if(dtset%mkmem/=0)cwavef(:,1:npw_k*my_nspinor)=&
&         cg(:,1+(iband-1)*npw_k*my_nspinor+icg:iband*npw_k*my_nspinor+icg)
         if(dtset%mkmem==0)cwavef(:,1:npw_k*my_nspinor)=&
&         cg_disk(:,1+(iband-1)*npw_k*my_nspinor:iband*npw_k*my_nspinor)

         choice=1-gs_hamk%usepaw ; signs=1 ; idir=0 ; nnlout=1 ; tim_nonlop=3
         paw_opt=gs_hamk%usepaw;cpopt=gs_hamk%usepaw-1
         call nonlop(atindx1,choice,cpopt,cwaveprj,gs_hamk%dimekb1,gs_hamk%dimekb2,dimffnl,dimffnl,&
&         gs_hamk%ekb,enlout,ffnl,ffnl,gmet,gprimd,idir,psps%indlmn,istwf_k,&
&         kg_k,kg_k,kpg_dum,kpg_dum,kpoint,kpoint,dum,psps%lmnmax,&
&         matblk,dtset%mgfft,mpi_enreg,psps%mpsang,&
&         psps%mpssoang,dtset%natom,nattyp,dtset%ngfft,nkpg,nkpg,dtset%nloalg,nnlout,npw_k,npw_k,&
&         my_nspinor,dtset%nspinor,psps%ntypat,0,paw_opt,gs_hamk%phkxred,gs_hamk%phkxred,ph1d,&
&         ph3d,ph3d,signs,nonlop_dum,nonlop_dum,tim_nonlop,ucvol,&
&         psps%useylm,cwavef,cwavef,use_gpu_cuda=dtset%use_gpu_cuda)

         if (psps%usepaw==0) enlk=enlk+occ_k(iband)*enlout(1)
         eeigk=eeigk+occ_k(iband)*eig_k(iband)

!        PAW: accumulate rhoij
         if (psps%usepaw==1) then
           if (mpi_enreg%paral_spin==1) then
             call cprj_gather_spin(cwaveprj,cwaveprj_gat,dtset%natom,1,my_nspinor,dtset%nspinor,&
&             mpi_enreg%comm_spin,ierr)
             call pawaccrhoij(atindx1,cplex,cwaveprj_gat,cwaveprj_gat,0,isppol,dtset%natom,&
&             dtset%nspinor,occ_k(iband),option_rhoij,pawrhoij,usetimerev,dtset%wtk(ikpt))
           else
             call pawaccrhoij(atindx1,cplex,cwaveprj,cwaveprj,0,isppol,dtset%natom,&
&             dtset%nspinor,occ_k(iband),option_rhoij,pawrhoij,usetimerev,dtset%wtk(ikpt))
           end if
         end if

!        End loop on bands
       end if
     end do

     if (psps%usepaw==0) energies%e_nonlocalpsp=energies%e_nonlocalpsp+dtset%wtk(ikpt)*enlk
     energies%e_eigenvalues=energies%e_eigenvalues+dtset%wtk(ikpt)*eeigk

!    DEBUG
!    write(std_out,*)' energy : after nonlop '
!    stop
!    ENDDEBUG

!    Compute residual of each band (for informative purposes)
     if(dtset%mkmem/=0)then
       call mkresi(cg,dimffnl,eig_k,ffnl,dtfil%filstat,&
&       gs_hamk,icg,ikpt,isppol,&
&       kg_k,kinpw,psps%lmnmax,matblk,mcg,dtset%mgfft,mpi_enreg,&
&       psps%mpsang,psps%mpssoang,dtset%natom,&
&       nband_k,npw_k,&
&       my_nspinor,psps%ntypat,nvloc,n4,n5,n6,&
&       dtset%paral_kgb,ph3d,dtset%prtvol,&
&       resid_k,psps%usepaw,vlocal)
     else if(dtset%mkmem==0)then
       call mkresi(cg_disk,dimffnl,eig_k,ffnl,dtfil%filstat,&
&       gs_hamk,icg,ikpt,isppol,&
&       kg_k,kinpw,psps%lmnmax,matblk,mcg_disk,dtset%mgfft,mpi_enreg,&
&       psps%mpsang,psps%mpssoang,dtset%natom,&
&       nband_k,npw_k,&
&       my_nspinor,psps%ntypat,nvloc,n4,n5,n6,&
&       dtset%paral_kgb,ph3d,dtset%prtvol,&
&       resid_k,psps%usepaw,vlocal)
     end if
     resid(1+bdtot_index : nband_k+bdtot_index) = resid_k(:)

     ABI_DEALLOCATE(eig_k)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(resid_k)

!    DEBUG
!    write(std_out,*)' isppol,ikpt',isppol,ikpt
!    write(std_out,*)resid(1+bdtot_index:nband_k+bdtot_index)
!    ENDDEBUG

     bdtot_index=bdtot_index+nband_k

     if (dtset%mkmem/=0) then
!      Handle case in which kg, cg, are kept in core
       icg=icg+npw_k*my_nspinor*nband_k
       ikg=ikg+npw_k
     end if

     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(kinpw)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(ylm_k)

!    End loops on isppol and ikpt
   end do
 end do

 if(mpi_enreg%paral_compil_kpt==1)then
   call leave_test()
   write(message,*) 'energy: loop on k-points and spins done in parallel'
   call wrtout(std_out,message,'COLL')
 end if

 ABI_DEALLOCATE(gs_hamk%atindx)
 ABI_DEALLOCATE(gs_hamk%atindx1)
 ABI_DEALLOCATE(gs_hamk%ekb)
 ABI_DEALLOCATE(gs_hamk%sij)
 ABI_DEALLOCATE(gs_hamk%gbound)
 ABI_DEALLOCATE(gs_hamk%indlmn)
 ABI_DEALLOCATE(gs_hamk%nattyp)
 ABI_DEALLOCATE(gs_hamk%phkxred)
 ABI_DEALLOCATE(gs_hamk%pspso)
 ABI_DEALLOCATE(gs_hamk%ph1d)
 ABI_DEALLOCATE(gs_hamk%xred)

 if(dtset%mkmem==0)  then
   ABI_DEALLOCATE(cg_disk)
 end if

!DEBUG
!write(std_out,*)' energy : after loop on kpts and spins '
!stop
!ENDDEBUG

 if(mpi_enreg%paral_compil_kpt==1)then
!  Accumulate enl eeig and ek on all proc.
   ABI_ALLOCATE(buffer,(3+dtset%mband*dtset%nkpt*dtset%nsppol))
   buffer(1)=energies%e_nonlocalpsp ; buffer(2)=energies%e_kinetic ; buffer(3)=energies%e_eigenvalues
   do iresid=1,dtset%mband*dtset%nkpt*dtset%nsppol
     buffer(iresid+3)=resid(iresid)
   end do
   call timab(48,1,tsec)
   call xsum_mpi(buffer,spaceComm,ierr)
   call timab(48,2,tsec)
   energies%e_nonlocalpsp=buffer(1) ; energies%e_kinetic=buffer(2) ; energies%e_eigenvalues=buffer(3)
   do iresid=1,dtset%mband*dtset%nkpt*dtset%nsppol
     resid(iresid)=buffer(iresid+3)
   end do
   ABI_DEALLOCATE(buffer)
!  Accumulate rhoij_
   if (psps%usepaw==1) then
     call timab(48,1,tsec)
     ABI_ALLOCATE(dimlmn,(dtset%natom))
     dimlmn(1:dtset%natom)=pawrhoij(1:dtset%natom)%cplex*pawrhoij(1:dtset%natom)%lmn2_size
     nsp2=pawrhoij(1)%nsppol
     if (pawrhoij(1)%nspden==4) nsp2=4
     bufdim=sum(dimlmn)*nsp2
     ABI_ALLOCATE(buffer,(bufdim))
     ABI_ALLOCATE(buffer2,(bufdim))
     ii=0
     do iatom=1,dtset%natom
       do isppol=1,nsp2
         buffer(ii+1:ii+dimlmn(iatom))=pawrhoij(iatom)%rhoij_(1:dimlmn(iatom),isppol)
         ii=ii+dimlmn(iatom)
       end do
     end do
     call xsum_mpi(buffer,buffer2,bufdim,spaceComm,ierr) !Build sum of everything
     ii=0
     do iatom=1,dtset%natom
       do isppol=1,nsp2
         pawrhoij(iatom)%rhoij_(1:dimlmn(iatom),isppol)=buffer2(ii+1:ii+dimlmn(iatom))
         ii=ii+dimlmn(iatom)
       end do
     end do
     ABI_DEALLOCATE(buffer)
     ABI_DEALLOCATE(buffer2)
     ABI_DEALLOCATE(dimlmn)
     call timab(48,2,tsec)
   end if
 end if

!Compute total (free) energy
 if (optene==0.or.optene==2) then
   if (psps%usepaw==0) then
     etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
&     energies%e_localpsp + energies%e_corepsp + energies%e_nonlocalpsp
   else
     etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
&     energies%e_localpsp + energies%e_corepsp + energies%e_paw
   end if
 else if (optene==1.or.optene==3) then
   if (psps%usepaw==0) then
     etotal = energies%e_eigenvalues - energies%e_hartree + energies%e_xc - &
&     energies%e_xcdc + energies%e_corepsp
   else
     etotal = energies%e_eigenvalues - energies%e_hartree + energies%e_xc - &
&     energies%e_xcdc + energies%e_corepsp + energies%e_pawdc
   end if
 end if
 etotal = etotal + energies%e_ewald
 if(dtset%occopt>=3 .and. dtset%occopt<=8) etotal=etotal-dtset%tsmear*energies%entropy

!Additional stuff for electron-positron
 if (dtset%positron/=0) then
   if (ipositron==0) then
     energies%e_electronpositron  =zero
     energies%edc_electronpositron=zero
   else
     energies%e_electronpositron  =electronpositron%e_hartree+electronpositron%e_xc
     energies%edc_electronpositron=electronpositron%e_hartree+electronpositron%e_xcdc
     if (psps%usepaw==1) then
       energies%e_electronpositron  =energies%e_electronpositron  +electronpositron%e_paw
       energies%edc_electronpositron=energies%edc_electronpositron+electronpositron%e_pawdc
     end if
   end if
   if (optene==0.or.optene==2) electronpositron%e0=etotal
   if (optene==1.or.optene==3) electronpositron%e0=etotal-energies%edc_electronpositron
   etotal=electronpositron%e0+energies%e0_electronpositron+energies%e_electronpositron
 end if

!Compute new charge density based on incoming wf
!Keep rhor and rhog intact for later use e.g. in stress. (? MT 08-12-2008: is that true now ?)
!=== Norm-conserving psps: simply compute rho from WFs
 paw_dmft%use_dmft=0 ! dmft not used here
 if (psps%usepaw==0) then
   tim_mkrho=3
   call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,&
&   npwarr,occ,paw_dmft,phnons,rhog,rhor,rprimd,tim_mkrho,ucvol,dtfil%unkg,wffnow,wfs, wvl)
   if(dtset%usekden==1)then
     call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,&
&     npwarr,occ,paw_dmft,phnons,taug,taur,rprimd,tim_mkrho,ucvol,dtfil%unkg,wffnow,wfs,wvl,option=1)
   end if
 else
!  === PAW case: add compensation charge density
   tim_mkrho=3;option=1;choice=1
   call symrhoij(choice,gprimd,psps%indlmn,indsym,0,psps%lmnmax,dtset%natom,dtset%natom,dtset%nsym,dtset%ntypat,option,&
&   pawang,dtset%pawprtvol,pawrhoij,rprimd,dtset%symafm,symrec,dtset%typat)
   do iatom=1,dtset%natom
     ABI_DEALLOCATE(pawrhoij(iatom)%rhoij_)
     pawrhoij(iatom)%use_rhoij_=0
   end do
   ider=0;izero=0;cplex_rf=1;ipert=0;idir=0;qpt(:)=zero
   call pawmknhat(compch_fft,cplex_rf,ider,idir,ipert,izero,gprimd,mpi_enreg,dtset%natom,dtset%natom,nfftf,ngfftf,&
&   0,dtset%nspden,dtset%ntypat,dtset%paral_kgb,pawang,pawfgrtab,rhodum,nhat,pawrhoij,pawrhoij,&
&   pawtab,qpt,rprimd,ucvol,xred)
   ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
   ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
   rhowfr(:,:)=zero
   call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,&
&   npwarr,occ,paw_dmft,phnons,rhowfg,rhowfr,rprimd,tim_mkrho,ucvol,dtfil%unkg,wffnow,wfs,wvl)
   call transgrid(1,mpi_enreg,dtset%nspden,+1,1,0,dtset%paral_kgb,pawfgr,rhowfg,rhodum,rhowfr,rhor)
   ABI_DEALLOCATE(rhowfr)
   ABI_DEALLOCATE(rhowfg)
   rhor(:,:)=rhor(:,:)+nhat(:,:)
   call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
 end if

 write(message, '(a,a,a,a)' )ch10, &
& ' energy: COMMENT -',ch10,&
& '  New density rho(r) made from input wfs'
 call wrtout(std_out,message,'COLL')

 call timab(59,2,tsec)

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(vlocal)
 if (psps%usepaw==1) then
   call cprj_free(cwaveprj)
   ABI_DEALLOCATE(cwaveprj)
   if (mpi_enreg%paral_spin==1) then
     call cprj_free(cwaveprj_gat)
     ABI_DEALLOCATE(cwaveprj_gat)
   else
     nullify(cwaveprj_gat)
   end if
 end if

end subroutine energy
!!***
