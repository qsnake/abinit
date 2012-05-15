!{\src2tex{textfont=tt}}
!!****f* ABINIT/forstr
!! NAME
!! forstr
!!
!! FUNCTION
!! Drives the computation of forces and/or stress tensor
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MB, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cg(2,mcg)=wavefunctions (may be read from disk instead of input)
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!       and electric field calculations (see initberry.f).
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryopt  = 4: electric field is on -> add the contribution of the
!!   |                - \Omega E.P term to the total energy
!!   |          /= 4: electric field is off
!!   |  from Etot(npw) data (at fixed geometry), used for making
!!   |  Pulay correction to stress tensor (hartree).  Should be <=0.
!!   | ecut=cut-off energy for plane wave basis sphere (Ha)
!!   | ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!   | effmass=effective mass for electrons (1. in common case)
!!   | efield = cartesian coordinates of the electric field in atomic units
!!   | ionmov=governs the movement of atoms (see help file)
!!   | iprcch=governs the mixed electronic-atomic part of the preconditioner
!!   | istwfk(nkpt)=input option parameter that describes the storage of wfs
!!   | kptns(3,nkpt)=reduced coordinates of k points in Brillouin zone
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem=maximum number of k points in core memory
!!   | mpw = maximum number of plane waves
!!   | natom=number of atoms in cell
!!   | nband(nkpt*nsppol)=number of bands to be included in summation at each k point
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!   | nkpt=number of k points in Brillouin zone
!!   | nloalg(5)=governs the choice of the algorithm for non-local operator.
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | pawprtvol=control print volume and debugging output for PAW
!!   | prtvol=integer controlling volume of printed output
!!   | symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!   | tfkinfunc=1 if use of Thomas-Fermi kinetic functional
!!   |          =2 if use of recursion method
!!   | typat(natom)=type integer for each atom in cell
!!   | wtk(nkpt)=weights associated with various k points
!!   | nsym=number of symmetries in space group
!!  energies <type(energies_type)>=all part of total energy.
!!   | e_localpsp(IN)=local psp energy (hartree)
!!   | e_hartree(IN)=Hartree part of total energy (hartree units)
!!   | e_corepsp(IN)=psp core-core energy
!!   | e_kinetic(IN)=kinetic energy part of total energy.
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  grewtn(3,natom)=d(Ewald)/d(xred) (hartree)
!!  gsqcut=cutoff value on G**2 for (large) sphere inside FFT box.
!!                       gsqcut=(boxcut**2)*ecut/(2._dp*(Pi**2)
!!  indsym(4,nsym,natom)=index showing transformation of atom labels
!!                       under symmetry operations (computed in symatm)
!!  kg(3,mpw*mkmem)=reduced (integer) coordinates of G vecs in basis sphere
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mgfftf= -PAW ONLY- maximum size of 1D FFTs for the fine grid
!!         (mgfftf=mgfft for norm-conserving potential runs)
!!  mpi_enreg=informations about MPI parallelization
!!  n3xccc=dimension of the xccc3d array (0 or nfftf).
!!  nattyp(ntypat)=number of atoms of each type
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  ngfftf(18)= -PAW ONLY- contain all needed information about 3D FFT for the fine grid
!!              (ngfftf=ngfft for norm-conserving potential runs)
!!  nhat(nfftf,nspden*psps%usepaw)= -PAW only- compensation density
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  npwarr(nkpt)=number of planewaves in basis and on boundary for each k
!!  ntypat=number of types of atoms
!!  nvresid(nfftf,nspden)=array for the residual of the density/potential
!!  occ(mband*nkpt*nsppol)=occupancies of bands at various k points
!!  optfor=1 if computation of forces is required
!!  optres=0 if the potential residual has to be used for forces corrections
!!        =1 if the density residual has to be used for forces corrections
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr(natom*usepaw) <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pel(3)=reduced coordinates of the electronic polarization (a. u.)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases
!!  ph1df(2,3*(2*mgfftf+1)*natom)=-PAW only- 1-dim structure factor phases for the fine grid
!!  pion(3)=reduced coordinates of the ionic polarization (a. u.)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum
!!  rhog(2,nfftf)=Fourier transform of charge density (bohr^-3)
!!  rhor(nfftf,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  strsxc(6)=xc correction to stress
!!  stress_needed=1 if computation of stress tensor is required
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  ucvol=unit cell volume in bohr**3
!!  unkg=unit number for (k+G) sphere data file
!!  unylm=unit number for Ylm(k) data (if used)
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vhartr(nfftf)=array for holding Hartree potential
!!  vpsp(nfftf)=array for holding local psp
!!  vxc(nfftf,nspden)=exchange-correlation potential (hartree) in real space
!!  wffnow=unit number for current wf disk file
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  ==== if (optfor==1) ====
!!   diffor=maximal absolute value of changes in the components of
!!          force between the input and the output.
!!   favg(3)=mean of the forces before correction for translational symmetry
!!   fcart(3,natom)=forces in cartesian coordinates (Ha/Bohr)
!!     at input, previous value of forces,
!!     at output, new value.
!!     Note : unlike fred, this array has been corrected by enforcing
!!     the translational symmetry, namely that the sum of force
!!     on all atoms is zero.
!!   forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!   fred(3,natom)=symmetrized grtn = d(etotal)/d(xred)
!!   gresid(3,natom)=forces due to the residual of the density/potential
!!   grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!   grxc(9+3*natom)=d(Exc)/d(xred) if core charges are used
!!   maxfor=maximal absolute value of the output array force.
!!   synlgr(3,natom)=symmetrized gradients of energy due to nonlocal contributions
!!  ==== if (stress_needed==1) ====
!!   strten(6)=components of the stress tensor (hartree/bohr^3) for the
!!    6 unique components of this symmetric 3x3 tensor:
!!    Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!  ===== if psps%usepaw==1
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!    (gradients of rhoij for each atom with respect to atomic positions are computed here)
!!
!! NOTES
!!  Be careful to the meaning of nfft (size of FFT grids):
!!   - In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!   - In case of PAW calculations:
!!     Two FFT grids are used; one with nfft points (coarse grid) for
!!     the computation of wave functions ; one with nfftf points
!!     (fine grid) for the computation of total density.
!!
!! PARENTS
!!      afterscfloop,setup_positron
!!
!! CHILDREN
!!      forces,forstrnps,leave_new,nres2vres,pawgrnl,stress,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine forstr(atindx1,cg,diffor,dtefield,dtset,eigen,electronpositron,energies,favg,fcart,&
&                 forold,fred,gresid,grewtn,grhf,grxc,gsqcut,indsym,&
&                 kg,kxc,maxfor,mcg,mgfftf,mpi_enreg,n3xccc,nattyp,&
&                 nfftf,ngfftf,nhat,nkxc,npwarr,&
&                 ntypat,nvresid,occ,optfor,optres,paw_ij,pawang,pawfgr,&
&                 pawfgrtab,pawrhoij,pawtab,pel,ph1d,ph1df,pion,psps,rhog,rhor,rprimd,stress_needed,&
&                 strsxc,strten,symrec,synlgr,ucvol,unkg,unylm,usexcnhat,vhartr,vpsp,&
&                 vxc,wffnow,wvl,xccc3d,xred,ylm,ylmgr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_efield
 use m_energies, only : energies_type
 use m_wffile
 use m_electronpositron, only : electronpositron_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'forstr'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_66_paw
 use interfaces_67_common, except_this_one => forstr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,mgfftf,n3xccc,nfftf,nkxc,ntypat,optfor,optres
 integer,intent(in) :: stress_needed,unkg,unylm,usexcnhat
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(inout) :: diffor,maxfor
 type(efield_type), intent(in) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(energies_type),intent(in) :: energies
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 integer,intent(in) :: atindx1(dtset%natom),indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),nattyp(ntypat),ngfftf(18)
 integer,intent(in) :: npwarr(dtset%nkpt),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: grewtn(3,dtset%natom),kxc(dtset%nfft,nkxc)
 real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol),pel(3)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
 real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*dtset%natom),pion(3)
 real(dp),intent(in) :: rhog(2,nfftf),rprimd(3,3),strsxc(6),vhartr(nfftf)
 real(dp),intent(in) :: vpsp(nfftf),vxc(nfftf,dtset%nspden)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: forold(3,dtset%natom)
 real(dp),intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
 real(dp),intent(inout) :: nvresid(nfftf,dtset%nspden),rhor(nfftf,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp),intent(out) :: favg(3),fcart(3,dtset%natom),fred(3,dtset%natom)
 real(dp),intent(out) :: gresid(3,dtset%natom),grhf(3,dtset%natom)
 real(dp),intent(out) :: grxc(3,dtset%natom),strten(6),synlgr(3,dtset%natom)
 type(paw_ij_type),intent(in) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: dimvcur,ifft,ispden,occopt_,optgr,optgr2,option,optnc,optstr
 character(len=500) :: message
!arrays
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: kinstr(6),nlstr(6),tsec(2)
 real(dp),allocatable :: dummy(:),grnl(:),vcurrent(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' forstr : enter '
!if(.true.)stop
!ENDDEBUG

 call timab(910,1,tsec)
 call timab(911,1,tsec)

!Do nothing if nothing is required
 if (optfor==0.and.stress_needed==0) return

!Test size of FFT grids (1 grid in norm-conserving, 2 grids in PAW)
 if ((psps%usepaw==1.and.pawfgr%nfft/=nfftf).or.(psps%usepaw==0.and.dtset%nfft/=nfftf)) then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' forstr :  BUG -',ch10,&
&   '  wrong values for nfft, nfftf !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if ((psps%usepaw==1.and.pawfgr%mgfft/=mgfftf).or.(psps%usepaw==0.and.dtset%mgfft/=mgfftf)) then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' forstr :  BUG -',ch10,&
&   '  wrong values for mgfft, mgfftf !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!==========================================================================
!Here compute terms common to forces and stresses
!==========================================================================

 if (optfor==1)  then
   ABI_ALLOCATE(grnl,(3*dtset%natom))
 end if

!Compute nonlocal pseudopotential parts of forces and stress tensor
!-involves summations over wavefunctions at all k points
 if ((dtset%tfkinfunc==1.or.dtset%tfkinfunc==2).and.stress_needed==1) then
   kinstr(1:3)=-two/three*energies%e_kinetic/ucvol ; kinstr(4:6)=zero
   nlstr(1:6)=zero
   write(std_out,*)'kinstr tf',kinstr(1)
 else
   occopt_=0 ! This means that occ are now fixed
   call forstrnps(atindx1,cg,dtefield,dtset%ecut,dtset%ecutsm,dtset%effmass,eigen,electronpositron,grnl,&
&   dtset%istwfk,kg,kinstr,nlstr,dtset%kptns,dtset%mband,mcg,dtset%mgfft,dtset%mkmem,&
&   mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,nattyp,dtset%nband,dtset%ngfft,dtset%nkpt,&
&   dtset%nloalg,npwarr,dtset%nspinor,dtset%nsppol,dtset%nsym,ntypat,occ,&
&   optfor,paw_ij,pawtab,ph1d,psps,rprimd,stress_needed,symrec,dtset%typat,&
&   unkg,unylm,dtset%use_gpu_cuda,wffnow,dtset%wtk,xred,ylm,ylmgr)
 end if

 call timab(911,2,tsec)
 call timab(912,1,tsec)

!PAW: add gradients due to Dij derivatives to non-local term
 if (psps%usepaw==1) then
   if (usexcnhat==0) then
     dimvcur=1
     ABI_ALLOCATE(vcurrent,(nfftf,dimvcur))
!    $OMP PARALLEL DO PRIVATE(ifft) &
!    $OMP&SHARED(vcurrent,nfftf,vhartr,vpsp)
     do ifft=1,nfftf
       vcurrent(ifft,1)=vhartr(ifft)+vpsp(ifft)
     end do
!    $OMP END PARALLEL DO
   else
     dimvcur=dtset%nspden
     ABI_ALLOCATE(vcurrent,(nfftf,dimvcur))
     do ispden=1,min(dtset%nspden,2)
!      $OMP PARALLEL DO PRIVATE(ifft) &
!      $OMP&SHARED(ispden,vcurrent,nfftf,vhartr,vpsp,vxc)
       do ifft=1,nfftf
         vcurrent(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)
       end do
!      $OMP END PARALLEL DO
     end do
     if (dtset%nspden==4) then
       do ispden=3,4
!        $OMP PARALLEL DO PRIVATE(ifft) &
!        $OMP&SHARED(ispden,vcurrent,nfftf,vhartr,vpsp,vxc)
         do ifft=1,nfftf
           vcurrent(ifft,ispden)=vxc(ifft,ispden)
         end do
!        $OMP END PARALLEL DO
       end do
     end if
   end if
   optgr=optfor;optgr2=0;optstr=stress_needed
   call pawgrnl(atindx1,dtset%nspden,dimvcur,dummy,1,grnl,gsqcut,mgfftf,mpi_enreg,dtset%natom,dtset%natom,&
&   nattyp,nfftf,ngfftf,nhat,nlstr,dtset%nspden,dtset%nsym,ntypat,optgr,optgr2,optstr,dtset%paral_kgb,&
&   pawang,pawfgrtab,pawrhoij,pawtab,ph1df,psps,k0,rprimd,symrec,dtset%typat,vcurrent,xred)
   ABI_DEALLOCATE(vcurrent)
 end if

 call timab(912,2,tsec)
 call timab(913,1,tsec)

!==========================================================================
!Here compute forces (if required)
!==========================================================================
 if (optfor==1) then
!  If residual is a density residual (and forces from residual asked),
!  has to convert it into a potential residualbefore calling forces routine
   if (optres==1 .and. dtset%usewvl==0.and.abs(dtset%iprcch)>=1 .and. &
&   abs(dtset%iprcch)<=6.and.abs(dtset%iprcch)/=5) then
     option=0; if (dtset%iprcch<0) option=1
     ABI_ALLOCATE(vcurrent,(nfftf,dtset%nspden))
     optnc=1;if (dtset%nspden==4.and.(abs(dtset%iprcch)==4.or.abs(dtset%iprcch)==6)) optnc=2
     call nres2vres(dtset,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,ngfftf,nhat,&
&     nkxc,nvresid,n3xccc,optnc,option,pawang,pawfgrtab,pawrhoij,pawtab,&
&     rhor,rprimd,psps%usepaw,usexcnhat,vcurrent,wvl,xccc3d,xred)
     call forces(atindx1,diffor,dtset,favg,fcart,forold,fred,gresid,grewtn,&
&     grhf,grnl,grxc,gsqcut,indsym,maxfor,mgfftf,&
&     mpi_enreg,psps%n1xccc,n3xccc,nattyp,&
&     nfftf,ngfftf,ntypat,pawtab,ph1df,psps,rhog,&
&     rhor,rprimd,symrec,synlgr,vcurrent,vxc,wvl,xred,&
&     electronpositron=electronpositron)
     ABI_DEALLOCATE(vcurrent)
   else
     call forces(atindx1,diffor,dtset,favg,fcart,forold,fred,gresid,grewtn,&
&     grhf,grnl,grxc,gsqcut,indsym,maxfor,mgfftf,&
&     mpi_enreg,psps%n1xccc,n3xccc,nattyp,&
&     nfftf,ngfftf,ntypat,pawtab,ph1df,psps,rhog,&
&     rhor,rprimd,symrec,synlgr,nvresid,vxc,wvl,xred,&
&     electronpositron=electronpositron)
   end if
 end if

 call timab(913,2,tsec)
 call timab(914,1,tsec)

!==========================================================================
!Here compute stress tensor (if required)
!==========================================================================
 if (stress_needed==1) then
   call stress(atindx1,dtset%berryopt,energies%e_localpsp,dtset%efield,&
&   energies%e_hartree,energies%e_corepsp,gsqcut,kinstr,mgfftf,&
&   mpi_enreg,psps%mqgrid_vl,psps%n1xccc,n3xccc,dtset%natom,nattyp,&
&   nfftf,ngfftf,nlstr,dtset%nspden,dtset%nsym,ntypat,dtset%paral_kgb,pawtab,pel,pion,ph1df,&
&   dtset%prtvol,psps%qgrid_vl,rhog,rprimd,strten,strsxc,symrec,dtset%typat,psps%usepaw,&
&   psps%vlspl,vxc,psps%xccc1d,xccc3d,psps%xcccrc,xred,psps%ziontypat,&
&   electronpositron=electronpositron)
 end if

!Memory deallocation
 if (optfor==1)  then
   ABI_DEALLOCATE(grnl)
 end if

!DEBUG
!write(std_out,*)' forstr : exit '
!if(.true.)stop
!ENDDEBUG

 call timab(914,2,tsec)
 call timab(910,2,tsec)

end subroutine forstr
!!***
