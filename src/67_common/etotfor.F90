!{\src2tex{textfont=tt}}
!!****f* ABINIT/etotfor
!! NAME
!! etotfor
!!
!! FUNCTION
!! This routine is called to compute the total energy and various parts of it.
!! The routine computes -if requested- the forces.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryopt  = 4: electric field is on -> add the contribution of the
!!   |                - \Omega E.P term to the total energy
!!   |          /= 4: electric field is off
!!   | bfield = cartesian coordinates of magnetic field in atomic units
!!   | efield = cartesian coordinates of the electric field in atomic units
!!   | iatfix(3,natom)=1 for frozen atom along some direction, 0 for unfrozen
!!   | ionmov=governs the movement of atoms (see help file)
!!   | iprcch=governs the mixed electronic-atomic part of the preconditioner
!!   | natom=number of atoms in cell.
!!   | nconeq=number of atomic constraint equations
!!   | nspden=number of spin-density components
!!   | nsym=number of symmetry elements in space group
!!   | occopt=option for occupancies
!!   | prtvol=integer controlling volume of printed output
!!   | tsmear=smearing energy or temperature (if metal)
!!   | typat(natom)=type integer for each atom in cell
!!   | wtatcon(3,natom,nconeq)=weights for atomic constraints
!!  efield_dot = reciprocal lattice coordinates of the electric field
!!  grewtn(3,natom)=grads of Ewald energy (hartree)
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!!  mag_cart = reduced coordinates of the orbital magnetization
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  ntypat=number of types of atoms in unit cell.
!!  nvresid(nfft,nspden)=potential or density residual
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  optene=option for the computation of total energy
!!         (-1=no computation; 0=direct scheme; 1=double-counting scheme)
!!  optforces=option for the computation of forces
!!  optres=0 if residual array (nvresid) contains the potential residual
!!        =1 if residual array (nvresid) contains the density residual
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pel = reduced coordinates of the electronic polarization
!!        (pel does not take into account the factor 1/ucvol)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase (structure factor) information.
!!  pion = reduced coordinates of the ionic polarization
!!        (pel does not take into account the factor 1/ucvol)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=array for Fourier transform of electron density
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vhartr(nfft)=array for holding Hartree potential
!!  vpsp(nfft)=array for holding local psp
!!  vxc(nfft,nspden)=array for holding XC potential
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  deltae=change in total energy between the previous and present SCF cycle
!!  etotal=total energy (hartree)
!!  ===== if optforces==1
!!   diffor=maximum absolute change in component of forces between present and previous SCF cycle.
!!   favg(3)=mean of fcart before correction for translational symmetry
!!   fcart(3,natom)=cartesian forces from fred (hartree/bohr)
!!   fred(3,natom)=symmetrized form of grtn (grads of Etot) (hartree)
!!   gresid(3,natom)=forces due to the residual of the density/potential
!!   grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!   grxc(3,natom)=d(Exc)/d(xred) derivatives (0 without core charges)
!!   maxfor=maximum absolute value of force
!!   synlgr(3,natom)=symmetrized form of grads of Enl (hartree)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  elast=previous value of the energy,
!!        needed to compute deltae, then updated.
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  energies <type(energies_type)>=all part of total energy.
!!   | entropy(IN)=entropy due to the occupation number smearing (if metal)
!!   | e_localpsp(IN)=local psp energy (hartree)
!!   | e_eigenvalues(IN)=Sum of the eigenvalues - Band energy (Hartree)
!!   | e_ewald(IN)=Ewald energy (hartree)
!!   | e_hartree(IN)=Hartree part of total energy (hartree units)
!!   | e_corepsp(IN)=psp core-core energy
!!   | e_kinetic(IN)=kinetic energy part of total energy.
!!   | e_nonlocalpsp(IN)=nonlocal pseudopotential part of total energy.
!!   | e_xc(IN)=exchange-correlation energy (hartree)
!!   | e_xcdc(IN)=exchange-correlation double-counting energy (hartree)
!!   | e_paw(IN)=PAW spherical part energy
!!   | e_pawdc(IN)=PAW spherical part double-counting energy
!!   | e_elecfield(OUT)=the term of the energy functional that depends explicitely
!!   |                  on the electric field:  enefield = -ucvol*E*P
!!   | e_magfield(OUT)=the term of the energy functional that depends explicitely
!!   |                  on the magnetic field:  e_magfield = -ucvol*E*P
!!   | e_entropy(OUT)=entropy energy due to the occupation number smearing (if metal)
!!   |                this value is %entropy * dtset%tsmear (hartree).
!!  ===== if optforces==1
!!   forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!   grnl(3*natom)=gradients of Etot due to nonlocal contributions
!!                 Input for norm-conserving psps, output for PAW
!!  ===== if psps%usepaw==1
!!   pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!    (gradients of rhoij for each atom with respect to atomic positions are computed here)
!!
!! NOTES
!!  In case of PAW calculations:
!!    All computations are done on the fine FFT grid.
!!    All variables (nfft,ngfft,mgfft) refer to this fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!  ! Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      forces,nres2vres,pawgrnl,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine etotfor(atindx1,deltae,diffor,dtset,efield_dot,&
&  elast,electronpositron,energies,&
&  etotal,favg,fcart,forold,fred,gresid,grewtn,grhf,grnl,&
&  grxc,gsqcut,indsym,kxc,mag_cart,maxfor,mgfft,mpi_enreg,nattyp,&
&  nfft,ngfft,nhat,nkxc,ntypat,nvresid,n1xccc,n3xccc,optene,optforces,optres,&
&  pawang,pawfgrtab,pawrhoij,pawtab,pel,ph1d,pion,psps,rhog,rhor,rprimd,symrec,synlgr,&
&  usepaw,usexcnhat,vhartr,vpsp,vxc,wvl,xccc3d,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_energies, only : energies_type
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'etotfor'
 use interfaces_18_timing
 use interfaces_66_paw
 use interfaces_67_common, except_this_one => etotfor
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,n1xccc,n3xccc,nfft,nkxc,ntypat,optene,optforces
 integer,intent(in) :: optres,usepaw,usexcnhat
 real(dp),intent(in) :: gsqcut
 real(dp),intent(inout) :: elast
 real(dp),intent(out) :: deltae,diffor,etotal,maxfor
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 integer,intent(in) :: atindx1(dtset%natom),indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: efield_dot(3),grewtn(3,dtset%natom),kxc(nfft,nkxc)
 real(dp),intent(in) :: mag_cart(3),pel(3),ph1d(2,3*(2*mgfft+1)*dtset%natom),pion(3)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: vhartr(nfft),vpsp(nfft),vxc(nfft,dtset%nspden)
 real(dp),intent(in) :: xccc3d(n3xccc)
 real(dp),intent(inout) :: forold(3,dtset%natom),grnl(3*dtset%natom)
 real(dp),intent(inout) :: nhat(nfft,dtset%nspden*psps%usepaw)
 real(dp),intent(inout) :: nvresid(nfft,dtset%nspden),xred(3,dtset%natom)
 real(dp),intent(out) :: favg(3),fcart(3,dtset%natom),fred(3,dtset%natom)
 real(dp),intent(out) :: gresid(3,dtset%natom),grhf(3,dtset%natom)
 real(dp),intent(out) :: grxc(3,dtset%natom),synlgr(3,dtset%natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: dimnhat,dimwork,ifft,ipositron,ispden,optgr,optgr2,option,optnc,optstr
 real(dp) :: emag,epel,epion
!arrays
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: tsec(2)
 real(dp),allocatable :: dummy(:),nhat_dum(:,:),work(:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)' etotfor : enter'
!ENDDEBUG

 call timab(80,1,tsec)

 ipositron=electronpositron_calctype(electronpositron)

 if (optene>-1) then

!  When the finite-temperature VG broadening scheme is used,
!  the total entropy contribution "tsmear*entropy" has a meaning,
!  and gather the two last terms of Eq.8 of the VG paper
!  Warning : might have to be changed for fixed moment calculations
   if(dtset%occopt>=3 .and. dtset%occopt<=8) then
     if (abs(dtset%tphysel) < tol10) then
       energies%e_entropy = - dtset%tsmear * energies%entropy
     else
       energies%e_entropy = - dtset%tphysel * energies%entropy
     end if
   else
     energies%e_entropy = zero
   end if

!  Turn it into an electric enthalpy, by adding both ionic and electronic contributions
   energies%e_elecfield = zero
   if (dtset%berryopt == 4 .and. ipositron/=1) then
!    First, ionic contribution (epion):
     epion = -1_dp*(efield_dot(1)*pion(1) + efield_dot(2)*pion(2) + &
&     efield_dot(3)*pion(3))
     energies%e_elecfield = energies%e_elecfield + epion
!    Now, electronic contribution (epel):
     epel = -1_dp*(efield_dot(1)*pel(1) + efield_dot(2)*pel(2) + &
&     efield_dot(3)*pel(3))
     energies%e_elecfield = energies%e_elecfield + epel
   end if

!  Turn it into an magnetic enthalpy, by adding orbital electronic contribution
   energies%e_magfield = zero
   if (dtset%berryopt == 5 .and. ipositron/=1) then
     emag = dot_product(mag_cart,dtset%bfield)
     energies%e_magfield = emag
   end if


!  Compute total (free)- energy by direct scheme
   if (optene==0) then
     etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
&     energies%e_localpsp + energies%e_corepsp + &
&     energies%e_entropy + energies%e_elecfield + energies%e_magfield
     etotal = etotal + energies%e_ewald
     if (usepaw==0) etotal = etotal + energies%e_nonlocalpsp
     if (usepaw/=0) etotal = etotal + energies%e_paw
     if (dtset%usewvl == 1) etotal = etotal - energies%e_vxc
   end if

!  Compute total (free) energy by double-counting scheme
   if (optene==1) then
     etotal = energies%e_eigenvalues - energies%e_hartree + energies%e_xc - &
&     energies%e_xcdc + energies%e_corepsp + &
&     energies%e_entropy + energies%e_elecfield + energies%e_magfield
     etotal = etotal + energies%e_ewald
     if (usepaw/=0) etotal = etotal + energies%e_pawdc
   end if

!  Additional stuff for electron-positron
   if (dtset%positron/=0) then
     if (ipositron==0) then
       energies%e_electronpositron  =zero
       energies%edc_electronpositron=zero
     else
       energies%e_electronpositron  =electronpositron%e_hartree+electronpositron%e_xc
       energies%edc_electronpositron=electronpositron%e_hartree+electronpositron%e_xcdc
       if (usepaw==1) then
         energies%e_electronpositron  =energies%e_electronpositron  +electronpositron%e_paw
         energies%edc_electronpositron=energies%edc_electronpositron+electronpositron%e_pawdc
       end if
     end if
     if (optene==0) electronpositron%e0=etotal
     if (optene==1) electronpositron%e0=etotal-energies%edc_electronpositron
     etotal=electronpositron%e0+energies%e0_electronpositron+energies%e_electronpositron
   end if

!  Compute energy residual
   deltae=etotal-elast
   elast=etotal
!  DEBUG
!  write(std_out,*) 'eeig-ehart+enxc-enxcdc+eew+eii+eent+enefield+epawdc',eeig,ehart,enxc,enxcdc,eew,eii,eent,enefield,epawdc
!  ENDDEBUG

 end if !optene/=-1

 call timab(80,2,tsec)

!------Compute forces-----------------------------------------------------

 if (optforces==1) then

!  PAW: add gradients due to Dij derivatives to non-local term
   if (usepaw==1) then
     if (usexcnhat==0) then
       dimwork=1
       ABI_ALLOCATE(work,(nfft,dimwork))
!      $OMP PARALLEL DO PRIVATE(ifft) &
!      $OMP&SHARED(nfft,vhartr,work,vpsp)
       do ifft=1,nfft
         work(ifft,1)=vhartr(ifft)+vpsp(ifft)
       end do
!      $OMP END PARALLEL DO
     else
       dimwork=dtset%nspden
       ABI_ALLOCATE(work,(nfft,dimwork))
       do ispden=1,min(dtset%nspden,2)
!        $OMP PARALLEL DO PRIVATE(ifft) &
!        $OMP&SHARED(ispden,nfft,vhartr,work,vpsp,vxc)
         do ifft=1,nfft
           work(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)
         end do
!        $OMP END PARALLEL DO
       end do
       if(dtset%nspden==4)then
         do ispden=3,4
!          $OMP PARALLEL DO PRIVATE(ifft) &
!          $OMP&SHARED(ispden,nfft,work,vxc)
           do ifft=1,nfft
             work(ifft,ispden)=vxc(ifft,ispden)
           end do
!          $OMP END PARALLEL DO
         end do
       end if
     end if
     dimnhat=0
     optgr=1
     optgr2=0
     optstr=0
     ABI_ALLOCATE(nhat_dum,(1,0))
     call pawgrnl(atindx1,dimnhat,dimwork,dummy,1,grnl,gsqcut,mgfft,mpi_enreg,dtset%natom,dtset%natom,&
&     nattyp,nfft,ngfft,nhat_dum,dummy,dtset%nspden,dtset%nsym,ntypat,optgr,optgr2,optstr,dtset%paral_kgb,&
&     pawang,pawfgrtab,pawrhoij,pawtab,ph1d,psps,k0,rprimd,symrec,dtset%typat,work,xred)
     ABI_DEALLOCATE(nhat_dum)
     ABI_DEALLOCATE(work)
   end if

!  If residual is a density residual (and forces from residual asked),
!  has to convert it into a potential residualbefore calling forces routine
   if (optres==1 .and. dtset%usewvl==0.and.abs(dtset%iprcch)>=1 .and. &
&   abs(dtset%iprcch)<=6.and.abs(dtset%iprcch)/=5) then
     option=0; if (dtset%iprcch<0) option=1
     ABI_ALLOCATE(work,(nfft,dtset%nspden))
     optnc=1;if (dtset%nspden==4.and.(abs(dtset%iprcch)==4.or.abs(dtset%iprcch)==6)) optnc=2
     call nres2vres(dtset,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,nhat,&
&     nkxc,nvresid,n3xccc,optnc,option,pawang,pawfgrtab,pawrhoij,pawtab,&
&     rhor,rprimd,usepaw,usexcnhat,work,wvl,xccc3d,xred)
     call forces(atindx1,diffor,dtset,favg,fcart,forold,fred,gresid,grewtn,&
&     grhf,grnl,grxc,gsqcut,indsym,maxfor,mgfft,mpi_enreg,&
&     n1xccc,n3xccc,nattyp,nfft,ngfft,ntypat,&
&     pawtab,ph1d,psps,rhog,rhor,rprimd,symrec,synlgr,work,vxc,wvl,xred,&
&     electronpositron=electronpositron)
     ABI_DEALLOCATE(work)
   else
     call forces(atindx1,diffor,dtset,favg,fcart,forold,fred,gresid,grewtn,&
&     grhf,grnl,grxc,gsqcut,indsym,maxfor,mgfft,mpi_enreg,&
&     n1xccc,n3xccc,nattyp,nfft,ngfft,ntypat,&
&     pawtab,ph1d,psps,rhog,rhor,rprimd,symrec,synlgr,nvresid,vxc,wvl,xred,&
&     electronpositron=electronpositron)
   end if

!  Returned fred are full symmetrized gradients of Etotal
!  wrt reduced coordinates xred, d(Etotal)/d(xred)
!  Forces are contained in array fcart

 else   ! if optforces==0
   fcart=zero
   fred=zero
   favg=zero
   diffor=zero
   gresid=zero
   grhf=zero
   maxfor=zero
   synlgr=zero
 end if

 call timab(80,2,tsec)

!DEBUG
!write(std_out,*)' etotfor : exit '
!write(std_out,*) etotal
!call leave_new("COLL")
!stop
!ENDDEBUG

end subroutine etotfor
!!***
