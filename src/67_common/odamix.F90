!{\src2tex{textfont=tt}}
!!****f* ABINIT/odamix
!! NAME
!! odamix
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
!! For the initials of contributors, see ~abinit/doc/developers/contributors.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryopt  = 4: electric field is on -> add the contribution of the
!!   |                - \Omega E.P term to the total energy
!!   | berryopt  = 5: magnetic field is on -> add the contribution of the
!!   |                - \Omega B.M term to the total energy
!!   |          /= 4: electric field is off
!!   |          /= 5: magnetic field is off
!!   | bfield = cartesian coordinates of the magnetic field in atomic units
!!   | efield = cartesian coordinates of the electric field in atomic units
!!   | iatfix(3,natom)=1 for frozen atom along some direction, 0 for unfrozen
!!   | ionmov=governs the movement of atoms (see help file)
!!   | natom=number of atoms in cell.
!!   | nconeq=number of atomic constraint equations
!!   | nspden=number of spin-density components
!!   | nsym=number of symmetry elements in space group
!!   | occopt=option for occupancies
!!   | prtvol=integer controlling volume of printed output
!!   | tsmear=smearing energy or temperature (if metal)
!!   | typat(natom)=type integer for each atom in cell
!!   | wtatcon(3,natom,nconeq)=weights for atomic constraints
!!   | xclevel= XC functional level
!!  efield_dot = reciprocal lattice coordinates of the electric field
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  grewtn(3,natom)=grads of Ewald energy (hartree)
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  mag_cart = orbital magnetization in cartesian coordinates
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  ntypat=number of types of atoms in unit cell.
!!  nvresid(nfft,nspden)=potential or density residual
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  optres=0 if residual array (nvresid) contains the potential residual
!!        =1 if residual array (nvresid) contains the density residual
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  paw_an
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrad
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pel = reduced coordinates of the electronic polarization
!!        (pel does not take into account the factor 1/ucvol)
!!  pion = reduced coordinates of the ionic polarization
!!        (pel does not take into account the factor 1/ucvol)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=array for Fourier transform of electron density
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density 
!!  ucvol = unit cell volume (Bohr**3)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vhartr(nfft)=array for holding Hartree potential
!!  vpsp(nfft)=array for holding local psp
!!  vxc(nfft,nspden)=array for holding XC potential
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  deltae=change in total energy
!!         between the previous and present SCF cycle
!!  etotal=total energy (hartree)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  elast=previous value of the energy,
!!        needed to compute deltae, then updated.
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
!!   |                  on the magnetic field:  enmagfield = -ucvol*B*M
!!   | e_entropy(OUT)=entropy energy due to the occupation number smearing (if metal)
!!   |                this value is %entropy * dtset%tsmear (hartree).
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
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
!!      dotprod_vn,fourdp,leave_new,pawdenpot,pawmknhat,rhohxc,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine odamix(deltae,dtset,efield_dot,elast,energies,etotal,&
&          gprimd,gsqcut,kxc,mag_cart,mpi_enreg,nfft,ngfft,nhat,&
&          nkxc,ntypat,nvresid,n3xccc,optres,paw_ij,&
&          paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,pel,&
&          pion,psps,rhog,rhor,rprimd,strsxc,taug,taur,ucvol,usepaw,&
&          usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,xccc3d,xred,&
&          vxctau) ! optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_energies, only : energies_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'odamix'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_53_ffts
 use interfaces_53_spacepar
 use interfaces_56_xc
 use interfaces_66_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n3xccc,nfft,nkxc,ntypat,optres
 integer,intent(in) :: usepaw,usexcnhat
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(inout) :: elast
 real(dp),intent(out) :: deltae,etotal,vxcavg
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(energies_type),intent(inout) :: energies
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: efield_dot(3),gprimd(3,3),mag_cart(3),pel(3)
 real(dp),intent(in) :: pion(3),rprimd(3,3),vpsp(nfft),xred(3,dtset%natom)
 real(dp),intent(inout) :: kxc(nfft,nkxc),nhat(nfft,dtset%nspden*usepaw)
 real(dp),intent(inout) :: nvresid(nfft,dtset%nspden),rhog(2,nfft),taug(2,nfft*dtset%usekden)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden),taur(nfft,dtset%nspden*dtset%usekden),vhartr(nfft)
 real(dp),intent(inout) :: vtrial(nfft,dtset%nspden),vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc)
 real(dp),intent(out) :: strsxc(6)
 real(dp),intent(inout),optional :: vxctau(nfft,dtset%nspden*dtset%usekden,4)
 type(paw_an_type),intent(inout) :: paw_an(dtset%natom)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: cplex,iatom,ider,idir,ifft,ipert,irhoij,ispden,itypat,izero
 integer :: jrhoij,klmn,klmn1,kmix,nfftot,nhatgrdim,nselect,nzlmopt
 integer :: option,optxc
 integer :: nk3xc
 real(dp) :: alphaopt,compch_fft,compch_sph,doti,e1t10,e_ksnm1,e_xcdc_vxctau
 real(dp) :: emag,enonlocal,epel,epion,fp0,gammp1,ro_dlt,ucvol_local
 logical :: with_vxctau
 character(len=500) :: message
!arrays
 real(dp) :: qpt(3),tsec(2)
 real(dp),allocatable :: nhatgr(:,:,:),rhoijtmp(:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)' odamix : enter'
!ENDDEBUG

!Check that usekden is not 0 if want to use vxctau
 with_vxctau = .false.
 if (present(vxctau) .and. dtset%usekden /= 0) with_vxctau = .true.

!to be adjusted for the call to rhohxc
 nk3xc=1

 call timab(80,1,tsec)
!faire un test sur optres=1, usewvl=0, nspden=1,nhatgrdim
 if(optres/=1)then
   write(message,'(a1,a,a1,a,i1,a)') ch10,' scfcv : exit ',&
&   ch10,'  optres=',optres,', not allowed in oda => stop '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(dtset%usewvl/=0)then
   write(message,'(a1,a,a1,a,i1,a)') ch10,' scfcv : exit ',&
&   ch10,'  usewvl=',dtset%usewvl,', not allowed in oda => stop '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(dtset%nspden/=1)then
   write(message,'(a1,a,a1,a,i1,a)') ch10,' scfcv : exit ',&
&   ch10,'  nspden=',dtset%nspden,', not allowed in oda => stop '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(paw_ij(1)%has_dijhat==0)then
   write(message,'(a1,a,a1,a,i1,a)') ch10,' odamix: BUG -',ch10,&
&   '   dijhat variable must be allocated ! '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(paw_ij(1)%cplex==2)then
   write(message,'(a1,a,a1,a,i1,a)') ch10,' odamix: ERROR -',ch10,&
&   '   complex dij not allowed ! '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!calcul de f'(0)= Eband_new-EH_old-E_xcdc_old-Ek_old-E_loc_old-E_nonloc_old
!!!!!!!!!!! stockage de l'energie precedente E(rho_tild_n)

 fp0=energies%e_eigenvalues-energies%h0-two*energies%e_hartree-energies%e_xcdc
 if (usepaw==1) then
   do iatom=1,dtset%natom
     itypat=dtset%typat(iatom)
     do ispden=1,pawrhoij(iatom)%nspden
       jrhoij=1
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(irhoij)
         ro_dlt=pawrhoij(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
         e1t10=e1t10+ro_dlt*(paw_ij(iatom)%dij(klmn,ispden)-paw_ij(iatom)%dijhat(klmn,ispden))
         jrhoij=jrhoij+pawrhoij(iatom)%cplex
       end do
       klmn1=1
       do klmn=1,pawrhoij(iatom)%lmn2_size
         ro_dlt=-pawrhoij(iatom)%rhoijres(klmn1,ispden)*pawtab(itypat)%dltij(klmn)
         e1t10=e1t10+ro_dlt*(paw_ij(iatom)%dij(klmn,ispden)-paw_ij(iatom)%dijhat(klmn,ispden))
         klmn1=klmn1+pawrhoij(iatom)%cplex
       end do
     end do
     if (paw_ij(iatom)%ndij>=2.and.pawrhoij(iatom)%nspden==1) then
       jrhoij=1
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(irhoij)
         ro_dlt=pawrhoij(iatom)%rhoijp(jrhoij,1)*pawtab(itypat)%dltij(klmn)
         e1t10=e1t10+ro_dlt*(paw_ij(iatom)%dij(klmn,2)-paw_ij(iatom)%dijhat(klmn,2))
         jrhoij=jrhoij+pawrhoij(iatom)%cplex
       end do
       klmn1=1
       do klmn=1,pawrhoij(iatom)%lmn2_size
         ro_dlt=-pawrhoij(iatom)%rhoijres(klmn1,1)*pawtab(itypat)%dltij(klmn)
         e1t10=e1t10+ro_dlt*(paw_ij(iatom)%dij(klmn,2)-paw_ij(iatom)%dijhat(klmn,2))
         klmn1=klmn1+pawrhoij(iatom)%cplex
       end do
       e1t10=half*e1t10
     end if
   end do

   fp0=fp0-e1t10
 end if
 e_ksnm1=etotal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Calcul des quantités qui dépendent de rho_n+1

!PAW: eventually recompute compensation density (and gradients)
 nhatgrdim=0
 if (usepaw==1) then
   ider=-1;if (dtset%xclevel==2.or.usexcnhat==0) ider=0
   if (dtset%xclevel==2.and.usexcnhat==1) ider=ider+2
   if (ider>0) then
     nhatgrdim=1
     ABI_ALLOCATE(nhatgr,(nfft,dtset%nspden,3))
   end if
   if (ider>=0) then
     ider=0;izero=0;cplex=1;ipert=0;idir=0;qpt(:)=zero
     call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,mpi_enreg,dtset%natom,dtset%natom,&
     nfft,ngfft,nhatgrdim,dtset%nspden,ntypat,dtset%paral_kgb,pawang,pawfgrtab,&
&     nhatgr,nhat,pawrhoij,pawrhoij,pawtab,qpt,rprimd,ucvol,xred)
   end if
 end if

!------Compute Hartree and xc potentials----------------------------------
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
!Compute xc potential (separate up and down if spin-polarized)
 optxc=1
 if(with_vxctau)then
   call rhohxc(dtset,energies%e_xc,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,&
&   nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,optxc,rhog,rhor,rprimd,strsxc,&
&   usexcnhat,vhartr,vxc,vxcavg,xccc3d,taug=taug,taur=taur,vxctau=vxctau)
 else  
   call rhohxc(dtset,energies%e_xc,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,&
&   nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,optxc,rhog,rhor,rprimd,strsxc,&
&   usexcnhat,vhartr,vxc,vxcavg,xccc3d,taug=taug,taur=taur)
 end if

!------Compute parts of total energy depending on potentials--------

 ucvol_local = ucvol


!Compute Hartree energy energies%e_hartree
 call dotprod_vn(1,rhor,energies%e_hartree,doti,mpi_enreg,nfft,nfftot,1,1,vhartr,ucvol_local)
 energies%e_hartree=half*energies%e_hartree


!Compute local psp energy energies%e_localpsp
 call dotprod_vn(1,rhor,energies%e_localpsp,doti,mpi_enreg,nfft,nfftot,1,1,vpsp,ucvol_local)

!Compute double-counting XC energy energies%e_xcdc
 call dotprod_vn(1,rhor,energies%e_xcdc,doti,mpi_enreg,nfft,nfftot,dtset%nspden,1,vxc,ucvol_local)
 if (with_vxctau) then
   call dotprod_vn(1,taur,e_xcdc_vxctau,doti,mpi_enreg,nfft,nfftot,dtset%nspden,1,vxctau(:,:,1),ucvol_local)
   energies%e_xcdc=energies%e_xcdc+e_xcdc_vxctau
 end if

 if (usepaw/=0) then

   nzlmopt=dtset%pawnzlm; option=2
   do iatom=1,dtset%natom
     ABI_ALLOCATE(paw_ij(iatom)%dijhartree,(pawtab(dtset%typat(iatom))%lmn2_size))
     paw_ij(iatom)%has_dijhartree=1
   end do
   call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,0,dtset%ixc,mpi_enreg,dtset%natom,dtset%natom,dtset%nspden,ntypat,&
&   nzlmopt,option,dtset%paral_kgb,paw_an,paw_an,paw_ij,pawang,dtset%pawprtvol,pawrad,pawrhoij,dtset%pawspnorb,&
&   pawtab,dtset%pawxcdev,dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,psps%znuclpsp)
   do iatom=1,dtset%natom
     ABI_DEALLOCATE(paw_ij(iatom)%dijhartree)
     paw_ij(iatom)%has_dijhartree=0
   end do
 end if


!When the finite-temperature VG broadening scheme is used,
!the total entropy contribution "tsmear*entropy" has a meaning,
!and gather the two last terms of Eq.8 of the VG paper
!Warning : might have to be changed for fixed moment calculations
 if(dtset%occopt>=3 .and. dtset%occopt<=8) then
   energies%e_entropy = - dtset%tsmear * energies%entropy
 else
   energies%e_entropy = zero
 end if
!Turn it into an electric enthalpy, by adding both ionic and electronic contributions
 energies%e_elecfield = zero
 if (dtset%berryopt == 4) then
!  First, ionic contribution (epion):
   epion = -1_dp*(efield_dot(1)*pion(1) + efield_dot(2)*pion(2) + &
&   efield_dot(3)*pion(3))
   energies%e_elecfield = energies%e_elecfield + epion
!  Now, electronic contribution (epel):
   epel = -1_dp*(efield_dot(1)*pel(1) + efield_dot(2)*pel(2) + &
&   efield_dot(3)*pel(3))
   energies%e_elecfield = energies%e_elecfield + epel
 end if
 
 energies%e_magfield = zero
 if (dtset%berryopt == 5) then
   emag = dot_product(mag_cart,dtset%bfield)
   energies%e_magfield = emag
 end if

 etotal = energies%e_kinetic+ energies%e_hartree + energies%e_xc + &
& energies%e_localpsp + energies%e_corepsp + &
& energies%e_entropy + energies%e_elecfield + energies%e_magfield
!etotal = energies%e_eigenvalues - energies%e_hartree + energies%e_xc - &
!& energies%e_xcdc + energies%e_corepsp + &
!& energies%e_entropy + energies%e_elecfield
 etotal = etotal + energies%e_ewald
 if (usepaw==0) then
   etotal = etotal + energies%e_nonlocalpsp
 else
   etotal = etotal + energies%e_paw
 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! now, compute mixed densities



 gammp1=etotal-e_ksnm1-fp0
 if (fp0>0.d0) then
   write(std_out,*) "fp0 est positif"
!  stop
 end if
 write(std_out,*) "fp0 ",fp0
 alphaopt=-fp0/two/gammp1

 if (alphaopt>one.or.alphaopt<0.d0) alphaopt=one
 if (abs(energies%h0)<=tol10) alphaopt=one
 write(std_out,*) " alphaopt",alphaopt

 if (usepaw==0) then
   enonlocal =energies%e_nonlocalpsp
 else
   enonlocal = zero
 end if

 energies%h0=(one-alphaopt)*energies%h0+alphaopt*(energies%e_kinetic+energies%e_localpsp+enonlocal)
 rhor= rhor+(alphaopt-one)*nvresid
 call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)


 if (usepaw==1) then
   do iatom=1,dtset%natom

     ABI_ALLOCATE(rhoijtmp,(pawrhoij(iatom)%cplex*pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
     rhoijtmp=zero
     if (pawrhoij(iatom)%cplex==1) then
       if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
         do ispden=1,pawrhoij(iatom)%nspden
           do irhoij=1,pawrhoij(iatom)%nrhoijsel
             klmn=pawrhoij(iatom)%rhoijselect(irhoij)
             rhoijtmp(klmn,ispden)=pawrhoij(iatom)%rhoijp(irhoij,ispden)
           end do
         end do
       end if
       do ispden=1,pawrhoij(iatom)%nspden
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           klmn=pawrhoij(iatom)%kpawmix(kmix)
           rhoijtmp(klmn,ispden)=rhoijtmp(klmn,ispden)+(alphaopt-one)*pawrhoij(iatom)%rhoijres(klmn,ispden)
         end do
       end do
       nselect=0
       do klmn=1,pawrhoij(iatom)%lmn2_size
         if (any(abs(rhoijtmp(klmn,:))>tol10)) then
           nselect=nselect+1
           pawrhoij(iatom)%rhoijselect(nselect)=klmn
           do ispden=1,pawrhoij(iatom)%nspden
             pawrhoij(iatom)%rhoijp(nselect,ispden)=rhoijtmp(klmn,ispden)
           end do
         end if
       end do
       pawrhoij(iatom)%nrhoijsel=nselect
       ABI_DEALLOCATE(rhoijtmp)
     else
       if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
         jrhoij=1
         do ispden=1,pawrhoij(iatom)%nspden
           do irhoij=1,pawrhoij(iatom)%nrhoijsel
             klmn=2*pawrhoij(iatom)%rhoijselect(irhoij)-1
             rhoijtmp(klmn-1:klmn+1,ispden)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+1,ispden)
             jrhoij=jrhoij+2
           end do
         end do
       end if
       do ispden=1,pawrhoij(iatom)%nspden
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           klmn=2*pawrhoij(iatom)%kpawmix(kmix)-1
           rhoijtmp(klmn:klmn+1,ispden)=rhoijtmp(klmn:klmn+1,ispden) &
&           +(alphaopt-one)*pawrhoij(iatom)%rhoijres(klmn:klmn+1,ispden)
         end do
       end do
       nselect=0
       do klmn=1,pawrhoij(iatom)%lmn2_size
         if (any(abs(rhoijtmp(2*klmn-1:2*klmn,:))>tol10)) then
           nselect=nselect+1
           pawrhoij(iatom)%rhoijselect(nselect)=klmn
           do ispden=1,pawrhoij(iatom)%nspden
             pawrhoij(iatom)%rhoijp(2*nselect-1:2*nselect,ispden)=rhoijtmp(2*klmn-1:2*klmn,ispden)
           end do
         end if
       end do
       pawrhoij(iatom)%nrhoijsel=nselect
       ABI_DEALLOCATE(rhoijtmp)
     end if

   end do
 end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Calcul des quantites qui dependent de rho_tilde_n+1 (rho apres mixing)

 if (usepaw==1) then
   if (ider>=0) then
     izero=0;cplex=1;ipert=0;idir=0;qpt(:)=zero
     call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,mpi_enreg,dtset%natom,dtset%natom,&
&     nfft,ngfft,nhatgrdim,dtset%nspden,ntypat,dtset%paral_kgb,pawang,pawfgrtab,nhatgr,&
&     nhat,pawrhoij,pawrhoij,pawtab,qpt,rprimd,ucvol,xred)
   end if
 end if

!------Compute Hartree and xc potentials----------------------------------

!Compute xc potential (separate up and down if spin-polarized)
 optxc=1;if (nkxc>0) optxc=2
 if(with_vxctau)then
   call rhohxc(dtset,energies%e_xc,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,&
&   nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,optxc,rhog,rhor,rprimd,strsxc,&
&   usexcnhat,vhartr,vxc,vxcavg,xccc3d,taug=taug,taur=taur,vxctau=vxctau)
 else
   call rhohxc(dtset,energies%e_xc,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,&
&   nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,optxc,rhog,rhor,rprimd,strsxc,&
&   usexcnhat,vhartr,vxc,vxcavg,xccc3d,taug=taug,taur=taur)
 end if

 if (nhatgrdim>0)  then
   ABI_DEALLOCATE(nhatgr)
 end if
!------Compute parts of total energy depending on potentials--------

 ucvol_local = ucvol


!Compute Hartree energy energies%e_hartree
 call dotprod_vn(1,rhor,energies%e_hartree,doti,mpi_enreg,nfft,nfftot,1,1,vhartr,ucvol_local)
 energies%e_hartree=half*energies%e_hartree
!Compute double-counting XC energy energies%e_xcdc
 call dotprod_vn(1,rhor,energies%e_xcdc,doti,mpi_enreg,nfft,nfftot,dtset%nspden,1,vxc,ucvol_local)


 etotal=energies%h0+energies%e_hartree+energies%e_xc+energies%e_corepsp + &
& energies%e_entropy + energies%e_elecfield + energies%e_magfield
 etotal = etotal + energies%e_ewald
 if (usepaw==1) then
   do iatom=1,dtset%natom
     ABI_ALLOCATE(paw_ij(iatom)%dijhartree,(pawtab(dtset%typat(iatom))%lmn2_size))
     paw_ij(iatom)%has_dijhartree=1
   end do
   call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,0,dtset%ixc,mpi_enreg,dtset%natom,dtset%natom,dtset%nspden,ntypat,&
&   nzlmopt,option,dtset%paral_kgb,paw_an,paw_an,paw_ij,pawang,dtset%pawprtvol,pawrad,pawrhoij,dtset%pawspnorb,&
&   pawtab,dtset%pawxcdev,dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,psps%znuclpsp)
   do iatom=1,dtset%natom
     ABI_DEALLOCATE(paw_ij(iatom)%dijhartree)
     paw_ij(iatom)%has_dijhartree=0
   end do
   etotal=etotal+energies%e_paw
 end if
!Compute energy residual
 deltae=etotal-elast
 elast=etotal

 do ispden=1,min(dtset%nspden,2)
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(ispden,nfft,dtset%nspden,vhartr,vnew,vpsp,vxc)
   do ifft=1,nfft
     vtrial(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)
   end do
!  $OMP END PARALLEL DO
 end do
 if(dtset%nspden==4) vtrial(:,3:4)=vxc(:,3:4)

!write(std_out,*) 'coucou',energies%h0,energies%e_kinetic+energies%e_localpsp+enonlocal

!DEBUG
!write(std_out,*) 'eeig-ehart+enxc-enxcdc+eew+eii+eent+enefield+epawdc',eeig,ehart,enxc,enxcdc,eew,eii,eent,enefield,epawdc
!ENDEBUG
 call timab(80,2,tsec)


 call timab(80,2,tsec)

!DEBUG
!write(std_out,*)' odamix: exit '
!write(std_out,*) etotal
!call leave_new("COLL")
!stop
!ENDDEBUG

end subroutine odamix
!!***

