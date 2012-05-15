!{\src2tex{textfont=tt}}
!!****f* ABINIT/setvtr
!!
!! NAME
!! setvtr
!!
!! FUNCTION
!! Set up the trial potential and some energy terms
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG, GMR, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(dtset%natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | ikhxc=exchange-correlation kernel treatment parameter
!!   |       if =0,1 no xc kernel, =2 spin-averaged (LDA) kernel
!!   | iprcch=govern the choice of preconditioner for the SCF cycle
!!   | iscf=determines the way the SCF cycle is handled
!!   | natom=number of atoms in cell.
!!   | nspden=number of spin-density components
!!   | qprtrb(3)= integer wavevector of possible perturbing potential
!!   |            in basis of reciprocal lattice translations
!!   | typat(natom)=type integer for each atom in cell
!!   | vprtrb(2)=complex amplitude of possible perturbing potential; if nonzero,
!!   |  perturbing potential is added of the form
!!   |  V(G)=(vprtrb(1)+I*vprtrb(2))/2 at the values G=qprtrb and
!!   |  (vprtrb(1)-I*vprtrb(2))/2 at G=-qprtrb (integers)
!!   |  for each type of atom, from psp (used in norm-conserving only)
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2) (sphere for density and potential)
!!  istep=step number in the main loop of scfcv
!!  mgfft=maximum size of 1D FFTs
!!  moved_rhor=1 if the density was moved just before
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!!  nhatgr(nfft,nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  nkxc=second dimension of the array kxc
!!  ntypat=number of types of atoms in unit cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  optene=>0 if some additional energies have to be computed
!!  pawtab(ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase (structure factor) information.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=Fourier transform of electron density
!!  rhor(nfft,nspden)=electron density in electrons/bohr**3.
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  ucvol = unit cell volume (bohr^3)
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  grewtn(3,natom)=grads of Ewald energy (hartree)
!!  kxc(nfft,nkxc)=exchange-correlation kernel, will be computed if nkxc/=0 .
!!                 see routine rhohxc for a more complete description
!!  strsxc(6)=xc contribution to stress tensor (hartree/bohr^3)
!!  vxcavg=mean of the vxc potential
!!  energies <type(energies_type)>=all part of total energy.
!!   | e_xc=exchange-correlation energy (hartree)
!!  ==== if optene==2 or 4
!!   | e_localpsp=local psp energy (hartree)
!!  ==== if dtset%icoulomb == 0
!!   | e_ewald=Ewald energy (hartree)
!!  ==== if optene>=1
!!   | e_hartree=Hartree part of total energy (hartree units)
!!  ==== if optene==3 or 4
!!   | e_xcdc=exchange-correlation double-counting energy (hartree)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  moved_atm_inside=1 if the atomic positions were moved inside the SCF loop.
!!  vhartr(nfft)=Hartree potential (Hartree)
!!  vpsp(nfft)=local psp (Hartree)
!!  vtrial(nfft,nspden)= trial potential (Hartree)
!!  vxc(nfft,nspden)= xc potential (Hartree)
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
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
!!      bethe_salpeter,scfcv,screening,sigma
!!
!! CHILDREN
!!      atm2fft,dotprod_vn,ewald,ionion_realspace,jellium,mkcore,mklocl
!!      psolver_rhohxc,rhohxc,rhohxcpositron,timab,wrtout,wvl_newvtr
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setvtr(atindx1,dtset,energies,gmet,gprimd,grewtn,gsqcut,&
&  istep,kxc,mgfft,moved_atm_inside,moved_rhor,mpi_enreg,&
&  nattyp,nfft,ngfft,nhat,nhatgr,nhatgrdim,nkxc,ntypat,n1xccc,n3xccc,&
&  optene,pawtab,ph1d,psps,rhog,rhor,rmet,rprimd,strsxc,&
&  ucvol,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,wvl,xccc3d,xred,&
&  electronpositron,taug,taur,vxctau) ! optional argument

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
#define ABI_FUNC 'setvtr'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_53_spacepar
 use interfaces_56_xc
 use interfaces_62_poisson
 use interfaces_65_psp
 use interfaces_67_common, except_this_one => setvtr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mgfft,n1xccc,n3xccc,nfft,nhatgrdim,nkxc,ntypat
 integer,intent(in) :: optene,usexcnhat
 integer,intent(inout) :: moved_atm_inside,moved_rhor
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(out) :: vxcavg
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer,optional :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 integer,intent(in) :: atindx1(dtset%natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: nhat(nfft,dtset%nspden*psps%usepaw)
 real(dp),intent(in) :: nhatgr(nfft,dtset%nspden,3*nhatgrdim),rhog(2,nfft)
 real(dp),intent(in),optional :: taug(2,nfft*dtset%usekden)
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3)
 real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden),vhartr(nfft),vpsp(nfft)
 real(dp),intent(inout),optional :: taur(nfft,dtset%nspden*dtset%usekden)
 real(dp),intent(inout) :: vtrial(nfft,dtset%nspden),vxc(nfft,dtset%nspden)
 real(dp),intent(out),optional :: vxctau(nfft,dtset%nspden*dtset%usekden,4)
 real(dp),intent(inout) :: xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp),intent(out) :: grewtn(3,dtset%natom),kxc(nfft,nkxc),strsxc(6)
 type(pawtab_type),intent(in) :: pawtab(ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: nk3xc
 integer :: iatom,ifft,ipositron,ispden,n1,n2,n3,nfftot,nn,offset
!integer :: jj,kk
 integer :: optatm,optdyfr,optgr,option,optn,optn2,optstr,optv
 real(dp) :: doti,e_xcdc_vxctau,ebb,ebn,ucvol_local
!real(dp) :: sum
 logical :: with_vxctau
 character(len=500) :: message
!arrays
 real(dp),parameter :: identity(1:4)=(/1._dp,1._dp,0._dp,0._dp/)
 real(dp),parameter :: qphon(3)=(/0._dp,0._dp,0._dp/)
 real(dp) :: dummy6(6),rhodum(1),tsec(2)
 real(dp),allocatable :: dummy_in(:),dummy_out(:),dyfr_dum(:,:,:),gr_dum(:,:),grtn(:,:)
 real(dp),allocatable :: rhojellg(:,:),rhojellr(:),rhowk(:,:),vjell(:)

! *********************************************************************

!DEBUG
!write(std_out,*)' setvtr : enter '
!sum=0.0
!do jj=1,dtset%nspden
!do kk=1,nfft
!sum=sum+abs(rhor(kk,jj))
!end do
!end do
!write(std_out,*) 'SUM(rhor)',sum
!write(std_out,*)' n1xccc=',n1xccc
!write(std_out,*)' moved_atm_inside=',moved_atm_inside
!write(std_out,*)' istep=',istep
!write(std_out,*)' moved_rhor=',moved_rhor
!stop
!ENDDEBUG

!Check that usekden is not 0 if want to use vxctau
 with_vxctau = .false.
 if (present(vxctau) .and. present(taur) .and. dtset%usekden /= 0) with_vxctau = .true.

 call timab(91,1,tsec)

!Get size of FFT grid
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)

!Test electron-positron case
 ipositron=0;if (present(electronpositron)) ipositron=electronpositron_calctype(electronpositron)

!Get Ewald energy and Ewald forces
!--------------------------------------------------------------
 call timab(5,1,tsec)
 if (ipositron/=1) then
   if ((dtset%icoulomb == 0) .or. (dtset%icoulomb == 2)) then
!    Periodic system, need to compute energy and forces due to replica and
!    to correct the shift in potential calculation.
     call ewald(energies%e_ewald,gmet,grewtn,dtset%natom,ntypat,rmet,dtset%typat,ucvol,xred,psps%ziontypat)
   else if (dtset%icoulomb == 1) then
!    In a non periodic system (real space computation), the G=0 divergence
!    doesn't occur and ewald is not needed. Only the ion/ion interaction
!    energy is relevant and used as Ewald energy and gradient.
     call ionion_realSpace(dtset, energies%e_ewald, grewtn, rprimd, xred, psps%ziontypat)
   end if
 else
   energies%e_ewald=zero
   grewtn=zero
 end if
 call timab(5,2,tsec)

!PAW: compute Vloc and core charge together in reciprocal space
!--------------------------------------------------------------
 if (psps%usepaw==1) then

   call timab(552,1,tsec)

   optatm=1;optdyfr=0;optgr=0;optstr=0;optv=1;optn=n3xccc/nfft;optn2=1
   call atm2fft(atindx1,xccc3d,vpsp,dummy_out,dummy_out,energies%e_localpsp,dummy_in,gmet,gprimd,&
&   dummy_out,dummy_out,gsqcut,mgfft,mpi_enreg,psps%mqgrid_vl,&
&   dtset%natom,nattyp,nfft,ngfft,ntypat,&
&   optatm,optdyfr,optgr,optn,optn2,optstr,optv,dtset%paral_kgb,&
&   pawtab,ph1d,psps%qgrid_vl,dtset%qprtrb,dummy_in,dummy6,dummy6,&
&   ucvol,psps%usepaw,dummy_in,dtset%vprtrb,psps%vlspl)
   call timab(552,2,tsec)

 else

!  Norm-conserving: compute Vloc in reciprocal space
!  and core charge in real space
!  --------------------------------------------------------------

!  Compute local ionic pseudopotential vpsp
   option=1
   ABI_ALLOCATE(dyfr_dum,(3,3,dtset%natom))
   ABI_ALLOCATE(gr_dum,(3,dtset%natom))
   call mklocl(dtset,dyfr_dum,energies%e_localpsp,gmet,gprimd,&
&   gr_dum,gsqcut,dummy6,mgfft,mpi_enreg,dtset%natom,nattyp,&
&   nfft,ngfft,dtset%nspden,ntypat,option,ph1d,psps,&
&   dtset%qprtrb,rhodum,rhor,rprimd,ucvol,dtset%vprtrb,vpsp,wvl,xred)

!  Compute 3D core electron density xccc3d
   if (n1xccc/=0) then
     call timab(91,2,tsec)
     call timab(92,1,tsec)
     n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
     call mkcore(dummy6,dyfr_dum,gr_dum,mpi_enreg,dtset%natom,nfft,dtset%nspden,ntypat,&
&     n1,n1xccc,n2,n3,option,rprimd,dtset%typat,ucvol,vxc,psps%xcccrc,psps%xccc1d,xccc3d,xred)
     call timab(92,2,tsec)
     call timab(91,1,tsec)
   end if
   ABI_DEALLOCATE(dyfr_dum)
   ABI_DEALLOCATE(gr_dum)

 end if  ! PAW or NC

!Adds the jellium potential to the local part of ionic potential
 if (dtset%jellslab/=0) then
   ABI_ALLOCATE(vjell,(nfft))
   ABI_ALLOCATE(rhojellg,(2,nfft))
   ABI_ALLOCATE(rhojellr,(nfft))
   option=1
   call jellium(gmet,gsqcut,mpi_enreg,nfft,ngfft,dtset%nspden,option,dtset%paral_kgb,&
&   dtset%slabwsrad,rhojellg,rhojellr,rprimd,vjell,dtset%slabzbeg,dtset%slabzend)
!  Compute background-background energy
   call dotprod_vn(1,rhojellr,ebb,doti,mpi_enreg,nfft,nfftot,1,1,vjell,ucvol)
   ebb=half*ebb
!  Compute electrostatic energy between background and nuclei before adding vjell to vpsp
   call dotprod_vn(1,rhojellr,ebn,doti,mpi_enreg,nfft,nfftot,1,1,vpsp,ucvol)
!  Update e_ewald with ebb and ebn
   energies%e_ewald=energies%e_ewald+ebb+ebn
!  Compute gradient of ebn wrt tn
   if (psps%usepaw==1) then
     write(message, '(a,a,a,a,a,a)' )ch10,&
&     ' setvtr : WARNING -',ch10,&
&     '  The computation of forces due to jellium background',ch10,&
&     '  has to be verified in the PAW formalism.'
     call wrtout(std_out,message,'COLL')
     ABI_ALLOCATE(grtn,(3,dtset%natom))
     optatm=0;optdyfr=0;optgr=1;optstr=0;optv=1;optn=0;optn2=1
     call atm2fft(atindx1,dummy_out,vpsp,dummy_out,dummy_out,energies%e_localpsp,dummy_in,gmet,gprimd,&
&     dummy_out,grtn,gsqcut,mgfft,mpi_enreg,psps%mqgrid_vl,&
&     dtset%natom,nattyp,nfft,ngfft,ntypat,&
&     optatm,optdyfr,optgr,optn,optn2,optstr,optv,dtset%paral_kgb,&
&     pawtab,ph1d,psps%qgrid_vl,dtset%qprtrb,rhojellg,dummy6,dummy6,&
&     ucvol,psps%usepaw,dummy_in,dtset%vprtrb,psps%vlspl)
!    Update grewtn with gradient of ebn wrt tn
     do iatom=1,dtset%natom
       grewtn(1:3,iatom)=grewtn(1:3,iatom)+grtn(1:3,iatom)
     end do
     ABI_DEALLOCATE(grtn)
   else ! of usepaw==1
     option=2
     ABI_ALLOCATE(dyfr_dum,(3,3,dtset%natom))
     ABI_ALLOCATE(grtn,(3,dtset%natom))
     call mklocl(dtset,dyfr_dum,energies%e_localpsp,gmet,gprimd,&
&     grtn,gsqcut,dummy6,mgfft,mpi_enreg,dtset%natom,nattyp,&
&     nfft,ngfft,1,ntypat,option,ph1d,psps,dtset%qprtrb,rhojellg,&
&     rhojellr,rprimd,ucvol,dtset%vprtrb,vpsp,wvl,xred)
!    Update grewtn with gradient of ebn wrt tn (reestablish order of atoms)
     do iatom=1,dtset%natom
       grewtn(1:3,atindx1(iatom))=grewtn(1:3,atindx1(iatom))+grtn(1:3,iatom)
     end do
     ABI_DEALLOCATE(dyfr_dum)
     ABI_DEALLOCATE(grtn)
   end if ! of usepaw==1
   vpsp(:)=vpsp(:)+vjell(:)
   ABI_DEALLOCATE(vjell)
   ABI_DEALLOCATE(rhojellg)
   ABI_DEALLOCATE(rhojellr)
 end if

!Additional stuff for electron-positron calculation
!Compute the electronic/positronic local (Hartree) potential
 if (ipositron==1) vpsp=-vpsp

!If we are at the initialisation, or
!if the atom positions has changed and the non-linear core correction
!is included, or the rhor has changed, one needs to compute the xc stuff.
!One needs also to compute the Hartree stuff if the density changed,
!or at initialisation.
!--------------------------------------------------------------

!DEBUG
 write(std_out,*)' setvtr : istep,n1xccc,moved_rhor=',istep,n1xccc,moved_rhor
!ENDDEBUG

 if(istep==1 .or. n1xccc/=0 .or. moved_rhor==1 .or. dtset%positron<0) then

   option=0
   if(istep==1 .or. moved_rhor==1 .or. dtset%positron<0) option=1
   if (nkxc>0) option=2
   if (dtset%xclevel==2.and.(nkxc==3-2*mod(dtset%nspden,2))) option=12
   if(dtset%iscf==-1) option=-2
   if(dtset%ikhxc==2) then
     write(std_out,*)' %setvtr: CALL rhohxc with option=2'         !MF!DEBUGLINE
     write(std_out,*)' %        computing kxc = (kxc++ + kxc+-)/2' !MF!DEBUGLINE
   end if

   if (ipositron/=1) then
     if (dtset%icoulomb == 0) then
!      Use the periodic solver to compute Hxc
       nk3xc=1
       if (ipositron==0) then
         if(with_vxctau)then
           call rhohxc(dtset,energies%e_xc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfft,ngfft,&
&           nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,&
&           option,rhog,rhor,rprimd,strsxc,&
&           usexcnhat,vhartr,vxc,vxcavg,xccc3d,taug=taug,taur=taur,vxctau=vxctau)
         else
           call rhohxc(dtset,energies%e_xc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfft,ngfft,&
&           nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,&
&           option,rhog,rhor,rprimd,strsxc,&
&           usexcnhat,vhartr,vxc,vxcavg,xccc3d,taug=taug,taur=taur)
         end if
       else if (ipositron==2) then
         if(with_vxctau)then
           call rhohxc(dtset,energies%e_xc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfft,ngfft,&
&           nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,&
&           option,rhog,rhor,rprimd,strsxc,&
&           usexcnhat,vhartr,vxc,vxcavg,xccc3d,electronpositron=electronpositron,taug=taug,taur=taur,vxctau=vxctau)
         else
           call rhohxc(dtset,energies%e_xc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfft,ngfft,&
&           nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,&
&           option,rhog,rhor,rprimd,strsxc,&
&           usexcnhat,vhartr,vxc,vxcavg,xccc3d,electronpositron=electronpositron,taug=taug,taur=taur)
         end if
       end if
     else
!      Use the free boundary solver
       call PSolver_rhohxc(dtset, energies%e_hartree, energies%e_xc, energies%e_vxc, &
&       mpi_enreg, rhor, rprimd, vhartr, vxc, vxcavg, wvl)
     end if
   else
     energies%e_xc=zero
     call rhohxcpositron(electronpositron,gprimd,kxc,mpi_enreg,nfft,ngfft,nhat,nkxc,dtset%nspden,n3xccc,&
&     dtset%paral_kgb,rhor,strsxc,ucvol,usexcnhat,psps%usepaw,vhartr,vxc,vxcavg,xccc3d,dtset%xc_denpos)
   end if
   if (ipositron/=0) then
     if (optene>=1) then
       call dotprod_vn(1,rhor,electronpositron%e_hartree,doti,mpi_enreg,&
&       nfft,nfftot,1,1,electronpositron%vha_ep,ucvol)
     end if
     vhartr=vhartr+electronpositron%vha_ep
   end if
 end if

!Compute the trial potential
!-------------------------------------------------------------
 if (dtset%usewvl == 0) then
!  Now, compute trial Hxc potential. Local psp potential will be added later.
   if(moved_atm_inside==0 .or.dtset%iscf>=10) then

!    Compute starting Hxc potential.
!    Multiply by identity, should not change anything if nspden /= 4
     do ispden=1,dtset%nspden
       vtrial(:,ispden)=vhartr(:)*identity(ispden)+vxc(:,ispden)
     end do

   else

!    One should be here only when moved_atm_inside==1
!    The (H)xc now added corrects the previous one.
     if(dtset%iprcch==1)then
!      xc was substracted off. This should be rationalized later
       do ispden=1,dtset%nspden
         vtrial(:,ispden)=vtrial(:,ispden)+vxc(:,ispden)
       end do
     else if(abs(dtset%iprcch)==2 .or.abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6)then
!      Hxc was substracted off. This should be rationalized later
       do ispden=1,dtset%nspden
         vtrial(:,ispden)=vtrial(:,ispden)+vhartr(:)*identity(ispden)+vxc(:,ispden)
       end do
     end if
   end if

!  Adds the local part of the potential
   if ((moved_atm_inside==0).or.(dtset%iprcch/=3)) then
     do ispden=1,min(2,dtset%nspden)
       do ifft=1,nfft
         vtrial(ifft,ispden)=vtrial(ifft,ispden)+vpsp(ifft)
       end do
     end do
   end if

 else

   call wvl_newvtr(dtset, mpi_enreg, nn, offset, vhartr, vpsp, vtrial, vxc, wvl)

 end if

!Add the zeeman field to vtrial
!Noncollinear:
 if(dtset%nspden==4)then
   do ispden=2,dtset%nspden
     do ifft=1,nfft
       vtrial(ifft,ispden) = vtrial(ifft,ispden) + dtset%zeemanfield(ispden-1)
     end do !ifft
   end do !ispden
 else if(dtset%nspden==2)then
   do ifft=1,nfft
     vtrial(ifft,2) = vtrial(ifft,2) + dtset%zeemanfield(3) ! ispden=2 is rho_up only
   end do !ifft
 end if

!Compute parts of total energy depending on potentials
!--------------------------------------------------------------
 if (dtset%usewvl == 0) then
   ucvol_local = ucvol
 else
!  We need to tune the volume when wavelets are used because, not
!  all FFT points are used.
   ucvol_local = (half * dtset%wvl_hgrid) ** 3 * nfftot
 end if

!DEBUG
!write(std_out,*)' setvtr : will compute ehart, optene,dtset%icoulomb=',optene,dtset%icoulomb
!ENDDEBUG

 if (optene>=1 .and. dtset%icoulomb == 0) then

!  Compute Hartree energy ehart
!  Already available in the Psolver case through psolver_rhohxc().
   if (ipositron/=1) then
     call dotprod_vn(1,rhor,energies%e_hartree,doti,mpi_enreg,nfft,nfftot,1,1,vhartr,ucvol_local)
     if (ipositron==0) energies%e_hartree = half * energies%e_hartree
     if (ipositron==2) energies%e_hartree = half *(energies%e_hartree-electronpositron%e_hartree)
   else
     energies%e_hartree=zero
   end if

 end if

 if (optene==2.or.optene==4) then

!  Compute local psp energy eei
   call dotprod_vn(1,rhor,energies%e_localpsp,doti,mpi_enreg,nfft,nfftot,1,1,vpsp,ucvol_local)

 end if

 if (optene==3.or.optene==4) then

!  Compute double-counting XC energy enxcdc
   if (ipositron/=1) then
     if (dtset%usepaw==0.or.usexcnhat/=0) then
       call dotprod_vn(1,rhor,energies%e_xcdc,doti,mpi_enreg,nfft,nfftot,dtset%nspden,1,vxc,ucvol_local)
       if(with_vxctau)then
         call dotprod_vn(1,taur,e_xcdc_vxctau,doti,mpi_enreg,nfft,nfftot,dtset%nspden,1,vxctau(:,:,1),ucvol_local)
         energies%e_xcdc=energies%e_xcdc+e_xcdc_vxctau
       end if
     else
       ABI_ALLOCATE(rhowk,(nfft,dtset%nspden))
       rhowk=rhor-nhat
       call dotprod_vn(1,rhowk,energies%e_xcdc,doti,mpi_enreg,nfft,nfftot,dtset%nspden,1,vxc,ucvol_local)
       ABI_DEALLOCATE(rhowk)
     end if
     if (ipositron==2) energies%e_xcdc=energies%e_xcdc-electronpositron%e_xcdc
   else
     energies%e_xcdc=zero
   end if

 end if

!--------------------------------------------------------------

!The initialisation for the new atomic positions has been done
 moved_atm_inside=0

 call timab(91,2,tsec)

!DEBUG
!write(std_out,*)' setvtr : exit '
!stop
!ENDDEBUG

end subroutine setvtr
!!***
