!{\src2tex{textfont=tt}}
!!****f* ABINIT/forces
!!
!! NAME
!! forces
!!
!! FUNCTION
!! Assemble gradients of various total energy terms with respect
!! to reduced coordinates, including possible symmetrization,
!! in order to produce forces.
!!     fcart(i,iat) = d(Etot)/(d(r(i,iat)))
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryopt  = 4: electric field is on -> add the contribution of the
!!   |                               - \Omega E.P term to the total energy
!!   |          /= 4: electric field is off
!!   | efield = cartesian coordinates of the electric field in atomic units
!!   | iatfix(3,natom)=1 for frozen atom along specified direction, 0 for unfrozen
!!   | ionmov=governs the movement of atoms (see help file)
!!   | iprcch=governs the mixed electronic-atomic part of the preconditioner
!!   | natom=number of atoms in cell
!!   | nconeq=number of atomic constraint equations
!!   | nspden=number of spin-density components
!!   | nsym=number of symmetries in space group
!!   | prtvol=integer controlling volume of printed output
!!   | typat(natom)=type integer for each atom in cell
!!   | wtatcon(3,natom,nconeq)=weights for atomic constraints
!!  grewtn(3,natom)=d(Ewald)/d(xred) (hartree)
!!  grnl(3*natom)=gradients of Etot due to nonlocal contributions
!!  gsqcut=cutoff value on G**2 for (large) sphere inside FFT box.
!!                       gsqcut=(boxcut**2)*ecut/(2._dp*(Pi**2)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms
!!  pawtab(ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase (structure factor) array
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=Fourier transform of charge density (bohr^-3)
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  vresid(nfft,nspden)=potential residual (if non-collinear magn., only trace of it)
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree) in real space
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)=previous reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  diffor=maximal absolute value of changes in the components of
!!         force between the input and the output.
!!  favg(3)=mean of the forces before correction for translational symmetry
!!  forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!  fred(3,natom)=symmetrized grtn = d(etotal)/d(xred)
!!  gresid(3,natom)=forces due to the residual of the density/potential
!!  grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!  grxc(9+3*natom)=d(Exc)/d(xred) if core charges are used
!!  maxfor=maximal absolute value of the output array force.
!!  synlgr(3,natom)=symmetrized d(enl)/d(xred)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  fcart(3,natom)=forces in cartesian coordinates (Ha/Bohr)
!!    Note : unlike fred, this array has been corrected by enforcing
!!    the translational symmetry, namely that the sum of force
!!    on all atoms is zero.
!!
!! NOTES
!! * Symmetrization of gradients with respect to reduced
!!   coordinates xred is conducted according to the expression
!!   [d(e)/d(t(n,a))]_symmetrized = (1/Nsym) Sum(S) symrec(n,m,S)*
!!                [d(e)/d(t(m,b))]_unsymmetrized
!!   where t(m,b)= (symrel^-1)(m,n)*(t(n,a)-tnons(n)) and tnons
!!   is a possible nonsymmorphic translation.  The label "b" here
!!   refers to the atom which gets rotated into "a" under symmetry "S".
!!   symrel is the symmetry matrix in real space, which is the inverse
!!   transpose of symrec.  symrec is the symmetry matrix in reciprocal
!!   space.  sym_cartesian = R * symrel * R^-1 = G * symrec * G^-1
!!   where the columns of R and G are the dimensional primitive translations
!!   in real and reciprocal space respectively.
!! * Note the use of "symrec" in the symmetrization expression above.
!!
!! PARENTS
!!      afterscfloop,etotfor,forstr
!!
!! CHILDREN
!!      atm2fft,constrf,fourdp,fred2fcart,fresid,fresidrsp,metric,mkcore,mklocl
!!      sygrad,timab,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine forces(atindx1,diffor,dtset,favg,fcart,forold,fred,gresid,grewtn,&
&                  grhf,grnl,grxc,gsqcut,indsym,&
&                  maxfor,mgfft,mpi_enreg,n1xccc,n3xccc,&
&                  nattyp,nfft,ngfft,ntypat,&
&                  pawtab,ph1d,psps,rhog,rhor,rprimd,symrec,synlgr,&
&                  vresid,vxc,wvl,xred,&
&                  electronpositron) ! optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'forces'
 use interfaces_18_timing
 use interfaces_42_geometry
 use interfaces_45_geomoptim
 use interfaces_53_ffts
 use interfaces_56_xc
 use interfaces_65_psp
 use interfaces_67_common, except_this_one => forces
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,n1xccc,n3xccc,nfft,ntypat
 real(dp),intent(in) :: gsqcut
 real(dp),intent(out) :: diffor,maxfor
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer,optional :: electronpositron
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 integer,intent(in) :: atindx1(dtset%natom),indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: grewtn(3,dtset%natom),grnl(3*dtset%natom)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: fcart(3,dtset%natom),forold(3,dtset%natom)
 real(dp),intent(inout) :: vresid(nfft,dtset%nspden),xred(3,dtset%natom)
 real(dp),intent(out) :: favg(3),fred(3,dtset%natom),gresid(3,dtset%natom)
 real(dp),intent(out) :: grhf(3,dtset%natom),grxc(3,dtset%natom)
 real(dp),intent(out) :: synlgr(3,dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: iatom,indx,ipositron,itypat,mu,optatm,optdyfr,optgr,option,optn
 integer :: optn2,optstr,optv
 real(dp) :: eei_dum,ucvol
!arrays
 integer :: qprtrb_dum(3)
 real(dp) :: dummy6(6),fioncart(3),gmet(3,3),gprimd(3,3),rmet(3,3),tsec(2)
 real(dp) :: vprtrb_dum(2)
 real(dp),allocatable :: dummy_in(:),dummy_out(:),dyfrlo_dum(:,:,:),dyfrx2_dum(:,:,:),fin(:,:)
 real(dp),allocatable :: fionred(:,:),grl(:,:),grnl_tmp(:,:),grtn(:,:)
 real(dp),allocatable :: grtn_indx(:,:),v_dum(:),vxctotg(:,:)
 real(dp),allocatable :: xccc3d_dum(:)

! *************************************************************************

!DEBUG
!write(std_out,*)' forces: enter '
!stop
!ENDDEBUG

 call timab(69,1,tsec)

!Save input value of forces
 ABI_ALLOCATE(fin,(3,dtset%natom))
 fin(:,:)=fcart(:,:)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!=======================================================================
!========= Local pseudopotential and core charge contributions =========
!=======================================================================

 ABI_ALLOCATE(grl,(3,dtset%natom))

!PAW: compute local psp and core charge contribs together
!in reciprocal space
!-----------------------------------------------------------------------
 if (psps%usepaw==1) then

   call timab(550,1,tsec)
   if (n3xccc>0) then
     ABI_ALLOCATE(v_dum,(nfft))
     ABI_ALLOCATE(vxctotg,(2,nfft))
     v_dum(:)=vxc(:,1);if (dtset%nspden>=2) v_dum(:)=0.5_dp*(v_dum(:)+vxc(:,2))
     call fourdp(1,vxctotg,v_dum,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
     call zerosym(vxctotg,2,mpi_enreg,ngfft(1),ngfft(2),ngfft(3))
     ABI_DEALLOCATE(v_dum)
   end if
   optatm=0;optdyfr=0;optgr=1;optstr=0;optv=1;optn=n3xccc/nfft;optn2=1
   call atm2fft(atindx1,dummy_out,dummy_out,dummy_out,dummy_out,eei_dum,dummy_in,gmet,gprimd,&
&   grxc,grl,gsqcut,mgfft,mpi_enreg,psps%mqgrid_vl,&
&   dtset%natom,nattyp,nfft,ngfft,ntypat,&
&   optatm,optdyfr,optgr,optn,optn2,optstr,optv,dtset%paral_kgb,&
&   pawtab,ph1d,psps%qgrid_vl,qprtrb_dum,rhog,dummy6,dummy6,&
&   ucvol,psps%usepaw,vxctotg,vprtrb_dum,psps%vlspl)
   if (n3xccc>0)  then
     ABI_DEALLOCATE(vxctotg)
   end if
   if (n3xccc==0) grxc=zero
   call timab(550,2,tsec)
 else

!  Norm-conserving: compute local psp contribution in reciprocal space
!  and core charge contribution in real space
!  -----------------------------------------------------------------------
   option=2
   ABI_ALLOCATE(dyfrlo_dum,(3,3,dtset%natom))
   ABI_ALLOCATE(grtn_indx,(3,dtset%natom))
   ABI_ALLOCATE(v_dum,(nfft))
   call mklocl(dtset,dyfrlo_dum,eei_dum,gmet,gprimd,grtn_indx,gsqcut,dummy6,mgfft,&
&   mpi_enreg,dtset%natom,nattyp,nfft,ngfft,dtset%nspden,ntypat,option,ph1d,psps,&
&   qprtrb_dum,rhog,rhor,rprimd,ucvol,vprtrb_dum,v_dum,wvl,xred)

   do iatom=1,dtset%natom
!    Has to use the indexing array atindx1
     grl(1:3,atindx1(iatom))=grtn_indx(1:3,iatom)
   end do
   ABI_DEALLOCATE(dyfrlo_dum)
   ABI_DEALLOCATE(grtn_indx)
   ABI_DEALLOCATE(v_dum)
!  If gradients are computed in real space, we need to symetrise
!  the system before summing.
!  Rshaltaf: I changed the following line to include surfaces BC
   if (dtset%icoulomb == 1 .or. dtset%icoulomb == 2) then
     ABI_ALLOCATE(grnl_tmp,(3,dtset%natom))
     call sygrad(grnl_tmp,dtset%natom,grl,dtset%nsym,symrec,indsym)
     grl(:, :) = grnl_tmp(:, :)
     ABI_DEALLOCATE(grnl_tmp)
   end if

   if (n3xccc>0) then
     call timab(53,1,tsec)
     ABI_ALLOCATE(dyfrx2_dum,(3,3,dtset%natom))
     ABI_ALLOCATE(xccc3d_dum,(n3xccc))
     call mkcore(dummy6,dyfrx2_dum,grxc,mpi_enreg,dtset%natom,nfft,dtset%nspden,ntypat,ngfft(1),n1xccc,ngfft(2),&
&     ngfft(3),option,rprimd,dtset%typat,ucvol,vxc,psps%xcccrc,psps%xccc1d,xccc3d_dum,xred)
     ABI_DEALLOCATE(dyfrx2_dum)
     ABI_DEALLOCATE(xccc3d_dum)
     call timab(53,2,tsec)
   else
     grxc(:,:)=zero
   end if
 end if

!=======================================================================
!===================== Nonlocal contributions ==========================
!=======================================================================

!Only has to apply symmetries
 ABI_ALLOCATE(grnl_tmp,(3,dtset%natom))
 do iatom=1,dtset%natom
   indx=3*(iatom-1);grnl_tmp(1:3,atindx1(iatom))=grnl(indx+1:indx+3)
 end do
 if (dtset%usewvl == 0) then
   call sygrad(synlgr,dtset%natom,grnl_tmp,dtset%nsym,symrec,indsym)
 else
   synlgr = grnl_tmp
 end if
 ABI_DEALLOCATE(grnl_tmp)

!=======================================================================
!============ Density/potential residual contributions =================
!=======================================================================

 if (dtset%usewvl==0.and.abs(dtset%iprcch)>=1.and.abs(dtset%iprcch)<=3) then
   call fresid(dtset,gresid,mpi_enreg,nfft,ngfft,ntypat,1,&
&   pawtab,rhor,rprimd,ucvol,vresid,xred,xred,psps%znuclpsp)
 else if (dtset%usewvl==0.and.(abs(dtset%iprcch)==4.or.abs(dtset%iprcch)==6)) then
   call fresidrsp(atindx1,dtset,gmet,gprimd,gresid,gsqcut,mgfft,&
&   mpi_enreg,psps%mqgrid_vl,nattyp,nfft,ngfft,ntypat,pawtab,ph1d,&
&   psps%qgrid_vl,ucvol,psps%usepaw,vresid,psps%zionpsp,psps%znuclpsp)
 else
   gresid(:,:)=zero
 end if

!=======================================================================
!======================= Other contributions ===========================
!=======================================================================

!Ewald energy contribution to forces as already been computed in "ewald"

!Potential residual contribution to forces as already been computed (forstr)

!Add Berry phase contributions (berryopt == 4)
!(compute the electric field force on the ion cores)
 if (dtset%berryopt==4) then
   ABI_ALLOCATE(fionred,(3,dtset%natom))
   fionred(:,:)=zero
   do iatom=1,dtset%natom
     itypat=dtset%typat(iatom)
     fioncart(:)=psps%ziontypat(itypat)*dtset%efield(:)
     do mu=1,3
       fionred(mu,iatom)=rprimd(1,mu)*fioncart(1) &
&       +rprimd(2,mu)*fioncart(2) &
&       +rprimd(3,mu)*fioncart(3)
     end do
   end do
 end if

!This was incorrect coding. Bug found by Jiawang Hong
!if (dtset%berryopt==4) then
!allocate(fionred(3,dtset%natom));fionred(:,:)=zero
!iatom = 0
!do itypat=1,ntypat
!do iattyp=1,nattyp(itypat)
!iatom=iatom+1
!fioncart(:)=psps%ziontypat(itypat)*dtset%efield(:)
!do mu=1,3
!fionred(mu,iatom)=rprimd(1,mu)*fioncart(1) &
!&         +rprimd(2,mu)*fioncart(2) &
!&         +rprimd(3,mu)*fioncart(3)
!end do
!end do
!end do
!end if

!=======================================================================
!======= Assemble the various contributions to the forces ==============
!=======================================================================

!write(std_out,*) "#grl", grl(:,:)
!write(std_out,*) "#gnl", grnl(:)
!write(std_out,*) "#gre", grewtn(:,:)
!write(std_out,*) "#syn", synlgr(:,:)
!write(std_out,*) "#grx", grxc(:,:)
!write(std_out,*) "#res", gresid(:,:)

!Collect grads of etot wrt reduced coordinates
!This gives non-symmetrized Hellman-Feynman reduced gradients
 ABI_ALLOCATE(grtn,(3,dtset%natom))
 grtn(:,:)=grl(:,:)+grewtn(:,:)+synlgr(:,:)+grxc(:,:)
 if (dtset%berryopt==4) grtn(:,:)=grtn(:,:)-fionred(:,:)

!write(std_out,*) "####### Gradients before sym #########"
!write(std_out,*) "#grt", grtn

!Symmetrize explicitly for given space group and store in grhf :
 call sygrad(grhf,dtset%natom,grtn,dtset%nsym,symrec,indsym)

!If residual forces are too large, there must be a problem: cancel them !
 if (dtset%usewvl==0.and.abs(dtset%iprcch)>0.and.abs(dtset%iprcch)/=5) then
   do iatom=1,dtset%natom
     do mu=1,3
       if (abs(gresid(mu,iatom))>10000._dp*abs(grtn(mu,iatom))) gresid(mu,iatom)=zero
     end do
   end do
 end if

!Add residual potential correction
 grtn(:,:)=grtn(:,:)+gresid(:,:)

!Additional stuff for electron-positron
 ipositron=0
 if (present(electronpositron)) then
   if (associated(electronpositron)) then
     if (associated(electronpositron%fred_ep)) ipositron=electronpositron_calctype(electronpositron)
   end if
 end if
 if (abs(ipositron)==1) then
   grtn(:,:)=grtn(:,:)-grxc(:,:)-grewtn(:,:)-gresid(:,:)-two*grl(:,:)
   grl(:,:)=-grl(:,:);grxc(:,:)=zero;gresid(:,:)=zero
   if (dtset%berryopt==4) grtn(:,:)=grtn(:,:)+fionred(:,:)
   if (dtset%berryopt==4) fionred(:,:)=zero
 end if
 if (ipositron>0) grtn(:,:)=grtn(:,:)+electronpositron%fred_ep(:,:)

!Symmetrize all grads explicitly for given space group:
 if (dtset%usewvl == 0) then
   call sygrad(fred,dtset%natom,grtn,dtset%nsym,symrec,indsym)
 else
   fred = grtn
 end if

!Conversion to cartesian coordinates (bohr) AND
!Subtract off average force from each force component
!to avoid spurious drifting of atoms across cell.
 call fred2fcart(favg,fcart,fred,gprimd,dtset%jellslab,dtset%natom)

!Compute maximal force and maximal difference
 maxfor=zero;diffor=zero
 do iatom=1,dtset%natom
   do mu=1,3
     if (dtset%iatfix(mu,iatom) /= 1) then
       maxfor=max(maxfor,abs(fcart(mu,iatom)))
       diffor=max(diffor,abs(fcart(mu,iatom)-fin(mu,iatom)))
     else if (dtset%ionmov==4 .or. dtset%ionmov==5) then
!      Make the force vanish on fixed atoms when ionmov=4 or 5
!      This is because fixing of atom cannot be imposed at the
!      level of a routine similar to brdmin or moldyn for these options.
       fcart(mu,iatom)=zero
     end if
   end do
 end do

!Apply any generalized constraints to the forces
 if (dtset%nconeq>0) call constrf(diffor,fcart,forold,fred,dtset%iatfix,dtset%ionmov,maxfor,&
& dtset%natom,dtset%nconeq,dtset%prtvol,rprimd,dtset%wtatcon,xred)

!=======================================================================
!Memory deallocations
 ABI_DEALLOCATE(grl)
 ABI_DEALLOCATE(grtn)
 ABI_DEALLOCATE(fin)
 if (dtset%berryopt==4)  then
   ABI_DEALLOCATE(fionred)
 end if

 call timab(69,2,tsec)

!DEBUG
!write(std_out,*)' forces: exit '
!stop
!ENDDEBUG

end subroutine forces
!!***
