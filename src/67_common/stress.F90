!{\src2tex{textfont=tt}}
!!****f* ABINIT/stress
!!
!! NAME
!! stress
!!
!! FUNCTION
!! Compute the stress tensor
!! strten(i,j) = (1/ucvol)*d(Etot)/(d(eps(i,j)))
!! where Etot is energy per unit cell, ucvol is the unstrained unit cell
!! volume, r(i,iat) is the ith position of atom iat,
!! and eps(i,j) is an infinitesimal strain which maps each
!! point r to r(i) -> r(i) + Sum(j) [eps(i,j)*r(j)].
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
!!  berryopt  = 4: electric field is on -> add the contribution of the
!!                 - \Omega E.P term to the total energy
!!           /= 4: electric field is off
!!   from Etot(npw) data (at fixed geometry), used for making
!!   Pulay correction to stress tensor (hartree).  Should be <=0.
!!  eei=local pseudopotential part of Etot (hartree)
!!  efield = cartesian coordinates of the electric field in atomic units
!!  ehart=Hartree energy (hartree)
!!  eii=pseudoion core correction energy part of Etot (hartree)
!!  gsqcut=cutoff value on G**2 for (large) sphere inside FFT box.
!!                       gsqcut=(boxcut**2)*ecut/(2._dp*(Pi**2)
!!  kinstr(6)=kinetic energy part of stress tensor
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mqgrid=dimensioned number of q grid points for local psp spline
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  natom=number of atoms in cell
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nlstr(6)=nonlocal part of stress tensor
!!  nspden=number of spin-density components
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pel(3)=reduced coordinates of the electronic polarization (a. u.)
!!  pion(3)=reduced coordinates of the ionic polarization (a. u.)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase (structure factor) array
!!  prtvol=integer controlling volume of printed output
!!  qgrid(mqgrid)=q point array for local psp spline fits
!!  rhog(2,nfft)=Fourier transform of charge density (bohr^-3)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  strsxc(6)=xc correction to stress
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  typat(natom)=type integer for each atom in cell
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vlspl(mqgrid,2,ntypat)=local psp spline
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree) in real space
!!  xccc1d(n1xccc*(1-usepaw),6,ntypat)=1D core charge function and five derivatives,
!!                          for each type of atom, from psp (used in Norm-conserving only)
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  zion(ntypat)=valence charge of each type of atom
!!
!! OUTPUT
!!  strten(6)=components of the stress tensor (hartree/bohr^3) for the
!!    6 unique components of this symmetric 3x3 tensor:
!!    Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
!!    The diagonal components of the returned stress tensor are
!!    CORRECTED for the Pulay stress.
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!
!! NOTES
!! * Concerning the stress tensor:
!!   See O. H. Nielsen and R. M. Martin, PRB 32, 3792 (1985).
!!   Note that first term in equation (2) should have minus sign
!!   (for kinetic energy contribution to stress tensor).
!!   Normalizations in this code differ somewhat from those employed
!!   by Nielsen and Martin.
!!   For the stress tensor contribution from the nonlocal Kleinman-Bylander
!!   separable pseudopotential, see D. M. Bylander, L. Kleinman, and
!!   S. Lee, PRB 42, 1394 (1990).
!!   Again normalization conventions differ somewhat.
!!   See Doug Allan s notes starting page 795 (13 Jan 1992).
!! * This subroutine calls different subroutines to compute the stress
!!   tensor contributions from the following parts of the total energy:
!!   (1) kinetic energy, (2) exchange-correlation energy,
!!   (3) Hartree energy, (4) local pseudopotential energy,
!!   (5) pseudoion core correction energy, (6) nonlocal pseudopotential energy,
!!   (7) Ewald energy.
!!
!! PARENTS
!!      forstr
!!
!! CHILDREN
!!      atm2fft,ewald2,fourdp,metric,mkcore,mklocl_recipspace,stresssym,strhar
!!      timab,wrtout,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine stress(atindx1,berryopt,eei,efield,ehart,eii,gsqcut,kinstr,&
&                  mgfft,mpi_enreg,mqgrid,n1xccc,n3xccc,natom,nattyp,&
&                  nfft,ngfft,nlstr,nspden,nsym,ntypat,paral_kgb,pawtab,pel,pion,ph1d,&
&                  prtvol,qgrid,rhog,rprimd,strten,strsxc,symrec,typat,usepaw,vlspl,&
&                  vxc,xccc1d,xccc3d,xcccrc,xred,zion,&
&                  electronpositron) ! optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'stress'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_42_geometry
 use interfaces_53_ffts
 use interfaces_56_xc
 use interfaces_65_psp
 use interfaces_67_common, except_this_one => stress
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: berryopt,mgfft,mqgrid,n1xccc,n3xccc,natom,nfft,nspden
 integer,intent(in) :: nsym,ntypat,paral_kgb,prtvol,usepaw
 real(dp),intent(in) :: eei,ehart,eii,gsqcut
 type(electronpositron_type),pointer,optional :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat),ngfft(18),symrec(3,3,nsym)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: efield(3),kinstr(6),nlstr(6),pel(3)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),pion(3),qgrid(mqgrid)
 real(dp),intent(in) :: rhog(2,nfft),rprimd(3,3),strsxc(6)
 real(dp),intent(in) :: vlspl(mqgrid,2,ntypat),vxc(nfft,nspden)
 real(dp),intent(in) :: xccc1d(n1xccc*(1-usepaw),6,ntypat),xcccrc(ntypat)
 real(dp),intent(in) :: xred(3,natom),zion(ntypat)
 real(dp),intent(inout) :: xccc3d(n3xccc)
 real(dp),intent(out) :: strten(6)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: ii,ipositron,mu,optatm,optdyfr,optgr,option,optn,optn2,optstr,optv
 real(dp),parameter :: tol=1.0d-15
 real(dp) :: strsii,ucvol
 character(len=500) :: message
!arrays
 integer :: qprtrb_dum(3)
 real(dp) :: berrystr(6),berrystr_tmp(3,3),corstr(6),ewestr(6),gmet(3,3)
 real(dp) :: gprimd(3,3),harstr(6),lpsstr(6),rmet(3,3),tsec(2),uncorr(3)
 real(dp) :: vprtrb_dum(2)
 real(dp),allocatable :: dummy(:),dummy_in(:),dummy_out(:),&
& dyfr_dum(:,:,:),gr_dum(:,:),rhog_ep(:,:),v_dum(:)
 real(dp),allocatable :: vxctotg(:,:)
 character(len=10) :: EPName(1:2)=(/"Electronic","Positronic"/)

! *************************************************************************

!DEBUG
!write(std_out,*)' stress: enter '
!stop
!ENDDEBUG

 call timab(37,1,tsec)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!=======================================================================
!========= Local pseudopotential and core charge contributions =========
!=======================================================================

 if (usepaw==1) then

!  PAW: compute local psp and core charge contribs together in reciprocal space
   call timab(551,1,tsec)
   if (n3xccc>0) then
     ABI_ALLOCATE(v_dum,(nfft))
     ABI_ALLOCATE(vxctotg,(2,nfft))
     v_dum(:)=vxc(:,1);if (nspden>=2) v_dum(:)=0.5_dp*(v_dum(:)+vxc(:,2))
     call fourdp(1,vxctotg,v_dum,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     call zerosym(vxctotg,2,mpi_enreg,ngfft(1),ngfft(2),ngfft(3))
     ABI_DEALLOCATE(v_dum)
   end if
   optatm=0;optdyfr=0;optgr=0;optstr=1;optv=1;optn=n3xccc/nfft;optn2=1
   call atm2fft(atindx1,dummy_out,dummy_out,dummy_out,dummy_out,&
&   eei,dummy_in,gmet,gprimd,dummy_out,dummy_out,&
&   gsqcut,mgfft,mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,&
&   optatm,optdyfr,optgr,optn,optn2,optstr,optv,paral_kgb,&
&   pawtab,ph1d,qgrid,qprtrb_dum,rhog,corstr,lpsstr,&
&   ucvol,usepaw,vxctotg,vprtrb_dum,vlspl)
   if (n3xccc>0)  then
     ABI_DEALLOCATE(vxctotg)
   end if
   if (n3xccc==0) corstr=zero
   call timab(551,2,tsec)

 else

!  Norm-conserving: compute local psp contribution in reciprocal space
!  and core charge contribution in real space
   option=3
   ABI_ALLOCATE(dyfr_dum,(3,3,natom))
   ABI_ALLOCATE(gr_dum,(3,natom))
   ABI_ALLOCATE(v_dum,(nfft))
   call mklocl_recipspace(dyfr_dum,eei,gmet,gprimd,gr_dum,gsqcut,lpsstr,mgfft,&
&   mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,option,paral_kgb,ph1d,qgrid,&
&   qprtrb_dum,rhog,ucvol,vlspl,vprtrb_dum,v_dum)
   if (n3xccc>0) then
     call timab(55,1,tsec)
     call mkcore(corstr,dyfr_dum,gr_dum,mpi_enreg,natom,nfft,nspden,ntypat,ngfft(1),&
&     n1xccc,ngfft(2),ngfft(3),option,rprimd,typat,ucvol,vxc,&
&     xcccrc,xccc1d,xccc3d,xred)
     call timab(55,2,tsec)
   else
     corstr(:)=zero
   end if
   ABI_DEALLOCATE(dyfr_dum)
   ABI_DEALLOCATE(gr_dum)
   ABI_DEALLOCATE(v_dum)
 end if

!=======================================================================
!======================= Hartree energy contribution ===================
!=======================================================================

 call strhar(ehart,gprimd,gsqcut,harstr,mpi_enreg,nfft,ngfft,rhog,ucvol)

!=======================================================================
!======================= Ewald contribution ============================
!=======================================================================

 call timab(38,1,tsec)
 call ewald2(gmet,natom,ntypat,rmet,rprimd,ewestr,typat,ucvol,xred,zion)
 call timab(38,2,tsec)

!=======================================================================
!=================== Berry phase contribution ==========================
!=======================================================================

 if (berryopt==4) then
   berrystr_tmp(:,:) = zero
!  Diagonal:
   do mu = 1, 3
     do ii = 1, 3
       berrystr_tmp(mu,mu) = berrystr_tmp(mu,mu) - &
&       efield(mu)*rprimd(mu,ii)*(pel(ii) + pion(ii))/ucvol
     end do
   end do
!  Off-diagonal (symmetrized before adding it to strten):
   do ii = 1, 3
     berrystr_tmp(3,2) = berrystr_tmp(3,2) &
&     - efield(3)*rprimd(2,ii)*(pel(ii) + pion(ii))/ucvol
     berrystr_tmp(2,3) = berrystr_tmp(2,3) &
&     - efield(2)*rprimd(3,ii)*(pel(ii) + pion(ii))/ucvol
     berrystr_tmp(3,1) = berrystr_tmp(3,1) &
&     - efield(3)*rprimd(1,ii)*(pel(ii) + pion(ii))/ucvol
     berrystr_tmp(1,3) = berrystr_tmp(1,3) &
&     - efield(1)*rprimd(3,ii)*(pel(ii) + pion(ii))/ucvol
     berrystr_tmp(2,1) = berrystr_tmp(2,1) &
&     - efield(2)*rprimd(1,ii)*(pel(ii) + pion(ii))/ucvol
     berrystr_tmp(1,2) = berrystr_tmp(1,2) &
&     - efield(1)*rprimd(2,ii)*(pel(ii) + pion(ii))/ucvol
   end do
   berrystr(1) = berrystr_tmp(1,1)
   berrystr(2) = berrystr_tmp(2,2)
   berrystr(3) = berrystr_tmp(3,3)
   berrystr(4) = (berrystr_tmp(3,2) + berrystr_tmp(2,3))/two
   berrystr(5) = (berrystr_tmp(3,1) + berrystr_tmp(1,3))/two
   berrystr(6) = (berrystr_tmp(2,1) + berrystr_tmp(1,2))/two
 end if

!=======================================================================
!================= Other (trivial) contributions =======================
!=======================================================================

!Nonlocal part of stress has already been computed
!(in forstrnps(norm-conserving) or pawgrnl(PAW))

!Kinetic part of stress has already been computed
!(in forstrnps)

!XC part of stress tensor has already been computed in "strsxc"

!ii part of stress (diagonal) is trivial
 strsii=-eii/ucvol

!=======================================================================
!===== Assemble the various contributions to the stress tensor =========
!=======================================================================
!In cartesian coordinates (symmetric storage)

 strten(:)=harstr(:)+lpsstr(:)+kinstr(:)+nlstr(:)+ewestr(:)+corstr(:)+strsxc(:)
 if (berryopt==4) strten(:)=strten(:)+berrystr(:)

!Additional stuff for electron-positron
 ipositron=0
 if (present(electronpositron)) then
   if (associated(electronpositron)) then
     if (associated(electronpositron%stress_ep)) ipositron=electronpositron_calctype(electronpositron)
   end if
 end if
 if (abs(ipositron)==1) then
   strten(:)=strten(:)-harstr(:)-ewestr(:)-corstr(:)-lpsstr(:)
   harstr(:)=zero;ewestr(:)=zero;corstr(:)=zero;strsii=zero
   lpsstr(:)=-lpsstr(:);lpsstr(1:3)=lpsstr(1:3)-two*eei/ucvol
   strten(:)=strten(:)+lpsstr(:)
   if (berryopt==4) strten(:)=strten(:)-berrystr(:)
   if (berryopt==4) berrystr(:)=zero
 end if
 if (abs(ipositron)==2) then
   ABI_ALLOCATE(rhog_ep,(2,nfft))
   ABI_ALLOCATE(dummy,(6))
   call fourdp(1,rhog_ep,electronpositron%rhor_ep,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   rhog_ep=-rhog_ep
   call strhar(electronpositron%e_hartree,gprimd,gsqcut,dummy,mpi_enreg,nfft,ngfft,rhog_ep,ucvol)
   strten(:)=strten(:)+dummy(:);harstr(:)=harstr(:)+dummy(:)
   ABI_DEALLOCATE(rhog_ep)
   ABI_DEALLOCATE(dummy)
 end if
 if (ipositron>0) strten(:)=strten(:)+electronpositron%stress_ep(:)

!Symmetrize resulting tensor if nsym>1
 if (nsym>1) call stresssym(gprimd,nsym,strten,symrec)

!Set to zero very small values of stress
 do mu=1,6
   if (abs(strten(mu))<tol) strten(mu)=zero
 end do

!Include diagonal terms, save uncorrected stress for output
 do mu=1,3
   uncorr(mu)=strten(mu)+strsii
   strten(mu)=uncorr(mu)
 end do

!=======================================================================
!================ Print out info about stress tensor ===================
!=======================================================================
 if (prtvol>=10.and.ipositron>=0) then
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' )&
&     ' stress: component',mu,' of hartree stress is',harstr(mu)
     call wrtout(std_out,message,'COLL')
   end do
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' )&
&     ' stress: component',mu,' of loc psp stress is',lpsstr(mu)
     call wrtout(std_out,message,'COLL')
   end do
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' )&
&     ' stress: component',mu,&
&     ' of kinetic stress is',kinstr(mu)
     call wrtout(std_out,message,'COLL')
   end do
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' )&
&     ' stress: component',mu,' of nonlocal ps stress is',nlstr(mu)
     call wrtout(std_out,message,'COLL')
   end do
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' )&
&     ' stress: component',mu,' of     core xc stress is',corstr(mu)
     call wrtout(std_out,message,'COLL')
   end do
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' )&
&     ' stress: component',mu,&
&     ' of Ewald energ stress is',ewestr(mu)
     call wrtout(std_out,message,'COLL')
   end do
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' ) &
&     ' stress: component',mu,' of xc stress is',strsxc(mu)
     call wrtout(std_out,message,'COLL')
   end do
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   write(message, '(a,1p,e22.12)' ) &
&   ' stress: ii (diagonal) part is',strsii
   call wrtout(std_out,message,'COLL')
   if (berryopt==4) then
     write(message, '(a)' ) ' '
     call wrtout(std_out,message,'COLL')
     do mu = 1, 3
       do ii = 1, 3
         write(message, '(a,i2,i2,a,1p,e22.12)' )&
&         ' stress: component',mu,ii,' of unsymmetrized Berry phase stress is',&
&         berrystr_tmp(mu,ii)
         call wrtout(std_out,message,'COLL')
       end do
     end do
   end if
   if (ipositron/=0) then
     write(message, '(a)' ) ' '
     call wrtout(std_out,message,'COLL')
     do mu=1,6
       write(message, '(a,i5,3a,1p,e22.12)' ) &
&       ' stress: component',mu,' of ',EPName(abs(ipositron)), &
&       ' stress is',electronpositron%stress_ep(mu)
       call wrtout(std_out,message,'COLL')
     end do
   end if

 end if ! prtvol
 if (ipositron>=0) then
   write(message, '(a,a)' )ch10,&
&   ' Cartesian components of stress tensor (hartree/bohr^3)'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(1 1)=',strten(1),'  sigma(3 2)=',strten(4)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(2 2)=',strten(2),'  sigma(3 1)=',strten(5)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(3 3)=',strten(3),'  sigma(2 1)=',strten(6)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a)' ) ' '
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

 call timab(37,2,tsec)

!DEBUG
!write(std_out,*)' stress : exit '
!if(.true.)stop
!ENDDEBUG

end subroutine stress
!!***
