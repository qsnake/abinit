!{\src2tex{textfont=tt}}
!!****f* ABINIT/getcprj
!! NAME
!! getcprj
!!
!! FUNCTION
!!  Compute <Proj_i|Cnk> for one wave function |Cnk> expressed in reciprocal space.
!!  Compute also derivatives of <Proj_i|Cnk>.
!!  |Proj_i> are non-local projectors (for each atom and each l,m,n)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  choice=chooses possible output:
!!    In addition to projected wave function:
!!    choice=1 => nothing else
!!          =2 => 1st gradients with respect to atomic position(s)
!!          =3 => 1st gradients with respect to strain(s)
!!          =23=> 1st gradients with respect to atm. pos. and strain(s)
!!          =4 => 2nd derivatives with respect to atomic pos.
!!          =24=> 1st and 2nd derivatives with respect to atomic pos.
!!          =5 => 1st gradients with respect to k wavevector
!!          =6 => 2nd derivatives with respect to strain and atm. pos.
!!  cpopt=1 if <Proj_i|Cnk> are already in memory; see below (side effects).
!!  cwavef(2,nspinor*npw_k)=input cmplx wavefunction coefficients <G|Cnk>
!!  dimekb1,dimekb2=dimensions of ekb (useful here only in non-PAW case)
!!  dimffnl=second dimension of ffnl
!!  ekb(dimekb1,dimekb2)= Kleinman-Bylander energies (hartree) (useful here only in non-PAW case)
!!  ffnl(npw_k,dimffnl,lmnmax,ntypat)=nonlocal form factors to be used for the application of the nl operator
!!  idir=direction of the derivative, i.e. dir. of - atom to be moved  in the case choice=2
!!                                                 - strain component  in the case choice=3
!!                                                 - k point direction in the case choice=5
!!       Compatible only with choice=2,3,5; if idir=0, all derivatives are computed
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=lmn
!!  istwf_k=option parameter that describes the storage of wfs
!!  kg_k(3,npw_k)=reduced planewave coordinates
!!  kpg(npw_k,npk)=(k+G) components and related data
!!  kpoint(3)=k point in terms of recip. translations
!!  lmnmax=max. number of (l,m,n) components over all types of atoms
!!  matblk=dimension of the array ph3d
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nattyp(ntypat)=number of atoms of each type
!!  ngfft(18)=contain all needed information about 3D FFT, see ~ABINIT/Infos/vargs.htm#ngfft
!!  nkpg=second size of array kpg
!!  nloalg(5)=governs the choice of the algorithm for nonlocal operator
!!  npw_k=number of planewaves for given k point
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  ntypat=number of types of atoms in unit cell
!!  phkxred(2,natom)=phase factors exp(2 pi kpoint.xred)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  ph3d(2,npw_k,matblk)=3D structure factors, for each atom and plane wave
!!   Note : is declared intent(inout), because it can be computed inside ph1d3d, and some compiler declare
!!   a conflict with the normal intent(in) of the present argument
!!  ucvol= unit cell volume
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  useylm=governs the way the nonlocal operator is to be applied
!!
!! SIDE EFFECTS
!!  cwaveprj(natom,nspinor) <type(cprj_type)>=projected input wave function <Proj_i|Cnk> with all NL projectors
!!                                (and derivatives)
!!                                if cpopt=1 the projected scalars have been already been computed and
!!                                           only derivatives are computed here
!!                                if cpopt=0 the projected scalars and derivatives are computed here
!!
!! TODO
!!  Spin-orbit
!!
!! PARENTS
!!      accrho3,berry_linemin,cgwf,ctocprj,debug_tools,ks_ddiago,m_wfs
!!
!! CHILDREN
!!      mkkpg,opernla_ylm,ph1d3d
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine getcprj(choice,cpopt,cwavef,cwaveprj,dimekb1,dimekb2,dimffnl,ekb,ffnl,&
&                   idir,indlmn,istwf_k,kg_k,kpg,kpoint,lmnmax,matblk,mgfft,mpi_enreg,&
&                   natom,nattyp,ngfft,nkpg,nloalg,npw_k,nspinor,ntypat,&
&                   phkxred,ph1d,ph3d,ucvol,usepaw,useylm)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getcprj'
 use interfaces_65_nonlocal, except_this_one => getcprj
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,cpopt,dimekb1,dimekb2,dimffnl,idir,istwf_k,lmnmax
 integer,intent(in) :: matblk,mgfft,natom,nkpg,npw_k,nspinor,ntypat,usepaw
 integer,intent(in) :: useylm
 real(dp) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),kg_k(3,npw_k),nattyp(ntypat)
 integer,intent(in) :: ngfft(18),nloalg(5)
 real(dp),intent(in) :: cwavef(2,npw_k*nspinor)
 real(dp),intent(in) :: ekb(dimekb1,dimekb2),ffnl(npw_k,dimffnl,lmnmax,ntypat),kpg(npw_k,nkpg)
 real(dp),intent(in) :: kpoint(3),ph1d(2,3*(2*mgfft+1)*natom),phkxred(2,natom)
 real(dp),intent(inout) :: ph3d(2,npw_k,matblk)
 type(cprj_type),intent(out) :: cwaveprj(natom,nspinor)

!Local variables-------------------------------
!scalars
 integer :: choice_,cplex,ia,ia1,ia2,ia3,ia4,iatm,ilmn,ispinor,itypat
 integer :: mincat,nd2gxdt,ndgxdt,nincat,nkpg_,nlmn,signs
 logical :: testnl
!arrays
 integer,allocatable :: indlmn_typ(:,:)
 real(dp),allocatable :: d2gxdt(:,:,:,:,:),dgxdt(:,:,:,:,:),ffnl_typ(:,:,:)
 real(dp),allocatable :: gx(:,:,:,:),kpg_(:,:)

! *********************************************************************

 DBG_ENTER('COLL')

!Nothing to do in that case
 if (cpopt==1.and.choice==1) return

!Not available for useylm=0
 if (useylm==0) then
   MSG_ERROR('  Not available for useylm=0 !')
 end if

!Error on bad choice
 if ((choice<1.or.choice>6).and.choice/=23.and.choice/=24) then
   MSG_BUG('  Does not presently support this choice !')
 end if

!Error on bad idir
 if (idir>0.and.choice/=2.and.choice/=3.and.choice/=5) then
   MSG_BUG('  Does not support idir>0 for this choice')
 end if

!Test: sizes of kpgin/kpgout
 if (nkpg>0.and. &
& ( (choice==2.and.nkpg<3) .or. &
& ((choice==4.or.choice==24).and.nkpg<9) .or. &
& ((choice==6.or.choice==3.or.choice==23).and.nkpg<3) )) then
   MSG_BUG('  Incorrect size for nkpg array !')
 end if

!Define dimensions of projected scalars
 ndgxdt=0;nd2gxdt=0
 if (idir==0) then
   if (choice==2) ndgxdt=3
   if (choice==3) ndgxdt=6
   if (choice==23) ndgxdt=9
   if (choice==4) nd2gxdt=6
   if (choice==24) then
     ndgxdt=3;nd2gxdt=6
   end if
   if (choice==5) ndgxdt=3
   if (choice==6) then
     ndgxdt=9;nd2gxdt=54
   end if
 else
   ndgxdt=1
 end if

!Eventually re-compute (k+G) vectors (and related data)
 nkpg_=0
 if (nkpg==0) then
   if (choice==4.or.choice==24) nkpg_=9
   if (choice==2.or.choice==3.or.choice==23) nkpg_=3
   if (nkpg_>0) then
     ABI_ALLOCATE(kpg_,(npw_k,nkpg_))
     call mkkpg(kg_k,kpg_,kpoint,nkpg_,npw_k)
   end if
 end if

!Some other dims
 mincat=min(nloalg(4),maxval(nattyp))
 cplex=2;if (istwf_k>1) cplex=1
 choice_=choice;if (cpopt==1) choice_=-choice
 signs=1;if (idir>0) signs=2

!Loop over atom types
 ia1=1;iatm=0
 do itypat=1,ntypat
   ia2=ia1+nattyp(itypat)-1;if (ia2<ia1) cycle
   nlmn=count(indlmn(3,:,itypat)>0)

!  Test on local part
   testnl=(usepaw>0)
   if (usepaw==0) then
     do ilmn=1,nlmn
       if (abs(ekb(indlmn(5,ilmn,itypat),itypat))>tol10) testnl=.true.
     end do
   end if
   if (.not.testnl) then
     if (cpopt>=0) then
       do ispinor=1,nspinor
         do ia=1,nattyp(itypat)
           if (cpopt==0) cwaveprj(iatm+ia,ispinor)%cp(:,1:nlmn)=zero
           if (choice>1) cwaveprj(iatm+ia,ispinor)%dcp(:,:,1:nlmn)=zero
         end do
       end do
     end if
     iatm=iatm+nattyp(itypat)
   else

!    Retrieve some data for this type of atom
     ABI_ALLOCATE(indlmn_typ,(6,nlmn))
     ABI_ALLOCATE(ffnl_typ,(npw_k,dimffnl,nlmn))
     indlmn_typ(:,1:nlmn)=indlmn(:,1:nlmn,itypat)
     ffnl_typ(:,:,1:nlmn)=ffnl(:,:,1:nlmn,itypat)

!    Loop on blocks of atoms inside type
     do ia3=ia1,ia2,mincat
       ia4=min(ia2,ia3+mincat-1);nincat=ia4-ia3+1

!      Prepare the phase factors if they were not already computed
       if (nloalg(1)<=0) call ph1d3d(ia3,ia4,kg_k,matblk,natom,&
&       npw_k,ngfft(1),ngfft(2),ngfft(3),phkxred,ph1d,ph3d)

!      Allocate memory for projected scalars
       ABI_ALLOCATE(gx,(cplex,nlmn,nincat,nspinor))
       ABI_ALLOCATE(dgxdt,(cplex,ndgxdt,nlmn,nincat,nspinor))
       ABI_ALLOCATE(d2gxdt,(cplex,nd2gxdt,nlmn,nincat,nspinor))

!      Retrieve eventually <p_i|c> coeffs
       if (cpopt==1) then
         do ispinor=1,nspinor
           do ia=1,nincat
             gx(1:cplex,1:nlmn,ia,ispinor)=cwaveprj(iatm+ia,ispinor)%cp(1:cplex,1:nlmn)
           end do
         end do
       end if

!      Compute <p_i|c> scalars (and derivatives) for this block of atoms
       if (nkpg_>0) then
         call opernla_ylm(choice_,cplex,dimffnl,d2gxdt,dgxdt,ffnl_typ,gx,ia3,idir,indlmn_typ,&
&         istwf_k,kpg_,matblk,mpi_enreg,nd2gxdt,ndgxdt,nincat,nkpg_,nlmn,&
&         nloalg,npw_k,nspinor,ph3d,signs,ucvol,cwavef)
       else
         call opernla_ylm(choice_,cplex,dimffnl,d2gxdt,dgxdt,ffnl_typ,gx,ia3,idir,indlmn_typ,&
&         istwf_k,kpg ,matblk,mpi_enreg,nd2gxdt,ndgxdt,nincat,nkpg ,nlmn,&
&         nloalg,npw_k,nspinor,ph3d,signs,ucvol,cwavef)
       end if

!      Transfer result to output variable cwaveprj
       if (cpopt==0) then
         do ispinor=1,nspinor
           do ia=1,nincat
             cwaveprj(iatm+ia,ispinor)%nlmn=nlmn
             cwaveprj(iatm+ia,ispinor)%cp(1:cplex,1:nlmn)=gx(1:cplex,1:nlmn,ia,ispinor)
             if (cplex==1) cwaveprj(iatm+ia,ispinor)%cp(2,1:nlmn)=zero
           end do
         end do
       end if
       if (cpopt>=0.and.choice>1) then
         do ispinor=1,nspinor
           do ia=1,nincat
             cwaveprj(iatm+ia,ispinor)%ncpgr=ndgxdt+nd2gxdt
             if (ndgxdt>0) cwaveprj(iatm+ia,ispinor)%dcp(1:cplex,1:ndgxdt,1:nlmn)=&
&             dgxdt(1:cplex,1:ndgxdt,1:nlmn,ia,ispinor)
             if (nd2gxdt>0)cwaveprj(iatm+ia,ispinor)%dcp(1:cplex,ndgxdt+1:ndgxdt+nd2gxdt,1:nlmn)=&
&             d2gxdt(1:cplex,1:nd2gxdt,1:nlmn,ia,ispinor)
             if (cplex==1) cwaveprj(iatm+ia,ispinor)%dcp(2,1:ndgxdt+nd2gxdt,1:nlmn)=zero
           end do
         end do
       end if

!      End loop inside block of atoms
       iatm=iatm+nincat
       ABI_DEALLOCATE(gx)
       ABI_DEALLOCATE(dgxdt)
       ABI_DEALLOCATE(d2gxdt)
     end do
!    End IF nonlocal part exists
     ABI_DEALLOCATE(indlmn_typ)
     ABI_DEALLOCATE(ffnl_typ)
   end if
!  End loop over atom types
   ia1=ia2+1
 end do

 DBG_EXIT('COLL')

 end subroutine getcprj
!!***
