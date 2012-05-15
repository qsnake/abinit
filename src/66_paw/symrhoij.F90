!{\src2tex{textfont=tt}}
!!****f* ABINIT/symrhoij
!! NAME
!! symrhoij
!!
!! FUNCTION
!! Symmetrize rhoij quantities (augmentation occupancies) and/or gradients
!! Compute also rhoij residuals (new-old values of rhoij and gradients)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  choice=select then type of rhoij gradients to symmetrize.
!!         choice=1 => no gradient
!!         choice=2 => gradient with respect to atomic position(s)
!!               =3 => a gradient with respect to strain(s)
!!               =4 => 2nd gradient with respect to atomic position(s)
!!               =23=> a gradient with respect to atm. pos. and strain(s)
!!               =24=> 1st and 2nd gradient with respect to atomic position(s)
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  indlmn(6,lmnmax,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn (for each atom type)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  ipert=index of perturbation if pawrhoij is a pertubed rhoij
!!        no meaning for ground-state calculations (should be 0)
!!  lmnmax=maximum number of PAW radial wavefunctions
!!  natom=number of atoms in cell
!!  nrhoij=size of pawrhoij (generally nrhoij=natom or 1)
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of types of atoms in unit cell.
!!  optrhoij= 1 if rhoij quantities have to be symmetrized
!!  pawrhoij(natom)%cplex=1 if rhoij are REAL, 2 if they are COMPLEX
!!  pawrhoij(natom)%lmn_size=number of (l,m,n) elements for the paw basis
!!  pawrhoij(natom)%nspden=number of spin-density components
!!  pawrhoij(natom)%nsppol=number of independant spin-density components
!!  pawrhoij(natom)%rhoij_(cplex*lmn2_size,nspden)=non-symetrized paw rhoij quantities
!!  pawang <type(pawang_type)>=angular mesh discretization and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!            Note: if pawprtvol=-10001, nothing is printed out
!!  rprimd(3,3)=real space primitive translations.
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!                   reciprocal space primitive translations
!!  typat(natom)=type for each atom
!!
!! OUTPUT
!!  if (optrhoij==1)
!!    pawrhoij(natom)%nrhoijsel=number of non-zero values of rhoij
!!    pawrhoij(natom)%rhoijp(cplex*lmn2_size,nspden)=symetrized paw rhoij quantities in PACKED STORAGE (only non-zero values)
!!    pawrhoij(natom)%rhoijres(cplex*lmn2_size,nspden)=paw rhoij quantities residuals (new values - old values)
!!    pawrhoij(natom)%rhoijselect(lmn2_size)=select the non-zero values of rhoij
!!
!! SIDE EFFECTS
!!  if (pawrhoij(:)%ngrhoij>0) (equivalent to choice>1)
!!    At input:
!!    pawrhoij(natom)%grhoij(ngrhoij,cplex*lmn2_size,nspden)=non-symetrized gradients of rhoij
!!    At output:
!!    pawrhoij(natom)%grhoij(ngrhoij,cplex*lmn2_size,nspden)=symetrized gradients of rhoij
!!
!! PARENTS
!!      dyfnl3,energy,paw_qpscgw,pawmkrho
!!
!! CHILDREN
!!      print_ij,symredcart,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine symrhoij(choice,gprimd,indlmn,indsym,ipert,lmnmax,natom,nrhoij,nsym,ntypat,optrhoij,&
&                   pawang,pawprtvol,pawrhoij,rprimd,symafm,symrec,typat)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symrhoij'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: choice,ipert,lmnmax,natom,nrhoij,nsym,ntypat,optrhoij,pawprtvol
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),indsym(4,nsym,natom)
 integer,intent(in) :: symafm(nsym),symrec(3,3,nsym),typat(natom)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 type(pawrhoij_type),intent(inout) :: pawrhoij(nrhoij)

!Local variables ---------------------------------------
 character(len=8),parameter :: dspin(6)=(/"up      ","down    ","dens (n)","magn (x)","magn (y)","magn (z)"/)
!scalars
 integer :: at_indx,cplex,cplex_eff,iafm,iatm,iatom,idum1,idum2,il,il0,ilmn,iln,iln0,ilpm,indexi
 integer :: indexii,indexj,indexjj,indexjj0,indexk,indexk1,iplex,irhoij,irot,ishift2,ishift3
 integer :: ishift4,ispden,itypat,j0lmn,jj,jl,jl0,jlmn,jln,jln0,jlpm,jrhoij,jspden,klmn,klmn1,kspden
 integer :: lmn_size,mi,mj,mu,mua,mub,mushift,natinc,ngrhoij,nselect,nselect1,nspinor,nu,nushift
 logical,parameter :: afm_noncoll=.true.  ! TRUE if antiferro symmetries are used with non-collinear magnetism
 real(dp) :: factafm,syma,zarot2
 logical :: antiferro,noncoll,use_afm,use_res
 character(len=8) :: pertstrg
 character(len=500) :: message
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer :: nsym_used(2)
 integer,allocatable :: idum(:)
 real(dp) :: ro(2),sumrho(2,2),sum1(2),xsym(3)
 real(dp),allocatable :: rotgr(:,:,:),rotmag(:,:),rotmaggr(:,:,:)
 real(dp),allocatable :: sumgr(:,:),summag(:,:),summaggr(:,:,:),symrec_cart(:,:,:),work1(:,:,:)
 real(dp),allocatable :: factsym(:)
 type(coeff3_type),allocatable :: tmp_grhoij(:)

! *********************************************************************

 DBG_ENTER("COLL")

!Symetrization occurs only when nsym>1
 if (nsym>1) then

!  Test: consistency between choice and ngrhoij
   ngrhoij=pawrhoij(1)%ngrhoij
   if ((choice==1.and.ngrhoij/=0) .or.(choice==2.and.ngrhoij/=3).or. &
&   (choice==3.and.ngrhoij/=6).or.(choice==23.and.ngrhoij/=9).or. &
&   (choice==4.and.ngrhoij/=6).or.(choice==24.and.ngrhoij/=9) ) then
     message='  Inconsistency between variables choice and ngrhoij !'
     MSG_BUG(message)
   end if

!  Symetrization of gradients not compatible with nspden=4
   if (choice>2.and.pawrhoij(1)%nspden==4) then
     message='  For the time being, choice>2 is not compatible with nspden=4 !'
     MSG_BUG(message)
   end if

!  Symetry matrixes must be in memory
   if (pawang%nsym==0) then
     message='  pawang%zarot must be allocated !'
     MSG_BUG(message)
   end if

!  Antiferro case ?
   antiferro=(pawrhoij(1)%nspden==2.and.pawrhoij(1)%nsppol==1)
!  Non-collinear case
   noncoll=(pawrhoij(1)%nspden==4)
!  Do we use antiferro symmetries ?
   use_afm=((antiferro).or.(noncoll.and.afm_noncoll))

!  Several inits/allocations
   cplex_eff=1
   if (ipert>0.or.antiferro.or.noncoll) cplex_eff=pawrhoij(1)%cplex ! Does not symmetrize imaginary part for GS calculations
   ABI_ALLOCATE(factsym,(cplex_eff))
   if (noncoll.and.optrhoij==1)  then
     ABI_ALLOCATE(summag,(cplex_eff,3))
     ABI_ALLOCATE(rotmag,(cplex_eff,3))
   end if
   if (noncoll) then
     ABI_ALLOCATE(symrec_cart,(3,3,nsym))
     do irot=1,nsym
       call symredcart(gprimd,rprimd,symrec_cart(:,:,irot),symrec(:,:,irot))
     end do
   end if
   ishift2=0;ishift3=0;ishift4=0
   if (choice>1) then
     ABI_ALLOCATE(sumgr,(cplex_eff,ngrhoij))
     if (choice>2)  then
       ABI_ALLOCATE(work1,(cplex_eff,3,3))
     end if
     if (antiferro) then
       ABI_ALLOCATE(rotgr,(cplex_eff,ngrhoij,2))
     else
       ABI_ALLOCATE(rotgr,(cplex_eff,ngrhoij,1))
     end if
     if (noncoll) then
       ABI_ALLOCATE(summaggr,(cplex_eff,ngrhoij,3))
       ABI_ALLOCATE(rotmaggr,(cplex_eff,ngrhoij,3))
     end if
     if (choice==23) ishift2=6
     if (choice==24) ishift4=3
!    Have to make a temporary copy of grhoij
     ABI_ALLOCATE(tmp_grhoij,(nrhoij))
     do iatm=1,nrhoij
       idum1=pawrhoij(iatm)%cplex*pawrhoij(iatm)%lmn2_size;idum2=pawrhoij(iatm)%nspden
       ABI_ALLOCATE(tmp_grhoij(iatm)%value,(ngrhoij,idum1,idum2))
       tmp_grhoij(iatm)%value(1:ngrhoij,1:idum1,1:idum2)=pawrhoij(iatm)%grhoij(1:ngrhoij,1:idum1,1:idum2)
     end do
   end if


!  Loops over atoms and spin components
!  ------------------------------------
   do iatm=1,nrhoij
     iatom=iatm;if (nrhoij==1.and.ipert>0.and.ipert<=natom) iatom=ipert
     itypat=typat(iatom)
     lmn_size=pawrhoij(iatm)%lmn_size
     nspinor=pawrhoij(iatm)%nspinor
     cplex=pawrhoij(iatm)%cplex
     cplex_eff=1;if (ipert>0.or.antiferro.or.noncoll) cplex_eff=cplex
     use_res=(pawrhoij(iatm)%use_rhoijres>0)

     nselect=0
     do ispden=1,pawrhoij(iatm)%nsppol
       jspden=min(3-ispden,pawrhoij(iatm)%nsppol)

!      Store old -rhoij in residual
       if (optrhoij==1.and.use_res) then
         pawrhoij(iatm)%rhoijres(:,ispden)=zero
         if (cplex==1) then
           do irhoij=1,pawrhoij(iatm)%nrhoijsel
             klmn=pawrhoij(iatm)%rhoijselect(irhoij)
             pawrhoij(iatm)%rhoijres(klmn,ispden)=-pawrhoij(iatm)%rhoijp(irhoij,ispden)
           end do
         else
           do irhoij=1,pawrhoij(iatm)%nrhoijsel
             klmn1=2*pawrhoij(iatm)%rhoijselect(irhoij);jrhoij=2*irhoij
             pawrhoij(iatm)%rhoijres(klmn1-1:klmn1,ispden)=-pawrhoij(iatm)%rhoijp(jrhoij-1:jrhoij,ispden)
           end do
         end if
         if (noncoll) then
           pawrhoij(iatm)%rhoijres(:,2:4)=zero
           if (cplex==1) then
             do mu=2,4
               do irhoij=1,pawrhoij(iatm)%nrhoijsel
                 klmn=pawrhoij(iatm)%rhoijselect(irhoij)
                 pawrhoij(iatm)%rhoijres(klmn,mu)=-pawrhoij(iatm)%rhoijp(irhoij,mu)
               end do
             end do
           else
             do mu=2,4
               do irhoij=1,pawrhoij(iatm)%nrhoijsel
                 klmn1=2*pawrhoij(iatm)%rhoijselect(irhoij);jrhoij=2*irhoij
                 pawrhoij(iatm)%rhoijres(klmn1-1:klmn1,mu)=-pawrhoij(iatm)%rhoijp(jrhoij-1:jrhoij,mu)
               end do
             end do
           end if
         end if
         if (antiferro) then
           pawrhoij(iatm)%rhoijres(:,2)=zero
           if (cplex==1) then
             do irhoij=1,pawrhoij(iatm)%nrhoijsel
               klmn=pawrhoij(iatm)%rhoijselect(irhoij)
               pawrhoij(iatm)%rhoijres(klmn,2)=-pawrhoij(iatm)%rhoijp(irhoij,2)
             end do
           else
             do irhoij=1,pawrhoij(iatm)%nrhoijsel
               klmn1=2*pawrhoij(iatm)%rhoijselect(irhoij);jrhoij=2*irhoij
               pawrhoij(iatm)%rhoijres(klmn1-1:klmn1,2)=-pawrhoij(iatm)%rhoijp(jrhoij-1:jrhoij,2)
             end do
           end if
         end if
       end if

!      Loops over (il,im) and (jl,jm)
!      ------------------------------
       jl0=-1;jln0=-1;indexj=1
       do jlmn=1,lmn_size
         jl=indlmn(1,jlmn,itypat)
         jlpm=1+jl+indlmn(2,jlmn,itypat)
         jln=indlmn(5,jlmn,itypat)
         if (jln/=jln0) indexj=indexj+2*jl0+1
         j0lmn=jlmn*(jlmn-1)/2
         il0=-1;iln0=-1;indexi=1
         do ilmn=1,jlmn
           il=indlmn(1,ilmn,itypat)
           ilpm=1+il+indlmn(2,ilmn,itypat)
           iln=indlmn(5,ilmn,itypat)
           if (iln/=iln0) indexi=indexi+2*il0+1
           klmn=j0lmn+ilmn;klmn1=cplex*klmn

           nsym_used(:)=0
           sumrho(1:cplex_eff,:)=zero
           if (noncoll.and.optrhoij==1) rotmag(:,:)=zero
           if (choice>1) rotgr(:,:,:)=zero
           if (choice>1.and.noncoll) rotmaggr(:,:,:)=zero

!          Loop over symmetries
!          --------------------
           do irot=1,nsym

             if ((symafm(irot)/=1).and.(.not.use_afm)) cycle
             kspden=ispden;if (symafm(irot)==-1) kspden=jspden
             iafm=1;if ((antiferro).and.(symafm(irot)==-1)) iafm=2
             factafm=dble(symafm(irot))

             nsym_used(iafm)=nsym_used(iafm)+1
             at_indx=min(indsym(4,irot,iatom),nrhoij)

!            Accumulate values over (mi,mj)
!            ------------------------------
             if (noncoll) summag(:,:)=zero
             if (choice>1) sumgr(:,:)=zero
             if (choice>1.and.noncoll) summaggr(:,:,:)=zero
             do mj=1,2*jl+1
               indexjj=indexj+mj;indexjj0=indexjj*(indexjj-1)/2
               do mi=1,2*il+1
                 factsym(:)=one
                 indexii=indexi+mi
                 if (indexii<=indexjj) then
                   indexk=indexjj0+indexii
                   if(cplex_eff==2.and.nspinor==2) factsym(cplex_eff)=one
                 else
                   indexk=indexii*(indexii-1)/2+indexjj
                   if(cplex_eff==2.and.nspinor==2) factsym(cplex_eff)=-one
                 end if
!                Be careful: use here R_rel^-1 in term of spherical harmonics
!                which is tR_rec in term of spherical harmonics
!                so, use transpose[zarot]
                 zarot2=pawang%zarot(mi,ilpm,il+1,irot)*pawang%zarot(mj,jlpm,jl+1,irot)
!                zarot2=pawang%zarot(ilpm,mi,il+1,irot)*pawang%zarot(jlpm,mj,jl+1,irot)

!                Rotate rhoij
                 if (optrhoij==1) then
                   if (cplex==1) then
                     sumrho(1,iafm)=sumrho(1,iafm)+zarot2*pawrhoij(at_indx)%rhoij_(indexk,kspden)
                   else
                     indexk1=2*(indexk-1)
                     sumrho(1,iafm)=sumrho(1,iafm) &
&                     +factsym(1)*zarot2*pawrhoij(at_indx)%rhoij_(indexk1+1,kspden)
                     if(cplex_eff==2) sumrho(cplex_eff,iafm)=sumrho(cplex_eff,iafm) &
&                     +factsym(cplex_eff)*factafm*zarot2*pawrhoij(at_indx)%rhoij_(indexk1+cplex_eff,kspden)
                   end if

!                  If non-collinear case, rotate rhoij magnetization
                   if (noncoll) then
                     if (cplex==1) then
                       do mu=1,3
                         summag(1,mu)=summag(1,mu)+zarot2*factafm*pawrhoij(at_indx)%rhoij_(indexk,1+mu)
                       end do
                     else
                       indexk1=2*(indexk-1)
                       do mu=1,3
!                        summag(1:cplex_eff,mu)=summag(1:cplex_eff,mu) &
!                        &             +zarot2*factsym(1:cplex_eff)*factafm*pawrhoij(at_indx)%rhoij_(indexk1+1:indexk1+cplex_eff,1+mu)
                         summag(1,mu)=summag(1,mu) &
&                         +zarot2*factsym(1)*factafm*pawrhoij(at_indx)%rhoij_(indexk1+1,1+mu)
                         if(cplex_eff==2) summag(cplex_eff,mu)=summag(cplex_eff,mu)&
&                         +zarot2*factsym(cplex_eff)*pawrhoij(at_indx)%rhoij_(indexk1+cplex_eff,1+mu)
                       end do
                     end if
                   end if
                 end if

!                Rotate gradients of rhoij
                 if (choice>1) then
                   if (cplex==1) then
                     do mu=1,ngrhoij
                       sumgr(1,mu)=sumgr(1,mu)+zarot2*tmp_grhoij(at_indx)%value(mu,indexk,kspden)
                     end do
                     if (noncoll) then
                       do mu=1,3
                         do nu=1,ngrhoij
                           summaggr(1,nu,mu)=summaggr(1,nu,mu)+zarot2*factafm*tmp_grhoij(at_indx)%value(nu,indexk,1+mu)
                         end do
                       end do
                     end if
                   else
                     indexk1=2*(indexk-1)
                     do mu=1,ngrhoij
                       sumgr(1:cplex_eff,mu)=sumgr(1:cplex_eff,mu) &
&                       +zarot2*factsym(1:cplex_eff)*tmp_grhoij(at_indx)%value(mu,indexk1+1:indexk1+cplex_eff,kspden)
                     end do
                     if (noncoll) then
                       do mu=1,3
                         do nu=1,ngrhoij
                           summaggr(1:cplex_eff,nu,mu)=summaggr(1:cplex_eff,nu,mu) &
&                           +zarot2*factsym(1:cplex_eff)*factafm*tmp_grhoij(at_indx)%value(nu,indexk1+1:indexk1+cplex_eff,1+mu)
                         end do
                       end do
                     end if
                   end if
                 end if

               end do
             end do

!            Rotate vector fields in real space (forces, magnetization, etc...)
!            Should use symrel^1 but use transpose[symrec] instead
!            ---------------------------------
!            ===== Rhoij magnetization ====
             if (noncoll.and.optrhoij==1) then
               do nu=1,3
                 do mu=1,3
                   rotmag(1:cplex_eff,mu)=rotmag(1:cplex_eff,mu)+symrec_cart(mu,nu,irot)*summag(1:cplex_eff,nu)
                 end do
               end do
             end if
!            ===== Derivatives vs atomic positions ====
             if (choice==2.or.choice==23.or.choice==24) then
               do nu=1,3
                 nushift=nu+ishift2
                 do mu=1,3
                   mushift=mu+ishift2
                   rotgr(1:cplex_eff,mushift,iafm)=&
&                   rotgr(1:cplex_eff,mushift,iafm)+dble(symrec(mu,nu,irot))*sumgr(1:cplex_eff,nushift)
                 end do
               end do
               if (noncoll) then
                 do mub=1,3 ! Loop on magnetization components
                   do mua=1,3 ! Loop on gradients
                     mushift=mua+ishift2
                     sum1(:)=zero;xsym(1:3)=dble(symrec(mua,1:3,irot))
                     do nu=1,3
                       syma=symrec_cart(mub,nu,irot)
                       sum1(1:cplex_eff)=sum1(1:cplex_eff)+syma*(summaggr(1:cplex_eff,ishift2+1,nu)*xsym(1) &
&                       +summaggr(1:cplex_eff,ishift2+2,nu)*xsym(2) &
&                       +summaggr(1:cplex_eff,ishift2+3,nu)*xsym(3))
                     end do
                     rotmaggr(1:cplex_eff,mushift,mub)=rotmaggr(1:cplex_eff,mushift,mub)+sum1(1:cplex_eff)
                   end do
                 end do
               end if
             end if
!            ===== Derivatives vs strain ====
             if (choice==3.or.choice==23) then
               work1(1:cplex_eff,1,1)=sumgr(1:cplex_eff,1+ishift3);work1(1:cplex_eff,2,2)=sumgr(1:cplex_eff,2+ishift3)
               work1(1:cplex_eff,3,3)=sumgr(1:cplex_eff,3+ishift3);work1(1:cplex_eff,2,3)=sumgr(1:cplex_eff,4+ishift3)
               work1(1:cplex_eff,1,3)=sumgr(1:cplex_eff,5+ishift3);work1(1:cplex_eff,1,2)=sumgr(1:cplex_eff,6+ishift3)
               work1(1:cplex_eff,3,1)=work1(1:cplex_eff,1,3);work1(1:cplex_eff,3,2)=work1(1:cplex_eff,2,3)
               work1(1:cplex_eff,2,1)=work1(1:cplex_eff,1,2)
               do mu=1,6
                 mushift=mu+ishift3
                 mua=alpha(mu);mub=beta(mu)
                 sum1(:)=zero;xsym(1:3)=dble(symrec(mub,1:3,irot))
                 do nu=1,3
                   syma=dble(symrec(mua,nu,irot))
                   sum1(1:cplex_eff)=sum1(1:cplex_eff)+syma*(work1(1:cplex_eff,nu,1)*xsym(1) &
&                   +work1(1:cplex_eff,nu,2)*xsym(2) &
&                   +work1(1:cplex_eff,nu,3)*xsym(3))
                 end do
                 rotgr(1:cplex_eff,mushift,iafm)=rotgr(1:cplex_eff,mushift,iafm)+sum1(1:cplex_eff)
               end do
             end if
!            ===== Second derivatives vs atomic positions ====
             if (choice==4.or.choice==24) then
               work1(1:cplex_eff,1,1)=sumgr(1:cplex_eff,1+ishift4);work1(1:cplex_eff,2,2)=sumgr(1:cplex_eff,2+ishift4)
               work1(1:cplex_eff,3,3)=sumgr(1:cplex_eff,3+ishift4);work1(1:cplex_eff,2,3)=sumgr(1:cplex_eff,4+ishift4)
               work1(1:cplex_eff,1,3)=sumgr(1:cplex_eff,5+ishift4);work1(1:cplex_eff,1,2)=sumgr(1:cplex_eff,6+ishift4)
               work1(1:cplex_eff,3,1)=work1(1:cplex_eff,1,3);work1(1:cplex_eff,3,2)=work1(1:cplex_eff,2,3)
               work1(1:cplex_eff,2,1)=work1(1:cplex_eff,1,2)
               do mu=1,6
                 mushift=mu+ishift4
                 mua=alpha(mu);mub=beta(mu)
                 sum1(:)=zero
                 xsym(1:3)=dble(symrec(mub,1:3,irot))
                 do nu=1,3
                   syma=dble(symrec(mua,nu,irot))
                   sum1(1:cplex_eff)=sum1(1:cplex_eff)+syma*(work1(1:cplex_eff,nu,1)*xsym(1) &
&                   +work1(1:cplex_eff,nu,2)*xsym(2) &
&                   +work1(1:cplex_eff,nu,3)*xsym(3))
                 end do
                 rotgr(1:cplex_eff,mushift,iafm)=rotgr(1:cplex_eff,mushift,iafm)+sum1(1:cplex_eff)
               end do
             end if

           end do ! End loop over symmetries

!          Store average result (over symmetries)
!          --------------------------------------
           if (optrhoij==1) then

!            Mean value for rhoij
             if (cplex==1) then
               ro(1)=sumrho(1,1)/nsym_used(1)
               if (abs(ro(1))>tol10) then
                 pawrhoij(iatm)%rhoijp(klmn,ispden)=ro(1)
                 if (use_res) pawrhoij(iatm)%rhoijres(klmn,ispden)=pawrhoij(iatm)%rhoijres(klmn,ispden)+ro(1)
               else
                 pawrhoij(iatm)%rhoijp(klmn,ispden)=zero
               end if
             else
               ro(1)=sumrho(1,1)/nsym_used(1)
               if (cplex_eff==2) then
                 ro(2)=sumrho(2,1)/nsym_used(1)
               else
                 ro(2)=pawrhoij(iatm)%rhoij_(klmn1,ispden)
               end if
               if (any(abs(ro(1:2))>tol10)) then
                 pawrhoij(iatm)%rhoijp(klmn1-1,ispden)=ro(1)
                 pawrhoij(iatm)%rhoijp(klmn1  ,ispden)=ro(2)
                 if (use_res) then
                   pawrhoij(iatm)%rhoijres(klmn1-1,ispden)=pawrhoij(iatm)%rhoijres(klmn1-1,ispden)+ro(1)
                   pawrhoij(iatm)%rhoijres(klmn1  ,ispden)=pawrhoij(iatm)%rhoijres(klmn1  ,ispden)+ro(2)
                 end if
               else
                 pawrhoij(iatm)%rhoijp(klmn1-1,ispden)=zero
                 pawrhoij(iatm)%rhoijp(klmn1  ,ispden)=zero
               end if
             end if

!            Non-collinear case: mean value for rhoij magnetization
             if (noncoll) then
!              Select on-zero elements
               if (cplex==1) then
                 do mu=2,4
                   ro(1)=rotmag(1,mu-1)/nsym_used(1)
                   if (abs(ro(1))>tol10) then
                     pawrhoij(iatm)%rhoijp(klmn,mu)=ro(1)
                     if (use_res) pawrhoij(iatm)%rhoijres(klmn,mu)=pawrhoij(iatm)%rhoijres(klmn,mu)+ro(1)
                   else
                     pawrhoij(iatm)%rhoijp(klmn,mu)=zero
                   end if
                 end do
               else
                 do mu=2,4
                   ro(1)=rotmag(1,mu-1)/nsym_used(1)
                   if (cplex_eff==2) then
                     ro(2)=rotmag(2,mu-1)/nsym_used(1)
                   else
                     ro(2)=pawrhoij(iatm)%rhoij_(klmn1,mu)
                   end if
                   if (any(abs(ro(1:2))>tol10)) then
                     pawrhoij(iatm)%rhoijp(klmn1-1,mu)=ro(1)
                     pawrhoij(iatm)%rhoijp(klmn1  ,mu)=ro(2)
                     if (use_res) then
                       pawrhoij(iatm)%rhoijres(klmn1-1,mu)=pawrhoij(iatm)%rhoijres(klmn1-1,mu)+ro(1)
                       pawrhoij(iatm)%rhoijres(klmn1  ,mu)=pawrhoij(iatm)%rhoijres(klmn1  ,mu)+ro(2)
                     end if
                   else
                     pawrhoij(iatm)%rhoijp(klmn1-1,mu)=zero
                     pawrhoij(iatm)%rhoijp(klmn1  ,mu)=zero
                   end if
                 end do
               end if
             end if

!            Antiferro case: mean value for down component
             if (antiferro.and.nsym_used(2)>0) then
               if (cplex==1) then
                 ro(1)=sumrho(1,2)/nsym_used(2)
                 if (abs(ro(1))>tol10) then
                   pawrhoij(iatm)%rhoijp(klmn,2)=ro(1)
                   if (use_res) pawrhoij(iatm)%rhoijres(klmn,2)=pawrhoij(iatm)%rhoijres(klmn,2)+ro(1)
                 else
                   pawrhoij(iatm)%rhoijp(klmn,2)=zero
                 end if
               else
                 ro(1:cplex_eff)=sumrho(1:cplex_eff,2)/nsym_used(2)
                 if (any(abs(ro(1:2))>tol10)) then
                   pawrhoij(iatm)%rhoijp(klmn1-1,2)=ro(1)
                   pawrhoij(iatm)%rhoijp(klmn1  ,2)=ro(2)
                   if (use_res) then
                     pawrhoij(iatm)%rhoijres(klmn1-1,2)=pawrhoij(iatm)%rhoijres(klmn1-1,2)+ro(1)
                     pawrhoij(iatm)%rhoijres(klmn1  ,2)=pawrhoij(iatm)%rhoijres(klmn1  ,2)+ro(2)
                   end if
                 else
                   pawrhoij(iatm)%rhoijp(klmn1-1,2)=zero
                   pawrhoij(iatm)%rhoijp(klmn1  ,2)=zero
                 end if
               end if
             end if

!            Select non-zero elements of rhoij
             if (ispden==pawrhoij(iatm)%nsppol) then
               if (cplex==1) then
                 if (any(abs(pawrhoij(iatm)%rhoijp(klmn,:))>tol10)) then
                   nselect=nselect+1
                   pawrhoij(iatm)%rhoijselect(nselect)=klmn
                   do jj=1,pawrhoij(iatm)%nspden
                     pawrhoij(iatm)%rhoijp(nselect,jj)=pawrhoij(iatm)%rhoijp(klmn,jj)
                   end do
                 end if
               else
                 if (any(abs(pawrhoij(iatm)%rhoijp(klmn1-1:klmn1,:))>tol10)) then
                   nselect=nselect+1;nselect1=2*nselect
                   pawrhoij(iatm)%rhoijselect(nselect)=klmn
                   do jj=1,pawrhoij(iatm)%nspden
                     pawrhoij(iatm)%rhoijp(nselect1-1,jj)=pawrhoij(iatm)%rhoijp(klmn1-1,jj)
                     pawrhoij(iatm)%rhoijp(nselect1  ,jj)=pawrhoij(iatm)%rhoijp(klmn1  ,jj)
                   end do
                 end if
               end if
             end if

           end if ! optrhoij==1

!          Store average result (over symmetries) for gradients
           if (choice>1) then
             do iplex=1,cplex_eff
               do mu=1,ngrhoij
                 pawrhoij(iatm)%grhoij(mu,klmn1+iplex-cplex,ispden)=rotgr(iplex,mu,1)/nsym_used(1)
               end do
             end do
             if (antiferro.and.nsym_used(2)>0) then
               do iplex=1,cplex_eff
                 do mu=1,ngrhoij
                   pawrhoij(iatm)%grhoij(mu,klmn1+iplex-cplex,2)=rotgr(iplex,mu,2)/nsym_used(2)
                 end do
               end do
             end if
             if (noncoll) then
               do nu=1,3
                 do iplex=1,cplex_eff
                   do mu=1,ngrhoij
                     pawrhoij(iatm)%grhoij(mu,klmn1+iplex-cplex,1+nu)=rotmaggr(iplex,mu,nu)/nsym_used(1)
                   end do
                 end do
               end do
             end if
           end if

           il0=il;iln0=iln  ! End loops over (il,im) and (jl,jm)
         end do
         jl0=jl;jln0=jln
       end do

     end do  ! End loop over ispden

!    Store number of non-zero values of rhoij
     if (optrhoij==1) pawrhoij(iatm)%nrhoijsel=nselect

   end do ! End loop over iatm

   if (noncoll.and.optrhoij==1)  then
     ABI_DEALLOCATE(summag)
     ABI_DEALLOCATE(rotmag)
   end if
   ABI_DEALLOCATE(factsym)
   if (noncoll)  then
     ABI_DEALLOCATE(symrec_cart)
   end if
   if (choice>1) then
     do iatm=1,nrhoij
       ABI_DEALLOCATE(tmp_grhoij(iatm)%value)
     end do
     ABI_DEALLOCATE(tmp_grhoij)
     ABI_DEALLOCATE(sumgr)
     ABI_DEALLOCATE(rotgr)
     if (choice>2)  then
       ABI_DEALLOCATE(work1)
     end if
     if (noncoll)  then
       ABI_DEALLOCATE(summaggr)
       ABI_DEALLOCATE(rotmaggr)
     end if
   end if

 else  ! nsym>1

!  *********************************************************************
!  If nsym==1, only copy rhoij_ into rhoij
!  also has to fill rhoijselect array

   if(pawrhoij(1)%nspden==2.and.pawrhoij(1)%nsppol==1) then
     message=' In the antiferromagnetic case, nsym cannot be 1'
     MSG_BUG(message)
   end if
   if (optrhoij==1) then
     do iatm=1,nrhoij
       cplex=pawrhoij(iatm)%cplex
       use_res=(pawrhoij(iatm)%use_rhoijres>0)
       if (use_res) then
         pawrhoij(iatm)%rhoijres(:,:)=zero
         if (cplex==1) then
           do ispden=1,pawrhoij(iatm)%nspden
             do irhoij=1,pawrhoij(iatm)%nrhoijsel
               klmn=pawrhoij(iatm)%rhoijselect(irhoij)
               pawrhoij(iatm)%rhoijres(klmn,ispden)=-pawrhoij(iatm)%rhoijp(irhoij,ispden)
             end do
           end do
         else
           do ispden=1,pawrhoij(iatm)%nspden
             do irhoij=1,pawrhoij(iatm)%nrhoijsel
               klmn1=2*pawrhoij(iatm)%rhoijselect(irhoij);jrhoij=2*irhoij
               pawrhoij(iatm)%rhoijres(klmn1-1,ispden)=-pawrhoij(iatm)%rhoijp(jrhoij-1,ispden)
               pawrhoij(iatm)%rhoijres(klmn1  ,ispden)=-pawrhoij(iatm)%rhoijp(jrhoij  ,ispden)
             end do
           end do
         end if
       end if
       nselect=0
       if (cplex==1) then
         do klmn=1,pawrhoij(iatm)%lmn2_size
           if (any(abs(pawrhoij(iatm)%rhoij_(klmn,:))>tol10)) then
             nselect=nselect+1
             pawrhoij(iatm)%rhoijselect(nselect)=klmn
             do jj=1,pawrhoij(iatm)%nspden
               ro(1)=pawrhoij(iatm)%rhoij_(klmn,jj)
               pawrhoij(iatm)%rhoijp(nselect,jj)=ro(1)
               if (use_res) pawrhoij(iatm)%rhoijres(klmn,jj)=pawrhoij(iatm)%rhoijres(klmn,jj)+ro(1)
             end do
           end if
         end do
       else
         do klmn=1,pawrhoij(iatm)%lmn2_size
           klmn1=2*klmn
           if (any(abs(pawrhoij(iatm)%rhoij_(klmn1-1:klmn1,:))>tol10)) then
             nselect=nselect+1;nselect1=2*nselect
             pawrhoij(iatm)%rhoijselect(nselect)=klmn
             do jj=1,pawrhoij(iatm)%nspden
               ro(1)=pawrhoij(iatm)%rhoij_(klmn1-1,jj)
               ro(2)=pawrhoij(iatm)%rhoij_(klmn1  ,jj)
               pawrhoij(iatm)%rhoijp(nselect1-1,jj)=ro(1)
               pawrhoij(iatm)%rhoijp(nselect1  ,jj)=ro(2)
               if (use_res) then
                 pawrhoij(iatm)%rhoijres(klmn1-1,jj)=pawrhoij(iatm)%rhoijres(klmn1-1,jj)+ro(1)
                 pawrhoij(iatm)%rhoijres(klmn1  ,jj)=pawrhoij(iatm)%rhoijres(klmn1  ,jj)+ro(2)
               end if
             end do
           end if
         end do
       end if
       pawrhoij(iatm)%nrhoijsel=nselect
     end do
   end if

 end if

!*********************************************************************
!Printing of Rhoij

 if (optrhoij==1.and.pawprtvol/=-10001) then
   pertstrg="RHOIJ";if (ipert>0) pertstrg="RHOIJ(1)"
   natinc=1;if(nrhoij>1.and.pawprtvol>=0) natinc=nrhoij-1
   do iatm=1,nrhoij,natinc
     iatom=iatm;if (nrhoij==1.and.ipert>0.and.ipert<=natom) iatom=ipert
     idum1=2;if (pawrhoij(iatm)%cplex==2.and.pawrhoij(iatm)%nspinor==1) idum1=1
     if (abs(pawprtvol)>=1) then
       write(message, '(6a,i3,a)') ch10," PAW TEST:",ch10,&
&       ' ====== Values of ',trim(pertstrg),' in symrhoij (iatom=',iatom,') ======'
       call wrtout(std_out,message,'COLL')
     end if
     do ispden=1,pawrhoij(iatm)%nspden
       if (abs(pawprtvol)>=1.and.pawrhoij(iatm)%nspden/=1) then
         write(message, '(3a)') '   Component ',trim(dspin(ispden+2*(pawrhoij(iatm)%nspden/4))),':'
       else if (pawrhoij(iatm)%nspden/=1) then
         if (pawrhoij(iatm)%nspden/=4) write(message, '(4a,i3,a,i1,a)') ch10,&
&         ' *********** ',trim(pertstrg),' (atom ',iatom,', ispden=',ispden,') **********'
         if (pawrhoij(iatm)%nspden==4) write(message, '(4a,i3,3a)') ch10,&
&         ' *********** ',trim(pertstrg),' (atom ',iatom,' - ',&
&         trim(dspin(ispden+2*(pawrhoij(iatm)%nspden/4))),') **********'
       else
         write(message, '(4a,i3,a)') ch10,&
&         ' *********** ',trim(pertstrg),' (atom ',iatom,') **********'
       end if
       call wrtout(std_out,message,'COLL')
       call print_ij(pawrhoij(iatm)%rhoijp(:,ispden),&
&       pawrhoij(iatm)%nrhoijsel,&
&       pawrhoij(iatm)%cplex,&
&       pawrhoij(iatm)%lmn_size,1,-1,idum,1,pawprtvol,&
&       pawrhoij(iatm)%rhoijselect(:),&
&       10.d0*dble(3-2*(ispden+ipert)),1,&
&       opt_sym=idum1)
     end do
   end do
   message=''
   call wrtout(std_out,message,'COLL')
 end if

 DBG_EXIT("COLL")

end subroutine symrhoij
!!***
