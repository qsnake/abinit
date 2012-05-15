!{\src2tex{textfont=tt}}
!!****f* ABINIT/symdij
!! NAME
!! symdij
!!
!! FUNCTION
!! Symmetrize dij quantities (PAW psp strengths)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  indlmn(6,lmnmax,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn (for each atom type)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  ipert=index of perturbation if pawrhoij is a pertubed rhoij
!!        no meaning for ground-state calculations (should be 0)
!!  lmnmax=maximum number of PAW radial wavefunctions
!!  natom=number of atoms in cell
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of types of atoms in unit cell.
!!  option_dij=choose which part of Dij has to be symetrized (which paw_ij(:)%dijxxx):
!!             0: total dij (dij)
!!             1: dij due to compensation charge (dijhat)
!!             2: dij due to +U (dijU)
!!             3: dij XC (dijxc)
!!             4: dij XC valence only (dijxc_val)
!!             5: dij spin-orbit (dijso)
!!             6: dij, RF frozen part (dijfr)
!!  paw_ij(natom)%cplex_dij=1 if dij are REAL, 2 if they are COMPLEX
!!  paw_ij(natom)%lmn_size=number of (l,m,n) elements for the paw basis
!!  paw_ij(natom)%nspden=number of spin-density components
!!  paw_ij(natom)%nsppol=number of independant spin-density components
!!  paw_ij(natom)%dij(lmn2_size,nspden)=non-symetrized paw dij quantities
!!  pawang <type(pawang_type)>=angular mesh discretization and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  rprimd(3,3)=real space primitive translations.
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!                   reciprocal space primitive translations
!!  typat(natom)=type for each atom
!!
!! SIDE EFFECTS
!!    paw_ij(natom)%dij???(cplex_dij*lmn2_size,nspden)=symetrized dij quantities as output
!!
!! PARENTS
!!      bethe_salpeter,paw_mknewh0,respfn,scfcv,scfcv3,screening,sigma,symdij
!!
!! CHILDREN
!!      symdij
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine symdij(gprimd,indlmn,indsym,ipert,lmnmax,natom,nsym,ntypat,option_dij,&
&                 paw_ij,pawang,pawprtvol,rprimd,symafm,symrec,typat)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symdij'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipert,lmnmax,natom,nsym,ntypat,option_dij,pawprtvol
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),indsym(4,nsym,natom)
 integer,intent(in) :: symafm(nsym),symrec(3,3,nsym),typat(natom)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)

!Local variables ---------------------------------------
!scalars
 integer :: at_indx,cplex,cplex_dij,iafm,iatom,il,il0,ilmn,iln,iln0,ilpm,indexi,indexii,indexj,i1,i2,i3,i4
 integer :: indexjj,indexjj0,indexk,indexkc,iplex,irot,ispden,itypat,j0lmn,jl,jl0
 integer :: jlmn,jln,jln0,jlpm,jspden,klmn,klmnc,kspden,lmn_size,mi,mj,mu,natinc,ndij0,ndij1,nu,optsym
 logical,parameter :: afm_noncoll=.true.  ! TRUE if antiferro symmetries are used with non-collinear magnetism
 logical,parameter :: lsymnew=.false.     ! TRUE if an altern. algorithm is used (with another representation)
 logical :: antiferro,noncoll,use_afm
 real(dp) :: factafm,zarot2
 character(len=6) :: pertstrg
 character(len=500) :: message
!arrays
 integer :: symm(3,3),symrel_conv(3,3),symm2(3,3)
 real(dp) :: spinrot(4)
 complex(dpc)::Rspinrot(2,2),dijt(2,2)
 complex(dpc)::dijt2(2,2)
 integer :: nsym_used(2)
 integer,allocatable :: idum(:)
 real(dp) :: sumdij(2,2)
 real(dp),allocatable :: dijnew(:,:),dijtemp(:,:),dijtmp(:,:),rotmag(:,:),summag(:,:)
 real(dp),allocatable :: symrec_cart(:,:,:),work(:,:),sumrhoso(:,:),factsym(:),dijprint(:,:)
 character(len=7),parameter :: dspin(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)
 type(coeff2_type),allocatable :: tmp_dij(:)

! *********************************************************************

 DBG_ENTER("COLL")

!Tests of compatibility:
 if ((option_dij==1.and.paw_ij(1)%has_dijhat==0).or.&
& (option_dij==2.and.paw_ij(1)%has_dijU==0).or.&
& (option_dij==3.and.paw_ij(1)%has_dijxc==0).or.&
& (option_dij==4.and.paw_ij(1)%has_dijxc_val==0).or.&
& (option_dij==5.and.paw_ij(1)%has_dijso==0).or.&
& (option_dij==6.and.paw_ij(1)%has_dijfr==0)) then
   MSG_BUG(' Incompatibilty between option_dij and allocation of Dij !')
 end if

 if(ipert==natom+3.or.ipert==natom+4) then
   MSG_ERROR(' ipert>natom not yet implemented !')
 end if

!Symetrization occurs only when nsym>1
 if (nsym>1.and.ipert/=natom+1.and.ipert/=natom+5) then

   if (pawang%nsym==0) then
     MSG_BUG(' pawang%zarot must be allocated !')
   end if
   cplex_dij=paw_ij(1)%cplex_dij
!  Antiferro case ?
   antiferro=(paw_ij(1)%nspden==2.and.paw_ij(1)%nsppol==1.and.paw_ij(1)%ndij/=4)
!  Non-collinear case
   noncoll=(paw_ij(1)%ndij==4)
   if (noncoll.and.paw_ij(1)%cplex_dij/=2) then
     message='  cplex_dij must be 2 with ndij=4 !'
     MSG_BUG(message)
   end if
   if (noncoll) then
     ABI_ALLOCATE(summag,(cplex_dij,3))
     ABI_ALLOCATE(rotmag,(cplex_dij,3))
     ABI_ALLOCATE(work,(cplex_dij,3))
     ABI_ALLOCATE(sumrhoso,(cplex_dij,4))
     ABI_ALLOCATE(symrec_cart,(3,3,nsym))
     do irot=1,nsym
       call symredcart(gprimd,rprimd,symrec_cart(:,:,irot),symrec(:,:,irot))
     end do
   end if
!  Do we use antiferro symmetries ?
   use_afm=((antiferro).or.(noncoll.and.afm_noncoll))

!  Have to make a temporary copy of dij
   ABI_ALLOCATE(tmp_dij,(natom))
   do iatom=1,natom
     ABI_ALLOCATE(tmp_dij(iatom)%value,(paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size,paw_ij(iatom)%ndij))
     ABI_ALLOCATE(dijtmp,(paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size,paw_ij(iatom)%ndij))
     if (option_dij==0) then
       dijtmp(:,:)=paw_ij(iatom)%dij(:,:)
!      If spin-orbit, substract spin-orbit contribution !This is obsolete now !
!      if (paw_ij(iatom)%has_dijso==2) dijtmp(:,:)=dijtmp(:,:)-paw_ij(iatom)%dijso(:,:)
     else if (option_dij==1) then
       dijtmp(:,:)=paw_ij(iatom)%dijhat(:,:)
     else if (option_dij==2) then
       dijtmp(:,:)=paw_ij(iatom)%dijU(:,:)
     else if (option_dij==3) then
       dijtmp(:,:)=paw_ij(iatom)%dijxc(:,:)
     else if (option_dij==4) then
       dijtmp(:,:)=paw_ij(iatom)%dijxc_val(:,:)
     else if (option_dij==5) then
       dijtmp(:,:)=paw_ij(iatom)%dijso(:,:)
     else if (option_dij==6) then
       dijtmp(:,:)=paw_ij(iatom)%dijfr(:,:)
     end if
     if (noncoll) then  ! Has to translate Dij^{alpha,beta} into (Dij, Dij magnetic field) format
       if(lsymnew) then
         tmp_dij(iatom)%value(:,:)=dijtmp(:,:)
       else
         tmp_dij(iatom)%value(:,1)=dijtmp(:,1)+dijtmp(:,2)
         tmp_dij(iatom)%value(:,2)=dijtmp(:,3)+dijtmp(:,4)
         tmp_dij(iatom)%value(:,4)=dijtmp(:,1)-dijtmp(:,2)
         do klmn=1,paw_ij(iatom)%lmn2_size
           tmp_dij(iatom)%value(2*klmn-1,3)=-dijtmp(2*klmn  ,3)+dijtmp(2*klmn  ,4)
           tmp_dij(iatom)%value(2*klmn  ,3)= dijtmp(2*klmn-1,3)-dijtmp(2*klmn-1,4)
         end do
       end if
     else
       tmp_dij(iatom)%value(:,:)=dijtmp(:,:)
     end if
     ABI_DEALLOCATE(dijtmp)
   end do

   ndij1=1
   if (antiferro) ndij1=2
   if (noncoll)   ndij1=4
   ndij0=ndij1-1
   ABI_ALLOCATE(dijnew,(cplex_dij,ndij1))
   ABI_ALLOCATE(factsym,(cplex_dij))

!  Loops over atoms and spin components
   do iatom=1,natom
     ABI_ALLOCATE(dijtemp,(paw_ij(iatom)%cplex_dij,paw_ij(iatom)%ndij))
     ABI_ALLOCATE(dijprint,(paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size,paw_ij(iatom)%ndij))
     itypat=typat(iatom)
     lmn_size=paw_ij(iatom)%lmn_size
     cplex_dij=paw_ij(iatom)%cplex_dij
     cplex=paw_ij(iatom)%cplex

     do ispden=1,paw_ij(iatom)%nsppol
       jspden=min(3-ispden,paw_ij(iatom)%nsppol)

!      Loops over (il,im) and (jl,jm)
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
           klmn=j0lmn+ilmn;klmnc=cplex_dij*(klmn-1)

           nsym_used(:)=0
           sumdij(:,:)=zero
           if (noncoll) rotmag(:,:)=zero
           if (noncoll) sumrhoso(:,:)=zero

!          Loop over symmetries
           do irot=1,nsym
             symm(:,:)=symrec(:,:,irot)
             if(lsymnew) then
               do mi=1,3
                 do mj=1,3
                   symm2(mi,mj)=symrec(mj,mi,irot)
                 end do
               end do
               call mati3inv(symm,symrel_conv)
               call getspinrot(rprimd,spinrot,symrel_conv)
!              Rspinrot(1,1)=cmplx(spinrot(1),spinrot(4))
!              Rspinrot(1,2)=cmplx(spinrot(3),spinrot(2))
!              Rspinrot(2,1)=cmplx(-spinrot(3),spinrot(2))
!              Rspinrot(2,2)=cmplx(spinrot(1),-spinrot(4))
               Rspinrot(1,1)=cmplx(spinrot(1),-spinrot(4))
               Rspinrot(1,2)=cmplx(-spinrot(3),-spinrot(2))
               Rspinrot(2,1)=cmplx(spinrot(3),-spinrot(2))
               Rspinrot(2,2)=cmplx(spinrot(1),spinrot(4))
             end if

             if ((symafm(irot)/=1).and.(.not.use_afm)) cycle
             kspden=ispden;if (symafm(irot)==-1) kspden=jspden
             iafm=1;if ((antiferro).and.(symafm(irot)==-1)) iafm=2
             factafm=dble(symafm(irot))

             nsym_used(iafm)=nsym_used(iafm)+1
             at_indx=indsym(4,irot,iatom)
             if (noncoll) summag(:,:)=zero

!            Accumulate values over (mi,mj) and symmetries
             do mj=1,2*jl+1
               indexjj=indexj+mj;indexjj0=indexjj*(indexjj-1)/2
               do mi=1,2*il+1
                 indexii=indexi+mi
                 factsym(:)=one
                 if (indexii<=indexjj) then
                   indexk=indexjj0+indexii
                   if(cplex_dij==2.and.cplex==1) factsym(cplex_dij)=one
                 else
                   indexk=indexii*(indexii-1)/2+indexjj
                   if(cplex_dij==2.and.cplex==1) factsym(cplex_dij)=-one
                 end if
                 indexkc=cplex_dij*(indexk-1)
                 if (noncoll.and.lsymnew) then
                   do iplex=1,cplex_dij
                     if(factafm>zero) then
                       dijtemp(iplex,1)=factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,1)
                       dijtemp(iplex,2)=factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,2)
                     else
                       dijtemp(iplex,1)=factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,2)
                       dijtemp(iplex,2)=factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,1)
                     end if
                     if(factsym(2)<zero) then ! to be changed if symafm
                       dijtemp(iplex,3)=factafm*factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,4)
                       dijtemp(iplex,4)=factafm*factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,3)
                     else
                       dijtemp(iplex,3)=factafm*factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,3)
                       dijtemp(iplex,4)=factafm*factsym(iplex)*tmp_dij(at_indx)%value(indexkc+iplex,4)
                     end if
                   end do
                 end if
!                Be careful: use here R_rel^-1 in term of spherical harmonics
!                which is tR_rec in term of spherical harmonics
!                so, use transpose[zarot]....  however, we use here zarot (??)
                 zarot2=pawang%zarot(mi,ilpm,il+1,irot)*pawang%zarot(mj,jlpm,jl+1,irot)
!                zarot2=pawang%zarot(ilpm,mi,il+1,irot)*pawang%zarot(jlpm,mj,jl+1,irot)

                 if((.not.noncoll).or.(.not.lsymnew)) then
                   sumdij(1,iafm)=sumdij(1,iafm)+factsym(1)*zarot2*tmp_dij(at_indx)%value(indexkc+1,kspden)
                   if(cplex_dij==2)&
&                   sumdij(cplex_dij,iafm)=sumdij(cplex_dij,iafm)+&
&                   factsym(cplex_dij)*factafm*zarot2*tmp_dij(at_indx)%value(indexkc+cplex_dij,kspden)
                 end if

                 if (noncoll.and.(.not.lsymnew)) then
                   do mu=1,3
                     summag(1,mu)=summag(1,mu)+factsym(1)*factafm*zarot2*tmp_dij(at_indx)%value(indexkc+1,1+mu)
                     if(cplex_dij==2)&
&                     summag(cplex_dij,mu)=summag(cplex_dij,mu)+&
&                     factsym(cplex_dij)*zarot2*tmp_dij(at_indx)%value(indexkc+cplex_dij,1+mu)
                   end do
                 end if
                 if (noncoll.and.(lsymnew)) then
                   dijt(1,1)=cmplx(dijtemp(1,1),dijtemp(2,1))
                   dijt(2,2)=cmplx(dijtemp(1,2),dijtemp(2,2))
                   dijt(1,2)=cmplx(dijtemp(1,3),dijtemp(2,3))
                   dijt(2,1)=cmplx(dijtemp(1,4),dijtemp(2,4))
                   dijt2(:,:)=czero
                   do i1=1,2
                     do i4=1,2
                       do i2=1,2
                         do i3=1,2
                           dijt2(i1,i4)=dijt2(i1,i4)+Rspinrot(i1,i2)*dijt(i2,i3)*conjg(Rspinrot(i4,i3))
                         end do
                       end do
                     end do
                   end do
                 end if

                 if (noncoll.and.(lsymnew)) then
                   do mu=1,4
                     if(mu==1) then
                       i1=1;i4=1
                     else if(mu==2) then
                       i1=2;i4=2
                     else if(mu==3) then
                       i1=1;i4=2
                     else if(mu==4) then
                       i1=2;i4=1
                     end if
                     sumrhoso(1,mu)=sumrhoso(1,mu)+zarot2*real(dijt2(i1,i4))
                     sumrhoso(2,mu)=sumrhoso(2,mu)+zarot2*imag(dijt2(i1,i4))
                   end do
                 end if

               end do ! mi
             end do ! mj

!            If non-collinear case, rotate Dij magnetization
             if (noncoll.and.(.not.lsymnew)) then
!              Should use symrel^1 but use transpose[symrec] instead
               do nu=1,3
                 do mu=1,3
                   rotmag(1:cplex_dij,mu)=rotmag(1:cplex_dij,mu)+symrec_cart(mu,nu,irot)*summag(1:cplex_dij,nu) !we need the transpose ?
                 end do
               end do
             end if

           end do ! End loop over symmetries

           if((.not.noncoll).or.(.not.lsymnew)) then
!            Store new value of dij
             do iplex=1,cplex_dij
               dijnew(iplex,1)=sumdij(iplex,1)/nsym_used(1)
               if (abs(dijnew(iplex,1))<=tol10) dijnew(iplex,1)=zero
             end do

!            Antiferromagnetic case: has to fill up "down" component of dij
             if (antiferro.and.nsym_used(2)>0) then
               do iplex=1,cplex_dij
                 dijnew(iplex,2)=sumdij(iplex,2)/nsym_used(2)
                 if (abs(dijnew(iplex,2))<=tol10) dijnew(iplex,2)=zero
               end do
             end if
           else if (noncoll.and.(lsymnew)) then
             do mu=1,4
               do iplex=1,cplex_dij
                 dijnew(iplex,mu)=sumrhoso(iplex,mu)/nsym_used(1)
                 if (abs(dijnew(iplex,mu))<=tol10) dijnew(iplex,mu)=zero
               end do
             end do
           end if

!          Non-collinear case: store new values of Dij magnetization
           if (noncoll.and.(.not.lsymnew)) then
!            Select on-zero elements
             do mu=1,3
               do iplex=1,cplex_dij
                 rotmag(iplex,mu)=rotmag(iplex,mu)/nsym_used(1)
                 if (abs(rotmag(iplex,mu))<=tol10) rotmag(iplex,mu)=zero
               end do
             end do
!            Transfer back to Dij^{alpha,beta}
             if(.not.lsymnew) then
               dijnew(1,1)=half*(dijnew(1,1)+rotmag(1,3))
               dijnew(2,1)=half*(dijnew(2,1)+rotmag(2,3))
               dijnew(1,2)=      dijnew(1,1)-rotmag(1,3)
               dijnew(2,2)=      dijnew(2,1)-rotmag(2,3)
               dijnew(1,3)=half*(rotmag(1,1)+rotmag(2,2))
               dijnew(2,3)=half*(rotmag(2,1)-rotmag(1,2))
               dijnew(1,4)=half*(rotmag(1,1)-rotmag(2,2))
               dijnew(2,4)=half*(rotmag(2,1)+rotmag(1,2))
             end if
           end if

!          Transfer new value of Dij in suitable pointer
           if (option_dij==0) then
             paw_ij(iatom)%dij(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==1) then
             paw_ij(iatom)%dijhat(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==2) then
             paw_ij(iatom)%dijU(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==3) then
             paw_ij(iatom)%dijxc(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==4) then
             paw_ij(iatom)%dijxc_val(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==5) then
             paw_ij(iatom)%dijso(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           else if (option_dij==6) then
             paw_ij(iatom)%dijfr(klmnc+1:klmnc+cplex_dij,ispden:ispden+ndij0)=dijnew(1:cplex_dij,1:ndij1)
           end if

           il0=il;iln0=iln  ! End loops over (il,im) and (jl,jm)
         end do
         jl0=jl;jln0=jln
       end do

       if(lsymnew.and.(abs(pawprtvol)>=3)) then
         dijprint(:,ispden:ispden+ndij0)=paw_ij(iatom)%dij(:,ispden:ispden+ndij0)
         write(message,'(2a,i4)') ch10,"Dij after sym in upup dndn updn dnup representation",ispden
         call wrtout(std_out,message,'COLL')
         do mu=1,ndij1
           write(message,'(2i4)') mu,ndij1
           call wrtout(std_out,message,'COLL')
           call print_ij(dijprint(:,mu),paw_ij(iatom)%lmn2_size,&
&           paw_ij(iatom)%cplex_dij,paw_ij(iatom)%lmn_size,1,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*ispden),1)
         end do
       end if

     end do ! ispden
     ABI_DEALLOCATE(dijtemp)
     ABI_DEALLOCATE(dijprint)
   end do ! iatom

   ABI_DEALLOCATE(dijnew)
   ABI_DEALLOCATE(factsym)
   if (noncoll)  then
     ABI_DEALLOCATE(summag)
     ABI_DEALLOCATE(rotmag)
     ABI_DEALLOCATE(symrec_cart)
     ABI_DEALLOCATE(work)
     ABI_DEALLOCATE(sumrhoso)
   end if
   do iatom=1,natom
     ABI_DEALLOCATE(tmp_dij(iatom)%value)
   end do
   ABI_DEALLOCATE(tmp_dij)

!  If spin-orbit, add again spin-orbit contribution  (Obsolete !!!)
!  if (option_dij==0.and.paw_ij(1)%has_dijso==2) then
!  do iatom=1,natom
!  paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)+paw_ij(iatom)%dijso(:,:)
!  end do
!  end if

 else if (ipert/=natom+1.and.ipert/=natom+5) then  ! nsym>1

!  *********************************************************************
!  If nsym==1, only cut small components of dij

   if(paw_ij(1)%nspden==2.and.paw_ij(1)%nsppol==1) then
     MSG_BUG(' In the antiferromagnetic case, nsym cannot be 1')
   end if
   do iatom=1,natom
     do ispden=1,paw_ij(iatom)%ndij

       if (option_dij==0) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dij(klmn,ispden))<=tol10) paw_ij(iatom)%dij(klmn,ispden)=zero
         end do
       else if (option_dij==1) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijhat(klmn,ispden))<=tol10) paw_ij(iatom)%dijhat(klmn,ispden)=zero
         end do
       else if (option_dij==2) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijU(klmn,ispden))<=tol10) paw_ij(iatom)%dijU(klmn,ispden)=zero
         end do
       else if (option_dij==3) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijxc(klmn,ispden))<=tol10) paw_ij(iatom)%dijxc(klmn,ispden)=zero
         end do
       else if (option_dij==4) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijxc_val(klmn,ispden))<=tol10) paw_ij(iatom)%dijxc_val(klmn,ispden)=zero
         end do
       else if (option_dij==5) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijso(klmn,ispden))<=tol10) paw_ij(iatom)%dijso(klmn,ispden)=zero
         end do
       else if (option_dij==6) then
         do klmn=1,paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij
           if (abs(paw_ij(iatom)%dijfr(klmn,ispden))<=tol10) paw_ij(iatom)%dijfr(klmn,ispden)=zero
         end do
       end if

     end do
   end do

 end if  ! nsym>1

!*********************************************************************
!Printing of Dij

 if (abs(pawprtvol)>=1.and.option_dij==0.and.ipert/=natom+1.and.ipert/=natom+5) then
   pertstrg="DIJ";if (ipert>0) pertstrg="DIJ(1)"
   natinc=1;if(natom>1.and.pawprtvol>=0) natinc=natom-1
   do iatom=1,natom,natinc
     write(message, '(6a,i3,a)') ch10," PAW TEST:",ch10,&
&     ' ====== Values of ',trim(pertstrg),' in symdij (iatom=',iatom,') (Hartree) ======'
     call wrtout(std_out,message,'COLL')
     optsym=2;if (paw_ij(iatom)%cplex_dij==2.and.ipert>0) optsym=1
     do ispden=1,paw_ij(iatom)%ndij
       if (paw_ij(iatom)%ndij==1) then
         write(message, '(4a,i3,a)') ch10,&
&         ' *********** ',trim(pertstrg),' (atom ',iatom,') **********'
       else
         write(message, '(4a,i3,3a)') ch10,&
&         ' *********** ',trim(pertstrg),' (atom ',iatom,', Component ', &
&         trim(dspin(ispden+2*(paw_ij(iatom)%ndij/4))),') **********'
       end if
       call wrtout(std_out,message,'COLL')
       if (paw_ij(iatom)%ndij/=4.or.ispden<=2) then
         call print_ij(paw_ij(iatom)%dij(:,ispden),paw_ij(iatom)%lmn2_size,&
&         paw_ij(iatom)%cplex_dij,paw_ij(iatom)%lmn_size,1,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*ispden),1,&
&         opt_sym=optsym)
       else
         if (ipert==0) then
           call print_ij(paw_ij(iatom)%dij(:,ispden),paw_ij(iatom)%lmn2_size,&
&           paw_ij(iatom)%cplex_dij,paw_ij(iatom)%lmn_size,1,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*ispden),1,&
&           asym_ij=paw_ij(iatom)%dij(:,7-ispden))
         else
           call print_ij(paw_ij(iatom)%dij(:,ispden),paw_ij(iatom)%lmn2_size,&
&           paw_ij(iatom)%cplex_dij,paw_ij(iatom)%lmn_size,1,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*ispden),1,&
&           opt_sym=optsym)
         end if
       end if
     end do
   end do
   call wrtout(std_out,"",'COLL')
 end if

 DBG_EXIT("COLL")

end subroutine symdij
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/symdij_all
!! NAME
!! symdij_all
!!
!! FUNCTION
!! Symmetrize all the dij quantities that have been computed. 
!!
!! INPUTS
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  indlmn(6,lmnmax,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn (for each atom type)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  ipert=index of perturbation if pawrhoij is a pertubed rhoij
!!        no meaning for ground-state calculations (should be 0)
!!  lmnmax=maximum number of PAW radial wavefunctions
!!  natom=number of atoms in cell
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of types of atoms in unit cell.
!!  paw_ij(natom)%cplex_dij=1 if dij are REAL, 2 if they are COMPLEX
!!  paw_ij(natom)%lmn_size=number of (l,m,n) elements for the paw basis
!!  paw_ij(natom)%nspden=number of spin-density components
!!  paw_ij(natom)%nsppol=number of independant spin-density components
!!  paw_ij(natom)%dij(lmn2_size,nspden)=non-symetrized paw dij quantities
!!  pawang <type(pawang_type)>=angular mesh discretization and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  rprimd(3,3)=real space primitive translations.
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!                   reciprocal space primitive translations
!!  typat(natom)=type for each atom
!!
!! SIDE EFFECTS
!!    paw_ij(natom)%dij???(cplex_dij*lmn2_size,nspden)=symetrized dij quantities as output
!!
!! PARENTS
!!      paw_mknewh0,screening,sigma
!!
!! CHILDREN
!!      symdij
!!
!! SOURCE

subroutine symdij_all(gprimd,indlmn,indsym,ipert,lmnmax,natom,nsym,ntypat,&
&                     paw_ij,pawang,pawprtvol,rprimd,symafm,symrec,typat)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symdij_all'
 use interfaces_66_paw, except_this_one => symdij_all
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipert,lmnmax,natom,nsym,ntypat,pawprtvol
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),indsym(4,nsym,natom)
 integer,intent(in) :: symafm(nsym),symrec(3,3,nsym),typat(natom)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)

!Local variables ---------------------------------------
!scalars
 integer,parameter :: MAX_NOPTS=9
 integer :: ii,option_dij,nopt
 !character(len=500) :: message
!arrays
 integer :: options(MAX_NOPTS)

! *********************************************************************

 nopt = 0
 if (ANY(paw_ij(:)%has_dij==2)) then   
   nopt = nopt + 1
   options(nopt) = 0
 end if

 if (ANY(paw_ij(:)%has_dijhat==2)) then   
   nopt = nopt + 1
   options(nopt) = 1
 end if

 if (ANY(paw_ij(:)%has_dijU==2))   then
   nopt = nopt + 1
   options(nopt) = 2
 end if

 if (ANY(paw_ij(:)%has_dijxc==2)) then 
   nopt = nopt + 1
   options(nopt) = 3
 end if

 if (ANY(paw_ij(:)%has_dijxc_val==2)) then 
   nopt = nopt + 1
   options(nopt) = 4
 end if

 if (ANY(paw_ij(:)%has_dijso==2)) then 
   nopt = nopt + 1
   options(nopt) = 5
 end if

 if (ANY(paw_ij(:)%has_dijfr==2)) then
   nopt = nopt + 1
   options(nopt) = 6
 end if

!FIXME  Dij_hartree and dij_exech_pot are not symmetrized, 

 if (ANY(paw_ij(:)%has_dijhartree==2)) then
   nopt = nopt + 1
   options(nopt) = 7
 end if

 if (ANY(paw_ij(:)%has_exexch_pot==2)) then
   nopt = nopt + 1
   options(nopt) = 8
   MSG_ERROR("Not coded")
 end if

 do ii=1,nopt
   option_dij = options(ii)
   call symdij(gprimd,indlmn,indsym,ipert,lmnmax,natom,nsym,ntypat,option_dij,&
&   paw_ij,pawang,pawprtvol,rprimd,symafm,symrec,typat)
 end do

end subroutine symdij_all
!!***
