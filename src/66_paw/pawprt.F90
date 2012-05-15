!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawprt
!! NAME
!! pawprt
!!
!! FUNCTION
!! Print out data concerning PAW formalism
!! (pseudopotential strength, augmentation occupancies...)
!! To be called at the end of the SCF cycle
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | enunit=parameter determining units of output energies
!!   | kptopt=option for the generation of k points
!!   | natom=number of atoms in cell
!!   | ntypat = number of atom types
!!   | pawprtvol= printing volume
!!   | pawspnorb=flag: 1 if spin-orbit coupling is activated
!!   | typat(natom)=type of each atom
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  indlmn(6,i,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all types of psps
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  (only printing)
!!
!! PARENTS
!!      bethe_salpeter,outscfcv,screening,sigma
!!
!! CHILDREN
!!      mat_mlms2jmj,mat_slm2ylm,print_ij,setnoccmmp,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawprt(dtset,indlmn,lmnmax,paw_ij,pawrhoij,pawtab,&
&                 electronpositron) ! optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_parameters
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype,EP_POSITRON
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawprt'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_66_paw, except_this_one => pawprt
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmnmax
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer,optional :: electronpositron
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,dtset%ntypat)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)


!Local variables-------------------------------
!scalars
 integer,parameter :: natmax=2
 integer :: iat,iatom,iatom1,im1,im2,ipositron,ispden,itypat,ll,llp,natprt
 integer :: nspden,nsppol,optsym,unt,irhoij,jrhoij
 real(dp) :: mnorm,mx,my,mz,ntot,valmx,localm
 logical :: useexexch,usepawu
 type(pawang_type):: pawang_dum
 character(len=7),parameter :: dspin1(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)
 character(len=8),parameter :: dspin2(6)=(/"up      ","down    ","dens (n)","magn (x)","magn (y)","magn (z)"/)
 character(len=9),parameter :: dspin3(6)=(/"up       ","down     ","up-up    ","down-down","Re[up-dn]","Im[up-dn]"/)
 character(len=500) :: message0,message
!arrays
 integer :: idum(1)
 integer,allocatable :: idum1(:),idum3(:,:,:),jatom(:)
 real(dp),allocatable :: rdum2(:,:),rdum4(:,:,:,:),rhoijs(:,:)
 complex(dpc),allocatable :: noccmmp_ylm(:,:,:),noccmmp_jmj(:,:),noccmmp_slm(:,:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Initializations
 natprt=natmax;if (dtset%natom==1) natprt=1
 iatom1=dtset%natom;if (dtset%pawprtvol<0) iatom1=2
 if (dtset%pawprtvol<0) natprt=dtset%natom
 ABI_ALLOCATE(jatom,(natprt))
 if (natprt==1) then
   jatom(1)=1
 else if (natprt==2) then
   jatom(1)=1;jatom(2)=dtset%natom
 else if (natprt==dtset%natom) then
   do iat=1,dtset%natom
     jatom(iat)=iat
   end do
 else
   message="  invalid value of natprt !"
   MSG_BUG(message)
 end if

 usepawu=(count(pawtab(:)%usepawu>0)>0)
 useexexch=(count(pawtab(:)%useexexch>0)>0)

 ipositron=0
 if (present(electronpositron)) then
   if (associated(electronpositron)) ipositron=electronpositron%calctype
 end if

 write(message, '(2a)' ) ch10,&
& ' ==== Results concerning PAW augmentation regions ===='
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 message=' '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

!Print out pseudopotential strength
!----------------------------------
 do unt=1,2
   if ((unt==1).and.(dtset%enunit==0.or.dtset%enunit==2)) then
     write(message, '(a)' ) ' Total pseudopotential strength Dij (hartree):'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   else if ((unt==2).and.(dtset%enunit==1.or.dtset%enunit==2)) then
     write(message, '(a)' ) ' Total pseudopotential strength Dij (eV):'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
   if (ipositron>0) then
     if (((unt==1).and.(dtset%enunit==0.or.dtset%enunit==2)).or.&
&     ((unt==2).and.(dtset%enunit==1.or.dtset%enunit==2))) then
       if (electronpositron%has_pos_ham==0) then
         write(message, '(a)' ) ' -Note: these are the electronic Dij'
       else
         write(message, '(a)' ) ' -Note: these are the positronic Dij'
       end if
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if
   end if
   if (((unt==1).and.(dtset%enunit==0.or.dtset%enunit==2)).or.&
&   ((unt==2).and.(dtset%enunit==1.or.dtset%enunit==2))) then
     do iat=1,natprt
       iatom=jatom(iat)
       nspden=paw_ij(iatom)%ndij
       optsym=2;if (paw_ij(iatom)%cplex_dij==2.and.dtset%nspinor==1) optsym=1
       do ispden=1,nspden
         valmx=100._dp;if (ispden==1) valmx=-1._dp
         message='' ; message0=''
         if (dtset%natom>1.or.nspden>1) write(message0, '(a,i3)' ) ' Atom #',iatom
         if (nspden==1) write(message, '(a)' ) trim(message0)
         if (nspden==2) write(message, '(2a,i1)' ) trim(message0),' - Spin component ',ispden
         if (nspden==4) write(message, '(3a)' )    trim(message0),' - Component ',trim(dspin1(ispden+2*(nspden/4)))
         if (dtset%natom>1.or.nspden>1) then
           call wrtout(ab_out,message,'COLL')
           call wrtout(std_out,message,'COLL')
         end if
         if (nspden/=4.or.ispden<=2) then
           call print_ij(paw_ij(iatom)%dij(:,ispden),paw_ij(iatom)%lmn2_size,paw_ij(iatom)%cplex_dij,&
&           paw_ij(iatom)%lmn_size,2,-1,idum,0,dtset%pawprtvol,idum,valmx,unt,opt_sym=optsym)
         else
           call print_ij(paw_ij(iatom)%dij(:,ispden),paw_ij(iatom)%lmn2_size,paw_ij(iatom)%cplex_dij,&
&           paw_ij(iatom)%lmn_size,2,-1,idum,0,dtset%pawprtvol,idum,valmx,unt,&
&           asym_ij=paw_ij(iatom)%dij(:,7-ispden),opt_sym=optsym)
         end if
       end do
     end do
   end if
   message=' '
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end do

!Print out SYMMETRIZED occupancies of the partial waves
!------------------------------------------------------
 write(message, '(a)' )  ' Augmentation waves occupancies Rhoij:'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 if (ipositron>0) then
   if (electronpositron%particle==EP_POSITRON) then
     write(message, '(a)' ) ' -Note: these are the electronic Rhoij'
   else
     write(message, '(a)' ) ' -Note: these are the positronic Rhoij'
   end if
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if
 if (dtset%pawspnorb>0.and.pawrhoij(1)%cplex==1.and.dtset%kptopt/=1.and.dtset%kptopt/=2) then
   write(message, '(6a)' ) ' pawprt: - WARNING:',ch10,&
&   '       Spin-orbit coupling is activated but only real part of Rhoij occupancies',ch10,&
&   '       has been computed; they could have an imaginary part (not printed here).'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if
 do iat=1,natprt
   iatom=jatom(iat)
   nspden=pawrhoij(iatom)%nspden
   optsym=2;if (pawrhoij(iatom)%cplex==2.and.dtset%nspinor==1) optsym=1
   do ispden=1,nspden
     valmx=25._dp;if (ispden==1) valmx=-1._dp
     message='' ; message0=''
     if (dtset%natom>1.or.nspden>1) write(message0, '(a,i3)' ) ' Atom #',iatom
     if (nspden==1) write(message, '(a)' ) trim(message0)
     if (nspden==2) write(message, '(2a,i1)' ) trim(message0),' - Spin component ',ispden
     if (nspden==4) write(message, '(3a)' )    trim(message0),' - Component ',dspin2(ispden+2*(nspden/4))
     if (dtset%natom>1.or.nspden>1) then
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if
     call print_ij(pawrhoij(iatom)%rhoijp(:,ispden),pawrhoij(iatom)%nrhoijsel,&
&     pawrhoij(iatom)%cplex,pawrhoij(iatom)%lmn_size,2,-1,idum,1,dtset%pawprtvol,&
&     pawrhoij(iatom)%rhoijselect(:),valmx,1,opt_sym=optsym)
   end do
 end do
 message=' '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

!PAW+U or local exact-exchange: print out +U components of occupancies
!-------------------------------------------------------------------------------
 if ((usepawu.or.useexexch).and.ipositron/=1) then
   if(useexexch) write(message, '(a)' ) &
&   ' "Local exact-exchange" part of augmentation waves occupancies Rhoij:'
   if(usepawu) write(message, '(a)' ) &
&   ' "PAW+U" part of augmentation waves occupancies Rhoij:'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   valmx=-1._dp
   do iatom=1,dtset%natom
     nspden=pawrhoij(iatom)%nspden
     itypat=dtset%typat(iatom)
     ll=-1;llp=-1
     if (pawtab(itypat)%usepawu>0) ll=pawtab(itypat)%lpawu
     if (pawtab(itypat)%useexexch>0) llp=pawtab(itypat)%lexexch
     if (ll/=llp.and.ll/=-1.and.llp/=-1) stop "pawprt: lpawu/=lexexch forbidden !"
     ll=max(ll,llp)
     if (ll>=0) then
       optsym=2;if (pawrhoij(iatom)%cplex==2.and.dtset%nspinor==1) optsym=1
       do ispden=1,nspden
         message='' ; message0=''
         write(message0, '(a,i3,a,i1,a)') ' Atom #',iatom,' - L=',ll,' ONLY'
         if (nspden==1) write(message, '(a)' ) trim(message0)
         if (nspden==2) write(message, '(2a,i1)' ) trim(message0),' - Spin component ',ispden
         if (nspden==4) write(message, '(3a)' )    trim(message0),' - Component ',dspin2(ispden+2*(nspden/4))
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
         call print_ij(pawrhoij(iatom)%rhoijp(:,ispden),pawrhoij(iatom)%nrhoijsel,&
&         pawrhoij(iatom)%cplex,pawrhoij(iatom)%lmn_size,2,ll,indlmn(1,1:pawtab(itypat)%lmn_size,itypat),&
&         1,dtset%pawprtvol,pawrhoij(iatom)%rhoijselect(:),valmx,1,opt_sym=optsym)
       end do
!      for tests and debug :
       if(paw_ij(iatom)%ndij==4.and.(3.eq.4)) then
         ABI_ALLOCATE(rhoijs,(paw_ij(iatom)%lmn2_size*paw_ij(iatom)%cplex_dij,paw_ij(iatom)%ndij))
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           jrhoij=paw_ij(iatom)%cplex_dij*(irhoij-1)+1
           rhoijs(jrhoij,1)=pawrhoij(iatom)%rhoijp(jrhoij,  1)+pawrhoij(iatom)%rhoijp(jrhoij  ,4)
           rhoijs(jrhoij+1,1)=pawrhoij(iatom)%rhoijp(jrhoij+1,1)+pawrhoij(iatom)%rhoijp(jrhoij+1,4)
           rhoijs(jrhoij,2)=pawrhoij(iatom)%rhoijp(jrhoij,  1)-pawrhoij(iatom)%rhoijp(jrhoij  ,4)
           rhoijs(jrhoij+1,2)=pawrhoij(iatom)%rhoijp(jrhoij+1,1)-pawrhoij(iatom)%rhoijp(jrhoij+1,4)
           rhoijs(jrhoij,3)=pawrhoij(iatom)%rhoijp(jrhoij  ,2)+pawrhoij(iatom)%rhoijp(jrhoij+1,3)
           rhoijs(jrhoij+1,3)=pawrhoij(iatom)%rhoijp(jrhoij+1,2)-pawrhoij(iatom)%rhoijp(jrhoij  ,3)
           rhoijs(jrhoij,4)=pawrhoij(iatom)%rhoijp(jrhoij  ,2)-pawrhoij(iatom)%rhoijp(jrhoij+1,3)
           rhoijs(jrhoij+1,4)=pawrhoij(iatom)%rhoijp(jrhoij+1,2)+pawrhoij(iatom)%rhoijp(jrhoij  ,3)
         end do
         do ispden=1,nspden
           write(message, '(3a)' )    trim(message0),' - Component ',dspin1(ispden+2*(nspden/4))
           call print_ij(rhoijs(:,ispden),pawrhoij(iatom)%nrhoijsel,&
&           pawrhoij(iatom)%cplex,pawrhoij(iatom)%lmn_size,2,ll,indlmn(1,1:pawtab(itypat)%lmn_size,itypat),&
&           1,dtset%pawprtvol,pawrhoij(iatom)%rhoijselect(:),valmx,1,opt_sym=optsym)
           write(std_out,*)  "warning,  half of the array is not correct"
         end do
         ABI_DEALLOCATE(rhoijs)
       end if
     end if
   end do
   message=' '
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!PAW+U: print out occupations for correlated orbitals
!----------------------------------------------------

 if (usepawu.and.ipositron/=1) then
   write(message, '(3a)' ) &
&   ' ---------- LDA+U DATA --------------------------------------------------- ',ch10
   call wrtout(std_out,  message,'COLL')
   call wrtout(ab_out,  message,'COLL')
   do iatom=1,dtset%natom
     itypat=dtset%typat(iatom);ll=pawtab(itypat)%lpawu
     nspden=paw_ij(iatom)%nspden
     if ((ll>=0).and.(pawtab(itypat)%usepawu>0)) then
       write(message,fmt='(a,i5,a,i4,a)') " ====== For Atom ", iatom,&
&       ", occupations for correlated orbitals. lpawu =",ll,ch10
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')
       if(nspden==2) then
         do ispden=1,nspden
           write(message,fmt='(a,i4,a,i3,a,f10.5)') " Atom", iatom,&
&           ". Occ. for lpawu and for spin",ispden," =",paw_ij(iatom)%nocctot(ispden)
           call wrtout(std_out,message,'COLL')
           call wrtout(ab_out,message,'COLL')
         end do
         localm=paw_ij(iatom)%nocctot(2)-paw_ij(iatom)%nocctot(1)
         write(message,fmt='(a,i4,a,2x,f12.6)') " => On atom",iatom,&
&         ",  local Mag. for lpawu is  ",localm
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,  message,'COLL')
       end if
       if(nspden==4) then
         ntot=paw_ij(iatom)%nocctot(1)
         mx=paw_ij(iatom)%nocctot(2)
         my=paw_ij(iatom)%nocctot(3)
         mz=paw_ij(iatom)%nocctot(4)
         mnorm=sqrt(mx*mx+my*my+mz*mz)
         write(message,'(a,i4,a,2x,e15.8)') " => On atom",iatom,", for  lpawu, local Mag. x is  ",mx
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,  message,'COLL')
         write(message,'(14x,a,2x,e15.8)') "               local Mag. y is  ",my
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,  message,'COLL')
         write(message,'(14x,a,2x,e15.8)') "               local Mag. z is  ",mz
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,  message,'COLL')
         write(message,'(14x,a,2x,e15.8)') "               norm of Mag. is  ",mnorm
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,  message,'COLL')
         write(message,fmt='(8x,a,2x,f10.5)') " (along mag axis)    occ. for majority spin is = ",half*(ntot+mnorm)
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,message,'COLL')
         write(message,fmt='(8x,a,2x,f10.5)') " (along mag axis)    occ. for minority spin is = ",half*(ntot-mnorm)
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,message,'COLL')
       end if
       write(message,'(3a)') ch10," == Occupation matrix for correlated orbitals:",ch10
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,  message,'COLL')
       do ispden=1,nspden
         if (nspden==1.and.(paw_ij(iatom)%cplex_dij==1)) write(message,fmt='(a)')   " Up component only..."
         if (nspden==2) write(message,fmt='(a,i3)')" Occupation matrix for spin",ispden
         if (nspden==4.or.(paw_ij(iatom)%cplex_dij==2))&
&         write(message,fmt='(2a)')  " Occupation matrix for component ",trim(dspin1(ispden+2*(nspden/4)))
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,  message,'COLL')
         do im1=1,ll*2+1
           if(paw_ij(iatom)%cplex_dij==1)&
&           write(message,'(12(1x,9(1x,f10.5)))') (paw_ij(iatom)%noccmmp(1,im1,im2,ispden),im2=1,ll*2+1)
           if(paw_ij(iatom)%cplex_dij==2)&
&           write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
&           (paw_ij(iatom)%noccmmp(:,im1,im2,ispden),im2=1,ll*2+1)
           call wrtout(std_out,message,'COLL')
           call wrtout(ab_out,message,'COLL')
         end do
         write(message, '(2a)' ) ch10,' '
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,message,'COLL')
       end do
!      Transformation matrices: real->complex spherical harmonics
       if(paw_ij(iatom)%ndij==4) then
         ABI_ALLOCATE(noccmmp_ylm,(2*ll+1,2*ll+1,paw_ij(iatom)%ndij))
         noccmmp_ylm=czero
         ABI_ALLOCATE(noccmmp_slm,(2*ll+1,2*ll+1,paw_ij(iatom)%ndij))
         noccmmp_slm=czero
!        go from real notation for complex noccmmp to complex notation in noccmmp_slm
         noccmmp_slm(:,:,:)=cmplx(paw_ij(iatom)%noccmmp(1,:,:,:)&
&         ,paw_ij(iatom)%noccmmp(2,:,:,:))
         call mat_slm2ylm(ll,noccmmp_slm,noccmmp_ylm,paw_ij(iatom)%ndij,1,1,dtset%pawprtvol) ! optspin=1 because up spin are first

         do ispden=1,paw_ij(iatom)%ndij
           write(message,'(3a)') ch10,&
&           "== Occupation matrix in the complex harmonics basis for component ",&
&           trim(dspin1(ispden+2*(paw_ij(iatom)%ndij/4)))
           call wrtout(std_out,message,'COLL')
           call wrtout(ab_out,message,'COLL')
           do im1=1,ll*2+1
             write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))') (noccmmp_ylm(im1,im2,ispden),im2=1,ll*2+1)
             call wrtout(std_out,message,'COLL')
             call wrtout(ab_out,message,'COLL')
           end do
         end do
         write(message,'(a)') ch10
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,message,'COLL')
         if (dtset%pawspnorb>0) then
           ABI_ALLOCATE(noccmmp_jmj,(2*(2*ll+1),2*(2*ll+1)))
           noccmmp_jmj=czero
           call mat_mlms2jmj(ll,noccmmp_ylm,noccmmp_jmj,paw_ij(iatom)%ndij,1,1,dtset%pawprtvol) !  optspin=1: up spin are first
           write(message,'(3a)') ch10,"== Occupation matrix in the J (= L-1/2, L+1/2) and M_J basis"
           call wrtout(std_out,message,'COLL')
           call wrtout(ab_out,message,'COLL')
           do im1=1,2*(ll*2+1)
             write(message,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))') (noccmmp_jmj(im1,im2),im2=1,2*(ll*2+1))
             call wrtout(std_out,message,'COLL')
             call wrtout(ab_out,message,'COLL')
           end do
           write(message,'(a)') ch10
           call wrtout(std_out,message,'COLL')
           call wrtout(ab_out,message,'COLL')
           ABI_DEALLOCATE(noccmmp_jmj)
         end if ! pawspnorb
         ABI_DEALLOCATE(noccmmp_ylm)
         ABI_DEALLOCATE(noccmmp_slm)
       end if ! ndij==4
     end if ! ((ll>=0).and.(pawtab(itypat)%usepawu>0))
   end do
 end if

!Exact exchange: print out occupations for correlated orbitals
!-------------------------------------------------------------
 if (useexexch.and.ipositron/=1) then
   write(message, '(3a)' ) &
&   ' ---------- Exact Exchange --------------------------------------------------- ',ch10
   call wrtout(ab_out,message,'COLL')
   nspden=paw_ij(1)%nspden
   nsppol=paw_ij(1)%nsppol
   do iatom=1,dtset%natom
     itypat=dtset%typat(iatom);ll=pawtab(itypat)%lexexch
     if (ll>=0.and.pawtab(itypat)%useexexch>0) then
       ABI_ALLOCATE(paw_ij(iatom)%noccmmp,(paw_ij(iatom)%cplex_dij,2*ll+1,2*ll+1,nspden))
       ABI_ALLOCATE(paw_ij(iatom)%nocctot,(nspden))
     end if
   end do
   call setnoccmmp(1,0,rdum4,0,0,idum3,dtset%natom,0,1,nsppol,0,dtset%ntypat,&
&   paw_ij,pawang_dum,dtset%pawprtvol,pawrhoij,pawtab,rdum2,idum1,1,0)
   do iatom=1,dtset%natom
     itypat=dtset%typat(iatom);ll=pawtab(itypat)%lexexch
     if ((ll>=0).and.(pawtab(itypat)%useexexch>0)) then
       write(message,fmt='(a,i5,a,i4,a)') " ====== For Atom",iatom,&
&       ", occupations for correlated orbitals. l =",ll,ch10
       call wrtout(ab_out,message,'COLL')
       do ispden=1,nspden
         if (nspden==1) write(message,fmt='(a)')   " Up component only..."
         if (nspden==2) write(message,fmt='(a,i3)')" Occupation matrix for spin",ispden
         if (nspden==4) write(message,fmt='(2a)')  " Occupation matrix for component ",trim(dspin2(ispden+2*(nspden/4)))
         call wrtout(ab_out,message,'COLL')
         do im1=1,ll*2+1
           if(paw_ij(iatom)%cplex_dij==1)&
&           write(message,'(12(1x,9(1x,f10.5)))') (paw_ij(iatom)%noccmmp(1,im1,im2,ispden),im2=1,ll*2+1)
           if(paw_ij(iatom)%cplex_dij==2)&
&           write(message,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))') (paw_ij(iatom)%noccmmp(:,im1,im2,ispden),im2=1,ll*2+1)
           call wrtout(ab_out,message,'COLL')
         end do
         write(message, '(a)' ) ' '
         call wrtout(ab_out,message,'COLL')
       end do
       ABI_DEALLOCATE(paw_ij(iatom)%noccmmp)
       ABI_DEALLOCATE(paw_ij(iatom)%nocctot)
     end if
   end do
 end if

 message=' '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 ABI_DEALLOCATE(jatom)

 DBG_EXIT("COLL")

end subroutine pawprt

!!***
