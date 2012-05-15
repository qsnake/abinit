!{\src2tex{textfont=tt}}
!!****f* ABINIT/symrhg
!! NAME
!! symrhg
!!
!! FUNCTION
!! From rho(r), generate rho(G), symmetrize it, and
!! come back to the real space for a symmetrized rho(r).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cplex=1 if rhor is real, 2 if rhor is complex
!! gprimd(3,3)=dimensional reciprocal space primitive translations
!! irrzon(nfft,2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!! mpi_enreg=informations about MPI parallelization
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nspden=number of spin-density components
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nsym=number of symmetry elements.
!! phnons(2,nfft,(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!! rprimd(3,3)=dimensional real space primitive translations
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!
!! OUTPUT
!! rhog(2,nfft)=symmetrized rho(G) (total) electron density in G space
!!
!! SIDE EFFECTS
!! Input/Output
!! rhor(cplex*nfft,nspden)=array for electron density in electrons/bohr**3.
!! Input, but also output, if symmetrization is applied.
!! Also output if nspden > 1 (change spin components)
!!
!! NOTES
!! When using spin-polarization (nspden==2),
!! put total density in first half of rhor array and spin up in second half
!! If (nspden=2 and nsppol=2) the density is transformed as  (up,down) => (up+down,up)
!! If (nspden=2 and nsppol=1) anti-ferromagnetic symmetry operations
!!  must be used, such as to transform (2*up) => (up+down,up)
!! In spin-polarized, and if there is no symmetry to be
!! applied on the system, only the total density is generated in G space
!!
!! PARENTS
!!      crho,mkrho,mkrho3,nstpaw3,rhofermi3,suscep_dyn,suscep_kxc_dyn
!!      suscep_stat,vtorho,vtorho3,vtorhorec,vtorhotf,wfd_mkrho
!!
!! CHILDREN
!!      fourdp,leave_new,matr3inv,symredcart,timab,wrtout,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine symrhg(cplex,gprimd,irrzon,mpi_enreg,nfft,nfftot,ngfft,nspden,nsppol,nsym,paral_kgb,&
&                 phnons,rhog,rhor,rprimd,symafm,symrel)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symrhg'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nfftot,nspden,nsppol,nsym,paral_kgb
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: irrzon(nfftot**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4)),ngfft(18)
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),phnons(2,nfftot**(1-1/nsym),(nspden/nsppol)-3*(nspden/4)),rprimd(3,3)
 real(dp),intent(inout) :: rhor(cplex*nfft,nspden)
 real(dp),intent(out) :: rhog(2,nfft)

!Local variables-------------------------------
!scalars
 integer :: ier,imagn,ind,ind2,indsy,ispden,isym,iup,izone,izone_max,j,j1,j2,j3,jsym
 integer :: k1,k2,k3,l1,l2,l3,map,me_fft
 integer :: n1,n2,n3,nd2,nproc_fft,nspden_eff,nsym_used,numpt,nup,old_paral_level
 integer :: r2,rep,spaceComm
 logical,parameter :: afm_noncoll=.true.  ! TRUE if antiferro symmetries are used in non-collinear magnetism
 real(dp) :: magxsu1,magxsu2,magysu1,magysu2,magzsu1,magzsu2,mxi,mxr,myi,myr,mzi,mzr,phi,phr,rhosu1,rhosu2
 character(len=500) :: message
!arrays
 integer,allocatable :: isymg(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: magngx(:,:),magngy(:,:),magngz(:,:)
 real(dp),allocatable :: rhosu1_arr(:),rhosu2_arr(:),work(:)
 real(dp),allocatable :: symafm_used(:),symrec_cart(:,:,:),symrel_cart(:,:,:)
!No abirules
!Statment function
 map(j1,n1)=mod(n1+mod(j1,n1),n1)
!integer :: ifft,ipw

!*************************************************************************
!
!Note the timing channel 17 excludes the
!different Fourier transforms

 ABI_ALLOCATE(work,(cplex*nfft))

!Special treatment for spin-polarized case
 if(nspden==2 .and. nsppol==2) then
   call timab(17,1,tsec)
!  When nspden=2 and nsppol=2, put total density in first half
!  of rhor array and spin up in second half  (up,down) => (up+down,up)
   call timab(17,1,tsec)
   work(:)=rhor(:,1)               ! up => work
   rhor(:,1)=rhor(:,1)+rhor(:,2)   ! up+down
   rhor(:,2)=work(:)               ! work => up
   call timab(17,2,tsec)
 end if

!Special treatment for antiferromagnetism case
 if(nspden==2 .and. nsppol==1) then
   call timab(17,1,tsec)
!  When nspden=2 and nsppol=1, (2*up) => (2*up,up)
!  Indeed, what was delivered to the present routine is a "total" density,
!  obtained from occupation numbers varying between 0 and 2,
!  but for spin up only potential.
   rhor(:,2)=half*rhor(:,1)
   call timab(17,2,tsec)
 end if

!Special treatment for non-collinear magnetism case
 if(nspden==4) then
   call timab(17,1,tsec)
   rhor(:,1)=rhor(:,1)+rhor(:,4)     !nup+ndown
   rhor(:,2)=rhor(:,2)-rhor(:,1)     !mx (n+mx-n)
   rhor(:,3)=rhor(:,3)-rhor(:,1)     !my (n+my-n)
   rhor(:,4)=rhor(:,1)-two*rhor(:,4) !mz=n-2ndown
   call timab(17,2,tsec)
 end if

!DEBUG
!write(std_out,*)'  ispden  ifft   rhor(ifft,ispden)'
!do ifft=1,nfft,123
!do ispden=1,nspden
!if(cplex==1)write(std_out,*)ispden,ifft,rhor(ifft,ispden)
!if(cplex==2)write(std_out,*)ispden,ifft,rhor(2*ifft-1,ispden),rhor(2*ifft,ispden)
!end do
!if(cplex==1)write(std_out,*)3,ifft,rhor(ifft,1)-rhor(ifft,2)
!if(cplex==2)write(std_out,*)3,ifft,rhor(2*ifft-1,1)-rhor(2*ifft-1,2),&
!&  rhor(2*ifft,1)-rhor(2*ifft,2)
!end do
!write(std_out,*)' symrhg : leave'
!ENDDEBUG

 if(nsym==1)then

   if(nspden==2 .and. nsppol==1) then
!    There must be at least one anti-ferromagnetic operation
     write(message,'(a,a,a)') ' symrhg : BUG -',ch10,&
&     ' In the antiferromagnetic case, nsym cannot be 1'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  If not using symmetry, still want total density in G space rho(G).
!  Fourier transform (incl normalization) to get rho(G)
   work(:)=rhor(:,1)
   call fourdp(cplex,rhog,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
 else

!  Treat either full density, spin-up density or magnetization
!  Note the decrease of ispden to the value 1, in order to finish
!  with rhog of the total density (and not the spin-up density or magnetization)
   nspden_eff=nspden;if (nspden==4) nspden_eff=1
   do ispden=nspden_eff,1,-1

!    Prepare the density to be symmetrized, in the reciprocal space
     if(nspden==1 .or. nsppol==2 .or. (nspden==4.and.(.not.afm_noncoll)))then
       imagn=1
       nsym_used=0
       do isym=1,nsym
         if(symafm(isym)==1)nsym_used=nsym_used+1
!        DEBUG
!        write(std_out,*)' symrhg : isym,symafm(isym)',isym,symafm(isym)
!        ENDDEBUG
       end do
     else if(nspden==2 .and. nsppol==1)then   ! antiferromagnetic case
       imagn=ispden
       nsym_used=nsym/ispden
     else if (nspden==4) then
       imagn=1
       nsym_used=nsym/ispden
     end if

!    DEBUG
!    write(std_out,*)' symrhg : nsym_used=',nsym_used
!    ENDDEBUG

!    rhor -fft-> rhog    (rhog is used as work space)
!    Note : it should be possible to reuse rhog in the antiferromagnetic case
!    this would avoid one FFT
     work(:)=rhor(:,ispden)
     call fourdp(cplex,rhog,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     if (nspden==4) then
       ABI_ALLOCATE(magngx,(2,nfft))
       ABI_ALLOCATE(magngy,(2,nfft))
       ABI_ALLOCATE(magngz,(2,nfft))
       work(:)=rhor(:,2)
       call fourdp(cplex,magngx,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
       work(:)=rhor(:,3)
       call fourdp(cplex,magngy,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
       work(:)=rhor(:,4)
       call fourdp(cplex,magngz,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     end if

!    Begins the timing here only , to exclude FFTs
     call timab(17,1,tsec)

     n1=ngfft(1);n2=ngfft(2);n3=ngfft(3);nproc_fft=ngfft(10);me_fft=ngfft(11);nd2=n2/nproc_fft

!    DEBUG
!    write(std_out,*)' symrhg : fourier space density'
!    ind=0
!    do j3=1,n3
!    do j2=1,n2
!    if(((j2-1)/nd2)==me_fft) then
!    do j1=1,n1
!    ind=ind+1
!    write(std_out,'(5i4,2es16.6)') j1-1,j2-1,j3-1,ind,n1*(n2*(j3-1)+j2-1)+j1,rhog(:,ind)
!    end do
!    end if
!    end do
!    end do
!    ENDDEBUG

!    DEBUG
!    phnons(2,:,1)=zero
!    write(std_out,*)' symrhg : density before symmetrization, phnons,irrzon'
!    do ipw=1,nfft
!    j=ipw-1;j1=modulo(j,n1);r2=modulo(j/n1,nd2);j3=j/(n1*nd2);j2=me_fft*nd2+r2
!    ind=n1*(n2*j3+j2)+j1+1 !this is ind in the full array proc
!    write(std_out,'(6i4,4es16.6,2i6)' )ipw,j1,j2,j3,r2,ind,rhog(:,ipw) !,phnons(:,ipw,1),irrzon(ipw,:,1)
!    write(std_out,* )ipw,rhog(:,ipw) !,phnons(:,ipw,1),irrzon(ipw,:,1)
!    write(std_out,'8i4')j+1,j1,j2,j3,r2,ind
!    end do
!    end do
!    ENDDEBUG

!    The following is only valid for total, up or dn density
!    -------------------------------------------------------

!    Get maxvalue of izone
     izone_max=count(irrzon(:,2,imagn)>0)
     ABI_ALLOCATE(rhosu1_arr,(izone_max))
     ABI_ALLOCATE(rhosu2_arr,(izone_max))

     numpt=0
     do izone=1,nfftot

!      Get repetition number
       rep=irrzon(izone,2,imagn)
       if(rep==0)exit

!      Compute number of unique points in this symm class:
       nup=nsym_used/rep

!      Accumulate charge over equivalent points
       rhosu1=zero
       rhosu2=zero
       do iup=1,nup
         ind=irrzon(iup+numpt,1,imagn)
         j=ind-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2);r2=modulo(j2,nd2)
         if(modulo(j/n1,n2)/nd2==me_fft) then ! this ind is to be treated by me_fft
           ind=n1*(nd2*j3+r2)+j1+1 !this is ind in the current proc
           rhosu1=rhosu1+rhog(1,ind)*phnons(1,iup+numpt,imagn)&
&           -rhog(2,ind)*phnons(2,iup+numpt,imagn)
           rhosu2=rhosu2+rhog(2,ind)*phnons(1,iup+numpt,imagn)&
&           +rhog(1,ind)*phnons(2,iup+numpt,imagn)
         end if

       end do
       rhosu1=rhosu1/dble(nup)
       rhosu2=rhosu2/dble(nup)
       rhosu1_arr(izone)=rhosu1
       rhosu2_arr(izone)=rhosu2
!      Keep index of how many points have been considered:
       numpt=numpt+nup

!      End loop over izone
     end do

!    Reduction in case of FFT parallelization
     if(mpi_enreg%mode_para=='b')then
       old_paral_level=mpi_enreg%paral_level
       mpi_enreg%paral_level=3
       spaceComm=mpi_enreg%comm_fft
       call xsum_mpi(rhosu1_arr,spaceComm,ier)
       call xsum_mpi(rhosu2_arr,spaceComm,ier)
       mpi_enreg%paral_level=old_paral_level
     end if

!    Now symmetrize the density
     numpt=0
     do izone=1,nfftot

!      Get repetition number
       rep=irrzon(izone,2,imagn)
       if(rep==0)exit

!      Compute number of unique points in this symm class:
       nup=nsym_used/rep

!      Define symmetrized rho(G) at equivalent points:
       do iup=1,nup
         ind=irrzon(iup+numpt,1,imagn)
!        decompose ind-1=n1(n2 j3+ j2)+j1
         j=ind-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2);r2=modulo(j2,nd2)
         if(modulo(j/n1,n2)/nd2==me_fft) then ! this ind is to be treated by me_fft
!          ind in the proc ind-1=n1(nd2 j3+ r2)+j1
           ind=n1*(nd2*j3+r2)+j1+1 !this is ind in the current proc
           rhog(1,ind)=rhosu1_arr(izone)*phnons(1,iup+numpt,imagn)&
&           +rhosu2_arr(izone)*phnons(2,iup+numpt,imagn)
           rhog(2,ind)=rhosu2_arr(izone)*phnons(1,iup+numpt,imagn)&
&           -rhosu1_arr(izone)*phnons(2,iup+numpt,imagn)
         end if
       end do

!      Keep index of how many points have been considered:
       numpt=numpt+nup

!      End loop over izone
     end do
     ABI_DEALLOCATE(rhosu1_arr)
     ABI_DEALLOCATE(rhosu2_arr)

!    The following is only valid for magnetization
!    ---------------------------------------------
     if (nspden==4) then

!      Transfer symmetries in cartesian coordinates
!      Compute symmetries in reciprocal space in cartesian coordinates
       ABI_ALLOCATE(symrec_cart,(3,3,nsym_used))
       ABI_ALLOCATE(symrel_cart,(3,3,nsym_used))
       ABI_ALLOCATE(symafm_used,(nsym_used))
       jsym=0
       do isym=1,nsym
         if (symafm(isym)/=1.and.(.not.afm_noncoll)) cycle
         jsym=jsym+1
         symafm_used(jsym)=dble(symafm(isym))
         call symredcart(rprimd,gprimd,symrel_cart(:,:,jsym),symrel(:,:,isym))
         call matr3inv(symrel_cart(:,:,jsym),symrec_cart(:,:,jsym))
       end do

       numpt=count(irrzon(:,1,imagn)>0)
       ABI_ALLOCATE(isymg,(numpt))
       isymg=0
       ABI_ALLOCATE(rhosu1_arr,(3*izone_max))
       ABI_ALLOCATE(rhosu2_arr,(3*izone_max))

!      Accumulate magnetization over equivalent points
!      Use all symmetries (not only those linking different g points)
!      Use Inverse[Transpose[symrel]]=symrec
       numpt=0
       do izone=1,izone_max
         magxsu1=zero;magxsu2=zero
         magysu1=zero;magysu2=zero
         magzsu1=zero;magzsu2=zero
         ind=irrzon(1+numpt,1,1)
         rep=irrzon(izone,2,1)
         nup=nsym_used/rep
         j=ind-1;l1=modulo(j,n1);l2=modulo(j/n1,n2);l3=j/(n1*n2)
         jsym=0
         do isym=1,nsym
           if (symafm(isym)/=1.and.(.not.afm_noncoll)) cycle
           jsym=jsym+1
           j1=symrel(1,1,isym)*l1+symrel(2,1,isym)*l2+symrel(3,1,isym)*l3
           j2=symrel(1,2,isym)*l1+symrel(2,2,isym)*l2+symrel(3,2,isym)*l3
           j3=symrel(1,3,isym)*l1+symrel(2,3,isym)*l2+symrel(3,3,isym)*l3
           k1=map(j1,n1);k2=map(j2,n2);k3=map(j3,n3)
           indsy=1+k1+n1*(k2+n2*k3)
           ind2=-1;iup=numpt
           do while (ind2/=indsy.and.iup<numpt+nup)
             iup=iup+1;ind2=irrzon(iup,1,1)
           end do
           if (ind2/=indsy) stop "ERROR (1) in symrhg !"
           if (isymg(iup)==0) isymg(iup)=jsym
           if(modulo((indsy-1)/n1,n2)/nd2==me_fft) then  ! this is indsy is to be treated by me_fft
             indsy=n1*(nd2*k3+modulo(k2,nd2))+k1+1        ! this is indsy in the current proc
             phr=phnons(1,iup,imagn);if (rep==1) phr=phr*symafm_used(jsym) !if rep==2, symafm is already included in phnons
             phi=phnons(2,iup,imagn);if (rep==1) phi=phi*symafm_used(jsym) !(see irrzg.F90)
             mxr=symrel_cart(1,1,jsym)*magngx(1,indsy)+symrel_cart(1,2,jsym)*magngy(1,indsy)+symrel_cart(1,3,jsym)*magngz(1,indsy)
             mxi=symrel_cart(1,1,jsym)*magngx(2,indsy)+symrel_cart(1,2,jsym)*magngy(2,indsy)+symrel_cart(1,3,jsym)*magngz(2,indsy)
             myr=symrel_cart(2,1,jsym)*magngx(1,indsy)+symrel_cart(2,2,jsym)*magngy(1,indsy)+symrel_cart(2,3,jsym)*magngz(1,indsy)
             myi=symrel_cart(2,1,jsym)*magngx(2,indsy)+symrel_cart(2,2,jsym)*magngy(2,indsy)+symrel_cart(2,3,jsym)*magngz(2,indsy)
             mzr=symrel_cart(3,1,jsym)*magngx(1,indsy)+symrel_cart(3,2,jsym)*magngy(1,indsy)+symrel_cart(3,3,jsym)*magngz(1,indsy)
             mzi=symrel_cart(3,1,jsym)*magngx(2,indsy)+symrel_cart(3,2,jsym)*magngy(2,indsy)+symrel_cart(3,3,jsym)*magngz(2,indsy)
             magxsu1=magxsu1+mxr*phr-mxi*phi;magxsu2=magxsu2+mxi*phr+mxr*phi
             magysu1=magysu1+myr*phr-myi*phi;magysu2=magysu2+myi*phr+myr*phi
             magzsu1=magzsu1+mzr*phr-mzi*phi;magzsu2=magzsu2+mzi*phr+mzr*phi
           end if
         end do
         rhosu1_arr(3*izone-2)=magxsu1/dble(nsym_used)
         rhosu1_arr(3*izone-1)=magysu1/dble(nsym_used)
         rhosu1_arr(3*izone  )=magzsu1/dble(nsym_used)
         rhosu2_arr(3*izone-2)=magxsu2/dble(nsym_used)
         rhosu2_arr(3*izone-1)=magysu2/dble(nsym_used)
         rhosu2_arr(3*izone  )=magzsu2/dble(nsym_used)
         numpt=numpt+nup
       end do

!      Reduction in case of FFT parallelization
       if(mpi_enreg%mode_para=='b')then
         old_paral_level=mpi_enreg%paral_level
         mpi_enreg%paral_level=3
         spaceComm=mpi_enreg%comm_fft
         call xsum_mpi(rhosu1_arr,spaceComm,ier)
         call xsum_mpi(rhosu2_arr,spaceComm,ier)
         mpi_enreg%paral_level=old_paral_level
       end if

!      Now symmetrize the magnetization at equivalent points
!      Use Transpose[symrel]
       numpt=0
       do izone=1,izone_max
         rep=irrzon(izone,2,imagn)
         nup=nsym_used/rep
         do iup=1,nup
           ind=irrzon(iup+numpt,1,imagn)
           j=ind-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2)
           if(modulo(j/n1,n2)/nd2==me_fft) then  ! this ind is to be treated by me_fft
             ind=n1*(nd2*j3+modulo(j2,nd2))+j1+1  ! this is ind in the current proc
             jsym=isymg(iup+numpt);if (jsym==0) stop "ERROR (2) in symrhg !"
             magxsu1=rhosu1_arr(3*izone-2);magxsu2=rhosu2_arr(3*izone-2)
             magysu1=rhosu1_arr(3*izone-1);magysu2=rhosu2_arr(3*izone-1)
             magzsu1=rhosu1_arr(3*izone  );magzsu2=rhosu2_arr(3*izone  )
             phr=phnons(1,iup,imagn);if (rep==1) phr=phr*symafm_used(jsym) !if rep==2, symafm is already included in phnons
             phi=phnons(2,iup,imagn);if (rep==1) phi=phi*symafm_used(jsym) !(see irrzg.F90)
             mxr=symrec_cart(1,1,jsym)*magxsu1+symrec_cart(2,1,jsym)*magysu1+symrec_cart(3,1,jsym)*magzsu1
             mxi=symrec_cart(1,1,jsym)*magxsu2+symrec_cart(2,1,jsym)*magysu2+symrec_cart(3,1,jsym)*magzsu2
             myr=symrec_cart(1,2,jsym)*magxsu1+symrec_cart(2,2,jsym)*magysu1+symrec_cart(3,2,jsym)*magzsu1
             myi=symrec_cart(1,2,jsym)*magxsu2+symrec_cart(2,2,jsym)*magysu2+symrec_cart(3,2,jsym)*magzsu2
             mzr=symrec_cart(1,3,jsym)*magxsu1+symrec_cart(2,3,jsym)*magysu1+symrec_cart(3,3,jsym)*magzsu1
             mzi=symrec_cart(1,3,jsym)*magxsu2+symrec_cart(2,3,jsym)*magysu2+symrec_cart(3,3,jsym)*magzsu2
             magngx(1,ind)=mxr*phr+mxi*phi
             magngx(2,ind)=mxi*phr-mxr*phi
             magngy(1,ind)=myr*phr+myi*phi
             magngy(2,ind)=myi*phr-myr*phi
             magngz(1,ind)=mzr*phr+mzi*phi
             magngz(2,ind)=mzi*phr-mzr*phi
           end if
         end do
         numpt=numpt+nup
       end do
       ABI_DEALLOCATE(isymg)
       ABI_DEALLOCATE(rhosu1_arr)
       ABI_DEALLOCATE(rhosu2_arr)
       ABI_DEALLOCATE(symrec_cart)
       ABI_DEALLOCATE(symrel_cart)
       ABI_DEALLOCATE(symafm_used)

     end if ! nspden==4

!    DEBUG
!    write(std_out,*)' symrhg : density after symmetrization, phnons,irrzon'
!    do ipw=1,nfft
!    if(abs(rhog(1,ipw))<1.0d-14)rhog(1,ipw)=0.0_dp
!    if(abs(rhog(2,ipw))<1.0d-14)rhog(2,ipw)=0.0_dp
!    write(std_out,*)ipw,rhog(:,ipw) !,phnons(:,ipw,1),irrzon(ipw,:,1)
!    end do
!    ENDDEBUG

     call timab(17,2,tsec)

!    Pull out full or spin up density, now symmetrized
     call fourdp(cplex,rhog,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     rhor(:,ispden)=work(:)
     if (nspden==4) then
       call fourdp(cplex,magngx,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
       rhor(:,2)=work(:)
       call fourdp(cplex,magngy,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
       rhor(:,3)=work(:)
       call fourdp(cplex,magngz,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
       rhor(:,4)=work(:)
       ABI_DEALLOCATE(magngx)
       ABI_DEALLOCATE(magngy)
       ABI_DEALLOCATE(magngz)
     end if

   end do ! ispden

!  End on the condition nsym==1
 end if

!DEBUG
!write(std_out,*)'  ispden  ifft   rhor(ifft,ispden)'
!do ifft=1,nfft,123
!do ispden=1,nspden
!if(cplex==1)write(std_out,*)ispden,ifft,rhor(ifft,ispden)
!if(cplex==2)write(std_out,*)ispden,ifft,rhor(2*ifft-1,ispden),rhor(2*ifft,ispden)
!end do
!if(cplex==1)write(std_out,*)3,ifft,rhor(ifft,1)-rhor(ifft,2)
!if(cplex==2)write(std_out,*)3,ifft,rhor(2*ifft-1,1)-rhor(2*ifft-1,2),&
!&  rhor(2*ifft,1)-rhor(2*ifft,2)
!end do
!write(std_out,*)' symrhg : leave'
!ENDDEBUG

 ABI_DEALLOCATE(work)

end subroutine symrhg
!!***
