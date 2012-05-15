!{\src2tex{textfont=tt}}
!!****f* ABINIT/wfd_mkrho
!! NAME
!! wfd_mkrho
!!
!! FUNCTION
!! Calculate the charge density on the fine FFT grid in real space.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MG,GMR, VO, LR, RWG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ngfftf(18)=array containing all the information for the "fine" FFT.
!!  Cryst<Crystal_structure> Info on the crystalline structure
!!     %nsym=number of symmetry operations.
!!     %ucvol=unit cell volume.
!!  optcalc=option for calculation. If =0 (default value) perform calculation
!!    of electronic density. If =1, perform calculation of kinetic energy density.
!!    In both cases, the result is returned in rhor.
!!  Psps<type(pseudopotential_type)>=variables related to pseudopotentials
!!  nfftf=Total number of points on the fine FFT grid (for this processor)
!!  Kmesh<bz_mesh_type>= Info on the k-sampling:
!!     %nibz=number of irreducible k-points.
!!     %nbz=number of k-points in the full Brillouin zone.
!!     %wt(nibz)=irreducible k-points weights.
!!     %timrev=2 if time-reversal symmetry can be used, 1 otherwise.
!!  Wfd<wfs_descriptor)=datatype gathering info on the wavefunctions.
!!    %npwwfn=Number of plane waves used to describe the wave functions.
!!    %nspinor=number of spinorial components.
!!    %nsppol=1 for unpolarized, 2 for spin-polarized calculations.
!!    %nspden=number of spin-density components.
!! [optcalc]=Optional option used to calculated the kinetic energy density. Defaults to 0.
!!
!! OUTPUT
!!  rhor(nfftf,nspden)=The density in the real space on the fine FFT grid.
!!   If nsppol==2, total charge in first half, spin-up component in second half.
!!   If optcalc==1 (optional argument, default value is 0), then rhor will actually
!!   contain kinetic energy density (taur) instead of electronic density.
!!
!! NOTES
!! In the case of PAW calculations:
!!    All computations are done on the fine FFT grid.
!!    All variables (nfftf,ngfftf,mgfftf) refer to this fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!    Developers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!! In the case of norm-conserving calculations:
!!    The mesh is the usual augmented FFT grid to treat correctly the convolution.
!!
!! PARENTS
!!      bethe_salpeter,screening,sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wfd_mkrho(Wfd,Cryst,Psps,Kmesh,Bstr,ngfftf,nfftf,rhor,&
&                    optcalc) ! optional arguments

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_iterators

 use m_io_tools,  only : get_unit
 use m_crystal,   only : crystal_structure
 use m_bz_mesh,   only : bz_mesh_type
 use m_wfs,       only : wfs_descriptor, wfd_get_ur, wfd_update_bkstab, wfd_iterator_bks, wfd_change_ngfft, fft_onewfn

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_mkrho'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_56_recipspace
 use interfaces_67_common
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf
 integer,intent(in),optional :: optcalc
 type(Bandstructure_type),intent(in) :: Bstr
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Crystal_structure),intent(in) :: Cryst
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(out) :: rhor(nfftf,Wfd%nspden)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp=5
 integer :: cplex,ib,ib_iter,ierr,ik,ir,is,n1,n2,n3,nfftotf,unt
 integer :: alpha,nalpha,ipw!,ipwsp
 integer :: myoptcalc
 real(dp) :: kpt_cart,kg_k_cart,gp2pi1,gp2pi2,gp2pi3,cwftmp
 character(len=100) :: frmt
 character(len=500) :: msg
 character(len=fnlen) :: filnam
!arrays
 integer,allocatable :: irrzon(:,:,:)
 real(dp),allocatable :: phnons(:,:,:),rhog(:,:),rhor_down(:),rhor_mx(:),rhor_my(:)
 real(dp),allocatable :: cwavef(:,:)
 real(dp),pointer :: occfact(:,:,:)
 complex(dpc),allocatable :: wfr_x(:),wfr_y(:)
 complex(gwpc),allocatable :: gradug(:),work(:)
 complex(gwpc),allocatable,target :: wfr(:)
 complex(gwpc),pointer :: cwavef1(:),cwavef2(:)
 type(iter2_t) :: Iter_bks

!*************************************************************************

 DBG_ENTER("COLL")

 ! === Consistency check ===
 if (Wfd%nspden==1.and.Wfd%nspinor==2) then
   MSG_ERROR('nspden==1 and nspinor==2 not implemented')
 end if
                                                                                     
 ABI_CHECK(Wfd%nsppol==Bstr%nsppol,"mismatch in nspppol")

 occfact => Bstr%occ(1:Bstr%mband,1:Bstr%nkpt,1:Bstr%nsppol)

 if ( ANY(ngfftf(1:3) /= Wfd%ngfft(1:3)) ) call wfd_change_ngfft(Wfd,Cryst,Psps,ngfftf) 
 !
 ! === Calculate IBZ contribution to the charge density ===
 ABI_ALLOCATE(wfr,(nfftf*Wfd%nspinor))

 if (Wfd%nspinor==2) then
   ABI_ALLOCATE(wfr_x,(nfftf))
   ABI_ALLOCATE(wfr_y,(nfftf))
   if (Wfd%nspden==4) then
     ABI_ALLOCATE(rhor_down,(nfftf))
     ABI_ALLOCATE(rhor_mx,(nfftf))
     ABI_ALLOCATE(rhor_my,(nfftf))
     rhor_down=zero
     rhor_mx  =zero
     rhor_my  =zero
   else !TODO
     MSG_ERROR('nspden/=4 and nspinor=2 not implemented') 
   end if
 end if

 ! Update the (b,k,s) distribution table.
 call wfd_update_bkstab(Wfd) 

 ! ==== Calculate the unsymmetrized density ====
 rhor=zero
 Iter_bks = wfd_iterator_bks(Wfd,bks_mask=ABS(Bstr%occ)>=tol8)

 myoptcalc=0
 if(present(optcalc))myoptcalc=optcalc
 nalpha=1
 if(myoptcalc==1)nalpha=3   

 do alpha=1,nalpha  
   do is=1,Wfd%nsppol
     do ik=1,Wfd%nkibz
       do ib_iter=1,iter_len(Iter_bks,ik,is)
         ib = yield(Iter_bks,ib_iter,ik,is)

         call wfd_get_ur(Wfd,ib,ik,is,wfr)

         cwavef1 => wfr(1:nfftf)
         if(myoptcalc==1) then
           ABI_ALLOCATE(gradug,(Wfd%Kdata(ik)%npw))
           ABI_ALLOCATE(cwavef,(2,Wfd%Kdata(ik)%npw))
           ABI_ALLOCATE(work,(nfftf))
           cwavef(1,:)= REAL(Wfd%Wave(ib,ik,is)%ug(:))
           cwavef(2,:)=AIMAG(Wfd%Wave(ib,ik,is)%ug(:))
!          Multiplication by 2pi i (k+G)_alpha
           gp2pi1=Cryst%gprimd(alpha,1)*two_pi ; gp2pi2=Cryst%gprimd(alpha,2)*two_pi ; gp2pi3=Cryst%gprimd(alpha,3)*two_pi
           kpt_cart=gp2pi1*Wfd%kibz(1,ik)+gp2pi2*Wfd%kibz(2,ik)+gp2pi3*Wfd%kibz(3,ik)
           do ipw=1,Wfd%Kdata(ik)%npw
             kg_k_cart=gp2pi1*Wfd%Kdata(ik)%kg_k(1,ipw)+gp2pi2*Wfd%Kdata(ik)%kg_k(2,ipw)+gp2pi3*Wfd%Kdata(ik)%kg_k(3,ipw)+kpt_cart
!             ipwsp=ipw!+(ispinor-1)*Wfd%Kdata(ik)%npw
             cwftmp=-cwavef(2,ipw)*kg_k_cart
             cwavef(2,ipw)=cwavef(1,ipw)*kg_k_cart
             cwavef(1,ipw)=cwftmp
           end do
           gradug(:)=CMPLX(cwavef(1,:),cwavef(2,:),gwpc)
           call fft_onewfn(Wfd%paral_kgb,Wfd%istwfk(ik),Wfd%nspinor,Wfd%npwarr(ik),nfftf,Wfd%mgfft,Wfd%ngfft,&
           &               gradug,work,Wfd%Kdata(ik)%igfft0,Wfd%Kdata(ik)%kg_k,Wfd%Kdata(ik)%gbound,tim_fourdp,Wfd%MPI_enreg)
           cwavef1(:)=work(:)
           ABI_DEALLOCATE(work)
           ABI_DEALLOCATE(cwavef)
           ABI_DEALLOCATE(gradug)
         end if        

         do ir=1,nfftf
           rhor(ir,is)=rhor(ir,is)+occfact(ib,ik,is)*CONJG(cwavef1(ir))*cwavef1(ir)*Kmesh%wt(ik)/Cryst%ucvol
         end do

         if (Wfd%nspinor==2) then
           cwavef2 => wfr(1+nfftf:2*nfftf)
           wfr_x(:)=cwavef1(:)+cwavef2(:)       ! $(\Psi^{1}+\Psi^{2})$
           wfr_y(:)=cwavef1(:)-j_dpc*cwavef2(:) ! $(\Psi^{1}-i\Psi^{2})$
           do ir=1,nfftf
             rhor_down(ir)=rhor_down(ir)+occfact(ib,ik,is)*CONJG(cwavef2(ir))*cwavef2(ir)*Kmesh%wt(ik)/Cryst%ucvol
             rhor_mx  (ir)=rhor_mx  (ir)+occfact(ib,ik,is)*CONJG(wfr_x  (ir))*wfr_x  (ir)*Kmesh%wt(ik)/Cryst%ucvol
             rhor_my  (ir)=rhor_my  (ir)+occfact(ib,ik,is)*CONJG(wfr_y  (ir))*wfr_y  (ir)*Kmesh%wt(ik)/Cryst%ucvol
           end do
         end if

       end do
     end do
   end do

 end do ! enddo alpha  

 if(myoptcalc==1)rhor(:,:)=half*rhor(:,:) ! convention for taur = 1/2 Sum_i |grad phi_i|^2

 call iter_destroy(Iter_bks)
 call xsum_mpi(rhor,Wfd%comm,ierr)
 !
 ! === Symmetrization in G-space implementing also the AFM case ===
 n1=ngfftf(1)
 n2=ngfftf(2)
 n3=ngfftf(3)
 nfftotf=n1*n2*n3

 ABI_ALLOCATE(irrzon,(nfftotf**(1-1/Cryst%nsym),2,(Wfd%nspden/Wfd%nsppol)-3*(Wfd%nspden/4)))
 ABI_ALLOCATE(phnons,(2,nfftotf,(Wfd%nspden/Wfd%nsppol)-3*(Wfd%nspden/4)))

 if (Cryst%nsym/=1) then
   call irrzg(irrzon,Wfd%nspden,Wfd%nsppol,Cryst%nsym,n1,n2,n3,phnons,Cryst%symafm,Cryst%symrel,Cryst%tnons)
 end if

 cplex=1
 ABI_ALLOCATE(rhog,(2,cplex*nfftf))

 call symrhg(cplex,Cryst%gprimd,irrzon,Wfd%MPI_enreg,nfftf,nfftotf,ngfftf,Wfd%nspden,Wfd%nsppol,&
&  Cryst%nsym,Wfd%paral_kgb,phnons,rhog,rhor,Cryst%rprimd,Cryst%symafm,Cryst%symrel)

 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(phnons)
 ABI_DEALLOCATE(irrzon)

 write(msg,'(a,f9.4)')' planewave contribution to nelect: ',SUM(rhor(:,1))*Cryst%ucvol/nfftf
 call wrtout(std_out,msg,'COLL')

 if (Wfd%nspden==4) then
   write(msg,'(a,3f9.4)')&
&     ' mx, my, mz: ',SUM(rhor(:,2))*Cryst%ucvol/nfftf,SUM(rhor(:,3))*Cryst%ucvol/nfftf,SUM(rhor(:,4))*Cryst%ucvol/nfftf
   call wrtout(std_out,msg,'COLL')
 end if

 ABI_DEALLOCATE(wfr)

 if (Wfd%nspinor==2) then
   ABI_DEALLOCATE(wfr_x)
   ABI_DEALLOCATE(wfr_y)
   if (Wfd%nspden==4)  then
     ABI_DEALLOCATE(rhor_down)
     ABI_DEALLOCATE(rhor_mx)
     ABI_DEALLOCATE(rhor_my)
   end if
 end if

 if (.FALSE..and.Wfd%my_rank==Wfd%master) then
   filnam='__rhor__.dat'
   call isfile(filnam,'new')
   unt=get_unit(); open(unit=unt,file=filnam)
   write(frmt,*)'(2x,',Wfd%nspden,'(1x,f8.3))'
   do ir=1,nfftf
     write(unt,frmt)(rhor(ir,:))
   end do
   close(unt)
 end if

 DBG_EXIT("COLL")

end subroutine wfd_mkrho
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/test_charge
!! NAME
!! test_charge
!!
!! FUNCTION
!!  Reports info on the electronic charge as well as Drude plasma frequency.
!!  Mainly used in the GW part.
!!
!! INPUTS
!!  nelectron_exp=Expected total number of electrons (used to normalize the charge)
!!
!! OUTPUT
!!
!! PARENTS
!!      bethe_salpeter,mrgscr,screening,sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine test_charge(nfftf,nelectron_exp,nspden,rhor,ucvol,&
& usepaw,usexcnhat,usefinegrid,compch_sph,compch_fft,omegaplasma)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'test_charge'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,nspden,usefinegrid,usepaw,usexcnhat
 real(dp),intent(in) :: compch_fft,compch_sph,ucvol,nelectron_exp
 real(dp),intent(out) :: omegaplasma
!arrays
 real(dp),intent(inout) :: rhor(nfftf,nspden)

!Local variables ------------------------------
!scalars
 real(dp) :: nelectron_tot,nelectron_fft 
 real(dp) :: nelectron_pw,nelectron_sph,rhoav,rs,nratio
 character(len=500) :: msg

!*************************************************************************

! ABI_UNUSED(usexcnhat)
if (usexcnhat==0)then
end if

 ! === For PAW output of compensation charges ===
 if (usepaw==1) then
!if (usepaw==1.and.usexcnhat>0) then ! TODO I still dont understand this if!
   write(msg,'(4a)')ch10,' PAW TEST:',ch10,' ==== Compensation charge inside spheres ============'
   if (compch_sph<greatest_real.and.compch_fft<greatest_real) &
&    write(msg,'(3a)')TRIM(msg),ch10,' The following values must be close...'
   if (compch_sph<greatest_real) &
&    write(msg,'(3a,f22.15)')TRIM(msg),ch10,' Compensation charge over spherical meshes = ',compch_sph
   if (compch_fft<greatest_real) then
     if (usefinegrid==1) then
       write(msg,'(3a,f22.15)')TRIM(msg),ch10,' Compensation charge over fine fft grid    = ',compch_fft
     else
       write(msg,'(3a,f22.15)')TRIM(msg),ch10,' Compensation charge over fft grid         = ',compch_fft
     end if
   end if
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a)')ch10
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if !PAW

 nelectron_pw =SUM(rhor(:,1))*ucvol/nfftf 
 nelectron_tot=nelectron_pw
 nratio       =nelectron_exp/nelectron_tot

 if (usepaw==1) then
   nelectron_sph=nelectron_pw+compch_sph
   nelectron_fft=nelectron_pw+compch_fft
   nelectron_tot=nelectron_sph
   nratio=(nelectron_exp-nelectron_sph)/nelectron_pw
 end if

 rhoav=nelectron_tot/ucvol ; rs=(three/(four_pi*rhoav))**third
 if (usepaw==0) then
  write(msg,'(2(a,f9.4))')&
&   ' Number of electrons calculated from density = ',nelectron_tot,'; Expected = ',nelectron_exp
 else
   write(msg,'(2(a,f9.4),a)')&
&   ' Total number of electrons per unit cell = ',nelectron_sph,' (Spherical mesh), ',nelectron_fft,' (FFT mesh)'
 end if
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

!$write(msg,'(a,f9.4)')' Renormalizing smooth charge density using nratio = ',nratio
!! rhor(:,:)=nratio*rhor(:,:)

 write(msg,'(a,f9.6)')' average of density, n = ',rhoav
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,f9.4)')' r_s = ',rs
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 omegaplasma=SQRT(four_pi*rhoav)
 write(msg,'(a,f9.4,2a)')' omega_plasma = ',omegaplasma*Ha_eV,' [eV]',ch10
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')

end subroutine test_charge
!!***
