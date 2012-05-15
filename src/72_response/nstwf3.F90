!{\src2tex{textfont=tt}}
!!****f* ABINIT/nstwf3
!! NAME
!! nstwf3
!!
!! FUNCTION
!! This routine computes the non-local contribution to the
!! 2DTE matrix elements, in the non-stationary formulation
!! Only for norm-conserving pseudopotentials (no PAW)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG,AR,MB,MVer,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions at k
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF wavefunctions at k,q.
!!  ddkfil(3)=unit numbers for the three possible ddk files for ipert1
!!       equal to 0 if no dot file is available for this direction
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig_k(mband*nsppol)=GS eigenvalues at k (hartree)
!!  eig1_k(2*nsppol*mband**2)=matrix of first-order eigenvalues (hartree)
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  icg=shift to be applied on the location of data in the array cg
!!  icg1=shift to be applied on the location of data in the array cg1
!!  idir=direction of the current perturbation
!!  ikpt=number of the k-point
!!  ipert=type of the perturbation
!!  isppol=1 for unpolarized, 2 for spin-polarized
!!  istwf_k=option parameter that describes the storage of wfs
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kg1_k(3,npw1_k)=reduced planewave coordinates at k+q, with RF k points
!!  kpt(3)=reduced coordinates of k points.
!!  mpert =maximum number of ipert
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  nband_k=number of bands at this k point for that spin polarization
!!  npw_k=number of plane waves at this k point
!!  npw1_k=number of plane waves at this k+q point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  wffddk(3)=struct info for for the three possible dot files for ipert1
!!       equal to 0 if no dot file is available for this direction
!!  wffnow=struct info for INPUT 1st-order wf file
!!  wfftgs=struct info for ground-state wf disk file
!!  wtk_k=weight assigned to the k point.
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k+q point
!!
!! OUTPUT
!!  d2bbb_k(2,3,mband,mband*prtbbb)=band by band decomposition of the second
!!   order derivatives, for the present k point, and perturbation idir, ipert
!!  d2nl_k(2,3,mpert)=non-local contributions to
!!   non-stationary 2DTE, for the present k point, and perturbation idir, ipert
!!
!! PARENTS
!!      nstdy3
!!
!! CHILDREN
!!      dotprod_g,gaugetransfo,getgh1c,mkffnl,mkkpg,ph1d3d,timab,wffreaddatarec
!!      wffreadnpwrec,wffreadskiprec,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine nstwf3(cg,cg1,ddkfil,dtfil,dtset,d2bbb_k,d2nl_k,eig_k,eig1_k,gs_hamkq,&
&                 icg,icg1,idir,ikpt,ipert,isppol,istwf_k,kg_k,kg1_k,kpt,mpert,&
&                 mpi_enreg,mpsang,mpw,mpw1,nband_k,npw_k,npw1_k,nsppol,&
&                 occ_k,psps,rmet,wffddk,wffnow,wfftgs,wtk_k,ylm,ylm1)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nstwf3'
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_spacepar
 use interfaces_59_io_mpi
 use interfaces_65_nonlocal
 use interfaces_72_response, except_this_one => nstwf3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,icg1,idir,ikpt,ipert,isppol,istwf_k
 integer,intent(in) :: mpert,mpsang,mpw,mpw1,nsppol
 integer,intent(inout) :: nband_k,npw1_k,npw_k
 real(dp),intent(in) :: wtk_k
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),intent(in) :: gs_hamkq
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow,wfftgs
!arrays
 integer,intent(in) :: ddkfil(3),kg1_k(3,npw1_k)
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*dtset%mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*dtset%nspinor*dtset%mband*dtset%mk1mem*nsppol)
 real(dp),intent(in) :: eig_k(dtset%mband*nsppol),kpt(3),occ_k(nband_k),rmet(3,3)
 real(dp),intent(in) :: ylm(npw_k,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(npw1_k,mpsang*mpsang*psps%useylm)
 real(dp),intent(inout) :: eig1_k(2*nsppol*dtset%mband**2)
 real(dp),intent(out) :: d2bbb_k(2,3,dtset%mband,dtset%mband*dtset%prtbbb),d2nl_k(2,3,mpert)
 type(wffile_type),intent(inout) :: wffddk(3)

!Local variables-------------------------------
!scalars
 integer :: berryopt,cplex,dime1kb,dimffnl,dimffnl1,dimphkxred
 integer :: iband,ider,idir1,ierr,ilmn,ipert1,ipw
 integer :: itypat,jband,matblk,me,n1,n2,n3,nband_kocc,nkpg,nkpg1
 integer :: npw_disk,nsp,optlocal,optnl,opt_gvnl1,sij_opt,tim_getgh1c,usecprj,usee1kb,usevnl
 real(dp) :: aa,arg,dot1i,dot1r,dot2i,dot2r,dot_ndiagi,dot_ndiagr,doti,dotr
 real(dp) :: lambda
 character(len=500) :: msg
!arrays
 integer :: dum_dimcprj(1),pspso_typ(1)
 integer,allocatable :: indlmn_typ(:,:,:)
 real(dp) :: dum_grad_berry(1,1),dum_gvnl1(1,1),dum_gs1(1,1)
 real(dp) :: dum_vlocal1(1,1),dum_wfraug(1,1),dum_ylmgr(1,3,1)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cg_k(:,:),cwave0(:,:),cwavef(:,:),cwavef_da(:,:)
 real(dp),allocatable :: cwavef_db(:,:),dkinpw(:),e1kbfr(:,:,:),eig2_k(:)
 real(dp),allocatable :: ekb_typ(:,:,:),ffnl(:,:,:,:)
 real(dp),allocatable :: ffnl1(:,:,:,:),ffnlk(:,:,:,:),ffnlkq(:,:,:,:)
 real(dp),allocatable :: gvnl1(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:)
 real(dp),allocatable :: ph3d(:,:,:),phkxred(:,:),sij_typ(:,:)
 type(cprj_type),allocatable :: dum_cwaveprj(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Not valid for PAW
 if (psps%usepaw==1) then
   msg='  This routine cannot be used for PAW (use pawnst3 instead) !'
   MSG_BUG(msg)
 end if

!Keep track of total time spent in nstwf3
 call timab(102,1,tsec)
 tim_getgh1c=2

!Define me
 call xme_init(mpi_enreg,me)

!Useful data for non-local operator
 dimffnl=1;dimffnl1=1
 usee1kb=0;dime1kb=0
 if (ipert/=dtset%natom+1) then
   ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
   ABI_ALLOCATE(ffnlk,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
   ABI_ALLOCATE(ffnlkq,(npw1_k,dimffnl,psps%lmnmax,1))
   ABI_ALLOCATE(ffnl1,(npw1_k,dimffnl1,psps%lmnmax,psps%ntypat))
   ABI_ALLOCATE(dum_cwaveprj,(0,0))
   ABI_ALLOCATE(ekb_typ,(gs_hamkq%dimekb1,1,dtset%nspinor**2))
   ABI_ALLOCATE(sij_typ,(gs_hamkq%dimekb1,usee1kb))
   ABI_ALLOCATE(e1kbfr,(dime1kb,gs_hamkq%dimekb2,usee1kb))
 end if

!Additional allocations
 dimphkxred=0
 if (ipert/=dtset%natom+1) then
   ABI_ALLOCATE(indlmn_typ,(6,psps%lmnmax,1))
   dimphkxred=1
   ABI_ALLOCATE(phkxred,(2,dimphkxred))
   ABI_ALLOCATE(dkinpw,(npw_k))
   ABI_ALLOCATE(kinpw1,(npw1_k))
   kinpw1=zero
 end if
 ABI_ALLOCATE(gvnl1,(2,npw1_k*dtset%nspinor))
 ABI_ALLOCATE(eig2_k,(2*nsppol*dtset%mband**2))
 ABI_ALLOCATE(cwave0,(2,npw_k*dtset%nspinor))
 ABI_ALLOCATE(cwavef,(2,npw1_k*dtset%nspinor))

!Compute (k+G) vectorsand  nonlocal form factors ffnl at (k+G)
 nkpg=0
 if (ipert/=dtset%natom+1) then
   nkpg=3*gs_hamkq%nloalg(5)
   ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
   if (nkpg>0) call mkkpg(kg_k,kpg_k,kpt,nkpg,npw_k)
   ider=0
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,gs_hamkq%gmet,&
&   gs_hamkq%gprimd,ider,ider,psps%indlmn,kg_k,kpg_k,kpt,psps%lmnmax,&
&   psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,psps%ntypat,psps%pspso,&
&   psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm,dum_ylmgr)
 end if


!Compute (k+q+G) vectors and nonlocal form factors ffnl1 at (k+q+G)
 nkpg1=0
 if (ipert/=dtset%natom+1) then
   nkpg1=3*gs_hamkq%nloalg(5)
   ABI_ALLOCATE(kpg1_k,(npw1_k,nkpg1))
   if (nkpg1>0) call mkkpg(kg1_k,kpg1_k,gs_hamkq%kpoint,nkpg1,npw1_k)
   ider=0
   call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1,psps%ffspl,gs_hamkq%gmet,&
&   gs_hamkq%gprimd,ider,ider,psps%indlmn,kg1_k,kpg1_k,gs_hamkq%kpoint,&
&   psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,npw1_k,psps%ntypat,&
&   psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1,dum_ylmgr)
 end if

!Compute ph3d (note: use npw1_k)
 if (ipert/=dtset%natom+1) then
   if(gs_hamkq%nloalg(1)<=0)then ! Here, allocation as well as precomputation
     matblk=gs_hamkq%nloalg(4)
     ABI_ALLOCATE(ph3d,(2,npw1_k,matblk))
   else                         ! Here, allocation as well as precomputation
     matblk=dtset%natom
     ABI_ALLOCATE(ph3d,(2,npw1_k,matblk))
     n1=gs_hamkq%ngfft(1);n2=gs_hamkq%ngfft(2);n3=gs_hamkq%ngfft(3)
     call ph1d3d(1,dtset%natom,kg1_k,matblk,dtset%natom,npw1_k,n1,n2,n3,gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d)
   end if
 end if

!Take care of the npw and kg records
!NOTE : one should be able to modify the rwwf routine to take care
!of the band parallelism, which is not the case yet ...
 do idir1=1,3
   if (ddkfil(idir1)/=0)then
!    Read npw record
     nsp=dtset%nspinor
     call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw_disk,nsp,wffddk(idir1))
     if (npw_k /= npw_disk) then
       write(unit=msg,fmt='(a,i3,a,i5,a,i3,a,a,i5,a,a,i5)')&
&       ' For isppol = ',isppol,', ikpt = ',ikpt,' and idir = ',idir,ch10,&
&       ' the number of plane waves in the ddk file is equal to', npw_disk,ch10,&
&       ' while it should be ',npw_k
       MSG_BUG(msg)
     end if
!    Skip k+G record
     call WffReadSkipRec(ierr,1,wffddk(idir1))
   end if
 end do
 if (dtset%mkmem==0) then
   nsp=dtset%nspinor
   call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw_k,nsp,wfftgs)
!  Skip k+G and eigenvalue records in wfftgs
   call WffReadSkipRec(ierr,2,wfftgs)
 end if
 if (dtset%mk1mem==0) then
   nsp=dtset%nspinor
   call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw1_k,nsp,wffnow)
!  Skip k+G record
   call WffReadSkipRec(ierr,1,wffnow)
 end if

 if (ipert==dtset%natom+1) then
   nband_kocc = 0
   do iband = 1,nband_k
     if (abs(occ_k(iband)) > tol8) nband_kocc = nband_kocc + 1
   end do
 end if

 if(dtset%prtbbb==1)then
   ABI_ALLOCATE(cwavef_da,(2,npw1_k*dtset%nspinor))
   ABI_ALLOCATE(cwavef_db,(2,npw1_k*dtset%nspinor))
   ABI_ALLOCATE(cg_k,(2,npw_k*dtset%nspinor*nband_k))
   if ((ipert == dtset%natom + 1).or.(ipert <= dtset%natom).or. &
&   (ipert == dtset%natom + 2).or.(ipert == dtset%natom + 5)) then
     if (dtset%mkmem /= 0) then
       cg_k(:,:) = cg(:,1+icg:icg+nband_k*npw_k*dtset%nspinor)
     else
       do iband=1,nband_k
         call WffReadDataRec(cwave0,ierr,2,npw_k*dtset%nspinor,wfftgs)
         cg_k(:,(iband-1)*npw_k*dtset%nspinor+1:iband*npw_k*dtset%nspinor)=cwave0(:,:)
       end do
     end if
   end if
   d2bbb_k(:,:,:,:) = zero
 end if

!Loop over bands
 do iband=1,nband_k

   if(mpi_enreg%paral_compil_kpt==1)then
     if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= me) then
       do idir1=1,3
!        Skip the eigenvalue and the wf records of this band
         if (ddkfil(idir1) /= 0) then
           call WffReadSkipRec(ierr,2,wffddk(idir1))
         end if
       end do
       if(dtset%mkmem==0)then
         if(dtset%prtbbb==0 .or. ipert==dtset%natom+2)then
           call WffReadSkipRec(ierr,1,wfftgs)
         end if
       end if
       if(dtset%mk1mem==0)then
         call WffReadSkipRec(ierr,2,wffnow)
       end if
       cycle
     end if
   end if ! mpi_enreg%paral_compil_kpt==1

!  Read ground-state wavefunctions
   if (dtset%prtbbb==0 .or. ipert==dtset%natom+2 .or. ipert==dtset%natom+5) then
     if(dtset%mkmem/=0) then
       cwave0(:,:)=cg(:,1+(iband-1)*npw_k*dtset%nspinor+icg:iband*npw_k*dtset%nspinor+icg)
     else
       call timab(286,1,tsec)
       call WffReadDataRec(cwave0,ierr,2,npw_k*dtset%nspinor,wfftgs)
       call timab(286,2,tsec)
     end if
   else    ! prtbbb==1 and ipert<=natom , already in cg_k
     cwave0(:,:)=cg_k(:,1+(iband-1)*npw_k*dtset%nspinor:iband*npw_k*dtset%nspinor)
   end if

!  Read first-order wavefunctions
   if(dtset%mk1mem/=0)then
     cwavef(:,:)=cg1(:,1+(iband-1)*npw1_k*dtset%nspinor+icg1:iband*npw1_k*dtset%nspinor+icg1)
   else
     call timab(286,1,tsec)
     call WffReadDataRec(eig1_k(1+(iband-1)*2*nband_k:2*iband*nband_k),ierr,2*nband_k,wffnow)
     call WffReadDataRec(cwavef,ierr,2,npw1_k*dtset%nspinor,wffnow)
     call timab(286,2,tsec)
   end if

!  In case non ddk perturbation
   if (ipert /= dtset%natom + 1) then

     do ipert1=1,mpert

       if( ipert1<=dtset%natom .or. ipert1==dtset%natom+2 )then

         if( ipert1<=dtset%natom )then
!          Transfer the data relative to the displaced atom: ekb, indlmn, pspso, ffspl
           itypat=gs_hamkq%typat(ipert1)
           ekb_typ(:,1,1)=psps%ekb(:,itypat)
           if (dtset%nspinor==2) then
             ekb_typ(:,1,2)=psps%ekb(:,itypat)
             ekb_typ(:,1,3:4)=zero
           end if
           indlmn_typ(:,:,1)=psps%indlmn(:,:,itypat)
           pspso_typ(1)=psps%pspso(itypat)
!          Copy the part needed for the displaced atom, in ffnlkq.
           do ilmn=1,psps%lmnmax
             ffnlkq(:,1:dimffnl,ilmn,1)=ffnl1(:,1:dimffnl,ilmn,itypat)
             ffnlk (:,1:dimffnl,ilmn,1)=ffnl (:,1:dimffnl,ilmn,itypat)
           end do
         end if ! ipert1

         if (((ipert <= dtset%natom).or.(ipert == dtset%natom + 2)) &
&         .and.(ipert1 == dtset%natom+2).and. dtset%prtbbb==1) then
           call gaugetransfo(cg_k,cwavef,cwavef_db,eig_k,eig1_k,iband,nband_k, &
&           dtset%mband,npw_k,npw1_k,dtset%nspinor,nsppol,occ_k)
           cwavef(:,:) = cwavef_db(:,:)
         end if

!        Define the direction along which to move the atom :
!        the polarisation (ipert1,idir1) is refered as j1.
         do idir1=1,3
           if (ipert1<=dtset%natom.or.(ipert1==dtset%natom+2.and.ddkfil(idir1)/=0)) then

!            Get |Vnon-locj^(1)|u0> :
!            First-order non-local, applied to zero-order wavefunction
!            This routine gives MINUS the non-local contribution

!            ==== Atomic displ. perturbation
             if( ipert1<=dtset%natom )then
!              Compute here phkxred for kpt
               arg=two_pi*(kpt(1)*gs_hamkq%xred(1,ipert1)+kpt(2)*gs_hamkq%xred(2,ipert1) &
&               +kpt(3)*gs_hamkq%xred(3,ipert1))
               phkxred(1,1)=cos(arg) ; phkxred(2,1)=sin(arg)
               lambda=eig_k((isppol-1)*nband_k+iband)
               berryopt=1;cplex=1;optlocal=0;optnl=1;usecprj=0;usevnl=0;opt_gvnl1=0;sij_opt=0
               call getgh1c(berryopt,cplex,cwave0,dum_cwaveprj,dum_dimcprj,gs_hamkq%dimekb1,dime1kb,&
&               dimffnl,dimffnl,dimffnl1,dimphkxred,dkinpw,ekb_typ,e1kbfr,e1kbfr,ffnlk,ffnlkq,&
&               ffnl1,dtfil%filstat,gs_hamkq%gbound,gvnl1,dum_grad_berry,dum_gs1,gs_hamkq,dum_gvnl1,&
&               idir1,indlmn_typ,ipert1,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,lambda,psps%lmnmax,&
&               matblk,dtset%mgfft,mpi_enreg,psps%mpsang,psps%mpssoang,dtset%natom,nkpg,nkpg1,&
&               npw_k,npw1_k,dtset%nspinor,gs_hamkq%ntypat,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,&
&               optlocal,optnl,opt_gvnl1,dtset%paral_kgb,ph3d,phkxred,dtset%prtvol,sij_opt,&
&               sij_typ,tim_getgh1c,usecprj,usee1kb,usevnl,dum_vlocal1,dum_wfraug)

!              ==== Electric field perturbation
             else if( ipert1==dtset%natom+2 )then
               call WffReadDataRec(eig2_k(1+(iband-1)*2*nband_k:2*iband*nband_k),ierr,2*nband_k,wffddk(idir1))
               call WffReadDataRec(gvnl1,ierr,2,npw1_k*dtset%nspinor,wffddk(idir1))
!              In case of band-by-band,
!              construct the first-order wavefunctions in the diagonal gauge
               if (((ipert <= dtset%natom).or.(ipert == dtset%natom + 2)).and.(dtset%prtbbb==1)) then
                 call gaugetransfo(cg_k,gvnl1,cwavef_da,eig_k,eig2_k,iband,nband_k, &
&                 dtset%mband,npw_k,npw1_k,dtset%nspinor,nsppol,occ_k)
                 gvnl1(:,:) = cwavef_da(:,:)
               end if
!              Multiplication by -i
               do ipw=1,npw1_k*dtset%nspinor
                 aa=gvnl1(1,ipw)
                 gvnl1(1,ipw)=gvnl1(2,ipw)
                 gvnl1(2,ipw)=-aa
               end do
             end if

!            MVeithen 021212 :
!            1) Case ipert1 = natom + 2 and ipert = natom + 2:
!            the second derivative of the energy with respect to an electric
!            field is computed from Eq. (38) of X. Gonze, PRB 55 ,10355 (1997).
!            The evaluation of this formula needs the operator $i \frac{d}{dk}.
!            2) Case ipert1 = natom + 2 and ipert < natom:
!            the computation of the Born effective charge tensor uses
!            the operator $-i \frac{d}{dk}.
             if (ipert==dtset%natom+2) gvnl1(:,:) = -gvnl1(:,:)

!            <G|Vnl1|Cnk> is contained in gvnl1
!            construct the matrix element (<uj2|vj1|u0>)complex conjug.
!            and add it to the 2nd-order matrix
!            XG030513 : use dotprod_g, for future parallelisation
             call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1_k*dtset%nspinor,2,cwavef,gvnl1)
             d2nl_k(1,idir1,ipert1)=d2nl_k(1,idir1,ipert1)+wtk_k*occ_k(iband)*two*dotr
             d2nl_k(2,idir1,ipert1)=d2nl_k(2,idir1,ipert1)-wtk_k*occ_k(iband)*two*doti

!            Band by band decomposition of the Born effective charges
!            calculated from a phonon perturbation
             if(dtset%prtbbb==1)then
               d2bbb_k(1,idir1,iband,iband) = wtk_k*occ_k(iband)*two*dotr
               d2bbb_k(2,idir1,iband,iband) = -one*wtk_k*occ_k(iband)*two*doti
             end if

           end if
         end do
       end if
     end do
   end if     ! ipert /= natom +1

!  Compute the localization tensor

   if (ipert==dtset%natom+1) then

     ipert1=dtset%natom+1
     if(dtset%prtbbb==1)then
       call gaugetransfo(cg_k,cwavef,cwavef_db,eig_k,eig1_k,iband,nband_k, &
&       dtset%mband,npw_k,npw1_k,dtset%nspinor,nsppol,occ_k)
       cwavef(:,:) = cwavef_db(:,:)
     end if

     do idir1 = 1,3
       eig2_k(:) = zero
       gvnl1(:,:) = zero
       if (idir == idir1) then
         if (ddkfil(idir1) /= 0) then
           call WffReadSkipRec(ierr,2,wffddk(idir1))
         end if
         gvnl1(:,:) = cwavef(:,:)
         eig2_k(:) = eig1_k(:)
       else
         if (ddkfil(idir1) /= 0) then
           call WffReadDataRec(eig2_k(1+(iband-1)*2*nband_k:2*iband*nband_k),ierr,2*nband_k,wffddk(idir1))
           call WffReadDataRec(gvnl1,ierr,2,npw1_k*dtset%nspinor,wffddk(idir1))

           if(dtset%prtbbb==1)then
             call gaugetransfo(cg_k,gvnl1,cwavef_da,eig_k,eig2_k,iband,nband_k, &
&             dtset%mband,npw_k,npw1_k,dtset%nspinor,nsppol,occ_k)
             gvnl1(:,:) = cwavef_da(:,:)
           end if

         end if    !ddkfil(idir1)
       end if    !idir == idir1

!      <G|du/dqa> is contained in gvnl1 and <G|du/dqb> in cwavef
!      construct the matrix elements <du/dqa|du/dqb> -> dot
!      <u|du/dqa> -> dot1
!      <du/dqb|u> -> dot2
!      and add them to the 2nd-order matrix

!      XG030513 : use dotprod_g, for future parallelisation
       call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1_k*dtset%nspinor,2,gvnl1,cwavef)
       d2nl_k(1,idir1,ipert1)=d2nl_k(1,idir1,ipert1)+wtk_k*occ_k(iband)*dotr/(nband_kocc*two)
       d2nl_k(2,idir1,ipert1)=d2nl_k(2,idir1,ipert1)+wtk_k*occ_k(iband)*doti/(nband_kocc*two)


!      XG 020216 : Marek, could you check the next forty lines
!      In the parallel gauge, dot1 and dot2 vanishes
       if(dtset%prtbbb==1)then
         d2bbb_k(1,idir1,iband,iband)=d2bbb_k(1,idir1,iband,iband)+dotr
         d2bbb_k(2,idir1,iband,iband)=d2bbb_k(2,idir1,iband,iband)+doti
         dot_ndiagr=zero ; dot_ndiagi=zero
         do jband = 1,nband_k              !compute dot1 and dot2
           if (abs(occ_k(jband)) > tol8) then
             dot1r=zero ; dot1i=zero
             dot2r=zero ; dot2i=zero
             cwave0(:,:)=cg_k(:,1+(jband-1)*npw_k*dtset%nspinor:jband*npw_k*dtset%nspinor)
!            XG030513 : use dotprod_g, for future parallelisation
             call dotprod_g(dot1r,dot1i,istwf_k,mpi_enreg,npw1_k*dtset%nspinor,2,cwave0,gvnl1)
             call dotprod_g(dot2r,dot2i,istwf_k,mpi_enreg,npw1_k*dtset%nspinor,2,cwavef,cwave0)
             dot_ndiagr = dot_ndiagr + dot1r*dot2r - dot1i*dot2i
             dot_ndiagi = dot_ndiagi + dot1r*dot2i + dot1i*dot2r
             d2bbb_k(1,idir1,iband,jband) = d2bbb_k(1,idir1,iband,jband) - &
&             (dot1r*dot2r - dot1i*dot2i)
             d2bbb_k(2,idir1,iband,jband) = d2bbb_k(2,idir1,iband,jband) - &
&             (dot1r*dot2i + dot1i*dot2r)
           end if  ! occ_k
         end do !jband
         d2bbb_k(:,idir1,iband,:)=d2bbb_k(:,idir1,iband,:)*wtk_k*occ_k(iband)*half
         d2nl_k(1,idir1,ipert1)= &
&         d2nl_k(1,idir1,ipert1)-wtk_k*occ_k(iband)*dot_ndiagr/(nband_kocc*two)
         d2nl_k(2,idir1,ipert1)=&
&         d2nl_k(2,idir1,ipert1)-wtk_k*occ_k(iband)*dot_ndiagi/(nband_kocc*two)
       end if ! prtbbb==1

     end do  ! idir1
   end if   ! Compute localization tensor, ipert=natom+1

!  End loop over bands
 end do

!Final deallocations
 ABI_DEALLOCATE(cwave0)
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(eig2_k)
 ABI_DEALLOCATE(gvnl1)
 if (ipert/=dtset%natom+1) then
   ABI_DEALLOCATE(ffnl)
   ABI_DEALLOCATE(ffnl1)
   ABI_DEALLOCATE(ffnlkq)
   ABI_DEALLOCATE(ffnlk)
   ABI_DEALLOCATE(kpg_k)
   ABI_DEALLOCATE(kpg1_k)
   ABI_DEALLOCATE(phkxred)
   ABI_DEALLOCATE(ph3d)
   ABI_DEALLOCATE(ekb_typ)
   ABI_DEALLOCATE(sij_typ)
   ABI_DEALLOCATE(e1kbfr)
   ABI_DEALLOCATE(indlmn_typ)
   ABI_DEALLOCATE(dkinpw)
   ABI_DEALLOCATE(kinpw1)
   ABI_DEALLOCATE(dum_cwaveprj)
 end if
 if(dtset%prtbbb==1)  then
   ABI_DEALLOCATE(cg_k)
   ABI_DEALLOCATE(cwavef_da)
   ABI_DEALLOCATE(cwavef_db)
 end if

 call timab(102,2,tsec)

 DBG_EXIT("COLL")

end subroutine nstwf3
!!***
