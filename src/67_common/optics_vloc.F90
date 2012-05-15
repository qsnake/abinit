!{\src2tex{textfont=tt}}
!!****f* ABINIT/optics_vloc
!! NAME
!! optics_vloc
!!
!! FUNCTION
!! Compute matrix elements need for optical conductivity in a LOCAL potential
!! and store them in a file.
!! Matrix elements = <Phi_i|Nabla|Phi_j>
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (SM,VR,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cg(2,mcg)=planewave coefficients of wavefunctions.
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  fildata= name of the output file
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum dimensioned size of npw.
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!
!! OUTPUT
!!  (only writing in a file)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      clsopn,hdr_io,hdr_skip,rdnpw,rwwf,timab,wffclose,wffopen,xcomm_init
!!      xcomm_world,xdefineoff,xexch_mpi,xme_init,xsum_master
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine optics_vloc(cg,dtfil,dtset,eigen0,gprimd,hdr,kg,mband,mcg,mkmem,mpi_enreg,mpw,&
&                       nkpt,npwarr,nsppol,wffnow)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'optics_vloc'
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mkmem,mpw,nkpt,nsppol
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg
 type(hdr_type),intent(inout) :: hdr
 type(wffile_type),intent(inout) :: wffnow
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
 real(dp),intent(inout) :: cg(2,mcg)

!Local variables-------------------------------
!scalars
 integer :: accesswff,bdtot_index,cplex,etiq,formeig,fformopt,ib,icg,ierr,ikg,ikpt
 integer :: ipw,isppol,istwf_k,iwavef,jb,jwavef
 integer :: mcg_disk,me,me_kpt,my_nspinor,nband_k,npw_k,nspinor0,old_paral_level,sender
 integer :: spaceComm_band,spaceComm_bandfftspin,spaceComm_fft,spaceComm_k,spaceComm_w
 integer :: tim_rwwf
 logical :: mykpt
 real(dp) :: cgnm1,cgnm2
 character(len=500) :: msg
!arrays
 integer :: tmp_shape(3)
 integer,allocatable :: kg_dum(:,:),kg_k(:,:)
 real(dp) :: kpoint(3),tsec(2)
 real(dp),allocatable :: cg_disk(:,:),eig_dum(:),eig0_k(:),kpg_k(:,:)
 real(dp),allocatable :: occ_dum(:),psinablapsi(:,:,:,:),tnm(:,:,:,:)
 type(wffile_type) :: wff1
! ************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 if(dtset%paral_kgb==1.and.mkmem==0)then
   msg='  Not compatible with mkmem=0 and band-fft parallelism !'
   MSG_ERROR(msg)
 end if

!----------------------------------------------------------------------------------
!2- Computation of <psi_n|-i.nabla|psi_m> for each k
!----------------------------------------------------------------------------------

!Init parallelism
 call xcomm_world(mpi_enreg,spaceComm_w,myrank=me)
 call xcomm_init(mpi_enreg,spaceComm_k,spaceComm_bandfft=mpi_enreg%comm_kpt)
 call xme_init(mpi_enreg,me_kpt)
 if (mpi_enreg%mode_para=='b') then
   spaceComm_fft=mpi_enreg%comm_fft
   spaceComm_band=mpi_enreg%comm_band
   spaceComm_bandfftspin=mpi_enreg%commcart_3d
 else
   spaceComm_band=0;spaceComm_fft=0;spaceComm_bandfftspin=0
 end if
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)

!Prepare temporary files if mkmem==0
!WF file
 if (mkmem==0) then
   formeig=0;mcg_disk=mpw*my_nspinor*mband
   call clsopn(wffnow)
   call hdr_skip(wffnow,ierr)
   call xdefineOff(formeig,wffnow,mpi_enreg,dtset%nband,npwarr,dtset%nspinor,nsppol,nkpt)
   ABI_ALLOCATE(cg_disk,(2,mcg_disk))
 end if

!Initialize main variables
 ABI_ALLOCATE(psinablapsi,(2,3,mband,mband))
 psinablapsi=zero

 accesswff=-1
 fformopt=612
 call WffOpen(accesswff,spaceComm_k,dtfil%fnameabo_app_opt,ierr,wff1,0,me,123)
 call hdr_io(fformopt,hdr,2,wff1)

!LOOP OVER SPINS
 icg=0
 do isppol=1,nsppol

!  LOOP OVER k POINTS
   ikg=0
   if (dtset%mkmem==0) rewind dtfil%unkg
   do ikpt=1,nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
     etiq=ikpt+(isppol-1)*nkpt
     if (me==0) then
       ABI_ALLOCATE(eig0_k,(nband_k))
       eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
     end if

!    Old FFT parallelism: define FFT communicator for this k-point
     if (mpi_enreg%paral_compil_fft==1.and.mpi_enreg%mode_para/='b') then
       mpi_enreg%num_group_fft=ikpt+(isppol-1)*nkpt
       old_paral_level=mpi_enreg%paral_level;mpi_enreg%paral_level=3
       call xcomm_init(mpi_enreg,spaceComm_fft)
       mpi_enreg%paral_level=old_paral_level
     end if

!    Select my k-points
     mykpt=.true.
     if(mpi_enreg%paral_compil_kpt==1)then
       mykpt=(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me_kpt))==0)
     end if
     if (mykpt) then

!      Allocations depending on k-point
       kpoint(:)=dtset%kptns(:,ikpt)
       istwf_k=dtset%istwfk(ikpt)
       npw_k=npwarr(ikpt)
       cplex=2;if (istwf_k>1) cplex=1
       ABI_ALLOCATE(kg_k,(3,npw_k))
       ABI_ALLOCATE(kpg_k,(npw_k*dtset%nspinor,3))

!      Extract G-vectors for this k-point according to mkmem
       if (mkmem==0) then
         nspinor0=dtset%nspinor
         call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor0,0,dtfil%unkg)
         read (dtfil%unkg) kg_k(1:3,1:npw_k)
       else
         kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
         ikg=ikg+npw_k
       end if

!      Calculation of k+G in cartesian coordinates
       do ipw=1,npw_k
         kpg_k(ipw,1)=(kpoint(1)+kg_k(1,ipw))*gprimd(1,1)&
&         +(kpoint(2)+kg_k(2,ipw))*gprimd(1,2)&
&         +(kpoint(3)+kg_k(3,ipw))*gprimd(1,3)
         kpg_k(ipw,2)=(kpoint(1)+kg_k(1,ipw))*gprimd(2,1)&
&         +(kpoint(2)+kg_k(2,ipw))*gprimd(2,2)&
&         +(kpoint(3)+kg_k(3,ipw))*gprimd(2,3)
         kpg_k(ipw,3)=(kpoint(1)+kg_k(1,ipw))*gprimd(3,1)&
&         +(kpoint(2)+kg_k(2,ipw))*gprimd(3,2)&
&         +(kpoint(3)+kg_k(3,ipw))*gprimd(3,3)
       end do !ipw
       kpg_k=two_pi*kpg_k
       if (dtset%nspinor==2) kpg_k(npw_k+1:2*npw_k,1:3)=kpg_k(1:npw_k,1:3)
       ABI_DEALLOCATE(kg_k)

!      Read the wavefunction block if mkmem=0
       if (mkmem==0) then
         tim_rwwf=1
         ABI_ALLOCATE(eig_dum,(dtset%mband))
         ABI_ALLOCATE(kg_dum,(3,0))
         ABI_ALLOCATE(occ_dum,(dtset%mband))
         call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,dtset%mband,mcg_disk,mpi_enreg,&
&         nband_k,nband_k,npw_k,my_nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
         ABI_DEALLOCATE(eig_dum)
         ABI_DEALLOCATE(kg_dum)
         ABI_DEALLOCATE(occ_dum)
       end if

!      2-A Computation of <psi_tild_n|-i.nabla|psi_tild_m>
!      ----------------------------------------------------------------------------------
!      Computation of (C_nk^*)*C_mk*(k+g) in cartesian coordinates

       ABI_ALLOCATE(tnm,(2,3,nband_k,nband_k))
       tnm=zero

!      Loops on bands
       do jb=1,nband_k
         jwavef=(jb-1)*npw_k*my_nspinor+icg
         if (mpi_enreg%paral_compil_kpt==1.and.mpi_enreg%mode_para/='b') then
           tmp_shape = shape(mpi_enreg%proc_distrb)
           if (ikpt > tmp_shape(1)) then
             msg='  ikpt out of bounds '
             MSG_BUG(msg)
           end if
           if (abs(mpi_enreg%proc_distrb(ikpt,jb,isppol)-me_kpt)/=0) cycle
         end if
         do ib=1,jb
           iwavef=(ib-1)*npw_k*my_nspinor+icg

!          Computation of (C_nk^*)*C_mk*(k+g) in cartesian coordinates
           if (mkmem==0) then
             if (cplex==1) then
               do ipw=1,npw_k*my_nspinor
                 cgnm1=cg_disk(1,ipw+iwavef)*cg_disk(1,ipw+jwavef)
                 tnm(1,1:3,ib,jb)=tnm(1,1:3,ib,jb)+cgnm1*kpg_k(ipw,1:3)
               end do
             else
               do ipw=1,npw_k*my_nspinor
                 cgnm1=cg_disk(1,ipw+iwavef)*cg_disk(1,ipw+jwavef)+cg_disk(2,ipw+iwavef)*cg_disk(2,ipw+jwavef)
                 cgnm2=cg_disk(1,ipw+iwavef)*cg_disk(2,ipw+jwavef)-cg_disk(2,ipw+iwavef)*cg_disk(1,ipw+jwavef)
                 tnm(1,1:3,ib,jb)=tnm(1,1:3,ib,jb)+cgnm1*kpg_k(ipw,1:3)
                 tnm(2,1:3,ib,jb)=tnm(2,1:3,ib,jb)+cgnm2*kpg_k(ipw,1:3)
               end do
             end if
           else
             if (cplex==1) then
               do ipw=1,npw_k*my_nspinor
                 cgnm1=cg(1,ipw+iwavef)*cg(1,ipw+jwavef)
                 tnm(1,1:3,ib,jb)=tnm(1,1:3,ib,jb)+cgnm1*kpg_k(ipw,1:3)
               end do
             else
               do ipw=1,npw_k*my_nspinor
                 cgnm1=cg(1,ipw+iwavef)*cg(1,ipw+jwavef)+cg(2,ipw+iwavef)*cg(2,ipw+jwavef)
                 cgnm2=cg(1,ipw+iwavef)*cg(2,ipw+jwavef)-cg(2,ipw+iwavef)*cg(1,ipw+jwavef)
                 tnm(1,1:3,ib,jb)=tnm(1,1:3,ib,jb)+cgnm1*kpg_k(ipw,1:3)
                 tnm(2,1:3,ib,jb)=tnm(2,1:3,ib,jb)+cgnm2*kpg_k(ipw,1:3)
               end do
             end if
           end if

!          Second half of the (n,m) matrix
           if (ib/=jb) then
             tnm(1,1:3,jb,ib)= tnm(1,1:3,ib,jb)
             tnm(2,1:3,jb,ib)=-tnm(2,1:3,ib,jb)
           end if

         end do ! ib
       end do ! jb

!      Reduction in case of parallelism
       if (mpi_enreg%paral_compil_fft == 1) then
         call timab(48,1,tsec)
         if (mpi_enreg%mode_para=='b') then
           call xsum_master(tnm,0,spaceComm_bandfftspin,ierr)
         else
           call xsum_master(tnm,0,spaceComm_fft,ierr)
         end if
         call timab(48,2,tsec)
       end if

       psinablapsi(:,:,:,:)=tnm(:,:,:,:)

       ABI_DEALLOCATE(tnm)

       if (mkmem/=0) then
         icg = icg + npw_k*my_nspinor*nband_k
       end if

       ABI_DEALLOCATE(kpg_k)

       if (me==0) then
         write(123)(eig0_k(ib),ib=1,nband_k)
         write(123)((psinablapsi(1:2,1,ib,jb),ib=1,nband_k),jb=1,nband_k)
         write(123)((psinablapsi(1:2,2,ib,jb),ib=1,nband_k),jb=1,nband_k)
         write(123)((psinablapsi(1:2,3,ib,jb),ib=1,nband_k),jb=1,nband_k)
       elseif (mpi_enreg%me_band==0.and.mpi_enreg%me_fft==0) then
         call xexch_mpi(psinablapsi,etiq,me_kpt,psinablapsi,0,spaceComm_k,ierr)
       end if

     elseif (me==0) then
       sender=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol))
       call xexch_mpi(psinablapsi,etiq,sender,psinablapsi,0,spaceComm_k,ierr)
       write(123)(eig0_k(ib),ib=1,nband_k)
       write(123)((psinablapsi(1:2,1,ib,jb),ib=1,nband_k),jb=1,nband_k)
       write(123)((psinablapsi(1:2,2,ib,jb),ib=1,nband_k),jb=1,nband_k)
       write(123)((psinablapsi(1:2,3,ib,jb),ib=1,nband_k),jb=1,nband_k)
     end if ! mykpt

     bdtot_index=bdtot_index+nband_k
     if (me==0)  then
       ABI_DEALLOCATE(eig0_k)
     end if
!    End loop on spin,kpt
   end do ! ikpt
 end do !isppol

!Close file
 call WffClose(wff1,ierr)

!Datastructures deallocations
 if (mkmem==0)  then
   ABI_DEALLOCATE(cg_disk)
 end if
 ABI_DEALLOCATE(psinablapsi)

 DBG_EXIT("COLL")

 end subroutine optics_vloc
!!***
