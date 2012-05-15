!{\src2tex{textfont=tt}}
!!****f* ABINIT/optics_paw
!! NAME
!! optics_paw
!!
!! FUNCTION
!! Compute matrix elements need for optical conductivity (in the PAW context) and store them in a file
!!  Matrix elements = <Phi_i|Nabla|Phi_j>
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (SM,VR,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mcg)=planewave coefficients of wavefunctions.
!!  cprj(natom,mcprj)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  dimcprj(natom)=array of dimensions of array cprj (not ordered)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  indlmn(6,lmnmax,ntypat)= array giving l,m,n,lm,ln,s for i=lmn
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  lmnmax=max. number of (l,m,n) numbers over all types of atom
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang =1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!     %mesh_size=Dimension of radial mesh
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  wffnow=struct infos for wf disk file
!!
!! OUTPUT
!!  (only writing in a file)
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      clsopn,cprj_alloc,cprj_diskinit_r,cprj_free,cprj_get,deducer0,hdr_io
!!      hdr_skip,int_ang,leave_new,nderiv_gen,rdnpw,rwwf,simp_gen,timab
!!      wffclose,wffopen,wrtout,xcomm_init,xcomm_world,xdefineoff,xexch_mpi
!!      xme_init,xsum_master,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine optics_paw(atindx1,cg,cprj,dimcprj,dtfil,dtset,eigen0,gprimd,hdr,indlmn,kg,lmnmax,&
&               mband,mcg,mcprj,mkmem,mpi_enreg,mpsang,mpw,natom,nkpt,npwarr,nsppol,&
&               pawrad,pawtab,wffnow)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_wffile

 use m_radmesh,   only : nderiv_gen, deducer0, simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'optics_paw'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_66_paw, except_this_one => optics_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmnmax,mband,mcg,mcprj,mkmem,mpsang,mpw,natom,nkpt,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(wffile_type),intent(inout) :: wffnow
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom),indlmn(6,lmnmax,dtset%ntypat)
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt)
 real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(inout) :: cg(2,mcg)
 type(cprj_type) :: cprj(natom,mcprj)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: accesswff,basis_size,bdtot_index,cplex,etiq,fformopt,formeig,iatom,ib,ibg,ibsp
 integer :: icg,ierr,ij_size,ikg,ikpt,il,ilm,ilmn,iln
 integer :: iorder_cprj,ipw,ispinor,isppol,istwf_k,itypat,iwavef
 integer :: jb,jbsp,jl,jlm,jlmn,jln,jwavef,lmn_size
 integer :: mcg_disk,me,me_kpt,mesh_size,my_nspinor,nband_k,npw_k,nspinor0,old_paral_level,sender
 integer :: spaceComm_band,spaceComm_bandfftspin,spaceComm_fft,spaceComm_k,spaceComm_spin,spaceComm_w
 integer :: tim_rwwf
 logical :: mykpt
 real(dp) :: cgnm1,cgnm2,cpnm1,cpnm2,intg
 character(len=500) :: message
!arrays
 integer :: tmp_shape(3)
 integer,allocatable :: indlmn_(:,:),kg_dum(:,:),kg_k(:,:)
 real(dp) :: ang_phipphj(mpsang**2,mpsang**2,8),kpoint(3)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cg_disk(:,:),dphi(:),dtphi(:),eig_dum(:)
 real(dp),allocatable ::ff(:),int1(:,:),int2(:,:)
 real(dp),allocatable :: kpg_k(:,:),occ_dum(:),phidphj(:,:),psinablapsi(:,:,:,:)
 real(dp),allocatable :: rad(:),tnm(:,:,:,:),tphidtphj(:,:)
 type(coeff3_type), allocatable :: phipphj(:)
 type(cprj_type),allocatable :: cprj_k(:,:)
 type(wffile_type) :: wff1

! ************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 if(dtset%paral_kgb==1.and.mkmem==0)then
   write(message, '(4a)' )ch10,&
&   ' optics_paw :  -',ch10,&
&   '  Not compatible with mkmem=0 and band-fft parallelism !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!----------------------------------------------------------------------------------
!1- Computation of phipphj=<phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j>
!----------------------------------------------------------------------------------

!1-A Integration of the angular part : all angular integrals have been
!computed outside Abinit and tabulated for each (l,m) value
!----------------------------------------------------------------------------------

 call int_ang(ang_phipphj,mpsang)

 ABI_ALLOCATE(phipphj,(dtset%ntypat))

!loop on atoms type
 do itypat=1,dtset%ntypat

   mesh_size=pawrad(itypat)%mesh_size
   lmn_size=pawtab(itypat)%lmn_size
   basis_size=pawtab(itypat)%basis_size
   ij_size=lmn_size*lmn_size

   ABI_ALLOCATE(indlmn_,(6,lmnmax))
   ABI_ALLOCATE(ff,(mesh_size))
   ABI_ALLOCATE(rad,(mesh_size))
   ABI_ALLOCATE(int2,(lmn_size,lmn_size))
   ABI_ALLOCATE(int1,(lmn_size,lmn_size))
   ABI_ALLOCATE(dphi,(mesh_size))
   ABI_ALLOCATE(dtphi,(mesh_size))
   ABI_ALLOCATE(phidphj,(mesh_size,ij_size))
   ABI_ALLOCATE(tphidtphj,(mesh_size,ij_size))
   ABI_ALLOCATE(phipphj(itypat)%value,(3,lmn_size,lmn_size))

   indlmn_(:,:)=indlmn(:,:,itypat)
   rad(1:mesh_size)=pawrad(itypat)%rad(1:mesh_size)

!  1-B  Computation of int1=\int phi phj /r dr - \int tphi tphj /r dr
!  ----------------------------------------------------------------------------------
   do jln=1,basis_size
     do iln=1,basis_size
       ff(2:mesh_size)=(pawtab(itypat)%phi(2:mesh_size,iln)*pawtab(itypat)%phi(2:mesh_size,jln)&
&       -pawtab(itypat)%tphi(2:mesh_size,iln)*pawtab(itypat)%tphi(2:mesh_size,jln))/rad(2:mesh_size)
       call deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int1(iln,jln)=intg
     end do
   end do

!  1-C Computation of int2=\int phi/r d/dr(phj/r) r^2dr - \int tphi/r d/dr(tphj/r)r^2 dr
!  ----------------------------------------------------------------------------------
   do jln=1,basis_size

     ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,jln)
     call nderiv_gen(dphi,ff,1,pawrad(itypat))
     ff(1:mesh_size)=pawtab(itypat)%tphi(1:mesh_size,jln)
     call nderiv_gen(dtphi,ff,1,pawrad(itypat))

     do iln=1,basis_size
       ff(2:mesh_size)=pawtab(itypat)%phi(2:mesh_size,iln)*dphi(2:mesh_size) &
&       -pawtab(itypat)%phi (2:mesh_size,iln)*pawtab(itypat)%phi(2:mesh_size,jln)/ &
&       rad(2:mesh_size)-(pawtab(itypat)%tphi(2:mesh_size,iln)*dtphi(2:mesh_size) &
&       -pawtab(itypat)%tphi (2:mesh_size,iln)*pawtab(itypat)%tphi(2:mesh_size,jln)/rad(2:mesh_size))
       call deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       int2(iln,jln)=intg
     end do
   end do

!  1-D Integration of the radial part
!  ----------------------------------------------------------------------------------
   do jlmn=1,lmn_size
     jlm=indlmn_(4,jlmn)
     jl=indlmn_(5,jlmn)
     do ilmn=1,lmn_size
       ilm=indlmn_(4,ilmn)
       il=indlmn_(5,ilmn)
       phipphj(itypat)%value(1,ilmn,jlmn)= int2(il,jl)*ang_phipphj(ilm,jlm,1)&
&       + int1(il,jl)*(ang_phipphj(ilm,jlm,2)+ang_phipphj(ilm,jlm,3))
       phipphj(itypat)%value(2,ilmn,jlmn)= int2(il,jl)*ang_phipphj(ilm,jlm,4)&
&       + int1(il,jl)*(ang_phipphj(ilm,jlm,5)+ang_phipphj(ilm,jlm,6))
       phipphj(itypat)%value(3,ilmn,jlmn)= int2(il,jl)*ang_phipphj(ilm,jlm,7)&
&       + int1(il,jl)*ang_phipphj(ilm,jlm,8)
     end do
   end do

   ABI_DEALLOCATE(indlmn_)
   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(rad)
   ABI_DEALLOCATE(int2)
   ABI_DEALLOCATE(int1)
   ABI_DEALLOCATE(dphi)
   ABI_DEALLOCATE(dtphi)
   ABI_DEALLOCATE(phidphj)
   ABI_DEALLOCATE(tphidtphj)

!  end loop on atoms type
 end do

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
   spaceComm_spin=mpi_enreg%comm_spin
   spaceComm_bandfftspin=mpi_enreg%commcart_3d
 else
   spaceComm_band=0;spaceComm_fft=0;spaceComm_spin=xmpi_self;spaceComm_bandfftspin=0
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
!PAW file
 iorder_cprj=0
 call cprj_diskinit_r(atindx1,natom,iorder_cprj,mkmem,natom,0,dimcprj,my_nspinor,dtfil%unpaw)

!Initialize main variables
 ABI_ALLOCATE(psinablapsi,(2,3,mband,mband))
 psinablapsi=zero

!open _OPT file for proc 0
 accesswff=-1
 fformopt=610
 call WffOpen(accesswff,spaceComm_k,dtfil%fnameabo_app_opt,ierr,wff1,0,me,123)
 call hdr_io(fformopt,hdr,2,wff1)
 if (me==0) then
   write(123)(eigen0(ib),ib=1,mband*nkpt*nsppol)
 end if

!LOOP OVER SPINS
 ibg=0;icg=0
 bdtot_index=0
 do isppol=1,nsppol

!  LOOP OVER k POINTS
   ikg=0
   if (dtset%mkmem==0) rewind dtfil%unkg
   do ikpt=1,nkpt

     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
     etiq=ikpt+(isppol-1)*nkpt

!    Old FFT parallelism: define FFT communicator for this k-point
     if (mpi_enreg%paral_compil_fft==1.and.mpi_enreg%mode_para/='b') then
       mpi_enreg%num_group_fft=ikpt+(isppol-1)*nkpt
       old_paral_level=mpi_enreg%paral_level;mpi_enreg%paral_level=3
       call xcomm_init(mpi_enreg,spaceComm_fft)
       mpi_enreg%paral_level=old_paral_level
     end if

!    Select k-points for current proc
     mykpt=.true.
     if (mpi_enreg%paral_compil_kpt==1) then
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
       ABI_ALLOCATE(cprj_k,(natom,my_nspinor*nband_k))
       call cprj_alloc(cprj_k,0,dimcprj)

!      Extract cprj for this k-point according to mkmem
       call cprj_get(atindx1,cprj_k,cprj,natom,1,ibg,ikpt,iorder_cprj,isppol,&
&       mband,mkmem,mpi_enreg,natom,nband_k,nband_k,my_nspinor,nsppol,dtfil%unpaw)

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
             write(message, '(4a)' )ch10,&
&             ' optics_paw :  -',ch10,' error: ikpt out of bounds '
             call wrtout(std_out,message,'COLL')
             call leave_new('COLL')
           end if
!          MJV 6/12/2008: looks like mpi_enreg may not be completely initialized here
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

!      ATTEMPT: USING BLAS (not successful)
!      allocate(dg(2,npw_k*my_nspinor,nband_k),tnm_tmp(2,nband_k,nband_k))
!      do ii=1,3
!      do ib=1,nband_k
!      iwavef=icg+(ib-1)*npw_k*my_nspinor
!      do ipw=1,my_nspinor*npw_k
!      dg(:,ipw,nband_k)=cg(:,ipw+iwavef)*kpg_k(ipw,ii)
!      end do
!      end do
!      call zgemm('C','N',nband_k,nband_k,npw_k*my_nspinor,one,dg,npw_k*my_nspinor,&
!      &                cg(:,1+icg:npw_k*my_nspinor*nband_k+icg),npw_k*my_nspinor,zero,tnm_tmp,nband_k)
!      tnm(:,ii,:,:)=tnm_tmp(:,:,:)
!      deallocate(dg,tnm_tmp)
!      end do

!      Reduction in case of parallelism
       if (mpi_enreg%paral_compil_fft == 1) then
         call timab(48,1,tsec)
         if (mpi_enreg%mode_para=='b') then
           call xsum_master(tnm,0,spaceComm_bandfftspin,ierr)
         else
           call xsum_master(tnm,0,spaceComm_fft,ierr)
!          In that case parallelisation on nspinor not implemented yet;put this line just in case
           call xsum_master(tnm,0,mpi_enreg%comm_spin,ierr)
         end if
         call timab(48,2,tsec)
       end if

       psinablapsi(:,:,:,:)=tnm(:,:,:,:)

!      2-B Computation of <psi_n|p_i><p_j|psi_m>(<phi_i|-i.nabla|phi_j>-<tphi_i|-i.nabla|tphi_j>)
!      ----------------------------------------------------------------------------------

       tnm=zero

!      Loops on bands
       do jb=1,nband_k

         if (mpi_enreg%mode_para=='b') then
           if (mod(jb-1,mpi_enreg%nproc_band)/=mpi_enreg%me_band) cycle
         elseif (mpi_enreg%paral_compil_kpt==1) then
           if (abs(mpi_enreg%proc_distrb(ikpt,jb,isppol)-me_kpt)/=0) cycle
         end if

         do ib=1,jb
           ibsp=(ib-1)*my_nspinor;jbsp=(jb-1)*my_nspinor

           if (cplex==1) then
             do ispinor=1,my_nspinor
               ibsp=ibsp+1;jbsp=jbsp+1
               do iatom=1,natom
                 itypat=dtset%typat(iatom)
                 lmn_size=pawtab(itypat)%lmn_size
                 do jlmn=1,lmn_size
                   do ilmn=1,lmn_size
                     cpnm1=cprj_k(iatom,ibsp)%cp(1,ilmn)*cprj_k(iatom,jbsp)%cp(1,jlmn)
                     tnm(2,:,ib,jb)=tnm(2,:,ib,jb)+cpnm1*phipphj(itypat)%value(:,ilmn,jlmn)
                   end do !ilmn
                 end do !jlmn
               end do !iatom
             end do !ispinor
           else
             do ispinor=1,my_nspinor
               ibsp=ibsp+1;jbsp=jbsp+1
               do iatom=1,natom
                 itypat=dtset%typat(iatom)
                 lmn_size=pawtab(itypat)%lmn_size
                 do jlmn=1,lmn_size
                   do ilmn=1,lmn_size
                     cpnm1=(cprj_k(iatom,ibsp)%cp(1,ilmn)*cprj_k(iatom,jbsp)%cp(1,jlmn) &
&                     +cprj_k(iatom,ibsp)%cp(2,ilmn)*cprj_k(iatom,jbsp)%cp(2,jlmn))
                     cpnm2=(cprj_k(iatom,ibsp)%cp(1,ilmn)*cprj_k(iatom,jbsp)%cp(2,jlmn) &
&                     -cprj_k(iatom,ibsp)%cp(2,ilmn)*cprj_k(iatom,jbsp)%cp(1,jlmn))
                     tnm(1,:,ib,jb)=tnm(1,:,ib,jb)+cpnm2*phipphj(itypat)%value(:,ilmn,jlmn)
                     tnm(2,:,ib,jb)=tnm(2,:,ib,jb)-cpnm1*phipphj(itypat)%value(:,ilmn,jlmn)
                   end do !ilmn
                 end do !jlmn
               end do !iatom
             end do !ispinor
           end if

!          Second half of the (n,m) matrix
           if (ib/=jb) then
             tnm(1,1:3,jb,ib)=-tnm(1,1:3,ib,jb)
             tnm(2,1:3,jb,ib)= tnm(2,1:3,ib,jb)
           end if

!          End loops on bands
         end do ! jb
       end do !ib

!      Reduction in case of parallelism
       if (mpi_enreg%mode_para=='b') then
         call timab(48,1,tsec)
         call xsum_mpi(tnm,mpi_enreg%comm_bandspin,ierr)
         call timab(48,2,tsec)
       else
!        In that case parallelisation on nspinor not implemented yet;put this line just in case
         call xsum_mpi(tnm,spacecomm_spin,ierr)
       end if

       psinablapsi(:,:,:,:)=psinablapsi(:,:,:,:)+tnm(:,:,:,:)
       ABI_DEALLOCATE(tnm)

       if (mkmem/=0) then
         ibg = ibg +       my_nspinor*nband_k
         icg = icg + npw_k*my_nspinor*nband_k
       end if

       call cprj_free(cprj_k)
       ABI_DEALLOCATE(kpg_k)
       ABI_DEALLOCATE(cprj_k)

       if (me==0) then
         write(123)((psinablapsi(1:2,1,ib,jb),ib=1,nband_k),jb=1,nband_k)
         write(123)((psinablapsi(1:2,2,ib,jb),ib=1,nband_k),jb=1,nband_k)
         write(123)((psinablapsi(1:2,3,ib,jb),ib=1,nband_k),jb=1,nband_k)
       elseif (mpi_enreg%me_band==0.and.mpi_enreg%me_fft==0) then
         call xexch_mpi(psinablapsi,etiq,me_kpt,psinablapsi,0,spaceComm_k,ierr)
       end if

     elseif (me==0) then
       sender=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol))
       call xexch_mpi(psinablapsi,etiq,sender,psinablapsi,0,spaceComm_k,ierr)
       write(123)((psinablapsi(1:2,1,ib,jb),ib=1,nband_k),jb=1,nband_k)
       write(123)((psinablapsi(1:2,2,ib,jb),ib=1,nband_k),jb=1,nband_k)
       write(123)((psinablapsi(1:2,3,ib,jb),ib=1,nband_k),jb=1,nband_k)
     end if ! mykpt

     bdtot_index=bdtot_index+nband_k
!    End loop on spin,kpt
   end do ! ikpt
 end do !isppol

!Close file
 call WffClose(wff1,ierr)

!Datastructures deallocations
 if (mkmem==0)  then
   ABI_DEALLOCATE(cg_disk)
 end if
 do itypat=1,dtset%ntypat
   ABI_DEALLOCATE(phipphj(itypat)%value)
 end do
 ABI_DEALLOCATE(phipphj)
 ABI_DEALLOCATE(psinablapsi)

 DBG_EXIT("COLL")

 end subroutine optics_paw
!!***
