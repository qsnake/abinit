!{\src2tex{textfont=tt}}
!!****f* ABINIT/partial_dos_fractions
!! NAME
!! partial_dos_fractions
!!
!! FUNCTION
!! calculate partial DOS fractions to feed to the tetrahedron method
!!  1 : project states on angular momenta
!!  2 : should be able to choose certain atoms or atom types, slabs of space...
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MVer,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mcg)=planewave coefficients of wavefunctions
!!  dtfil = structured datatype for disk files: units etc...
!!  dtset     structured datatype, from which one uses :
!!   exchn2n3d=if 1, n2 and n3 are exchanged
!!   kpt(3,nkpt)  =irreducible kpoints
!!   kptrlatt(3,3)=lattice vectors for full kpoint grid
!!   mband        =maximum number of bands
!!   mkmem        =number of kpoints in memory
!!   natom        =number of atoms in total
!!   natsph       =number of atoms for which the spherical decomposition must be done
!!   nband        =number of electronic bands for each kpoint
!!   nkpt         =number of irreducible kpoints
!!   nshiftk      =number of kpoint grid shifts
!!   nspinor      =number of spinor components
!!   nsppol       =1 or 2 spin polarization channels
!!   nsym         =number of symmetries
!!   shiftk(3,nshiftk)=kpoint shifts
!!   symrel(3,3,nsym)=symmetry matrices in real space
!!  hdr= header of the wavefunction file (contains many informations)
!!  mbesslang=maximum angular momentum for Bessel function expansion
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=informations about MPI parallelization
!!  m_dos_flag=option for the m-contributions to the partial DOS
!!  ndosfraction=natsph*mbesslang
!!  partial_dos= option for this routine - only 1 is supported at present
!!  wffnow = eventual disk file for mkmem = 0
!!
!! OUTPUT
!!  dos_fractions(ikpt,iband,isppol,ndosfraction) = percentage of s, p, d..
!!    character on each atom for the wavefunction # ikpt,iband, isppol
!!  == if m_dos_flag==1
!!  dos_fractions_m(ikpt,iband,isppol,ndosfraction*mbesslang) = percentage of s, p, d..
!!    character on each atom for the wavefunction # ikpt,iband, isppol (m-resolved)
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      clsopn,dens_in_sph,getkpgnorm,getph,hdr_skip,init_bess_spl,initylmg
!!      kpgio,metric,ph1d3d,recip_ylm,rwwf,sort_dp,splint,xcomm_init,xredxcart
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine partial_dos_fractions(cg,dos_fractions,dos_fractions_m,dtfil,&
&           dtset,hdr,mbesslang,mcg,mpi_enreg,m_dos_flag,ndosfraction,partial_dos,wffnow)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_splines
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'partial_dos_fractions'
 use interfaces_28_numeric_noabirule
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
 use interfaces_59_io_mpi
 use interfaces_62_occeig
 use interfaces_65_nonlocal
 use interfaces_67_common, except_this_one => partial_dos_fractions
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: m_dos_flag,mbesslang,mcg,ndosfraction,partial_dos
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(wffile_type),intent(inout) :: wffnow
!arrays
 real(dp),intent(inout) :: cg(2,mcg)
 real(dp),intent(out) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
 real(dp),intent(out) :: dos_fractions_m(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang*m_dos_flag)

!Local variables-------------------------------
!scalars
 integer :: cg1kptshft,cgshift,ia,iatom,iband,ierr,ii,ikpt,ilang,ioffkg
 integer :: ioffylm,iout,ipw,ispinor,isppol,ixint,mbess,mcg_disk,me
 integer :: mgfft,my_nspinor,n1,n2,n3,natsph,nband_dum,nfit,npw_k,nradintmax,oldkpt,prtsphere
 integer :: spaceComm,tim_rwwf,unkg_dum,unylm
 real(dp) :: arg,bessarg,bessargmax,bessint_delta,kpgmax,rmax,ucvol
 character(len=4) :: mode_paral
 character(len=fnlen) :: kgnam
 type(MPI_type) :: mpi_enreg_dummy
!arrays
 integer :: atindx(dtset%natom),iindex(dtset%mpw),kg_dum(3,dtset%mpw)
 integer :: kg_k(3,dtset%mpw),npwarr_trivial(1)
 integer,allocatable :: iatsph(:),kg(:,:),npwarr1(:),npwtot(:),nradint(:)
 real(dp) :: cmax(dtset%natom),gmet(3,3)
 real(dp) :: eig_dum(dtset%mband),gprimd(3,3),kpoint(3),occ_dum(dtset%mband)
 real(dp) :: phkxred(2,dtset%natom),rmet(3,3),xcart(3,dtset%natom)
 real(dp) :: xfit(dtset%mpw),yfit(dtset%mpw)
 real(dp) :: ylm_k(dtset%mpw,mbesslang*mbesslang),ylmgr_dum(1)
 real(dp),allocatable :: bess_fit(:,:,:),bess_spl(:,:),bess_spl_der(:,:)
 real(dp),allocatable :: cg_1band(:,:),cg_1kpt(:,:),kpgnorm(:),ph1d(:,:)
 real(dp),allocatable :: ph3d(:,:,:),ratsph(:),rint(:),sum_1atom_1ll(:,:)
 real(dp),allocatable :: sum_1atom_1lm(:,:),x_bess(:),ylm(:,:)

!*************************************************************************

!for the moment, only support projection on angular momenta
 if (partial_dos /= 1) then
   write(std_out,*) 'Error : partial_dos_fractions only supports angular '
   write(std_out,*) ' momentum projection for the moment. return to outscfcv'
   write(std_out,*) ' partial_dos = ', partial_dos
   return
 end if

!impose all kpoints have same number of bands
 do isppol=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
     if (dtset%nband((isppol-1)*dtset%nkpt + ikpt) /= dtset%mband) then
       write(std_out,*) 'Error : partial_dos_fractions wants same number of',&
&       ' bands at each kpoint'
       write(std_out,*) ' isppol, ikpt = ', isppol,ikpt, &
&       dtset%nband((isppol-1)*dtset%nkpt + ikpt), dtset%mband
       write(std_out,*) ' all nband = ', dtset%nband
       return
     end if
   end do
 end do

 if (dtset%nspinor==2) then
   write(std_out,'(4a)') ch10,&
&   ' partial_dos_fractions: ERROR -',ch10,&
&   '  Cannot compute projected DOS when nspinor=2 ! '
   return
 end if

!initialize atindx
 do iatom=1,dtset%natom
   atindx(iatom) = iatom
 end do

!initialize dos_fractions
 dos_fractions(:,:,:,:) = zero
 if (m_dos_flag==1) dos_fractions_m(:,:,:,:) = zero

!initialize mpi_enreg_dummy
 mpi_enreg_dummy%paral_compil_kpt = 0
 mpi_enreg_dummy%me = 0

 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_kpt)
 me=xcomm_rank(spaceComm)

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)
 mcg_disk=dtset%mpw*my_nspinor*dtset%mband

 if (dtset%mkmem==0) then
   call clsopn(wffnow)
   call hdr_skip(wffnow,ierr)
!  Should use xdefineOff for MPI I/O
!  Define offsets, in case of MPI I/O
!  formeig=0
!  call xdefineOff(formeig,wffnow,mpi_enreg,dtset%nband,hdr%npwarr,dtset%nspinor !  &
!  &             ,dtset%nsppol,dtset%nkpt)
 end if

!##############################################################
!FIRST CASE : project on angular momenta to get dos parts
!##############################################################

 if (partial_dos == 1) then

   natsph = dtset%natsph
   ABI_ALLOCATE(iatsph,(natsph))
   ABI_ALLOCATE(ratsph,(dtset%natom))
   ABI_ALLOCATE(nradint,(natsph))
!  
!  FIXME: normally this min(size) is not useful: why should dtset%iatsph have a size different from natsph? Because of other datasets?
!  
   iatsph(1:min(natsph,size(dtset%iatsph)))=dtset%iatsph(1:min(natsph,size(dtset%iatsph)))
   do ii=1,dtset%natom
     ratsph(ii)=dtset%ratsph(dtset%typat(ii))
   end do

!  init bessel function integral for recip_ylm
!  max ang mom + 1
   ABI_ALLOCATE(sum_1atom_1ll,(mbesslang,natsph))
   ABI_ALLOCATE(sum_1atom_1lm,(mbesslang**2,natsph))
   rmax=zero
   do ii=1,natsph
     rmax=max(rmax,ratsph(iatsph(ii)))
   end do

   bessint_delta = 0.1_dp
   kpgmax = sqrt(dtset%ecut)
   rmax = zero ; bessargmax = zero ; nradintmax = 0
   do ii=1,natsph
     rmax=max(rmax,ratsph(iatsph(ii)))
     bessarg=ratsph(iatsph(ii))*two_pi*kpgmax
     bessargmax=max(bessargmax,bessarg)
     nradint(ii) = int (bessarg / bessint_delta) + 1
     nradintmax=max(nradintmax,nradint(ii))
   end do
   write(std_out,*) ' partial_dos_fractions :  rmax = ', rmax
!  use same number of grid points to calculate Bessel function
!  and to do the integration later on r
   mbess = nradintmax
!  make sure bessargmax is a multiple of bessint_delta
   bessargmax = bessint_delta*mbess

   ABI_ALLOCATE(bess_spl,(mbess,mbesslang))
   ABI_ALLOCATE(bess_spl_der,(mbess,mbesslang))
   ABI_ALLOCATE(x_bess,(nradintmax))
   ABI_ALLOCATE(rint,(nradintmax))
   ABI_ALLOCATE(bess_fit,(dtset%mpw,nradintmax,mbesslang))

!  
!  initialize general Bessel function array on uniform grid
!  x_bess, from 0 to (2 \pi |k+G|_{max} |r_{max}|)
!  
!  call init_bess_spl(mbess,bessargmax,bessint_delta,mbesslang,
   call init_bess_spl(mbess,bessint_delta,mbesslang,&
&   bess_spl,bess_spl_der,x_bess)

!  get xcart from xred and rprim which are in dtset
   call xredxcart(dtset%natom,1,hdr%rprimd,xcart,hdr%xred)

!  get recip space metric
!  if iout<0, the output of metric will not be printed
   iout=ab_out
   call metric(gmet,gprimd,iout,rmet,hdr%rprimd,ucvol)

!  get kg matrix of the positions of G vectors in recip space
   mgfft=maxval(dtset%ngfft(1:3))
!  write(std_out,*) 'DEBUG : about to allocate npwarr1 kg'

!  kg contains G vectors only for kpoints used by this processor
   mode_paral='PERS'
   ABI_ALLOCATE(kg,(3,dtset%mpw*dtset%mkmem))

!  If mkmem /= 0 fill kg array using kpgio
   if (dtset%mkmem /= 0) then
     ABI_ALLOCATE(npwarr1,(dtset%nkpt))
     ABI_ALLOCATE(npwtot,(dtset%nkpt))
     kg(:,:) = 0
!    kgnam is dummy here
     call kpgio(dtset%ecut,dtset%exchn2n3d,gmet,dtset%istwfk,kg,kgnam,dtset%kpt,&
&     dtset%mkmem,dtset%nband,dtset%nkpt,&
&     mode_paral,mpi_enreg,dtset%mpw,npwarr1,npwtot,dtset%nsppol,unkg_dum)
     ABI_DEALLOCATE(npwarr1)
     ABI_DEALLOCATE(npwtot)
   end if

!  
!  for each electronic state, get corresponding wavefunction and project on Ylm
!  real(dp) :: cg(2,dtset%mpw*my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
!  
   n1 = dtset%ngfft(1); n2 = dtset%ngfft(2); n3 = dtset%ngfft(3)

   ABI_ALLOCATE(ph1d,(2,(2*n1+1+2*n2+1+2*n3+1)*dtset%natom))
   call getph(atindx,dtset%natom,n1,n2,n3,ph1d,hdr%xred)

!  kpgnorm contains norms only for kpoints used by this processor
   if (dtset%mkmem /= 0) then
     ABI_ALLOCATE(kpgnorm,(dtset%mpw*dtset%mkmem))
!    ... or all the kpoints if mkmem==0
   else if (dtset%mkmem == 0) then
     ABI_ALLOCATE(kpgnorm,(dtset%mpw*dtset%nkpt))
   end if

!  
!  Now get Ylm factors: returns "real Ylms", which are real (+m) and
!  imaginary (-m) parts of actual complex Ylm. Yl-m = Ylm*
!  
!  Single call to initylmg for all kg (all mkmem are in memory)
!  in this call dtfil%unkg and unylm should not be used
   if (dtset%mkmem/=0) then
     ABI_ALLOCATE(ylm,(dtset%mpw*dtset%mkmem,mbesslang*mbesslang))
     call initylmg(gprimd,kg,dtset%kpt,dtset%mkmem,mpi_enreg,mbesslang,&
&     dtset%mpw,dtset%nband,dtset%nkpt,&
&     hdr%npwarr,dtset%nsppol,0,hdr%rprimd,dtfil%unkg,unylm,ylm,ylmgr_dum)
   else
     ABI_ALLOCATE(ylm,(dtset%mpw*dtset%nkpt,mbesslang*mbesslang))
   end if

   kpgnorm (:) = zero
   ioffkg=0
   ioffylm=0
   if (dtset%mkmem==0) rewind (unit=dtfil%unkg)
   do ikpt=1,dtset%nkpt
     if(mpi_enreg%paral_compil_kpt==1)then
       if (mpi_enreg%proc_distrb(ikpt,1,1)/=me) then
!        in case mkmem==0 the ylm array has full nkpt size
         if (dtset%mkmem==0) ioffylm=ioffylm+hdr%npwarr(ikpt)
         cycle
       end if
     end if

     kg_k(:,:) = 0
     if (dtset%mkmem==0) then
       read(dtfil%unkg) npw_k
       read(dtfil%unkg)
       read(dtfil%unkg) ((kg_k(ii,ipw),ii=1,3),ipw=1,npw_k)

       npwarr_trivial(1) = npw_k
!      In this case ylm still have to be calculated from disk data for kg
!      dtfil%unkg and unylm should still not be used
!      and dummy mpi_enreg forces sequential-like execution
       call initylmg(gprimd,kg_k,dtset%kpt(:,ikpt),1,mpi_enreg_dummy,mbesslang,&
&       dtset%mpw,dtset%nband,1,&
&       npwarr_trivial,dtset%nsppol,0,hdr%rprimd,dtfil%unkg,&
&       unylm,ylm_k,ylmgr_dum)
!      works with sequential mkmem0
!      ylm(ioffkg+1:ioffkg+npw_k,:) = ylm_k(1:npw_k,:)
       ylm(ioffylm+1:ioffylm+npw_k,:) = ylm_k(1:npw_k,:)

     else
       npw_k = hdr%npwarr(ikpt)
       kg_k(:,1:min(size(kg_k,2),(npw_k))) = kg(:,ioffkg+1:ioffkg+npw_k) ! check compatibility of dimensions (PMA)
     end if

     call getkpgnorm(gprimd,dtset%kpt(:,ikpt),kg_k(:,1:npw_k),&
&     kpgnorm(ioffylm+1:ioffylm+dtset%mpw),hdr%npwarr(ikpt))

     ioffkg=ioffkg+hdr%npwarr(ikpt)
     ioffylm=ioffylm+hdr%npwarr(ikpt)
   end do !ikpt

!  In case of parallel mkmem 0 the ylm need to be collected
   if (dtset%mkmem == 0) then
     call xsum_mpi(ylm,spaceComm,ierr)
     call xsum_mpi(kpgnorm,spaceComm,ierr)
   end if

   ABI_ALLOCATE(cg_1kpt,(2,mcg_disk))

   cgshift = 0
   oldkpt = 0

   do isppol=1,dtset%nsppol
     ioffkg = 0
     ioffylm = 0
!    kg array is the same for both sppol
     if (dtset%mkmem==0) rewind (unit=dtfil%unkg)

     do ikpt=1,dtset%nkpt
       if(mpi_enreg%paral_compil_kpt==1)then
         if (mpi_enreg%proc_distrb(ikpt,1,isppol)/=mpi_enreg%me) then
           if (dtset%mkmem==0) ioffylm=ioffylm+hdr%npwarr(ikpt)
           cycle
         end if
       end if

       npw_k = hdr%npwarr(ikpt)
       ABI_ALLOCATE(cg_1band,(2,npw_k))
       kpoint(:) = dtset%kpt(:,ikpt)

!      for each kpoint set up the phase factors, ylm factors
       if (dtset%mkmem==0) then
         read(dtfil%unkg) npw_k
         read(dtfil%unkg)
         read(dtfil%unkg) ((kg_k(ii,ipw),ii=1,3),ipw=1,npw_k)
       else
         npw_k = hdr%npwarr(ikpt)
         kg_k(:,1:min(size(kg_k,2),(npw_k))) = kg(:,ioffkg+1:ioffkg+npw_k) !check dimensions compatibility
       end if

       do ilang=1,mbesslang*mbesslang
         do ipw=1,npw_k
           ylm_k(ipw,ilang) = ylm(ioffylm+ipw,ilang)
         end do
       end do

!      make phkred for all atoms
       do ia=1,dtset%natom
         arg=two_pi*( kpoint(1)*hdr%xred(1,ia) &
&         + kpoint(2)*hdr%xred(2,ia) &
&         + kpoint(3)*hdr%xred(3,ia) )
         phkxred(1,ia)=cos(arg)
         phkxred(2,ia)=sin(arg)
       end do

       ABI_ALLOCATE(ph3d,(2,npw_k,dtset%natom))

!      need simple 3d phases for dens_in_sph
       call ph1d3d(1,dtset%natom,kg_k(:,1:npw_k),&
&       dtset%natom,dtset%natom,npw_k,n1,n2,n3,&
&       phkxred,ph1d,ph3d)
!      phases exp (2 pi i G.x_tau) are now in ph3d

       if (dtset%mkmem == 0) then
         tim_rwwf = 0
         call rwwf(cg_1kpt,eig_dum,0,0,0,ikpt,isppol,kg_dum(:,:),dtset%mband,mcg_disk,&
&         mpi_enreg,dtset%mband,nband_dum,npw_k,my_nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
       else
         cg_1kpt(:,:) = cg(:,cgshift+1:cgshift+mcg_disk)
       end if

       if (dtset%prtvol > 3) write(std_out,*) 'get dens in sphere for ikpt,isppol = ', ikpt,isppol
       cg1kptshft = 0
       do iband=1,dtset%mband
         if (dtset%prtvol > 3) write(std_out,*) 'get dens in sphere for iband = ', iband
         do ispinor=1,my_nspinor

           cg_1band(:,:) = cg_1kpt(:,cg1kptshft+1:cg1kptshft+npw_k)

           call dens_in_sph(cmax,cg_1band,&
&           gmet,dtset%istwfk(ikpt),&
&           kg_k(:,1:npw_k),dtset%natom,dtset%ngfft,mpi_enreg,&
&           npw_k,dtset%paral_kgb,ph1d,ratsph,ucvol)
           cg1kptshft = cg1kptshft + npw_k
         end do
       end do

!      get full phases for the following
       call ph1d3d(1,dtset%natom,kg_k(:,1:npw_k),&
&       dtset%natom,dtset%natom,npw_k,n1,n2,n3,&
&       phkxred,ph1d,ph3d)
!      phases exp (2 pi i (k+G).x_tau) are now in ph3d

!      get Bessel function factors on array of |k+G|*r distances
!      since we need many r distances and have a large number of different
!      |k+G|, get j_l on uniform grid (above, in array gen_besj),
!      and spline it for each kpt Gvector set.
       nfit = npw_k
       do ixint=1,nradintmax
         rint(ixint) = (ixint-1)*rmax / (nradintmax-1)
         do ipw=1,npw_k
           xfit(ipw) = two_pi * kpgnorm(ipw+ioffylm) * rint(ixint)
           iindex(ipw) = ipw
         end do

         call sort_dp(npw_k,xfit,iindex,tol14)
         do ilang=1,mbesslang
           call splint(mbess,x_bess,bess_spl(:,ilang),bess_spl_der(:,ilang),&
&           nfit,xfit,yfit)
!          re-order results for different G vectors
           do ipw=1,npw_k
             bess_fit(iindex(ipw),ixint,ilang) = yfit(ipw)
           end do
         end do
       end do

       cg1kptshft=0
       do iband=1,dtset%mband
         do ispinor=1,my_nspinor

           cg_1band(:,:) = cg_1kpt(:,cg1kptshft+1:cg1kptshft+npw_k)

           prtsphere=0  ! do not output all the band by band details for projections.
           call recip_ylm (bess_fit,&
&           cg_1band,&
&           iatsph,dtset%istwfk(ikpt),&
&           nradint,nradintmax,mbesslang,mpi_enreg,dtset%mpw,dtset%natom,&
&           natsph,npw_k,dtset%ntypat,&
&           ph3d,prtsphere,rint,ratsph,sum_1atom_1ll,sum_1atom_1lm,&
&           dtset%typat,ucvol,&
&           ylm_k,hdr%znuclpsp)

           do iatom=1,natsph
             do ilang=1,mbesslang
               dos_fractions(ikpt,iband,isppol,mbesslang*(iatom-1) + ilang) &
&               = sum_1atom_1ll(ilang,iatom)
             end do
           end do

           if (m_dos_flag==1) then
             do iatom=1,natsph
               do ilang=1,mbesslang**2
                 dos_fractions_m(ikpt,iband,isppol,mbesslang**2*(iatom-1) + ilang) &
&                 = sum_1atom_1lm(ilang,iatom)
               end do
             end do
           end if

           cg1kptshft=cg1kptshft + npw_k

         end do ! end spinor
       end do ! end band
!      cgshift=cgshift + cg1kptshft
       cgshift=cgshift + dtset%mband*my_nspinor*npw_k

       ABI_DEALLOCATE(ph3d)
       ioffkg = ioffkg + npw_k
       ioffylm = ioffylm + npw_k

       ABI_DEALLOCATE(cg_1band)
     end do ! end kpt
   end do ! end sppol
   ABI_DEALLOCATE(cg_1kpt)

!  Gather all contributions from different processors
   call xsum_mpi(dos_fractions,spaceComm,ierr)
   if (m_dos_flag==1) call xsum_mpi(dos_fractions_m,spaceComm,ierr)
   if (mpi_enreg%paral_spin == 1)then
     call xsum_mpi(dos_fractions,mpi_enreg%comm_spin,ierr)
     if (m_dos_flag==1) call xsum_mpi(dos_fractions_m,mpi_enreg%comm_spin,ierr)
   end if

   ABI_DEALLOCATE(iatsph)
   ABI_DEALLOCATE(nradint)
   ABI_DEALLOCATE(ratsph)
   ABI_DEALLOCATE(sum_1atom_1ll)
   ABI_DEALLOCATE(sum_1atom_1lm)
   ABI_DEALLOCATE(bess_spl)
   ABI_DEALLOCATE(bess_spl_der)
   ABI_DEALLOCATE(x_bess)
   ABI_DEALLOCATE(bess_fit)
   ABI_DEALLOCATE(rint)
   ABI_DEALLOCATE(kg)
   ABI_DEALLOCATE(ph1d)
   ABI_DEALLOCATE(kpgnorm)
   ABI_DEALLOCATE(ylm)

 else
   write(std_out,*) ' partial_dos_fractions: only partial_dos==1 is coded '
 end if


end subroutine partial_dos_fractions
!!***
