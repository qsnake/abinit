!{\src2tex{textfont=tt}}
!!****f* ABINIT/gstateimg
!! NAME
!! gstateimg
!!
!! FUNCTION
!! Routine for conducting DFT calculations for a set of (dynamical) images
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG, AR, GG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  codvsn=code version
!!  cpui=initial CPU time
!!  filnam(5)=character strings giving file names
!!  filstat=character strings giving name of status file
!!  nimage=number of images of the cell (treated by current proc)
!!  === Optional arguments (needed when nimage>1) ===
!!    filnam(5)=character strings giving file names
!!    filstat=character strings giving name of status file
!!    idtset=index of the dataset
!!    jdtset(0:ndtset)=actual index of the datasets
!!    ndtset=number of datasets
!!
!! OUTPUT
!!  etotal_img=total energy, for each image
!!  fcart_img(3,natom,nimage)=forces, in cartesian coordinates, for each image
!!  fred_img(3,natom,nimage)=forces, in reduced coordinates, for each image
!!  npwtot(nkpt) = total number of plane waves at each k point
!!  strten_img(6,nimage)=stress tensor, for each image
!!
!! SIDE EFFECTS
!!  acell_img(3,nimage)=unit cell length scales (bohr), for each image
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | mband =maximum number of bands (IN)
!!   | mgfft =maximum single fft dimension (IN)
!!   | mkmem =maximum number of k points which can fit in core memory (IN)
!!   | mpw   =maximum number of planewaves in basis sphere (large number) (IN)
!!   | natom =number of atoms in unit cell (IN)
!!   | nfft  =(effective) number of FFT grid points (for this processor) (IN)
!!   | nkpt  =number of k points (IN)
!!   | nspden=number of spin-density components (IN)
!!   | nsppol=number of channels for spin-polarization (1 or 2) (IN)
!!   | nsym  =number of symmetry elements in space group
!!  iexit= exit flag
!!  mpi_enreg=MPI-parallelisation information (some already initialized,
!!            some others to be initialized here)
!!  occ_img(mband*nkpt*nsppol,nimage) = occupation number for each band and k, for each image
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   Before entering the first time in gstateimg, a significant part of
!!   psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,
!!     ntypat,n1xccc,usepaw,useylm, and the arrays dimensioned to npsp
!!   All the remaining components of psps are to be initialized in the call
!!   to pspini .
!!   The next time the code enters gstateimg, psps might be identical to the
!!   one of the previous dtset, in which case, no reinitialisation is scheduled
!!   in pspini.f .
!!  rprim_img(3,3,nimage)=dimensionless real space primitive translations, for each image
!!  vel_img(3,natom,nimage)=value of velocity,for each image
!!  xred_img(3,natom,nimage) = reduced atomic coordinates, for each image
!!
!! NOTES
!! USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!! In case of wavelets:
!! --------------------
!!    - Only the usual FFT grid (defined by wvl_crmult) is used.
!!      It is defined by nfft, ngfft, mgfft, ... This is strictly not
!!      an FFT grid since its dimensions are not suited for FFTs. They are
!!      defined by wvl_setngfft().
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! TODO
!! Not yet possible to use restartxf in parallel when localrdwf==0
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      abi_io_redirect,appdig,destroy_results_img,destroy_scf_history
!!      dtfil_init1,flush_unit,gstate,init_results_img,int2char4,leave_test
!!      nullify_scf_history,predictimg,prtimg,status,timab,wrtout,wrtout_myproc
!!      xbarrier_mpi,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine gstateimg(acell_img,codvsn,cpui,dtfil,dtset,etotal_img,fcart_img,&
&                    fred_img,iexit,mpi_enreg,nimage,npwtot,occ_img,&
&                    pawang,pawrad,pawtab,psps,rprim_img,strten_img,vel_img,xred_img,&
&                    filnam,filstat,idtset,jdtset,ndtset) ! optional arguments

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_scftypes
 use defs_wvltypes
 use defs_parameters
 use defs_rectypes
 use m_pimd
 use m_xmpi
 use m_errors
 use m_wffile
 use m_rec
 use m_results_img
 use m_scf_history
 use m_io_tools,         only : flush_unit
 use m_paw_dmft,         only : init_sc_dmft,destroy_sc_dmft,print_sc_dmft,paw_dmft_type
 use m_paw_toolbox,      only : init_pawfgr
 use m_electronpositron, only : electronpositron_type,init_electronpositron,destroy_electronpositron, &
&                               electronpositron_calctype
 use m_header,           only : hdr_init, hdr_clean
 use m_ebands,           only : bstruct_init, bstruct_clean
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gstateimg'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_27_toolbox_oop
 use interfaces_32_util
 use interfaces_45_geomoptim
 use interfaces_51_manage_mpi
 use interfaces_67_common
 use interfaces_95_drive, except_this_one => gstateimg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nimage
 integer,optional,intent(in) :: idtset,ndtset
 integer,intent(inout) :: iexit
 real(dp),intent(in) :: cpui
 character(len=6),intent(in) :: codvsn
 character(len=fnlen),optional,intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),target,intent(inout) :: dtfil
 type(dataset_type),target,intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,optional,intent(in) :: jdtset(:)
 integer,intent(out) :: npwtot(dtset%nkpt)
 character(len=fnlen),optional,intent(in) :: filnam(:)
 real(dp), intent(out) :: etotal_img(nimage),fcart_img(3,dtset%natom,nimage)
 real(dp), intent(out) :: fred_img(3,dtset%natom,nimage),strten_img(6,nimage)
 real(dp),intent(inout) :: acell_img(3,nimage),occ_img(dtset%mband*dtset%nkpt*dtset%nsppol,nimage)
 real(dp),intent(inout) :: rprim_img(3,3,nimage),vel_img(3,dtset%natom,nimage)
 real(dp),intent(inout) :: xred_img(3,dtset%natom,nimage)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!Define file format for different type of files. Presently,
!only one file format is supported for each type of files, but this might
!change soon ...
!2   for wavefunction file, new format (version 2.0 and after)    (fform)   NOT USED
!52  for density rho(r)       (fformr)
!102 for potential V(r) file. (fformv)  NOT USED
!scalars
 integer,parameter :: formeig=0,level=100,ndtpawuj=0,response=0
 integer :: history_size,idynimage,ierr,ii,iimage,ios,itimimage,itimimage_eff
 integer :: last_itimimage,ndynimage,nocc,ntimimage_eff,slideimg
 logical :: eof,ex,compute_all_images,compute_static_images
 real(dp) :: delta_energy
 character(len=4) :: appen,tag
 character(len=500) :: msg
 character(len=fnlen) :: line
 type(pimd_type) :: pimd_param
!arrays
 integer,allocatable :: list_dynimage(:),scf_initialized(:)
 character(len=60),parameter :: imagealgo_str(0:13)=(/ &
&   'IMAGE COPY                                                  ',& ! 0
&   'IMAGE STEEPEST DESCENT                                      ',& ! 1
&   'STRING METHOD                                               ',& ! 2
&   'METADYNAMICS                                                ',& ! 3
&   'GENETIC ALGORITHM                                           ',& ! 4
&   '                                                            ',& ! 5
&   '                                                            ',& ! 6
&   '                                                            ',& ! 7
&   '                                                            ',& ! 8
&   'PATH-INTEGRAL MOLECULAR DYNAMICS (LANGEVIN)                 ',& ! 9
&   '                                                            ',& ! 10
&   '                                                            ',& ! 11
&   '                                                            ',& ! 12
&   'PATH-INTEGRAL MOLECULAR DYNAMICS (CHAIN OF THERMOSTATS)     '/) ! 13
 character(len=fnlen),allocatable :: fillog(:),filout(:)
 real(dp) :: acell(3),rprim(3,3),tsec(2)
 real(dp),allocatable :: occ(:),vel(:,:),xred(:,:)
 type(results_img_type),pointer :: results_img_timimage(:,:)
 type(scf_history_type),allocatable :: scf_history(:)

! ***********************************************************************

 call timab(700,1,tsec)
 call timab(703,3,tsec)

 DBG_ENTER("COLL")

!Arguments check
 if (dtset%nimage>1) then
   if ((.not.present(filnam)).or.(.not.present(filnam)).or.(.not.present(idtset)).or.&
&   (.not.present(ndtset)).or.(.not.present(jdtset))) then
     write(msg,'()') &
&     '  When nimage>1, all the following argument should be present:',ch10,&
&     '        filnam, filstat, idtset, ndtset, jdtset  !'
     MSG_BUG(msg)
   end if
 end if

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Set flag for the effective computation (once) of static images
!For the time being only set when parallelization is activated
!Note: if you modify this flag, do not forget to change it in outvars and outvar1
 compute_static_images=(dtset%istatimg>0)

!Prepare the allocations, by computing flags and dimensions
 slideimg=0;if(dtset%imgmov==9 .or. dtset%imgmov==13)slideimg=1
 ntimimage_eff=dtset%ntimimage;if(slideimg==1)ntimimage_eff=2
 nocc=dtset%mband*dtset%nkpt*dtset%nsppol

!Allocations
 ABI_ALLOCATE(occ,(nocc))
 ABI_ALLOCATE(vel,(3,dtset%natom))
 ABI_ALLOCATE(xred,(3,dtset%natom))
 ABI_ALLOCATE(results_img_timimage,(nimage,ntimimage_eff))
 ABI_ALLOCATE(list_dynimage,(dtset%ndynimage))
 do itimimage=1,ntimimage_eff
   call init_results_img(dtset%natom,results_img_timimage(:,itimimage))
   do iimage=1,nimage
     results_img_timimage(iimage,itimimage)%acell(:)  =acell_img(:,iimage)
     results_img_timimage(iimage,itimimage)%rprim(:,:)=rprim_img(:,:,iimage)
     results_img_timimage(iimage,itimimage)%xred(:,:) =xred_img(:,:,iimage)
     results_img_timimage(iimage,itimimage)%vel(:,:)  =vel_img(:,:,iimage)
   end do
 end do
 ndynimage=0
 do iimage=1,nimage
   ii=mpi_enreg%index_img(iimage)
   if (dtset%dynimage(ii)==1) then
     ndynimage=ndynimage+1
     list_dynimage(ndynimage)=iimage
   end if
 end do

!Management of SCF history (density/WF predictions from one time step to another)
 ABI_ALLOCATE(scf_history,(nimage))
 ABI_ALLOCATE(scf_initialized,(nimage))
 scf_initialized=0
 history_size=-1
 if (dtset%ntimimage<=1) then
   if (dtset%usewvl==0.and.dtset%ionmov>0.and. &
&   (abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6)) history_size=2
 else
   if (abs(dtset%iprcch)==2.or.abs(dtset%iprcch)==3) history_size=0
   if (dtset%usewvl==0.and.(abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6)) history_size=2
 end if
 do iimage=1,nimage
   call nullify_scf_history(scf_history(iimage))
   scf_history(iimage)%history_size=history_size
 end do

!PIMD only: fill in the data structure pimd_param
 if((dtset%imgmov==9).or.(dtset%imgmov==13))then
   pimd_param%irandom    =>dtset%irandom
   pimd_param%mdtemp     =>dtset%mdtemp
   pimd_param%vis        =>dtset%vis
   pimd_param%dtion      =>dtset%dtion
   pimd_param%optcell    =>dtset%optcell
   pimd_param%nnos       =>dtset%nnos
   pimd_param%ntypat     =>dtset%ntypat
   pimd_param%pimass     =>dtset%pimass
   pimd_param%pitransform=>dtset%pitransform
   pimd_param%bmass      =>dtset%bmass
   pimd_param%friction   =>dtset%friction
   pimd_param%strtarget  =>dtset%strtarget
   pimd_param%amu        =>dtset%amu
   pimd_param%qmass      =>dtset%qmass
   pimd_param%typat      =>dtset%typat
 end if

!If outputs are redirected (for each image),
!define names of the files
 if ((dtset%nimage>1)) then
   ABI_ALLOCATE(filout,(dtset%nimage))
   ABI_ALLOCATE(fillog,(dtset%nimage))
   call int2char4(mpi_enreg%me,tag)
   do ii=1,dtset%nimage
     call appdig(ii,'',appen)
     filout(ii)=trim(filnam(2))//'_IMG'//trim(appen)
     if (mpi_enreg%me_one_img==0) then
       fillog(ii)=trim(filnam(5))//'_IMG'//trim(appen)//"_LOG"
     else
       fillog(ii)=trim(filnam(5))//'_IMG'//trim(appen)//"_LOG_"//tag
     end if
   end do
   if (mpi_enreg%me==0) then
     do ii=1,dtset%nimage
       inquire(file=filout(ii),iostat=ios,exist=ex)
       if (ios/=0) ex=.false.
       if (ex) then
         open(unit=tmp_unit,file=filout(ii),iostat=ios)
         if (ios==0) close(unit=tmp_unit,status='delete')
       end if
     end do
   end if
   call xbarrier_mpi(mpi_enreg%comm_img)
!  MT june 2011: are the following lines useful ?
!  do ii=1,dtset%nimage
!  inquire(file=fillog(ii),iostat=ios,exist=ex)
!  if (ios/=0) ex=.false.
!  if (ex) then
!  open(unit=tmp_unit,file=fillog(ii),iostat=ios)
!  if (ios==0) close(unit=tmp_unit,status='delete')
!  end if
!  end do
 end if

 last_itimimage=1

 call timab(703,2,tsec)

!-----------------------------------------------------------------------------------------
!Big loop on the propagation of all images
 do itimimage=1,dtset%ntimimage


   call timab(704,1,tsec)

   itimimage_eff=itimimage;if(slideimg==1) itimimage_eff=1
   compute_all_images=(compute_static_images.and.itimimage==1)

!  Print title for time step
   if (dtset%nimage>1.or.dtset%ntimimage>1) then
     if (dtset%prtvolimg<2) then
       msg=ch10;if (itimimage >1) write(msg,'(2a)') ch10,ch10
       write(msg,'(5a)') trim(msg),&
&       '================================================================================',&
&       ch10,' ',trim(imagealgo_str(dtset%imgmov))
     else
       msg='';if (itimimage >1) msg=ch10
       write(msg,'(5a)') trim(msg),&
&       '--------------------------------------------------------------------------------',&
&       ch10,' ',trim(imagealgo_str(dtset%imgmov))
     end if
     if (dtset%ntimimage==1) write(msg,'(2a)')    trim(msg),' FOR 1 TIME STEP'
     if (dtset%ntimimage >1) write(msg,'(2a,i5)') trim(msg),' - TIME STEP ',itimimage
     if (dtset%prtvolimg<2) then
       write(msg,'(3a)') trim(msg),ch10,&
&       '================================================================================'
     end if
     call wrtout(ab_out ,msg,'COLL')
     call wrtout(std_out,msg,'PERS')
   end if

   call timab(704,2,tsec)

!  Loop on the dynamical images
   idynimage=0
   do iimage=1,nimage

     call timab(705,1,tsec)

     ii=mpi_enreg%index_img(iimage)
     if (dtset%dynimage(ii)==1) idynimage=idynimage+1

!    Compute static image only at first time step
     if (dtset%dynimage(ii)==1.or.compute_all_images) then

!      Change file names according to image index (if nimage>1)
       if (dtset%nimage>1) then
         call dtfil_init1(dtfil,dtset,filnam,filstat,idtset,jdtset,mpi_enreg,ndtset,&
&         image_index=ii)
         if (itimimage>1) then
           dtfil%ireadwf=0;dtfil%ireadden=0;dtfil%ireadkden=0
         end if
       end if

!      Redefine output units
       if (dtset%nimage>1) then
         if (mpi_enreg%paral_img==1) then
           call abi_io_redirect(new_io_comm=mpi_enreg%comm_one_img,&
&           new_leave_comm=mpi_enreg%comm_one_img)
         end if
         if (dtset%prtvolimg>0) then
           if (do_write_log) then
             call abi_io_redirect(new_ab_out=100)
             if (mpi_enreg%me_one_img==0) open(unit=ab_out,file=NULL_FILE,status='unknown')
           else
             call abi_io_redirect(new_ab_out=std_out)
           end if
         else if (mpi_enreg%paral_img==1) then
           call abi_io_redirect(new_ab_out=100)
           if (mpi_enreg%me_one_img==0) open(unit=ab_out,file=filout(ii),status='unknown')
         end if
         if (mpi_enreg%paral_img==1.and.mpi_enreg%me_one_img==0.and.do_write_log) then
           call abi_io_redirect(new_std_out=200)
           open(unit=std_out,file=fillog(ii),status='unknown')
         end if
       end if

!      Print title for image
       if (dtset%nimage>1.and.(dtset%prtvolimg==0.or.do_write_log)) then
         if (ii==1) write(msg,'(a)' ) ch10
         if (ii >1) write(msg,'(2a)') ch10,ch10
         write(msg,'(6a,i4,a,i4,3a)') trim(msg),&
&         '--------------------------------------------------------------------------------',ch10,&
&         ' ',trim(imagealgo_str(dtset%imgmov)),' - CELL # ',ii,'/',dtset%nimage,ch10,&
&         '--------------------------------------------------------------------------------',ch10
         if (dtset%prtvolimg==0) then
           call wrtout(ab_out ,msg,'COLL')
         end if
         if (do_write_log) then
           call wrtout(std_out,msg,'PERS')
         end if
       end if

       acell(:)  =results_img_timimage(iimage,itimimage_eff)%acell(:)
       rprim(:,:)=results_img_timimage(iimage,itimimage_eff)%rprim(:,:)
       vel(:,:)  =results_img_timimage(iimage,itimimage_eff)%vel(:,:)
       xred(:,:) =results_img_timimage(iimage,itimimage_eff)%xred(:,:)
       occ(:)    =occ_img(:,iimage)

       call timab(705,2,tsec)

       call gstate(acell,codvsn,cpui,dtfil,dtset,iexit,scf_initialized(iimage),&
&       mpi_enreg,npwtot,occ,pawang,pawrad,pawtab,psps,&
&       results_img_timimage(iimage,itimimage_eff)%results_gs,&
&       rprim,scf_history(iimage),vel,xred)

       call timab(706,1,tsec)

       if (dtset%dynimage(ii)==1) then
         results_img_timimage(iimage,itimimage_eff)%acell(:)  =acell(:)
         results_img_timimage(iimage,itimimage_eff)%rprim(:,:)=rprim(:,:)
         results_img_timimage(iimage,itimimage_eff)%vel(:,:)  =vel(:,:)
         results_img_timimage(iimage,itimimage_eff)%xred(:,:) =xred(:,:)
         occ_img(:,iimage)=occ(:)
       end if

!      Close output units ; restore defaults
       if (dtset%nimage>1) then
         if (dtset%prtvolimg>0) then
           if (do_write_log.and.mpi_enreg%me_one_img==0) close(unit=ab_out)
         else if (mpi_enreg%paral_img==1) then
           if (mpi_enreg%me_one_img==0) close(unit=ab_out)
         end if
         if (mpi_enreg%paral_img==1.and.mpi_enreg%me_one_img==0.and.do_write_log) close(unit=std_out)
         call abi_io_redirect(new_ab_out=ab_out_default,new_std_out=std_out_default,&
&         new_io_comm=xmpi_world,new_leave_comm=xmpi_world)
       end if

       call timab(706,2,tsec)

     end if
   end do ! iimage

   if(mpi_enreg%paral_img==1)then
     call timab(702,1,tsec)
     call leave_test()
     call timab(702,2,tsec)
   end if

   call timab(707,1,tsec)

!  Output when images are used
   if (dtset%nimage>1) then
!    === 1st option: reduced outputs ===
     if (dtset%prtvolimg>0) then
       call prtimg(dtset%dynimage,imagealgo_str(dtset%imgmov),dtset%imgmov,ab_out,&
&       mpi_enreg,nimage,dtset%nimage,compute_all_images,dtset%prtvolimg,&
&       results_img_timimage(:,itimimage_eff))

!      === 2nd option: full gathered outputs ===
     else if (mpi_enreg%paral_img==1) then
       if (mpi_enreg%me==0) then
         do ii=1,dtset%nimage
           open(unit=tmp_unit,file=filout(ii),status='old',iostat=ios)
           if (ios==0) then
             if (dtset%dynimage(ii)==1.or.compute_all_images) then
               eof=.false.
               do while (.not.eof)
                 read(tmp_unit,fmt='(a)',err=111,end=111,iostat=ios) line
                 if (ios==0) write(ab_out,'(a)') trim(line)
                 goto 112
                 111              eof=.true.
                 112              continue
               end do
             end if
             close(unit=tmp_unit,status='delete')
           end if
         end do
         call flush_unit(ab_out)
       end if
       call xbarrier_mpi(mpi_enreg%comm_img)
     end if
   end if

!  Manage log files when images are used
   if (dtset%nimage>1) then
     if (mpi_enreg%paral_img==1.and.do_write_log) then
       if (mpi_enreg%me==0) then
         do ii=1,dtset%nimage
           open(unit=tmp_unit,file=fillog(ii),status='old',iostat=ios)
           if (ios==0) then
             eof=.false.
             do while (.not.eof)
               read(tmp_unit,fmt='(a)',err=113,end=113,iostat=ios) line
               if (ios==0) write(std_out,'(a)') trim(line)
               goto 114
               113            eof=.true.
               114            continue
             end do
             close(unit=tmp_unit,status='delete')
           end if
         end do
!        else if (mpi_enreg%me_one_img==0) then
!        do iimage=1,nimage
!        ii=mpi_enreg%index_img(iimage)
!        open(unit=tmp_unit,file=fillog(ii),status='old',iostat=ios)
!        if (ios==0) then
!        eof=.false.
!        do while (.not.eof)
!        read(tmp_unit,fmt='(a)',err=115,end=115,iostat=ios) line
!        if (ios==0) write(std_out,'(a)') trim(line)
!        goto 116
!        115            eof=.true.
!        116            continue
!        end do
!        close(unit=tmp_unit,status='delete')
!        end if
!        end do
       end if
     end if
     call xbarrier_mpi(mpi_enreg%comm_img)
   end if

!  TESTS WHETHER ONE CONTINUES THE LOOP
!  Here we calculate the change in energy, and exit if delta_energy < tolimg
   delta_energy=zero
   if (slideimg/=1 .and. itimimage>1) then
     do idynimage=1,ndynimage
       iimage=list_dynimage(idynimage)
       delta_energy=delta_energy &
&       +abs(results_img_timimage(iimage,itimimage  )%results_gs%etotal &
&       -results_img_timimage(iimage,itimimage-1)%results_gs%etotal)
     end do
     if (mpi_enreg%paral_img==1) then
       call xsum_mpi(delta_energy,mpi_enreg%comm_img,ierr)
     end if
     delta_energy=delta_energy/dtset%ndynimage
     if (delta_energy<dtset%tolimg) then
       if (dtset%prtvolimg<2) then
         write(msg,'(5a,i5,6a,es11.3,a,es11.3,2a)') ch10,ch10,&
&         '================================================================================',ch10,&
&         ' At time step ',itimimage,ch10,&
&         ' ',trim(imagealgo_str(dtset%imgmov)),' has reached energy convergence',ch10,&
&         ' with Average[Abs(Etotal(t)-Etotal(t-dt))]=',delta_energy,'<tolimg=',dtset%tolimg,ch10,&
&         '================================================================================'
       else
         write(msg,'(4a,i5,6a,es11.3,a,es11.3)') ch10,&
&         '--------------------------------------------------------------------------------',ch10,&
&         ' At time step ',itimimage,ch10,&
&         ' ',trim(imagealgo_str(dtset%imgmov)),' has reached energy convergence',ch10,&
&         ' with Average[Abs(Etotal(t)-Etotal(t-dt))]=',delta_energy,'<tolimg=',dtset%tolimg
       end if
       call wrtout(ab_out ,msg,'COLL')
       call wrtout(std_out,msg,'COLL')
       call timab(707,2,tsec)
       exit   ! exit itimimage
     end if
   end if

!  In any case, stop at the maximal value of itimimage
   last_itimimage=itimimage_eff
   if(itimimage==dtset%ntimimage) then
     if (dtset%imgmov/=9.and.dtset%imgmov/=13) exit
   end if

!  Predict the next value of the images
   if (dtset%nimage>1) then
     call predictimg(delta_energy,dtset%fxcartfactor,dtset%iatfix,imagealgo_str(dtset%imgmov),&
&     dtset%imgmov,itimimage,list_dynimage,mpi_enreg,dtset%natom,&
&     ndynimage,nimage,dtset%nimage,dtset%ntimimage,pimd_param,&
&     dtset%prtvolimg,results_img_timimage)
   end if

   call timab(707,2,tsec)

 end do ! itimimage
!-----------------------------------------------------------------------------------------

 call timab(708,1,tsec)

!Copy the results of the computation in the appropriate arguments of the routine
 do iimage=1,nimage
   ii=mpi_enreg%index_img(iimage)
   if (dtset%dynimage(ii)==1) then
     acell_img(:,iimage)   =results_img_timimage(iimage,last_itimimage)%acell(:)
     rprim_img(:,:,iimage) =results_img_timimage(iimage,last_itimimage)%rprim(:,:)
     vel_img(:,:,iimage)   =results_img_timimage(iimage,last_itimimage)%vel(:,:)
     xred_img(:,:,iimage)  =results_img_timimage(iimage,last_itimimage)%xred(:,:)
     etotal_img(iimage)    =results_img_timimage(iimage,last_itimimage)%results_gs%etotal
     fcart_img(:,:,iimage) =results_img_timimage(iimage,last_itimimage)%results_gs%fcart(:,:)
     fred_img(:,:,iimage)  =results_img_timimage(iimage,last_itimimage)%results_gs%fred(:,:)
     strten_img(:,iimage)  =results_img_timimage(iimage,last_itimimage)%results_gs%strten(:)
   else if (compute_static_images) then
     etotal_img(iimage)    =results_img_timimage(iimage,1)%results_gs%etotal
     fcart_img(:,:,iimage) =results_img_timimage(iimage,1)%results_gs%fcart(:,:)
     fred_img(:,:,iimage)  =results_img_timimage(iimage,1)%results_gs%fred(:,:)
     strten_img(:,iimage)  =results_img_timimage(iimage,1)%results_gs%strten(:)
   end if
 end do

!When parallelizattion over images is activated, has to sum number of warnings
!and comments written in log file
 if (mpi_enreg%paral_img==1.and.mpi_enreg%comm_one_img==0) then
   call wrtout_myproc(std_out,msg,mpi_comm=mpi_enreg%comm_img)
 end if

!Final deallocations

 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(vel)
 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(list_dynimage)
 if ((dtset%nimage>1))  then
   ABI_DEALLOCATE(filout)
   ABI_DEALLOCATE(fillog)
 end if
 do itimimage=1,ntimimage_eff
   call destroy_results_img(results_img_timimage(:,itimimage))
 end do
 ABI_DEALLOCATE(results_img_timimage)
 do iimage=1,nimage
   call destroy_scf_history(scf_history(iimage))
 end do
 ABI_DEALLOCATE(scf_history)
 ABI_DEALLOCATE(scf_initialized)

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(708,2,tsec)
 call timab(700,2,tsec)

 DBG_EXIT("COLL")

end subroutine gstateimg
!!***
