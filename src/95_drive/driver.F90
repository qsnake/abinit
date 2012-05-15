!{\src2tex{textfont=tt}}
!!****f* ABINIT/driver
!! NAME
!! driver
!!
!! FUNCTION
!! Driver for ground state, response function, susceptibility, screening
!! and sigma calculations. The present routine drives the following operations.
!! An outer loop allows computation related to different data sets.
!! For each data set, either a GS calculation, a RF calculation,
!! a SUS calculation, a SCR calculation or a SIGMA calculation is made.
!! In both cases, the input variables are transferred in the proper variables,
!! selected big arrays are allocated, then the gstate, respfn or suscep subroutines are called.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG,MKV,MM,MT,FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! codvsn= code version
!! cpui=initial CPU time
!! filnam(5)=character strings giving file names
!! filstat=character strings giving name of status file
!! mpi_enreg=informations about MPI parallelization
!! ndtset=number of datasets
!! ndtset_alloc=number of datasets, corrected for allocation of at
!!               least one data set.
!! npsp=number of pseudopotentials
!! pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  Input/Output
!! results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!   Default values are set up in the calling routine
!! dtsets(0:ndtset_alloc)=<type datasets_type>
!!   intput: all input variables initialized from the input file.
!!   output: the effective set of variables used in the different datasets.
!!           Some variables, indeed, might have been redefined in one of the children.
!!           of this routine.
!!
!! NOTES
!! The array filnam is used for the name of input and output files,
!! and roots for generic input, output or temporary files.
!! Pseudopotential file names are set in pspini and pspatm,
!! using another name. The name filstat will be needed beyond gstate to check
!! the appearance of the "exit" flag, to make a hasty exit, as well as
!! in order to output the status of the computation.
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      bethe_salpeter,chkdilatmx,chkexi,clnmpi_img,destroy_results_out
!!      destroy_results_respfn,destroy_timer,dtfil_init1,dtfil_init_img
!!      dtsetcopy,dtsetfree,echo_xc_name,fftw3_cleanup,fftw3_set_nthreads
!!      find_getdtset,gather_results_out,gstateimg,gw_driver
!!      init_results_respfn,init_timer,initmpi_img,int2char4
!!      libxc_functionals_end,libxc_functionals_init,memocc,mkrdim,nonlinear
!!      pawalloc,print_timer,psps_free,psps_init_from_dtset,psps_init_global
!!      rdm,respfn,screening,sigma,status,suscep,timab,wrtout,wvl_timing
!!      xc_init,xcast_mpi,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine driver(codvsn,cpui,dtsets,filnam,filstat,&
&                 mpi_enreg,ndtset,ndtset_alloc,npsp,pspheads,results_out)

 use defs_basis
 use defs_parameters
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_results_out
 use m_results_respfn
 use m_xmpi
 use m_profiling
#ifdef HAVE_TIMER_ABINIT
 use m_timer
#endif
#if defined HAVE_DFT_LIBXC
 use libxc_functionals
#endif
#if defined HAVE_DFT_BIGDFT
 use BigDFT_API, only: xc_init, XC_MIXED, XC_ABINIT, wvl_timing => timing
#endif

 use m_fftw3,   only : fftw3_gain_wisdom, fftw3_cleanup, fftw3_set_nthreads

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'driver'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_27_toolbox_oop
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_abiutil
 use interfaces_56_xc
 use interfaces_59_io_mpi
 use interfaces_65_psp
 use interfaces_66_paw
 use interfaces_77_suscep
 use interfaces_93_rdm
 use interfaces_95_drive, except_this_one => driver
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndtset,ndtset_alloc,npsp
 real(dp),intent(in) :: cpui
 character(len=6),intent(in) :: codvsn
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 character(len=fnlen),intent(in) :: filnam(5)
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)
 type(results_out_type),target,intent(inout) :: results_out(0:ndtset_alloc)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=2
 integer,save :: paw_size_old=-1
 integer :: idtset,ierr,iexit,iget_cell,iget_occ,iget_vel,iget_xcart,iget_xred
 integer :: ii,iimage,iimage_get,ireadwf,jdtset
 integer :: jdtset_status!,mband
 integer :: mdtset,mtypalch,mu,mxnimage !,mdtsetmpsang,mgfft,mk1mem,mkmem,mkqmem,mpssoang,mpw
 integer :: nimage,openexit,paw_size,prtvol !,natom,nfft,nkpt,nspden,nsppol,nsym
 real(dp) :: etotal
 character(len=500) :: message
#ifdef HAVE_TIMER_ABINIT
 character(len=4) :: tag
 character(len=fnlen) :: timer_fname
#endif
 logical :: converged,results_gathered,test_img,use_results_all
 type(dataset_type) :: dtset
 type(datafiles_type) :: dtfil
 type(pawang_type) :: pawang
 type(pseudopotential_type) :: psps
 type(results_respfn_type) :: results_respfn
!!! type(vardims_type) :: abidims
!arrays
 integer :: mkmems(3)
 integer,allocatable :: jdtset_(:),npwtot(:)
 real(dp) :: acell(3),rprim(3,3),rprimd(3,3),rprimdget(3,3),strten(6),tsec(2)
 real(dp),allocatable :: acell_img(:,:),rprim_img(:,:,:)
 real(dp),allocatable :: fcart(:,:),fred(:,:)
 real(dp),allocatable :: fcart_img(:,:,:),fred_img(:,:,:)
 real(dp),allocatable :: etotal_img(:),strten_img(:,:)
 real(dp),allocatable :: miximage(:,:)
 real(dp),allocatable :: occ(:),vel(:,:),xcart(:,:),xred(:,:),xredget(:,:)
 real(dp),allocatable :: occ_img(:,:),vel_img(:,:,:),xred_img(:,:,:)
 type(pawrad_type),allocatable :: pawrad(:)
 type(pawtab_type),allocatable :: pawtab(:)
 type(results_out_type),pointer :: results_out_all(:)

!******************************************************************

 DBG_ENTER("COLL")

 ireadwf=0

 call timab(640,1,tsec)
 call timab(641,1,tsec)
 call status(0,filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 prtvol=dtsets(1)%prtvol
 if(prtvol==-level)then
   write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&   ' driver : enter , debug mode '
   call wrtout(std_out,message,'COLL')
 end if

 mdtset=9999

 if(ndtset>mdtset)then
   write(message,'(a,i2,a,i5,a)')&
&   '  The maximal allowed ndtset is ',mdtset,' while the input value is ',ndtset,'.'
   MSG_BUG(message)
 end if

 mtypalch=dtsets(1)%ntypalch
 do ii=1,ndtset_alloc
   mtypalch=max(dtsets(ii)%ntypalch,mtypalch)
 end do
 call psps_init_global(mtypalch, npsp, psps, pspheads)

 ABI_ALLOCATE(jdtset_,(0:ndtset))
 if(ndtset/=0)then
   jdtset_(:)=dtsets(0:ndtset)%jdtset
 else
   jdtset_(0)=0
 end if

 do idtset=1,ndtset_alloc
   call init_results_respfn(dtsets,ndtset_alloc,results_respfn)
 end do

 call timab(641,2,tsec)

!*********************************************************************
!Big loop on datasets

!Do loop on idtset (allocate statements are present)
 do idtset=1,ndtset_alloc

   call timab(642,1,tsec)

#ifdef HAVE_TIMER_ABINIT
#if 0
   call init_timer()
#endif
#endif

   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=1

   if(ndtset>=2)then
     jdtset_status=jdtset
   else
     jdtset_status=0
   end if

   call status(jdtset_status,filstat,iexit,level,'loop jdtset   ')

   if(jdtset>=100)then
     write(message,'(a,80a,a,a,i4,a,64a,a)') ch10,&
&     ('=',mu=1,80),ch10,&
&     '== DATASET ',jdtset,' ',('=',mu=1,64),ch10
   else
     write(message,'(a,80a,a,a,i2,a,66a,a)') ch10,&
&     ('=',mu=1,80),ch10,&
&     '== DATASET ',jdtset,' ',('=',mu=1,66),ch10
   end if
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'PERS')     ! PERS is choosen to make debugging easier

!  Copy input variables into a local dtset.
   call dtsetCopy(dtset, dtsets(idtset))

!  Set other values
   dtset%jdtset = jdtset
   dtset%ndtset = ndtset

!  Copy input values
   mkmems(1)        = dtset%mkmem
   mkmems(2)        = dtset%mkqmem
   mkmems(3)        = dtset%mk1mem

   mpi_enreg%paral_compil_respfn=dtset%paral_rf
   mpi_enreg%ngroup_respfn=dtset%ngroup_rf

   ABI_ALLOCATE(occ,(dtset%mband*dtset%nkpt*dtset%nsppol))
   ABI_ALLOCATE(vel,(3,dtset%natom))
   ABI_ALLOCATE(xred,(3,dtset%natom))

!  Initialize MPI data for the parallelism over images
   call initmpi_img(dtset,mpi_enreg,-1)
   nimage=mpi_enreg%nimage

!  Retrieve arrays concerned by the parallelism over images
   ABI_ALLOCATE(acell_img,(3,nimage))
   ABI_ALLOCATE(occ_img,(dtset%mband*dtset%nkpt*dtset%nsppol,nimage))
   ABI_ALLOCATE(rprim_img,(3,3,nimage))
   ABI_ALLOCATE(vel_img,(3,dtset%natom,nimage))
   ABI_ALLOCATE(xred_img,(3,dtset%natom,nimage))
   do iimage=1,nimage
     ii=mpi_enreg%index_img(iimage)
     acell_img(:  ,iimage) = dtset%acell_orig(:,ii)
     rprim_img(:,:,iimage) = dtset%rprim_orig(:,:,ii)
     vel_img  (:,:,iimage) = dtset%vel_orig(:,1:dtset%natom,ii)
     xred_img (:,:,iimage) = dtset%xred_orig(:,1:dtset%natom,ii)
!    Note that occ is not supposed to depend on the image in the input file.
!    However, the results will depend on the image. And occ can be used
!    to reinitialize the next dataset, or the next timimage.
     occ_img  (:  ,iimage) = dtset%occ_orig(1:dtset%mband*dtset%nkpt*dtset%nsppol)
   end do

!  ****************************************************************************
!  Treat the file names (get variables)

   call status(jdtset_status,filstat,iexit,level,'filenames     ')

   if (dtset%optdriver==RUNL_GSTATE.and.dtset%nimage>1) then
     call dtfil_init_img(dtfil,dtset,dtsets,idtset,jdtset_,ndtset,ndtset_alloc)
!    Call to dtfil_init1 is postponed in the loop over images
   else
     call dtfil_init1(dtfil,dtset,filnam,filstat,idtset,jdtset_,mpi_enreg,ndtset)
   end if

!  ****************************************************************************
!  Treat other get variables

!  If multi dataset mode, and already the second dataset,
!  treatment of other get variables.
   if( ndtset>1 .and. idtset>1 )then

     call status(jdtset_status,filstat,iexit,level,'get variables ')

!    Check if parallelization over images is activated
     mxnimage=maxval(dtsets(1:ndtset_alloc)%nimage)
     ABI_ALLOCATE(miximage,(mxnimage,mxnimage))
     test_img=(mxnimage/=1.and.maxval(dtsets(:)%npimage)>1)
     use_results_all=.false.
     if (test_img.and.mpi_enreg%me_one_img==0) then
       use_results_all=.true.
       ABI_ALLOCATE(results_out_all,(0:ndtset_alloc))
     else
       results_out_all => results_out
     end if
     results_gathered=.false.


     call find_getdtset(dtsets,dtset%getocc,'getocc',idtset,iget_occ,miximage,mxnimage,ndtset_alloc)
     if(iget_occ/=0)then
!      Gather contributions to results_out from images, if needed
       if (test_img.and.mpi_enreg%me_one_img==0.and.(.not.results_gathered)) then
         call gather_results_out(dtsets,results_out,results_out_all,use_results_all, &
&         allgather=.true.,only_one_per_img=.true.)
         results_gathered=.true.
       end if
       if ((.not.test_img).or.mpi_enreg%me_one_img==0) then
         do iimage=1,nimage
           ii=mpi_enreg%index_img(iimage)
           occ_img(:,iimage)=zero
           do iimage_get=1,dtsets(iget_occ)%nimage
             occ_img(:,iimage)=occ_img(:,iimage)+miximage(ii,iimage_get) &
&             *results_out_all(iget_occ)%occ(1:dtset%mband*dtset%nkpt*dtset%nsppol,iimage_get)
           end do
         end do
       end if
     end if

!    Getcell has to be treated BEFORE getxcart
!    since acell and rprim will be used
     call find_getdtset(dtsets,dtset%getcell,'getcell',idtset,iget_cell,miximage,mxnimage,ndtset_alloc)
     if(iget_cell/=0)then
!      Gather contributions to results_out from images, if needed
       if (test_img.and.mpi_enreg%me_one_img==0.and.(.not.results_gathered)) then
         call gather_results_out(dtsets,results_out,results_out_all,use_results_all, &
&         allgather=.true.,only_one_per_img=.true.)
         results_gathered=.true.
       end if
       if ((.not.test_img).or.mpi_enreg%me_one_img==0) then
         do iimage=1,nimage
           ii=mpi_enreg%index_img(iimage)
           acell_img(:,iimage)=zero
           rprim_img(:,:,iimage)=zero
           do iimage_get=1,dtsets(iget_cell)%nimage
             acell_img(:,iimage)=acell_img(:,iimage)+&
&             miximage(ii,iimage_get)*results_out_all(iget_cell)%acell(:,iimage_get)
             rprim_img(:,:,iimage)=rprim_img(:,:,iimage)+&
&             miximage(ii,iimage_get)*results_out_all(iget_cell)%rprim(:,:,iimage_get)
!            Check that the new acell and rprim are consistent with the input dilatmx
             call mkrdim(acell_img(:,iimage),rprim_img(:,:,iimage),rprimd)
             call chkdilatmx(dtset%dilatmx,rprimd,dtset%rprimd_orig(1:3,1:3,iimage))
           end do
         end do
       end if
     end if

     call find_getdtset(dtsets,dtset%getxred,'getxred',idtset,iget_xred,miximage,mxnimage,ndtset_alloc)
     if(iget_xred/=0)then
!      Gather contributions to results_out from images, if needed
       if (test_img.and.mpi_enreg%me_one_img==0.and.(.not.results_gathered)) then
         call gather_results_out(dtsets,results_out,results_out_all,use_results_all, &
&         allgather=.true.,only_one_per_img=.true.)
         results_gathered=.true.
       end if
       if ((.not.test_img).or.mpi_enreg%me_one_img==0) then
         do iimage=1,nimage
           ii=mpi_enreg%index_img(iimage)
           xred_img(:,:,iimage)=zero
           do iimage_get=1,dtsets(iget_xred)%nimage
             xred_img(:,:,iimage)=xred_img(:,:,iimage)+&
&             miximage(ii,iimage_get)*results_out_all(iget_xred)%xred(:,1:dtset%natom,iimage_get)
           end do
         end do
       end if
     end if

     call find_getdtset(dtsets,dtset%getxcart,'getxcart',idtset,iget_xcart,miximage,mxnimage,ndtset_alloc)
     if(iget_xcart/=0)then
!      Gather contributions to results_out from images, if needed
       if (test_img.and.mpi_enreg%me_one_img==0.and.(.not.results_gathered)) then
         call gather_results_out(dtsets,results_out,results_out_all,use_results_all, &
&         allgather=.true.,only_one_per_img=.true.)
         results_gathered=.true.
       end if
       if ((.not.test_img).or.mpi_enreg%me_one_img==0) then
         ABI_ALLOCATE(xcart,(3,dtset%natom))
         ABI_ALLOCATE(xredget,(3,dtset%natom))
         do iimage=1,nimage
           ii=mpi_enreg%index_img(iimage)
           xred_img(:,:,iimage)=zero
           do iimage_get=1,dtsets(iget_xcart)%nimage
!            Compute xcart of the previous dataset
             rprimdget(:,:)=results_out_all(iget_xcart)%rprimd(:,:,iimage_get)
             xredget (:,:) =results_out_all(iget_xcart)%xred(:,1:dtset%natom,iimage_get)
             call xredxcart(dtset%natom,1,rprimdget,xcart,xredget)
!            xcart from previous dataset is computed. Now, produce xred for the new dataset,
!            with the new acell and rprim ...
             call mkrdim(acell_img(:,iimage),rprim_img(:,:,iimage),rprimd)
             call xredxcart(dtset%natom,-1,rprimd,xcart,xredget(:,:))
             xred_img(:,:,iimage)=xred_img(:,:,iimage)+miximage(ii,iimage_get)*xredget(:,:)
           end do
         end do
         ABI_DEALLOCATE(xcart)
         ABI_DEALLOCATE(xredget)
       end if
     end if

     call find_getdtset(dtsets,dtset%getvel,'getvel',idtset,iget_vel,miximage,mxnimage,ndtset_alloc)
     if(iget_vel/=0)then
!      Gather contributions to results_out from images, if needed
       if (test_img.and.mpi_enreg%me_one_img==0.and.(.not.results_gathered)) then
         call gather_results_out(dtsets,results_out,results_out_all,use_results_all, &
&         allgather=.true.,only_one_per_img=.true.)
         results_gathered=.true.
       end if
       if ((.not.test_img).or.mpi_enreg%me_one_img==0) then
         do iimage=1,nimage
           ii=mpi_enreg%index_img(iimage)
           vel_img(:,:,iimage)=zero
           do iimage_get=1,dtsets(iget_vel)%nimage
             vel_img(:,:,iimage)=vel_img(:,:,iimage)+&
&             miximage(ii,iimage_get)*results_out_all(iget_vel)%vel(:,1:dtset%natom,iimage_get)
           end do
         end do
       end if
     end if

!    In the case of parallelization over images, has to distribute data
     if (test_img) then
       if (iget_occ/=0) then
         call xcast_mpi(occ_img,0,mpi_enreg%comm_one_img,ierr)
       end if
       if (iget_cell/=0) then
         call xcast_mpi(acell_img,0,mpi_enreg%comm_one_img,ierr)
         call xcast_mpi(rprim_img,0,mpi_enreg%comm_one_img,ierr)
       end if
       if (iget_vel/=0) then
         call xcast_mpi(vel_img,0,mpi_enreg%comm_one_img,ierr)
       end if
       if (iget_xred/=0.or.iget_xcart/=0) then
         call xcast_mpi(xred_img,0,mpi_enreg%comm_one_img,ierr)
       end if
     end if

!    Clean memory
     ABI_DEALLOCATE(miximage)
     if (test_img.and.mpi_enreg%me==0) then
       if (results_gathered) then
         call destroy_results_out(results_out_all)
       end if
       ABI_DEALLOCATE(results_out_all)
     else
       nullify(results_out_all)
     end if

   end if

!  ****************************************************************************
!  Treat the pseudopotentials : initialize the psps variable

   call status(jdtset_status,filstat,iexit,level,'init psps     ')
   call psps_init_from_dtset(dtset, idtset, psps, pspheads)

!  ****************************************************************************
!  PAW allocations.

   call status(jdtset_status,filstat,iexit,level,'PAW allocs    ')

!  The correct dimension of pawrad/tab is ntypat.
!  In case of paw, no alchemical psp is allowed, so npsp=ntypat
!  However, in case of alchemical psps, pawrad/tab(ipsp) is invoked in
!  pspini. So, in order to avoid any problem, declare pawrad/tab
!  at paw_size=max(ntypat,npsp).
   paw_size=0;if (psps%usepaw==1) paw_size=max(dtset%ntypat,psps%npsp)
   if (paw_size/=paw_size_old) then
     if(idtset/=1) then
       call pawalloc(dtset,idtset,psps%mpsang,psps%mqgrid_vl,npsp,2,paw_size,paw_size_old,&
&       pawang,pawrad,pawtab,pspheads)
       ABI_DEALLOCATE(pawrad)
       ABI_DEALLOCATE(pawtab)
     end if
     ABI_ALLOCATE(pawrad,(paw_size))
     ABI_ALLOCATE(pawtab,(paw_size))
   end if
   call pawalloc(dtset,idtset,psps%mpsang,psps%mqgrid_vl,npsp,1,paw_size,paw_size_old,&
&   pawang,pawrad,pawtab,pspheads)
   paw_size_old=paw_size

!  ****************************************************************************
!  WVL allocations.

   call status(jdtset_status,filstat,iexit,level,'WVL allocs    ')

!  Set up mpi informations from the dataset
   if (dtset%usewvl == 1) then
!    WVL - data distribution
     ABI_ALLOCATE(mpi_enreg%nscatterarr,(0:mpi_enreg%nproc - 1, 4))
     ABI_ALLOCATE(mpi_enreg%ngatherarr,(0:mpi_enreg%nproc - 1, 2))
#if defined HAVE_DFT_BIGDFT
!    WVL - debugging
     call wvl_timing(mpi_enreg%nproc,trim(dtfil%filnam_ds(4)) // "_TIME.prc",'IN')
#endif
   end if

!  ****************************************************************************
!  ETSF I/O

   call status(jdtset_status,filstat,iexit,level,'ETSF IO       ')

!  DC20101117: currently unused code.
!  !!  Fill in abidims structure
!  !!   abidims%mband    = dtset%mband
!  !!   abidims%mproj    = 1                 ! FIXME
!  !!   abidims%mpsang   = psps%mpsang
!  !!   abidims%mpw      = dtset%mpw
!  !!   abidims%natom    = dtset%natom
!  !!   abidims%natsph   = dtset%natsph
!  !!   abidims%nberry   = dtset%nberry
!  !!   abidims%nconeq   = dtset%nconeq
!  !!   abidims%nfft     = dtset%nfft
!  !!   abidims%nfreqsus = dtset%nfreqsus
!  !!   abidims%ngrid1   = dtset%ngfft(1)
!  !!   abidims%ngrid2   = dtset%ngfft(2)
!  !!   abidims%ngrid3   = dtset%ngfft(3)
!  !!   abidims%nkpt     = dtset%nkpt
!  !!   abidims%nkptgw   = dtset%nkptgw
!  !!   abidims%npsp     = npsp
!  !!   abidims%npspalch = dtset%npspalch
!  !!   abidims%nqptdm   = dtset%nqptdm
!  !!   abidims%nshiftk  = dtset%nshiftk
!  !!   abidims%nspden   = dtset%nspden
!  !!   abidims%nspinor  = dtset%nspinor
!  !!   abidims%nsppol   = dtset%nsppol
!  !!   abidims%nsym     = dtset%nsym
!  !!   abidims%ntypat   = dtset%ntypat
!  !!   abidims%ntypalch = dtset%ntypalch
!  !!   abidims%wfs_dim1 = -1
!  !!   abidims%wfs_dim2 = -1
!  !!   abidims%npw_tiny = -1


!  ****************************************************************************
!  At this stage, all the data needed for the treatment of one dataset
!  have been transferred from multi-dataset arrays.

   iexit=0

!  Smaller integer arrays :
   ABI_ALLOCATE(npwtot,(dtset%nkpt))
   ABI_ALLOCATE(fcart,(3,dtset%natom))
   ABI_ALLOCATE(fred,(3,dtset%natom))

!  Initialize these to zero (needed when hasty exit)
   etotal=zero
   strten(:)=zero
   fcart(:,:)=zero ; fred(:,:)=zero

   if(dtset%optdriver/=RUNL_GSTATE)then
     acell(:)=acell_img(:,1)
     occ(:)=occ_img(:,1)
     rprim(:,:)=rprim_img(:,:,1)
     vel(:,:)=vel_img(:,:,1)
     xred(:,:)=xred_img(:,:,1)
   end if

   call echo_xc_name (dtset%ixc)

#if defined HAVE_DFT_LIBXC
   if (dtset%ixc<0) then
     call libxc_functionals_init(dtset%ixc,dtset%nspden)
   end if
#endif
#if defined HAVE_DFT_BIGDFT
   if (dtset%usewvl > 0) then
     if (dtset%ixc < 0) then
       call xc_init(dtset%ixc, XC_MIXED, dtset%nsppol)
     else
       call xc_init(dtset%ixc, XC_ABINIT, dtset%nsppol)
     end if
   end if
#endif


!  #ifdef HAVE_FFTW3
   if (dtset%ngfft(7)==312) then
     call fftw3_set_nthreads()
!    $call fftw3_gain_wisdom(dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),1) ! ndat set to 1
   end if
!  #endif

   call timab(642,2,tsec)

   select case(dtset%optdriver)

     case(RUNL_GSTATE)

       call status(jdtset_status,filstat,iexit,level,'call gstateimg')

       ABI_ALLOCATE(fcart_img,(3,dtset%natom,nimage))
       ABI_ALLOCATE(fred_img,(3,dtset%natom,nimage))
       ABI_ALLOCATE(etotal_img,(nimage))
       ABI_ALLOCATE(strten_img,(6,nimage))

       call gstateimg(acell_img,codvsn,cpui,dtfil,dtset,etotal_img,fcart_img,&
&       fred_img,iexit,mpi_enreg,nimage,npwtot,occ_img,&
&       pawang,pawrad,pawtab,psps,rprim_img,strten_img,vel_img,xred_img,&
&       filnam,filstat,idtset,jdtset_,ndtset)

       call status(jdtset_status,filstat,iexit,level,'after gstateimg')

     case(RUNL_RESPFN)

       call status(jdtset_status,filstat,iexit,level,'call respfn   ')

       mpi_enreg%paral_level=2
       call respfn(codvsn,cpui,dtfil,dtset,etotal,iexit,mkmems,mpi_enreg,&
&       npwtot,occ,pawang,pawrad,pawtab,psps,results_respfn,xred)

       call status(jdtset_status,filstat,iexit,level,'after respfn  ')

     case(RUNL_SUSCEP)

       call status(jdtset_status,filstat,iexit,level,'call suscep   ')

       mpi_enreg%paral_level=2
       call suscep(dtfil,dtset,iexit,&
&       dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,dtset%nfft,dtset%nkpt,&
&       dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%nsym,occ,xred)

       call status(jdtset_status,filstat,iexit,level,'after suscep  ')

     case(RUNL_SCREENING)

       call screening(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim)

     case(RUNL_SIGMA)

       call sigma(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim,xred,converged)

     case (RUNL_SCGW)

       call gw_driver(idtset,jdtset_,ndtset,acell,codvsn,filnam,dtfil,dtset,&
&       pawang,pawrad,pawtab,psps,rprim,xred)

     case(RUNL_NONLINEAR)

       call status(jdtset_status,filstat,iexit,level,'call nonlinear   ')

       mpi_enreg%paral_level=2
       call nonlinear(codvsn,dtfil,dtset,etotal,iexit,&
&       dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,dtset%nfft,dtset%nkpt,npwtot,dtset%nspden,&
&       dtset%nspinor,dtset%nsppol,dtset%nsym,occ,pawrad,pawtab,psps,xred)

       call status(jdtset_status,filstat,iexit,level,'after nonlinear  ')

     case(6)

       write(message, '(a,i12,a,a,a,a)' )&
&       '  The optdriver value 6 has been disabled since ABINITv6.0.',ch10,&
&       '  Action : modify optdriver in the input file.'
       MSG_ERROR(message)

     case (RUNL_BSE) ! Bethe-Salpeter
       call status(jdtset_status,filstat,iexit,level,'call bethe_salpeter')

       call bethe_salpeter(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim,xred)

       call status(jdtset_status,filstat,iexit,level,'after bethe_salpeter')

     case(RUNL_RDM)

       call rdm(acell,dtfil,dtset,pawtab,mpi_enreg,rprim)

       case default ! Bad value for optdriver
       write(message,'(a,i12,4a)')&
&       '  Unknown value for the variable optdriver: ',dtset%optdriver,ch10,&
&       '  This is not allowed. ',ch10,&
&       '  Action : modify optdriver in the input file.'
       MSG_ERROR(message)
   end select

   call timab(643,1,tsec)
!  ****************************************************************************

!  Transfer of multi dataset outputs from temporaries :
!  acell, xred, occ rprim, and vel might be modified from their
!  input values
!  etotal, fcart, fred, and strten have been computed
!  npwtot was already computed before, but is stored only now

   if(dtset%optdriver==RUNL_GSTATE)then
     do iimage=1,nimage
       results_out(idtset)%etotal(iimage)                 =etotal_img(iimage)
!      results_out(idtset)%deltae(iimage)                 =etotal_img(iimage)
!      results_out(idtset)%diffor(iimage)                 =etotal_img(iimage)
!      results_out(idtset)%residm(iimage)                 =etotal_img(iimage)
!      results_out(idtset)%res2  (iimage)                 =etotal_img(iimage)
       results_out(idtset)%acell(:,iimage)                =acell_img(:,iimage)
       results_out(idtset)%rprim(:,:,iimage)              =rprim_img(:,:,iimage)
       call mkrdim(acell_img(:,iimage),rprim_img(:,:,iimage),rprimd)
       results_out(idtset)%rprimd(:,:,iimage)             =rprimd(:,:)
       results_out(idtset)%strten(:,iimage)                =strten_img(:,iimage)
       results_out(idtset)%fcart(1:3,1:dtset%natom,iimage)=fcart_img(:,:,iimage)
       results_out(idtset)%fred(1:3,1:dtset%natom,iimage) =fred_img(:,:,iimage)
       results_out(idtset)%npwtot(1:dtset%nkpt,iimage)    =npwtot(1:dtset%nkpt)
       results_out(idtset)%occ(1:dtset%mband*dtset%nkpt*dtset%nsppol,iimage)=&
&       occ_img(1:dtset%mband*dtset%nkpt*dtset%nsppol,iimage)
       results_out(idtset)%vel(:,1:dtset%natom,iimage)    =vel_img(:,:,iimage)
       results_out(idtset)%xred(:,1:dtset%natom,iimage)   =xred_img(:,:,iimage)
     end do
     ABI_DEALLOCATE(etotal_img)
     ABI_DEALLOCATE(fcart_img)
     ABI_DEALLOCATE(fred_img)
     ABI_DEALLOCATE(strten_img)
   else
     results_out(idtset)%acell(:,1)                =acell(:)
     results_out(idtset)%etotal(1)                 =etotal
     results_out(idtset)%rprim(:,:,1)              =rprim(:,:)
     call mkrdim(acell,rprim,rprimd)
     results_out(idtset)%rprimd(:,:,1)             =rprimd(:,:)
     results_out(idtset)%npwtot(1:dtset%nkpt,1)    =npwtot(1:dtset%nkpt)
     results_out(idtset)%occ(1:dtset%mband*dtset%nkpt*dtset%nsppol,1)=&
&     occ(1:dtset%mband*dtset%nkpt*dtset%nsppol)
     results_out(idtset)%xred(:,1:dtset%natom,1)   =xred(:,:)
   end if
   ABI_DEALLOCATE(fcart)
   ABI_DEALLOCATE(fred)
   ABI_DEALLOCATE(acell_img)
   ABI_DEALLOCATE(occ_img)
   ABI_DEALLOCATE(rprim_img)
   ABI_DEALLOCATE(vel_img)
   ABI_DEALLOCATE(xred_img)

!  Clean MPI data (images parallelization)
   call clnmpi_img(mpi_enreg)

#if defined HAVE_DFT_LIBXC
   if (dtset%ixc<0) then
     call libxc_functionals_end()
   end if
#endif

!  MG TODO Update dtset(idtset) since some values might have been changed.
!  Several automatic tests will change due to this change though.
#if 0
   call dtsetFree(dtsets(idtset))
   call dtsetCopy(dtsets(idtset),dtset)
#endif

   call dtsetFree(dtset)

   ABI_DEALLOCATE(occ)
   ABI_DEALLOCATE(vel)
   ABI_DEALLOCATE(xred)
   ABI_DEALLOCATE(npwtot)

#ifdef HAVE_TIMER_ABINIT
   call int2char4(mpi_enreg%me,tag)
   timer_fname=TRIM(Dtfil%filnam_ds(4))//"_P-"//TRIM(tag)//'_TIMER'
#if 0 
   call print_timer(timer_fname)
   call destroy_timer()
#endif
#endif

!  Check whether exiting was required by the user.
!  If found then beat a hasty exit from time steps
   openexit=1 ; if(dtset%chkexit==0) openexit=0

   call chkexi(zero,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)

   call timab(643,2,tsec)

   if (iexit/=0)exit

 end do ! idtset (allocate statements are present - an exit statement is present)

!*********************************************************************

 call timab(644,1,tsec)

!PSP deallocation
 call psps_free(psps)

!PAW deallocation
 call pawalloc(dtset,idtset,psps%mpsang,psps%mqgrid_vl,npsp,3,paw_size,paw_size_old,&
& pawang,pawrad,pawtab,pspheads)
 ABI_DEALLOCATE(pawrad)
 ABI_DEALLOCATE(pawtab)

 ABI_DEALLOCATE(jdtset_)

!Results_respfn deallocation
 call destroy_results_respfn(results_respfn)

!#ifdef HAVE_FFTW3
 if (dtset%ngfft(7)==312) then
   call fftw3_cleanup()
 end if
!#endif

#if defined HAVE_DFT_BIGDFT
 if (dtset%usewvl == 1) then
!  WVL - debugging
   call memocc(0,mpi_enreg%me,'count','stop')
   call wvl_timing(mpi_enreg%me,'WFN_OPT','PR')
   call wvl_timing(mpi_enreg%me,'              ','RE')
 end if
#endif

 call status(0,filstat,iexit,level,'exit          ')

 call timab(644,2,tsec)
 call timab(640,2,tsec)

 DBG_EXIT("COLL")

end subroutine driver
!!***
