!{\src2tex{textfont=tt}}
!!****p* ABINIT/mrgscr
!! NAME
!! mrgscr
!!
!! FUNCTION
!! This code reads partial (SCR|SUSC) files for different q points creating a single file that
!! can be used to perform a sigma calculation.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (RS, MG, MS)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! If the number of SCR files to be merged is equal to 1, the program checks
!! the integrity of the file reporting the list of missing q-points.
!! Note that the list of required q-points depends on the k-mesh
!! used during the calculation of the KSS file. We assume indeed that the same k-mesh
!! is used during the calculation of the matrix elements of sigma.
!!
!! INPUTS
!!  (Main program)
!!
!! OUTPUT
!!  Only checking and writing
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,copy_scrhdr,cqratio,cspint,decompose_epsm1
!!      destroy_bz_mesh_type,destroy_crystal,destroy_epsilonm1_results
!!      destroy_gpairs_type,destroy_gsphere,destroy_mpi_enreg,find_qmesh,fourdp
!!      free_scrhdr,get_ppm_eigenvalues,get_rhor,getem1_from_ppm_one_ggp,getng
!!      herald,im_screening,init_crystal_from_hdr,init_er_from_file
!!      init_gpairs_type,init_gsphere,init_kmesh,initmpi_seq,int2char4
!!      interpolate_w,isfile,merge_scrhdr,metric,mkdump_er,mpi_comm_rank
!!      mpi_comm_size,my_gwannier,nullify_epsilonm1_results,nullify_gpairs_type
!!      nullify_mpi_enreg,ppm_free,ppm_init,print_bz_mesh
!!      print_epsilonm1_results,print_scrhdr,prompt
!!      re_and_im_screening_with_phase,re_screening,read_pole_screening
!!      read_screening,remove_phase,scr_hdr_io,sequential_fitting,setup_ppmodel
!!      sort_dp,test_charge,vcoul_free,vcoul_init,write_pole_screening
!!      write_screening,wrtout,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program mrgscr

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_build_info
 use m_errors
 use m_gwdefs,         only : GW_TOLQ, GW_TOLQ0
 use m_io_tools,       only : prompt, get_unit, file_exist
 use m_numeric_tools,  only : iseven
 use m_header,         only : hdr_get_nelect_byocc
 use m_geometry,       only : normv
 use m_crystal,        only : destroy_crystal, crystal_structure
 use m_crystal_io,     only : init_crystal_from_hdr
 use m_gsphere,        only : init_gsphere, destroy_gsphere, nullify_Gpairs_type, destroy_Gpairs_type,&
&                             init_Gpairs_type, gvectors_type, gpairs_type
 use m_bz_mesh,        only : bz_mesh_type, find_qmesh, init_kmesh, print_BZ_mesh, destroy_bz_mesh_type
 use m_vcoul,          only : vcoul_t, vcoul_init, vcoul_free
 use m_io_screening,   only : scr_hdr_io, print_scrhdr, copy_scrhdr, merge_scrhdr, read_screening, &
&                             write_screening, write_pole_screening,free_scrhdr, ScrHdr_type, &
&                             read_pole_screening
 use m_ppmodel,        only : ppm_init, ppm_free, setup_ppmodel, getem1_from_PPm_one_ggp, &
&                             get_PPm_eigenvalues, ppmodel_type, cqratio
 use m_model_screening,     only : init_peaks_from_grid,im_screening,re_screening,remove_phase, &
&                                  init_peaks_even_dist,init_single_peak,sequential_fitting, &
&                                  re_and_im_screening_with_phase
 use m_levenberg_marquardt, only : lmdif1,lm_fit_print_info
 use m_gwannier,       only : my_gwannier
 use m_screening,      only : mkdump_er, destroy_epsilonm1_results, print_epsilonm1_results,decompose_epsm1,&
&                             nullify_epsilonm1_results, init_er_from_file, Epsilonm1_results, interpolate_w
 use m_gwdefs,         only : GW_Q0_DEFAULT
 use m_fft_mesh,       only : g2ifft

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mrgscr'
 use interfaces_14_hidewrite
 use interfaces_27_toolbox_oop
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_69_wfdesc
 use interfaces_70_gw
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,paral_kgb0=0
 integer :: accesswff,fform1,fform_merge
 integer :: ifile,ifound,ierr,ii,ios,iqiA,iqibz,iqf,istat !,idx1,idx2
 integer :: nfiles,nomega4m,npwe4m,rdwr,timrev,prtvol=0,unt_out,unitem1
 integer :: unt_dump,idx,ig1,ig2,iomega,ppmodel,npwe_asked,mqmem,io
 integer :: unt_dump2
 integer :: id_required,ikxc,approx_type,option_test,dim_kxcg
 integer :: usexcnhat,usefinegrid
 integer :: mgfft,nqlwl,nfft,igmax,spaceComm
 integer :: nq_selected,nqibzA,nomega_asked,unt_tmp,kptopt
 integer :: choice,indx_imfreq_file,nfreq_tot,nfreqre,nfreqim,nfreqc!,npairs
 integer :: ifrq,npwe4mI,npwe4mJ,indx_q0,imax,new_nkibz,old_nshiftq
 integer :: ig1_start,ig1_end,ig2_start,ig2_end,gmgp_idx
 integer :: orig_npwe,new_nfreqre
#ifdef HAVE_ALGO_LEVMAR
 integer :: npoles
#endif
 real(dp) :: ucvol,boxcutmin,ecut,drude_plsmf,compch_fft,compch_sph
 real(dp) :: nelectron_exp,freqremax,freqremin,eps_diff,eps_norm,eps_ppm_norm
 real(dp) :: value1,value2,factor,GN_drude_plsmf,domegareal,freqimmax
 real(dp) :: dwre,dwim
 real(gwp) :: phase
#ifdef HAVE_ALGO_LEVMAR
 real(gwp) :: rerot,imrot,retemp,imtemp
#endif
 complex(gwpc),pointer :: vc_sqrt(:)
 logical :: ltest,is_sus,is_scr,is_scr_p,skip,same_freqs,calc_epsilon,only_diag
 character(len=1) :: ans
 character(len=4) :: tagq
 character(len=24) :: codename
 character(len=500) :: msg
 character(len=fnlen) :: fname_out,fname,fname_dump,fname_rho,prefix,fname_eigen,intp_fname
 character(len=fnlen) :: fname_dump2
 type(ScrHdr_type) :: Hscr_recovery
 type(ScrHdr_type),pointer :: Hscr0
 type(ScrHdr_type),target :: Hscr_merge
 type(MPI_type) :: MPI_enreg_seq
 type(BZ_mesh_type) :: Kmesh,Qmesh
 type(Crystal_structure) :: Cryst
 type(Gpairs_type) :: Gpairs_q
 type(Gvectors_type)  :: Gsphere
 type(PPmodel_type) :: PPm
 type(Epsilonm1_results) :: Er
 type(vcoul_t) :: Vcp
 type(Dataset_type) :: Dtset
 type(Datafiles_type) :: Dtfil
 type(MPI_type) :: mpi_enreg
!arrays
 integer :: ngfft(18) !,ipvt(100)
 integer,allocatable :: merge_table(:,:),foundq(:),freq_indx(:,:),ifile_indx(:)
 integer,allocatable :: pos_indx(:),i_temp(:),i2_temp(:,:)
!integer,pointer :: ip2fp(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),qdiff(3),rmet(3,3),qtmp(3)
 real(dp)  ,allocatable :: qlwl(:,:),real_omega(:)
 real(dp)  ,allocatable :: rhor(:,:),rhog(:,:),nhat(:,:)
 real(dp)  ,allocatable :: work(:),ftab(:),ysp(:,:),eint(:),qratio(:,:)
 real(gwp) ,allocatable :: epsm1_pole(:,:,:,:)
#ifdef HAVE_ALGO_LEVMAR
 real(gwp) ,allocatable :: coeff(:),startcoeff(:)
 real(gwp) ,allocatable :: refval(:),imfval(:),refval2(:),imfval2(:)
#endif
 complex(gwpc),allocatable :: epsm1(:,:,:,:),epsm1_temp(:,:,:,:),kxcg(:,:)
 complex(gwpc),allocatable :: fval(:)
 complex(dpc),allocatable :: omega_storage(:),omega_new(:),omega(:)
 complex(dpc),allocatable :: em1_ppm(:),epsm1_eigen(:,:),ppm_eigen(:,:)
 complex(dpc),allocatable :: rhoggp(:,:)
 character(len=fnlen),allocatable :: filenames(:)
 type(ScrHdr_type),target,allocatable :: Hscr_file(:)

 integer :: old_qptrlatt(3,3)
 real(dp),allocatable :: new_kibz(:,:),old_shiftq(:,:)

! *************************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

 call xmpi_init()

!Initialize MPI : one should write a separate routine -init_mpi_enreg- for doing that.
!Default for sequential use
 call nullify_mpi_enreg(mpi_enreg)
 mpi_enreg%world_comm=0
 mpi_enreg%world_group=0
 mpi_enreg%me=0
 mpi_enreg%nproc=1
 mpi_enreg%nproc_atom=1
 mpi_enreg%num_group_fft = 0 ! in some cases not initialized but referenced in xdef_comm.F90
 mpi_enreg%paral_compil=0
 mpi_enreg%paral_compil_fft=0
 mpi_enreg%paral_compil_mpio=0
 mpi_enreg%mode_para="n"
 mpi_enreg%flag_ind_kg_mpi_to_seq = 0
!MG080916 If we want to avoid MPI preprocessing options, %proc_distr should be always allocated and
!set to mpi_enreg%me. In such a way we can safely test its value inside loops parallelized over k-points
!For the time being, do not remove this line since it is needed in outkss.F90.
 nullify(mpi_enreg%proc_distrb)
 nullify(mpi_enreg%bandfft_kpt,mpi_enreg%tab_kpt_distrib)

!Initialize MPI
#if defined HAVE_MPI

 mpi_enreg%world_comm=xmpi_world
 mpi_enreg%world_group=MPI_GROUP_NULL
 call MPI_COMM_RANK(xmpi_world,mpi_enreg%me,ierr)
 call MPI_COMM_SIZE(xmpi_world,mpi_enreg%nproc,ierr)
!write(std_out,*)' mrgscr : nproc,me=',mpi_enreg%nproc,mpi_enreg%mempi_enreg%paral_compil=1
 mpi_enreg%paral_compil_respfn=0
 mpi_enreg%paral_level=2
#endif

!Signal MPI I/O compilation has been activated
#if defined HAVE_MPI_IO
 mpi_enreg%paral_compil_mpio=1
!if(mpi_enreg%paral_compil==0)then
!write(msg,'(3a)')&
!&   ' In order to use MPI_IO, you must compile with the MPI flag ',ch10,&
!&   ' Action : recompile your code with different CPP flags.'
!MSG_ERROR(msg)
!end if
#endif

 spaceComm=mpi_enreg%world_comm
!Initialize paral_compil_kpt, actually always equal to paral_compil
!(paral_compil_kpt should be suppressed after big cleaning)
 mpi_enreg%paral_compil_kpt=0
 if(mpi_enreg%paral_compil==1) mpi_enreg%paral_compil_kpt=1

 is_sus=.FALSE.; is_scr=.FALSE.; is_scr_p=.FALSE.

!Other values of mpi_enreg are dataset dependent, and should NOT be initialized
!inside mrgscr.F90.


!* Init fake MPI type with values for sequential case.
 call initmpi_seq(MPI_enreg_seq)

!=== Write greetings, and read the number of files ===
 codename='MRGSCR'//REPEAT(' ',18)
 call herald(codename,abinit_version,std_out)
!YP: calling dump_config() makes tests fail => commented
!call dump_config()

 call prompt(' Enter the number of files to merge: ',nfiles)
!ABI_CHECK(nfiles>0,'nfiles must be >0')
 if (nfiles<0 ) then
   call my_GWannier()
 end if

 ABI_ALLOCATE(filenames,(nfiles))
 ABI_ALLOCATE(Hscr_file,(nfiles))

 if (nfiles==1) then
   call prompt(' Enter the name of the file to be analyzed: ',filenames(1))
   write(msg,'(7a)')ch10,&
&   ' Running single-file mode:',ch10,&
&   ' Checking the integrity of file: ',TRIM(filenames(1)),ch10,&
&   ' reporting the list of q-points that are missing. '
   call wrtout(std_out,msg,'COLL')

 else if (nfiles>1) then ! Read name of files to be merged and check for existence.
   call prompt(' Enter the prefix for the final output file: ',fname_out)

   do ifile=1,nfiles
     write(msg,'(a,i4)')' Enter the name for the partial screening file no.',ifile
     call prompt(msg,filenames(ifile))

     if (.not.file_exist(filenames(ifile))) then
       write(msg,'(3a)')' File ',TRIM(filenames(ifile)),' does not exist. '
       MSG_ERROR(msg)
     end if
   end do
 end if
!
!=== Read the header of each file ===
 do ifile=1,nfiles

   unitem1=get_unit()
   open(unit=unitem1,file=filenames(ifile),status='old',form='unformatted',iostat=ios)
   ABI_CHECK(ios==0,' Opening '//TRIM(filenames(ifile)))

!  TODO this should initialized by the user
   accesswff=IO_MODE_FORTRAN
   rdwr=1

   call scr_hdr_io(fform1,rdwr,unitem1,spaceComm,master,accesswff,Hscr_file(ifile))

   if (.not.((fform1==1002).or.(fform1==1102).or.fform1==2002)) then
     write(msg,'(a,i5)')'Error while reading ScrHdr, fform/=1002,2002 or 1102, ',fform1
     MSG_ERROR(msg)
   end if
   if (fform1==1102) is_sus=.TRUE.
   if (fform1==1002) is_scr=.TRUE.
   if (fform1==2002) is_scr_p=.TRUE.
   if (is_sus.AND.is_scr) then
     MSG_ERROR('Error, files seem to be of mixed type (_SUS *and* _SCR)')
   end if

!  rdwr=4; call scr_hdr_io(fform,rdwr,std_out,spaceComm,master,accesswff,Hscr_file(ifile))
   close(unitem1)

   call print_ScrHdr(Hscr_file(ifile),unit=std_out,prtvol=1)

   if (ifile==1) then
     call metric(gmet,gprimd,-1,rmet,Hscr_file(ifile)%Hdr%rprimd,ucvol)
   end if
 end do !ifile
 if (nfiles>1) then
!  
!  Put the correct ending on the output file
   if (is_scr) fname_out=TRIM(fname_out)//'_SCR'
   if (is_sus) fname_out=TRIM(fname_out)//'_SUS'
   if (is_scr_p) fname_out=TRIM(fname_out)//'_SCR'
   call isfile(fname_out,'new')
 end if
!============================
!=== Merge multiple files ===
!============================
 if (nfiles>1) then

   if (is_scr_p) then
     write(std_out,'(2(a))') ch10,' Option not available for pole-fits! Exiting ...'
     goto 100
   end if

!  Check what kind of merging is to be performed
   write(std_out,'(2(a))') ch10,' Do you want to merge q-points        (= 1) ?'
   write(std_out,'(a)')         '  or do you want to merge frequencies (= 2) ?'
   read(std_in,*)choice

   select case(choice)
     case(1)

       write(std_out,'(3a)') ch10,' 1 => merging q-points',ch10

!      * Merge the headers creating the full list of q-points.
       call merge_ScrHdr(Hscr_file,Hscr_merge)
       call print_ScrHdr(Hscr_merge,header='Header of the final file',unit=std_out,prtvol=1)

!      For each q to be merged, save the index of the file where q is stored as well as its sequential index.
!      Useful to do the merge point-by-point thus avoiding the allocation of the entire epsm1 array.
       ABI_ALLOCATE(merge_table,(Hscr_merge%nqibz,2))
       do iqibz=1,Hscr_merge%nqibz
         ifound=0
         fl: do ifile=1,nfiles
           do iqf=1,Hscr_file(ifile)%nqibz
             qdiff(:)=Hscr_merge%qibz(:,iqibz)-Hscr_file(ifile)%qibz(:,iqf)
             if (normv(qdiff,gmet,'G')<GW_TOLQ) then
               merge_table(iqibz,1)=ifile
               merge_table(iqibz,2)=iqf
               ifound=ifound+1
               write(msg,'(a,3f12.6,2a)')' q-point :',Hscr_merge%qibz(:,iqibz),' will be taken from ',TRIM(filenames(ifile))
               call wrtout(std_out,msg,'COLL')
               EXIT fl
             end if
           end do
         end do fl
!        Check if point has been found, multiple points not allowed.
         ABI_CHECK(ifound==1,'ifound/=1')
       end do

       unt_out=get_unit()
       open(unit=unt_out,file=fname_out,status='new',form='unformatted',iostat=ios)
       if (ios/=0) then
         write(msg,'(3a)')' Opening file ',TRIM(fname_out),' as new-unformatted'
         MSG_ERROR(msg)
       end if

!      * Write the header.
!      TODO Use new format but first of all fix problem with hdr_check
       rdwr=2
       if (is_scr) fform_merge=1002
       if (is_sus) fform_merge=1102
       call scr_hdr_io(fform_merge,rdwr,unt_out,spaceComm,master,accesswff,Hscr_merge)

       npwe4m   = Hscr_merge%npwe
       nomega4m = Hscr_merge%nomega

       ABI_ALLOCATE(epsm1,(npwe4m,npwe4m,nomega4m,1))
       istat = ABI_ALLOC_STAT
       ABI_CHECK(istat==0,'out of memory in epsm1')

       do iqibz=1,Hscr_merge%nqibz
         ifile=merge_table(iqibz,1)
         iqiA =merge_table(iqibz,2)
         fname=filenames(ifile)
         call read_screening(fname,npwe4m,1,nomega4m,epsm1,accesswff,spaceComm,iqiA=iqiA)
         call write_screening(unt_out,accesswff,npwe4m,nomega4m,epsm1)
       end do

       ABI_DEALLOCATE(epsm1)
       ABI_DEALLOCATE(merge_table)
       close(unt_out)
       write(msg,'(3a)')ch10,' ==== Files have been merged successfully === ',ch10
       call wrtout(std_out,msg,'COLL')

     case(2) ! Merge frequencies

       write(std_out,'(3a)') ch10,' 2 => merging frequency grids',ch10
       write(std_out,'(2a)') ' WARNING: Advanced user option, consistency in fform,',&
&       ' etc. will not be checked.'

!      Check that q-point sets are the same
       do ifile=1,nfiles
         do ii=1,nfiles
           if (ii==ifile) CYCLE
           if (Hscr_file(ifile)%nqibz /= Hscr_file(ii)%nqibz) then
             MSG_ERROR(' One or more files do not have the same number of q-points!')
           end if
           do iqibz=1,Hscr_file(1)%nqibz
             if (ABS(SUM(Hscr_file(ifile)%qibz(:,iqibz)-Hscr_file(ii)%qibz(:,iqibz))) > tol6) then
               MSG_ERROR(' Q-point set differs between one or more files!')
             end if
           end do
         end do
       end do

!      Find out which file, if any, holds the imaginary frequencies to merge
       indx_imfreq_file=0
       write(std_out,'(2a)') ch10,' Enter index of file to take imaginary frequencies from:'
       write(std_out,'(a)')       '  (enter 0 if no imaginary frequencies are to be merged)'
       do ifile=1,nfiles
         write(std_out,'(a,I0,2a)')    '(',ifile,') => ',TRIM(filenames(ifile))
       end do
       read(std_in,*)indx_imfreq_file

       write(std_out,'(2a)') ch10,' Enter freqremax [eV] for the merged file (Enter 0 to use all freq. found):'
       read(std_in,*)freqremax
       freqremax = freqremax/Ha_eV
       if (freqremax<tol16) freqremax = HUGE(freqremax)

!      nfreq_tot here is the total *possible* number of freq.
       nfreq_tot=0
       do ifile=1,nfiles
         nfreq_tot = nfreq_tot + Hscr_file(ifile)%nomega
       end do
       ABI_ALLOCATE(omega_storage,(nfreq_tot))
       ABI_ALLOCATE(freq_indx,(nfreq_tot,nfiles))
       ABI_ALLOCATE(ifile_indx,(nfreq_tot))
       omega_storage = CMPLX(-one,-one); freq_indx = 0; ifile_indx = 0

!      Calculate the total number of real freq and store
       nfreqre = 0
       do ifile=1,nfiles
         do ifrq=1,Hscr_file(ifile)%nomega
           skip = .FALSE.
!          Check whether to skip this point
           if (AIMAG(Hscr_file(ifile)%omega(ifrq)) > tol16) skip = .TRUE.
           if (REAL(Hscr_file(ifile)%omega(ifrq)) > freqremax) skip = .TRUE.
!          Check for repetition or non-monotonic points
           if (nfreqre>1) then
             do ii=1,nfreqre
               if (ABS(REAL(Hscr_file(ifile)%omega(ifrq)) - REAL(omega_storage(ii))) < tol6) skip = .TRUE.
             end do
           end if
           if (skip) CYCLE
           nfreqre = nfreqre + 1
!          Store (complex) frequency and index
           omega_storage(nfreqre) = Hscr_file(ifile)%omega(ifrq)
           ifile_indx(nfreqre) = ifile
           freq_indx(nfreqre,ifile) = ifrq
           write(std_out,'(a,F12.6,a,I0)') ' Found frequency: ',REAL(omega_storage(nfreqre)),&
&           ' number: ',nfreqre
           write(std_out,'(2(a,I0))') ' in file: ',ifile,' local index: ',ifrq
         end do
       end do

!      Sort real frequencies
       ABI_ALLOCATE(real_omega,(nfreqre))
       ABI_ALLOCATE(pos_indx,(nfreqre))
       ABI_ALLOCATE(i_temp,(nfreqre))
       ABI_ALLOCATE(i2_temp,(nfreqre,nfiles))
       real_omega(1:nfreqre) = REAL(omega_storage(1:nfreqre)) ! Copy real frequencies to temp. sorting array
       do ii=1,nfreqre ! Set up indexing array
         pos_indx(ii) = ii
       end do
!      Sort frequencies while keeping track of index
       call sort_dp(nfreqre,real_omega,pos_indx,tol16)
       i_temp(1:nfreqre) = ifile_indx(1:nfreqre)
       i2_temp(1:nfreqre,1:nfiles) = freq_indx(1:nfreqre,1:nfiles)
!      Copy sorted frequencies plus file and frequency index
       do ii=1,nfreqre
         omega_storage(ii) = CMPLX(real_omega(ii),zero)
         ifile_indx(ii) = i_temp(pos_indx(ii))
         freq_indx(ii,1:nfiles) = i2_temp(pos_indx(ii),1:nfiles)
       end do
       ABI_DEALLOCATE(real_omega)
       ABI_DEALLOCATE(pos_indx)
       ABI_DEALLOCATE(i_temp)
       ABI_DEALLOCATE(i2_temp)

!      Check imaginary frequencies and store them
       nfreqim = 0
       if (indx_imfreq_file/=0) then
         ifile = indx_imfreq_file
         do ifrq=1,Hscr_file(ifile)%nomega
           if (REAL(Hscr_file(ifile)%omega(ifrq)) > tol8) CYCLE
           if (AIMAG(Hscr_file(ifile)%omega(ifrq)) < tol8) CYCLE
           nfreqim = nfreqim + 1
           omega_storage(nfreqre+nfreqim) = Hscr_file(ifile)%omega(ifrq)
           ifile_indx(nfreqre+nfreqim) = ifile
           freq_indx(nfreqre+nfreqim,ifile) = ifrq
           write(std_out,'(a,F12.6,a,I0)') ' Found imaginary frequency: ',AIMAG(omega_storage(nfreqre+nfreqim)),&
&           ' number: ',nfreqim
           write(std_out,'(2(a,I0))') ' in file: ',ifile,' local index: ',ifrq
         end do
       end if
       nfreq_tot = nfreqre + nfreqim ! Here nfreq_tot becomes the *true* number of freq
       write(std_out,'(2a,I0,a)') ch10,' Merging ',nfreq_tot,' frequencies.'
       write(std_out,'(2(a,I0),2a)') ' ',nfreqre,' real, and ',nfreqim,' imaginary.',ch10

!      First copy an old header
       ios = 1; if (ios==indx_imfreq_file) ios = 2
       if (is_sus) then ! Hack to avoid bug for SUS files
         ii = Hscr_file(ios)%headform
         Hscr_file(ios)%headform = 56
       end if
       call copy_ScrHdr(Hscr_file(ios),Hscr_merge)
       if (is_sus) then ! Hack to avoid bug for SUS files
         Hscr_file(ios)%headform = ii
         Hscr_merge%headform = ii
       end if


!      Then modify entries for new frequency grid
       Hscr_merge%nomega = nfreq_tot
       ABI_DEALLOCATE(Hscr_merge%omega)
       ABI_ALLOCATE(Hscr_merge%omega,(nfreq_tot))
       Hscr_merge%omega(1:nfreq_tot) = omega_storage(1:nfreq_tot)
!      Check if q=0 is part of the set and reallocate and set lwing and uwing
       do iqibz=1,Hscr_merge%nqibz
         if (Hscr_merge%qibz(1,iqibz)<tol6.AND.&
&         Hscr_merge%qibz(2,iqibz)<tol6.AND.&
&         Hscr_merge%qibz(3,iqibz)<tol6) then
           indx_q0 = iqibz
         end if
         npwe4mI = Hscr_merge%npwe*Hscr_merge%nI
         npwe4mJ = Hscr_merge%npwe*Hscr_merge%nJ
         ABI_DEALLOCATE(Hscr_merge%uwing)
         ABI_DEALLOCATE(Hscr_merge%lwing)
         ABI_ALLOCATE(Hscr_merge%uwing,(npwe4mJ,nfreq_tot,Hscr_merge%nqlwl))
         ABI_ALLOCATE(Hscr_merge%lwing,(npwe4mI,nfreq_tot,Hscr_merge%nqlwl))
         do ifrq=1,nfreq_tot
           Hscr_merge%uwing(:,ifrq,:)=&
&           Hscr_file(ifile_indx(ifrq))%uwing(:,freq_indx(ifrq,ifile_indx(ifrq)),:)
           Hscr_merge%lwing(:,ifrq,:)=&
&           Hscr_file(ifile_indx(ifrq))%lwing(:,freq_indx(ifrq,ifile_indx(ifrq)),:)
         end do
       end do

!      Print new header for info
       call print_ScrHdr(Hscr_merge,header='Header of the final file',unit=std_out,prtvol=1)

       unt_out=get_unit()
       open(unit=unt_out,file=fname_out,status='new',form='unformatted',iostat=ios)
       if (ios/=0) then
         write(msg,'(3a)')' Opening file ',TRIM(fname_out),' as new-unformatted'
         MSG_ERROR(msg)
       end if

!      * Write the header.
!      TODO Use new format but first of all fix problem with hdr_check
       rdwr=2
       if (is_scr) fform_merge=1002
       if (is_sus) fform_merge=1102
       call scr_hdr_io(fform_merge,rdwr,unt_out,spaceComm,master,accesswff,Hscr_merge)

       npwe4mI = Hscr_merge%npwe*Hscr_merge%nI
       npwe4mJ = Hscr_merge%npwe*Hscr_merge%nJ
       nomega4m = Hscr_merge%nomega
       write(std_out,'(a,f16.1,a)')' Memory required for merged matrix   = ',&
&       2.0*gwpc*npwe4mI*npwe4mJ*nomega4m*b2Mb," [Mb]."
       ABI_ALLOCATE(epsm1,(npwe4mI,npwe4mJ,nomega4m,1))
       istat = ABI_ALLOC_STAT
       ABI_CHECK(istat==0,'out of memory in epsm1')

       do iqibz=1,Hscr_merge%nqibz
         do ifile=1,nfiles

!          allocate temporary array
           npwe4mI = Hscr_file(ifile)%npwe*Hscr_file(ifile)%nI
           npwe4mJ = Hscr_file(ifile)%npwe*Hscr_file(ifile)%nJ
           nomega4m = Hscr_file(ifile)%nomega
           write(std_out,'(a,f16.1,a)')' Memory required for temporary matrix= ',&
&           2.0*gwpc*npwe4mI*npwe4mJ*nomega4m*b2Mb,' [Mb].'
           ABI_ALLOCATE(epsm1_temp,(npwe4mI,npwe4mJ,nomega4m,1))
           istat = ABI_ALLOC_STAT
           ABI_CHECK(istat==0,'out of memory in epsm1_temp')

!          read screening
           fname=filenames(ifile)
           call read_screening(fname,npwe4mI,1,nomega4m,epsm1_temp,accesswff,spaceComm,iqiA=iqibz)

!          Copy matrices for relevant frequencies
           do ifrq=1,nfreqre
             if (ifile_indx(ifrq)==ifile) then
               epsm1(:,:,ifrq,1)=epsm1_temp(:,:,freq_indx(ifrq,ifile),1)
             end if
           end do

!          Add imaginary part if needed
           if (indx_imfreq_file==ifile) then
             do ifrq=nfreqre+1,Hscr_merge%nomega
               epsm1(:,:,ifrq,1)=epsm1_temp(:,:,freq_indx(ifrq,ifile),1)
             end do
           end if

           ABI_DEALLOCATE(epsm1_temp)

         end do !ifile

         npwe4mI = Hscr_merge%npwe*Hscr_merge%nI
         nomega4m = Hscr_merge%nomega
         call write_screening(unt_out,accesswff,npwe4mI,nomega4m,epsm1)

       end do ! iqibz

       ABI_DEALLOCATE(epsm1)
       close(unt_out)
       write(msg,'(3a)')ch10,' ==== Files have been merged successfully === ',ch10
       call wrtout(std_out,msg,'COLL')

       ABI_DEALLOCATE(omega_storage)
       ABI_DEALLOCATE(freq_indx)
       ABI_DEALLOCATE(ifile_indx)

       case default

       write(std_out,*) ' Invalid choice! Exiting...'
       goto 100

   end select

 end if ! nfiles>1

!=== Now check if the list of q-points is complete ===
!* Here we assume that the k-mesh reported in the header is the same as
!that used during the sigma calculation.
 write(msg,'(3a)') ch10,' Checking if the list of q-points is complete. ',ch10
 call wrtout(std_out,msg,'COLL')

 Hscr0 => Hscr_file(1)
 fname =filenames(1)
 if (nfiles>1) then
   Hscr0 => Hscr_merge
   fname =fname_out
 end if

 timrev=2 ! This should be read from kptopt
 call init_crystal_from_hdr(Cryst,HScr0%Hdr,timrev,remove_inv=.FALSE.)

 kptopt=1
 call init_kmesh(Kmesh,Cryst,HScr0%Hdr%nkpt,Hscr0%Hdr%kptns,kptopt)

 call print_BZ_mesh(Kmesh,"K-mesh for the wavefunctions",prtvol=prtvol)

 call find_qmesh(Qmesh,Cryst,Kmesh)

 call print_BZ_mesh(Qmesh,"Q-mesh for the screening function",prtvol=prtvol)

 ABI_ALLOCATE(foundq,(Qmesh%nibz))
 foundq(:)=0
 do iqibz=1,Qmesh%nibz
   do iqf=1,Hscr0%nqibz
     qdiff(:)=Qmesh%ibz(:,iqibz)-Hscr0%qibz(:,iqf)
     if (normv(qdiff,gmet,'G')<GW_TOLQ) foundq(iqibz)=foundq(iqibz)+1
   end do
 end do

 if (ANY(foundq(:)==0)) then
   write(msg,'(6a)')ch10,&
&   ' File ',TRIM(fname),' is not complete ',ch10,&
&   ' The following q-points are missing :'
   call wrtout(std_out,msg,'COLL')
   ii=0
   do iqibz=1,Qmesh%nibz
     if (foundq(iqibz)==0) then
       ii=ii+1
       write(msg,'(i3,a,3f12.6)')ii,') ',Qmesh%ibz(:,iqibz)
       call wrtout(std_out,msg,'COLL')
     end if
   end do
 end if

 if (ANY(foundq(:)>1)) then
   write(msg,'(6a)')ch10,&
&   ' File ',TRIM(fname),' is overcomplete ',ch10,&
&   ' The following q-points are present more than once :'
   call wrtout(std_out,msg,'COLL')
   ii=0
   do iqibz=1,Qmesh%nibz
     if (foundq(iqibz)>1) then
       ii=ii+1
       write(msg,'(i3,a,3f12.6)')ii,') ',Qmesh%ibz(:,iqibz)
       call wrtout(std_out,msg,'COLL')
     end if
   end do
 end if

 if (ALL(foundq(:)==1)) then
   write(msg,'(5a)')ch10,&
&   ' File ',TRIM(fname),' contains a complete list of q-points ',ch10
   call wrtout(std_out,msg,'COLL')
 end if

!=====================
!=== Recovery mode ===
!=====================
 if (nfiles==1) then

   write(std_out,'(2(a))') ch10,' Do you want to recover a subset of q-points    (= 1) ?'
   write(std_out,'(a)')         '  or extract the contents of the file           (= 2) ?'
   write(std_out,'(a)')         '  or create dielectric function (SCR file)'
   write(std_out,'(a)')         '    and/or extract plasmon-pole parameters      (= 3) ?'
   write(std_out,'(a)')         '  or remove real frequencies                    (= 4) ?'
   write(std_out,'(a)')         '  or remove imaginary frequencies               (= 5) ?'
   write(std_out,'(a)')         '  or calculate a model screening                (= 6) ?'
!  write(std_out,'(a)')         '  or interpolate a new real freq. grid          (= 7) ?'
!  write(std_out,'(a)')         '  or interpolate the screening in k-space       (= 8) ?'
   read(std_in,*)choice

   select case(choice)
     case(1) ! Recover subset of q-points --------------------------------------------------

       if (is_scr_p) then
         write(std_out,'(2(a))') ch10,' Option not available for pole-fits! Exiting ...'
         goto 100
       end if

       write(std_out,'(a)') ' 1 => Recovering subset of q-points'
       call prompt(' Enter the number of q-points to be extracted: ',nq_selected)
       ltest = (nq_selected>0 .and. nq_selected<=Hscr0%nqibz)
       ABI_CHECK(ltest,"Wrong number of q-points")

       call prompt(' Enter the name of the final output file: ',fname_out)

       unt_tmp = get_unit()
       open(unit=unt_tmp,file=fname_out,status='new',form='unformatted',iostat=ios)
       msg = ' Opening file '//TRIM(fname_out)//' as new-unformatted'
       ABI_CHECK(ios==0,msg)

!      Initialize new header with correct dimensions.
       if (is_sus) then ! Hack to avoid bug for SUS files
         ii = Hscr_file(1)%headform
         Hscr_file(1)%headform = 56
       end if
       call copy_ScrHdr(Hscr0,Hscr_recovery)
       if (is_sus) then ! Hack to avoid bug for SUS files
         Hscr_file(1)%headform = ii
         Hscr_merge%headform = ii
       end if

       call print_ScrHdr(Hscr_recovery,header="Header of the new SCR file",unit=std_out,prtvol=1)

!      Had to change dimensions and arrays associated to nqibz.
       Hscr_recovery%nqibz = nq_selected
       ABI_DEALLOCATE(Hscr_recovery%qibz)
       ABI_ALLOCATE(Hscr_recovery%qibz,(3,nq_selected))
       Hscr_recovery%qibz = Hscr0%qibz(:,1:nq_selected)

!      Write the new header of the recovered file.
!      rdwr=2; fform_chi0=1102 ! Use the new format
!      TODO Use new format but first of all fix problem with hdr_check
       rdwr=2
       if (is_scr) fform1=1002
       if (is_sus) fform1=1102
       call scr_hdr_io(fform1,rdwr,unt_tmp,spaceComm,master,accesswff,Hscr_recovery)

       nqibzA=1; nomega_asked=Hscr0%nomega; npwe_asked=Hscr0%npwe

       ABI_ALLOCATE(epsm1,(npwe_asked,npwe_asked,nomega_asked,1))
       istat = ABI_ALLOC_STAT
       ABI_CHECK(istat==0,'out of memory in epsm1')

       do iqiA=1,Hscr_recovery%nqibz
         call read_screening(fname,npwe_asked,nqibzA,nomega_asked,epsm1,accesswff,spaceComm,iqiA=iqiA)
         call write_screening(unt_tmp,accesswff,npwe_asked,nomega_asked,epsm1)
       end do

       call free_scrhdr(Hscr_recovery)
       ABI_DEALLOCATE(epsm1)
       close(unt_tmp)

       call wrtout(std_out,"Recovery completed",'COLL')

     case(2) ! Analyse file ----------------------------------------------------------------
       
       write(std_out,'(a)') ' 2 => Extraction of file contents' 

       select case(Hscr0%fform)
         case(1002,1102)
           
!          === Initialize the G-sphere ===
           call init_gsphere(Gsphere,.FALSE.,Cryst,Hscr0%npwe,gvec=Hscr0%gvec)

           call nullify_Gpairs_type(Gpairs_q)

           ABI_ALLOCATE(epsm1,(Hscr0%npwe,Hscr0%npwe,Hscr0%nomega,1))
           istat = ABI_ALLOC_STAT
           if (istat/=0) STOP 'out of memory in epsm1'

!          Give option to output epsilon instead of chi0
           calc_epsilon = .FALSE.
           if (is_sus) then
             write(std_out,'(2a)') ch10,&
&             ' You have provided a chi_0 file for analysis. Would you like to output'
             write(std_out,'(2a)',advance='no') ' the dielectric function epsilon_GG'' ',&
&             '= delta_GG'' - v_G*chi0_GG''[Y/N] ? '
             read(std_in,*)ans
             if (ans=='Y'.or.ans=='y') then
!              Initialise Coulomb terms
               if (Er%Hscr%nqlwl==0) then
                 nqlwl=1
                 ABI_ALLOCATE(qlwl,(3,nqlwl))
                 qlwl(:,1)= GW_Q0_DEFAULT
               else
                 nqlwl=Er%Hscr%nqlwl
                 ABI_ALLOCATE(qlwl,(3,nqlwl))
                 qlwl(:,:)=Er%Hscr%qlwl(:,1:nqlwl)
               end if

               Dtset%icutcoul=3; Dtset%rcut=zero
               Dtset%vcutgeo=(/zero,zero,zero/);
               Dtset%boxcenter=(/zero,zero,zero/)

               write(std_out,'(2a)',advance='no') ch10,&
&               ' Was a Coulomb cutoff technique used [Y/N] ? '
               read(std_in,*)ans
               if (ans=='Y'.or.ans=='y') then
                 write(std_out,'(2a)',advance='no') ' Enter icutcoul: '
                 read(std_in,*)Dtset%icutcoul
                 write(std_out,'(2a)',advance='no') ' Enter vcutgeo: '
                 read(std_in,*)Dtset%vcutgeo
                 write(std_out,'(2a)',advance='no') ' Enter boxcenter: '
                 read(std_in,*)Dtset%boxcenter
               end if

               call vcoul_init(Vcp,Gsphere,Qmesh,Kmesh,Dtset%rcut,Dtset%icutcoul,&
&               Dtset%vcutgeo,Hscr0%npwe,nqlwl,qlwl,Cryst%rprimd,ngfft,spaceComm)
               ABI_DEALLOCATE(qlwl)

               calc_epsilon = .TRUE.
             end if
           end if

           ig1 = 0; ig2 = 0
           write(std_out,'(2(a),I0,a)',advance='NO') ch10,' Enter the starting index for G (1 - ',Hscr0%npwe,' ): '
           read(std_in,*)ig1_start
           if (ig1_start<1.OR.ig1_start>Hscr0%npwe) then
             MSG_ERROR(' Starting index out of bounds')
           end if
           write(std_out,'(a,I0,a,I0,a)',advance='NO')    ' Enter the ending index for G ( ',ig1_start,' - ',Hscr0%npwe,' ): '
           read(std_in,*)ig1_end
           if (ig1_end<ig1_start.OR.ig1_end>Hscr0%npwe) then
             MSG_ERROR(' Ending index out of bounds')
           end if
           write(std_out,'(a,I0,a)',advance='NO')         ' Enter the starting index for G'' (1 - ',Hscr0%npwe,' ): '
           read(std_in,*)ig2_start
           if (ig2_start<1.OR.ig2_start>Hscr0%npwe) then
             MSG_ERROR(' Starting index out of bounds')
           end if
           write(std_out,'(a,I0,a,I0,a)',advance='NO')    ' Enter the ending index for G'' ( ',ig2_start,' - ',Hscr0%npwe,' ): '
           read(std_in,*)ig2_end
           if (ig2_end<ig2_start.OR.ig2_end>Hscr0%npwe) then
             MSG_ERROR(' Ending index out of bounds')
           end if

           only_diag = .FALSE.
           write(std_out,'(a)',advance='no') ' Would you like to output only the diagonal [Y/N] ? '
           read(std_in,*)ans
           if (ans=='Y'.or.ans=='y') only_diag = .TRUE.

           do iqibz=1,Hscr0%nqibz

!            === Find the independent set of G-Gp pairs for this q-point. ===
!            * In the long wavelength limit we set q==0, because we still can use symmetries for the Body.
             qtmp(:)=Hscr0%qibz(:,iqibz); if (normv(qtmp,Cryst%gmet,'G')<GW_TOLQ0) qtmp(:)=zero
             call init_Gpairs_type(Gpairs_q,qtmp,Gsphere,Cryst)

             call read_screening(fname,Hscr0%npwe,1,Hscr0%nomega,epsm1,accesswff,spaceComm,iqiA=iqibz)

             if (calc_epsilon) then ! Calculate epsilon
               do iomega=1,Hscr0%nomega
                 if (iqibz==1) then
                   if (nqlwl>1) then
                     MSG_ERROR('nqlwl>1 not coded yet!')
                   end if
                   vc_sqrt => Vcp%vcqlwl_sqrt(:,iqibz)  ! Use Coulomb term for q-->0
                 else
                   vc_sqrt => Vcp%vc_sqrt(:,iqibz)
                 end if
                 do ig2=ig2_start,ig2_end
                   do ig1=ig1_start,ig1_end
                     epsm1(ig1,ig2,iomega,1) = -(vc_sqrt(ig1)**2)*epsm1(ig1,ig2,iomega,1)
                   end do ! ig1
                   epsm1(ig2,ig2,iomega,1) = one + epsm1(ig2,ig2,iomega,1)
                 end do ! ig2
               end do ! iomega
             end if ! Do we calculate epsilon

!            Find out the total number of frequencies along real/imaginary axes
!            and possibly in the z-plane
             nfreqre=0; nfreqim=0; nfreqc=0;
             do iomega=1,Hscr0%nomega
               if (ABS(REAL(Hscr0%omega(iomega)))<tol8.AND.&
&               ABS(AIMAG(Hscr0%omega(iomega)))<tol8) nfreqre = nfreqre + 1 
               if (ABS(REAL(Hscr0%omega(iomega)))>tol8.AND.&
&               ABS(AIMAG(Hscr0%omega(iomega)))<tol8) nfreqre = nfreqre + 1 
               if (ABS(REAL(Hscr0%omega(iomega)))<tol8.AND.&
&               ABS(AIMAG(Hscr0%omega(iomega)))>tol8) nfreqim = nfreqim + 1
             end do
             if (Hscr0%nomega-nfreqre-nfreqim/=0) then
               write(std_out,'(/,a)') ' WARNING: There are frequencies in the full complex plane.'
               write(std_out,'(a)')   '          The _SCR or _SUS file might not be suitable'
               write(std_out,'(a,/)') '          for self-energy calculations.'
               nfreqc = Hscr0%nomega-nfreqre-nfreqim
             end if
             write(std_out,'(2a,I0,a)') ch10,' Found ',Hscr0%nomega,' frequencies.'
             write(std_out,'(2(a,I0),2a)') ' ',nfreqre,' real, and ',nfreqim,' imaginary.',ch10
             if (nfreqc>0) then
               write(std_out,'(a,I0)') ' There is a grid in the complex plane with ',nfreqc
               write(std_out,'(2a)')   '  extra frequencies in the list.',ch10
             end if

!            Get Q index for name
             call int2char4(iqibz,tagq)

             if (nfreqre>0) then ! Output real frequency axis
               if (calc_epsilon) then
                 fname_dump=TRIM(fname)//'_EPS_Q'//TRIM(tagq); unt_dump=get_unit()
               else
                 fname_dump=TRIM(fname)//'_Q'//TRIM(tagq); unt_dump=get_unit()
               end if
               open(file=fname_dump,unit=unt_dump,status='replace',form='formatted',iostat=ios)
!              Check if number of pairs are reduced by symmetry
!              if (Gpairs_q%niggp<Gpairs_q%ng**2) then ! We have reduced pairs
!              idx=0
!              do igpair=1,Gpairs_q%niggp
!              if (npairs/=0) then
!              if (igpair>npairs.and.igpair<=Gpairs_q%niggp-npairs) CYCLE ! By default only first and last 50 pairs are printed
!              end if
!              ig1 = Gpairs_q%ip2fp(1,igpair)
!              ig2 = Gpairs_q%ip2fp(2,igpair)
!              idx=idx+1
!              write(unt_dump,'(a,i4,a,i8,/,a,3f12.6,/,a,3i6,a,3i6,/,a,/)')&
!              &             '# index= ',idx,'    pair number = ',igpair,&
!              &             '# q = ',Hscr0%qibz(:,iqibz),&
!              &             '# G = ',Hscr0%gvec(:,ig1),'  G''= ',Hscr0%gvec(:,ig2),&
!              &             '#   omega [eV]           Re             Im '
!              do iomega=1,Hscr0%nomega
!              if (ABS(AIMAG(Hscr0%omega(iomega)))>tol8) EXIT !only real frequencies
!              write(unt_dump,'(f8.2,4x,2es16.8)') REAL(Hscr0%omega(iomega))*Ha_eV,epsm1(ig1,ig2,iomega,1)
!              end do
!              write(unt_dump,*)
!              write(unt_dump,*)
!              end do
!              else ! We do not have reduced pairs
!              idx1=0
               do ig1=ig1_start,ig1_end
!                idx1=idx1+1
!                if (idx1>npairs.and.idx1<((Gpairs_q%ng**2)-npairs)) CYCLE
!                idx2=0
                 do ig2=ig2_start,ig2_end
                   if (only_diag.AND.ig1/=ig2) CYCLE
!                  idx2=idx2+1
!                  if (idx2>npairs.and.idx2<((Gpairs_q%ng**2)-npairs)) CYCLE
                   write(unt_dump,'(2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,/,a,/)')&
&                   '# ig1= ',ig1,'    ig2= ',ig2,&
&                   '# q = ',Hscr0%qibz(:,iqibz),&
&                   '# G = ',Hscr0%gvec(:,ig1),'  G''= ',Hscr0%gvec(:,ig2),&
&                   '#   omega [eV]           Re             Im '
                   do iomega=1,nfreqre
                     write(unt_dump,'(f8.2,4x,2es16.8)') REAL(Hscr0%omega(iomega))*Ha_eV,&
&                     REAL(epsm1(ig1,ig2,iomega,1)),AIMAG(epsm1(ig1,ig2,iomega,1))
                   end do
                   write(unt_dump,*)
                   write(unt_dump,*)
                 end do !ig2
               end do !ig1
!              end if
               close(unt_dump)
             end if ! Output real frequency axis

             if (nfreqim>0) then ! output imaginary frequency axis
               if (calc_epsilon) then
                 fname_dump=TRIM(fname)//'_EPS_Imfrq_Q'//TRIM(tagq); unt_dump=get_unit()
               else
                 fname_dump=TRIM(fname)//'_Imfrq_Q'//TRIM(tagq); unt_dump=get_unit()
               end if
               open(file=fname_dump,unit=unt_dump,status='replace',form='formatted',iostat=ios)
!              idx=0
!              if (Gpairs_q%niggp<Gpairs_q%ng**2) then ! We have reduced pairs
!              do igpair=1,Gpairs_q%niggp
!              if (npairs/=0) then
!              if (igpair>npairs.and.igpair<=Gpairs_q%niggp-npairs) CYCLE ! default npairs=50
!              end if
!              ig1 = Gpairs_q%ip2fp(1,igpair)
!              ig2 = Gpairs_q%ip2fp(2,igpair)
!              idx=idx+1
!              write(unt_dump,'(a,i4,a,i8,/,a,3f12.6,/,a,3i6,a,3i6,/,a,/)')&
!              &             '# index= ',idx,'    pair number = ',igpair,&
!              &             '# q = ',Hscr0%qibz(:,iqibz),&
!              &             '# G = ',Hscr0%gvec(:,ig1),'  G''= ',Hscr0%gvec(:,ig2),&
!              &             '#   iomega [eV]           Re             Im '
!              do iomega=1,Hscr0%nomega
!              if (ABS(REAL(Hscr0%omega(iomega)))>tol8) CYCLE !only imaginary frequencies
!              write(unt_dump,'(f8.2,4x,2es16.8)') AIMAG(Hscr0%omega(iomega))*Ha_eV,epsm1(ig1,ig2,iomega,1)
!              end do
!              write(unt_dump,*)
!              write(unt_dump,*)
!              end do
!              else
!              idx1=0
               do ig1=ig1_start,ig1_end
!                idx1=idx1+1
!                if (idx1>npairs.and.idx1<((Gpairs_q%ng**2)-npairs)) CYCLE
!                idx2=0
                 do ig2=ig2_start,ig2_end
                   if (only_diag.AND.ig1/=ig2) CYCLE
!                  idx2=idx2+1
!                  if (idx2>npairs.and.idx2<((Gpairs_q%ng**2)-npairs)) CYCLE
                   write(unt_dump,'(a,i4,2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,/,a,/)')&
&                   '# index= ',idx,'    ig1= ',ig1,'    ig2= ',ig2,&
&                   '# q = ',Hscr0%qibz(:,iqibz),&
&                   '# G = ',Hscr0%gvec(:,ig1),'  G''= ',Hscr0%gvec(:,ig2),&
&                   '#   omega [eV]           Re             Im '
                   do iomega=nfreqre+1,nfreqre+nfreqim
                     write(unt_dump,'(f8.2,4x,2es16.8)') AIMAG(Hscr0%omega(iomega))*Ha_eV,epsm1(ig1,ig2,iomega,1)
                   end do
                   write(unt_dump,*)
                   write(unt_dump,*)
                 end do !ig2
               end do !ig1
!              end if
               close(unt_dump)
             end if ! Check for imaginary frequencies

!            Check for complex plane values
             if (nfreqc>0) then
               if (calc_epsilon) then
                 fname_dump=TRIM(fname)//'_EPS_ZPLANE_Q'//TRIM(tagq); unt_dump=get_unit()
               else
                 fname_dump=TRIM(fname)//'_ZPLANE_Q'//TRIM(tagq); unt_dump=get_unit()
               end if
               open(file=fname_dump,unit=unt_dump,status='replace',form='formatted',iostat=ios)
               do ig1=ig1_start,ig1_end
                 do ig2=ig2_start,ig2_end
                   if (only_diag.AND.ig1/=ig2) CYCLE
                   write(unt_dump,'(a,i4,2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,/,a,/)')&
&                   '# index= ',idx,'    ig1= ',ig1,'    ig2= ',ig2,&
&                   '# q = ',Hscr0%qibz(:,iqibz),&
&                   '# G = ',Hscr0%gvec(:,ig1),'  G''= ',Hscr0%gvec(:,ig2),&
&                   '#   omega [eV]           Re             Im '
                   do iomega=1,nfreqre
                     write(unt_dump,'(2(f8.2),4x,2es16.8)') REAL(Hscr0%omega(iomega))*Ha_eV,&
&                     AIMAG(Hscr0%omega(iomega))*Ha_eV,epsm1(ig1,ig2,iomega,1)
                   end do
                   write(unt_dump,*)
                   do ios=1,nfreqim
                     do iomega=1,nfreqre
                       if (iomega==1) then
                         io = nfreqre + ios
                       else
                         io = nfreqre + nfreqim + (ios-1)*(nfreqre-1) + (iomega-1)
                       end if
                       write(unt_dump,'(2(f8.2),4x,2es16.8)') REAL(Hscr0%omega(io))*Ha_eV,&
&                       AIMAG(Hscr0%omega(io))*Ha_eV,epsm1(ig1,ig2,io,1)
                     end do
                     write(unt_dump,*)
                   end do
                   write(unt_dump,*)
                   write(unt_dump,*)
                 end do !ig2
               end do !ig1
               close(unt_dump) 
             end if ! Check for complex plane freqs

           end do !iqibz

           ABI_DEALLOCATE(epsm1)

           call destroy_Gpairs_type(Gpairs_q)
           call destroy_gsphere(Gsphere)

         case(2002) ! We are dealing with a pole-fit screening

!          === Initialize the G-sphere ===
           call init_gsphere(Gsphere,.FALSE.,Cryst,Hscr0%npwe,gvec=Hscr0%gvec)

           call nullify_Gpairs_type(Gpairs_q)

!          Specify grid
           write(std_out,'(2a)') ch10,' Enter new number of gridpoints along real axis:'
           read(std_in,*)nfreqre
           write(std_out,'(2a)') ch10,' Enter freqremin [in eV] for the grid:'
           read(std_in,*)freqremin
           freqremin = freqremin/Ha_eV
           write(std_out,'(2a)') ch10,' Enter freqremax [in eV] for the grid:'
           read(std_in,*)freqremax
           freqremax = freqremax/Ha_eV
           write(std_out,'(2a)') ch10,' Enter new number of gridpoints along imaginary axis:'
           read(std_in,*)nfreqim
           write(std_out,'(2a)') ch10,' Enter freqimmax [in eV] for the grid:'
           read(std_in,*)freqimmax
           freqimmax = freqimmax/Ha_eV

           if (nfreqre<1) then
             write(std_out,'(2(a))') ch10,' ERROR - One or more points need to be selected! Exiting ...'
             goto 100
           end if
           if (freqremin>freqremax) then
             write(std_out,'(2(a))') ch10,' ERROR - freqremin > freqremax! Exiting ...'
             goto 100
           end if
           if (freqremin<zero.OR.freqremax<zero.OR.freqimmax<zero) then
             write(std_out,'(2(a))') ch10,' ERROR - freqre(im)max or freqremin negative! Exiting ...'
             goto 100
           end if

           ABI_ALLOCATE(omega,(nfreqre*nfreqim))
           idx = 1
           dwre = (freqremax-freqremin)/real(nfreqre-1,dp)
           dwim = freqimmax/real(nfreqim-1,dp)
           do io=1,nfreqim
             do iomega=1,nfreqre
               omega(idx) = CMPLX(freqremin+(iomega-1)*dwre,(io-1)*dwim)
               idx = idx + 1
             end do
           end do

           ig1 = 0; ig2 = 0
           write(std_out,'(2(a),I0,a)',advance='NO') ch10,' Enter the starting index for G (1 - ',Hscr0%npwe,' ): '
           read(std_in,*)ig1_start
           if (ig1_start<1.OR.ig1_start>Hscr0%npwe) then
             MSG_ERROR(' Starting index out of bounds')
           end if
           write(std_out,'(a,I0,a,I0,a)',advance='NO')    ' Enter the ending index for G ( ',ig1_start,' - ',Hscr0%npwe,' ): '
           read(std_in,*)ig1_end
           if (ig1_end<ig1_start.OR.ig1_end>Hscr0%npwe) then
             MSG_ERROR(' Ending index out of bounds')
           end if
           write(std_out,'(a,I0,a)',advance='NO')         ' Enter the starting index for G'' (1 - ',Hscr0%npwe,' ): '
           read(std_in,*)ig2_start
           if (ig2_start<1.OR.ig2_start>Hscr0%npwe) then
             MSG_ERROR(' Starting index out of bounds')
           end if
           write(std_out,'(a,I0,a,I0,a)',advance='NO')    ' Enter the ending index for G'' ( ',ig2_start,' - ',Hscr0%npwe,' ): '
           read(std_in,*)ig2_end
           if (ig2_end<ig2_start.OR.ig2_end>Hscr0%npwe) then
             MSG_ERROR(' Ending index out of bounds')
           end if

           only_diag = .FALSE.
           write(std_out,'(a)',advance='no') ' Would you like to output only the diagonal [Y/N] ? '
           read(std_in,*)ans
           if (ans=='Y'.or.ans=='y') only_diag = .TRUE.

           ABI_ALLOCATE(epsm1_pole,(Hscr0%npwe,Hscr0%npwe,Hscr0%ncoeff,1))
           ABI_ALLOCATE(fval,(nfreqre*nfreqim))
           fval = CMPLX(0.0_gwp,0.0_gwp)

           do iqibz=1,Hscr0%nqibz
!            Get Q index for name
             call int2char4(iqibz,tagq)
             fname_dump=TRIM(fname)//'_POLE_EPS_ZPLANE_Q'//TRIM(tagq); unt_dump=get_unit()
             open(file=fname_dump,unit=unt_dump,status='replace',form='formatted',iostat=ios)
             qtmp(:)=Hscr0%qibz(:,iqibz); if (normv(qtmp,Cryst%gmet,'G')<GW_TOLQ0) qtmp(:)=zero
             call init_Gpairs_type(Gpairs_q,qtmp,Gsphere,Cryst)

             call read_pole_screening(fname,Hscr0%npwe,1,Hscr0%ncoeff,epsm1_pole,accesswff,&
&             spaceComm,iqiA=iqibz)
             
             do ig1=ig1_start,ig1_end
               do ig2=ig2_start,ig2_end

                 if (only_diag.AND.ig1/=ig2) CYCLE

                 call re_and_im_screening_with_phase(omega,fval,nfreqre*nfreqim, &
&                 epsm1_pole(ig1,ig2,:,1),Hscr0%ncoeff)
                 
                 if (ig1==ig2) fval(:) = fval(:) + CMPLX(1.0_gwp,0.0_gwp)

                 write(unt_dump,'(a,i4,2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,/,a,/)')&
&                 '# index= ',idx,'    ig1= ',ig1,'    ig2= ',ig2,&
&                 '# q = ',Hscr0%qibz(:,iqibz),&
&                 '# G = ',Hscr0%gvec(:,ig1),'  G''= ',Hscr0%gvec(:,ig2),&
&                 '#   omega [eV]           Re             Im '

                 do ios=1,nfreqim
                   do iomega=1,nfreqre
                     io = iomega + nfreqre*(ios-1)
                     write(unt_dump,'(2(f8.2),4x,2es16.8)') REAL(omega(io))*Ha_eV,&
&                     AIMAG(omega(io))*Ha_eV,fval(io)
                   end do
                   write(unt_dump,*)
                 end do
                 write(unt_dump,*)
                 write(unt_dump,*)
               end do ! ig2
             end do ! ig2
             close(unt_dump)

           end do ! iqibz
           
           ABI_DEALLOCATE(omega)
           ABI_DEALLOCATE(epsm1_pole)
           ABI_DEALLOCATE(fval)
           call destroy_Gpairs_type(Gpairs_q)
           call destroy_gsphere(Gsphere)

       end select ! Select case fform

     case(3) ! Extract dielectric function and plasmon-pole stuff --------------------------

       if (is_scr_p) then
         write(std_out,'(2(a))') ch10,' Option not available for pole-fits! Exiting ...'
         goto 100
       end if

       write(std_out,'(a)') ' 3 => Calculation of dielectric function and plasmon-pole model'

!      Check that a _SUS file was given
!      if (is_scr) then
!      MSG_ERROR(' You have to give a _SUS file as input for epsilon^(-1) and plasmon-pole model extraction')
!      end if

       npwe_asked=Hscr0%npwe; mqmem=Hscr0%nqibz
       call nullify_epsilonm1_results(Er)
       call init_Er_from_file(Er,fname,mqmem,npwe_asked,accesswff,spaceComm)

!      === Initialize the G-sphere ===
       call init_gsphere(Gsphere,.FALSE.,Cryst,Hscr0%npwe,gvec=Hscr0%gvec)
       call nullify_Gpairs_type(Gpairs_q)

       boxcutmin=two; igmax=Gsphere%shlim(Gsphere%nsh)
       ecut=Er%Hscr%Hdr%ecutdg

       call getng(boxcutmin,ecut,Gsphere%gmet,MPI_enreg_seq%me_fft,&
&       mgfft,nfft,ngfft,MPI_enreg_seq%nproc_fft,Cryst%nsym,MPI_enreg_seq%fft_option_lob,&
&       MPI_enreg_seq%paral_fft,Cryst%symrel)

!      I am using standard valued, it would be better to call indefo
!      ngfft(1:3)=Er%Hscr%Hdr%ngfft(1:3)
       ngfft(7)=112
       ngfft(8)=256     ! Was optimized for my PII 450MHz
       nfft = PRODUCT(ngfft(1:3))

       Dtset%icutcoul=3; Dtset%rcut=zero
       Dtset%vcutgeo=(/zero,zero,zero/); Dtset%boxcenter=(/zero,zero,zero/)

       if (Er%Hscr%nqlwl==0) then
         nqlwl=1
         ABI_ALLOCATE(qlwl,(3,nqlwl))
         qlwl(:,1)= GW_Q0_DEFAULT
       else
         nqlwl=Er%Hscr%nqlwl
         ABI_ALLOCATE(qlwl,(3,nqlwl))
         qlwl(:,:)=Er%Hscr%qlwl(:,1:nqlwl)
       end if

       call vcoul_init(Vcp,Gsphere,Qmesh,Kmesh,Dtset%rcut,Dtset%icutcoul,Dtset%vcutgeo,Hscr0%npwe,nqlwl,&
&       qlwl,Cryst%rprimd,ngfft,spaceComm)
       ABI_DEALLOCATE(qlwl)

!      === Get the density from an external file ===
!      * If meshes are not the same, do an FFT interpolation to have rhor on ngfft.

       call prompt(' Enter name for external DEN (or PAWDEN) file: ',fname_rho)

       ABI_ALLOCATE(rhor,(nfft,Hscr0%Hdr%nspden))
       if (Hscr0%Hdr%usepaw==1) then
         call get_rhor(fname_rho,accesswff,Hscr0%Hdr%nspden,nfft,ngfft,paral_kgb0, &
&         MPI_enreg_seq,rhor,get_pawden=.TRUE.)
       else
         call get_rhor(fname_rho,accesswff,Hscr0%Hdr%nspden,nfft,ngfft,paral_kgb0,MPI_enreg_seq,rhor)
       end if

       ABI_ALLOCATE(rhog,(2,nfft))
       call fourdp(1,rhog,rhor(:,1),-1,MPI_enreg_seq,nfft,ngfft,paral_kgb0,0)

       ABI_ALLOCATE(nhat,(nfft,Hscr0%Hdr%nspden*Hscr0%Hdr%usepaw))
       compch_sph=greatest_real; compch_fft=greatest_real
!      ABI_CHECK(Hscr0%Hdr%usepaw==0,'PAW not implemented')
       usexcnhat=0; usefinegrid=0

       nelectron_exp = hdr_get_nelect_byocc(Hscr0%Hdr)

       call test_charge(nfft,nelectron_exp,Hscr0%Hdr%nspden,rhor,Cryst%ucvol,&
&       Hscr0%Hdr%usepaw,usexcnhat,usefinegrid,compch_sph,compch_fft,drude_plsmf)
       GN_drude_plsmf = drude_plsmf

!      * Read and in case make Epsilon^{-1} according the the options specified
       id_required=4; ikxc=0; approx_type=0; option_test=0; dim_kxcg=0
       ABI_ALLOCATE(kxcg,(nfft,dim_kxcg))

       call prompt(' Enter prefix for output files: ',prefix)

!      TODO get rid of Dtfil
       fname_dump=TRIM(prefix)//'_SCR'
       Dtfil%filnam_ds(4)=prefix

       orig_npwe = Er%npwe
       write(std_out,'(2a,I0)') ch10,' Number of plane waves is: ',Er%npwe
       write(std_out,'(a)',advance='no') ' Would you like to change it [Y/N] ?'
       read(std_in,*) ans
       if (ans=='Y'.or.ans=='y') then
         write(std_out,'(a)',advance='no') ' Enter new no. of plane waves (0 means use old value): '
         read(std_in,*) ii
         if (ii>0.or.ii<=Er%npwe) Er%npwe = ii
         if (ii<0.or.ii>Er%npwe) then
           MSG_ERROR(' Wrong value for no. of plane waves!')
         end if
       end if

       if (is_scr) Er%mqmem=1
       if (is_sus) Er%mqmem=0
       call mkdump_Er(Er,Vcp,Er%npwe,Gsphere%gvec,dim_kxcg,kxcg,id_required,approx_type,ikxc,option_test,&
&       fname_dump,accesswff,nfft,ngfft,spaceComm)
       Er%mqmem=1

       call print_epsilonm1_results(Er)

       write(std_out,'(2a)',advance='no') ch10,&
&       ' Would you like to calculate the eigenvalues of eps^{-1}_GG''(omega) [Y/N] ? '
       read(std_in,*) ans
       if (ans=='Y'.or.ans=='y') then

         ABI_ALLOCATE(epsm1_eigen,(Er%npwe,Er%nomega))
         imax = 10
         if (Er%npwe < imax) imax = Er%npwe
         do iqibz=1,Er%nqibz
           call int2char4(iqibz,tagq)
           fname_eigen=TRIM(prefix)//'_EM1_EIG_Q'//TRIM(tagq)
           unt_dump=get_unit()
           open(file=fname_eigen,unit=unt_dump,status='replace',form='formatted',iostat=ios)
           call decompose_epsm1(Er,iqibz,epsm1_eigen)
           write(unt_dump,'(a)')       '# First (max 10) eigenvalues of eps^{-1}(omega)'
           write(unt_dump,'(a,3f12.6)')'# q = ',Hscr0%qibz(:,iqibz)
           write(unt_dump,'(a)')       '# REAL omega [eV]  REAL(eigen(esp^-1(1,w)))  AIMAG(eigen(esp^-1(1,w))  ...'
           do iomega=1,Er%nomega_r
             write(unt_dump,'(21(es16.8))')REAL(Er%omega(iomega))*Ha_eV,&
&             (REAL(epsm1_eigen(ii,iomega)),ii=1,imax),(AIMAG(epsm1_eigen(ii,iomega)),ii=1,imax)
           end do
           close(unt_dump)
         end do
         ABI_DEALLOCATE(epsm1_eigen)
         ABI_DEALLOCATE(kxcg)

       end if ! Calculate eigenvalues

!      === Analyze the PPmodel ===
       write(std_out,'(2a)') ch10,' Would you like to analyse plasmon-pole models [Y/N] ? '
       read(std_in,*)ans

       if (ans=='Y'.or.ans=='y') then

         write(std_out,'(2a,f6.2,a)') ch10,' Plasma frequency for GN PPM is: ',GN_drude_plsmf*Ha_eV, ' eV'
         write(std_out,'(a)',advance='no') ' Would you like to change it [Y/N] ?'
         read(std_in,*) ans
         if (ans=='Y'.or.ans=='y') then
           write(std_out,'(2a)',advance='no') ch10,' Enter plasma frequency [eV]: '
           read(std_in,*) GN_drude_plsmf
           GN_drude_plsmf = GN_drude_plsmf/Ha_eV
         end if

         write(std_out,'(2a)') ch10,' Would you like to calculate the plasmon-pole model'
         write(std_out,'(a)',advance='no')       '        eigenvalues of eps^{-1}_GG''(omega) [Y/N] ? '
         read(std_in,*) ans

         if (ans=='Y'.or.ans=='y') then

           ABI_ALLOCATE(ppm_eigen,(PPm%npwc,Er%nomega))
           imax = 10; if (Er%npwe < imax) imax = Er%npwe
           do iqibz=1,Er%nqibz

             do ppmodel=1,2

               call int2char4(iqibz,tagq)
               if (ppmodel==1) fname_dump=TRIM(prefix)//'_PPM_GN_EM1_EIG_Q'//TRIM(tagq)
               if (ppmodel==2) fname_dump=TRIM(prefix)//'_PPM_HL_EM1_EIG_Q'//TRIM(tagq)
               if (ppmodel==3) fname_dump=TRIM(prefix)//'_PPM_vdLH_EM1_EIG_Q'//TRIM(tagq)
               if (ppmodel==4) fname_dump=TRIM(prefix)//'_PPM_EF_EM1_EIG_Q'//TRIM(tagq)
               unt_dump=get_unit()
               open(file=fname_eigen,unit=unt_dump,status='new',form='formatted',iostat=ios)

               call ppm_free(PPm)
               if (ppmodel==1) then
                 call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,ppmodel,GN_drude_plsmf)
               else
                 call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,ppmodel,drude_plsmf)
               end if
               call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,nfft,Gsphere%gvec,ngfft,rhor(:,1),iqibz)

               call get_PPm_eigenvalues(PPm,iqibz,Er%Hscr%zcut,Er%nomega,Er%omega,Vcp,ppm_eigen)

               write(unt_dump,'(a)')       '# First (max 10) eigenvalues of eps^{-1}(omega) from Plasmon-pole model'
               write(unt_dump,'(a,3f12.6)')'# q = ',Hscr0%qibz(:,iqibz)
               select case(ppmodel)
                 case(1)
                   write(unt_dump,'(a)')     '# ppmodel = 1 : Godby - Needs'
                 case(2)
                   write(unt_dump,'(a)')     '# ppmodel = 2 : Hybertsen - Louie'
                 case(3)
                   write(unt_dump,'(a)')     '# ppmodel = 3 : von der Linden - Horsch'
                 case(4)
                   write(unt_dump,'(a)')     '# ppmodel = 4 : Engel - Farid'
               end select
               write(unt_dump,'(a)')       '# REAL omega [eV]  REAL(eigen(ppm_eps^-1(1,w)))  AIMAG(eigen(ppm_eps^-1(1,w))  ...'
               do iomega=1,Er%nomega_r
                 write(unt_dump,'(21(es16.8))')REAL(Er%omega(iomega))*Ha_eV,&
&                 (REAL(ppm_eigen(ii,iomega)),ii=1,imax),(AIMAG(ppm_eigen(ii,iomega)),ii=1,imax)
               end do
               close(unt_dump)

             end do !ppmodel

           end do ! iqibz
           ABI_DEALLOCATE(ppm_eigen)

         end if ! Calculate PPM eigenvalues

!        Optionally output eps^{-1}_GG''(w) for a given set of GG' and gridpoints
         write(std_out,'(2a)',advance='no') ch10,' Would you like to extract eps^{-1}_GG''(omega) for the PPM [Y/N] ?'
         read(std_in,*) ans

         if (ans=='Y'.or.ans=='y') then
!          Reconstruct e^{-1}_GG'(w) according to PPmodel for statistical analysis.
           write(std_out,'(a)') ' Enter the number of frequency points in the'
           write(std_out,'(a)') '  interval 0 - freqremax (0 means same as input file ): '
           read(std_in,*) nfreqre
           if (nfreqre==0) then
             nfreqre   = Er%nomega_r
             nfreqim   = Er%nomega_i
             nfreq_tot = Er%nomega
             freqremax = REAL(Er%omega(Er%nomega_r))
             ABI_ALLOCATE(omega,(nfreq_tot))
             omega(:) = Er%omega(:)
             same_freqs = .TRUE.
           else
             write(std_out,'(a)') ' Enter the value of freqremax (in eV): '
             read(std_in,*) freqremax
             nfreqim   = Er%nomega_i
             nfreq_tot = nfreqre+nfreqim
             ABI_ALLOCATE(omega,(nfreqre+Er%nomega_i))
             do iomega=1,nfreqre
               omega(iomega) =  CMPLX((freqremax/REAL((nfreqre-1)))*(iomega-1),zero)
             end do
             omega(nfreqre+1:nfreq_tot) = Er%omega(Er%nomega_r+1:Er%nomega)
             same_freqs = .FALSE.
           end if ! frequencies

           do iqibz=1,Er%nqibz

             qtmp(:)=Er%qibz(:,iqibz); if (normv(qtmp,Cryst%gmet,'G')<GW_TOLQ0) qtmp(:)=zero
             call init_Gpairs_type(Gpairs_q,qtmp,Gsphere,Cryst)
             call int2char4(iqibz,tagq)

             do ppmodel=1,2

               call ppm_free(PPm)
               if (ppmodel==1) then
                 call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,ppmodel,GN_drude_plsmf)
               else
                 call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,ppmodel,drude_plsmf)
               end if
               call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,&
&               nfft,Gsphere%gvec,ngfft,rhor(:,1),iqibz)

!              Prepare file for data on real omega axis
               if (ppmodel==1) fname_dump=TRIM(prefix)//'_PPM_w_GN_Q'//TRIM(tagq)
               if (ppmodel==2) fname_dump=TRIM(prefix)//'_PPM_w_HL_Q'//TRIM(tagq)
               if (ppmodel==3) fname_dump=TRIM(prefix)//'_PPM_w_vdLH_Q'//TRIM(tagq)
               if (ppmodel==4) fname_dump=TRIM(prefix)//'_PPM_w_EF_Q'//TRIM(tagq)
               unt_dump=get_unit()
               open(file=fname_dump,unit=unt_dump,status='replace',form='formatted',iostat=ios)
!              Prepare file for data on imaginary omega axis
               if (ppmodel==1) fname_dump2=TRIM(prefix)//'_PPM_iw_GN_Q'//TRIM(tagq)
               if (ppmodel==2) fname_dump2=TRIM(prefix)//'_PPM_iw_HL_Q'//TRIM(tagq)
               if (ppmodel==3) fname_dump2=TRIM(prefix)//'_PPM_iw_vdLH_Q'//TRIM(tagq)
               if (ppmodel==4) fname_dump2=TRIM(prefix)//'_PPM_iw_EF_Q'//TRIM(tagq)
               unt_dump2=get_unit()
               open(file=fname_dump2,unit=unt_dump2,status='replace',form='formatted',iostat=ios)

               ABI_ALLOCATE(em1_ppm,(nfreq_tot))

               ig1 = 0; ig2 = 0
               write(std_out,'(3a,I0,a,I0)') ch10,' Enter indices for G and G''.',&
&               'Entering 0 exits the loop. iqibz = ',iqibz,' ppmodel = ',ppmodel

               do
                 write(std_out,'(2(a),I0,a)',advance='NO') ch10,' Enter index for G (1 - ',Hscr0%npwe,' ): '
                 read(std_in,*)ig1
                 if (ig1==0) EXIT
                 if (ig1<0.OR.ig1>Er%npwe) MSG_ERROR(' index out of bounds')
                 write(std_out,'(2(a),I0,a)',advance='NO') ch10,' Enter index for G'' (1 - ',Hscr0%npwe,' ): '
                 read(std_in,*)ig2
                 if (ig2==0) EXIT
                 if (ig2<0.OR.ig2>Er%npwe) MSG_ERROR(' index out of bounds')

!                Generate the PPM representation of epsilon^-1
                 call getem1_from_PPm_one_ggp(PPm,iqibz,Er%Hscr%zcut,nfreq_tot,omega,Vcp,&
&                 em1_ppm,ig1,ig2)

                 write(unt_dump,'(a,I1)') '# epsilon^-1_GG''(omega) from ppmodel = ',ppmodel 
                 write(unt_dump,'(2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,&
&                 /,a,3F9.4,a,3F9.4,a,/a,f9.4,a,f9.4,a,/,a,/)')&
&                 '# ig1= ',ig1,'    ig2= ',ig2,&
&                 '# q = ',Er%qibz(:,iqibz),&
&                 '# G = ',Er%gvec(:,ig1),'  G''= ',Er%gvec(:,ig2),&
&                 '# G = (',MATMUL(two_pi*Cryst%gmet,Er%gvec(:,ig1)),&
&                 ')  G''= (',MATMUL(two_pi*Cryst%gmet,Er%gvec(:,ig2)),')',&
&                 '# 1/2|G|^2 =',half*normv(Er%gvec(:,ig1),Cryst%gmet,'G')**2,&
&                 ' Ha 1/2|G''|^2 =',half*normv(Er%gvec(:,ig2),Cryst%gmet,'G')**2,' Ha',&
&                 '#   omega [eV]           Re             Im '
                 write(unt_dump2,'(a,I1)') '# epsilon^-1_GG''(iomega) from ppmodel = ',ppmodel 
                 write(unt_dump2,'(2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,&
&                 /,a,3F9.4,a,3F9.4,a,/a,f9.4,a,f9.4,a,/,a,/)')&
&                 '# ig1= ',ig1,'    ig2= ',ig2,&
&                 '# q = ',Er%qibz(:,iqibz),&
&                 '# G = ',Er%gvec(:,ig1),'  G''= ',Er%gvec(:,ig2),&
&                 '# G = (',MATMUL(two_pi*Cryst%gmet,Er%gvec(:,ig1)),&
&                 ')  G''= (',MATMUL(two_pi*Cryst%gmet,Er%gvec(:,ig2)),')',&
&                 '# 1/2|G|^2 =',half*normv(Er%gvec(:,ig1),Cryst%gmet,'G')**2,&
&                 ' Ha 1/2|G''|^2 =',half*normv(Er%gvec(:,ig2),Cryst%gmet,'G')**2,' Ha',&
&                 '#   iomega [eV]           Re             Im '

                 do iomega=1,nfreqre
                   if (same_freqs) then
                     write(unt_dump,'(f8.2,4x,4es16.8)') REAL(omega(iomega))*Ha_eV,em1_ppm(iomega),&
&                     Er%epsm1(ig1,ig2,iomega,iqibz)
                   else
                     write(unt_dump,'(f8.2,4x,2es16.8)') REAL(omega(iomega))*Ha_eV,em1_ppm(iomega)
                   end if
                 end do
!                First output the iomega = 0 point
                 write(unt_dump2,'(f8.2,4x,4es16.8)') AIMAG(omega(1))*Ha_eV,em1_ppm(1),&
&                 Er%epsm1(ig1,ig2,1,iqibz)
!                Then the rest
                 do iomega=nfreqre+1,nfreq_tot
                   write(unt_dump2,'(f8.2,4x,4es16.8)') AIMAG(omega(iomega))*Ha_eV,em1_ppm(iomega),&
&                   Er%epsm1(ig1,ig2,iomega,iqibz)
                 end do
                 write(unt_dump,*)
                 write(unt_dump,*)
                 write(unt_dump2,*)
                 write(unt_dump2,*)
               end do ! Empty
               ABI_DEALLOCATE(em1_ppm)
               close(unt_dump); close(unt_dump2)

             end do ! ppmodel

           end do ! iqibz

           ABI_DEALLOCATE(omega)

         end if ! Output epsilon for PPM

!        Optionally statistics for all PPMs
         write(std_out,'(2a)',advance='no') ch10,' Would you like to output statistics for all PPMs [Y/N] ?'
         read(std_in,*) ans

         if (ans=='Y'.or.ans=='y') then

           nfreqre   = Er%nomega_r
           nfreq_tot = Er%nomega
           freqremax = REAL(Er%omega(Er%nomega_r))
           ABI_ALLOCATE(real_omega,(nfreqre))
           real_omega(:) = REAL(Er%omega(1:nfreqre))

           do iqibz=1,Er%nqibz

             do ppmodel=1,2

               qtmp(:)=Er%qibz(:,iqibz)
               if (normv(qtmp,Cryst%gmet,'G')<GW_TOLQ0) qtmp(:)=zero
               call init_Gpairs_type(Gpairs_q,qtmp,Gsphere,Cryst)

               call ppm_free(PPm)
               if (ppmodel==1) then
                 call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,ppmodel,GN_drude_plsmf)
               else
                 call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,ppmodel,drude_plsmf)
               end if
               call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,&
&               nfft,Gsphere%gvec,ngfft,rhor(:,1),iqibz)

!              Prepare ratios and density for the f-sum rule
               ABI_ALLOCATE(qratio,(orig_npwe,orig_npwe))
               istat = ABI_ALLOC_STAT
               ABI_CHECK(istat==0,"out-of-memory in qratio")
               ABI_ALLOCATE(rhoggp,(Er%npwe,Er%npwe))
               istat = ABI_ALLOC_STAT
               ABI_CHECK(istat==0,"out-of-memory in rhoggp")
               call cqratio(orig_npwe,Gsphere%gvec,qtmp,Cryst%gmet,Cryst%gprimd,qratio)
!              Arrange n(G-G')->n(G,G')
               ierr=0
               do ig1=1,Er%npwe
                 do ig2=1,Er%npwe
                   gmgp_idx = g2ifft(Gsphere%gvec(:,ig1)-Gsphere%gvec(:,ig2),ngfft)
                   if (gmgp_idx/=0) then
                     rhoggp(ig1,ig2)=CMPLX(rhog(1,gmgp_idx),rhog(2,gmgp_idx))
                   else
                     ierr=ierr+1
                     rhoggp(ig1,ig2)=czero
                   end if
                 end do
               end do
               if (ierr/=0) then
                 write(std_out,'(a,i0,a)')&
&                 ' Found ',ierr,' G1-G2 vectors falling outside the FFT box. '
               end if

!              Prepare files
               call int2char4(iqibz,tagq)
               if (ppmodel==1) fname_dump=TRIM(prefix)//'_norms_GN_Q'//TRIM(tagq)
               if (ppmodel==2) fname_dump=TRIM(prefix)//'_norms_HL_Q'//TRIM(tagq)
               unt_dump=get_unit()
               open(file=fname_dump,unit=unt_dump,status='replace',form='formatted',iostat=ios)
               write(unt_dump,'(a)') '# Various norms integrated through spline interpolation'
               write(unt_dump,'(a)') '# over all frequencies in the input file,'
               write(unt_dump,'(a)') '# for all G and G'' vectors.'
               write(unt_dump,'(a,I0)') '#               ppmodel: ',ppmodel
               write(unt_dump,'(a,I0)') '# Number of frequencies: ',nfreqre
               write(unt_dump,'(a,f12.6)') '# Maximum frequency    : ',freqremax
               write(unt_dump,'(a)') '# Columns:'
               write(unt_dump,'(2a)') '#  ig1      ig2   |eps-eps_PPM|/|eps|',&
&               '   |eps-eps_PPM|    |eps|     |eps_PPM|            G                  G'''
               if (ppmodel==1) fname_dump2=TRIM(prefix)//'_f_sumrule_GN_Q'//TRIM(tagq)
               if (ppmodel==2) fname_dump2=TRIM(prefix)//'_f_sumrule_HL_Q'//TRIM(tagq)
               unt_dump2=get_unit()
               open(file=fname_dump2,unit=unt_dump2,status='replace',form='formatted',iostat=ios)
               write(unt_dump2,'(a)') '# The fulfillment of the f-sum rule: I(epsilon) ='
               write(unt_dump2,'(a)') '#   int_0^{inf}{omega*Im[epsilon_G,G''(omega)]}/C_qGG'''
               write(unt_dump2,'(a)') '# C_qGG'' = '
               write(unt_dump2,'(a)') '#   -Pi/2*omega_p^2*(q+G)*(q+G'')/|q+G|^2*n(G-G'')/n(0)'
               write(unt_dump2,'(a)') '# for all G and G'' vectors.'
               write(unt_dump2,'(a,I0)') '#               ppmodel: ',ppmodel
               write(unt_dump2,'(a,I0)') '# Number of frequencies: ',nfreqre
               write(unt_dump2,'(a,f12.6)') '# Maximum frequency    : ',freqremax
               write(unt_dump2,'(a)') '# Columns:'
               write(unt_dump2,'(3a)') '#  ig1      ig2   I(epsilon)',&
&               '   I(eps_PPM)   Re[n(G-G'')]    Im[n(G-G'')]    qratio      I1*C_qGG''',&
&               ' Re[Omegatwsq] Im[Omegatwsq]   Re[omegatw]   Im[omegatw]    |G|    1/2|G|^2'

               ABI_ALLOCATE(em1_ppm,(nfreq_tot))
               ABI_ALLOCATE(ftab,(nfreqre))
               ABI_ALLOCATE(ysp,(3,nfreqre))
               ABI_ALLOCATE(work,(nfreqre))
               ABI_ALLOCATE(eint,(nfreqre))
               do ig1=1,Er%npwe
                 write(std_out,'(2(a,I0))') ' ig1= ',ig1, ' of ',Er%npwe
                 do ig2=1,Er%npwe
!                  ig2 = ig1

                   call getem1_from_PPm_one_ggp(PPm,iqibz,Er%Hscr%zcut,nfreq_tot,Er%omega,Vcp,&
&                   em1_ppm,ig1,ig2)

!                  Calculate norms in real
                   eps_diff=0; eps_norm=0; eps_ppm_norm=0
                   ftab(1:nfreqre) = ABS(Er%epsm1(ig1,ig2,1:nfreqre,iqibz)-em1_ppm(1:nfreqre))
                   call cspint(ftab,real_omega,nfreqre,real_omega(1),&
&                   real_omega(nfreqre),ysp,eint,work,eps_diff)
                   ftab(1:nfreqre) = ABS(Er%epsm1(ig1,ig2,1:nfreqre,iqibz))
                   call cspint(ftab,real_omega,nfreqre,real_omega(1),&
&                   real_omega(nfreqre),ysp,eint,work,eps_norm)
                   ftab(1:nfreqre) = ABS(em1_ppm(1:nfreqre))
                   call cspint(ftab,real_omega,nfreqre,real_omega(1),&
&                   real_omega(nfreqre),ysp,eint,work,eps_ppm_norm)
                   write(unt_dump,'(2i6,f12.4,3es14.4,6i4)') ig1,ig2,eps_diff/eps_norm,eps_diff,&
&                   eps_norm,eps_ppm_norm,Er%gvec(:,ig1),Er%gvec(:,ig2)

!                  Evaluate the f-sum rule
                   if (ig1==ig2) then 
                     ftab(1:nfreqre) = real_omega(1:nfreqre)*AIMAG(Er%epsm1(ig1,ig2,1:nfreqre,iqibz))
                   else ! Dephase first - HERE epsm1 is changed!
                     call remove_phase(epsm1(ig1,ig2,:,1),Hscr_file(1)%nomega,phase)
                     ftab(1:nfreqre) = real_omega(1:nfreqre)*AIMAG(Er%epsm1(ig1,ig2,1:nfreqre,iqibz))
                   end if

                   call cspint(ftab,real_omega,nfreqre,real_omega(1),&
&                   real_omega(nfreqre),ysp,eint,work,eps_diff)
                   
                   if (ig1==ig2) then
                     factor = -two*pi*pi*REAL(rhoggp(ig1,ig2))*qratio(ig1,ig2)
                   else
                     rhoggp(ig1,ig2) = CMPLX(COS(phase),-SIN(phase))*rhoggp(ig1,ig2)
                     factor = -two*pi*pi*REAL(rhoggp(ig1,ig2))*qratio(ig1,ig2)
                   end if

                   if (ABS(qratio(ig1,ig2))>zero) then
                     value1 = eps_diff/factor
                     if (ppmodel==1) then
                       value2 = -pi*half*(REAL(PPm%bigomegatwsq(iqibz)%value(ig1,ig2))&
&                       /(REAL(PPm%omegatw(iqibz)%value(ig1,ig2))))&
&                       /factor*(2*sqrt(pi*rhoggp(1,1)))
                     else
                       value2 = -pi*half*(SQRT(REAL(PPm%bigomegatwsq(iqibz)%value(ig1,ig2))))&
&                       /factor*(2*sqrt(pi*rhoggp(1,1)))
                     end if
                   else
                     value1 = zero
                     value2 = zero
                   end if

                   write(unt_dump2,'(2i6,12es14.4)') ig1,ig2,value1,value2,&
&                   REAL(rhoggp(ig1,ig2)),AIMAG(rhoggp(ig1,ig2)),qratio(ig1,ig2),&
&                   eps_diff,REAL(PPm%bigomegatwsq(iqibz)%value(ig1,ig2)),&
&                   AIMAG(PPm%bigomegatwsq(iqibz)%value(ig1,ig2)),&
&                   REAL(PPm%omegatw(iqibz)%value(ig1,ig2)),&
&                   AIMAG(PPm%omegatw(iqibz)%value(ig1,ig2)),&
&                   normv(Er%gvec(:,ig1),Cryst%gmet,'G'),&
&                   half*normv(Er%gvec(:,ig1),Cryst%gmet,'G')**2

                 end do !ig2
               end do !ig1

               ABI_DEALLOCATE(em1_ppm)
               ABI_DEALLOCATE(ftab)
               ABI_DEALLOCATE(ysp)
               ABI_DEALLOCATE(work)
               ABI_DEALLOCATE(eint)
               ABI_DEALLOCATE(qratio)
               ABI_DEALLOCATE(rhoggp)
               close(unt_dump); close(unt_dump2)

             end do ! ppmodel

           end do ! iqibz

           ABI_DEALLOCATE(real_omega)

!          write(std_out,'(2a)') ch10,' Would you like to check the fulfillment of the f-sum rule ?'
!          read(std_in,*) ans
!          if (ans=='Y'.or.ans=='y') then

!          end if ! Check f-sum rule

         end if ! Output statistics

         call ppm_free(PPm)

       end if ! If ppmodel>0

       ABI_DEALLOCATE(rhor)
       ABI_DEALLOCATE(rhog)
       ABI_DEALLOCATE(nhat)

       call vcoul_free(Vcp)
       call destroy_Epsilonm1_results(Er)
       call destroy_Gpairs_type(Gpairs_q)
       call destroy_gsphere(Gsphere)

     case(4) ! Remove frequencies ----------------------------------------------------------

       if (is_scr_p) then
         write(std_out,'(2(a))') ch10,' Option not available for pole-fits! Exiting ...'
         goto 100
       end if

!      Calculate the total number of real freq
       nfreqre = 0
       nfreqim = 0
       do ifrq=1,Hscr_file(1)%nomega
!        If frequency is not imaginary, count.
         if (AIMAG(Hscr_file(1)%omega(ifrq))<tol8) nfreqre = nfreqre + 1
         if (REAL(Hscr_file(1)%omega(ifrq))<tol8.and.AIMAG(Hscr_file(1)%omega(ifrq))>tol8)  nfreqim = nfreqim + 1
       end do ! ifrq

!      Test for no real frequencies
       if (nfreqre==0) then
         write(std_out,'(2(a))') ch10,' No real frequencies in file! Exiting ...'
         goto 100
       end if

       nfreq_tot = nfreqre + nfreqim ! Here nfreq_tot becomes the *true* number of freq
       write(std_out,'(2a,I0,a)') ch10,' Found ',nfreq_tot,' frequencies.'
       write(std_out,'(2(a,I0),2a)') ' ',nfreqre,' real, and ',nfreqim,' imaginary.',ch10


       ABI_ALLOCATE(omega_storage,(nfreq_tot))
       ABI_ALLOCATE(freq_indx,(nfreq_tot,1))
       omega_storage = CMPLX(-one,-one); freq_indx = 0;

       write(std_out,'(2(a))') ch10,' Do you want to remove every other real frequency  (= 1) ?'
       write(std_out,'(a)')         '  or specify for each real frequency individually  (= 2) ?'
       write(std_out,'(a)')         '  or remove ALL real frequencies                   (= 3) ?'
       read(std_in,*)choice

       select case(choice)

         case(1) ! Remove every other frequency

           write(std_out,'(2(a))') ch10,' Removing every other real frequency, i.e. every even one.'
           write(std_out,'(a)')         ' If the total number of frequencies is odd, the first and last one will be kept.'
           write(std_out,'(a)')         ' If the total number is even, the first one will still be in the final set.'

           ii=nfreqre
           nfreqre = 0
           do ifrq=1,ii
             if (.not.iseven(ifrq)) then
               nfreqre = nfreqre + 1
               omega_storage(nfreqre) = Hscr_file(1)%omega(ifrq)
               freq_indx(nfreqre,1) = ifrq
             end if
           end do ! ifrq
           write(std_out,'(2a,I0,a)') ch10,' ',nfreqre,' real frequencies will be kept.'

         case(2) ! Specify freq. individually

           ii=nfreqre
           nfreqre = 0
           do ifrq=1,ii
             write(std_out,'(a,f12.6,a)') ' Would you like to keep freq. at: ',REAL(Hscr_file(1)%omega(ifrq))*Ha_eV,' eV? [y/n]'
             read(std_in,*) ans
             if (ans=='Y'.or.ans=='y') then
               nfreqre = nfreqre + 1
               omega_storage(nfreqre) = Hscr_file(1)%omega(ifrq)
               freq_indx(nfreqre,1) = ifrq
             end if
           end do ! ifrq
           write(std_out,'(2a,I0,a)') ch10,' ',nfreqre,' real frequencies will be kept.'

         case(3) ! Remove all real freq.

           ii=nfreqre
           nfreqre = 0

           case default ! Bail if choice is wrong
           write(std_out,*) ' Invalid choice! Exiting...'
           goto 100
       end select

!      Add imaginary frequencies if any
       if (nfreqim>0) then
         nfreqim = 0
         do ifrq=1,Hscr_file(1)%nomega
           if (AIMAG(Hscr_file(1)%omega(ifrq)) > tol8) then
             nfreqim = nfreqim + 1
             omega_storage(nfreqre+nfreqim) = Hscr_file(1)%omega(ifrq)
             freq_indx(nfreqre+nfreqim,1) = ifrq
           end if
         end do
       end if

       nfreq_tot = nfreqre + nfreqim ! Here nfreq_tot becomes the *true* number of freq
       write(std_out,'(2a,I0,a)') ch10,' Finally, we have ',nfreq_tot,' frequencies.'
       write(std_out,'(2(a,I0),2a)') ' ',nfreqre,' real, and ',nfreqim,' imaginary.',ch10

!      First copy the old header
       if (is_sus) then ! Hack to avoid bug for SUS files
         ii = Hscr_file(1)%headform
         Hscr_file(1)%headform = 56
       end if
       call copy_ScrHdr(Hscr_file(1),Hscr_merge)
       if (is_sus) then ! Hack to avoid bug for SUS files
         Hscr_file(1)%headform = ii
         Hscr_merge%headform = ii
       end if

!      Then modify entries for new frequency grid
       Hscr_merge%nomega = nfreq_tot
       ABI_DEALLOCATE(Hscr_merge%omega)
       ABI_ALLOCATE(Hscr_merge%omega,(nfreq_tot))
       Hscr_merge%omega(1:nfreq_tot) = omega_storage(1:nfreq_tot)
!      Check if q=0 is part of the set and reallocate and set lwing and uwing
       do iqibz=1,Hscr_merge%nqibz
         if (Hscr_merge%qibz(1,iqibz)<tol6.AND.&
&         Hscr_merge%qibz(2,iqibz)<tol6.AND.&
&         Hscr_merge%qibz(3,iqibz)<tol6) then
           indx_q0 = iqibz
         end if
         npwe4mI = Hscr_merge%npwe*Hscr_merge%nI
         npwe4mJ = Hscr_merge%npwe*Hscr_merge%nJ
         ABI_DEALLOCATE(Hscr_merge%uwing)
         ABI_DEALLOCATE(Hscr_merge%lwing)
         ABI_ALLOCATE(Hscr_merge%uwing,(npwe4mJ,nfreq_tot,Hscr_merge%nqlwl))
         ABI_ALLOCATE(Hscr_merge%lwing,(npwe4mI,nfreq_tot,Hscr_merge%nqlwl))
         do ifrq=1,nfreq_tot
           Hscr_merge%uwing(:,ifrq,:)=&
&           Hscr_file(1)%uwing(:,freq_indx(ifrq,1),:)
           Hscr_merge%lwing(:,ifrq,:)=&
&           Hscr_file(1)%lwing(:,freq_indx(ifrq,1),:)
         end do
       end do

       call prompt(' Enter the full name of the final output file: ',fname_out)

!      Print new header for info
       call print_ScrHdr(Hscr_merge,header='Header of the final file',unit=std_out,prtvol=1)

       unt_out=get_unit()
       open(unit=unt_out,file=fname_out,status='new',form='unformatted',iostat=ios)
       ABI_CHECK(ios==0,"Opening "//TRIM(fname_out))

!      * Write the header.
!      TODO Use new format but first of all fix problem with hdr_check
       rdwr=2
       if (is_scr) fform_merge=1002
       if (is_sus) fform_merge=1102
       call scr_hdr_io(fform_merge,rdwr,unt_out,spaceComm,master,accesswff,Hscr_merge)

       npwe4mI = Hscr_merge%npwe*Hscr_merge%nI
       npwe4mJ = Hscr_merge%npwe*Hscr_merge%nJ
       nomega4m = Hscr_merge%nomega
       write(std_out,'(a,f16.1,a)')' Memory required for new matrix   = ',&
&       2.0*gwpc*npwe4mI*npwe4mJ*nomega4m*b2Mb," [Mb]."
       ABI_ALLOCATE(epsm1,(npwe4mI,npwe4mJ,nomega4m,1))
       istat = ABI_ALLOC_STAT
       ABI_CHECK(istat==0,'out of memory in epsm1')

       do iqibz=1,Hscr_merge%nqibz
         ifile=1

!        allocate temporary array
         npwe4mI = Hscr_file(ifile)%npwe*Hscr_file(ifile)%nI
         npwe4mJ = Hscr_file(ifile)%npwe*Hscr_file(ifile)%nJ
         nomega4m = Hscr_file(ifile)%nomega
         write(std_out,'(a,f16.1,a)')' Memory required for temporary matrix= ',&
         2.0*gwpc*npwe4mI*npwe4mJ*nomega4m*b2Mb,' [Mb].'
         ABI_ALLOCATE(epsm1_temp,(npwe4mI,npwe4mJ,nomega4m,1))
         istat = ABI_ALLOC_STAT
         ABI_CHECK(istat==0,'out of memory in epsm1_temp')

!        read screening
         fname=filenames(ifile)
         call read_screening(fname,npwe4mI,1,nomega4m,epsm1_temp,accesswff,spaceComm,iqiA=iqibz)

!        Copy matrices for relevant frequencies
         do ifrq=1,nfreqre
           epsm1(:,:,ifrq,1)=epsm1_temp(:,:,freq_indx(ifrq,ifile),1)
         end do

!        Add imaginary part if needed
         if (nfreqim>0) then
           do ifrq=nfreqre+1,Hscr_merge%nomega
             epsm1(:,:,ifrq,1)=epsm1_temp(:,:,freq_indx(ifrq,ifile),1)
           end do
         end if

         ABI_DEALLOCATE(epsm1_temp)

         npwe4mI = Hscr_merge%npwe*Hscr_merge%nI
         nomega4m = Hscr_merge%nomega
         call write_screening(unt_out,accesswff,npwe4mI,nomega4m,epsm1)

       end do ! iqibz

       ABI_DEALLOCATE(epsm1)
       close(unt_out)
       write(msg,'(3a)')ch10,' ==== Real frequencies have been removed successfully === ',ch10
       call wrtout(std_out,msg,'COLL')

       ABI_DEALLOCATE(omega_storage)
       ABI_DEALLOCATE(freq_indx)

     case(5) ! Remove imaginary frequencies ------------------------------------------------

       if (is_scr_p) then
         write(std_out,'(2(a))') ch10,' Option not available for pole-fits! Exiting ...'
         goto 100
       end if

!      Calculate the total number of freq
       nfreqre = 0
       nfreqim = 0
       do ifrq=1,Hscr_file(1)%nomega
!        If frequency is not imaginary, count.
         if (AIMAG(Hscr_file(1)%omega(ifrq))<tol8) nfreqre = nfreqre + 1
         if (REAL(Hscr_file(1)%omega(ifrq))<tol8.and.AIMAG(Hscr_file(1)%omega(ifrq))>tol8)  nfreqim = nfreqim + 1
       end do ! ifrq

!      Test for no real frequencies
       if (nfreqim==0) then
         write(std_out,'(2(a))') ch10,' No imaginary frequencies in file! Exiting ...'
         goto 100
       end if

       nfreq_tot = nfreqre + nfreqim ! Here nfreq_tot becomes the *true* number of freq
       write(std_out,'(2a,I0,a)') ch10,' Found ',nfreq_tot,' frequencies.'
       write(std_out,'(2(a,I0),2a)') ' ',nfreqre,' real, and ',nfreqim,' imaginary.',ch10


       ABI_ALLOCATE(omega_storage,(nfreq_tot))
       ABI_ALLOCATE(freq_indx,(nfreq_tot,1))
       omega_storage = CMPLX(-one,-one); freq_indx = 0;

!      Specify freq. individually
       ii=nfreq_tot
       nfreqim = 0
       do ifrq=nfreqre+1,ii
         write(std_out,'(a,f12.6,a)')&
&         ' Would you like to keep imaginary freq. at: ',AIMAG(Hscr_file(1)%omega(ifrq))*Ha_eV,' eV? [y/n]'
         read(std_in,*) ans
         if (ans=='Y'.or.ans=='y') then
           nfreqim = nfreqim + 1
           omega_storage(nfreqre+nfreqim) = Hscr_file(1)%omega(ifrq)
           freq_indx(nfreqre+nfreqim,1) = ifrq
         end if
       end do ! ifrq
       write(std_out,'(2a,I0,a)') ch10,' ',nfreqim,' imaginary frequencies will be kept.'

!      Add real frequencies if any
       if (nfreqre>0) then
         nfreqre = 0
         do ifrq=1,Hscr_file(1)%nomega
           if (AIMAG(Hscr_file(1)%omega(ifrq)) < tol8) then
             nfreqre = nfreqre + 1
             omega_storage(nfreqre) = Hscr_file(1)%omega(ifrq)
             freq_indx(nfreqre,1) = ifrq
           end if
         end do
       end if

       nfreq_tot = nfreqre + nfreqim ! Here nfreq_tot becomes the *true* number of freq
       write(std_out,'(2a,I0,a)') ch10,' Finally, we have ',nfreq_tot,' frequencies.'
       write(std_out,'(2(a,I0),2a)') ' ',nfreqre,' real, and ',nfreqim,' imaginary.',ch10

!      First copy the old header
       if (is_sus) then ! Hack to avoid bug for SUS files
         ii = Hscr_file(1)%headform
         Hscr_file(1)%headform = 56
       end if
       call copy_ScrHdr(Hscr_file(1),Hscr_merge)
       if (is_sus) then ! Hack to avoid bug for SUS files
         Hscr_file(1)%headform = ii
         Hscr_merge%headform = ii
       end if

!      Then modify entries for new frequency grid
       Hscr_merge%nomega = nfreq_tot
       ABI_DEALLOCATE(Hscr_merge%omega)
       ABI_ALLOCATE(Hscr_merge%omega,(nfreq_tot))
       Hscr_merge%omega(1:nfreq_tot) = omega_storage(1:nfreq_tot)
!      Check if q=0 is part of the set and reallocate and set lwing and uwing
       do iqibz=1,Hscr_merge%nqibz
         if (Hscr_merge%qibz(1,iqibz)<tol6.AND.&
&         Hscr_merge%qibz(2,iqibz)<tol6.AND.&
&         Hscr_merge%qibz(3,iqibz)<tol6) then
           indx_q0 = iqibz
         end if
         npwe4mI = Hscr_merge%npwe*Hscr_merge%nI
         npwe4mJ = Hscr_merge%npwe*Hscr_merge%nJ
         ABI_DEALLOCATE(Hscr_merge%uwing)
         ABI_DEALLOCATE(Hscr_merge%lwing)
         ABI_ALLOCATE(Hscr_merge%uwing,(npwe4mJ,nfreq_tot,Hscr_merge%nqlwl))
         ABI_ALLOCATE(Hscr_merge%lwing,(npwe4mI,nfreq_tot,Hscr_merge%nqlwl))
         do ifrq=1,nfreq_tot
           Hscr_merge%uwing(:,ifrq,:)=&
&           Hscr_file(1)%uwing(:,freq_indx(ifrq,1),:)
           Hscr_merge%lwing(:,ifrq,:)=&
&           Hscr_file(1)%lwing(:,freq_indx(ifrq,1),:)
         end do
       end do

       call prompt(' Enter the full name of the final output file: ',fname_out)

!      Print new header for info
       call print_ScrHdr(Hscr_merge,header='Header of the final file',unit=std_out,prtvol=1)

       unt_out=get_unit()
       open(unit=unt_out,file=fname_out,status='new',form='unformatted',iostat=ios)
       if (ios/=0) then
         write(msg,'(3a)')' Opening file ',TRIM(fname_out),' as new-unformatted'
         MSG_ERROR(msg)
       end if

!      * Write the header.
!      TODO Use new format but first of all fix problem with hdr_check
       rdwr=2
       if (is_scr) fform_merge=1002
       if (is_sus) fform_merge=1102
       call scr_hdr_io(fform_merge,rdwr,unt_out,spaceComm,master,accesswff,Hscr_merge)

       npwe4mI = Hscr_merge%npwe*Hscr_merge%nI
       npwe4mJ = Hscr_merge%npwe*Hscr_merge%nJ
       nomega4m = Hscr_merge%nomega
       write(std_out,'(a,f16.1,a)')' Memory required for new matrix   = ',&
&       2.0*gwpc*npwe4mI*npwe4mJ*nomega4m*b2Mb," [Mb]."
       ABI_ALLOCATE(epsm1,(npwe4mI,npwe4mJ,nomega4m,1))
       istat = ABI_ALLOC_STAT
       ABI_CHECK(istat==0,'out of memory in epsm1')

       do iqibz=1,Hscr_merge%nqibz
         ifile=1

!        allocate temporary array
         npwe4mI = Hscr_file(ifile)%npwe*Hscr_file(ifile)%nI
         npwe4mJ = Hscr_file(ifile)%npwe*Hscr_file(ifile)%nJ
         nomega4m = Hscr_file(ifile)%nomega
         write(std_out,'(a,f16.1,a)')' Memory required for temporary matrix= ',&
         2.0*gwpc*npwe4mI*npwe4mJ*nomega4m*b2Mb,' [Mb].'
         ABI_ALLOCATE(epsm1_temp,(npwe4mI,npwe4mJ,nomega4m,1))
         istat = ABI_ALLOC_STAT
         ABI_CHECK(istat==0,'out of memory in epsm1_temp')

!        read screening
         fname=filenames(ifile)
         call read_screening(fname,npwe4mI,1,nomega4m,epsm1_temp,accesswff,spaceComm,iqiA=iqibz)

!        Copy matrices for relevant frequencies
         do ifrq=1,nfreqre
           epsm1(:,:,ifrq,1)=epsm1_temp(:,:,freq_indx(ifrq,ifile),1)
         end do

!        Add imaginary part if needed
         if (nfreqim>0) then
           do ifrq=nfreqre+1,Hscr_merge%nomega
             epsm1(:,:,ifrq,1)=epsm1_temp(:,:,freq_indx(ifrq,ifile),1)
           end do
         end if

         ABI_DEALLOCATE(epsm1_temp)

         npwe4mI = Hscr_merge%npwe*Hscr_merge%nI
         nomega4m = Hscr_merge%nomega
         call write_screening(unt_out,accesswff,npwe4mI,nomega4m,epsm1)

       end do ! iqibz

       ABI_DEALLOCATE(epsm1)
       close(unt_out)
       write(msg,'(3a)')ch10,' ==== Imaginary frequencies have been removed successfully === ',ch10
       call wrtout(std_out,msg,'COLL')

       ABI_DEALLOCATE(omega_storage)
       ABI_DEALLOCATE(freq_indx)



     case(6) ! Model screening -------------------------------------------------------------

       if (is_scr_p) then
         write(std_out,'(2(a))') ch10,' This is an already pole-fit file! Exiting ...'
         goto 100
       end if

#ifdef HAVE_ALGO_LEVMAR

       write(std_out,'(a)') ' 6 => Fitting of a generalised Drude-Lorentz model'

!      === Initialize the G-sphere ===
       call init_gsphere(Gsphere,.FALSE.,Cryst,Hscr0%npwe,gvec=Hscr0%gvec)

       call nullify_Gpairs_type(Gpairs_q)

       write(std_out,'(2(a))') ch10,' Do you want to fit just one q,G,G'' matrix element?  (= 1) ?'
       write(std_out,'(a)')         '                             or ALL matrix elements  (= 2) ?'
       read(std_in,*)choice

       select case(choice)

         case(1) ! Fit just one element

           ABI_ALLOCATE(epsm1,(Hscr_file(1)%npwe,Hscr_file(1)%npwe,Hscr_file(1)%nomega,1))
           istat = ABI_ALLOC_STAT
           if (istat/=0) STOP 'out of memory in epsm1'

           ig1 = 0; ig2 = 0; iqibz = 0;
           write(std_out,'(2(a),I0,a)',advance='NO') ch10,'  Enter the index for q (1 - ',Hscr_file(1)%nqibz,' ): '
           read(std_in,*)iqibz
           if (iqibz<1.OR.iqibz>Hscr_file(1)%nqibz) MSG_ERROR(' q index out of bounds')
           write(std_out,'(a,I0,a)',advance='NO')         '  Enter the index for G (1 - ',Hscr_file(1)%npwe,' ): '
           read(std_in,*)ig1
           if (ig1<1.OR.ig1>Hscr_file(1)%npwe) MSG_ERROR(' G index out of bounds')
           write(std_out,'(a,I0,a)',advance='NO')         ' Enter the index for G'' (1 - ',Hscr_file(1)%npwe,' ): '
           read(std_in,*)ig2
           if (ig2<1.OR.ig2>Hscr_file(1)%npwe) MSG_ERROR(' G'' index out of bounds')

           write(std_out,'(2a,I0,a)',advance='NO')   ch10,' Enter the number of poles to fit: '
           read(std_in,*)npoles
           ABI_CHECK(npoles>0,' Number of poles should be positive and at least 1')

           fname=filenames(1)
           call read_screening(fname,Hscr_file(1)%npwe,1,Hscr_file(1)%nomega,epsm1,accesswff,spaceComm,iqiA=iqibz)

!          Calculate the total number of freq
           nfreqre = 0
           nfreqim = 0
           do ifrq=1,Hscr_file(1)%nomega
!            If frequency is not imaginary, count.
             if (AIMAG(Hscr_file(1)%omega(ifrq))<tol8) nfreqre = nfreqre + 1
             if (REAL(Hscr_file(1)%omega(ifrq))<tol8.and.AIMAG(Hscr_file(1)%omega(ifrq))>tol8)  nfreqim = nfreqim + 1
           end do ! ifrq

!          Test for grid
           if (nfreqre+nfreqim==Hscr_file(1)%nomega) then
             write(std_out,'(2(a))') ch10,' This file does not contain a grid in the complex plane! Exiting ...'
             goto 100
           end if
           
           nfreqc = Hscr_file(1)%nomega - nfreqre - nfreqim
           write(std_out,'(2a,I0,a)') ch10,' Found ',Hscr0%nomega,' frequencies.'
           write(std_out,'(2(a,I0),2a)') ' ',nfreqre,' real, and ',nfreqim,' imaginary.',ch10
           write(std_out,'(a,I0)') ' There is a grid in the complex plane with ',nfreqc
           write(std_out,'(2a)')   '  extra frequencies in the list.',ch10

!          Allocate
           ABI_ALLOCATE(coeff,(npoles*3+1))
           ABI_ALLOCATE(startcoeff,(npoles*3+1))
           ABI_ALLOCATE(refval,(Hscr_file(1)%nomega))
           ABI_ALLOCATE(imfval,(Hscr_file(1)%nomega))
           ABI_ALLOCATE(imfval2,(Hscr_file(1)%nomega))
           ABI_ALLOCATE(refval2,(Hscr_file(1)%nomega))
           coeff = 0.0_gwp; startcoeff = 0.0_gwp

!          Remove phase factor if ig1/=ig2
           if (ig1/=ig2) then
             call remove_phase(epsm1(ig1,ig2,:,1),Hscr_file(1)%nomega,phase)
             coeff(npoles*3+1)=phase; startcoeff(npoles*3+1)=phase
             write(std_out,'(a,f6.3,a)') ' ig1/=ig2, Removed phase factor: pi*(',phase/pi,')'
           end if

           imfval(:)   = AIMAG(epsm1(ig1,ig2,:,1))
           refval(:)   = REAL(epsm1(ig1,ig2,:,1))
           if (ig1==ig2.and.is_scr) refval = refval - one

           call sequential_fitting(Hscr_file(1)%omega(:),refval,imfval,&
&           Hscr_file(1)%nomega,nfreqre,coeff(1:npoles*3),npoles*3,11,startcoeff(1:npoles*3))

!          Calculate fit
           refval = zero; refval2 = zero; imfval = zero; imfval2 = zero
           call im_screening(Hscr_file(1)%omega,imfval,Hscr_file(1)%nomega,coeff(1:npoles*3),npoles*3)
           call re_screening(Hscr_file(1)%omega,refval,Hscr_file(1)%nomega,coeff(1:npoles*3),npoles*3)
           call im_screening(Hscr_file(1)%omega,imfval2,Hscr_file(1)%nomega,startcoeff(1:npoles*3),npoles*3)
           call re_screening(Hscr_file(1)%omega,refval2,Hscr_file(1)%nomega,startcoeff(1:npoles*3),npoles*3)

!          Restore phase factors
           if (ig1/=ig2) then
             write(std_out,'(a,f6.3,a)') ' ig1/=ig2, Restoring phase factor: pi*(',coeff(npoles*3+1)/pi,')'
             rerot = COS(coeff(npoles*3+1))
             imrot = SIN(coeff(npoles*3+1))
             do iomega=1,Hscr_file(1)%nomega
               retemp = REAL(epsm1(ig1,ig2,iomega,1))
               imtemp = AIMAG(epsm1(ig1,ig2,iomega,1))
               epsm1(ig1,ig2,iomega,1) = CMPLX(rerot*retemp-imrot*imtemp,&
&               rerot*imtemp+imrot*retemp)
               retemp = refval(iomega)
               imtemp = imfval(iomega)
               refval(iomega) =  rerot*retemp - imrot*imtemp
               imfval(iomega) =  rerot*imtemp + imrot*retemp
               retemp = refval2(iomega)
               imtemp = imfval2(iomega)
               refval2(iomega) = rerot*retemp - imrot*imtemp
               imfval2(iomega) = rerot*imtemp + imrot*retemp
             end do
           end if

           if (ig1==ig2.and.is_scr) refval = refval + one
           if (ig1==ig2.and.is_scr) refval2 = refval2 + one

           fname_dump='epsilon_model_fit_real_omega.dat'; unt_dump=get_unit()
           open(file=TRIM(fname_dump),unit=unt_dump,status='replace',form='formatted',iostat=ios)
           write(unt_dump,'(a,i4,2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,/,a,i0,/,a,100ES16.8)')&
&           '# index= ',idx,'    ig1= ',ig1,'    ig2= ',ig2,&
&           '# q = ',Hscr_file(1)%qibz(:,iqibz),&
&           '# G = ',Hscr_file(1)%gvec(:,ig1),'  G''= ',Hscr_file(1)%gvec(:,ig2),&
&           '# npoles = ',npoles, &
&           '# coeffs = ',coeff
           write(unt_dump,'(2a)')'#   omega (eV)   Re[Model]   Im[Model]    Re[epsilon]',&
&           '    Im[epsilon]   Re[initmodel]   Im[initmodel]'
           do iomega=1,nfreqre
             write(unt_dump,'(f8.2,4x,6ES16.8)') REAL(Hscr_file(1)%omega(iomega))*Ha_eV, &
&             refval(iomega),imfval(iomega), &
&             REAL(epsm1(ig1,ig2,iomega,1)),AIMAG(epsm1(ig1,ig2,iomega,1)), &
&             refval2(iomega),imfval2(iomega)
           end do
           close(unt_dump)

           fname_dump='epsilon_model_fit_imag_omega.dat'; unt_dump=get_unit()
           open(file=TRIM(fname_dump),unit=unt_dump,status='replace',form='formatted',iostat=ios)
           write(unt_dump,'(a,i4,2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,/,a,i0,/,a,100ES16.8)')&
&           '# index= ',idx,'    ig1= ',ig1,'    ig2= ',ig2,&
&           '# q = ',Hscr_file(1)%qibz(:,iqibz),&
&           '# G = ',Hscr_file(1)%gvec(:,ig1),'  G''= ',Hscr_file(1)%gvec(:,ig2),&
&           '# npoles = ',npoles, &
&           '# coeffs = ',coeff
           write(unt_dump,'(2a)')'#   iomega (eV)   Re[Model]   Im[Model]    Re[epsilon]',&
&           '    Im[epsilon]   Re[initmodel]   Im[initmodel]'
           write(unt_dump,'(f8.2,4x,6ES16.8)') AIMAG(Hscr_file(1)%omega(1))*Ha_eV, &
&           refval(1),imfval(1), &
&           REAL(epsm1(ig1,ig2,1,1)),AIMAG(epsm1(ig1,ig2,1,1)), &
&           refval2(1),imfval2(1)
           do iomega=nfreqre+1,nfreqre+nfreqim
             write(unt_dump,'(f8.2,4x,6ES16.8)') AIMAG(Hscr_file(1)%omega(iomega))*Ha_eV, &
&             refval(iomega),imfval(iomega), &
&             REAL(epsm1(ig1,ig2,iomega,1)),AIMAG(epsm1(ig1,ig2,iomega,1)), &
&             refval2(iomega),imfval2(iomega)
           end do
           close(unt_dump)

           fname_dump='epsilon_model_fit_zplane.dat'; unt_dump=get_unit()
           open(file=fname_dump,unit=unt_dump,status='replace',form='formatted',iostat=ios)
           write(unt_dump,'(a,i4,2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,/,a,/)')&
&           '# index= ',idx,'    ig1= ',ig1,'    ig2= ',ig2,&
&           '# q = ',Hscr0%qibz(:,iqibz),&
&           '# G = ',Hscr0%gvec(:,ig1),'  G''= ',Hscr0%gvec(:,ig2),&
&           '#   omega [eV]           Re             Im '
           do iomega=1,nfreqre
             write(unt_dump,'(2(f8.2),4x,6es16.8)') REAL(Hscr0%omega(iomega))*Ha_eV,&
&             AIMAG(Hscr0%omega(iomega))*Ha_eV,refval(iomega),imfval(iomega),epsm1(ig1,ig2,iomega,1),&
&             refval2(iomega),imfval2(iomega)
           end do
           write(unt_dump,*)
           do ios=1,nfreqim
             do iomega=1,nfreqre
               if (iomega==1) then
                 io = nfreqre + ios
               else
                 io = nfreqre + nfreqim + (ios-1)*(nfreqre-1) + (iomega-1)
               end if
               write(unt_dump,'(2(f8.2),4x,6es16.8)') REAL(Hscr0%omega(io))*Ha_eV,&
&               AIMAG(Hscr0%omega(io))*Ha_eV,refval(io),imfval(io),epsm1(ig1,ig2,io,1),&
&               refval2(io),imfval2(io)
             end do
             write(unt_dump,*)
           end do
           close(unt_dump) 

         case(2) ! Fit everything

           if (is_scr/=.true.) then
             write(std_out,'(2(a))') ch10,&
&             ' Complete fits can only be done for _SCR file for now. Exiting ...'
             goto 100
           end if

           call prompt(' Enter the full name of the final output file: ',fname_out)

           ABI_ALLOCATE(epsm1,(Hscr_file(1)%npwe,Hscr_file(1)%npwe,Hscr_file(1)%nomega,1))
           istat = ABI_ALLOC_STAT

           write(std_out,'(2a,I0,a)',advance='NO')   ch10,' Enter the number of poles to fit: '
           read(std_in,*)npoles
           ABI_CHECK(npoles>0,' Number of poles should be positive and at least 1')

!          Calculate the total number of freq
           nfreqre = 0
           nfreqim = 0
           do ifrq=1,Hscr_file(1)%nomega
!            If frequency is not imaginary, count.
             if (AIMAG(Hscr_file(1)%omega(ifrq))<tol8) nfreqre = nfreqre + 1
             if (REAL(Hscr_file(1)%omega(ifrq))<tol8.and.AIMAG(Hscr_file(1)%omega(ifrq))>tol8)&
&             nfreqim = nfreqim + 1
           end do ! ifrq
!          Test for grid
           if (nfreqre+nfreqim==Hscr_file(1)%nomega) then
             write(std_out,'(2(a))') ch10,&
&             ' This file does not contain a grid in the complex plane! Exiting ...'
             goto 100
           end if

!          Set up header
           call copy_ScrHdr(Hscr_file(1),Hscr_merge)
!          Then modify entries. It is no longer frequencies that are stored
!          for a pole-fitted file, but rather the pole coefficients (fform=2002)
           Hscr_merge%fform  = 2002
           Hscr_merge%npoles = npoles
           Hscr_merge%ncoeff = npoles*3+1
!          Check if q=0 is part of the set and reallocate and set lwing and uwing
           do iqibz=1,Hscr_merge%nqibz
             if (Hscr_merge%qibz(1,iqibz)<tol6.AND.&
&             Hscr_merge%qibz(2,iqibz)<tol6.AND.&
&             Hscr_merge%qibz(3,iqibz)<tol6) then
               indx_q0 = iqibz
             end if
             npwe4mI = Hscr_merge%npwe*Hscr_merge%nI
             npwe4mJ = Hscr_merge%npwe*Hscr_merge%nJ
             ABI_DEALLOCATE(Hscr_merge%uwing)
             ABI_DEALLOCATE(Hscr_merge%lwing)
             ABI_ALLOCATE(Hscr_merge%uwing,(npwe4mJ,nfreq_tot,Hscr_merge%nqlwl))
             ABI_ALLOCATE(Hscr_merge%lwing,(npwe4mI,nfreq_tot,Hscr_merge%nqlwl))
             Hscr_merge%uwing = czero
             Hscr_merge%lwing = czero
           end do

!          Print new header for info
           write(std_out,*) ''
           call print_ScrHdr(Hscr_merge,header='Header of the pole-fit file',unit=std_out,prtvol=1)

!          Write header and screening
           unt_out=get_unit()
           open(unit=unt_out,file=fname_out,status='replace',form='unformatted',iostat=ios)
           if (ios/=0) then
             write(msg,'(3a)')' Opening file ',TRIM(fname_out),' as replace-unformatted'
             MSG_ERROR(msg)
           end if
           rdwr=2; fform_merge=2002
           call scr_hdr_io(fform_merge,rdwr,unt_out,spaceComm,master,accesswff,Hscr_merge)

           do iqibz=1,Hscr_file(1)%nqibz

             ABI_ALLOCATE(epsm1_pole,(Hscr_file(1)%npwe,Hscr_file(1)%npwe,npoles*3+1,1))

             fname=filenames(1)
             call read_screening(fname,Hscr_file(1)%npwe,1,Hscr_file(1)%nomega,epsm1,accesswff,&
&             spaceComm,iqiA=iqibz)

#ifdef HAVE_GW_OPENMP
!!OMP *** OPENMP SECTION *** Added by MS
!            $OMP PARALLEL SHARED(Hscr_file,epsm1,nfreqre,epsm1_pole, &
!            $OMP                 is_scr,npoles,std_out,iqibz,one) &
!            $OMP          PRIVATE(ig1,ig2,refval,imfval,phase)
!!OMP write(std_out,'(a,i0)') ' Entering openmp loop. Number of threads: ',omp_get_num_threads()
#endif
             ABI_ALLOCATE(refval,(Hscr_file(1)%nomega))
             ABI_ALLOCATE(imfval,(Hscr_file(1)%nomega))
             epsm1_pole = 0.0_gwp; refval = 0.0_gwp; imfval = 0.0_gwp
#ifdef HAVE_GW_OPENMP
!            $OMP DO 
#endif
             do ig1=1,Hscr_file(1)%npwe
               do ig2=1,Hscr_file(1)%npwe
                 
!                Remove phase factor if ig1/=ig2
                 if (ig1/=ig2) then
                   call remove_phase(epsm1(ig1,ig2,:,1),Hscr_file(1)%nomega,phase)
                   epsm1_pole(ig1,ig2,npoles*3+1,1)=phase
                 end if

                 imfval(:)   = AIMAG(epsm1(ig1,ig2,:,1))
                 refval(:)   = REAL(epsm1(ig1,ig2,:,1))
                 if (ig1==ig2.and.is_scr) refval = refval - one

                 call sequential_fitting(Hscr_file(1)%omega(:),refval,imfval,&
&                 Hscr_file(1)%nomega,nfreqre,epsm1_pole(ig1,ig2,1:npoles*3,1),&
&                 npoles*3,0)

               end do ! ig2
               if (mod(ig1,10)==1)write(std_out,'(2(a,i0))') &
&               ' Done ig1 = ',ig1,' of ',Hscr_file(1)%npwe
             end do   ! ig1
             ABI_DEALLOCATE(refval)
             ABI_DEALLOCATE(imfval)

#ifdef HAVE_GW_OPENMP
!            $OMP END DO
!            $OMP END PARALLEL
#endif
             write(std_out,'(a,i0)') ' Done iqibz = ',iqibz 

!            Write the screening
             npwe4mI = Hscr_merge%npwe*Hscr_merge%nI
             call write_pole_screening(unt_out,accesswff,npwe4mI,Hscr_merge%ncoeff,epsm1_pole)

             ABI_DEALLOCATE(epsm1_pole)

           end do !iqibz

!          Clean up
           ABI_DEALLOCATE(epsm1)
           close(unt_out)
           write(msg,'(5(a))')ch10,' ==== Pole-fit saved to ',TRIM(fname_out),' === ',ch10
           call wrtout(std_out,msg,'COLL')

       end select ! What to fit for model func

#else
!      There is no levmar library in the linking
       write(std_out,'(2(a))') ch10,' MRGSCR has not been linked with the levmar library! Exiting ...'
       goto 100
#endif

     case (7) ! Interpolation for real frequency -------------------------------------
       MSG_ERROR("Still under development")
!      Perform b-spline interpolation onto a user defined real frequency grid

!      Calculate the total number of real freq
       nfreqre = 0
       nfreqim = 0
       do ifrq=1,Hscr_file(1)%nomega
!        If frequency is not imaginary, count.
         if (AIMAG(Hscr_file(1)%omega(ifrq))<tol8) nfreqre = nfreqre + 1
         if (REAL(Hscr_file(1)%omega(ifrq))<tol8.and.AIMAG(Hscr_file(1)%omega(ifrq))>tol8)  nfreqim = nfreqim + 1
       end do ! ifrq

!      Test for no real frequencies or grid file
       if (nfreqre==0) then
         write(std_out,'(2(a))') ch10,' ERROR - No real frequencies in file! Exiting ...'
         goto 100
       end if
       if ((nfreqre+nfreqim)/=Hscr_file(1)%nomega) then
         write(std_out,'(2(a))') ch10,' ERROR - This is a file with a full grid! Exiting ...'
         goto 100
       end if

       nfreq_tot = nfreqre + nfreqim ! Here nfreq_tot becomes the *true* number of freq
       write(std_out,'(2a,I0,a)') ch10,' Found ',nfreq_tot,' frequencies.'
       write(std_out,'(2(a,I0),2a)') ' ',nfreqre,' real, and ',nfreqim,' imaginary.',ch10


       write(std_out,'(2(a))') ch10,' Do you want to use an equidistant grid  (= 1) ?'
       write(std_out,'(a)')         '  or a tangent transform grid            (= 2) ?'
       read(std_in,*)choice

       select case(choice)
         
         case(1) ! Equidistant grid

           write(std_out,'(2a)') ch10,' Enter new number of gridpoints for the equidistant grid:'
           read(std_in,*)new_nfreqre
           write(std_out,'(2a)') ch10,' Enter freqremin [in eV] for the equidistant grid:'
           read(std_in,*)freqremin
           freqremin = freqremin/Ha_eV
           write(std_out,'(2a)') ch10,' Enter freqremax [in eV] for the equidistant grid:'
           read(std_in,*)freqremax
           freqremax = freqremax/Ha_eV

           if (new_nfreqre<1) then
             write(std_out,'(2(a))') ch10,' ERROR - One or more points need to be selected! Exiting ...'
             goto 100
           end if
           if (freqremin<freqremax) then
             write(std_out,'(2(a))') ch10,' ERROR - freqremin > freqremax! Exiting ...'
             goto 100
           end if
           if (freqremin<zero.OR.freqremax<zero) then
             write(std_out,'(2(a))') ch10,' ERROR - freqremax or freqremin negative! Exiting ...'
             goto 100
           end if

!          allocate and evaluate new frequency grid
           ABI_ALLOCATE(omega_new,(new_nfreqre+nfreqim))
           omega_new = CMPLX(-one,-one)
           domegareal=(freqremax-freqremin)/(new_nfreqre-1)
           do ifrq=2,new_nfreqre
             omega_new(ifrq) = CMPLX(freqremin+(ifrq-1)*domegareal,zero,kind=dpc)
           end do
           omega_new((new_nfreqre+1):(new_nfreqre+nfreqim))=Hscr_file(1)%omega((nfreqre+1):nfreq_tot)

         case(2) ! Tangent transform grid
           write(std_out,'(2(a))') ch10,' Tangent grid not implemented yet. Exiting ...'
           goto 100

           case default
           write(std_out,'(2(a))') ch10,' Wrong choice of frequency grid! Exiting ...'
           goto 100

       end select ! Type of grid

       



     case (8) ! Interpolation in q-space --------------------------------------------

       write(std_out,'(a)') ' 8 => Interpolating the matrix with linear interpolation'
       MSG_ERROR("Still under development")

       call prompt(' Enter the name of the output file: ',intp_fname)

       call prompt(' Enter the number of k-points: ',new_nkibz)
       ABI_ALLOCATE(new_kibz,(3,new_nkibz))
       call prompt(' Enter the list of new k-points: ',new_kibz)

       call prompt(' Enter old_qptrlatt: ',old_qptrlatt)
       call prompt(' Enter old_nshiftq: ',old_nshiftq)
       ABI_ALLOCATE(old_shiftq,(3,old_nshiftq))
       call prompt(' Enter old_shiftq: ',old_shiftq)

       call interpolate_w(filenames(1),intp_fname,old_qptrlatt,old_nshiftq,old_shiftq,new_nkibz,new_kibz,spacecomm)

       ABI_DEALLOCATE(new_kibz)
       ABI_DEALLOCATE(old_shiftq)


       case default ! Bail if choice is wrong

       write(std_out,*) ' Invalid choice! Exiting...'
       goto 100

   end select

 end if ! Single file mode

 100 continue

!=====================
!==== Free memory ====
!=====================
 call destroy_crystal(Cryst)
 call destroy_BZ_mesh_type(Kmesh)
 call destroy_BZ_mesh_type(Qmesh)

 nullify(Hscr0)
 if (nfiles>1) then
   call free_scrhdr(Hscr_merge)
 end if
 do ifile=1,nfiles
   call free_scrhdr(Hscr_file(ifile))
 end do
 ABI_DEALLOCATE(Hscr_file)
 ABI_DEALLOCATE(filenames)

 call destroy_mpi_enreg(mpi_enreg)
 call xmpi_end()

 end program mrgscr
!!***
