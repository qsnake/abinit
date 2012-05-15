!{\src2tex{textfont=tt}}
!!****f* ABINIT/outwf
!! NAME
!! outwf
!!
!! FUNCTION
!! Conduct output of a "wave-functions" file.
!!  - Compute the maximal residual
!!  - Then open a permanent file wff2 for final output of wf data
!!  - Create a new header for the file.
!!  - Write wave-functions (and energies)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, AR, MB, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg(2,mcg)=wavefunction array (storage if nkpt>1)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen( (2*mband)**response *mband*nkpt*nsppol)=
!!                  eigenvalues (hartree) for all bands at each k point
!!  filnam= character string giving the root to form the name of the
!!   output WFK or WFQ file if response==0, otherwise it is the filename.
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kptns(3,nkpt)=k points in terms of recip primitive translations
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mkmem=maximum number of k-points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum number of plane waves
!!  mxfh=last dimension of the xfhist array
!!  natom=number of atoms in unit cell
!!  nband=number of bands
!!  nkpt=number of k points
!!  npwarr(nkpt)=number of planewaves in basis and on boundary for each k
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nstep=desired number of electron iteration steps
!!  nxfh=actual number of (x,f) history pairs, see xfhist array.
!!  occ(mband*nkpt*nsppol)=occupations for all bands at each k point
!!  resid(mband*nkpt*nsppol)=squared residuals for each band and k point
!!   where resid(n,k)=|<C(n,k)|(H-e(n,k))|C(n,k)>|^2 for the ground state
!!  response: if == 0, GS wavefunctions , if == 1, RF wavefunctions
!!  unwff2=unit for output of wavefunction
!!  wffnow=structure information for current wavefunction (if nkpt>1)
!!  xfhist(3,natom+4,2,mxfh)=(x,f) history array,
!!                                 also includes rprim and stress
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! * The name of the file wff2 might be the same as that of the file wff1.
!! * The routine includes closing wffnow.
!!
!! PARENTS
!!      berryphase_new,gstate,loper3
!!
!! CHILDREN
!!      handle_ncerr,hdr_io,hdr_io_etsf,hdr_io_netcdf,hdr_skip,ini_wf_netcdf
!!      leave_new,leave_test,rwwf,timab,wffclose,wffdelete,wffkg,wffoffset
!!      wffopen,wrtout,wvl_write,xbarrier_mpi,xcomm_init,xcomm_self,xcomm_world
!!      xdefineoff,xexch_mpi,xmaster_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine outwf(cg,dtset,eigen,filnam,hdr,kg,kptns,mband,mcg,mkmem,&
 &                mpi_enreg,mpw,mxfh,natom,nband,nkpt,npwarr,&
 &                nsppol,nstep,nxfh,occ,resid,response,unwff2,&
 &                wffnow,wfs,wvl,xfhist)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_wffile
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outwf'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
 use interfaces_61_ionetcdf
 use interfaces_62_wvl_wfs
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif
!Arguments ------------------------------------
 integer, intent(in) :: mband,mcg,mkmem,mpw,mxfh,natom,nkpt,nsppol
 integer, intent(in) :: nstep,nxfh,response,unwff2
 character(len=fnlen), intent(in) :: filnam
 type(MPI_type), intent(inout) :: mpi_enreg
 type(dataset_type), intent(in) :: dtset
 type(hdr_type), intent(inout) :: hdr
 type(wffile_type), intent(inout) :: wffnow
 type(wvl_wf_type), intent(in) :: wfs
 type(wvl_internal_type), intent(in) :: wvl
 integer, intent(in) :: kg(3,mpw*mkmem),nband(nkpt*nsppol),npwarr(nkpt)
 real(dp), intent(inout) :: cg(2,mcg)
 real(dp), intent(in) :: eigen((2*mband)**response*mband*nkpt*nsppol),kptns(3,nkpt)
 real(dp), intent(in) :: occ(mband*nkpt*nsppol),resid(mband*nkpt*nsppol)
 real(dp), intent(in) :: xfhist(3,natom+4,2,mxfh)

!Local variables-------------------------------
 integer,parameter :: nkpt_max=50
 integer :: accesswff,action,band_index,fform,formeig,headform,iband,ibdkpt,icg
  integer :: icg0,ierr,ii,ikg,ikpt,isppol,ixfh,master,mcg_disk,me,me0,my_nspinor
 integer :: nband_disk,nband_k,nkpt_eff,nmaster,npw_k,option,optkg,rdwr,sender,source
 integer :: spaceComm,spaceComm_io,spacecomsender,spaceWorld,sread,sskip,tim_rwwf
 integer :: xfdim2
 real(dp) :: residk,residm,resims
 logical :: mydata,swrite,tmaster
 character(len=500) :: message
 type(wffile_type) :: wff2
 integer,allocatable :: kg_disk(:,:)
 real(dp) :: tsec(2)
!There are two cases of use of the cg_disk : if mkmem==0, or parallel treatment.
 real(dp),allocatable :: cg_disk(:,:),eig_dum(:),eig_k(:),occ_dum(:),occ_k(:)
!no_abirules
#if defined HAVE_MPI
           !Variables introduced for MPI version
           integer :: ipwnbd
#endif

#if defined HAVE_TRIO_NETCDF
           !netCDF variables
           integer :: ncid_hdr,ncerr
           integer :: nxfh_id, mxfh_id, xfdim2_id, dim2inout_id, dimr3_id,xfhist_id
           integer :: nxfh_tmp,mxfh_tmp,xfdim2_tmp,dim2inout_tmp

#endif

! *************************************************************************
!For readability of the source file, define a "me" variable
!also in the sequential case

!DEBUG
!write(std_out,*)' outwf : enter'
!write(std_out,*)' outwf : trim(filnam)=',trim(filnam)
!ENDDEBUG
 xfdim2 = natom+4

!Init me
 call xme_init(mpi_enreg,me)
 me0=me
!Define master
 call xmaster_init(mpi_enreg,master)
!Init mpi_comm
 call xcomm_world(mpi_enreg,spaceWorld)
 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart_4d)
 call xcomm_self(spaceComm_io)
 if (mpi_enreg%mode_para=='b' ) spaceComm_io= mpi_enreg%commcart_3d
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)
 tim_rwwf =0
 source = master
 sread = master
 tmaster=(master==me)
 swrite=tmaster
 sender=-1

!Compute mean square and maximum residual over all bands and k points and spins
!(disregard k point weights and occupation numbers here)
 band_index=sum(nband(1:nkpt*nsppol))
 resims=sum(resid(1:band_index))/dble(band_index)

!Find largest residual over bands, k points, and spins, except for nbdbuf highest bands
!Already AVAILABLE in hdr ?!
 ibdkpt=1
 residm=zero
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     nband_k=max(1,nband_k-dtset%nbdbuf)
     residm=max(residm,maxval(resid(ibdkpt:ibdkpt+nband_k-1)))
     ibdkpt=ibdkpt+nband_k
   end do
 end do

 write(message, '(a,1p,e12.4,a,e12.4)' ) &
& ' Mean square residual over all n,k,spin= ',resims,'; max=',residm
 call wrtout(ab_out,message,'COLL')

 band_index=0
 nkpt_eff=nkpt
 if( (dtset%prtvol==0 .or. dtset%prtvol==1) .and. nkpt_eff>nkpt_max ) nkpt_eff=nkpt_max

!Loop over spin again
 do isppol=1,nsppol
!  Give (squared) residuals for all bands at each k
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
!    Will not print all residuals when prtvol=0 or 1
     if(ikpt<=nkpt_eff)then
!      Find largest residual over all bands for given k point
       residk=maxval(resid(1+band_index:nband_k+band_index))
       write(message, '(1x,3f8.4,3x,i2,1p,e13.5,a)' ) &
&       kptns(1:3,ikpt),isppol,residk,' kpt; spin; max resid(k); each band:'
       call wrtout(ab_out,message,'COLL')
       do ii=0,(nband_k-1)/8
         write(message, '(1x,1p,8e9.2)' ) &
&         (resid(iband+band_index),iband=1+ii*8,min(nband_k,8+ii*8))
         call wrtout(ab_out,message,'COLL')
       end do
     else if(ikpt==nkpt_eff+1)then
       write(message, '(a,a)' ) &
&       ' outwf : prtvol=0 or 1, do not print more k-points.',ch10
       call wrtout(ab_out,message,'COLL')
     end if
     band_index=band_index+nband_k
   end do
 end do

!Will write the wavefunction file only when nstep>0
 if (nstep>0 .and. dtset%prtwf/=0) then

!  Only the master write the file, except if MPI I/O, but the
!  full wff dataset should be provided to WffOpen in this case
   accesswff=-1
   if(dtset%accesswff==1) then
     accesswff=1
#if defined HAVE_TRIO_NETCDF
   else if(dtset%accesswff==2) then
     accesswff=2
!    Create empty netCDF file
     ncerr = nf90_create(path=filnam, cmode=NF90_CLOBBER, ncid=ncid_hdr)
     call handle_ncerr(ncerr," create netcdf wavefunction file")
     ncerr = nf90_close(ncid_hdr)
     call handle_ncerr(ncerr," close netcdf wavefunction file")
#endif
#if defined HAVE_TRIO_ETSF_IO
   else if (dtset%accesswff == 3) then
     accesswff = 3
#endif
   end if
!  DEBUG
!  write(std_out,*) 'outwf : accesswff = ', accesswff
!  write(std_out,*) 'outwf : wffnow%accesswff = ', wffnow%accesswff,wffnow%kgwff
!  ENDDEBUG
   call WffOpen(accesswff,spaceComm,filnam,ierr,wff2,master,me0,unwff2,spaceComm_io)
!  Conduct wavefunction output to wff2
   write(message, '(a,a)' )&
&   ' outwf  : write wavefunction to file ',trim(filnam)
   call wrtout(std_out,message,'COLL')

   ABI_ALLOCATE(kg_disk,(3,mpw))

   mcg_disk=mpw*my_nspinor*mband
   formeig=0 ; if(response==1)formeig=1
   if(mkmem == 0 )then
     ABI_ALLOCATE(eig_dum,( (2*mband)**formeig * mband))
     ABI_ALLOCATE(occ_dum,(mband))
   end if
   ABI_ALLOCATE(eig_k,( (2*mband)**formeig * mband))
   ABI_ALLOCATE(occ_k,(mband))

#if defined HAVE_MPI
   call leave_test()
!  Compute mband and mpw
   if(mkmem/=0) ABI_ALLOCATE(cg_disk,(2,mcg_disk))
#endif

   if (mkmem==0) then

!    Skip wffnow header
     call hdr_skip(wffnow,ierr)

     ABI_ALLOCATE(cg_disk,(2,mcg_disk))
!    Define offsets, in case of MPI I/O
     call WffKg(wffnow,1)
     call xdefineOff(formeig,wffnow,mpi_enreg,nband,npwarr,dtset%nspinor,nsppol,nkpt)
   end if  !mkmem = 0


   band_index=0
   icg=0
   if((mpi_enreg%paralbd==0) .or. (mpi_enreg%paralbd>1))tim_rwwf=6
   if(mpi_enreg%paralbd==1)tim_rwwf=12

!  Write header info for new wf file
   rdwr=2
   if (dtset%usewvl == 0) then
     fform=2
   else
!    Use 200 as radical for naming file format
!    used by wavelets.
     fform = 200
   end if

   if (wff2%accesswff < 2) then
     call hdr_io(fform,hdr,rdwr,wff2)
     call WffKg(wff2,1)
#if defined HAVE_TRIO_NETCDF
   else if (wff2%accesswff == 2.and.tmaster) then

!    DEBUG
!    write(std_out,*) 'outwf : entering hdr_io_netcdf'
!    ENDDEBUG

     call hdr_io_netcdf(fform,hdr,rdwr,wff2)

     call ini_wf_netcdf(mpw,wff2%unwff,response)
#endif
#if defined HAVE_TRIO_ETSF_IO
   else if (wff2%accesswff == 3 .and. tmaster) then

!    DEBUG
!    write(std_out,*) 'outwf : entering hdr_io_etsf'
!    ENDDEBUG

     call hdr_io_etsf(fform, hdr, rdwr, wff2%unwff)
#endif
   end if

   do isppol=1,nsppol

     ikg=0

     do ikpt=1,nkpt
       nband_k=nband(ikpt+(isppol-1)*nkpt)
       npw_k=npwarr(ikpt)

!      Read the wavefunction block, without the eigenvalues
       if(mkmem==0)then
#if defined HAVE_MPI
         sread=-1
         if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))==0) sread=me
#endif
         if(sread==me)then

           headform=0 ; icg0=0 ; option=-2 ; optkg=1
           call rwwf(cg_disk,eig_dum,formeig,headform,&
&           icg0,ikpt,isppol,kg_disk,mband,mcg_disk,mpi_enreg,nband_k,nband_disk,&
&           npw_k,my_nspinor,occ_dum,option,optkg,tim_rwwf,wffnow)

           if(nband_k/=nband_disk)then
             write(message, '(a,a,a,a,i4,a,i6,a,a,a,i6,a)' ) ch10,&
&             ' outwf : BUG -',ch10,&
&             '  For k pt number',ikpt,' disk file has',nband_disk,' bands',ch10,&
&             '  but input file gave nband=',nband_k,'.'
             call wrtout(std_out,message,'PERS')
             call leave_new('PERS')
           end if ! nband check

         end if ! sread==me
       end if ! mkmem

#if defined HAVE_MPI
       if (dtset%usewvl == 0) then
         call xbarrier_mpi(spaceWorld)

!        Must transfer the wavefunctions to the master processor
!        Separate sections for paralbd=1 or other values ; might be merged
         if(mpi_enreg%paralbd==0 .or. mpi_enreg%paralbd>1)then
           nmaster=0
           source=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol))
           mydata=.false.
           if(source==me)mydata=.true.
           action=0
!          I am the master node, and I have the data in cg or cg_disk
           if((tmaster).and.(mydata))action=1
!          I am not the master, and I have the data => send to master
           if((.not.tmaster).and.(mydata))action=2
!          I am the master, and I receive the data
           if((tmaster).and.(.not.mydata))action=3

!          I have the data in cg or cg_disk ( MPI_IO case)
           if (accesswff==1  ) then
             action = 0
             sender=-1
             swrite=.false.
             if (mydata)then
               action=1
               swrite=.true.
               sender=me
             end if
           end if

!          I am the master node, and I have the data in cg or cg_disk
!          I have the data in cg or cg_disk ( MPI_IO case)
           if(action==1)then
             if(mkmem/=0)then
!              Copy from kg to kg_disk
               kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
!              Copy from cg to cg_disk
               do ipwnbd=1,nband_k*npw_k*my_nspinor
                 cg_disk(1,ipwnbd)=cg(1,ipwnbd+icg)
                 cg_disk(2,ipwnbd)=cg(2,ipwnbd+icg)
               end do
             end if
           end if


!          I am not the master, and I have the data => send to master
!          I am the master, and I receive the data
           if ( action==2.or.action==3) then
             call timab(48,1,tsec)
             if(mkmem/=0 .and. action==2)then
               call xexch_mpi(kg(:,1+ikg:npw_k+ikg),3*npw_k,source,kg_disk,nmaster,spaceWorld,ierr)
               call xexch_mpi(cg(:,icg+1:icg+nband_k*npw_k*my_nspinor),2*nband_k*npw_k*my_nspinor, &
&               source,cg_disk,nmaster,spaceWorld,ierr)
             else
               call xexch_mpi(kg_disk,3*npw_k,source,kg_disk,nmaster,spaceWorld,ierr)
               call xexch_mpi(cg_disk,2*nband_k*npw_k*my_nspinor,source,cg_disk,nmaster, &
&               spaceWorld,ierr)
             end if
             call timab(48,2,tsec)
           end if


         else if(mpi_enreg%paralbd==1)then
           nmaster=0
#if defined HAVE_MPI_IO
           sender=-1
           if( accesswff ==1 ) then
             nmaster=mpi_enreg%proc_distrb(ikpt,1,isppol)
             sender=nmaster
           end if
#endif

!          Note the loop over bands
           do iband=1,nband_k

!            The message passing related to kg is counted as one band
             action=0

!            I am the master node, and I have the data in cg or cg_disk
             if( mpi_enreg%proc_distrb(ikpt,iband,isppol)==nmaster .and. &
&             me==nmaster) then
               action=1
!              I am not the master, and I have the data => send to master
             elseif( mpi_enreg%proc_distrb(ikpt,iband,isppol)==me &
&               .and. me/=nmaster ) then
               action = 2
!              I am the master, and I receive the data
             elseif( mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me &
&               .and. me==nmaster ) then
               action=3
             end if

             if(action==1) then
!              I am the master node, and I have the data in cg or cg_disk
               if(mkmem/=0)then
!                Copy from kg to kg_disk
                 if(iband==1)kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
!                Copy from cg to cg_disk
                 do ipwnbd=1,npw_k*my_nspinor
                   cg_disk(1,(iband-1)*npw_k*my_nspinor+ipwnbd)= &
&                   cg(1,(iband-1)*npw_k*my_nspinor+ipwnbd+icg)
                   cg_disk(2,(iband-1)*npw_k*my_nspinor+ipwnbd)= &
&                   cg(2,(iband-1)*npw_k*my_nspinor+ipwnbd+icg)
                 end do
               end if
             end if  ! action=1

             if ( action==2.or.action==3) then
!              action=2 :  I am not the master, and I have the data => send to master
!              action=3 :  I am the master, and I receive the data
               call timab(48,1,tsec)
               if ( iband == 1 ) then
                 if ( mkmem/=0 .and. action==2) then
                   call xexch_mpi(kg(:,1+ikg:npw_k+ikg),3*npw_k,mpi_enreg%proc_distrb(ikpt,iband,isppol), &
&                   kg_disk,nmaster,spaceWorld,ierr)
                 else
                   call xexch_mpi(kg_disk,3*npw_k,mpi_enreg%proc_distrb(ikpt,iband,isppol),  &
&                   kg_disk,nmaster,spaceWorld,ierr)
                 end if
               end if       ! iband =1
               ipwnbd=(iband-1)*npw_k*my_nspinor
               if(mkmem/=0 .and. action==2)then
                 call xexch_mpi( cg(:,ipwnbd+icg+1:ipwnbd+icg+npw_k*my_nspinor),2*npw_k*my_nspinor &
&                 ,mpi_enreg%proc_distrb(ikpt,iband,isppol)                    &
&                 ,cg_disk(:,ipwnbd+1:ipwnbd+npw_k*my_nspinor),nmaster,spaceWorld,ierr)
               else
                 call xexch_mpi( cg_disk(:,ipwnbd+1:ipwnbd+npw_k*my_nspinor),2*npw_k*my_nspinor    &
&                 ,mpi_enreg%proc_distrb(ikpt,iband,isppol)                    &
&                 ,cg_disk(:,ipwnbd+1:ipwnbd+npw_k*my_nspinor),nmaster,spaceWorld,ierr)
               end if

               call timab(48,2,tsec)
             end if        ! action=2 or action=3



             if(accesswff ==1 ) then
!              I have the data in cg or cg_disk
               swrite=.false.
               if (nmaster == me) then
                 swrite=.true.
               end if
             end if

!            End of loop over bands
           end do

!          End of paralbd=1
         end if
       end if
#endif

!      Only the master will write to disk the final output wf file.
!      in MPI_IO case only swrite will write to disk the final output wf file.
       if(swrite) then
!        DEBUG
!        write(std_out,*) 'outwf : I am master and will write wf file'
!        ENDDEBUG
         if(formeig==0)then
           eig_k(1:nband_k)=eigen(1+band_index:nband_k+band_index)
           occ_k(1:nband_k)=occ(1+band_index:nband_k+band_index)
         else
           eig_k(1:2*nband_k*nband_k)=eigen(1+band_index:2*nband_k*nband_k+band_index)
         end if
         option=2
!        if (dtset%prtwf == 2 .and. mkmem/=0) option=4

         if (dtset%usewvl == 0) then
#if defined HAVE_MPI
           call rwwf(cg_disk,eig_k,formeig,0,0,ikpt,isppol,kg_disk,mband,mcg_disk,mpi_enreg, &
&           nband_k, nband_k,npw_k,my_nspinor,occ_k,option,1,tim_rwwf,wff2)
#elif !defined HAVE_MPI
           if(mkmem==0)then
             call rwwf(cg_disk,eig_k,formeig,0,0,ikpt,isppol,kg_disk,mband,mcg_disk,mpi_enreg, &
&             nband_k,nband_k, npw_k,my_nspinor,occ_k,2,1,tim_rwwf,wff2)
           else if(mkmem/=0)then
             kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
             call rwwf(cg,eig_k,formeig,0,icg,ikpt,isppol,kg_disk,mband,mcg,mpi_enreg,nband_k, &
&             nband_k, npw_k,my_nspinor,occ_k,option,1,tim_rwwf,wff2)
           end if
#endif
         else
           call wvl_write(dtset, eigen, mpi_enreg, option, hdr%rprimd, &
&           wff2, wfs, wvl, hdr%xred)
         end if
       end if

!      The wavefunctions for the present k point and spin are written
       if(response==0)band_index=band_index+nband_k
       if(response==1)band_index=band_index+2*nband_k*nband_k

       if (mkmem/=0) then

         sskip=1
#if defined HAVE_MPI
         if (dtset%usewvl == 0) then
           sskip=0
           if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))==0)sskip=1
         end if
#endif
         if(sskip==1)then
           icg=icg+npw_k*my_nspinor*nband_k
           ikg=ikg+npw_k
         end if

       end if !mkem/=0


#if defined HAVE_MPI_IO
       spacecomsender=spaceComm
       if (mpi_enreg%mode_para=='b') spacecomsender =mpi_enreg%comm_kpt
       call WffOffset(wff2,sender,spacecomsender,ierr)
#endif

     end do ! ikpt
   end do ! isppol
   ABI_DEALLOCATE(kg_disk)
   if(mkmem==0) ABI_DEALLOCATE(cg_disk)
#if defined HAVE_MPI
   if(mkmem/=0) ABI_DEALLOCATE(cg_disk)
#endif

   if(mkmem==0)  then
     ABI_DEALLOCATE(eig_dum)
     ABI_DEALLOCATE(occ_dum)
   end if
   ABI_DEALLOCATE(eig_k)
   ABI_DEALLOCATE(occ_k)

!  Write the (x,f) history

   if(me0==0 .and. nxfh>0 .and. response==0)then
     if (wff2%accesswff /= 2) then
#if defined HAVE_MPI_IO
       if(wff2%accesswff == 1 ) then
         close(unit=wff2%unwff)
!        the file is to be positioned at the terminal point
         open(unit=wff2%unwff,file=wff2%fname,form='unformatted',POSITION="APPEND")
       end if
#endif
       write(unit=wff2%unwff)nxfh
       do ixfh=1,nxfh
         write(unit=wff2%unwff)xfhist(:,:,:,ixfh)
       end do
#if defined HAVE_TRIO_NETCDF
     else if (wff2%accesswff == 2) then
       ncid_hdr=wff2%unwff

!      check if nxfh and xfhist are defined
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)

       if (ncerr /= NF90_NOERR) then
!        need to define everything
         ncerr = nf90_redef (ncid=ncid_hdr)
         call handle_ncerr(ncerr," outwf : going to define mode ")

         ncerr = nf90_def_dim(ncid=ncid_hdr,name="dim2inout",len=2,dimid=dim2inout_id)
         call handle_ncerr(ncerr," outwf : define dim2inout")
         ncerr = nf90_def_dim(ncid=ncid_hdr,name="mxfh",len=mxfh,dimid=mxfh_id)
         call handle_ncerr(ncerr," outwf : define mxfh")
         ncerr = nf90_def_dim(ncid=ncid_hdr,name="nxfh",len=nxfh,dimid=nxfh_id)
         call handle_ncerr(ncerr," outwf : define nxfh")
         ncerr = nf90_def_dim(ncid=ncid_hdr,name="xfdim2",len=xfdim2,dimid=xfdim2_id)
         call handle_ncerr(ncerr," outwf : define xfdim2")

         ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="dimr3",dimid=dimr3_id)
         call handle_ncerr(ncerr," outwf : inquire dimr3")

!        xfhist(3,natom+4,2,mxfh)
         ncerr = nf90_def_var(ncid=ncid_hdr,name="xfhist",xtype=NF90_DOUBLE,&
&         dimids=(/dimr3_id,xfdim2_id,dim2inout_id,mxfh_id/),varid=xfhist_id)
         call handle_ncerr(ncerr," outwf : define xfhist")

!        End define mode and go to data mode
         ncerr = nf90_enddef(ncid=ncid_hdr)
         call handle_ncerr(ncerr," outwf : enddef call ")
       else
!        check that the dimensions are correct
         ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
         call handle_ncerr(ncerr," outwf : inquire nxfh")
         ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
&         len=nxfh_tmp)
         call handle_ncerr(ncerr,"  outwf : get nxfh")
         ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="xfdim2",dimid=xfdim2_id)
         call handle_ncerr(ncerr," outwf : inquire xfdim2")
         ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=xfdim2_id,&
&         len=xfdim2_tmp)
         call handle_ncerr(ncerr,"  outwf : get xfdim2")
         ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="mxfh",dimid=mxfh_id)
         call handle_ncerr(ncerr," outwf : inquire mxfh")
         ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=mxfh_id,&
&         len=mxfh_tmp)
         call handle_ncerr(ncerr,"  outwf : get mxfh")
         ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="dim2inout",dimid=dim2inout_id)
         call handle_ncerr(ncerr," outwf : inquire dim2inout")
         ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=dim2inout_id,&
&         len=dim2inout_tmp)
         call handle_ncerr(ncerr,"  outwf : get dim2inout")

         ncerr = nf90_inq_varid(ncid=ncid_hdr,name="xfhist",varid=xfhist_id)
         call handle_ncerr(ncerr," outwf : inquire xfhist")

         if (mxfh_tmp /= mxfh .or. dim2inout_tmp /= 2 .or. xfdim2_tmp /= xfdim2) then
           stop
         end if

       end if

!      Now fill the data
       ncerr = nf90_put_var(ncid=ncid_hdr,varid=xfhist_id,values=xfhist,&
&       start=(/1,1,1,1/),count=(/3,xfdim2,2,nxfh/))
       call handle_ncerr(ncerr," outwf : fill xfhist")

!      end NETCDF definition ifdef
#endif
     end if
!    end accesswff if

   end if
!  end tmaster if

!  Close the wavefunction file (and do NOT delete it !)
   if (wff2%accesswff /= 2) then
     call WffClose(wff2,ierr)
#if defined HAVE_TRIO_NETCDF
   else if (wff2%accesswff == 2 .and. tmaster) then
     ncerr = nf90_close(wff2%unwff)
     call handle_ncerr(ncerr," close netcdf wavefunction file")
#endif
   end if
!  End condition of nstep>0
 end if

!Close the temporary data file, if any
 if (mkmem==0) then
   call WffDelete(wffnow,ierr)
 end if

!DEBUG
!write(std_out,*)' outwf : exit'
!ENDDEBUG

end subroutine outwf
!!***
