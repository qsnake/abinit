!{\src2tex{textfont=tt}}
!!****f* ABINIT/inwffil3
!! NAME
!! inwffil3
!!
!! FUNCTION
!! Reads eigenvalues from the wavefunction file.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (GMR,AR,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables in this dataset
!!  nkpt = number of k points
!!  istwfk(nkpt) = input option parameter that describes the storage of wfs
!!  mband = maximum number of bands
!!  nband(nkpt*nsppol) = number of bands at each k point, for each polarization
!!  npwarr(nkpt) = array holding npw for each k point, taking into account
!!   the effect of istwfk, and the spreading over processors
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  prtvol = controls print volume and debugging
!!  unwff1 = file unit number for the wavefunction file
!!  wffnm = name of the wavefunction file
!!
!! OUTPUT
!!  eigen(nband*nkpt*nsppol) = eigenvalues of the wavefunctions
!!  hdr <type(hdr_type)> = the header structured variable
!!  wff1 <type(wffile_type)> = structured info about the wavefunction file
!!
!! SIDE EFFECTS
!!  mpi_enreg = informations about MPI parallelization
!!
!! PARENTS
!!      suscep
!!
!! CHILDREN
!!      hdr_io,hdr_io_netcdf,leave_new,leave_test,timab,wffopen,wffreadeigk
!!      wffreadskipk,wrtout,xcomm_world,xdefineoff,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine inwffil3(dtset,eigen,hdr,istwfk,mband,mpi_enreg,nband,&
& nkpt,npwarr,nsppol,prtvol,wff1,unwff1,wffnm)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_wffile
#if defined HAVE_TRIO_NETCDF
  use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inwffil3'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkpt,nsppol,prtvol,unwff1
 character(len=fnlen),intent(in) :: wffnm
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(out) :: hdr
 type(wffile_type),intent(out) :: wff1
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt)
 real(dp),intent(out) :: eigen(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt_max=50,nspinor=1
 integer :: band_index,fform0,formeig,ierr,ikpt,ikpt0,ikptsp_old,isppol
 integer :: isppol0,istwf_k,master,me,nband_k,nkpt_eff,npw_k
 integer :: rdwr,spaceworld,tim_rwwf
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: buffer1(:),eig_k(:)

! *************************************************************************

 call timab(19,1,tsec)

 master=0
!Init me
 call xme_init(mpi_enreg,me)
!Init spaceworld
 call xcomm_world(mpi_enreg,spaceworld)

!Open the file in which the input wfs are stored
 call WffOpen(dtset%accesswff,spaceworld,wffnm,ierr,wff1,master,me,unwff1)

 write(message, '(a,a)' )&
& ' inwffil : will read wavefunctions from disk file ',trim(wffnm)
 call wrtout(ab_out,message,'COLL')

 rdwr=1
 if (dtset%accesswff /= IO_MODE_NETCDF) then
   call hdr_io(fform0,hdr,rdwr,wff1)
#if defined HAVE_TRIO_NETCDF
 else if (dtset%accesswff == IO_MODE_NETCDF) then
   call hdr_io_netcdf(fform0,hdr,rdwr,wff1)
 else if (dtset%accesswff == IO_MODE_ETSF) then
   write (std_out,*) "FIXME: ETSF I/O support in inwffil3"
#endif
 end if

!Define offsets, in case of MPI I/O
 formeig=0
 call xdefineOff(formeig,wff1,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)

 if(hdr%headform<40 .and. hdr%headform/=0)then
   write(message, '(7a,i5,3a)' ) ch10,&
&   ' inwffil3 : ERROR -',ch10,&
&   '  The header format of the file',wffnm,ch10,&
&   '  is headform=',hdr%headform,', that is, pre-v4.0.',ch10,&
&   '  Action : either generate a new file using a current ABINIT version, or use an old ABINIT version.'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

 nkpt_eff=nkpt
 if( (prtvol==0.or.prtvol==1) .and. nkpt_eff>nkpt_max ) nkpt_eff=nkpt_max

 band_index=0
 ikptsp_old=0

!Loop over spins
 do isppol=1,nsppol

!  Loop over k points
   do ikpt=1,nkpt

     nband_k=nband(ikpt+(isppol-1)*nkpt)
     npw_k=npwarr(ikpt)
     istwf_k=istwfk(ikpt)

     if(mpi_enreg%paral_compil_kpt==0)then
       if(ikpt<=nkpt_eff)then
         write(message, '(a,i6,a,i8,a,i4)' ) &
&         ' inwffil3: treating ',nband_k,' bands with npw=',npw_k,&
&         ' for ikpt=',ikpt
         call wrtout(ab_out,message,'COLL')
       else if(ikpt==nkpt_eff+1)then
         write(message, '(a)' ) &
&         ' inwffil : prtvol=0 or 1, do not print more k-points.'
         call wrtout(ab_out,message,'COLL')
       end if
     else if(mpi_enreg%paral_compil_kpt==1)then
       if(ikpt<=nkpt_eff)then
         write(message, '(a,i6,a,i8,a,i4,a,i4)' )&
&         ' inwffil3: treating ',nband_k,' bands with npw=',npw_k,&
&         ' for ikpt=',ikpt,' by node ',&
&         mpi_enreg%proc_distrb(ikpt,1,isppol)
         call wrtout(ab_out,message,'COLL')
       else if(ikpt==nkpt_eff+1)then
         write(message, '(16x,a)' ) &
&         ' inwffil : prtvol=0 or 1, do not print more k-points.'
         call wrtout(ab_out,message,'COLL')
       end if
       if(mpi_enreg%proc_distrb(ikpt,1,isppol)/=mpi_enreg%me)then
         eigen(1+band_index : nband_k+band_index) = zero
         band_index=band_index+nband_k
!        In the case this k point does not belong to me, cycle
         cycle
       end if
       if(ikpt<=nkpt_eff)then
         write(message, '(a,i6,a,i8,a,i4,a,i4)' ) &
&         ' inwffil3: treating ',nband_k,' bands with npw=',npw_k,&
&         ' for ikpt=',ikpt,' by node ',mpi_enreg%me
         call wrtout(std_out,message,'PERS')
       else if(ikpt==nkpt_eff+1)then
         write(message, '(a)' ) &
&         ' inwffil : prtvol=0 or 1, do not print more k-points.'
         call wrtout(std_out,message,'PERS')
       end if
     end if

!    Read the eigenvalues for this k point

     ABI_ALLOCATE(eig_k,(nband_k))

!    Skip wavefunctions for k-points not treated by this proc.
!    (from ikptsp_old+1 to ikpt+(isppol-1)*nkpt-1)
     if(ikptsp_old<ikpt+(isppol-1)*nkpt-1)then
       do ikpt0=ikptsp_old+1,ikpt+(isppol-1)*nkpt-1
         isppol0=1 ! For the time being, non spin-polarized systems only
         call WffReadSkipK(0,0,ikpt0,isppol0,mpi_enreg,wff1)
       end do
     end if

     tim_rwwf=0
!    DEBUG
!    write(std_out,*)' inwffil3 : will call WffReadEig '
!    ENDDEBUG
     call WffReadEigK(eig_k,formeig,hdr%headform,ikpt,isppol,nband_k,mpi_enreg,nband_k,tim_rwwf,wff1)

     eigen(1+band_index:nband_k+band_index)=eig_k(:)
     ABI_DEALLOCATE(eig_k)

     ikptsp_old=ikpt+(isppol-1)*nkpt
     band_index=band_index+nband_k

!    End of the k loop
   end do

!  End of spin loop
 end do

 if(mpi_enreg%paral_compil_kpt == 1) then
!  BEGIN TF_CHANGES
   call leave_test()
!  END TF_CHANGES
   write(message, '(a)' ) 'inwffil3: loop on k-points and spins done in parallel'
   call wrtout(std_out,message,'COLL')

!  Prepare transmission of eigen
   band_index=0
   ABI_ALLOCATE(buffer1,(mband*nkpt*nsppol))
   do isppol=1,nsppol
     do ikpt=1,nkpt
       nband_k=nband(ikpt+(isppol-1)*nkpt)
       buffer1(band_index+1:band_index+nband_k)=eigen(band_index+1:band_index+nband_k)
       band_index=band_index+nband_k
     end do
   end do

!  Build sum of everything
   call timab(48,1,tsec)
   call xsum_mpi(buffer1,eigen,band_index,spaceworld,ierr)
   call timab(48,2,tsec)
 end if

!****************************************************************************

 write(message,*)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')


 call timab(19,2,tsec)

end subroutine inwffil3
!!***
