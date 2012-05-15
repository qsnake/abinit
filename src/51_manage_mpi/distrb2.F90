!{\src2tex{textfont=tt}}
!!****f* ABINIT/distrb2
!! NAME
!!  distrb2
!!
!! FUNCTION
!!  This routine creates the tabs of repartition of processors
!!  for sharing the jobs on k-points, spins and bands.
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2012 ABINIT group (AR,XG,MB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mband = maximum number of bands
!!  nband(nkpt*nsppol) = number of bands per k point, for each spin
!!  nkpt = number of k-points
!!  nsppol = 1 for unpolarized, 2 for polarized
!!
!! SIDE EFFECTS
!!  mpi_enreg = informations about MPI parallelization
!!   mpi_enreg%proc_distrb(nkpt,mband,nsppol)=number of the processor
!!       that will treat each band in each k point.
!!   mpi_enreg%nproc_kpt is set
!!
!!
!! NOTES
!!  For the time being, the band parallelisation works only
!!  when the number of bands is identical for spin up and spin down
!!  at the same k point. The problem is the most clearly seen
!!  in the kpgio routine, where a different parallel repartition
!!  of k points for spin up and spin down would conflict with the
!!  present computation of k+G sphere, independent of the spin.
!!
!! PARENTS
!!      eig2tot,initmpi_gs,invars1,invars2m,loper3,nonlinear,respfn,suscep
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine distrb2(mband, nband, nkpt, nsppol, mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'distrb2'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: mband,nkpt,nsppol
 integer,intent(in) :: nband(nkpt*nsppol)
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: ind,ios,nband_k,proc_max,proc_min
 integer :: iiband,iikpt,iisppol,nbsteps,nproc,nstates
 integer :: kpt_distrb(nkpt)
 integer,save :: file_exist
 logical,save :: first=.true.
 character(len=500) :: message

!******************************************************************
!BEGIN EXECUTABLE SECTION
!DEBUG
!write(std_out,*)' distrb2: enter '
!write(std_out,*)' mpi_enreg%paralbd=',mpi_enreg%paralbd
!write(std_out,*)' mpi_enreg%paral_compil_respfn=',mpi_enreg%paral_compil_respfn
!write(std_out,*)' mpi_enreg%paral_compil_kpt=',mpi_enreg%paral_compil_kpt
!write(std_out,*)' mpi_enreg%paral_compil_fft=',mpi_enreg%paral_compil_fft
!write(std_out,*)' mpi_enreg%nproc_kpt=',mpi_enreg%nproc_kpt
!ENDDEBUG

!This is to take in to account the possible distribution
!over images of the cell
 if (mpi_enreg%paral_img==0) then
   nproc=mpi_enreg%nproc
 else
   nproc=mpi_enreg%nproc_one_img
 end if

 if (mpi_enreg%paral_compil_fft==0) then
   mpi_enreg%nproc_kpt = nproc
 end if

!Initialization of proc_distrb
 do iisppol=1,nsppol
   do iiband=1,mband
     do iikpt=1,nkpt
       mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=mpi_enreg%nproc_kpt-1
     end do
   end do
 end do
!That s all for an empty communication space
 if (nproc==0) return

!Testing section
 if (mpi_enreg%paralbd >= 1) then
!  Calculation of number of states for isppol

!  MT 2011-29-11: change nstates in order to make parallelization over bands
!  work even is nproc is not a dividor of nstates
   if (nproc<=nkpt*nsppol) then
     nstates = 0
     do iisppol=1,nsppol
       do iikpt=1,nkpt
         nstates= nstates + nband(iikpt+(iisppol-1)*nkpt)
       end do
     end do
   else
     nstates=nsppol*nkpt*mband
   end if

!  Tests
   if (mpi_enreg%nproc_kpt > nstates) then
!    Too much proc. with respect to nstates
     write(message, '(a,a,a,a,i4,a,i4,a,a)' ) ch10,&
&     ' distrb2: WARNING -',ch10,&
&     '  nproc_kpt=',mpi_enreg%nproc_kpt,' >= nstates=',nstates,ch10,&
&     '  The number of processors is larger than nstates. This is a waste.'
     call wrtout(std_out,message,'COLL')
   end if

 elseif (mpi_enreg%paralbd==0) then

!  Check if nkpt and nproc_kpt match
   if(mpi_enreg%nproc_kpt>nkpt*nsppol) then
!    Too much proc. with respect to nkpt
     write(message, '(a,a,a,a,i4,a,i4,a,i4,a,a)' ) ch10,&
&     ' distrb2: WARNING -',ch10,&
&     '  nproc_kpt=',mpi_enreg%nproc_kpt,' >= nkpt=',nkpt,'* nsppol=',nsppol,ch10,&
&     '  The number of processors is larger than nkpt. This is a waste.'
     call wrtout(std_out,message,'COLL')
   elseif(mod(nkpt*nsppol,mpi_enreg%nproc_kpt)/=0) then
!    nkpt not a multiple of nproc_kpt
     write(message, '(4a,i5,a,i5,3a)' ) ch10,&
&     ' distrb2: WARNING -',ch10,&
&     '  nkpt*nsppol (', nkpt*nsppol, ') is not a multiple of nproc_kpt (',&
&     mpi_enreg%nproc_kpt, ')', ch10,&
&     '  The k-point parallelisation is not efficient.'
     call wrtout(std_out,message,'COLL')
   end if

 end if ! mpi_enreg%paralbd

 if (mpi_enreg%paralbd /= 1 .and. mpi_enreg%paralbd /= 0) then
   if (mod(mpi_enreg%nproc_kpt,nkpt*nsppol) /= 0) then
     write(message, '(a,a,a,a,i4,a,i4,a,a,i4,a)' ) ch10,&
&     ' distrb: WARNING -',ch10,&
&     '  nproc_kpt=',mpi_enreg%nproc_kpt,' nkpt*nsppol=',nkpt*nsppol,ch10, &
&     '  The number of processors is not a mulptiple of nkpt*nsppol. '&
&     ,  mod(mpi_enreg%nproc_kpt,nkpt*nsppol),' proc is unusefull. This is a waste.'
     call wrtout(std_out,message,'COLL')
   end if
   if (mpi_enreg%nproc_kpt > (mpi_enreg%paralbd*nkpt*nsppol)) then
     write(message, '(a,a,a,a,i4,a,i4,a,a)' ) ch10,&
&     ' distrb2: WARNING -',ch10,&
&     '  nproc_kpt=',mpi_enreg%nproc_kpt,' > nkpt*nbdblock*nsppol', &
&     nkpt*mpi_enreg%paralbd*nsppol,ch10, &
&     '  The number of processors is larger than nkpt*nbdblock*nsppol. This is a waste.'
     call wrtout(std_out,message,'COLL')
   end if
 end if

 if(mpi_enreg%paral_compil_respfn == 1) then
   if(mpi_enreg%nproc_kpt<mpi_enreg%ngroup_respfn) then
     write(message, '(a,a,a,a,i4,a,i4,a,a)' ) ch10,&
&     ' distrb2: ERROR -',ch10,&
&     '  nproc_kpt=',mpi_enreg%nproc_kpt,' < ngroup_respfn=',mpi_enreg%ngroup_respfn,ch10,&
&     '  The number of processors is smaller than the number of computation groups.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if
!End of testing section

!Inquire whether there exist a file containing the processor distribution
 if (first) then
!  Case first time : test file to do
!  Open the file containing the k-point distribution
   open(unit=tmp_unit,file='kpt_distrb',form='formatted',status='old',iostat=ios)
   if(ios==0) then
!    'kpt_distrb' file exists
     file_exist=1
     close(tmp_unit)
   else
     file_exist=0
   end if
   first=.false.
 end if

!Initialize the processor distribution, either from a file, or from an algorithm
 if (file_exist == 1) then

   open(unit=tmp_unit,file='kpt_distrb',form='formatted',status='old',iostat=ios)
   rewind(unit=tmp_unit)
   if (mpi_enreg%paralbd >= 1) then
!    -> read bands distribution
     read(tmp_unit,*) mpi_enreg%proc_distrb
   else
     read(tmp_unit,*) kpt_distrb
   end if
   close(tmp_unit)
   proc_max=0
   proc_min=mpi_enreg%nproc_kpt
!  -> determine the range of proc. requested
   if (mpi_enreg%paralbd >= 1) then
     do iisppol=1,nsppol
       do iikpt=1,nkpt
         nband_k = nband(iikpt+(iisppol-1)*nkpt)
         proc_max=maxval(mpi_enreg%proc_distrb(iikpt,1:nband_k,iisppol))
         proc_min=minval(mpi_enreg%proc_distrb(iikpt,1:nband_k,iisppol))
       end do
     end do
   else
     proc_max=maxval(kpt_distrb(1:nkpt))
     proc_min=minval(kpt_distrb(1:nkpt))
!    -> fill the tab proc_distrb with kpt_distrb
     do iisppol=1,nsppol
       do iikpt=1,nkpt
         nband_k = nband(iikpt+(iisppol-1)*nkpt)
         do iiband=1,nband_k
           mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=kpt_distrb(iikpt)
         end do
       end do
     end do
   end if ! mpi_enreg%paralbd

   if(proc_max>(mpi_enreg%nproc_kpt-1)) then
!    Too much proc. requested
     write(message, '(a,a,a,a,a,a,i4,a,a,a)' ) ch10,&
&     ' distrb2: ERROR -',ch10,&
&     '  The number of processors mentioned in the kpt_distrb file',ch10,&
&     '  must be lower or equal to the actual number of processors =',&
&     mpi_enreg%nproc_kpt-1,ch10,&
&     '  Action : change the kpt_distrb file, or increase the',&
&     '  number of processors.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   if(proc_max/=(mpi_enreg%nproc_kpt-1)) then
!    Too few proc. used
     write(message, '(a,a,a,a,i4,a,a,a,i4,a,a,a)' ) ch10,&
&     ' distrb2: ERROR -',ch10,&
&     '  Only ',proc_max+1,' processors are used (from kpt_distrb file),',&
&     ch10,'  when',mpi_enreg%nproc_kpt,' processors are available.',ch10,&
&     '  Action : adjust number of processors and kpt_distrb file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   if(proc_min<0) then
     write(message, '(a,a,a,a,a,a)' ) ch10,&
&     ' distrb2: ERROR -',ch10,&
&     '  The number of processors must be bigger than 0 in kpt_distrb file.',&
&     ch10,' Action : modify kpt_distrb file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

 else

!  'kpt_distrb' file does not exist

   if (mpi_enreg%paralbd==1) then

     if(mpi_enreg%paral_compil_respfn == 1) then
       nbsteps = (nstates*mpi_enreg%ngroup_respfn) / mpi_enreg%nproc_kpt
       if (mod(nstates,mpi_enreg%nproc_kpt / mpi_enreg%ngroup_respfn) /=0) then
         nbsteps=nbsteps+1
       end if
     else ! mpi_enreg%paral_compil_respfn
       nbsteps=nstates/mpi_enreg%nproc_kpt
       if (mod(nstates,mpi_enreg%nproc_kpt) /=0) then
         nbsteps=nbsteps+1
       end if
     end if

!    XG060807 : new processor distribution, correct for nsppol=2
!    One relies on the fact that the number of bands for the same k point,
!    spin up and spin down, is equal
     ind=0
     do iikpt=1,nkpt
       nband_k = nband(iikpt)
       do iiband=1,nband_k
         mpi_enreg%proc_distrb(iikpt,iiband,1)=ind/nbsteps
         ind = ind + 1
         if(nsppol==2)then
           mpi_enreg%proc_distrb(iikpt,iiband,2)=mpi_enreg%nproc_kpt-mpi_enreg%proc_distrb(iikpt,iiband,1)-1
         end if
       end do
!      This ensure that each proc either treats all bands of given k-points or some bands for ONE k-point
       if (nproc>nkpt*nsppol.and.mod(ind,nbsteps)/=0) ind=((ind/nbsteps)+1)*nbsteps
     end do
!    XG060807 : OLD CODING
!    ind=0
!    do iisppol=1,nsppol
!    do iikpt=1,nkpt
!    nband_k = nband(iikpt+(iisppol-1)*nkpt)
!    do iiband=1,nband_k
!    mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=ind/nbsteps
!    ind = ind + 1
!    end do
!    end do
!    end do
!    XG060807 : END OF OLD CODING

   elseif (mpi_enreg%paralbd==0) then

     nbsteps=(nsppol*nkpt)/mpi_enreg%nproc_kpt;
     if (mod((nsppol*nkpt),mpi_enreg%nproc_kpt) /=0) then
       nbsteps=nbsteps+1
     end if

!    XG060807 : new processor distribution, correct for nsppol=2
     ind=0
     do iikpt=1,nkpt
       nband_k = nband(iikpt)
       do iiband=1,nband_k
         mpi_enreg%proc_distrb(iikpt,iiband,1)=ind/nbsteps
         if(nsppol==2)then
           mpi_enreg%proc_distrb(iikpt,iiband,2)=mpi_enreg%nproc_kpt-mpi_enreg%proc_distrb(iikpt,iiband,1)-1
         end if
       end do
       ind=ind + 1
     end do
!    XG060807 : OLD CODING
!    ind=0
!    do iisppol=1,nsppol
!    do iikpt=1,nkpt
!    nband_k = nband(iikpt+(iisppol-1)*nkpt)
!    do iiband=1,nband_k
!    Distribute k-points homogeneously
!    proc_distrb(iikpt,iiband,iisppol)=mod(iikpt-1,mpi_enreg%nproc_kpt)
!    mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=ind/nbsteps
!    DEBUG
!    write(std_out,fmt='(a8,i8,a,3i8,a,i8)')' proc ',mpi_enreg%me       &
!    &          ,' proc_distrb0( ', iikpt,iiband,iisppol,') proc : '          &
!    &          , mpi_enreg%proc_distrb(iikpt,iiband,iisppol)
!    ENDDEBUG
!    end do
!    ind=ind + 1
!    end do
!    end do
!    XG060807 : END OF OLD CODING

   else ! mpi_enreg%paralbd

     mpi_enreg%nproc_per_kpt=mpi_enreg%nproc_kpt / (nkpt*nsppol)
     if (mpi_enreg%nproc_per_kpt == 0 ) then
       write(message, '(a,a,a,a,a,i4,a,a)' ) ch10,&
&       ' distrb: ERROR -',ch10,&
&       ' The number of processors is not sufficient to activate the band parallelism',ch10,&
&       mpi_enreg%nproc_kpt,ch10,&
&       ' Action : increase the number of processors, or suppress band parallelism.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if (mpi_enreg%nproc_per_kpt > mpi_enreg%paralbd) then
       mpi_enreg%nproc_per_kpt=mpi_enreg%paralbd
     end if

!    XG060807 : new processor distribution, correct for nsppol=2
!    One relies on the fact that the number of bands for the same k point,
!    spin up and spin down, is equal
     do iikpt=1,nkpt
       nband_k = nband(iikpt)
       do iiband=1,nband_k
         mpi_enreg%proc_distrb(iikpt,iiband,1)= &
&         (iikpt-1)*mpi_enreg%nproc_per_kpt + &
&         mod((iiband-1),mpi_enreg%nproc_per_kpt)
         if(nsppol==2)then
           mpi_enreg%proc_distrb(iikpt,iiband,2)=mpi_enreg%nproc_kpt-mpi_enreg%proc_distrb(iikpt,iiband,1)-1
         end if
       end do
     end do
!    XG060807 : OLD CODING
!    do iisppol=1,nsppol
!    do iikpt=1,nkpt
!    nband_k = nband(iikpt+(iisppol-1)*nkpt)
!    do iiband=1,nband_k
!    mpi_enreg%proc_distrb(iikpt,iiband,iisppol)= &
!    &                     (iisppol-1)*mpi_enreg%nproc_per_kpt*nkpt + &
!    &                     (iikpt-1)*mpi_enreg%nproc_per_kpt + &
!    &                     mod((iiband-1),mpi_enreg%nproc_per_kpt)
!    end do
!    end do
!    end do
!    XG060807 : END OF OLD CODING

   end if ! mpi_enreg%paralbd

 end if ! file_exist

!DEBUG
!write(std_out,*)' distrb2: exit '
!ENDDEBUG

end subroutine distrb2
!!***
