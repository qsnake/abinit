!{\src2tex{textfont=tt}}
!!****f* ABINIT/xdefineOff
!! NAME
!!  xdefineOff
!!
!! FUNCTION
!!  In case of MPI I/O, defines the offset for each processor
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2012 ABINIT group (MB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  formeig option (format of the eigenvalues and occupations) :
!!   0 => ground-state format (initialisation of eigenvectors with
!!        random numbers, vector of eigenvalues, occupations are present)
!!   1 => respfn format (initialisation of eigenvectors with 0 s,
!!        hermitian matrix of eigenvalues)
!!  nkpt = number of k points
!!  nspinor = total number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  nband(nkpt*nsppol) = number of bands at each k point, for each polarization
!!  npwarr(nkpt) = number of planewaves at each k point
!!  mpi_enreg <type(MPI_type)> = informations about MPI parallelization
!!
!! OUTPUT
!!  (no output)
!!
!! SIDE EFFECTS
!!  wff <type(wffile_type)> =
!!
!! PARENTS
!!      ctocprj,dyfnl3,eltfrkin3,eltfrnl3,energy,forstrnps,inwffil,inwffil3
!!      ladielmt,lavnl,mkrho,mkrho3,nselt3,nstdy3,nstpaw3,optics_paw
!!      optics_vloc,outwf,pawmkaewf,prctfvw1,prctfvw2,rhofermi3,tddft,uderiv
!!      vtorho,vtorho3
!!
!! CHILDREN
!!      xcomm_world
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine xdefineOff(formeig,wff,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_wffile

#if defined HAVE_MPI2 && defined HAVE_MPI_IO
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xdefineOff'
 use interfaces_51_manage_mpi, except_this_one => xdefineOff
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 && defined HAVE_MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer, intent(in) ::  nsppol,nkpt,nspinor,formeig
 integer, intent(in) ::  nband(nkpt*nsppol),npwarr(nkpt)
 type(wffile_type),intent(inout) :: wff
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables-------------------------------
#if defined HAVE_MPI_IO
!scalars
 integer :: comm,iproc
 integer :: nband_k,npw_k,nproc,me,ipp
 integer :: nbrec,isppol,ikpt,nbint,nbreal,nbd,ippband
 integer :: nrecnpw,nreckg
 integer(kind=MPI_OFFSET_KIND) :: pos_start
!arrays
 integer(kind=MPI_OFFSET_KIND),allocatable  :: offproc(:)
#endif

! *************************************************************************
!nbOct_int octet number of int value
!nbOct_dp octet number of dp value
!nbOct_ch octet number of character value
!lght_recs length of record

 if(.false.)write(std_out,*)wff%me,mpi_enreg%nproc,formeig,nband,npwarr,nspinor,nkpt
#if defined HAVE_MPI_IO
 if(wff%accesswff==IO_MODE_MPI)then

   call xcomm_world(mpi_enreg,comm,myrank=me,mysize=nproc)
   pos_start=wff%offwff

   ABI_ALLOCATE(offproc,(0:nproc))
   offproc = 0
   nbrec =2
   nrecnpw=3+nbrec

   do isppol=1,nsppol
     do ikpt=1,nkpt
       nband_k=nband(ikpt+(isppol-1)*nkpt)
       npw_k=npwarr(ikpt)
       iproc=mpi_enreg%proc_distrb(ikpt,1,isppol)
       if (mpi_enreg%paralbd==1) iproc=mpi_enreg%proc_distrb(ikpt,1,isppol)
!      record kg
       nreckg=nbrec+ wff%kgwff*3*npw_k

!      Record npw,nspinor,nband, Record kg
       offproc(iproc) = offproc(iproc) + wff%nbOct_int*(nrecnpw+nreckg)

       if (formeig == 0) then
!        Records eigen,occ
         nbint=nbrec
         nbreal =  2 *nband_k
         offproc(iproc) = offproc(iproc) + (wff%nbOct_int*nbint+wff%nbOct_dp*nbreal)

!        Records cg
         offproc(iproc) = offproc(iproc) &
&         + (wff%nbOct_int*nbrec+wff%nbOct_dp*2*npw_k*nspinor)*nband_k

         ippband=iproc
         do nbd=1,nband_k
           ipp=mpi_enreg%proc_distrb(ikpt,nbd,isppol)
           if (ipp /= ippband ) then
             ippband=ipp
             offproc(ippband)=offproc(ippband)+ wff%nbOct_int*(nrecnpw+nreckg)
             offproc(ippband) = offproc(ippband) + (wff%nbOct_int*nbint &
&             +wff%nbOct_dp*nbreal)
             offproc(ippband) = offproc(ippband) + (wff%nbOct_int*nbrec &
&             + wff%nbOct_dp*2*npw_k*nspinor)*nband_k
           end if
         end do
       else if (formeig == 1) then
!        record eigen
         offproc(iproc) = offproc(iproc) + (wff%nbOct_int*2*nbrec  &
&         + wff%nbOct_dp*2*npw_k*nspinor &
&         + wff%nbOct_dp*2*nband_k)*nband_k
         ippband=iproc
         do nbd=1,nband_k
           ipp=mpi_enreg%proc_distrb(ikpt,nbd,isppol)
           if (ipp /= ippband) then
             ippband=ipp
             offproc(ippband)=offproc(ippband)+ wff%nbOct_int*(nrecnpw+nreckg)
             offproc(ippband) = offproc(ippband) + (wff%nbOct_int*2*nbrec  &
&             + wff%nbOct_dp*2*npw_k*nspinor &
&             + wff%nbOct_dp*2*nband_k)*nband_k
           end if
         end do
       end if   ! formeig
     end do ! ikpt

   end do ! isppol

!  pos_start=wff%offwff
!  wff%offwff = pos_start

   if (me/=0)then
     do iproc=0,me-1
       wff%offwff=wff%offwff+offproc(iproc)
     end do
   end if
   ABI_DEALLOCATE(offproc)

 end if ! accesswff
#endif

end subroutine xdefineOff
!!***
