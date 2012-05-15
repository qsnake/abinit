!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_veloc_tr
!!
!! NAME
!! get_veloc_tr
!!
!! FUNCTION
!!  calculate the (in) and (out) velocity factors for transport
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (JPC)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  elph_ds
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%k_fine%nkpt = number of kpts included in the FS integration
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%minFSband = index of the lowest FS band
!!    elph_ds%nqpt_full  = number of Q pts 
!!    elph_ds%nqptirred  = number of irreducible Q pts 
!!  to index the GS electronic states :
!!  nband =full number of bands 
!!  kphon_full2irr = mapping of full FS kpts to irreducible ones
!!   FSfullpqtofull = mapping of k + q to k
!!   FSirredtoGS = mapping of irreducible kpoints to GS set
!! OUTPUT
!! elph_tr_ds%eta_trout = prefactor for out gkk_qpt
!! elph_tr_ds%eta_trin = prefactor for in_gkk_qpt
!! elph_tr_ds%FSelecveloc_sq = avergae FS electronic velocity
!! SIDE EFFECTS
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      leave_new,read_el_veloc,wrtout
!!
!! NOTES
!!   
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine get_veloc_tr(elph_ds,mpi_enreg,nband,elph_tr_ds)

 use m_profiling

  use defs_datatypes
  use defs_abitypes
  use defs_elphon

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_veloc_tr'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_77_ddb, except_this_one => get_veloc_tr
!End of the abilint section

  implicit none


!Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nband
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  type(elph_type),intent(in) :: elph_ds
  type(elph_tr_type):: elph_tr_ds


!Local variables-------------------------------
  !scalars
  integer :: ikpt_fine
  integer :: ib1,fib1,isppol, ii
  real(dp) :: eta2
  !arrays
  real(dp) :: elvelock(3)
  character(len=500) :: message

! *********************************************************************

!check inputs
!TODO: should be done at earlier stage of initialization and checking
 if (elph_ds%ngkkband /= elph_ds%nFSband) then
   write (message,'(a)') 'Error: need to keep electron band dependency in memory for transport calculations'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 write(std_out,*)'reading of electronic velocities'

!
!FIXME: should have gkqwrite affect the storage of el_veloc, now that the trin trout variables have been eliminated
!
 ABI_ALLOCATE(elph_tr_ds%el_veloc,(elph_ds%k_fine%nkpt,nband,3,elph_ds%nsppol))
 call read_el_veloc(mpi_enreg,nband,elph_ds%k_fine%nkpt,elph_ds%k_fine%kpt,elph_ds%nsppol,elph_tr_ds)

!precalculate the Fermi speed modulus squared
 ABI_ALLOCATE(elph_tr_ds%FSelecveloc_sq,(3,elph_ds%nsppol))
 elph_tr_ds%FSelecveloc_sq = zero
 do isppol=1,elph_ds%nsppol
   do ikpt_fine=1,elph_ds%k_fine%nkpt
     do ib1=1,elph_ds%nFSband
       fib1=ib1+elph_ds%minFSband-1
       elvelock(:)=elph_tr_ds%el_veloc(ikpt_fine,fib1,:,isppol)
       do ii=1, 3
         eta2=elvelock(ii)*elvelock(ii)
         elph_tr_ds%FSelecveloc_sq(ii, isppol)=elph_tr_ds%FSelecveloc_sq(ii, isppol)&
&         +eta2*elph_ds%k_fine%wtk(ib1,ikpt_fine,isppol)
       end do
     end do
   end do
   elph_tr_ds%FSelecveloc_sq(:,isppol) = elph_tr_ds%FSelecveloc_sq(:,isppol)/elph_ds%k_fine%nkpt/elph_ds%n0(isppol)
!  for factor 1/elph_ds%n0(isppol) see eq 12 of Allen prb 17 3725: sum of v**2 over all k gives n0 times FSelecveloc_sq
 end do ! end isppol
 write (std_out,*) '  get_veloc_tr: FSelecveloc_sq ', elph_tr_ds%FSelecveloc_sq 

 write(std_out,*)'out of get_veloc_tr'

end subroutine get_veloc_tr
!!***
