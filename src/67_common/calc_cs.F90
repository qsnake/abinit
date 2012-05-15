!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_cs
!! NAME
!! calc_cs
!!
!! FUNCTION
!! calculation and output of chemical shielding tensor at each atomic site
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (SM,VR,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cg(2,mcg)=planewave coefficients of wavefunctions
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kptns(3,nkpt)=coordinates of k points (reduced units)
!!  mcg=size of wave-functions array (cg) =mpw**mband*mkmem*nsppol
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  occopt=option used for band occupation
!!  ucvol=unit cell volume (bohr**3)
!!  usepaw=option controlling use of PAW
!!  typat(natom)=atom type of each atom in the cell
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      leave_new,mkkpg,wrtout,xcomm_world
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine calc_cs(cg,dtefield,gprimd,kg,kptns,mcg,mkmem,mpi_enreg,mpw,natom,nkpt,npwarr,&
&                   occopt,usepaw)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_linalg_interfaces
 use m_xmpi
 use m_efield

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_cs'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,mkmem,mpw,natom,nkpt,occopt,usepaw
 type(efield_type),intent(in) :: dtefield
 type(MPI_type), intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mcg),gprimd(3,3),kptns(3,nkpt)

!Local variables-------------------------------
!scalars
 integer :: bpw,iatom,iband,icg,idir,ikpt,kpw,mu,nkpg,npw_ctr,npw_k,spaceComm
 real(dp) :: ldg2
 complex(dpc) :: bra_c,ket_c
 character(len=500) :: message
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: bra_kg(3),dg(3),ket_kg(3),kpoint(3)
 real(dp),allocatable :: cs_tensor(:,:,:),kpg_k(:,:)

! ************************************************************************

!DEBUG
!write(std_out,*)' calc_cs : enter'
!ENDDEBUG
!Compatibility tests
 if (usepaw /= 1) then
   write (message,'(4a)')' calc_cs : ERROR- ',ch10,&
&   ' usepaw /= 1 but CS calculation requires PAW ',ch10
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if (occopt > 1) then
   write (message,'(4a)')' calc_cs : ERROR- ',ch10,&
&   ' occopt > 1 but CS calculation requires occopt = 1 or 0 (no metals for now, sorry !)  ',ch10
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if (mkmem == 0) then
   write (message,'(4a)')' calc_cs : ERROR- ',ch10,&
&   ' mkmem == 0 but CS calculation requires mkmem > 0  ',ch10
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 write(message,'(a,a,a)' ) ch10,' Chemical Shielding Calculation ',ch10
 call wrtout(ab_out,message,'COLL')

 call xcomm_world(mpi_enreg,spaceComm)

 nkpg = 3
 npw_ctr = 0
 icg = 0
 ABI_ALLOCATE(cs_tensor,(2,3,natom))
!complex number (2)
!direction (3)
!on each atom (natom)
!note that this calculation is for whatever direction of mag field
!was imposed, in general will to run the finite mag field computation
!three times, one in each of x, y, z in order to compute the full
!3x3 cs tensor
 cs_tensor(:,:,:) = zero

 do ikpt = 1, nkpt
   kpoint(:) = kptns(:,ikpt) ! this may be inconsistent with kpoints used in update_mmat in kptopt /= 3 case
   npw_k = npwarr(ikpt)
   ABI_ALLOCATE(kg_k,(3,npw_k))
   ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
   kg_k(:,1:npw_k) = kg(:,npw_ctr+1:npw_ctr+npw_k)

   call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k) ! construct k+G at this k point

   do iband = 1, dtefield%nband_occ

     do kpw = 1, npw_k
       ket_c = cmplx(cg(1,icg+(iband-1)*npw_k+kpw),cg(2,icg+(iband-1)*npw_k+kpw))
       do mu=1,3
         ket_kg(mu)=dot_product(kpg_k(kpw,:),gprimd(mu,:))
       end do 
       do bpw = 1, npw_k
         bra_c = cmplx(cg(1,icg+(iband-1)*npw_k+bpw),cg(2,icg+(iband-1)*npw_k+bpw))
         do mu=1,3
           bra_kg(mu)=dot_product(kpg_k(bpw,:),gprimd(mu,:))
         end do 
         dg(:) = ket_kg(:) - bra_kg
         ldg2 = dot_product(dg,dg)
         if(ldg2 < tol8) cycle
         do idir = 1, 3
           do iatom = 1, natom

           end do ! end loop over atoms
         end do ! end loop over directions
       end do ! end loop over bra plane waves
     end do ! end loop over ket plane waves

   end do ! end loop over bands

   npw_ctr = npw_ctr + npw_k
   icg = icg + npw_k*dtefield%nband_occ
   ABI_DEALLOCATE(kg_k)
   ABI_DEALLOCATE(kpg_k)
 end do ! end loop over k points

!final deallocations
 ABI_DEALLOCATE(cs_tensor)
 end subroutine calc_cs
!!***
