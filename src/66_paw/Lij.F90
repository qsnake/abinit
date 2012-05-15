!{\src2tex{textfont=tt}}
!!****f* ABINIT/Lij
!! NAME
!! Lij
!!
!! FUNCTION
!! Routine which computes PAW onsite angular momentum expectation values
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ntypat :: number of types of atoms
!!  type(pseudopotential_type) :: psps Data on pseudopotentials
!!  type(pawrad_type) :: pawrad(ntypat) Data on PAW radial grid
!!  type(pawtab_type) :: pawtab(ntypat) Tabulated PAW data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  type(efield_type) :: dtefield. The onsite Lij values and Lij/r^3
!!                                 Lij/r^3 are stored in dtefield
!!
!! NOTES
!!
!! PARENTS
!!      initberry
!!
!! CHILDREN
!!      deducer0,simp_gen,slxyzs
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine Lij(dtefield,ntypat,pawrad,pawtab,psps)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_efield

 use m_radmesh,   only : simp_gen, deducer0

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Lij'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: ntypat
 type(efield_type),intent(inout) :: dtefield
 type(pseudopotential_type),intent(in) :: psps

!arrays
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables---------------------------
!scalars
 integer :: idir, itypat, ilmn,il,im,iln,ilm, jlmn,jl,jm,jlm,jln,j0lmn
 integer :: klmn,kln, mesh_size
 real(dp) :: intg,intgr3
 complex(dpc) :: lms
!arrays
 real(dp),allocatable :: ff(:)

! *************************************************************************

!loop over types of atoms in cell
 do itypat = 1, ntypat
   mesh_size=pawrad(itypat)%mesh_size
   ABI_ALLOCATE(ff,(mesh_size))

!  loop over basis state pairs for this type
   do jlmn=1,pawtab(itypat)%lmn_size
     jl= psps%indlmn(1,jlmn,itypat)
     jm=psps%indlmn(2,jlmn,itypat)
     jlm = psps%indlmn(4,jlmn,itypat)
     jln=psps%indlmn(5,jlmn,itypat)
     j0lmn=jlmn*(jlmn-1)/2
     do ilmn=1,jlmn
       il= psps%indlmn(1,ilmn,itypat)
       im=psps%indlmn(2,ilmn,itypat)
       iln=psps%indlmn(5,ilmn,itypat)
       ilm = psps%indlmn(4,ilmn,itypat)
       klmn=j0lmn+ilmn
       kln = pawtab(itypat)%indklmn(2,klmn)

!      Computation of <phi_i|phi_j>- <tphi_i|tphi_j> radial integral
!      this is NOT the same as the sij non-local overlap, because that also
!      involves an angular integral over the S_i*S_j spherical harmonics
       ff(2:mesh_size)=pawtab(itypat)%phiphj(2:mesh_size,kln)-&
&       pawtab(itypat)%tphitphj(2:mesh_size,kln)
       call deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))

       ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln)-&
&       pawtab(itypat)%tphitphj(2:mesh_size,kln))/&
&       pawrad(itypat)%rad(2:mesh_size)**3
       call deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intgr3,ff,pawrad(itypat))

       do idir = 1, 3

         call slxyzs(il,im,idir,jl,jm,lms)
         dtefield%Lij(1,klmn,itypat,idir)=intg*real(lms)
         dtefield%Lij(2,klmn,itypat,idir)=intg*aimag(lms)
         dtefield%Lijr3(1,klmn,itypat,idir)=intgr3*real(lms)
         dtefield%Lijr3(2,klmn,itypat,idir)=intgr3*aimag(lms)

       end do

     end do ! end loop over ilmn
   end do ! end loop over jlmn


   ABI_DEALLOCATE(ff)

 end do ! end loop over atom types

 dtefield%has_Lij = 2

 end subroutine Lij
!!***

