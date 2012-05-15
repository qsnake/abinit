!{\src2tex{textfont=tt}}
!!****f* ABINIT/jvec_to_B
!! NAME
!! jvec_to_B
!!
!! FUNCTION
!! compute B field and hence chemical shieldings due to a vector current field at atomic points
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! gcart(ngfft(1),ngfft(2),ngfft(3),3), G vectors in cartesian space
!! jvec(Bdir=1..3,jvec_i=1..3,nfft), current vector field on the grid for each of three external B-field directions 
!! natom, number of atoms in unit cell
!! nfft,ngfft(18), number of FFT points and details of FFT
!! paral_kgb, flag about parallel calculation used in call to fourdp
!! rprimd(3,3), conversion from crystal coordinates to cartesian coordinates
!!
!! OUTPUT
!! cshield(3,3,natom) 3x3 shielding tensor at each atom site due to this current field
!!
!! SIDE EFFECTS
!! xred(3,natom), location of atoms in crystal coordinates. It appears with intent(inout) here
!!                because routine xredxcart requires the ability to change xred, although we do not
!!                use that facility in the present call.
!!
!! NOTES
!! The vector field jvec, itself induced by an external magnetic field, induces a local magnetic field
!! at each atomic site. The induced field is related to the external field by $B_{ind} = -\sigma B_{ext}$,
!! where $\sigma$ is the chemical shielding tensor. Thus for an external field of unit strength, the
!! internal field is numerically equivalent to the shielding. Here the vector field induced by each direction
!! of the external field is passed as input and then the Biot-Savart Law in reciprocal space is applied
!! to calculate the induced field. The specific equation is
!! $\mathbf{B(R)}=\frac{2i}{c}\sum_{\mathbf{G}}e^{2\pi i\mathbf{G\cdot R}}\frac{\mathbf{j(G)\wedge G}{G^2}$
!! (see Yates, Pickard, Mauri, PRB 76, 024401 (2007), Eq. 45, but note that our form has one factor of
!! $2\pi$ less, because we define reciprocal space by $e^{2\pi i \mathbf{k\cdot r}$, NOT by
!! $e^{i\mathbf{k\cdot r}}$ like they do). 
!!
!! PARENTS
!!
!! CHILDREN
!!      acrossb,fourdp,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine jvec_to_B(cshield,gcart,jvec,natom,nfft,ngfft,paral_kgb,rprimd,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'jvec_to_B'
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,paral_kgb
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gcart(ngfft(1),ngfft(2),ngfft(3),3),jvec(3,3,nfft)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: xred(3,natom)
 real(dp),intent(out) :: cshield(3,3,natom)

!Local variables-------------------------------
!scalars
 integer :: cplex,iatom,idir,igfft1,igfft2,igfft3,index,isign,jdir
 integer :: tim_fourdp
 real(dp) :: cph,phase,sph,trace
 type(MPI_type) :: mpi_enreg
!arrays
 real(dp) :: BG(2,3),Bind(2,3,3,natom),fofg(2,nfft),fofr(nfft),gvec(3)
 real(dp) :: jvecg(2,3,nfft),jxg(3),ratom(3),xcart(3,natom)
!no_abirules
!

! ************************************************************************

!DEBUG
!write(std_out,*)' jvec_to_B : enter'
!ENDDEBUG

 call xredxcart(natom,1,rprimd,xcart,xred) ! get atomic locations in cartesian coords

 Bind(:,:,:,:) = zero
 do idir = 1, 3 ! loop over directions of external field
   do jdir = 1, 3 ! loop over components of the vector field
     tim_fourdp = 0 ! timing code, not using
     isign = -1 ! FT from R to G
     cplex = 1 ! fofr is real
!    extract the jdir component of the current vector field, produced from the idir field direction
     fofr(:) = jvec(idir,jdir,:)
     call fourdp(cplex,fofg,fofr,isign,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp) ! construct current in G space
     jvecg(1,jdir,:)=fofg(1,:); jvecg(2,jdir,:)=fofg(2,:);
   end do ! end loop on real space current components
!  now have j(G), the vector current field expressed in G space, for external B in the idir direction

   do igfft1 = 1, ngfft(1) ! looping over points in the G grid
     do igfft2 = 1, ngfft(2)
       do igfft3 = 1, ngfft(3)
         index = (igfft3-1)*ngfft(2)*ngfft(1) + (igfft2-1)*ngfft(1) + igfft1
         gvec(:) = gcart(igfft1,igfft2,igfft3,:) ! gvec is the vector in G space
         trace = dot_product(gvec,gvec)
         if (trace > zero) then ! avoid G = 0 point
           call acrossb(jvecg(2,:,index),gvec,jxg)
           BG(1,:) = -2.0*jxg(:)/(trace*Sp_Lt)
           call acrossb(jvecg(1,:,index),gvec,jxg)
           BG(2,:) = 2.0*jxg(:)/(trace*Sp_Lt)
           do iatom = 1, natom ! sum over atoms in unit cell
             ratom(:) = xcart(:,iatom) ! extract location of atom iatom
             phase = two_pi*dot_product(gvec,ratom) ! argument of $e^{2\pi i G\cdot R}$
             cph = cos(phase)
             sph = sin(phase)
             do jdir = 1, 3 ! loop over directions of induced field
               Bind(1,jdir,idir,iatom) = Bind(1,jdir,idir,iatom) + BG(1,jdir)*cph-BG(2,jdir)*sph
               Bind(2,jdir,idir,iatom) = Bind(2,jdir,idir,iatom) + BG(1,jdir)*sph+BG(2,jdir)*cph
             end do ! end loop over induced field directions
           end do ! end loop over atoms in cell
         end if ! end statement avoiding G = 0 point
       end do ! end loop over ngfft(3)
     end do ! end loop over ngfft(2)
   end do ! end loop over ngfft(1)
 end do ! end loop over directions of external field

 cshield(:,:,:) = -1.0D6*Bind(1,:,:,:) ! Bind = -sigma_cs * B, given here in ppm


!DEBUG
!write(std_out,*)' jvec_to_B : exit '
!stop
!ENDDEBUG

end subroutine jvec_to_B
!!***
