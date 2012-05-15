!{\src2tex{textfont=tt}}
!!****f* ABINIT/crho
!! NAME
!! crho
!!
!! FUNCTION
!! Calculate the charge density rho on the FFT grid.
!! In case of nsppol==2 calculate rho_up and rho_down
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  mpi_enreg= datatype gathering information on parallelism, variables used
!!   |gwpara= if 2 bands are spread btw processors
!!  nbnds = number of bands.
!!  nkibz = number of k-points in the irreducible Brillouin zone.
!!  nsym = number of symmetry operations.
!!  nfftot = total number of points on the FFT grid.
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(nkibz,nbnds,nsppol) = occupation numbers for nbnds bands at nkibz irreducible k-points, for each spin
!!  rprimd(3,3)=dimensional real space primitive translations
!!  ucvol = unit cell volume.
!!  wfr(nfftot,nbnds,nkibz,nsppol) = wavefunctions on the FFT grid for nbnds bands at nkibz irreducible k-points, for each spin
!!  wtk(nkibz) = irreducible k-points weights.
!!
!! OUTPUT
!!  omegaplasma = the plasma frequency.
!!  rho(nfftot,nsppol) = the density on the FFT grid.
!!   (total in first half and spin-up in second half if nsppol=2)
!!
!! PARENTS
!!      rdm
!!
!! CHILDREN
!!      initmpi_seq,irrzg,symrhg,wrtout,xcomm_init,xmaster_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine crho(paral_kgb,ngfft,gprimd,nbnds,nkibz,nsym,symrel,tnons,symafm,&
& nfftot,nspden,nsppol,occ,omegaplasma,rho,rprimd,ucvol,wfr,wtk,mpi_enreg,my_minb,my_maxb)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'crho'
 use interfaces_14_hidewrite
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
 use interfaces_67_common
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_maxb,my_minb,nbnds,nfftot,nkibz,nspden,nsppol
 integer,intent(in) :: nsym,paral_kgb
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: omegaplasma
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),symafm(nsym)
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),occ(nkibz,nbnds,nsppol),rprimd(3,3)
 real(dp),intent(in) :: tnons(3,nsym),wtk(nkibz)
 real(dp),intent(out) :: rho(nfftot,nsppol)
 complex(gwpc),intent(in) :: wfr(nfftot,my_minb:my_maxb,nkibz,nsppol)

!Local variables ------------------------------
!scalars
 integer :: cplex,ib,ik,ir,is,master,me,n1,n2,n3,spaceComm
 real(dp) :: rhoav,rs,tnepuc
 character(len=500) :: message
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer,allocatable :: irrzon(:,:,:)
 real(dp),allocatable :: phnons(:,:,:),rho2(:),rhog(:,:)

!*************************************************************************

#if defined DEBUG_MODE
 write(message,'(a)')' crho : enter '
 call wrtout(std_out,message,'PERS')
#endif

 write(message,'(2a)')ch10,' crho: calculating charge density...'
 call wrtout(std_out,message,'COLL')
!
!Initialize some MPI related variables
 call xcomm_init(mpi_enreg,spaceComm)
 call xme_init(mpi_enreg,me)
 call xmaster_init(mpi_enreg,master)
!
!Calculate IBZ contribution to the charge density
 ABI_ALLOCATE(rho2,(nfftot))
 rho(:,:)= zero
 do is=1,nsppol
   rho2(:)= zero
!  Loop over k-points in IBZ
   do ik=1,nkibz
!    Skip the higher bands if occupation is less than tol8
!    do while ((abs(occ(ik,ib,is))>tol8).and.(ib<=nbnds))
     do ib=1,nbnds
!      if (mpi_enreg%gwpara==2) then
!      if (mpi_enreg%proc_distrb(ik,ib,is)/=me) cycle
!      end if
       if (abs(occ(ik,ib,is))<tol8) cycle
       do ir=1,nfftot
         rho2(ir)= rho2(ir) + occ(ik,ib,is)*conjg(wfr(ir,ib,ik,is))*wfr(ir,ib,ik,is)*wtk(ik)/SUM(wtk)/ucvol
       end do !ir
     end do !ib
   end do !ik
!  we could sum rho outside the loop over is
!  if (mpi_enreg%gwpara==2) then
!  call xsum_mpi(rho2,spaceComm,ier)
!  end if

   rho(:,is) = rho2(:)

 end do !is
!
!Store the total charge in the first half
!if (nsppol==2) then
!rho2(:) = rho(:,1)
!rho(:,1)= rho(:,1)+rho(:,2)
!rho(:,2)= rho2(:)
!end if

!NEW symmetrization in G space implementing also the AFM case.
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 ABI_ALLOCATE(irrzon,(nfftot,2,(nspden/nsppol)-3*(nspden/4)))
 ABI_ALLOCATE(phnons,(2,nfftot,(nspden/nsppol)-3*(nspden/4)))

 call irrzg(irrzon,nspden,nsppol,nsym,n1,n2,n3,&
& phnons,symafm,symrel,tnons)

!* Fake MPI_type for sequential part
 call initmpi_seq(MPI_enreg_seq)

 cplex=1
 ABI_ALLOCATE(rhog,(2,cplex*nfftot))

 call symrhg(cplex,gprimd,irrzon,MPI_enreg_seq,nfftot,nfftot,ngfft,nspden,nsppol,&
& nsym,paral_kgb,phnons,rhog,rho,rprimd,symafm,symrel)

 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(phnons)
 ABI_DEALLOCATE(irrzon)
!
!Calculate total number of electrons as a check
 tnepuc=zero
 do ir=1,nfftot
   tnepuc=tnepuc+rho(ir,1)
 end do
 tnepuc=tnepuc*ucvol/nfftot ; rhoav=tnepuc/ucvol ; rs=(three/(four_pi*rhoav))**third

 write(message,'(a,f9.4)')' total number of electrons per unit cell = ',tnepuc
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 write(message,'(a,f9.6)')' average of density, n = ',rhoav
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 write(message,'(a,f9.4)')' r_s = ',rs
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 omegaplasma= sqrt(four_pi*rhoav)
 write(message,'(a,f9.4,2a)')' omega_plasma = ',omegaplasma*Ha_eV,' [eV]',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 ABI_DEALLOCATE(rho2)

#if defined DEBUG_MODE
 write(message,'(a)')' crho : exit '
 call wrtout(std_out,message,'PERS')
#endif

end subroutine crho
!!***
