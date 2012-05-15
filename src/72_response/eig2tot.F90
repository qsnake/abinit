!{\src2tex{textfont=tt}}
!!****f* ABINIT/eig2tot
!! NAME
!! eig2tot
!!
!! FUNCTION
!! This routine calculates the second-order eigenvalues.
!! The output eig2nkq is this quantity for the input k points.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (PB, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  bdeigrf = number of bands for which to calculate the second-order eigenvalues
!!  clflg(3,mpert)= Array on calculated perturbations for eig2rf
!!  dim_eig2nkq=1 if eig2nkq is to be computed
!!  cg1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol,3,mpert)= first-order wf in G space for each perturbation 
!!    The wavefunction is orthogonal to the active space.
!!  gh0c1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol,3,mpert)= matrix containing the vector:  <G|H(0)|psi(1)>, for each perturbation
!!  gh1c_pert(2,mpw1*nspinor*mband*mk1mem*nsppol,3,mpert))= matrix containing the vector:  <G|H(1)|n,k>, for each perturbation
!!    The wavefunction is orthogonal to the active space. 
!!  eigbrd(2,mband*nsppol,nkpt,3,natom,3,natom)=broadening factors for the electronic eigenvalues (optional)
!!  eigen0(nkpt_rbz*mband*nsppol)= 0-order eigenvalues at all K-points: <k,n'|H(0)|k,n'> (hartree)
!!  eigenq(nkpt_rbz*mband*nsppol)= 0-order eigenvalues at all shifted K-points: <k+Q,n'|H(0)|k+Q,n'> (hartree)
!!  eigen1(nkpt_rbz*2*nsppol*mband**2,3,mpert)= matrix of first-order: <k+Q,n'|H(1)|k,n> (hartree) (calculated in cgwf3)
!!  eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom*dim_eig2nkq)=second derivatives of the electronic eigenvalues
!!  elph2_imagden=imaginary part of the denominator of the sum-over-state expression for the electronic eigenenergy shift 
!!              due to second-order electron-phonon interation
!!  ieig2rf= integer for calculation type
!!  indsym(4,nsym,natom)= indirect indexing array for atom labels (not used yet, but will be used with symmetries)
!!  istwfk_pert(nkpt_rbz,3,mpert)= integer for choice of storage of wavefunction at each k point for each perturbation
!!  mband= maximum number of bands
!!  mk1mem= maximum number of k points which can fit in memory (RF data)  ; 0 if use disk
!!  mpert= maximum number of perturbations
!!  natom= number of atoms in the unit cell
!!  npert= number of phonon perturbations, without taking into account directions: natom 
!!  nsym= number of symmetries (not used yet)
!!  mpi_enreg= informations about MPI parallelization
!!  mpw1= maximum number of planewaves used to represent first-order wavefunctions
!!  nkpt_rbz= number of k-points for each perturbation
!!  npwar1(nkpt_rbz,mpert)= number of planewaves at k-point for first-order
!!  nspinor= number of spinorial components of the wavefunctions
!!  nsppol= 1 for unpolarized, 2 for spin-polarized
!!  smdelta= integer controling the calculation of electron lifetimes
!!  symq(4,2,nsym)= 1 if symmetry preserves present qpoint. From symq3 (not used yet)
!!  symrec(3,3,nsym)= 3x3 matrices of the group symmetries (reciprocal space) (not used yet)
!!  symrel(3,3,nsym)= array containing the symmetries in real space (not used yet)
!!  timrev= 1 if time-reversal preserves the q wavevector; 0 otherwise (not in use yet)
!!
!! OUTPUT
!!  eig2nkq(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)= diagonal part of the second-order eigenvalues: E^{(2),diag}_{k,q,j}
!!  eigbrd(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)= OPTIONAL, array containing the the electron lifetimes
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      distrb2,dotprod_g,leave_new,smeared_delta,timab,wrtout,xcomm_world
!!      xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine eig2tot(bdeigrf,clflg,cg1_pert,dim_eig2nkq,dim_eig2rf,eigen0,eigenq,eigen1,eig2nkq,elph2_imagden,esmear,&
&  gh0c1_pert,gh1c_pert,ieig2rf,istwfk_pert,mband,mk1mem,mpert,npert,mpi_enreg,mpw1,&
&  nkpt_rbz,npwar1,nspinor,nsppol,smdelta,eigbrd)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eig2tot'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_spacepar
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bdeigrf,dim_eig2nkq,dim_eig2rf,ieig2rf,mband,mk1mem,mpert,mpw1,nkpt_rbz
 integer,intent(in) :: npert,nspinor,nsppol,smdelta
 real(dp),intent(in) :: elph2_imagden,esmear
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: clflg(3,mpert)
 integer,intent(in) :: istwfk_pert(nkpt_rbz,3,mpert)
 integer,intent(in) :: npwar1(nkpt_rbz,mpert)
 real(dp),intent(in) :: cg1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
 real(dp),intent(in) :: gh0c1_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
 real(dp),intent(in) :: gh1c_pert(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf,3,mpert)
 real(dp),intent(in) :: eigen0(nkpt_rbz*mband*nsppol)
 real(dp),intent(in) :: eigen1(nkpt_rbz*2*nsppol*mband**2,3,mpert)
 real(dp),intent(in) :: eigenq(nkpt_rbz*mband*nsppol)
 real(dp),intent(out) :: eig2nkq(2,mband*nsppol,nkpt_rbz,3,npert,3,npert*dim_eig2nkq)
 real(dp),intent(out),optional :: eigbrd(2,mband*nsppol,nkpt_rbz,3,npert,3,npert)

!Local variables-------------------------------
!tolerance for non degenerated levels
!scalars
 integer :: band2tot_index,band_index,bandtot_index,iband,icg2,idir1,idir2
 integer :: ikpt,ipert1,ipert2,isppol,istwf_k,jband,npw1_k
!integer :: ipw
 integer :: master,me,spaceworld,ierr
!real(dp),parameter :: etol=1.0d-3
 real(dp),parameter :: etol=1.0d-6
!real(dp),parameter :: etol=zero 
 real(dp) :: ar,ai,deltae,den,dot2i,dot2r,dot3i,dot3r,doti,dotr,eig1_i1,eig1_i2
 real(dp) :: eig1_r1,eig1_r2,eig2_diai
 real(dp) :: eig2_diar,eigbrd_i,eigbrd_r
 character(len=500) :: message
 logical :: test_do_band
!arrays
 integer :: blk1flg(3,mpert,3,mpert)
 integer, allocatable :: nband_rbz(:),icg2_rbz(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwavef(:,:),cwavef2(:,:)
 real(dp) :: eigen(mband*nsppol),eigen_prime(mband*nsppol)
 real(dp),allocatable :: gh(:,:),gh1(:,:),ghc(:,:)
 real(dp),allocatable :: smdfun(:,:)

! *********************************************************************

!Init parallelism
 call xcomm_world(mpi_enreg,spaceworld)
 call xme_init(mpi_enreg,me)
 master =0

!DEBUG
!write(std_out,*)' eig2tot : enter '
!write(std_out,*)' mpw1=',mpw1
!write(std_out,*)' mband=',mband
!write(std_out,*)' nsppol=',nsppol
!write(std_out,*)' nkpt_rbz=',nkpt_rbz
!write(std_out,*)' npert=',npert
!ENDDEBUG

 call timab(148,1,tsec)

 if(nsppol==2)then
   write(message, '(4a)' )ch10,&
&   ' eig2tot : ERROR -',ch10,&
&   '  nsppol=2 is not yet treated. Sorry for this ... '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if        

 band2tot_index =0
 bandtot_index=0
 band_index=0

 if(ieig2rf == 1 .or. ieig2rf == 2) then
   eig2nkq(:,:,:,:,:,:,:) = zero
 end if
 if(present(eigbrd))then
   eigbrd(:,:,:,:,:,:,:) = zero
 end if
 blk1flg(:,:,:,:) = 0

 if(mpi_enreg%paral_compil_kpt==1) then
   ABI_ALLOCATE(mpi_enreg%proc_distrb,(nkpt_rbz,mband,nsppol))
   ABI_ALLOCATE(nband_rbz,(nkpt_rbz))
!  Assume the number of bands is the same for all k points.
   nband_rbz(:)=mband
   call distrb2(mband,nband_rbz,nkpt_rbz,nsppol,mpi_enreg)
 end if

 icg2=0
 ipert1=1 ; isppol=1 ! Suppose that the situation is the same for all perturbations
 ABI_ALLOCATE(icg2_rbz,(nkpt_rbz))
 do ikpt=1,nkpt_rbz
   icg2_rbz(ikpt)=icg2
   if(mpi_enreg%paral_compil_kpt==1)then
     if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:mband,isppol)-me))/=0) then
       cycle
     end if
   end if
   icg2 = icg2 + npwar1(ikpt,ipert1)*nspinor*mband !does not work with isppol
 end do

 do isppol=1,nsppol
   do ikpt =1,nkpt_rbz

     if(mpi_enreg%paral_compil_kpt==1)then
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:mband,isppol)-me))/=0) then
         band2tot_index = band2tot_index + 2*mband**2
         bandtot_index = bandtot_index + mband
         cycle
       end if
     end if

     if(smdelta >0) then   !broadening
       if(.not.allocated(smdfun))  then
         ABI_ALLOCATE(smdfun,(mband,mband))
       end if
       smdfun(:,:) = zero
       do iband=1,mband
         eigen(iband) = eigen0(iband+bandtot_index)
         eigen_prime(iband) =eigenq(iband+bandtot_index)
       end do
       if(esmear>tol6) call smeared_delta(eigen,eigen_prime,esmear,mband,smdelta,smdfun)
     end if
     icg2=icg2_rbz(ikpt)

     ipert1=1 ! Suppose all perturbations lead to the same number of planewaves
     npw1_k = npwar1(ikpt,ipert1)
     ABI_ALLOCATE(cwavef,(2,npw1_k*nspinor))
     ABI_ALLOCATE(cwavef2,(2,npw1_k*nspinor))
     ABI_ALLOCATE(gh,(2,npw1_k*nspinor))
     ABI_ALLOCATE(gh1,(2,npw1_k*nspinor))
     ABI_ALLOCATE(ghc,(2,npw1_k*nspinor))

     do iband=1,bdeigrf

!      If the k point and band belong to me, compute the contribution
       test_do_band=.true.
       if(mpi_enreg%paral_compil_kpt==1)then
         if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me)test_do_band=.false.
       end if
       if(test_do_band)then

         do ipert1=1,npert

           do idir1=1,3
             if(clflg(idir1,ipert1)==0)cycle
             istwf_k = istwfk_pert(ikpt,idir1,ipert1)

             do ipert2=1,npert
               do idir2=1,3
                 if(clflg(idir2,ipert2)==0)cycle
                 blk1flg(idir1,ipert1,idir2,ipert2)=1

                 eig2_diar = zero ; eig2_diai = zero ; eigbrd_r = zero ; eigbrd_i = zero

                 do jband=1,mband
                   eig1_r1 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                   eig1_r2 = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index,idir2,ipert2)
                   eig1_i1 = eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir1,ipert1)
                   eig1_i2 = - eigen1(2*jband+(iband-1)*2*mband+band2tot_index,idir2,ipert2) !the negative sign is from the CC
                   deltae= eigenq(jband+bandtot_index)-eigen0(iband+bandtot_index)
                   ar=eig1_r1*eig1_r2-eig1_i1*eig1_i2
                   ai=eig1_r1*eig1_i2+eig1_i1*eig1_r2
                   
!                  Sum over all active space to retrieve the diagonal gauge
                   if(ieig2rf == 1 .or. ieig2rf ==2 ) then
                     if(abs(deltae)>etol) then

                       den=-one/(deltae**2+elph2_imagden**2)

!                      The following should be the most general implementation of the presence of elph2_imagden
!                      eig2_diar=eig2_diar+(ar*deltae+ai*elph2_imagden)*den
!                      eig2_diai=eig2_diai+(ai*deltae-ar*elph2_imagden)*den
!                      This gives back the implementation without elph2_imagden
!                      eig2_diar=eig2_diar+ar*deltae*den
!                      eig2_diai=eig2_diai+ai*deltae*den
!                      This is what Samuel had implemented
!                      eig2_diar=eig2_diar+ar*deltae*den
!                      eig2_diai=eig2_diai+ai*elph2_imagden*den
!                      Other possibility : throw away the broadening part, that is actually treated separately.
                       eig2_diar=eig2_diar+ar*deltae*den
                       eig2_diai=eig2_diai+ai*deltae*den

                     end if ! abs(deltae)>etol
                   end if ! ieig2rf==1 or 2

                   if(present(eigbrd))then
                     if(smdelta >0) then   !broadening
                       eigbrd_r = eigbrd_r + ar*smdfun(iband,jband)
                       eigbrd_i = eigbrd_i + ai*smdfun(iband,jband)
                     end if
                   end if

                 end do !jband

!                Add the contribution of non-active bands, if DFPT calculation
                 if(ieig2rf == 1) then

                   dotr=zero ; doti=zero
                   dot2r=zero ; dot2i=zero
                   dot3r=zero ; dot3i=zero


                   cwavef(:,:) = cg1_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir2,ipert2)
                   cwavef2(:,:)= cg1_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir1,ipert1)
                   gh1(:,:)    = gh1c_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir1,ipert1)
                   gh(:,:)     = gh1c_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir2,ipert2)
                   ghc(:,:)    = gh0c1_pert(:,1+(iband-1)*npw1_k*nspinor+icg2:iband*npw1_k*nspinor+icg2,idir1,ipert1)

!                  The first two dotprod corresponds to:  <Psi(1)|H(1)|Psi(0)> + cc.
!                  They are calculated using wavefunctions <Psi(1)| that are orthogonal to the active space.
                   call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1_k*nspinor,2,cwavef,gh1)
                   call dotprod_g(dot2r,dot2i,istwf_k,mpi_enreg,npw1_k*nspinor,2,gh,cwavef2)

!                  This dotprod corresponds to : <Psi(1)|H(0)- E(0)|Psi(1)>
!                  It is calculated using wavefunctions that are orthogonal to the active space.
!                  Should work for metals. (But adiabatic approximation is bad in this case...)
                   call dotprod_g(dot3r,dot3i,istwf_k,mpi_enreg,npw1_k*nspinor,2,cwavef,ghc)

                   eig2_diar= eig2_diar + dotr + dot2r + dot3r
                   eig2_diai= eig2_diai + doti + dot2i + dot3i

                 end if

!                Store the contribution
                 if(ieig2rf == 1 .or. ieig2rf == 2) then
                   eig2nkq(1,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eig2_diar
                   eig2nkq(2,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eig2_diai 
                 end if

                 if(present(eigbrd))then
                   if(smdelta >0) then   !broadening
                     eigbrd(1,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eigbrd_r
                     eigbrd(2,iband+band_index,ikpt,idir1,ipert1,idir2,ipert2) = eigbrd_i
                   end if
                 end if

               end do !idir2
             end do !ipert2
           end do  !idir1
         end do   !ipert1

       end if ! Selection of processor
       
     end do !iband

     ABI_DEALLOCATE(cwavef)
     ABI_DEALLOCATE(cwavef2)
     ABI_DEALLOCATE(gh)
     ABI_DEALLOCATE(gh1)
     ABI_DEALLOCATE(ghc)
     band2tot_index = band2tot_index + 2*mband**2
     bandtot_index = bandtot_index + mband

   end do    !ikpt
   band_index = band_index + mband
 end do !isppol

!Accumulate eig2nkq and/or eigbrd
 if(mpi_enreg%paral_compil_kpt==1) then
   if(ieig2rf == 1 .or. ieig2rf == 2) then
     call xsum_mpi(eig2nkq,spaceworld,ierr)
   end if
   if(present(eigbrd))then
     if(smdelta >0) then
       call xsum_mpi(eigbrd,spaceworld,ierr)
     end if
   end if
   ABI_DEALLOCATE(nband_rbz)
   ABI_DEALLOCATE(mpi_enreg%proc_distrb)
 end if 

 if(ieig2rf == 1 .or. ieig2rf == 2) then
   write(ab_out,'(a)')' Components of second-order derivatives of the electronic energy, EIGR2D.'
   write(ab_out,'(a)')' For automatic tests, printing the matrix for the first k-point, first band, first atom.'
   do idir1=1,3
     do idir2=1,3
       ar=eig2nkq(1,1,1,idir1,1,idir2,1) ; if(abs(ar)<tol10)ar=zero
       ai=eig2nkq(2,1,1,idir1,1,idir2,1) ; if(abs(ai)<tol10)ai=zero
       write (ab_out,'(4i4,2es20.10)') idir1,1,idir2,1,ar,ai 
     end do ! idir2
   end do ! idir1
 end if 

 if(present(eigbrd))then
   if(smdelta >0) then   !broadening
     write(ab_out,'(a)')' '
     write(ab_out,'(a)')' Components of second-order derivatives of the electronic energy, EIGI2D.'
     write(ab_out,'(a)')' For automatic tests, printing the matrix for the first k-point, first band, first atom.'
     do idir1=1,3
       do idir2=1,3
         ar=eigbrd(1,1,1,idir1,1,idir2,1) ; if(abs(ar)<tol10)ar=zero
         ai=eigbrd(2,1,1,idir1,1,idir2,1) ; if(abs(ai)<tol10)ai=zero
         write (ab_out,'(4i4,2es20.10)') idir1,1,idir2,1,ar,ai
       end do
     end do !nband
   end if
 end if

 if(allocated(smdfun))  then
   ABI_DEALLOCATE(smdfun)
 end if
 ABI_DEALLOCATE(icg2_rbz)

 call timab(148,2,tsec)

!DEBUG
!write(std_out,*)' eig2tot: exit'
!ENDDEBUG

end subroutine eig2tot
!!***

