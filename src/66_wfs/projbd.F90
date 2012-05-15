!{\src2tex{textfont=tt}}
!!****f* ABINIT/projbd
!!
!! NAME
!! projbd
!!
!! FUNCTION
!! Project out vector "direc" onto the bands contained in "cg".
!! if useoverlap==0
!!  New direc=direc-$sum_{j/=i} { <cg_{j}|direc>.|cg_{j}> }$
!! if useoverlap==1 (use of overlap matrix S)
!!  New direc=direc-$sum_{j/=i} { <cg_{j}|S|direc>.|cg_{j}> }$
!! (index i can be set to -1 to sum over all bands)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg(2,mcg)=wavefunction coefficients for ALL bands
!!  iband0=which particular band we are interested in
!!         ("i" in the above formula)
!!         Can be set to -1 to sum over all bands...
!!  icg=shift to be given to the location of the data in cg
!!  iscg=shift to be given to the location of the data in cg
!!  istwf_k=option parameter that describes the storage of wfs
!!  mcg=maximum size of second dimension of cg
!!  mpi_enreg=informations about MPI parallelization
!!  mscg=maximum size of second dimension of scg
!!  nband=number of bands
!!  npw=number of planewaves
!!  nspinor=number of spinorial components (on current proc)
!!  ortalg= choice of algorithm for the projection
!!   Note : negative values are used outside to select whether
!!   the perpendicular projection is done prior preconditioning.
!!  printopt= if 1, print intermediate dot products
!!  scg(2,mscg*useoverlap)=<G|S|band> for ALL bands,
!!                        where S is an overlap matrix
!!  scprod_io=0 if scprod array has to be computed; ! if it is input (already in memory)
!!  tim_projbd=timing code of the calling subroutine(can be set to 0 if not attributed)
!!  useoverlap=describe the overlap of wavefunctions:
!!               0: no overlap (S=Identity_matrix)
!!               1: wavefunctions are overlapping
!!
!! SIDE EFFECTS
!!  direc(2,npw)= input: vector to be orthogonalised with respect to cg (and S)
!!                output: vector that has been orthogonalized wrt cg (and S)
!!
!!  scprod(2,nband)=scalar_product
!!        if useoverlap==0: scalar_product_i=$<cg_{j}|direc_{i}>$
!!        if useoverlap==1: scalar_product_i=$<cg_{j}|S|direc_{i}>$
!!    if scprod_io=0, scprod is output
!!    if scprod_io=1, scprod is input
!!
!! NOTES
!!  XG030513: MPIWF Might have to be recoded for efficient paralellism
!!  Note for PAW: ref.= PRB 73, 235101 (2006), equations (71) and (72):
!!                in normal use, projbd applies P_c projector
!!                if cg and scg are inverted, projbd applies P_c+ projector
!! PARENTS
!!      cgwf,cgwf3,getdc1,nstpaw3
!!
!! CHILDREN
!!      timab,wrtout,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine projbd(cg,direc,iband0,icg,iscg,istwf_k,mcg,mpi_enreg,mscg,nband,&
&                 npw,nspinor,ortalg,printopt,scg,scprod,scprod_io,tim_projbd,useoverlap)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'projbd'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: iband0,icg,iscg,istwf_k,mcg,mscg,nband,npw,nspinor
 integer,intent(in) :: ortalg,printopt,scprod_io,tim_projbd,useoverlap
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: cg(2,mcg),scg(2,mscg*useoverlap)
 real(dp),intent(inout) :: direc(2,npw*nspinor)
 real(dp),intent(inout) :: scprod(2,nband)

!Local variables-------------------------------
!scalars
 integer :: iband,iband2,iblock,ierr,index1,index2,index3,index4,ipw,ipw1,isp
 integer :: nbandm,npw_sp,old_paral_level,spaceComm
 real(dp) :: ai,ai2,ai3,ai4,ar,ar2,ar3,ar4,cg_im,cg_im2,cg_im3,cg_im4,cg_re
 real(dp) :: cg_re2,cg_re3,cg_re4,direc_im,direc_re,scg_im,scg_im2,scg_im3
 real(dp) :: scg_im4,scg_re,scg_re2,scg_re3,scg_re4
 character(len=500) :: message
!arrays
 real(dp) :: buffer2(2),tsec(2)
 real(dp),allocatable :: atab(:)

! *************************************************************************
!

!DEBUG
!write(std_out,*)' projbd : enter '
!stop
!ENDDEBUG

 old_paral_level=mpi_enreg%paral_level
 mpi_enreg%paral_level=3
 call timab(210+tim_projbd,1,tsec)

 npw_sp=npw*nspinor

!Here the common coding
 if(ortalg==0 .or. ortalg==1 .or. ortalg==-1 )then

   nbandm=nband

   if(istwf_k==1)then

     if (scprod_io==0) then
       ABI_ALLOCATE(atab,(2*nbandm))
       atab=zero
       if (useoverlap==1) then
         do iband=1,nbandm
           ar=zero ; ai=zero
           index1=npw_sp*(iband-1)+iscg
!          $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ai,ar) &
!          $OMP&SHARED(scg,direc,index1,npw_sp)
           do ipw=1,npw_sp
             ar=ar+scg(1,index1+ipw)*direc(1,ipw)+scg(2,index1+ipw)*direc(2,ipw)
             ai=ai-scg(2,index1+ipw)*direc(1,ipw)+scg(1,index1+ipw)*direc(2,ipw)
           end do
!          $OMP END PARALLEL DO
           iband2=2*iband
           atab(iband2-1)=ar;atab(iband2)=ai
         end do ! Loop on iband
       else
         do iband=1,nbandm
           ar=zero ; ai=zero
           index1 =npw_sp*(iband-1)+icg
!          $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ai,ar) &
!          $OMP&SHARED(cg,direc,index1,npw_sp)
           do ipw=1,npw_sp
             ar=ar+cg(1,index1+ipw)*direc(1,ipw)+cg(2,index1+ipw)*direc(2,ipw)
             ai=ai-cg(2,index1+ipw)*direc(1,ipw)+cg(1,index1+ipw)*direc(2,ipw)
           end do
!          $OMP END PARALLEL DO
           iband2=2*iband
           atab(iband2-1)=ar;atab(iband2)=ai
         end do ! Loop on iband
       end if ! useoverlap
       if (mpi_enreg%paral_compil_fft==1) then
         call xcomm_init(mpi_enreg,spaceComm)
         call timab(48,1,tsec)
         call xsum_mpi(atab,spaceComm,ierr)
         call timab(48,2,tsec)
       end if
       scprod(:,1:nbandm)=reshape(atab(1:2*nbandm),(/2,nbandm/))
       ABI_DEALLOCATE(atab)
     end if

     do iband=1,nbandm
       if (iband==iband0) cycle
       index1=npw_sp*(iband-1)+icg
       ar=scprod(1,iband);ai=scprod(2,iband)
!      $OMP PARALLEL DO PRIVATE(ipw,cg_re,cg_im) &
!      $OMP&SHARED(ar,ai,cg,direc,index,npw_sp)
       do ipw=1,npw_sp
         cg_re=cg(1,index1+ipw) ; cg_im=cg(2,index1+ipw)
         direc(1,ipw)=direc(1,ipw)-ar*cg_re+ai*cg_im
         direc(2,ipw)=direc(2,ipw)-ar*cg_im-ai*cg_re
       end do
!      $OMP END PARALLEL DO
       if(printopt==1)then
         write(message,'(a,i3,2f14.6)') &
&         'projbd : called from cgwf ; iband,ar,ai=',iband,ar,ai
         call wrtout(std_out,message,'PERS')
       end if
     end do ! Loop on iband

   else if(istwf_k>=2)then

     if (scprod_io==0) then
       ABI_ALLOCATE(atab,(nbandm))
       atab=zero
       if (useoverlap==1) then
         do iband=1,nbandm
           index1=npw_sp*(iband-1)+iscg
           if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
             ar=half*scg(1,index1+1)*direc(1,1) ; ipw1=2
           else
             ar=zero ; ipw1=1
           end if
!          $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar) &
!          $OMP&SHARED(scg,direc,index1,ipw1,npw,nspinor)
           do isp=1,nspinor
             do ipw=ipw1+(isp-1)*npw,npw*isp
               ar=ar+scg(1,index1+ipw)*direc(1,ipw)+scg(2,index1+ipw)*direc(2,ipw)
             end do
           end do
!          $OMP END PARALLEL DO
           atab(iband)=two*ar
         end do ! Loop on iband
       else
         do iband=1,nbandm
           index1=npw_sp*(iband-1)+icg
           if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
             ar=half*cg(1,index1+1)*direc(1,1) ; ipw1=2
           else
             ar=zero ; ipw1=1
           end if
!          $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar) &
!          $OMP&SHARED(cg,direc,index1,ipw1,npw,nspinor)
           do isp=1,nspinor
             do ipw=ipw1+(isp-1)*npw,npw*isp
               ar=ar+cg(1,index1+ipw)*direc(1,ipw)+cg(2,index1+ipw)*direc(2,ipw)
             end do
           end do
!          $OMP END PARALLEL DO
           atab(iband)=two*ar
         end do ! Loop on iband
       end if ! useoverlap
       if (mpi_enreg%paral_compil_fft==1) then
         call xcomm_init(mpi_enreg,spaceComm)
         call timab(48,1,tsec)
         call xsum_mpi(atab,spaceComm,ierr)
         call timab(48,2,tsec)
       end if
       scprod(1,1:nbandm)=atab(1:nbandm)
       scprod(2,1:nbandm)=zero
       ABI_DEALLOCATE(atab)
     end if

     do iband=1,nbandm
       if (iband==iband0) cycle
       index1=npw_sp*(iband-1)+icg
       ar=scprod(1,iband)
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(ar,cg,direc,index1,npw_sp)
       do ipw=1,npw_sp
         direc(1,ipw)=direc(1,ipw)-ar*cg(1,index1+ipw)
         direc(2,ipw)=direc(2,ipw)-ar*cg(2,index1+ipw)
       end do
!      $OMP END PARALLEL DO
       if(printopt==1)then
         write(message,'(a,i3,f14.6)') &
&         'projbd : called from cgwf ; iband,ar=',iband,ar
         call wrtout(std_out,message,'PERS')
       end if
     end do ! Loop on iband

   end if ! Test on istwf_k

!  Here better use of the registers
 else

   if (ortalg==2 .or. ortalg==-2)then

     iblock=2;nbandm=(nband/iblock)*iblock

     if(nband>=iblock)then

       if(istwf_k==1)then

         if (scprod_io==0) then
           ABI_ALLOCATE(atab,(2*nbandm))
           atab=zero
           if (useoverlap==1) then
             do iband=1,nbandm,iblock
               ar=zero  ; ai=zero
               ar2=zero ; ai2=zero
               index1=npw_sp*(iband-1)+iscg
               index2=npw_sp* iband   +iscg
!              $OMP PARALLEL DO PRIVATE(scg_re,scg_re2,scg_im,scg_im2) &
!              $OMP&PRIVATE(direc_re,direc_im,ipw) &
!              $OMP&REDUCTION(+ :ai,ai2,ar,ar2) &
!              $OMP&SHARED(scg,direc,index1,index2,npw_sp)
               do ipw=1,npw_sp
                 direc_re=direc(1,ipw)     ; direc_im=direc(2,ipw)
                 scg_re =scg(1,index1+ipw) ; scg_im =scg(2,index1+ipw)
                 scg_re2=scg(1,index2+ipw) ; scg_im2=scg(2,index2+ipw)
                 ar=ar+scg_re*direc_re+scg_im*direc_im
                 ai=ai-scg_im*direc_re+scg_re*direc_im
                 ar2=ar2+scg_re2*direc_re+scg_im2*direc_im
                 ai2=ai2-scg_im2*direc_re+scg_re2*direc_im
               end do
!              $OMP END PARALLEL DO
               iband2=2*iband
               atab(iband2-1)=ar ;atab(iband2  )=ai
               atab(iband2+1)=ar2;atab(iband2+2)=ai2
             end do ! Loop on iband
           else
             do iband=1,nbandm,iblock
               ar=zero  ; ai=zero
               ar2=zero ; ai2=zero
               index1=npw_sp*(iband-1)+icg
               index2=npw_sp* iband   +icg
!              $OMP PARALLEL DO PRIVATE(cg_re,cg_re2,cg_im,cg_im2) &
!              $OMP&PRIVATE(direc_re,direc_im,ipw) &
!              $OMP&REDUCTION(+ :ai,ai2,ar,ar2) &
!              $OMP&SHARED(cg,direc,index1,index2,npw_sp)
               do ipw=1,npw_sp
                 direc_re=direc(1,ipw)   ; direc_im=direc(2,ipw)
                 cg_re=cg(1,index1+ipw)  ; cg_im=cg(2,index1+ipw)
                 cg_re2=cg(1,index2+ipw) ; cg_im2=cg(2,index2+ipw)
                 ar=ar+cg_re*direc_re+cg_im*direc_im
                 ai=ai-cg_im*direc_re+cg_re*direc_im
                 ar2=ar2+cg_re2*direc_re+cg_im2*direc_im
                 ai2=ai2-cg_im2*direc_re+cg_re2*direc_im
               end do
!              $OMP END PARALLEL DO
               iband2=2*iband
               atab(iband2-1)=ar ;atab(iband2  )=ai
               atab(iband2+1)=ar2;atab(iband2+2)=ai2
             end do ! Loop on iband
           end if ! useoverlap
           if (mpi_enreg%paral_compil_fft==1) then
             call xcomm_init(mpi_enreg,spaceComm)
             call timab(48,1,tsec)
             call xsum_mpi(atab,spaceComm,ierr)
             call timab(48,2,tsec)
           end if
           scprod(:,1:nbandm)=reshape(atab(1:2*nbandm),(/2,nbandm/))
           ABI_DEALLOCATE(atab)
         end if

         do iband=1,nbandm,iblock
           index1=npw_sp*(iband-1)+icg
           index2=npw_sp*iband+icg
           if (iband/=iband0) then
             ar=scprod(1,iband);ai=scprod(2,iband)
           else
             ar=zero;ai=zero
           end if
           if (iband+1/=iband0) then
             ar2=scprod(1,iband+1);ai2=scprod(2,iband+1)
           else
             ar2=zero;ai2=zero
           end if
!          $OMP PARALLEL DO PRIVATE(ipw,cg_re,cg_re2,cg_im,cg_im2) &
!          $OMP&SHARED(ar,ai,ar2,ai2,cg,direc,index1,index2,npw_sp)
           do ipw=1,npw_sp
             cg_re =cg(1,index1+ipw) ; cg_im =cg(2,index1+ipw)
             cg_re2=cg(1,index2+ipw) ; cg_im2=cg(2,index2+ipw)
             direc(1,ipw)=direc(1,ipw)-ar*cg_re+ai*cg_im-ar2*cg_re2+ai2*cg_im2
             direc(2,ipw)=direc(2,ipw)-ar*cg_im-ai*cg_re-ar2*cg_im2-ai2*cg_re2
           end do
!          $OMP END PARALLEL DO
         end do ! Loop on iband

       else if(istwf_k>=2)then

         if (scprod_io==0) then
           ABI_ALLOCATE(atab,(nbandm))
           atab=zero
           if (useoverlap==1) then
             do iband=1,nbandm,iblock
               index1=npw_sp*(iband-1)+iscg
               index2=npw_sp* iband   +iscg
               if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
                 ar =half*scg(1,index1 +1)*direc(1,1)
                 ar2=half*scg(1,index2+1)*direc(1,1)
                 ipw1=2
               else
                 ar=zero ; ar2=zero ; ipw1=1
               end if
!              $OMP PARALLEL DO PRIVATE(scg_re,scg_re2,scg_im,scg_im2) &
!              $OMP&PRIVATE(direc_re,direc_im,ipw) &
!              $OMP&REDUCTION(+ :ar,ar2) &
!              $OMP&SHARED(scg,direc,index1,index2,ipw1,npw,nspinor)
               do isp=1,nspinor
                 do ipw=ipw1+(isp-1)*npw,npw*isp
                   direc_re=direc(1,ipw)     ; direc_im=direc(2,ipw)
                   scg_re =scg(1,index1+ipw) ; scg_im =scg(2,index1+ipw)
                   scg_re2=scg(1,index2+ipw) ; scg_im2=scg(2,index2+ipw)
                   ar =ar +scg_re *direc_re+scg_im *direc_im
                   ar2=ar2+scg_re2*direc_re+scg_im2*direc_im
                 end do
               end do
!              $OMP END PARALLEL DO
               atab(iband  )=two*ar
               atab(iband+1)=two*ar2
             end do ! Loop on iband
           else
             do iband=1,nbandm,iblock
               index1=npw_sp*(iband-1)+icg
               index2=npw_sp* iband   +icg
               if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
                 ar =half*cg(1,index1+1)*direc(1,1)
                 ar2=half*cg(1,index2+1)*direc(1,1)
                 ipw1=2
               else
                 ar=zero ; ar2=zero ; ipw1=1
               end if
!              $OMP PARALLEL DO PRIVATE(cg_re,cg_re2,cg_im,cg_im2) &
!              $OMP&PRIVATE(direc_re,direc_im,ipw) &
!              $OMP&REDUCTION(+ :ar,ar2) &
!              $OMP&SHARED(cg,direc,index1,index2,ipw1,npw,nspinor)
               do isp=1,nspinor
                 do ipw=ipw1+(isp-1)*npw,npw*isp
                   direc_re=direc(1,ipw)   ; direc_im=direc(2,ipw)
                   cg_re=cg(1,index1+ipw)  ; cg_im=cg(2,index1+ipw)
                   cg_re2=cg(1,index2+ipw) ; cg_im2=cg(2,index2+ipw)
                   ar =ar+cg_re *direc_re +cg_im *direc_im
                   ar2=ar2+cg_re2*direc_re+cg_im2*direc_im
                 end do
               end do
!              $OMP END PARALLEL DO
               atab(iband  )=two*ar
               atab(iband+1)=two*ar2
             end do ! Loop on iband
           end if ! useoverlap
           if (mpi_enreg%paral_compil_fft==1) then
             call xcomm_init(mpi_enreg,spaceComm)
             call timab(48,1,tsec)
             call xsum_mpi(atab,spaceComm,ierr)
             call timab(48,2,tsec)
           end if
           scprod(1,1:nbandm)=atab(1:nbandm)
           scprod(2,1:nbandm)=zero
           ABI_DEALLOCATE(atab)
         end if

         do iband=1,nbandm,iblock
           index1=npw_sp*(iband-1)+icg
           index2=npw_sp*iband+icg
           if (iband/=iband0) then
             ar=scprod(1,iband)
           else
             ar=zero
           end if
           if (iband+1/=iband0) then
             ar2=scprod(1,iband+1)
           else
             ar2=zero
           end if
!          $OMP PARALLEL DO PRIVATE(cg_re,cg_re2,cg_im,cg_im2) &
!          $OMP&SHARED(ar,ar2,cg,direc,index1,index2,npw_sp)
           do ipw=1,npw_sp
             cg_re =cg(1,index1+ipw) ; cg_im =cg(2,index1+ipw)
             cg_re2=cg(1,index2+ipw) ; cg_im2=cg(2,index2+ipw)
             direc(1,ipw)=direc(1,ipw)-ar*cg_re-ar2*cg_re2
             direc(2,ipw)=direc(2,ipw)-ar*cg_im-ar2*cg_im2
           end do
!          $OMP END PARALLEL DO
         end do ! Loop on iband
       end if ! Test on istwf_k
     end if ! Test on nband

   else if (ortalg==3 .or. ortalg==-3)then

     iblock=3;nbandm=(nband/iblock)*iblock

     if(nband>=iblock)then

       if(istwf_k==1)then

         if (scprod_io==0) then
           ABI_ALLOCATE(atab,(2*nbandm))
           atab=zero
           if (useoverlap==1) then
             do iband=1,nbandm,iblock
               ar=zero  ; ai=zero
               ar2=zero ; ai2=zero
               ar3=zero ; ai3=zero
               index1=npw_sp*(iband-1)+iscg
               index2=npw_sp* iband   +iscg
               index3=npw_sp*(iband+1)+iscg
!              $OMP PARALLEL DO PRIVATE(scg_re,scg_re2,scg_re3,scg_im,scg_im2) &
!              $OMP&PRIVATE(scg_im3,direc_re,direc_im,ipw) &
!              $OMP&REDUCTION(+ :ai,ai2,ai3,ar,ar2,ar3) &
!              $OMP&SHARED(scg,direc,index1,index2,index3,npw_sp)
               do ipw=1,npw_sp
                 direc_re=direc(1,ipw)     ; direc_im=direc(2,ipw)
                 scg_re =scg(1,index1+ipw) ; scg_im =scg(2,index1+ipw)
                 scg_re2=scg(1,index2+ipw) ; scg_im2=scg(2,index2+ipw)
                 scg_re3=scg(1,index3+ipw) ; scg_im3=scg(2,index3+ipw)
                 ar=ar+scg_re*direc_re+scg_im*direc_im
                 ai=ai-scg_im*direc_re+scg_re*direc_im
                 ar2=ar2+scg_re2*direc_re+scg_im2*direc_im
                 ai2=ai2-scg_im2*direc_re+scg_re2*direc_im
                 ar3=ar3+scg_re3*direc_re+scg_im3*direc_im
                 ai3=ai3-scg_im3*direc_re+scg_re3*direc_im
               end do
!              $OMP END PARALLEL DO
               iband2=2*iband
               atab(iband2-1)=ar ;atab(iband2  )=ai
               atab(iband2+1)=ar2;atab(iband2+2)=ai2
               atab(iband2+3)=ar3;atab(iband2+4)=ai3
             end do ! Loop on iband
           else
             do iband=1,nbandm,iblock
               ar=zero  ; ai=zero
               ar2=zero ; ai2=zero
               ar3=zero ; ai3=zero
               index1=npw_sp*(iband-1)+icg
               index2=npw_sp*iband+icg
               index3=npw_sp*(iband+1)+icg
!              $OMP PARALLEL DO PRIVATE(cg_re,cg_re2,cg_re3,cg_im,cg_im2) &
!              $OMP&PRIVATE(cg_im3,direc_re,direc_im,ipw) &
!              $OMP&REDUCTION(+ :ai,ai2,ai3,ar,ar2,ar3) &
!              $OMP&SHARED(cg,direc,index1,index2,index3,npw_sp)
               do ipw=1,npw_sp
                 direc_re=direc(1,ipw)      ; direc_im=direc(2,ipw)
                 cg_re=cg(1,index1+ipw)  ; cg_im=cg(2,index1+ipw)
                 cg_re2=cg(1,index2+ipw) ; cg_im2=cg(2,index2+ipw)
                 cg_re3=cg(1,index3+ipw) ; cg_im3=cg(2,index3+ipw)
                 ar=ar+cg_re*direc_re+cg_im*direc_im
                 ai=ai-cg_im*direc_re+cg_re*direc_im
                 ar2=ar2+cg_re2*direc_re+cg_im2*direc_im
                 ai2=ai2-cg_im2*direc_re+cg_re2*direc_im
                 ar3=ar3+cg_re3*direc_re+cg_im3*direc_im
                 ai3=ai3-cg_im3*direc_re+cg_re3*direc_im
               end do
!              $OMP END PARALLEL DO
               iband2=2*iband
               atab(iband2-1)=ar ;atab(iband2  )=ai
               atab(iband2+1)=ar2;atab(iband2+2)=ai2
               atab(iband2+3)=ar3;atab(iband2+4)=ai3
             end do ! Loop on iband
           end if ! useoverlap
           if (mpi_enreg%paral_compil_fft==1) then
             call xcomm_init(mpi_enreg,spaceComm)
             call timab(48,1,tsec)
             call xsum_mpi(atab,spaceComm,ierr)
             call timab(48,2,tsec)
           end if
           scprod(:,1:nbandm)=reshape(atab(1:2*nbandm),(/2,nbandm/))
           ABI_DEALLOCATE(atab)
         end if

         do iband=1,nbandm,iblock
           index1=npw_sp*(iband-1)+icg
           index2=npw_sp* iband   +icg
           index3=npw_sp*(iband+1)+icg
           if (iband/=iband0) then
             ar=scprod(1,iband);ai=scprod(2,iband)
           else
             ar=zero;ai=zero
           end if
           if (iband+1/=iband0) then
             ar2=scprod(1,iband+1);ai2=scprod(2,iband+1)
           else
             ar2=zero;ai2=zero
           end if
           if (iband+2/=iband0) then
             ar3=scprod(1,iband+2);ai3=scprod(2,iband+2)
           else
             ar3=zero;ai3=zero
           end if
!          $OMP PARALLEL DO PRIVATE(cg_re,cg_re2,cg_re3,cg_im,cg_im2) &
!          $OMP&PRIVATE(cg_im3,ipw) &
!          $OMP&SHARED(ar,ai,ar2,ai2,ar3,ai3,cg,direc,index1,index2,index3,npw_sp)
           do ipw=1,npw_sp
             cg_re =cg(1,index1+ipw) ; cg_im =cg(2,index1+ipw)
             cg_re2=cg(1,index2+ipw) ; cg_im2=cg(2,index2+ipw)
             cg_re3=cg(1,index3+ipw) ; cg_im3=cg(2,index3+ipw)
             direc(1,ipw)=direc(1,ipw)-ar *cg_re +ai *cg_im-ar2*cg_re2+ai2*cg_im2 &
&             -ar3*cg_re3+ai3*cg_im3
             direc(2,ipw)=direc(2,ipw)-ar *cg_im -ai *cg_re-ar2*cg_im2-ai2*cg_re2 &
&             -ar3*cg_im3-ai3*cg_re3
           end do
!          $OMP END PARALLEL DO
         end do ! Loop on iband

       else if(istwf_k>=2)then

         if (scprod_io==0) then
           ABI_ALLOCATE(atab,(nbandm))
           atab=zero
           if (useoverlap==1) then
             do iband=1,nbandm,iblock
               index1=npw_sp*(iband-1)+iscg
               index2=npw_sp* iband   +iscg
               index3=npw_sp*(iband+1)+iscg
               if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
                 ar =half*scg(1,index1 +1)*direc(1,1)
                 ar2=half*scg(1,index2+1)*direc(1,1)
                 ar3=half*scg(1,index3+1)*direc(1,1)
                 ipw1=2
               else
                 ar=zero ; ar2=zero ; ar3=zero ; ipw1=1
               end if
!              $OMP PARALLEL DO PRIVATE(scg_re,scg_re2,scg_re3,scg_im,scg_im2) &
!              $OMP&PRIVATE(scg_im3,direc_re,direc_im,ipw) &
!              $OMP&REDUCTION(+ :ar,ar2,ar3) &
!              $OMP&SHARED(scg,direc,index1,index2,index3,npw,nspinor)
               do isp=1,nspinor
                 do ipw=ipw1+(isp-1)*npw,npw*isp
                   direc_re=direc(1,ipw)     ; direc_im=direc(2,ipw)
                   scg_re=scg(1,index1+ipw)  ; scg_im=scg(2,index1+ipw)
                   scg_re2=scg(1,index2+ipw) ; scg_im2=scg(2,index2+ipw)
                   scg_re3=scg(1,index3+ipw) ; scg_im3=scg(2,index3+ipw)
                   ar =ar +scg_re *direc_re+scg_im *direc_im
                   ar2=ar2+scg_re2*direc_re+scg_im2*direc_im
                   ar3=ar3+scg_re3*direc_re+scg_im3*direc_im
                 end do
               end do
!              $OMP END PARALLEL DO
               atab(iband  )=two*ar
               atab(iband+1)=two*ar2
               atab(iband+2)=two*ar3
             end do
           else
             do iband=1,nbandm,iblock
               index1=npw_sp*(iband-1)+icg
               index2=npw_sp* iband   +icg
               index3=npw_sp*(iband+1)+icg
               if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
                 ar =half*cg(1,index1+1)*direc(1,1)
                 ar2=half*cg(1,index2+1)*direc(1,1)
                 ar3=half*cg(1,index3+1)*direc(1,1)
                 ipw1=2
               else
                 ar=zero ; ar2=zero ; ar3=zero ; ipw1=1
               end if
!              $OMP PARALLEL DO PRIVATE(cg_re,cg_re2,cg_re3,cg_im,cg_im2) &
!              $OMP&PRIVATE(cg_im3,direc_re,direc_im,ipw) &
!              $OMP&REDUCTION(+ :ar,ar2,ar3) &
!              $OMP&SHARED(cg,direc,index1,index2,index3,npw,nspinor)
               do isp=1,nspinor
                 do ipw=ipw1+(isp-1)*npw,npw*isp
                   direc_re=direc(1,ipw)   ; direc_im=direc(2,ipw)
                   cg_re=cg(1,index1+ipw)  ; cg_im=cg(2,index1+ipw)
                   cg_re2=cg(1,index2+ipw) ; cg_im2=cg(2,index2+ipw)
                   cg_re3=cg(1,index3+ipw) ; cg_im3=cg(2,index3+ipw)
                   ar=ar+cg_re*direc_re+cg_im*direc_im
                   ar2=ar2+cg_re2*direc_re+cg_im2*direc_im
                   ar3=ar3+cg_re3*direc_re+cg_im3*direc_im
                 end do
               end do
!              $OMP END PARALLEL DO
               atab(iband  )=two*ar
               atab(iband+1)=two*ar2
               atab(iband+2)=two*ar3
             end do ! Loop on iband
           end if ! useoverlap
           if (mpi_enreg%paral_compil_fft==1) then
             call xcomm_init(mpi_enreg,spaceComm)
             call timab(48,1,tsec)
             call xsum_mpi(atab,spaceComm,ierr)
             call timab(48,2,tsec)
           end if
           scprod(1,1:nbandm)=atab(1:nbandm)
           scprod(2,1:nbandm)=zero
           ABI_DEALLOCATE(atab)
         end if

         do iband=1,nbandm,iblock
           index1=npw_sp*(iband-1)+icg
           index2=npw_sp* iband   +icg
           index3=npw_sp*(iband+1)+icg
           if (iband/=iband0) then
             ar=scprod(1,iband)
           else
             ar=zero
           end if
           if (iband+1/=iband0) then
             ar2=scprod(1,iband+1)
           else
             ar2=zero
           end if
           if (iband+2/=iband0) then
             ar3=scprod(1,iband+2)
           else
             ar3=zero
           end if
!          $OMP PARALLEL DO PRIVATE(cg_re,cg_re2,cg_re3,cg_im,cg_im2) &
!          $OMP&PRIVATE(cg_im3,direc_re,direc_im,ipw) &
!          $OMP&SHARED(ar,ar2,ar3,cg,direc,index1,index2,index3,npw_sp)
           do ipw=1,npw_sp
             cg_re =cg(1,index1+ipw) ; cg_im =cg(2,index1+ipw)
             cg_re2=cg(1,index2+ipw) ; cg_im2=cg(2,index2+ipw)
             cg_re3=cg(1,index3+ipw) ; cg_im3=cg(2,index3+ipw)
             direc(1,ipw)=direc(1,ipw)-ar*cg_re-ar2*cg_re2-ar3*cg_re3
             direc(2,ipw)=direc(2,ipw)-ar*cg_im-ar2*cg_im2-ar3*cg_im3
           end do
!          $OMP END PARALLEL DO
         end do ! Loop on iband
       end if ! Test on istwf_k
     end if ! Test on nband

   else if (ortalg==4 .or. ortalg==-4)then

     iblock=4;nbandm=(nband/iblock)*iblock

     if(nband>=iblock)then

       if(istwf_k==1)then

         if (scprod_io==0) then
           ABI_ALLOCATE(atab,(2*nbandm))
           atab=zero
           if (useoverlap==1) then
             do iband=1,nbandm,iblock
               ar=zero  ; ai=zero
               ar2=zero ; ai2=zero
               ar3=zero ; ai3=zero
               ar4=zero ; ai4=zero
               index1=npw_sp*(iband-1)+iscg
               index2=npw_sp* iband   +iscg
               index3=npw_sp*(iband+1)+iscg
               index4=npw_sp*(iband+2)+iscg
!              $OMP PARALLEL DO PRIVATE(scg_re,scg_re2,scg_re3,scg_re4,scg_im) &
!              $OMP&PRIVATE(scg_im2,scg_im3,scg_im4,direc_re,direc_im,ipw) &
!              $OMP&REDUCTION(+ :ai,ai2,ai3,ai4,ar,ar2,ar3,ar4) &
!              $OMP&SHARED(scg,direc) &
!              $OMP&SHARED(index1,index2,index3,index4,npw_sp)
               do ipw=1,npw_sp
                 direc_re=direc(1,ipw)     ; direc_im=direc(2,ipw)
                 scg_re=scg(1,index1+ipw)  ; scg_im=scg(2,index1+ipw)
                 scg_re2=scg(1,index2+ipw) ; scg_im2=scg(2,index2+ipw)
                 scg_re3=scg(1,index3+ipw) ; scg_im3=scg(2,index3+ipw)
                 scg_re4=scg(1,index4+ipw) ; scg_im4=scg(2,index4+ipw)
                 ar=ar+scg_re*direc_re+scg_im*direc_im
                 ai=ai-scg_im*direc_re+scg_re*direc_im
                 ar2=ar2+scg_re2*direc_re+scg_im2*direc_im
                 ai2=ai2-scg_im2*direc_re+scg_re2*direc_im
                 ar3=ar3+scg_re3*direc_re+scg_im3*direc_im
                 ai3=ai3-scg_im3*direc_re+scg_re3*direc_im
                 ar4=ar4+scg_re4*direc_re+scg_im4*direc_im
                 ai4=ai4-scg_im4*direc_re+scg_re4*direc_im
               end do
!              $OMP END PARALLEL DO
               iband2=2*iband
               atab(iband2-1)=ar ;atab(iband2  )=ai
               atab(iband2+1)=ar2;atab(iband2+2)=ai2
               atab(iband2+3)=ar3;atab(iband2+4)=ai3
               atab(iband2+5)=ar4;atab(iband2+6)=ai4
             end do ! Loop on iband
           else
             do iband=1,nbandm,iblock
               ar=zero  ; ai=zero
               ar2=zero ; ai2=zero
               ar3=zero ; ai3=zero
               ar4=zero ; ai4=zero
               index1=npw_sp*(iband-1)+icg
               index2=npw_sp* iband   +icg
               index3=npw_sp*(iband+1)+icg
               index4=npw_sp*(iband+2)+icg
!              $OMP PARALLEL DO PRIVATE(cg_re,cg_re2,cg_re3,cg_re4,cg_im) &
!              $OMP&PRIVATE(cg_im2,cg_im3,cg_im4,direc_re,direc_im,ipw) &
!              $OMP&REDUCTION(+ :ai,ai2,ai3,ai4,ar,ar2,ar3,ar4) &
!              $OMP&SHARED(cg,direc) &
!              $OMP&SHARED(index1,index2,index3,index4,npw_sp)
               do ipw=1,npw_sp
                 direc_re=direc(1,ipw)   ; direc_im=direc(2,ipw)
                 cg_re =cg(1,index1+ipw) ; cg_im =cg(2,index1+ipw)
                 cg_re2=cg(1,index2+ipw) ; cg_im2=cg(2,index2+ipw)
                 cg_re3=cg(1,index3+ipw) ; cg_im3=cg(2,index3+ipw)
                 cg_re4=cg(1,index4+ipw) ; cg_im4=cg(2,index4+ipw)
                 ar=ar+cg_re*direc_re+cg_im*direc_im
                 ai=ai-cg_im*direc_re+cg_re*direc_im
                 ar2=ar2+cg_re2*direc_re+cg_im2*direc_im
                 ai2=ai2-cg_im2*direc_re+cg_re2*direc_im
                 ar3=ar3+cg_re3*direc_re+cg_im3*direc_im
                 ai3=ai3-cg_im3*direc_re+cg_re3*direc_im
                 ar4=ar4+cg_re4*direc_re+cg_im4*direc_im
                 ai4=ai4-cg_im4*direc_re+cg_re4*direc_im
               end do
!              $OMP END PARALLEL DO
               iband2=2*iband
               atab(iband2-1)=ar ;atab(iband2  )=ai
               atab(iband2+1)=ar2;atab(iband2+2)=ai2
               atab(iband2+3)=ar3;atab(iband2+4)=ai3
               atab(iband2+5)=ar4;atab(iband2+6)=ai4
             end do ! Loop on iband
           end if ! useoverlap
           if (mpi_enreg%paral_compil_fft==1) then
             call xcomm_init(mpi_enreg,spaceComm)
             call timab(48,1,tsec)
             call xsum_mpi(atab,spaceComm,ierr)
             call timab(48,2,tsec)
           end if
           scprod(:,1:nbandm)=reshape(atab(1:2*nbandm),(/2,nbandm/))
           ABI_DEALLOCATE(atab)
         end if

         do iband=1,nbandm,iblock
           index1=npw_sp*(iband-1)+icg
           index2=npw_sp* iband   +icg
           index3=npw_sp*(iband+1)+icg
           index4=npw_sp*(iband+2)+icg
           if (iband/=iband0) then
             ar=scprod(1,iband);ai=scprod(2,iband)
           else
             ar=zero;ai=zero
           end if
           if (iband+1/=iband0) then
             ar2=scprod(1,iband+1);ai2=scprod(2,iband+1)
           else
             ar2=zero;ai2=zero
           end if
           if (iband+2/=iband0) then
             ar3=scprod(1,iband+2);ai3=scprod(2,iband+2)
           else
             ar3=zero;ai3=zero
           end if
           if (iband+3/=iband0) then
             ar4=scprod(1,iband+3);ai4=scprod(2,iband+3)
           else
             ar4=zero;ai4=zero
           end if
!          $OMP PARALLEL DO PRIVATE(cg_re,cg_re2,cg_re3,cg_re4,cg_im) &
!          $OMP&PRIVATE(cg_im2,cg_im3,cg_im4,ipw) &
!          $OMP&SHARED(ar,ai,ar2,ai2,ar3,ai3,ar4,ai4,cg,direc) &
!          $OMP&SHARED(index1,index2,index3,index4,npw_sp)
           do ipw=1,npw_sp
             cg_re =cg(1,index1+ipw) ; cg_im =cg(2,index1+ipw)
             cg_re2=cg(1,index2+ipw) ; cg_im2=cg(2,index2+ipw)
             cg_re3=cg(1,index3+ipw) ; cg_im3=cg(2,index3+ipw)
             cg_re4=cg(1,index4+ipw) ; cg_im4=cg(2,index4+ipw)
             direc(1,ipw)=direc(1,ipw)-ar *cg_re +ai *cg_im -ar2*cg_re2+ai2*cg_im2 &
&             -ar3*cg_re3+ai3*cg_im3-ar4*cg_re4+ai4*cg_im4
             direc(2,ipw)=direc(2,ipw)-ar *cg_im -ai *cg_re -ar2*cg_im2-ai2*cg_re2 &
&             -ar3*cg_im3-ai3*cg_re3-ar4*cg_im4-ai4*cg_re4
           end do
!          $OMP END PARALLEL DO
         end do ! Loop on iband

       else if(istwf_k>=2)then

         if (scprod_io==0) then
           ABI_ALLOCATE(atab,(nbandm))
           atab=zero
           if (useoverlap==1) then
             do iband=1,nbandm,iblock
               index1=npw_sp*(iband-1)+iscg
               index2=npw_sp* iband   +iscg
               index3=npw_sp*(iband+1)+iscg
               index4=npw_sp*(iband+2)+iscg
               if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
                 ar =half*scg(1,index1 +1)*direc(1,1)
                 ar2=half*scg(1,index2+1)*direc(1,1)
                 ar3=half*scg(1,index3+1)*direc(1,1)
                 ar4=half*scg(1,index4+1)*direc(1,1)
                 ipw1=2
               else
                 ar=zero ; ar2=zero ; ar3=zero ; ar4=zero ; ipw1=1
               end if
!              $OMP PARALLEL DO PRIVATE(scg_re,scg_re2,scg_re3,scg_re4,scg_im) &
!              $OMP&PRIVATE(scg_im2,scg_im3,scg_im4,direc_re,direc_im,ipw) &
!              $OMP&REDUCTION(+ :ar,ar2,ar3,ar4) &
!              $OMP&SHARED(scg,direc) &
!              $OMP&SHARED(index1,index2,index3,index4,npw,nspinor)
               do isp=1,nspinor
                 do ipw=ipw1+(isp-1)*npw,npw*isp
                   direc_re=direc(1,ipw)     ; direc_im=direc(2,ipw)
                   scg_re=scg(1,index1+ipw)  ; scg_im=scg(2,index1+ipw)
                   scg_re2=scg(1,index2+ipw) ; scg_im2=scg(2,index2+ipw)
                   scg_re3=scg(1,index3+ipw) ; scg_im3=scg(2,index3+ipw)
                   scg_re4=scg(1,index4+ipw) ; scg_im4=scg(2,index4+ipw)
                   ar=ar+scg_re*direc_re+scg_im*direc_im
                   ar2=ar2+scg_re2*direc_re+scg_im2*direc_im
                   ar3=ar3+scg_re3*direc_re+scg_im3*direc_im
                   ar4=ar4+scg_re4*direc_re+scg_im4*direc_im
                 end do
               end do
               atab(iband  )=two*ar
               atab(iband+1)=two*ar2
               atab(iband+2)=two*ar3
               atab(iband+3)=two*ar4
             end do ! Loop on iband
           else
             do iband=1,nbandm,iblock
               index1=npw_sp*(iband-1)+icg
               index2=npw_sp* iband   +icg
               index3=npw_sp*(iband+1)+icg
               index4=npw_sp*(iband+2)+icg
               if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
                 ar =half*cg(1,index1+1)*direc(1,1)
                 ar2=half*cg(1,index2+1)*direc(1,1)
                 ar3=half*cg(1,index3+1)*direc(1,1)
                 ar4=half*cg(1,index4+1)*direc(1,1)
                 ipw1=2
               else
                 ar=zero ; ar2=zero ; ar3=zero ; ar4=zero ; ipw1=1
               end if
!              $OMP PARALLEL DO PRIVATE(cg_re,cg_re2,cg_re3,cg_re4,cg_im) &
!              $OMP&PRIVATE(cg_im2,cg_im3,cg_im4,direc_re,direc_im,ipw) &
!              $OMP&REDUCTION(+ :ar,ar2,ar3,ar4) &
!              $OMP&SHARED(cg,direc) &
!              $OMP&SHARED(index1,index2,index3,index4,npw,nspinor)
               do isp=1,nspinor
                 do ipw=ipw1+(isp-1)*npw,npw*isp
                   direc_re=direc(1,ipw)   ; direc_im=direc(2,ipw)
                   cg_re=cg(1,index1+ipw)  ; cg_im=cg(2,index1+ipw)
                   cg_re2=cg(1,index2+ipw) ; cg_im2=cg(2,index2+ipw)
                   cg_re3=cg(1,index3+ipw) ; cg_im3=cg(2,index3+ipw)
                   cg_re4=cg(1,index4+ipw) ; cg_im4=cg(2,index4+ipw)
                   ar=ar+cg_re*direc_re+cg_im*direc_im
                   ar2=ar2+cg_re2*direc_re+cg_im2*direc_im
                   ar3=ar3+cg_re3*direc_re+cg_im3*direc_im
                   ar4=ar4+cg_re4*direc_re+cg_im4*direc_im
                 end do
               end do
!              $OMP END PARALLEL DO
               atab(iband  )=two*ar
               atab(iband+1)=two*ar2
               atab(iband+2)=two*ar3
               atab(iband+3)=two*ar4
             end do ! Loop on iband
           end if ! useoverlap
           if (mpi_enreg%paral_compil_fft==1) then
             call xcomm_init(mpi_enreg,spaceComm)
             call timab(48,1,tsec)
             call xsum_mpi(atab,spaceComm,ierr)
             call timab(48,2,tsec)
           end if
           scprod(1,1:nbandm)=atab(1:nbandm)
           scprod(2,1:nbandm)=zero
           ABI_DEALLOCATE(atab)
         end if

         do iband=1,nbandm,iblock
           index1=npw_sp*(iband-1)+icg
           index2=npw_sp* iband   +icg
           index3=npw_sp*(iband+1)+icg
           index4=npw_sp*(iband+2)+icg
           if (iband/=iband0) then
             ar=scprod(1,iband)
           else
             ar=zero
           end if
           if (iband+1/=iband0) then
             ar2=scprod(1,iband+1)
           else
             ar2=zero
           end if
           if (iband+2/=iband0) then
             ar3=scprod(1,iband+2)
           else
             ar3=zero
           end if
           if (iband+3/=iband0) then
             ar4=scprod(1,iband+3)
           else
             ar4=zero
           end if
!          $OMP PARALLEL DO PRIVATE(cg_re,cg_re2,cg_re3,cg_re4,cg_im) &
!          $OMP&PRIVATE(cg_im2,cg_im3,cg_im4,ipw) &
!          $OMP&SHARED(ar,ar2,ar3,ar4,cg,direc) &
!          $OMP&SHARED(index1,index2,index3,index4,npw_sp)
           do ipw=1,npw_sp
             cg_re=cg(1,index1+ipw)  ; cg_im=cg(2,index1+ipw)
             cg_re2=cg(1,index2+ipw) ; cg_im2=cg(2,index2+ipw)
             cg_re3=cg(1,index3+ipw) ; cg_im3=cg(2,index3+ipw)
             cg_re4=cg(1,index4+ipw) ; cg_im4=cg(2,index4+ipw)
             direc(1,ipw)=direc(1,ipw)-ar *cg_re -ar2*cg_re2 &
&             -ar3*cg_re3-ar4*cg_re4
             direc(2,ipw)=direc(2,ipw)-ar *cg_im -ar2*cg_im2 &
&             -ar3*cg_im3-ar4*cg_im4
           end do
!          $OMP END PARALLEL DO
         end do ! Loop on iband
       end if ! Test on istwf_k
     end if ! Test on nband
   end if ! Test on ortalg=2, 3 or 4

!  Need to treat the bands not yet treated
   if( nbandm /= nband )then

     if(istwf_k==1)then

       do iband=nbandm+1,nband

         if (scprod_io==0) then
           if (useoverlap==1) then
             ar=zero ; ai=zero
             index1=npw_sp*(iband-1)+iscg
!            $OMP PARALLEL DO PRIVATE(scg_re,scg_im,direc_re,direc_im,ipw) &
!            $OMP&REDUCTION(+ :ai,ar) &
!            $OMP&SHARED(scg,direc,index1,npw_sp)
             do ipw=1,npw_sp
               direc_re=direc(1,ipw)    ; direc_im=direc(2,ipw)
               scg_re=scg(1,index1+ipw) ; scg_im=scg(2,index1+ipw)
               ar=ar+scg_re*direc_re+scg_im*direc_im
               ai=ai-scg_im*direc_re+scg_re*direc_im
             end do
!            $OMP END PARALLEL DO
           else
             ar=zero ; ai=zero
             index1=npw_sp*(iband-1)+icg
!            $OMP PARALLEL DO PRIVATE(cg_re,cg_im,direc_re,direc_im,ipw) &
!            $OMP&REDUCTION(+ :ai,ar) &
!            $OMP&SHARED(cg,direc,index1,npw_sp)
             do ipw=1,npw_sp
               direc_re=direc(1,ipw)  ; direc_im=direc(2,ipw)
               cg_re=cg(1,index1+ipw) ; cg_im=cg(2,index1+ipw)
               ar=ar+cg_re*direc_re+cg_im*direc_im
               ai=ai-cg_im*direc_re+cg_re*direc_im
             end do
!            $OMP END PARALLEL DO
           end if ! useoverlap
           if (mpi_enreg%paral_compil_fft==1) then
             call xcomm_init(mpi_enreg,spaceComm)
             buffer2(1)=ar;buffer2(2)=ai
             call timab(48,1,tsec)
             call xsum_mpi(buffer2,spaceComm,ierr)
             call timab(48,2,tsec)
             ar=buffer2(1);ai=buffer2(2)
           end if
           scprod(1,iband)=ar  ; scprod(2,iband)=ai
         end if

         if (iband==iband0) cycle
         index1=npw_sp*(iband-1)+icg
!        $OMP PARALLEL DO PRIVATE(cg_re,cg_im,ipw) &
!        $OMP&SHARED(ai,ar,cg,direc,index1,npw_sp)
         do ipw=1,npw_sp
           cg_re=cg(1,index1+ipw) ; cg_im=cg(2,index1+ipw)
           direc(1,ipw)=direc(1,ipw)-ar*cg_re+ai*cg_im
           direc(2,ipw)=direc(2,ipw)-ar*cg_im-ai*cg_re
         end do
!        $OMP END PARALLEL DO
       end do ! Loop on iband

     else if(istwf_k>=2)then

       do iband=nbandm+1,nband

         if (scprod_io==0) then
           if (useoverlap==1) then
             index1=npw_sp*(iband-1)+iscg
             if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
               ar=half*scg(1,1+index1)*direc(1,1) ; ipw1=2
             else
               ar=zero ; ipw1=1
             end if
!            $OMP PARALLEL DO PRIVATE(scg_re,scg_im,direc_re,direc_im,ipw) &
!            $OMP&REDUCTION(+ :ar) &
!            $OMP&SHARED(scg,direc,index1,npw,nspinor)
             do isp=1,nspinor
               do ipw=ipw1+(isp-1)*npw,npw*isp
                 direc_re=direc(1,ipw)    ; direc_im=direc(2,ipw)
                 scg_re=scg(1,index1+ipw) ; scg_im=scg(2,index1+ipw)
                 ar=ar+scg_re*direc_re+scg_im*direc_im
               end do
             end do
!            $OMP END PARALLEL DO
             ar=two*ar
           else
             index1=npw_sp*(iband-1)+icg
             if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
               ar=half*cg(1,1+index1)*direc(1,1) ; ipw1=2
             else
               ar=zero ; ipw1=1
             end if
!            $OMP PARALLEL DO PRIVATE(cg_re,cg_im,direc_re,direc_im,ipw) &
!            $OMP&REDUCTION(+ :ar) &
!            $OMP&SHARED(cg,direc,index1,npw,nspinor)
             do isp=1,nspinor
               do ipw=ipw1+(isp-1)*npw,npw*isp
                 direc_re=direc(1,ipw)  ; direc_im=direc(2,ipw)
                 cg_re=cg(1,index1+ipw) ; cg_im=cg(2,index1+ipw)
                 ar=ar+cg_re*direc_re+cg_im*direc_im
               end do
             end do
!            $OMP END PARALLEL DO
             ar=two*ar
           end if ! useoverlap
           if (mpi_enreg%paral_compil_fft==1) then
             call xcomm_init(mpi_enreg,spaceComm)
             call timab(48,1,tsec)
             call xsum_mpi(ar,spaceComm,ierr)
             call timab(48,2,tsec)
           end if
           scprod(1,iband)=ar ; scprod(2,iband)=zero
         end if

         if (iband==iband0) cycle
         index1=npw_sp*(iband-1)+icg
!        $OMP PARALLEL DO PRIVATE(cg_re,cg_im,ipw) &
!        $OMP&SHARED(ar,cg,direc,index1,npw_sp)
         do ipw=1,npw_sp
           cg_re=cg(1,index1+ipw) ; cg_im=cg(2,index1+ipw)
           direc(1,ipw)=direc(1,ipw)-ar*cg_re
           direc(2,ipw)=direc(2,ipw)-ar*cg_im
         end do
!        $OMP END PARALLEL DO
       end do ! Loop on iband
     end if ! Test on istwf_k
   end if ! Test on nband

 end if ! on ortalg=0,1, 2 or 4

 mpi_enreg%paral_level=old_paral_level
 call timab(210+tim_projbd,2,tsec)

!DEBUG
!write(std_out,*)' projbd: debug, enter.'
!ENDDEBUG

end subroutine projbd

!ALTERNATIVE ROUTINE TO BE EVENTUALLY USED FOR PAW
!KEEP IT FOR THE TIME BEING
!==================================================
!!{\src2tex{textfont=tt}}
!!!****f* ABINIT/projbd1
!!!
!!! NAME
!!! projbd1
!!!
!!! FUNCTION
!!! Project out a band cg_{i} onto the other bands contained in cg_{j}.
!!!  direc=$sum_{j/=i} { <cg_{j}|S|cg_{i}>.|cg_{j}> }$
!!!    (s=overlap matrix) - only use for PAW
!!!
!!! COPYRIGHT
!!! Copyright (C) 2008-2012 ABINIT group (FJ,MT)
!!! This file is distributed under the terms of the
!!! GNU General Public License, see ~abinit/COPYING
!!! or http://www.gnu.org/copyleft/gpl.txt .
!!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!!
!!! INPUTS
!!!  cg(2,mcg)=wavefunction coefficients for ALL bands
!!!  iband0=which particular band we are interested in
!!!         ("i" in the above formula)
!!!         Can be set to -1 to sum over all bands...
!!!  icg=shift to be given to the location of the data in cg
!!!  istwf_k=option parameter that describes the storage of wfs
!!!  mcg=maximum size of second dimension of cg
!!!  mpi_enreg=informations about MPI parallelization
!!!  nband=number of bands
!!!  npw=number of planewaves
!!!  nspinor=number of spinorial components
!!!  printopt= if 1, print intermediate dot products
!!!  scg(2,npw)=auxilliary vector, S|cg_{i}>
!!!  tim_projbd=timing code of the calling subroutine(can be set to 0 if not attributed)
!!!
!!! OUTPUT
!!!  direc(2,npw)= $sum_{j} { <cg_{j}|S|cg_{i}>.|cg_{j}> }$
!!!
!!! PARENTS
!!!      cgwf3
!!!
!!! CHILDREN
!!!      timab,wrtout,xcomm_init,xsum_mpi
!!!
!!! SOURCE
!
!#if defined HAVE_CONFIG_H
!#include "config.h"
!#endif
!
!subroutine projbd1(cg,direc,iband0,icg,istwf_k,mcg,mpi_enreg,nband,&
!&                  npw,nspinor,printopt,scg,tim_projbd)
!
!use defs_basis
!use defs_datatypes
!use defs_abitypes
!
!!This section has been created automatically by the script Abilint (TD).
!!Do not modify the following lines by hand.
!use interfaces_51_manage_mpi
!use interfaces_12_hide_mpi
!!End of the abilint section
!
!implicit none
!
!!Arguments ------------------------------------
!!This type is defined in defs_mpi
!!scalars
!integer,intent(in) :: iband0,icg,istwf_k,mcg,nband,npw,nspinor,printopt
!integer,intent(in) :: tim_projbd
!type(MPI_type),intent(inout) :: mpi_enreg
!!arrays
!real(dp),intent(in) :: cg(2,mcg),scg(2,npw)
!real(dp),intent(out) :: direc(2,npw*nspinor)
!
!!Local variables-------------------------------
!!scalars
!integer :: iband,iblock,ierr,index1,ipw,ipw1,isp
!integer :: nbandm,npw_sp,old_paral_level,spaceComm
!real(dp) :: ai,ar,cg_im,cg_re,direc_im,direc_re,scg_im,scg_re
!character(len=500) :: message
!!arrays
!real(dp) :: tsec(2)
!real(dp),allocatable :: atab1(:),atab2(:,:)
!
!! *************************************************************************
!!
!
!!DEBUG
!!write(std_out,*)' projbd1 : enter '
!!stop
!!ENDDEBUG
!old_paral_level=mpi_enreg%paral_level
!mpi_enreg%paral_level=3
!call timab(210+tim_projbd,1,tsec)
!
!npw_sp=npw*nspinor
!nbandm=nband
!direc(:,:)=zero
!
!if(istwf_k==1)then
!
!allocate(atab2(2,nbandm));atab2=zero
!
!do iband=1,nbandm
!ar=zero ; ai=zero
!index1=npw_sp*(iband-1)+icg
!!  $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ai,ar) &
!!  $OMP&SHARED(cg,direc,index1,npw_sp)
!do ipw=1,npw_sp
!ar=ar+cg(1,index1+ipw)*scg(1,ipw)+cg(2,index1+ipw)*scg(2,ipw)
!ai=ai-cg(2,index1+ipw)*scg(1,ipw)+cg(1,index1+ipw)*scg(2,ipw)
!end do
!!  $OMP END PARALLEL DO
!atab2(1,iband)=ar;atab2(2,iband)=ai
!end do! Loop on iband
!
!if (mpi_enreg%paral_compil_fft==1) then
!call xcomm_init(mpi_enreg,spaceComm)
!call timab(48,1,tsec)
!call xsum_mpi(atab2,spaceComm,ierr)
!call timab(48,2,tsec)
!end if
!
!do iband=1,nbandm
!if (iband==iband0) cycle
!index1=npw_sp*(iband-1)+icg
!ar=atab2(1,iband);ai=atab2(2,iband)
!!  $OMP PARALLEL DO PRIVATE(ipw,cg_re,cg_im) &
!!  $OMP&SHARED(ar,ai,cg,direc,index1,npw_sp)
!do ipw=1,npw_sp
!cg_re=cg(1,index1+ipw) ; cg_im=cg(2,index1+ipw)
!direc(1,ipw)=direc(1,ipw)+ar*cg_re-ai*cg_im
!direc(2,ipw)=direc(2,ipw)+ar*cg_im+ai*cg_re
!end do
!!  $OMP END PARALLEL DO
!if(printopt==1)then
!write(message,'(a,i3,2f14.6)') &
!&    'projbd1 : called from cgwf3 ; iband,ar,ai=',iband,ar,ai
!call wrtout(std_out,message,'PERS')
!end if
!end do! Loop on iband
!deallocate(atab2)
!
!else if(istwf_k>=2)then
!
!allocate(atab1(nbandm));atab1=zero
!
!do iband=1,nbandm
!index1=npw_sp*(iband-1)+icg
!if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
!ar=half*cg(1,index1+1)*direc(1,1) ; ipw1=2
!else
!ar=zero ; ipw1=1
!end if
!!  $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar) &
!!  $OMP&SHARED(cg,direc,index1,ipw1,npw,nspinor)
!do isp=1,nspinor
!do ipw=ipw1+(isp-1)*npw,npw*isp
!ar=ar+cg(1,index1+ipw)*scg(1,ipw)+cg(2,index1+ipw)*scg(2,ipw)
!end do
!end do
!!  $OMP END PARALLEL DO
!atab1(iband)=two*ar
!end do! Loop on iband
!
!if (mpi_enreg%paral_compil_fft==1) then
!call xcomm_init(mpi_enreg,spaceComm)
!call timab(48,1,tsec)
!call xsum_mpi(atab1,spaceComm,ierr)
!call timab(48,2,tsec)
!end if
!
!do iband=1,nbandm
!if (iband==iband0) cycle
!index1=npw_sp*(iband-1)+icg
!ar=atab1(iband)
!!  $OMP PARALLEL DO PRIVATE(ipw) &
!!  $OMP&SHARED(ar,cg,direc,index1,npw_sp)
!do ipw=1,npw_sp
!direc(1,ipw)=direc(1,ipw)+ar*cg(1,index1+ipw)
!direc(2,ipw)=direc(2,ipw)+ar*cg(2,index1+ipw)
!end do
!!  $OMP END PARALLEL DO
!if(printopt==1)then
!write(message,'(a,i3,f14.6)') &
!&    'projbd1 : called from cgwf3 ; iband,ar=',iband,ar
!call wrtout(std_out,message,'PERS')
!end if
!end do! Loop on iband
!deallocate(atab1)
!
!end if! Test on istwf_k
!
!mpi_enreg%paral_level=old_paral_level
!call timab(210+tim_projbd,2,tsec)
!
!!DEBUG
!!write(std_out,*)' projbd1: debug, enter.'
!!ENDDEBUG
!
!end subroutine projbd1
!!***
