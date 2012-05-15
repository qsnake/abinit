!{\src2tex{textfont=tt}}
!!****f* ABINIT/approxdelta
!! NAME
!!  approxdelta
!!
!! FUNCTION
!!  Approximate the Dirac function using :
!!  method 1) a triangular funtion centered at the value egwdiff_re (Eq 17 of PRB 74, 035101 (2006) 
!!  method 2) a gaussian of witdth ep%spsmear expandended in Taylor serie 
!!  (at the moment only the 0-th moments) 
!!
!!  Subroutine needed to implement the calculation 
!!  of the polarizability using the spectral representation as proposed in :
!!  PRB 74, 035101 (2006) and PRB 61, 7172 (1999)
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nomegasf=number of frequencies in the grid for Im \chi_0
!!  omegasf(0:nomega+1)= frequencies (real)
!!  egwdiff_re = transition energy where the delta function is centered 
!!
!!  method= 1 : a triangular shaped function used to approximated the delta
!!          2 : gaussian approximation with standard deviation (smear)
!! smear= used only in case of method==2, defines the width of the gaussian
!!
!! OUTPUT
!!  wl = weight associated to omegal (last omega wich is smaller than egwdiff_re
!!  wr = weight associate to omegar  (first omega larger than egwdff_re
!!  iomegal= index in the array omegasf of the last frequency < egwdiff
!!  iomegar= index in the array omegasf of the first frequency > egwdiff
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!      flush_unit,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine approxdelta(nomegasf,omegasf,egwdiff_re,smear,iomegal,iomegar,wl,wr,spmeth)

 use m_profiling
    
 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'approxdelta'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomegasf,spmeth
 integer,intent(out) :: iomegal,iomegar
 real(dp),intent(in) :: egwdiff_re,smear  
 real(dp),intent(out) :: wl,wr
!arrays
 real(dp),intent(in) :: omegasf(nomegasf)

!Local variables-------------------------------
 integer :: io,iomega
 real(dp) :: omegal,omegar,deltal,deltar
 character(len=500) :: msg                 
! *************************************************************************
 
 iomega=-999
 do io=nomegasf,1,-1
   if (omegasf(io)<egwdiff_re) then
    iomega=io; EXIT
   end if 
 end do

 iomegal=iomega    ; omegal=omegasf(iomegal)
 iomegar=iomegal+1 ; omegar=omegasf(iomegar)

 SELECT CASE (spmeth) 

 CASE (1) ! Weights for triangular shaped function
   wr=  (egwdiff_re-omegal)/(omegar-omegal)
   wl= -(egwdiff_re-omegar)/(omegar-omegal)

 CASE (2) ! Weights for gaussian method (0-th moment)
   deltal=(egwdiff_re-omegal)/smear
   deltar=(omegar-egwdiff_re)/smear
   if (deltar>=deltal) then
     wl=EXP(-deltal*deltal)
     ! this value is used to avoid double counting and speed-up
     wr=huge(one)*1.d-10 
   else 
     wl=huge(one)*1.d-10
     wr=exp(-deltal*deltal) 
   end if 

 CASE DEFAULT
   write(msg,'(a,i4)')'Wrong value for spmeth = ',spmeth 
   MSG_BUG(msg)

 END SELECT 

end subroutine approxdelta
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/calc_kkweight
!! NAME
!!  calc_kkweight
!!
!! FUNCTION
!!  Calculate frequency dependent weights needed to perform the Hilbert transform 
!!
!!  Subroutine needed to implement the calculation 
!!  of the polarizability using the spectral representation as proposed in :
!!  PRB 74, 035101 (2006) and PRB 61, 7172 (1999)
!!  
!! INPUTS
!! nsp=number of frequencies where the imaginary part of the polarizability is evaluated
!! ne=number of frequencies for the polarizability (same as in epsilon^-1)
!! omegasp(nsp)=real frequencies for the imaginary part of the polarizability 
!! omegae(ne)= imaginary frequencies for the polarizability
!! delta=small imaginary part used to avoid poles, input variables
!!
!! OUTPUT
!! kkweight(nsp,ne)=frequency dependent weights (Eq A1 PRB 74, 035101 (2006)
!!
!! PARENTS
!!      spectral
!!
!! CHILDREN
!!      flush_unit,wrtout
!!
!! SOURCE
!!

subroutine calc_kkweight(ne,omegae,nsp,omegasp,delta,omegamax,kkw)

 use m_profiling
    
 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_kkweight'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ne,nsp
 real(dp),intent(in) :: delta,omegamax
!arrays
 real(dp),intent(in) :: omegasp(nsp)
 complex(dpc),intent(in) :: omegae(ne)
 complex(dpc),intent(out) :: kkw(nsp,ne)

!Local variables-------------------------------
!scalars
 integer :: isp,je
 real(dp) :: eta,xx1,xx2
 real(dp) :: den1,den2
 complex(dpc) :: c1,c2,wt
!************************************************************************

 DBG_ENTER("COLL")
 
 kkw(:,:)=czero

 do je=1,ne
   eta=delta
   wt=omegae(je)
   ! Not include shift at omega==0, what about metallic systems?
   if (abs(real(omegae(je)))<tol6 .and. abs(aimag(wt))<tol6) eta=tol12
   !  Not include shift along the imaginary axis
   if (abs(aimag(wt))>tol6) eta=zero
   do isp=1,nsp
     if (isp==1) then  
       ! Skip negative point, should check that this would not lead to spurious effects
       c1=czero
       den1=one
     else 
       xx1=omegasp(isp-1)
       xx2=omegasp(isp)
       den1= xx2-xx1 
       c1= -(wt-xx1+j_dpc*eta)*log( (wt-xx2+j_dpc*eta)/(wt-xx1+j_dpc*eta) )&
&          +(wt+xx1-j_dpc*eta)*log( (wt+xx2-j_dpc*eta)/(wt+xx1-j_dpc*eta) )
       c1= c1/den1
     end if 
     xx1=omegasp(isp)
     if (isp==nsp) then ! Skip last point should check that this would not lead to spurious effects 
       xx2=omegamax
     else
       xx2=omegasp(isp+1)
     end if
     den2=xx2-xx1
     c2=  (wt-xx2+j_dpc*eta)*log( (wt-xx2+j_dpc*eta)/(wt-xx1+j_dpc*eta) )&
&        -(wt+xx2-j_dpc*eta)*log( (wt+xx2-j_dpc*eta)/(wt+xx1-j_dpc*eta) )
     c2= c2/den2
     kkw(isp,je)=  c1/den1 + c2/den2
   end do
 end do 

 DBG_EXIT("COLL")

end subroutine calc_kkweight
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/setup_spectral
!! NAME
!!  setup_spectral
!!
!! FUNCTION
!! Calculation of \chi_o based on the spectral method as proposed in PRB 74, 035101 (2006) and PRB 61, 7172 (1999).
!! Setup of the real frequency mesh for $\Im\chi_o$ and of the frequency-dependent weights for 
!! Hilbert transform. Note that CPU time does not depend dramatically on nomegasf unlike memory.
!! spmeth defines the approximant for the delta function:
!!  ==1 : use Triangular approximant (Kresse method)
!!  ==2 : use Gaussian method, requiring smearing (Miyake method)
!!
!! INPUTS
!! nomegasf=number of points for the imaginary part of $\chi0(q,\omega)$ 
!! nomega=number of frequencies in $\chi0(q,\omega)$.
!! max_rest,min_res=max and min resonant transition energy (for this q-point)
!! my_max_rest,my_min_rest=max and min resonant transition energy treated by this processor
!! method=integer flag defining the type of frequency mesh used for $\Im chi0$
!!  | 0 for a linear mesh 
!!  | 1 for a mesh densified around omegaplasma
!! omegaplasma=frequency around which the mesh is densifies (usually Drude plasma frequency)
!!  used only in case of method==1
!! zcut=small imaginary shift to avoid pole in chi0
!!
!! OUTPUT
!!  kkweight(nomegasf,nomega)=Frequency dependent weight for Hilber transform.
!!  omegasf(nomegasf+1)=frequencies for imaginary part.
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!      flush_unit,wrtout
!!
!! SOURCE

subroutine setup_spectral(nomega,omega,nomegasf,omegasf,max_rest,min_rest,my_max_rest,my_min_rest,&
&  method,zcut,omegaplasma,my_wl,my_wr,kkweight)

 use m_profiling

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_spectral'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_70_gw, except_this_one => setup_spectral
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: method,nomega,nomegasf
 integer,intent(out) :: my_wl,my_wr
 real(dp),intent(in) :: max_rest,min_rest,omegaplasma,zcut
 real(dp),intent(in) :: my_max_rest,my_min_rest
!arrays
 real(dp),intent(out) :: omegasf(nomegasf)
 complex(dpc),intent(in) :: omega(nomega)       
 complex(dpc),intent(out) :: kkweight(nomegasf,nomega)

!Local variables-------------------------------
!scalars
 integer :: io,ii
 real(dp) :: nu_min,nu_max,nu1,nu2,dd,domegasf,wp,deltat
 character(len=500) :: msg                 
!arrays
 integer,allocatable :: insort(:)
!************************************************************************

 ! === The mesh has to enclose the entire range of transitions ===
 dd=(max_rest-min_rest)/(nomegasf-1)
 domegasf=(max_rest-min_rest+2*dd)/(nomegasf-1) 

 write(msg,'(4a,f8.3,3a,i5,2a,f8.5,a)')ch10,&
&  ' === Info on the real frequency mesh for spectral method === ',ch10,&
&  '  maximum frequency = ',max_rest*Ha_eV,' [eV]',ch10,&
&  '  nomegasf = ',nomegasf,ch10,&
&  '  domegasf = ',domegasf*Ha_eV,' [eV]'
 call wrtout(std_out,msg,'COLL') !; call wrtout(ab_out,msg,'COLL')

 if (min_rest<tol6) then  
   write(msg,'(a)')" system seems to be metallic"
   call wrtout(std_out,msg,'COLL') 
 end if
 !
 ! ======================================================
 ! === Setup of the w-mesh for the spectral function ====
 ! ======================================================
 SELECT CASE (method) 

 CASE (0) ! Linear mesh.
   call wrtout(std_out,' Using linear mesh for Im chi0','COLL')
   do io=1,nomegasf
     omegasf(io)=(io-1)*domegasf+min_rest-dd
   end do 

 CASE (1) ! Non-homogeneous mesh densified around omega_plasma, do not improve results ===
   ! WARNING_ this part has to be checked since I modified omegasf
   write(msg,'(a,f7.4,a)')' Using mesh densified around ',omegaplasma*Ha_eV,' [eV] '
   call wrtout(std_out,msg,'COLL')
   wp=omegaplasma ; deltat=max_rest-min_rest 
   nu_min=zero
   if (deltat<wp ) then 
     nu_max = wp/sqrt2 *   ATAN(sqrt2*deltat*wp/(-deltat**2+wp**2))
   else 
     nu_max = wp/sqrt2 * ( ATAN(sqrt2*deltat*wp/(-deltat**2+wp**2)) + pi)
   end if 
   domegasf=(nu_max-nu_min)/(nomegasf+1)
   !write(std_out,*)  -(wp/sqrt2) * atan(sqrt2*deltat*wp/(deltat**2-wp**2))
   omegasf(1)=zero ; omegasf(nomegasf+1)=deltat
   ii=0
   do io=2,nomegasf
     nu1=domegasf*(io-1) ; nu2=TAN(-sqrt2*nu1/wp)
     if (nu2<0) then 
       omegasf(io) = wp * (one - SQRT(1+2*nu2**2))/(sqrt2*nu2)
     else
       omegasf(io) = wp * (one + SQRT(1+2*nu2**2))/(sqrt2*nu2) 
     end if 
     if (omegasf(io)> deltat ) then  
       omegasf(io)= deltat-0.1*ii
       ii=ii+1
     end if 
     ! write(102,'(i4,2x,3(f9.4,2x))')io,nu1,nu2,ep%omegasf(io)*Ha_eV
   end do 

   ! Reorder frequencies in ascending order
   ABI_ALLOCATE(insort,(nomegasf+1))
   insort(:)=(/ (io,io=1,nomegasf+1) /)
   call sort_dp(nomegasf+1,omegasf,insort,tol14)
   ABI_DEALLOCATE(insort)

 CASE DEFAULT 
   write(msg,'(a,i4)')'Wrong value for method= ',method
   MSG_BUG(msg)
 END SELECT

 !write(std_out,*)omegasf(1)*Ha_eV,omegasf(nomegasf)*Ha_eV
 !
 ! === Find min and max index in omegasf treated by this processor ===
 my_wr=-999
 do io=1,nomegasf
   if (omegasf(io)>my_max_rest) then 
     my_wr=io; EXIT
   end if 
 end do
 if (my_wr==nomegasf+2) my_wr=nomegasf+1
 my_wl=-999
 do io=nomegasf,1,-1
   if (omegasf(io)< my_min_rest) then ! Check metals
     my_wl=io; EXIT
   end if 
 end do 

 write(msg,'(a,2i6)')' my_wl and my_wr:',my_wl,my_wr
 call wrtout(std_out,msg,'PERS')

 if (my_wl==-999 .or. my_wr==-999) then 
   write(msg,'(a,2i6)')' wrong value in my_wl and/or my_wr ',my_wl,my_wr
   MSG_ERROR(msg)
 end if 
 ! 
 ! Calculate weights for Hilbert transform.
 call calc_kkweight(nomega,omega,nomegasf,omegasf,zcut,max_rest,kkweight)

end subroutine setup_spectral
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/hilbert_transform
!! NAME
!!  hilbert_transform
!!
!! FUNCTION
!!
!! INPUTS
!! nomegasf=number of points for the imaginary part of $\chi0(q,\omega)$ 
!! nomega=number of frequencies in $\chi0(q,\omega)$.
!! max_rest,min_res=max and min resonant transition energy (for this q-point)
!! my_max_rest,my_min_rest=max and min resonant transition energy treated by this processor
!!
!! OUTPUT
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!      flush_unit,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine hilbert_transform(npwe,nomega,nomegasf,my_wl,my_wr,kkweight,sf_chi0,chi0,spmeth)

 use m_profiling

 use defs_basis
 use m_errors
 use m_io_tools, only : flush_unit
#if HAVE_GW_OPENMP
!$ use omp_lib, only: omp_get_thread_num, omp_get_num_threads
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hilbert_transform'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spmeth,nomega,nomegasf
 integer,intent(in) :: my_wl,my_wr,npwe
! real(dp),intent(in) :: max_rest,min_rest
! real(dp),intent(in) :: my_max_rest,my_min_rest
! real(dp),intent(in) :: max_fit_rest
!arrays
! real(dp),intent(in) :: omegasf(nomegasf)
 complex(dpc),intent(in) :: kkweight(nomegasf,nomega)
 complex(gwpc), intent(inout) :: sf_chi0(npwe,npwe,my_wl:my_wr)
 complex(gwpc), intent(inout) :: chi0(npwe,npwe,nomega)

!Local variables-------------------------------
!scalars
 integer :: ig2,io,iw
#if HAVE_GW_OPENMP
 integer :: ig1,local_std_out
#endif
 complex(gwpc) :: kkw
 character(len=500) :: msg                 
!arrays
!************************************************************************

#if HAVE_GW_OPENMP
 write(msg,'(2a,i3,a)')ch10,&
&    ' Performing Hilbert transform (with OpenMP) using method ',spmeth,' It might take some time...'
 call wrtout(std_out,msg,'COLL')
#else
 write(msg,'(2a,i3,a)')ch10,&
&    ' Performing Hilbert transform using method ',spmeth,' It might take some time...'
 call wrtout(std_out,msg,'COLL')
#endif
 call flush_unit(std_out) 
 
 ! First coding, The loop over ig1, ig2 could be optimised taking into account symmetries

#if HAVE_GW_OPENMP
 local_std_out = std_out
!!OMP *** OPENMP SECTION *** Added by MS
!$OMP PARALLEL SHARED(sf_chi0,chi0,kkweight,my_wl,my_wr,npwe,nomega,local_std_out) &
!$OMP          PRIVATE(io,iw,ig1,ig2,kkw)
!  !$  if (prtvol>9.and.omp_get_thread_num()==0) then
!  !$    write(local_std_out,'(a,i0)') ' Entering openmp loop. Number of threads: ',omp_get_num_threads()
!  !$    call flush_unit(local_std_out)
!  !$  end if
!$OMP DO SCHEDULE(DYNAMIC,1)
 do io=1,nomega
!  !$   if (prtvol>9) then
!  !$     if (mod(io,10)==0) write(local_std_out,'(3(a,i0))') &
!  !$     ' Processor ',omp_get_thread_num(),' doing ',io,' of ',nomega
!  !$     call flush_unit(local_std_out)
!  !$   end if
   do iw=my_wl,my_wr
     kkw = kkweight(iw,io)
     do ig2=1,npwe
       do ig1=1,npwe
        chi0(ig1,ig2,io) = chi0(ig1,ig2,io) + kkw*sf_chi0(ig1,ig2,iw)
       end do
     end do
   end do
 end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
#else
 do io=1,nomega
   do iw=my_wl,my_wr
     kkw = kkweight(iw,io)
     do ig2=1,npwe
#if defined HAVE_GW_DPC
      CALL zaxpy(npwe,kkw,sf_chi0(1:npwe,ig2,iw),1,chi0(1:npwe,ig2,io),1)
#else
      CALL caxpy(npwe,kkw,sf_chi0(1:npwe,ig2,iw),1,chi0(1:npwe,ig2,io),1)
#endif
     end do
   end do
 end do
#endif

end subroutine hilbert_transform
!!***

!!****f* ABINIT/hilbert_transform_headwings
!! NAME
!!  hilbert_transform_headwings
!!
!! FUNCTION
!!
!! INPUTS
!! nomegasf=number of points for the imaginary part of $\chi0(q,\omega)$ 
!! nomega=number of frequencies in $\chi0(q,\omega)$.
!! max_rest,min_res=max and min resonant transition energy (for this q-point)
!! my_max_rest,my_min_rest=max and min resonant transition energy treated by this processor
!!
!! OUTPUT
!!
!! PARENTS
!!      cchi0q0
!!
!! CHILDREN
!!      flush_unit,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine hilbert_transform_headwings(npwe,nomega,nomegasf,my_wl,my_wr,kkweight, &
& sf_lwing,sf_uwing,sf_head,chi0_lwing,chi0_uwing,chi0_head,spmeth)

 use m_profiling

 use defs_basis
 use m_errors
 use m_io_tools, only : flush_unit
! use m_blas,     only : xaxpy
#if HAVE_GW_OPENMP
!$ use omp_lib, only: omp_get_thread_num, omp_get_num_threads
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hilbert_transform_headwings'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spmeth,nomega,nomegasf
 integer,intent(in) :: my_wl,my_wr,npwe
! real(dp),intent(in) :: max_rest,min_rest
! real(dp),intent(in) :: my_max_rest,my_min_rest
! real(dp),intent(in) :: max_fit_rest
!arrays
! real(dp),intent(in) :: omegasf(nomegasf)
 complex(dpc),intent(in) :: kkweight(nomegasf,nomega)
 complex(dpc), intent(inout) :: sf_lwing(npwe,my_wl:my_wr,3)
 complex(dpc), intent(inout) :: sf_uwing(npwe,my_wl:my_wr,3)
 complex(dpc), intent(inout) :: sf_head(3,3,my_wl:my_wr)
 complex(dpc), intent(inout) :: chi0_lwing(npwe,nomega,3)
 complex(dpc), intent(inout) :: chi0_uwing(npwe,nomega,3)
 complex(dpc), intent(inout) :: chi0_head(3,3,nomega)

!Local variables-------------------------------
!scalars
 integer :: ig1,idir,io,iw
#if HAVE_GW_OPENMP
 integer :: local_std_out
#endif
 complex(gwpc) :: kkw 
 character(len=500) :: msg                 
!arrays
! complex(gwpc) :: sf_lwing_slice(my_wr-my_wl+1)
! complex(gwpc) :: sf_uwing_slice(my_wr-my_wl+1)
!************************************************************************

#if HAVE_GW_OPENMP
 write(msg,'(2a,i3,a)')ch10,&
&    ' Performing Hilbert transform (with OpenMP) using method ',spmeth,' It might take some time...'
 call wrtout(std_out,msg,'COLL')
#else
 write(msg,'(2a,i3,a)')ch10,&
&    ' Performing Hilbert transform using method ',spmeth,' It might take some time...'
 call wrtout(std_out,msg,'COLL')
#endif
 call flush_unit(std_out) 
 
 ! First coding, The loop over ig1, ig2 could be optimised taking into account symmetries

#if HAVE_GW_OPENMP
 local_std_out = std_out
!!OMP  *** OPENMP SECTION *** Added by MS
!$OMP PARALLEL SHARED(sf_lwing,sf_uwing,sf_head,chi0_lwing,chi0_uwing,chi0_head,kkweight,&
!$OMP                 my_wl,my_wr,npwe,nomega,local_std_out) &
!$OMP          PRIVATE(io,ig1,iw,kkw)
!  !$  if (prtvol>9.and.omp_get_thread_num()==0) then
!  !$    write(local_std_out,'(a,i0)') ' Entering openmp loop. Number of threads: ',omp_get_num_threads()
!  !$    call flush_unit(local_std_out)
!  !$  end if
 do idir=1,3
!$OMP DO SCHEDULE(DYNAMIC,1)
   do io=1,nomega
     do iw=my_wl,my_wr 
       kkw = kkweight(iw,io)
       do ig1=1,npwe
         chi0_lwing(ig1,io,idir) = chi0_lwing(ig1,io,idir) + kkw*sf_lwing(ig1,iw,idir)
         chi0_uwing(ig1,io,idir) = chi0_uwing(ig1,io,idir) + kkw*sf_uwing(ig1,iw,idir)
       end do
     end do
   end do
!$OMP END DO
 end do
 !
 ! Hilbert transform of the head.
!$OMP DO SCHEDULE(DYNAMIC,1)
 do io=1,nomega
     chi0_head(1,1,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(1,1,my_wl:my_wr))
     chi0_head(2,1,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(2,1,my_wl:my_wr))
     chi0_head(3,1,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(3,1,my_wl:my_wr))
     chi0_head(1,2,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(1,2,my_wl:my_wr))
     chi0_head(2,2,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(2,2,my_wl:my_wr))
     chi0_head(3,2,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(3,2,my_wl:my_wr))
     chi0_head(1,3,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(1,3,my_wl:my_wr))
     chi0_head(2,3,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(2,3,my_wl:my_wr))
     chi0_head(3,3,io) = SUM(kkweight(my_wl:my_wr,io)*sf_head(3,3,my_wl:my_wr))
 end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
#else
 ! * Hilbert transform for wings.
 ! Partial contributions to chi0 will be summed afterwards.
 do idir=1,3
   do io=1,nomega
     do iw=my_wl,my_wr
       kkw = kkweight(iw,io)
       do ig1=1,npwe
         chi0_lwing(ig1,io,idir) = chi0_lwing(ig1,io,idir) + kkw*sf_lwing(ig1,iw,idir)
         chi0_uwing(ig1,io,idir) = chi0_uwing(ig1,io,idir) + kkw*sf_uwing(ig1,iw,idir)
       end do  
     end do
   end do
 end do  !idir
 !
 ! Hilbert transform of the head.
 do io=1,nomega
   do iw=my_wl,my_wr
     kkw = kkweight(iw,io)
     chi0_head(1,1,io) = chi0_head(1,1,io) + kkw*sf_head(1,1,iw)
     chi0_head(2,1,io) = chi0_head(2,1,io) + kkw*sf_head(2,1,iw)
     chi0_head(3,1,io) = chi0_head(3,1,io) + kkw*sf_head(3,1,iw)
     chi0_head(1,2,io) = chi0_head(1,2,io) + kkw*sf_head(1,2,iw)
     chi0_head(2,2,io) = chi0_head(2,2,io) + kkw*sf_head(2,2,iw)
     chi0_head(3,2,io) = chi0_head(3,2,io) + kkw*sf_head(3,2,iw)
     chi0_head(1,3,io) = chi0_head(1,3,io) + kkw*sf_head(1,3,iw)
     chi0_head(2,3,io) = chi0_head(2,3,io) + kkw*sf_head(2,3,iw)
     chi0_head(3,3,io) = chi0_head(3,3,io) + kkw*sf_head(3,3,iw)
   end do
 end do
#endif

end subroutine hilbert_transform_headwings
!!***

!----------------------------------------------------------------------
