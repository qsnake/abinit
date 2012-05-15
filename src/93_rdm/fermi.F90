!{\src2tex{textfont=tt}}
!!****f* ABINIT/fermi
!! NAME
!! fermi
!!
!! FUNCTION
!! Calculate the Fermi level and occupation numbers
!! (This routine should be cleaned !!)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  en(nkibz,nbnds,nsppol)=KS energies, for different k points, bands and spins
!!  fixmom=if differ from -99.99d0, fix the magnetic moment (in Bohr magneton)
!!  Hdr=header of previously read file, containing many variables
!!  nbnds=number of bands
!!  nkibz=number of irred k-points
!!  nsppol=number of spins (should be changed to nsppol)
!!  wtk(nkibz)=weights for k points (input variable)
!!  occ(nkibz,nbnds,nsppol)= occupations numbers for each k-point, band and spin
!!
!! OUTPUT
!!  nbv(nsppol)=number of valence bands (for each spin)
!!  fermie=Fermi energy
!!
!! SIDE EFFECTS
!!  nel= at input, if nel==0, recalculate nelect from occupation numbers
!!  occ(nkibz,nbnds,nsppol)=occupation numbers, for each k point, each band and spin
!!
!! PARENTS
!!      rdm
!!
!! CHILDREN
!!      newocc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine fermi(hdr,nbnds,nkibz,fixmom,nsppol,wtk,en,occ,nel,nbv,fermie)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fermi'
 use interfaces_14_hidewrite
 use interfaces_62_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbnds,nkibz,nsppol
 integer,intent(inout) :: nel
 real(dp),intent(in) :: fixmom
 real(dp),intent(out) :: fermie
 type(Hdr_type),intent(in) :: Hdr
!arrays
 integer,intent(out) :: nbv(nsppol)
 real(dp),intent(in) :: en(nkibz,nbnds,nsppol),wtk(nkibz)
 real(dp),intent(inout) :: occ(nkibz,nbnds,nsppol)

!Local variables-------------------------------
!scalars
 integer :: ib,idx,ik,is
 real(dp),parameter :: eps5=5*1.d-2
 real(dp) :: cbot,diff,entropy,nelect,vtop
 character(len=500) :: msg
!arrays
 integer,allocatable :: nbandt(:)
 real(dp) :: condbottom(nsppol),n(nsppol),valencetop(nsppol),wtk_norm(nkibz)
 real(dp),allocatable :: doccdet(:),eigent(:),occt(:)

! *************************************************************************

!=== Get number of electrons for each channel ===
 wtk_norm(:)=wtk(:)/SUM(wtk(:))
 n(:)=zero
 do is=1,nsppol
   do ik=1,nkibz
     do ib=1,nbnds
       n(is)= n(is)+wtk_norm(ik)*occ(ik,ib,is)
     end do
   end do
 end do

 if (nel==0) then 
!  === Avoid problem in newocc, if fermi is called with nel==0 ===
   nel=NINT(SUM(n))
   nelect=nel
   write(msg,'(2a,i5)')ch10,' fermi : total number of electrons = ',nel
   call wrtout(std_out,msg,'COLL')
 else
   nelect=nel
   write(msg,'(2a,i5)')ch10,' fermi : total number of electrons = ',nel
   call wrtout(std_out,msg,'COLL')
!  Check the input value nel with the calculated one.
!  Warn if the difference is small, stop if too large 
   diff=nelect-SUM(n)
   write(msg,'(2a,f12.8)')ch10,&
&   ' fermi : input and calculated number of electrons differ by ',diff 
   call wrtout(std_out,msg,'COLL')
   if (ABS(diff)>eps5*nelect) then 
     write(msg,'(4a)')ch10,&
&     ' fermi : ERROR - ',ch10,&
&     '  Found too large difference in the number of electrons'
     call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
   end if
 end if  

 if (Hdr%nspinor==1) then
   if (nsppol==1) then 
     nbv(:)=(nel+1)/2
   else 
     write(msg,'(2a,5x,2(f8.4,2x))')&
&     ' number of electrons for each spin channel = ',ch10,(n(is),is=1,nsppol)
     call wrtout(std_out,msg,'COLL')
!    Take into account the case in which there are two bands with 
!    the same spin index close to the gap
     nbv(:)=NINT(n(:))
     write(msg,'(a,2i8)')' maximum occupied band index for spin up and down= ',nbv(:)
     call wrtout(std_out,msg,'COLL')
   end if
 else 
   nbv(:)=nel
 end if

 if (Hdr%occopt<3.or.Hdr%occopt>8) then
!  
!  === Semiconductor or Insulator ===
!  * Calculate valence index for each channel
   do is=1,nsppol
     valencetop(is)= -huge(0.0_dp)
     condbottom(is)=  huge(0.0_dp)
     do ik=1,nkibz
       do ib=1,nbv(is)
         if (valencetop(is)<en(ik,ib,is)) valencetop(is)=en(ik,ib,is)
       end do
       do ib=nbv(is)+1,nbnds
         if (condbottom(is)>en(ik,ib,is)) condbottom(is)=en(ik,ib,is)
       end do
     end do 
   end do 
   vtop=MAXVAL(valencetop)
   cbot=MINVAL(condbottom)
   write(msg,'(a,f6.2,2a,f6.2)')&
&   ' top of valence       [eV] ',vtop*Ha_eV,ch10,&
&   ' bottom of conduction [eV] ',cbot*Ha_eV
   call wrtout(std_out,msg,'COLL')
   if (nsppol==2) then 
     if (ABS(vtop-MINVAL(valencetop))>tol6) then 
       write(msg,'(a,i2)')' top of valence is spin ',MAXLOC(valencetop)
       call wrtout(std_out,msg,'COLL')
     end if
     if (ABS(cbot-MAXVAL(condbottom))>tol6) then 
       write(msg,'(a,i2)')' bottom of conduction is spin ',MINLOC(condbottom)
       call wrtout(std_out,msg,'COLL')
     end if
   end if
   fermie=(vtop+cbot)/2 
   if (ABS(cbot-vtop)<1.d-4) fermie= vtop !to avoid error on the last digit
 else 
!  
!  === For Metallic system. Recalculate occupations factors ===
   ABI_ALLOCATE(eigent,(nbnds*nkibz*nsppol))
   ABI_ALLOCATE(nbandt,(nkibz*nsppol))
   eigent(:)=zero ; nbandt(:)=nbnds
   idx=0
   do is=1,nsppol
     do ik=1,nkibz
       do ib=1,nbnds
         idx=idx+1
         eigent(idx)=en(ik,ib,is)
       end do
     end do
   end do
   write(msg,'(a,f9.5)')' fermi : metallic system, calling newocc with fixmom = ',fixmom
   call wrtout(std_out,msg,'COLL')
   ABI_ALLOCATE(occt,(nbnds*nkibz*nsppol))
   ABI_ALLOCATE(doccdet,(nbnds*nkibz*nsppol))

   call newocc(doccdet,eigent,entropy,fermie,fixmom,nbnds,nbandt,nelect,nkibz,&
&   Hdr%nspinor,nsppol,occt,Hdr%occopt,1,Hdr%stmbias,Hdr%tphysel,Hdr%tsmear,wtk_norm)
   
   idx=0
   do is=1,nsppol
     do ik=1,nkibz
       do ib=1,nbnds
         idx=idx+1
         if (ABS(occt(idx)-occ(ik,ib,is))>1.d-2) then
           write(msg,'(3a,i5,a,i5,a,i3,2a,f6.3)')&
&           ' fermi: WARNING - ',ch10,&
&           '  Occupation number of band ',ib,' for kpt ',ik,' spin ',is,ch10,&
           '  is changed with respect to KSS file : ',occt(idx)-occ(ik,ib,is)
           call wrtout(std_out,msg,'COLL')
         end if
         occ(ik,ib,is)=occt(idx)
       end do
     end do
   end do
   ABI_DEALLOCATE(eigent)
   ABI_DEALLOCATE(nbandt)
   ABI_DEALLOCATE(occt)
   ABI_DEALLOCATE(doccdet)
 end if

 write(msg,'(a,f6.2,a)')' Fermi energy         [eV] ',fermie*Ha_eV,ch10
 call wrtout(std_out,msg,'COLL')
 
end subroutine fermi
!!***
