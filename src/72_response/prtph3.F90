!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtph3
!!
!! NAME
!! prtph3
!!
!! FUNCTION
!! Print the phonon frequencies, on unit 6 as well as the printing
!! unit (except if the associated number -iout- is negative),
!! and for the latter, in Hartree, meV, Thz, Kelvin or cm-1.
!! If eivec==1,2, also print the eigenmodes : displacements
!! in cartesian coordinates.
!! If eivec==3, print the dynamic information to the lwf-formatted file.
!! If eivec==4, generate output files for band2eps (drawing tool for the
!!               phonon band structure

!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  displ(2,3*natom,3*natom)= contain
!!  the displacements of atoms in cartesian coordinates.
!!  The first index means either the real or the imaginary part,
!!  The second index runs on the direction and the atoms displaced
!!  The third index runs on the modes.
!!  eivec=(if eivec==0, the eigendisplacements are not printed,
!!    if eivec==1,2, the eigendisplacements are printed,
!!    if eivec==3, the eigendisplacements are printed in lwf-output format file)
!!    if eivec==4, files for band2eps
!!  enunit=units for output of the phonon frequencies :
!!    0=> Hartree and cm-1, 1=> eV and Thz, other=> Ha,Thz,eV,cm-1 and K
!!  iodyn=unit number for the LWF output (if negative, do not write)
!!  iout= unit for long print (if negative, the routine only print
!!        on unit 6, and in Hartree only).
!!  natom= number of atom
!!  phfreq(3*natom)= phonon frequencies in Hartree
!!  qphnrm=phonon wavevector normalisation factor
!!  qphon(3)=phonon wavevector
!!
!! OUTPUT
!!
!! NOTES
!! called by one processor only
!!
!! PARENTS
!!      anaddb,mkifc9,mkphbs,respfn
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prtph3(displ,eivec,enunit,iodyn,iout,natom,phfrq,qphnrm,qphon)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtph3'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: eivec,enunit,iodyn,iout,natom
 real(dp),intent(in) :: qphnrm
!arrays
 real(dp),intent(in) :: displ(2,3*natom,3*natom),phfrq(3*natom),qphon(3)

!Local variables -------------------------
!scalars
 integer :: i,idir,ii,imode,jj
 real(dp) :: tolerance
 character(len=500) :: message
!arrays
 real(dp) :: vecti(3),vectr(3)

! *********************************************************************

!DEBUG
!write(std_out,'(a)' )' prtph3 : enter '
!ENDDEBUG

!Check the value of eivec
 if(eivec<0.or.eivec>4)then
   write(message, '(a,a,a,i6,a,a)' )&
&   ' prtph3 : BUG -',ch10,&
&   '  In the calling subroutine, eivec is',eivec,ch10,&
&   '  but allowed values are between 0 and 4.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!write the phonon frequencies on unit 6
 write(std_out,'(a)' )' '
 write(std_out,'(a,a)' )' phonon wavelength (reduced coordinates) , ', &
& 'norm, and energies in hartree'
!The next format should be rewritten
 write(std_out,'(a,4f5.2)' )&
& ' ',(qphon(i),i=1,3),qphnrm
 do jj=1,3*natom,5
   if (3*natom-jj<5) then
     write(std_out,'(5es17.9)') (phfrq(ii),ii=jj,3*natom)
   else
     write(std_out,'(5es17.9)') (phfrq(ii),ii=jj,jj+4)
   end if
 end do
 write(std_out,'(a,es17.9)')' Zero Point Motion energy (sum of freqs/2)=',sum(phfrq(1:3*natom))/2

!DEBUG
!write(std_out,'(a,2es17.9)')' Prepared by dividing by 0.0009375, Ha, eV=',&
!&   sum(phfrq(1:3*natom))/2/0.0009375,&
!&   sum(phfrq(1:3*natom))/2/0.0009375*27.211
!ENDDEBUG

!Put the wavevector in nice format
 if(iout>=0)then
   write(iout, '(a)' )' '
   if(qphnrm/=0.0_dp)then
     write(iout, '(a,3f9.5)' )&
&     '  Phonon wavevector (reduced coordinates) :',&
&     (qphon(i)/qphnrm+tol10,i=1,3)
   else
     write(iout, '(a,/,a,3f9.5)' )&
&     '  Phonon at Gamma, with non-analyticity in the',&
&     '  direction (cartesian coordinates)',qphon(1:3)+tol10
   end if

!  Write it, in different units.
   if(enunit/=1)then
     write(iout, '(a)' )' Phonon energies in Hartree :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(iout, '(1x,5es14.6)') (phfrq(ii),ii=jj,3*natom)
       else
         write(iout, '(1x,5es14.6)') (phfrq(ii),ii=jj,jj+4)
       end if
     end do
   end if
   if(enunit/=0)then
     write(iout, '(a)' )' Phonon energies in meV     :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(iout, '("-",5es14.6)') (phfrq(ii)*Ha_eV*1.0d3,ii=jj,3*natom)
       else
         write(iout, '("-",5es14.6)') (phfrq(ii)*Ha_eV*1.0d3,ii=jj,jj+4)
       end if
     end do
   end if
   if(enunit/=1)then
     write(iout, '(a)' )' Phonon frequencies in cm-1    :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(iout, '("-",5es14.6)') (phfrq(ii)*Ha_cmm1,ii=jj,3*natom)
       else
         write(iout, '("-",5es14.6)') (phfrq(ii)*Ha_cmm1,ii=jj,jj+4)
       end if
     end do
   end if
   if(enunit/=0)then
     write(iout, '(a)' )' Phonon frequencies in Thz     :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(iout, '("-",5es14.6)') (phfrq(ii)*Ha_THz,ii=jj,3*natom)
       else
         write(iout, '("-",5es14.6)') (phfrq(ii)*Ha_THz,ii=jj,jj+4)
       end if
     end do
   end if
   if(enunit/=0.and.enunit/=1)then
     write(iout, '(a)' )' Phonon energies in Kelvin  :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(iout, '("-",5es14.6)') (phfrq(ii)/kb_HaK,ii=jj,3*natom)
       else
         write(iout, '("-",5es14.6)') (phfrq(ii)/kb_HaK,ii=jj,jj+4)
       end if
     end do
   end if
 end if

!Write the lwf-formatted file
!wavevector position
 if (eivec==3) then
   if (iodyn>0) then
     write(iodyn, '(a,a,3f9.5)' )&
&     '  Phonon_wavevector_(reduced_coordinates)',ch10,&
&     (qphon(i)/qphnrm+tol10,i=1,3)
     write(iodyn, '(a)' )' Phonon_frequencies_in_cm-1'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(iodyn, '(5f15.9)') (phfrq(ii)*Ha_cmm1,ii=jj,3*natom)
       else
         write(iodyn, '(5f15.9)') (phfrq(ii)*Ha_cmm1,ii=jj,jj+4)
       end if
     end do
   end if
 end if

!Take care of the eigendisplacements
 if(eivec==1 .or. eivec==2)then
   write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&   ' Eigendisplacements ',ch10,&
&   ' (will be given, for each mode : in cartesian coordinates',ch10,&
&   '   for each atom the real part of the displacement vector,',ch10,&
&   '   then the imaginary part of the displacement vector)'
   call wrtout(std_out,message,'COLL')
   if(iout>=0)then
     call wrtout(iout,message,'COLL')
   end if
   do imode=1,3*natom
     write(std_out,'(a,i4,a,es16.6)' )&
&     '  Mode number ',imode,'   Energy',phfrq(imode)
     if(iout>=0)then
       write(iout, '(a,i4,a,es16.6)' )&
&       '  Mode number ',imode,'   Energy',phfrq(imode)
     end if
     tolerance=1.0d-6
     if(abs(phfrq(imode))<1.0d-5)tolerance=2.0d-6
     if(phfrq(imode)<1.0d-5)then
       write(std_out,'(a)' )' Attention : low frequency mode.'
       write(std_out,'(a)' )'   (Could be unstable or acoustic mode)'
       if(iout>=0)then
         write(iout, '(a)' )' Attention : low frequency mode.'
         write(iout, '(a)' )'   (Could be unstable or acoustic mode)'
       end if
     end if
     do ii=1,natom
       do idir=1,3
         vectr(idir)=displ(1,idir+(ii-1)*3,imode)
         if(abs(vectr(idir))<tolerance)vectr(idir)=0.0_dp
         vecti(idir)=displ(2,idir+(ii-1)*3,imode)
         if(abs(vecti(idir))<tolerance)vecti(idir)=0.0_dp
       end do
       write(std_out,'(i4,3es16.8)' ) ii,vectr(:)
       write(std_out,'(4x,3es16.8)' )    vecti(:)
       if(iout>=0)then
         write(iout,'(a,i3,3es16.8)') ';',ii,vectr(:)
         write(iout,'(a,3x,3es16.8)') ';',   vecti(:)
       end if
     end do
   end do
 end if

!Write the lwf-formatted file (if iodyn>0)
!eigenvectors
 if (eivec==3) then
   if (iodyn>0) then
     do imode=1,3*natom
       write(iodyn, '(a,i4,a,f14.9)' )&
&       '  Mode_number ',imode,' Frequency',phfrq(imode)*Ha_cmm1
       do ii=1,natom
         do idir=1,3
           vectr(idir)=displ(1,idir+(ii-1)*3,imode)
           if(abs(vectr(idir))<tolerance)vectr(idir)=0.0_dp
           vecti(idir)=displ(2,idir+(ii-1)*3,imode)
           if(abs(vecti(idir))<tolerance)vecti(idir)=0.0_dp
         end do
         write(iodyn,'(i5,3f18.12)') ii,vectr(:)
         write(iodyn,'(i5,3f18.12)') ii,vecti(:)
       end do
     end do
   end if
 end if

!DEBUG
!write(std_out,'(a)' )' prtph3 : exit '
!ENDDEBUG

end subroutine prtph3
!!***
