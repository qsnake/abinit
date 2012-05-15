!{\src2tex{textfont=tt}}
!!****f* ABINIT/outvars9
!!
!! NAME
!! outvars9
!!
!! FUNCTION
!! Open input file for the anaddb code, then
!! echoes the input information.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG,JCC,CL,XW)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! anaddb_dtset= (derived datatype) contains all the input variables
!! nunit=unit number for input or output
!!
!! OUTPUT
!!  (only writing)
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine outvars9 (anaddb_dtset,nunit)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outvars9'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nunit
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset

!Local variables -------------------------
!Set routine version number here:
!scalars
 integer :: ii,iph1,iph2,iqpt,iqshft

!*********************************************************************

!Write the heading
 write(nunit,'(a,80a,a)') ch10,('=',ii=1,80),ch10
 write(nunit, '(a,a)' )&
& ' -outvars9: echo values of input variables ----------------------',ch10

!The flags
 if(anaddb_dtset%dieflag/=0 .or. anaddb_dtset%ifcflag/=0 .or. &
& anaddb_dtset%nlflag/=0 .or. anaddb_dtset%thmflag/=0 .or. &
& anaddb_dtset%elaflag/=0 .or. anaddb_dtset%elphflag/=0 .or. &
& anaddb_dtset%polflag/=0 .or. anaddb_dtset%instrflag/=0 .or. &
& anaddb_dtset%piezoflag/=0                                   &
& )then
   write(nunit,'(a)')' Flags :'
   if(anaddb_dtset%dieflag/=0)write(nunit,'(3x,a9,3i10)')'  dieflag',anaddb_dtset%dieflag
   if(anaddb_dtset%ifcflag/=0)write(nunit,'(3x,a9,3i10)')'  ifcflag',anaddb_dtset%ifcflag
   if(anaddb_dtset%nlflag/=0)write(nunit,'(3x,a9,3i10)')'   nlflag',anaddb_dtset%nlflag
   if(anaddb_dtset%thmflag/=0)write(nunit,'(3x,a9,3i10)')'  thmflag',anaddb_dtset%thmflag
   if(anaddb_dtset%elaflag/=0)write(nunit,'(3x,a9,3i10)')'  elaflag',anaddb_dtset%elaflag
   if(anaddb_dtset%elphflag/=0)write(nunit,'(3x,a9,3i10)')' elphflag',anaddb_dtset%elphflag
   if(anaddb_dtset%polflag/=0)write(nunit,'(3x,a9,3i10)')'  polflag',anaddb_dtset%polflag
   if(anaddb_dtset%instrflag/=0)write(nunit,'(3x,a9,3i10)')'instrflag',anaddb_dtset%instrflag
   if(anaddb_dtset%piezoflag/=0)write(nunit,'(3x,a9,3i10)')'piezoflag',anaddb_dtset%piezoflag
 end if

!Write the general information
 if( anaddb_dtset%rfmeth/=1 .or. &
& anaddb_dtset%enunit/=0 .or. &
& anaddb_dtset%eivec/=0 .or. &
& anaddb_dtset%asr/=0 .or. &
& anaddb_dtset%chneut/=0 .or. &
& anaddb_dtset%selectz/=0        )then
   write(nunit,'(a)')' Miscellaneous information :'
   if(anaddb_dtset%rfmeth/=1)write(nunit,'(3x,a9,3i10)')'   rfmeth',anaddb_dtset%rfmeth
   if(anaddb_dtset%enunit/=0)write(nunit,'(3x,a9,3i10)')'   enunit',anaddb_dtset%enunit
   if(anaddb_dtset%eivec/=0) write(nunit,'(3x,a9,3i10)')'    eivec',anaddb_dtset%eivec
   if(anaddb_dtset%asr/=0)   write(nunit,'(3x,a9,3i10)')'      asr',anaddb_dtset%asr
   if(anaddb_dtset%chneut/=0)write(nunit,'(3x,a9,3i10)')'   chneut',anaddb_dtset%chneut
   if(anaddb_dtset%selectz/=0)write(nunit,'(3x,a9,3i10)')'  selectz',anaddb_dtset%selectz
 end if

!Frequency information
 if(anaddb_dtset%dieflag==1)then
   write(nunit,'(a)')' Frequency information :'
   write(nunit,'(3x,a9,3i10)')'    nfreq',anaddb_dtset%nfreq
   write(nunit,'(3x,a9,7x,3es16.8)')'    frmin',anaddb_dtset%frmin
   write(nunit,'(3x,a9,7x,3es16.8)')'    frmax',anaddb_dtset%frmax
 end if

!For interatomic force constant information
 if(anaddb_dtset%ifcflag/=0)then
   write(nunit,'(a)')' Interatomic Force Constants Inputs :'
   write(nunit,'(3x,a9,3i10)')'   dipdip',anaddb_dtset%dipdip
   if(anaddb_dtset%nsphere/=0)write(nunit,'(3x,a9,3i10)')'  nsphere',anaddb_dtset%nsphere
   if(abs(anaddb_dtset%rifcsph)>tol10)write(nunit,'(3x,a9,E16.6)')'  nsphere',anaddb_dtset%rifcsph
   write(nunit,'(3x,a9,3i10)')'   ifcana',anaddb_dtset%ifcana
   write(nunit,'(3x,a9,3i10)')'   ifcout',anaddb_dtset%ifcout
   if(anaddb_dtset%natifc>=1)then
     write(nunit,'(3x,a9,3i10)')'   natifc',anaddb_dtset%natifc
     write(nunit,'(3x,a9,8i10)')'    atifc',(anaddb_dtset%atifc(ii),ii=1,anaddb_dtset%natifc)
   end if
   write(nunit,'(a)')' Description of grid 1 :'
   write(nunit,'(3x,a9,3i10)')'     brav',anaddb_dtset%brav
   write(nunit,'(3x,a9,3i10)')'    ngqpt',anaddb_dtset%ngqpt(1:3)
   write(nunit,'(3x,a9,3i10)')'   nqshft',anaddb_dtset%nqshft
   if (anaddb_dtset%nqshft/=0)then
     write(nunit,'(3x,a9)')'   q1shft'
     do iqshft=1,anaddb_dtset%nqshft
       write(nunit,'(19x,4es16.8)') (anaddb_dtset%q1shft(ii,iqshft),ii=1,3)
     end do
   end if
   if (anaddb_dtset%qrefine > 1) then
     write(nunit,'(3x,a9,i10)')'  qrefine', anaddb_dtset%qrefine
   end if
 end if

!Phonon density of states with gaussian method
 if(anaddb_dtset%prtdos/=0)then
   write(nunit,'(a)')' Phonon DOS information :'
   write(nunit,'(3x,a9,es16.8)')'dosdeltae',anaddb_dtset%dosdeltae
   write(nunit,'(3x,a9,es16.8)')' dossmear',anaddb_dtset%dossmear
   write(nunit,'(a)')' Description of grid 2 for Fourier interpolation :'
   write(nunit,'(3x,a9,3i10)')'   ng2qpt',anaddb_dtset%ng2qpt(1:3)
   write(nunit,'(3x,a9,3i10)')'   ngrids',anaddb_dtset%ngrids
   write(nunit,'(3x,a9,7x,3es16.8)')'   q2shft',anaddb_dtset%q2shft(1:3)
 end if

!Thermal information
 if(anaddb_dtset%thmflag/=0)then
   write(nunit,'(a)')' Thermal information :'
   write(nunit,'(3x,a9,3i10)')'    nchan',anaddb_dtset%nchan
   write(nunit,'(3x,a9,3i10)')'   nwchan',anaddb_dtset%nwchan
   write(nunit,'(3x,a9,7x,3es16.8)')'   dostol',anaddb_dtset%dostol
   write(nunit,'(3x,a9,7x,3es16.8)')'   thmtol',anaddb_dtset%thmtol
   write(nunit,'(3x,a9,3i10)')'  ntemper',anaddb_dtset%ntemper
   write(nunit,'(3x,a9,7x,3es16.8)')'temperinc',anaddb_dtset%temperinc
   write(nunit,'(3x,a9,7x,3es16.8)')'tempermin',anaddb_dtset%tempermin
   if (anaddb_dtset%iavfrq/=0) write(nunit,'(3x,a9,3i10)')'    iavfrq',anaddb_dtset%iavfrq
   write(nunit,'(a)')' Description of grid 2 :'
   write(nunit,'(3x,a9,3i10)')'   ng2qpt',anaddb_dtset%ng2qpt(1:3)
   write(nunit,'(3x,a9,3i10)')'   ngrids',anaddb_dtset%ngrids
   write(nunit,'(3x,a9,7x,3es16.8)')'   q2shft',anaddb_dtset%q2shft(1:3)
 end if

!Non-linear response information
 if (anaddb_dtset%nlflag /= 0) then
   write(nunit,'(a)')' Non-linear response information :'
   write(nunit,'(3x,a9,i10)') '   alphon',anaddb_dtset%alphon
   write(nunit,'(3x,a9,3i10)')'   prtmbm',anaddb_dtset%prtmbm
   write(nunit,'(3x,a9,3i10)')'  ramansr',anaddb_dtset%ramansr
 end if

!Structural relaxation at fixed polarization
 if (anaddb_dtset%polflag /= 0) then
   write(nunit,'(a)')' Relaxation at fixed polarization :'
   if (anaddb_dtset%relaxat == 1) then
     write(nunit,'(3x,a9,i10)') '  relaxat',anaddb_dtset%relaxat
   end if
   if (anaddb_dtset%relaxstr == 1) then
     write(nunit,'(a12,i10)') ' relaxstr',anaddb_dtset%relaxstr
   end if
 end if

!Elphon information
 if (anaddb_dtset%elphflag /= 0) then
   write(nunit,'(a)')' Elphon calculation will be carried out'
   write(nunit,'(a12,E16.6)') 'elphsmear', anaddb_dtset%elphsmear
   write(nunit,'(a12,E16.6)') 'a2fsmear', anaddb_dtset%a2fsmear
   write(nunit,'(a12,E16.6)') 'mustar', anaddb_dtset%mustar
   write(nunit,'(a12,i10)') 'nqpath', anaddb_dtset%nqpath
   write(nunit,'(a12)') 'qpath'
   do iqpt=1,anaddb_dtset%nqpath
     write(nunit,'(12x,3(E16.6,1x))') anaddb_dtset%qpath(:,iqpt)
   end do
   write(nunit,'(a12,i10)') 'telphint', anaddb_dtset%telphint
   if (anaddb_dtset%telphint == 0) then
     write(nunit,'(a)') ' Tetrahedron integration for elphon'
   else if (anaddb_dtset%telphint == 1) then
     write(nunit,'(a)') ' Smeared weight integration for elphon'
   else if (anaddb_dtset%telphint == 2) then
     write(nunit,'(a)') ' Band filtered integration for elphon'
   end if
   if (anaddb_dtset%elph_fermie /= 0) then
     write(nunit,'(a12,E16.6)')  'elph_fermie', anaddb_dtset%elph_fermie
   end if
   if (anaddb_dtset%ep_extrael /= 0) then
     write(nunit,'(a12,E16.6)')  'Elphon: extra electrons per unit cell = ', anaddb_dtset%ep_extrael
   end if

!  anaddb_dtset%telphint == 0 .or. anaddb_dtset%ep_alter_int_gam == 1 .or. &
!  &      anaddb_dtset%prtnest==1 .or. anaddb_dtset%prtnest==2) then
   if (sum(abs(anaddb_dtset%kptrlatt)) > 0) then
     write(nunit,'(a12,3(3(i3,1x),2x))' ) 'kptrlatt',&
&     reshape( anaddb_dtset%kptrlatt(:,:), (/9/) )
   end if

   if (sum(abs(anaddb_dtset%kptrlatt_fine)) > 0) then
     write(nunit,'(a12,3(3(i3,1x),2x))' ) 'kptrlatt_fine ',&
&     reshape( anaddb_dtset%kptrlatt_fine(:,:), (/9/) )
   end if

   if (anaddb_dtset%ep_keepbands == 1) then
     write(nunit, '(a)') ' Will keep band dependency in gkk in memory.'
     write(nunit, '(a)') ' WARNING: the memory requirements will be multiplied by nbands**2 !!!'
   end if

   if (anaddb_dtset%ep_scalprod == 1) then
     write(nunit, '(a)') ' scalar product will be performed when assembling the gamma matrices.'
     write(nunit, '(a)') ' WARNING: with this option you can not distinguish which '
     write(nunit, '(a)') '    linewidth comes from which phonon mode !!!'
   end if

   if (anaddb_dtset%prtfsurf == 1) then
     write(nunit, '(a)') ' Will output fermi surface in XCrysDen format'
   end if

   if (anaddb_dtset%prtnest == 1) then
     write(nunit, '(a)') ' Will output nesting factor'
   end if

   if (anaddb_dtset%ifltransport == 1) then
     write(nunit, '(a)') ' Will perform transport calculation in elphon to get'
     write(nunit, '(a,a)') ' resistivity and thermal conductivity as a function of T',ch10
     write(nunit, '(a,es16.6,a)' ) ' Minimum temperature for transport outputs: ', &
&     anaddb_dtset%tempermin, ' K'
     write(nunit, '(a,es16.6,a)' ) ' Maximum temperature for transport outputs: ', &
&     anaddb_dtset%tempermin+anaddb_dtset%temperinc*anaddb_dtset%ntemper, ' K'
     write(nunit, '(a,i6)' ) ' Number of temperature points for transport outputs: ', anaddb_dtset%ntemper
     write(nunit, '(a)' ) 
   end if

   if (anaddb_dtset%gkqwrite == 1) then
     write(nunit,'(a,a)' ) 'Gkk matrix elements on input grid of ',&
     'qpoints will be written to disk. File gkqfile must be absent.'
   end if
   if (anaddb_dtset%gkk_rptwrite == 1) then
     write(nunit,'(a,a)' ) 'Gkk matrix elements in real space ',&
     'will be written to disk. File gkk_rpt_file must be absent.'
   end if
   if (anaddb_dtset%gkk2write == 1) then
     write(nunit,'(a,a)' ) 'Full grid gkk matrix elements ',&
     'will be written to disk. File gkk2file must be absent.'
   end if
 end if

!List of vector 1  (reduced coordinates)
 if(anaddb_dtset%nph1l/=0)then
   write(nunit,'(a)')' First list of wavevector (reduced coord.) :'
   write(nunit,'(3x,a9,3i10)')'    nph1l',anaddb_dtset%nph1l
   write(nunit,'(3x,a9)')'    qph1l'
   do iph1=1,anaddb_dtset%nph1l
     write(nunit,'(19x,3es16.8,2x,es11.3)') &
&     (anaddb_dtset%qph1l(ii,iph1),ii=1,3),anaddb_dtset%qnrml1(iph1)
   end do
 end if

!DEBUG
!write(std_out,*)' outvars9 : 1 stop '
!stop
!ENDDEBUG

!List of vector 2  (cartesian coordinates)
 if(anaddb_dtset%nph2l/=0)then
   write(nunit,'(a)')' Second list of wavevector (cart. coord.) :'
   write(nunit,'(3x,a9,3i10)')'    nph2l',anaddb_dtset%nph2l
   write(nunit,'(3x,a9)')'    qph2l'
   do iph2=1,anaddb_dtset%nph2l
     write(nunit,'(19x,3es16.8,2x,es11.3)') &
&     (anaddb_dtset%qph2l(ii,iph2),ii=1,3),anaddb_dtset%qnrml2(iph2)
   end do
 end if

!DEBUG
!write(std_out,*)' outvars9 : 2 stop '
!stop
!ENDDEBUG


!phonon frozen in supercell
 if (abs(anaddb_dtset%freeze_displ) > tol10) then
   write(nunit,'(a)') 'Phonon displacements will be output, frozen into supercells'
   write(nunit,'(a,E20.10)') ' Chosen amplitude of frozen displacements = ', anaddb_dtset%freeze_displ
 end if

!atom projected bs files
 if (abs(anaddb_dtset%natprj_bs) > 0) then
   write(nunit,'(a)') 'Phonon band structure files, with atomic projections, will be output '
   write(nunit,'(a)') ' Chosen atoms for projection = '
   write(nunit,'(10I6)') anaddb_dtset%iatprj_bs
 end if


 write(nunit,'(a,80a,a)') ch10,('=',ii=1,80),ch10

!DEBUG
!write(std_out,*)' outvars9 : exit '
!stop
!ENDDEBUG

end subroutine outvars9
!!***
