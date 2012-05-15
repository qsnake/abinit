!{\src2tex{textfont=tt}}
!!****f* ABINIT/prt_cml2
!! NAME
!! prt_cml2
!!
!!
!! FUNCTION
!! Produce a CML (Chemical Markup Language) file
!! with crystalline cell description, symmetries,
!! and atomic coordinates in reduced or cartesian coordinates.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  coord_option=option for writing CML file (see below)
!!  fnameabo_app_cml_xml=name of the file _CML.xml
!!  natom=number of atoms in unit cell
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of types of atoms in unit cell.
!!  rprimd(3,3)=real space dimensional primitive translations (bohr)
!!  spgroup=symmetry space group number
!!  symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!  tnons(3,nsym)=reduced nonsymmorphic translations
!!  typat(natom)=type integer for each atom in cell
!!  xred(3,natom)=reduced coordinates of atoms
!!  znucl(ntypat)=real(dp), atomic number of atom type
!!
!! OUTPUT
!! data written in file fnameabo_app_cml_xml
!!   if coord_option=1, atomic positions are written
!!      in reduced format (xFract,yFract,zFract)
!!   if coord_option=2, atomic positions are written
!!      in cartesian format (x3,y3,z3)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      atmdata,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prt_cml2(coord_option,fnameabo_app_cml_xml,natom,nsym,ntypat,rprimd,spgroup,symrel,tnons,typat,xred,znucl)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prt_cml2'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: coord_option,natom,nsym,ntypat,spgroup
 character(len=fnlen),intent(in) :: fnameabo_app_cml_xml
!arrays
 integer,intent(in) :: symrel(3,3,nsym),typat(natom)
 real(dp),intent(in) :: rprimd(3,3),tnons(3,nsym),xred(3,natom),znucl(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iatom,ii,isym,mu,nu
 real(dp) :: amu,rcov
 character(len=2) :: string2,symbol
 character(len=20) :: string20,xstring20,ystring20,zstring20
 character(len=3) :: string3
 character(len=500) :: message
!arrays
 real(dp) :: angle(3),rmet(3,3)
 real(dp), allocatable :: xangst(:,:),xcart(:,:),xred_(:,:)
 character(len=3) :: string2array(3,3)

! *************************************************************************

!Initialize the file fnameabo_app_cml_xml
 write(message, '(a,a)' ) ' prt_cml2 : about to open file ',fnameabo_app_cml_xml
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 open (unit=tmp_unit,file=fnameabo_app_cml_xml,status='unknown',form='formatted')
 rewind(tmp_unit)

!Take care of the header, and initialize the <molecule> element
!FIXME: this xml name space is no longer available online and CML has been
!upgraded to version 2.4 then 3. They use xslt schemas now.
 write(message, '(3a)' ) &
& '<?xml version="1.0" encoding="iso-8859-1"?>',ch10,&
& '<molecule id="crystal1" xmlns="http://www.xml-cml.org/schema/cml2/core">'
 call wrtout(tmp_unit,message,'COLL')

!Compute real space metrics
 do ii=1,3
   rmet(ii,:)=rprimd(1,ii)*rprimd(1,:)+&
&   rprimd(2,ii)*rprimd(2,:)+&
&   rprimd(3,ii)*rprimd(3,:)
 end do

!Compute angles in degree
 angle(1)=acos(rmet(2,3)/sqrt(rmet(2,2)*rmet(3,3)))/two_pi*360.0d0
 angle(2)=acos(rmet(1,3)/sqrt(rmet(1,1)*rmet(3,3)))/two_pi*360.0d0
 angle(3)=acos(rmet(1,2)/sqrt(rmet(1,1)*rmet(2,2)))/two_pi*360.0d0

!Write the <crystal> element
 write(message, '(a)' ) ' <crystal>'
 call wrtout(tmp_unit,message,'COLL')
 write(string20, '(f20.12)')sqrt(rmet(1,1))*Bohr_Ang
 write(message, '(a,a,a)' )&
& '  <scalar title="a" units="angstrom">',&
& trim(adjustl(string20)),'</scalar>'
 call wrtout(tmp_unit,message,'COLL')
 write(string20, '(f20.12)')sqrt(rmet(2,2))*Bohr_Ang
 write(message, '(a,a,a)' )&
& '  <scalar title="b" units="angstrom">',&
& trim(adjustl(string20)),'</scalar>'
 call wrtout(tmp_unit,message,'COLL')
 write(string20, '(f20.12)')sqrt(rmet(3,3))*Bohr_Ang
 write(message, '(a,a,a)' )&
& '  <scalar title="c" units="angstrom">',&
& trim(adjustl(string20)),'</scalar>'
 call wrtout(tmp_unit,message,'COLL')
 write(string20, '(f20.12)')angle(1)
 write(message, '(a,a,a)' )&
& '  <scalar title="alpha" units="degrees">',&
& trim(adjustl(string20)),'</scalar>'
 call wrtout(tmp_unit,message,'COLL')
 write(string20, '(f20.12)')angle(2)
 write(message, '(a,a,a)' )&
& '  <scalar title="beta"  units="degrees">',&
& trim(adjustl(string20)),'</scalar>'
 call wrtout(tmp_unit,message,'COLL')
 write(string20, '(f20.12)')angle(3)
 write(message, '(a,a,a)' )&
& '  <scalar title="gamma" units="degrees">',&
& trim(adjustl(string20)),'</scalar>'
 call wrtout(tmp_unit,message,'COLL')
 write(message, '(a)' )' </crystal>'
 call wrtout(tmp_unit,message,'COLL')

 write(message, '(a)' )&
& ' <!-- "ITC" refers to the space group number in the International Tables for Crystallography -->'
 call wrtout(tmp_unit,message,'COLL')

!Write the <symmetry> element
 write(string3, '(i3)')spgroup
 write(message, '(4a)' )' <symmetry id="s1" ',&
& 'pointGroup="ITC:',trim(adjustl(string3)),'">'
 call wrtout(tmp_unit,message,'COLL')
 do isym=1,nsym
   do mu=1,3
     do nu=1,3
       write(string2,'(i2)')symrel(mu,nu,isym)
       string2array(mu,nu)=string2
     end do
   end do
   write(string3, '(i3)')isym
   write(message, '(3a)')&
&   '  <matrix id="symOp',trim(adjustl(string3)),'" rows="3" columns="4">'
   call wrtout(tmp_unit,message,'COLL')
   do mu=1,3
     write(string20,'(f20.12)')tnons(mu,isym)
     write(message, '(8a)' ) '    ',&
&     trim(adjustl(string2array(mu,1))),' ',&
&     trim(adjustl(string2array(mu,2))),' ',&
&     trim(adjustl(string2array(mu,3))),' ',&
&     trim(adjustl(string20))
     call wrtout(tmp_unit,message,'COLL')
   end do
   write(message, '(a)')'  </matrix>'
   call wrtout(tmp_unit,message,'COLL')
 end do
 write(message, '(a)' )' </symmetry>'
 call wrtout(tmp_unit,message,'COLL')

!Initialize the <atomArray> element
 write(message, '(a)' )' <atomArray>'
 call wrtout(tmp_unit,message,'COLL')

!Loop over all atoms
 if(coord_option==1)then
   do iatom=1,natom
     call atmdata(amu,rcov,symbol,znucl(typat(iatom)))
     write(xstring20,'(f20.12)')xred(1,iatom)
     write(ystring20,'(f20.12)')xred(2,iatom)
     write(zstring20,'(f20.12)')xred(3,iatom)
     write(string3, '(i3)')iatom
     write(message, '(14a)')'  <atom id="',trim(adjustl(string3)),&
&     '" elementType="',trim(adjustl(symbol)),'" ',&
&     'xFract="',trim(adjustl(xstring20)),'" ',&
&     'yFract="',trim(adjustl(ystring20)),'" ',&
&     'zFract="',trim(adjustl(zstring20)),'"/>'
     call wrtout(tmp_unit,message,'COLL')
   end do
 else if(coord_option==2)then
   ABI_ALLOCATE(xcart,(3,natom))
   ABI_ALLOCATE(xangst,(3,natom))
   ABI_ALLOCATE(xred_,(3,natom))
   xred_(:,:)=xred(:,:)
   call xredxcart(natom,1,rprimd,xcart,xred_)
   xangst(:,:)=xcart(:,:)*Bohr_Ang
   do iatom=1,natom
     call atmdata(amu,rcov,symbol,znucl(typat(iatom)))
     write(xstring20,'(f20.12)')xangst(1,iatom)
     write(ystring20,'(f20.12)')xangst(2,iatom)
     write(zstring20,'(f20.12)')xangst(3,iatom)
     write(string3, '(i3)')iatom
     write(message, '(14a)')'  <atom id="',trim(adjustl(string3)),&
&     '" elementType="',trim(adjustl(symbol)),'" ',&
&     'x3="',trim(adjustl(xstring20)),'" ',&
&     'y3="',trim(adjustl(ystring20)),'" ',&
&     'z3="',trim(adjustl(zstring20)),'"/>'
     call wrtout(tmp_unit,message,'COLL')
   end do
   ABI_DEALLOCATE(xcart)
   ABI_DEALLOCATE(xangst)
 end if

!Finalize the <atomArray> element
 write(message, '(a)' )' </atomArray>'
 call wrtout(tmp_unit,message,'COLL')

!Finalize the <molecule> element
 write(message, '(a)' ) '</molecule>'
 call wrtout(tmp_unit,message,'COLL')

 close(tmp_unit)

end subroutine prt_cml2
!!***
