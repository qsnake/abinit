!{\src2tex{textfont=tt}}
!!****f* ABINIT/inpspheads.F90
!! NAME
!! inpspheads
!!
!! FUNCTION
!! Read the pseudopotential header of each psp file, in order to initialize pspheads(1:npsp).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, FrD, AF, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  npsp=number of pseudopotentials
!!
!! OUTPUT
!!  pspheads(npsp)=<type pspheader_type>=all the important information from the
!!   pseudopotential file headers, as well as the psp file names
!!
!! PARENTS
!!      m_ab6_invars_f90
!!
!! CHILDREN
!!      close_xml_t,leave_new,open_xml_file,parse,psxml2ab,upfheader2abi,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine inpspheads(filnam,npsp,pspheads)

 use m_profiling

 use defs_basis
 use defs_datatypes
#if defined HAVE_TRIO_FOX
 use m_xml_pseudo_types
 use m_xml_pseudo
 use fox_sax
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inpspheads'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_45_psp_parser, except_this_one => inpspheads
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npsp
!arrays
 character(len=fnlen), intent(in) :: filnam(npsp)
 type(pspheader_type),intent(out) :: pspheads(npsp)

!Local variables-------------------------------
!In case a xc core correction is to be taken into account,
!the n1xccc value will be given by n1xccc_default. Otherwise it is set to 0.
!scalars
 integer,parameter :: n1xccc_default=2501
 integer :: ios
 integer :: idum,ii,ilmax,ipsp,lang,lmax,mmax,mpsang,n1xccc,nmesh
 integer :: pspcod,pspso,test_paw,usexml, useupf !,,pspxc
 real(dp) :: al,e990,e999,fchrg,qchrg,r1,rchrg,rr ! ,rp,rs
 character(len=3) :: testxc
 character(len=500) :: message
 character(len=70) :: testxml
 character(len=80) :: pspline
!arrays
 integer :: nproj(0:3),nprojso(1:3)
 integer,allocatable :: orb(:)
 real(dp) :: hdum(3)
!no_abirules
#if defined HAVE_TRIO_FOX
 integer :: il,maxn_pots
 type(xml_t)               :: fxml
 type(pseudo_t), pointer   :: psxml
#endif

!*************************************************************************

 test_paw=0

 do ipsp=1,npsp

   pspheads(ipsp)%filpsp=trim(filnam(ipsp))

!  Check if the file is written in XML
   usexml = 0
   open (unit=tmp_unit,file=filnam(ipsp),form='formatted',status='old',iostat=ios)
   if (ios/=0) then 
     call wrtout(std_out,"Error opening pseudopotential file"//TRIM(filnam(ipsp)),"PERS")
     call leave_new('PERS')
   end if
   rewind (unit=tmp_unit)
   read(tmp_unit,*) testxml
   if(testxml(1:5)=='<?xml')then
     usexml = 1
   else
     usexml = 0
   end if
   close (unit=tmp_unit)

!  Check if pseudopotential file is a Q-espresso UPF file
   useupf = 0
   open (unit=tmp_unit,file=filnam(ipsp),form='formatted',status='old',iostat=ios)
   if (ios/=0) then 
     call wrtout(std_out,"Error opening pseudopotential file"//TRIM(filnam(ipsp)),"PERS")
     call leave_new('PERS')
   end if
   rewind (unit=tmp_unit)
   read(tmp_unit,*) testxml ! just a string, no relation to xml.
   if(testxml(1:9)=='<PP_INFO>')then
     useupf = 1
   else
     useupf = 0
   end if
   close (unit=tmp_unit)

!  Read the header of the pseudopotential file
   if (usexml /= 1 .and. useupf /= 1) then
!    Open the psp file
     open (unit=tmp_unit,file=filnam(ipsp),form='formatted',status='old',iostat=ios)
     if (ios/=0) then 
       call wrtout(std_out,"Error opening pseudopotential file"//TRIM(filnam(ipsp)),"PERS")
       call leave_new('PERS')
     end if
     rewind (unit=tmp_unit)

!    Read the three first lines
     read (tmp_unit, '(a)' )pspheads(ipsp)%title
     read (tmp_unit,*)pspheads(ipsp)%znuclpsp,pspheads(ipsp)%zionpsp,pspheads(ipsp)%pspdat
     read (tmp_unit,*)pspheads(ipsp)%pspcod,pspheads(ipsp)%pspxc,&
&     pspheads(ipsp)%lmax,idum,mmax
     pspcod=pspheads(ipsp)%pspcod
     lmax=pspheads(ipsp)%lmax
     write(message,'(a,f5.1,a,i4,a,i4)')'  read the values zionpsp=',pspheads(ipsp)%zionpsp,' , pspcod=',pspcod,' , lmax=',lmax
     call wrtout(std_out,message,'PERS')

     nproj(0:3)=0 ; nprojso(1:3)=0

     pspheads(ipsp)%xccc=0
     pspheads(ipsp)%pspso=0

   else if( usexml == 1) then

#if defined HAVE_TRIO_FOX
     if(usexml==1)then
       write(message,'(a,a)')  &
&       '- inpspheads : Reading pseudopotential header in XML form from ', trim(filnam(ipsp))
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')

       call open_xml_file(fxml,filnam(ipsp),ios)
       if (ios/=0) then 
         call wrtout(std_out,"Error opening pseudopotential file"//TRIM(filnam(ipsp)),"PERS")
         call leave_new('PERS')
       end if

       call parse(fxml,pcdata_chunk,startElement_handler=begin_element,endElement_handler=end_element)

!      problem: this depends on pseudo which is defined in m_xml_pseudo, used from 65_psp
!      in order to include this call in the header verification, need to allocate a dummy pseudo or something
       psxml => pseudo

       call psxml2ab( psxml,                       &
&       pspheads(ipsp)%znuclpsp,     &
&       pspheads(ipsp)%zionpsp,      &
&       pspheads(ipsp)%pspcod,       &
&       pspheads(ipsp)%pspxc,        &
&       pspheads(ipsp)%lmax, 0 )

       pspcod = pspheads(ipsp)%pspcod
       lmax   = pspheads(ipsp)%lmax

       nproj(0:3) = 0 ; nprojso(1:3) = 0
       if( psxml%header%core_corrections .eq. "yes") then
         pspheads(ipsp)%xccc  = n1xccc_default
       else
         pspheads(ipsp)%xccc  = 0
       end if
       pspheads(ipsp)%pspso = 0

       call close_xml_t(fxml)

!      Deallocate the different parts of psxml
       maxn_pots=size(psxml%pswf)

       do il = 1, maxn_pots
         if(associated( psxml%pswf(il)%V%data))   ABI_DEALLOCATE (psxml%pswf(il)%V%data)
         if(associated( psxml%pot(il)%V%data))    ABI_DEALLOCATE (psxml%pot(il)%V%data)
       end do
       if(associated( psxml%core_charge%data))    ABI_DEALLOCATE (psxml%core_charge%data)
       if(associated (psxml%valence_charge%data)) ABI_DEALLOCATE (psxml%valence_charge%data)

     else
#endif
       write(std_out,'(4a)' ) ch10,&
&       ' inpspheads : In order to use XML pseudopotentials, you need to compile ABINIT',ch10,&
&       '  with the -DFOX preprocessing option, and also to compile the FoX library. Stop.'
       stop
#if defined HAVE_TRIO_FOX
     end if
#endif

   else if (useupf == 1) then
     pspheads(ipsp)%pspcod = 11

     pspheads(ipsp)%xccc  = n1xccc_default ! will be set to 0 if no nlcc
!    call upfoctheader2abi (filnam(ipsp),  &
     call upfheader2abi (filnam(ipsp),  &
&     pspheads(ipsp)%znuclpsp, &
&     pspheads(ipsp)%zionpsp,  &
&     pspheads(ipsp)%pspxc,    &
&     pspheads(ipsp)%lmax,     &
&     pspheads(ipsp)%xccc,     &
&     nproj, nprojso)

     pspcod = pspheads(ipsp)%pspcod
     lmax   = pspheads(ipsp)%lmax
!    FIXME : generalize for SO pseudos
     pspheads(ipsp)%pspso = 0
   end if

!  DEBUG
!  write(std_out,*) pspheads(ipsp)%znuclpsp
!  write(std_out,*) pspheads(ipsp)%zionpsp
!  write(std_out,*) pspheads(ipsp)%pspcod
!  write(std_out,*) pspheads(ipsp)%pspxc
!  write(std_out,*) pspheads(ipsp)%lmax
!  stop
!  ENDDEBUG

!  Initialize nproj, nprojso, pspso, as well as xccc, for each type of psp
   pspheads(ipsp)%GTHradii = zero

   if(pspcod==1 .or. pspcod==4)then

!    Teter format
     do ilmax=0,lmax
       read (tmp_unit,*) lang,e990,e999,nproj(ilmax)
       read (tmp_unit,*)
     end do
     read (tmp_unit,*) rchrg,fchrg,qchrg
     if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default

   else if(pspcod==2)then

!    GTH pseudopotentials
     read (tmp_unit,*) pspheads(ipsp)%GTHradii(0) !rloc
     read (tmp_unit,*) pspheads(ipsp)%GTHradii(1),hdum(1),hdum(2)
     if(abs(hdum(1))>1.d-9) nproj(0)=1
     if(abs(hdum(2))>1.d-9) nproj(0)=2
     read (tmp_unit,*) pspheads(ipsp)%GTHradii(2),hdum(3)
     if(abs(hdum(3))>1.d-9) nproj(1)=1

   else if(pspcod==3)then

!    HGH pseudopotentials
     read (tmp_unit,*) pspheads(ipsp)%GTHradii(0) !rloc
     do ilmax=0,lmax
       read (tmp_unit,*) pspheads(ipsp)%GTHradii(ilmax + 1),hdum(1),hdum(2),hdum(3)
       if (abs(hdum(1))>1.d-9)nproj(ilmax)=1
       if (abs(hdum(2))>1.d-9)nproj(ilmax)=2
       if (abs(hdum(3))>1.d-9)nproj(ilmax)=3
       if (ilmax>0.and.ilmax<3) then
         read (tmp_unit,*) hdum(1),hdum(2),hdum(3)
         if (abs(hdum(1))>1.d-9)nprojso(ilmax)=1
         if (abs(hdum(2))>1.d-9)nprojso(ilmax)=2
         if (abs(hdum(3))>1.d-9)nprojso(ilmax)=3
         if(nprojso(ilmax)>0)pspheads(ipsp)%pspso=2
       end if
       if (ilmax==3) then
         read (tmp_unit,*) hdum(1)
         if (abs(hdum(1))>1.d-9)nprojso(3)=1
         if(nprojso(3)>0)pspheads(ipsp)%pspso=2
       end if
     end do

   else if(pspcod==5)then

!    PHONEY pseudopotentials
!    read parameter for Hamman grid
     pspso=1
     read (tmp_unit,fmt=*,err=10,end=10) r1,al,pspso
     10 continue
     do ilmax=0,lmax
       read (tmp_unit,*) lang,e990,e999,nproj(ilmax)
       read (tmp_unit,*)
       if (ilmax>0.and.pspso/=1) then
         read (tmp_unit,*) lang,e990,e999,nprojso(ilmax)
         read (tmp_unit,*)
         pspheads(ipsp)%pspso=pspso
!        Meaning of pspso internally to ABINIT has been changed in v5.4
!        So : file must contain pspso 1 , but ABINIT will have pspso 0 .
         if(pspso==1)pspheads(ipsp)%pspso=0
       end if
     end do
     read (tmp_unit,*) rchrg,fchrg,qchrg
     if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default

   else if(pspcod==6)then

!    FHI pseudopotentials
     read (tmp_unit, '(a3)') testxc
!    Note : prior to version 2.2, this 4th line started with  4--  ,
!    and no core-correction was available.
     if(testxc/='4--')then
       backspace(tmp_unit)
       read (tmp_unit,*) rchrg,fchrg,qchrg
     else
       fchrg=0.0_dp
     end if
     if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default
!    XG020728 : Should take lloc into account ??
     do ilmax=0,lmax
       nproj(ilmax)=1
     end do

   else if(pspcod==7)then

!    PAW pseudopotentials
     test_paw=1;pspheads(ipsp)%pawheader%pawver=1
     read (tmp_unit,'(a80)') pspline;pspline=adjustl(pspline)
     if (pspline(1:3)=="paw".or.pspline(1:3)=="PAW") &
&     read(unit=pspline(4:80),fmt=*) pspheads(ipsp)%pawheader%pawver
     if (pspheads(ipsp)%pawheader%pawver==1) then   ! Compatibility with Abinit v4.2.x
       read (unit=pspline,fmt=*) pspheads(ipsp)%pawheader%basis_size,&
&       pspheads(ipsp)%pawheader%lmn_size
       ABI_ALLOCATE(orb,(pspheads(ipsp)%pawheader%basis_size))
       orb(:)=0
       read (tmp_unit,*) (orb(ii), ii=1,pspheads(ipsp)%pawheader%basis_size)
       read (tmp_unit,*);read (tmp_unit,*) pspheads(ipsp)%pawheader%rpaw
       pspheads(ipsp)%pawheader%rshp=pspheads(ipsp)%pawheader%rpaw
       read (tmp_unit,*) pspheads(ipsp)%pawheader%mesh_size
       read (tmp_unit,*) pspheads(ipsp)%pawheader%shape_type
       if (pspheads(ipsp)%pawheader%shape_type==3) pspheads(ipsp)%pawheader%shape_type=-1
     else
       read (tmp_unit,*) pspheads(ipsp)%pawheader%basis_size,&
&       pspheads(ipsp)%pawheader%lmn_size
       ABI_ALLOCATE(orb,(pspheads(ipsp)%pawheader%basis_size))
       orb(:)=0
       read (tmp_unit,*) (orb(ii), ii=1,pspheads(ipsp)%pawheader%basis_size)
       pspheads(ipsp)%pawheader%mesh_size=mmax
       read (tmp_unit,*) nmesh
       do ii=1,nmesh
         read(tmp_unit,*)
       end do
       read (tmp_unit,*) pspheads(ipsp)%pawheader%rpaw
       pspheads(ipsp)%pawheader%rshp=pspheads(ipsp)%pawheader%rpaw
       read (tmp_unit,'(a80)') pspline;pspline=adjustl(pspline);write(std_out,*) pspline
       read(unit=pspline,fmt=*) pspheads(ipsp)%pawheader%shape_type
       if (pspheads(ipsp)%pawheader%pawver==2.and.&
&       pspheads(ipsp)%pawheader%shape_type==3) pspheads(ipsp)%pawheader%shape_type=-1
       if (pspheads(ipsp)%pawheader%pawver>=3.and.pspheads(ipsp)%pawheader%shape_type==-1) then
         rr=zero;read(unit=pspline,fmt=*,err=20,end=20) ii,rr
         20   continue
         if (rr>=tol8) pspheads(ipsp)%pawheader%rshp=rr
       end if
     end if
     do ilmax=0,lmax
       do ii=1,pspheads(ipsp)%pawheader%basis_size
         if(orb(ii)==ilmax) nproj(ilmax)=nproj(ilmax)+1
       end do
     end do
     pspheads(ipsp)%pawheader%l_size=2*maxval(orb)+1
     pspheads(ipsp)%xccc=1  ! We suppose apriori that cc is used (but n1xccc is not used in PAW)
     ABI_DEALLOCATE(orb)

   else if(pspcod==8)then

!    DRH pseudopotentials
     read (tmp_unit,*) rchrg,fchrg,qchrg
     if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default
     read(tmp_unit,*) nproj(0:lmax)
     pspso=0
     pspheads(ipsp)%pspso=pspso

   else if(pspcod==9)then

     fchrg = 0.0_dp
     do ilmax = 0, lmax
       nproj(ilmax) = 1
     end do

   else if(pspcod==10)then

!    HGH pseudopotentials, full h/k matrices
     read (tmp_unit,*) pspheads(ipsp)%GTHradii(0) !rloc
     read (tmp_unit,*) idum; if(idum-1/=lmax) stop "in inpspheads: nnonloc-1/ = lmax"
     do ilmax=0,lmax
       read (tmp_unit,*) pspheads(ipsp)%GTHradii(ilmax + 1),nproj(ilmax),(hdum(idum),idum=1,nproj(ilmax))
       do idum=2,nproj(ilmax) !skip the rest of h_ij
         read (tmp_unit,*)
       end do
       if (ilmax==0) cycle
       nprojso(ilmax)=nproj(ilmax)
       if(nprojso(ilmax)>0)then
         pspheads(ipsp)%pspso=2
         do idum=1,nprojso(ilmax) !skip the rest of k_ij
           read (tmp_unit,*)
         end do
       end if
     end do

   else if (pspcod == 11) then
!    already done above
   else

     write(message, '(a,a,a,a,i4,a,a,a,a)' ) ch10,&
&     ' inpspheads : ERROR -',ch10,&
&     '  The pseudopotential code (pspcod) read from file is ',pspcod,ch10,&
&     '  This value is not allowed (should be between 1 and 10). ',ch10,&
&     '  Action : use a correct pseudopotential file.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')

   end if ! pspcod=...

!  Store in pspheads
   pspheads(ipsp)%nproj(0:3)=nproj(0:3)
   pspheads(ipsp)%nprojso(1:3)=nprojso(1:3)

   close(tmp_unit)

 end do ! ipsp=1,npsp

!Note that mpsang is the max of 1+lmax, with minimal value 1 (even for local psps, at present)
!mpsang=max(maxval(pspheads(1:npsp)%lmax)+1,1) ! Likely troubles with HP compiler
!n1xccc=maxval(pspheads(1:npsp)%xccc)
 mpsang=1
 n1xccc=pspheads(1)%xccc
 do ii=1,npsp
   mpsang=max(pspheads(ii)%lmax+1,mpsang)
   n1xccc=max(pspheads(ii)%xccc,n1xccc)
 end do

 write(message,'(2a,i4,a,i4,a)')ch10,&
& ' inpspheads : deduce mpsang  =',mpsang,', n1xccc  =',n1xccc,'.'
 call wrtout(std_out,message,'PERS')

!Test: if one psp is PAW, all must be
 if (test_paw==1) then
   do ipsp=1,npsp
     if (pspheads(ipsp)%pspcod/=7) then
       write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&       ' inpspheads : ERROR -',ch10,&
&       '  One pseudopotential is PAW (pspcod=7) !',ch10,&
&       '  All pseudopotentials must be PAW (this is not the case here) !',ch10,&
&       '  Action : use only PAW pseudopotential files.'
       call wrtout(std_out,message,'PERS')
       call leave_new('PERS')
     end if
   end do
 end if

end subroutine inpspheads
!!***
