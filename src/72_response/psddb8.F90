!{\src2tex{textfont=tt}}
!!****f* ABINIT/psddb8
!!
!! NAME
!! psddb8
!!
!! FUNCTION
!! Take care of the i/o of pseudopotentials for the
!! Derivative DataBase, and also the number of data blocks.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  choice=(1 => read), (2=> write)
!!  dimekb=dimension of ekb (contains Kleimann-Bylander energies)
!!         used only for norm-conserving pseudopotentials
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  nunit=unit number for the Derivative DataBase.
!!  ntypat=number of atom types
!!  pspso(ntypat)=For each type of psp, 1 if no spin-orbit component is taken
!!     into account, 2 if a spin-orbit component is used
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  vrsddb=Derivative Database version, for check of compatibility.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                                 or i=lmn  (if useylm=1)
!!  ekb(dimekb,ntypat)= (norm-conserving psps only) (Real) Kleinman-Bylander energies (hartree)
!!                      Presently the only characteristics of the psp
!!  fullinit=0 if the ekb are not available, at input as well as at output
!!  pawtab(ntypat*usepaw)= (PAW only) PAW datasets characteristics
!!                  Presently only pawtab%basis_size,pawtab%lmn_size,pawtab%shape_type
!!                  pawtab%rpaw,pawtab%rshp,pawtab%dij0  are used
!!  nblok=number of blocks
!!
!! NOTES
!! Only executed by one processor
!!
!! PARENTS
!!      gstate,loper3,mblktyp1,mblktyp5,nonlinear,rdddb9,respfn,thmeig
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,&
&          nblok,ntypat,nunit,pawtab,pspso,usepaw,useylm,vrsddb)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psddb8'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,dimekb,lmnmax,ntypat,nunit,usepaw,useylm
 integer,intent(in) :: vrsddb
 integer,intent(inout) :: fullinit,nblok
!arrays
 integer,intent(in) :: pspso(ntypat)
 integer,intent(inout) :: indlmn(6,lmnmax,ntypat)
 real(dp),intent(inout) :: ekb(dimekb,ntypat)
 type(pawtab_type),intent(inout) :: pawtab(ntypat*usepaw)

!Local variables -------------------------
!Set the version number
!scalars
 integer,parameter :: vrsio8=100401,vrsio8_old=010929,vrsio8_old_old=990527
 integer :: basis_size0,dimekb0,iekb,ii,ij,il,ilm,ilmn,iln,iln0,im,ios,iproj,iproj0,itypat,itypat0
 integer :: jekb,jlmn,jln,lmnmax0,lmn_size0,lmn2_size0,lpsang,nekb,nproj,npsang,pspso0,shape_type0
 integer :: usepaw0,vrspsp8
 real(dp) :: rpaw0,rshape0
 character(len=12) :: string
 character(len=500) :: message
!arrays
 integer,allocatable :: i1(:),i2(:),nprj(:),orbitals(:)
 real(dp),allocatable :: dij0(:),ekb0(:,:)

! *********************************************************************


!Check psddb8 version number (vrsio8) against DDB version number
!(vrsddb)
 if (vrsio8/=vrsddb) then
   write(message, '(a,a,a,a,i10,a,a,i10,a)' ) ch10,&
&   ' psddb8 : BUG -',ch10,&
&   '  the psddb8 DDB version number=',vrsio8,ch10,&
&   '    is not equal to the calling code DDB version number=',vrsddb,'.'
   call wrtout(std_out,message,'COLL')
!  call leave_new('COLL')
 end if

!Check the value of choice
 if (choice<=0.or.choice>=3) then
   write(message, '(a,a,a,a,a,a,i10,a)' ) ch10,&
&   ' psddb8: BUG -',ch10,&
&   '  The permitted values for choice are 1 or 2.',ch10,&
&   '  The calling routine asks ',choice,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!==================================================================================
!First option: read psp characteristic from file =================================
!==================================================================================
 if (choice==1) then

   read(nunit,*)
   read(nunit, '(a12)' )string
   fullinit=1 ; if (lmnmax>0) indlmn(:,:,:)=0

!  --------------------------------------------
!  -----  NEW FORMAT (NCPP+PAW) ---------------
!  --------------------------------------------
   if (string=='  Descriptio')then


!    ==============================
!    ======== Common data =========
!    ==============================
     read(nunit, '(32x,i6)' )vrspsp8
     if (vrspsp8==vrsio8_old.or.vrspsp8==vrsio8_old_old) then
       usepaw0=0
       read(nunit, '(10x,i3,14x,i3,11x,i3)', iostat=ios)dimekb0,lmnmax0,usepaw0
       if(ios/=0)then
         backspace(nunit)
         read (nunit, '(10x,i3,14x,i3)' )dimekb0,lmnmax0
         usepaw0=0
       end if
     else if (vrspsp8==vrsio8) then
       read(nunit, '(10x,i3)') usepaw0
       if (usepaw/=usepaw0) then
         write(message, '(4a,i1,a,i1,a)' ) ch10,&
&         ' psddb8: ERROR -',ch10,&
&         '  usepaw is announced to be ',usepaw,' but read usepaw is ',usepaw0,' !'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
       if (usepaw==0) then
         read (nunit, '(10x,i3,14x,i3)' )dimekb0,lmnmax0
       end if
     end if

!    ==============================
!    === Norm-conserving psps =====
!    ==============================
     if (usepaw==0) then
       ekb(:,:)=zero
       ABI_ALLOCATE(ekb0,(dimekb,dimekb))
       do itypat=1,ntypat
         read(nunit, '(13x,i4,9x,i3,8x,i4)' )itypat0,pspso0,nekb
!        Check the compatibility with the main code dimensioning
         if(nekb>dimekb)then
           write(message, '(a,a,a,a,i8,a,a,a,i3,a)' ) ch10,&
&           ' psddb8 : BUG -',ch10,&
&           '  ',nekb,' components of ekb are announced',ch10,&
&           '  but dimekb=',dimekb,'.'
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if
         read(nunit,*)
         ilmn=0;iproj0=0
         do iekb=1,nekb
           read(nunit, '(3i6,3x,8d15.7)' ) iln,lpsang,iproj,&
&           (ekb0(ii,iekb),ii=1,min(nekb,4))
           if(nekb>4)then
             do jekb=5,nekb,4
               read(nunit, '(21x,8d15.7)' )&
&               (ekb0(ii,iekb),ii=jekb,min(nekb,jekb+3))
             end do
           end if
           if (lpsang==0.and.iproj>iproj0) iproj0=iproj
           if (useylm==1) then
             do im=-lpsang,lpsang
               ilmn=ilmn+1
               indlmn(1,ilmn,itypat)=lpsang
               indlmn(2,ilmn,itypat)=im
               indlmn(3,ilmn,itypat)=iproj
               indlmn(4,ilmn,itypat)=lpsang**2+lpsang+1+im
               indlmn(5,ilmn,itypat)=iln
               indlmn(6,ilmn,itypat)=1
               if (pspso0/=1.and.iln>(nekb-iproj0)/2) indlmn(6,ilmn,itypat)=2
             end do
           else
             ilmn=ilmn+1
             indlmn(1,ilmn,itypat)=lpsang
             indlmn(2,ilmn,itypat)=lpsang
             indlmn(3,ilmn,itypat)=iproj
             indlmn(4,ilmn,itypat)=lpsang**2+lpsang+1
             indlmn(5,ilmn,itypat)=iln
             indlmn(6,ilmn,itypat)=1
             if (pspso0/=1.and.iln>(nekb-iproj0)/2) indlmn(6,ilmn,itypat)=2
           end if
!          For the time being, only diagonal ekb are treated in abinit v3
           ekb(iekb,itypat)=ekb0(iekb,iekb)
!          For non-diagonal ekb, one could use:
!          do jekb=iekb to nekb
!          ekb(jekb+iekb*(iekb-1)/2,itypat)=ekb0(jekb,iekb)
!          end do
         end do
       end do
       ABI_DEALLOCATE(ekb0)

!      ==============================
!      ============ PAW =============
!      ==============================
     else
       do itypat=1,ntypat
         read(nunit, '(12x,i4,12x,i3,12x,i5)' )itypat0,basis_size0,lmn_size0
         lmn2_size0=lmn_size0*(lmn_size0+1)/2
         ABI_ALLOCATE(orbitals,(basis_size0))
         read(nunit, '(20x,50i2)' ) orbitals(1:basis_size0)
         read(nunit, '(11x,f6.3,13x,i2,11x,f6.3)' ) rpaw0,shape_type0,rshape0
         read(nunit,'(24x,i3)') nekb
         read(nunit,*)
         ABI_ALLOCATE(dij0,(nekb))
         ABI_ALLOCATE(i1,(nekb))
         ABI_ALLOCATE(i2,(nekb))
         do ii=1,nekb,4
           read(nunit,'(3x,4(1x,i4,1x,i4,1x,d12.5))') (i1(ij),i2(ij),dij0(ij),ij=ii,min(ii+3,nekb))
         end do
         if (lmn_size0>lmnmax) then
           write(message, '(4a,i5,3a,i5,a)' ) ch10,&
&           ' psddb8 : BUG -',ch10,&
&           '  max. value of ',lmnmax,' for lmn_size is announced',ch10,&
&           '  but ',lmn_size0,' is read.'
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if
         if (associated(pawtab(itypat)%dij0)) then
           if (lmn_size0>pawtab(itypat)%lmn_size) then
             write(message, '(4a,i5,3a,i5,a)' ) ch10,&
&             ' psddb8 : BUG -',ch10,&
&             '  lmn_size=,',pawtab(itypat)%lmn_size,' is announced',ch10,&
&             '  but ',lmn_size0,' is read.'
             call wrtout(std_out,message,'COLL')
             call leave_new('COLL')
           end if
         end if
         ABI_ALLOCATE(nprj,(0:maxval(orbitals)))
         ilmn=0;nprj=0
         do iln=1,basis_size0
           il=orbitals(iln)
           nprj(il)=nprj(il)+1
           do ilm=1,2*il+1
             indlmn(1,ilmn+ilm,itypat)=il
             indlmn(2,ilmn+ilm,itypat)=ilm-(il+1)
             indlmn(3,ilmn+ilm,itypat)=nprj(il)
             indlmn(4,ilmn+ilm,itypat)=il*il+ilm
             indlmn(5,ilmn+ilm,itypat)=iln
             indlmn(6,ilmn+ilm,itypat)=1
           end do
           ilmn=ilmn+2*il+1
         end do
         pawtab(itypat)%basis_size=basis_size0
         pawtab(itypat)%lmn_size  =lmn_size0
         pawtab(itypat)%lmn2_size =lmn2_size0
         pawtab(itypat)%shape_type=shape_type0
         pawtab(itypat)%rpaw      =rpaw0
         pawtab(itypat)%rshp      =rshape0
         if (.not.associated(pawtab(itypat)%dij0))  then
           ABI_ALLOCATE(pawtab(itypat)%dij0,(lmn2_size0))
         end if
         pawtab(itypat)%dij0(1:lmn2_size0)=zero
         do ii=1,nekb
           ij=i1(ii)+i2(ii)*(i2(ii)-1)/2
           pawtab(itypat)%dij0(ij)=dij0(ii)
         end do
         ABI_DEALLOCATE(nprj)
         ABI_DEALLOCATE(orbitals)
         ABI_DEALLOCATE(dij0)
         ABI_DEALLOCATE(i1)
         ABI_DEALLOCATE(i2)
       end do

     end if ! NCPP or PAW

!    --------------------------------------------
!    -----  OLD FORMAT (NCPP only) --------------
!    --------------------------------------------
   else if (string==' Description')then
     if (usepaw==1) stop 'BUG: old DDB pspformat not compatible with PAW 1'

     read (nunit, '(10x,i3,10x,i3)' )nproj,npsang
     nekb=nproj*npsang
!    Check the compatibility with the main code dimensioning
     if(nekb>dimekb)then
       write(message, '(a,a,a,a,i8,a,a,a,i3,a)' ) ch10,&
&       ' psddb8 : BUG -',ch10,&
&       '  ',nekb,' components of ekb are announced',ch10,&
&       '  but the maximum is dimekb=',dimekb,'.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if(useylm/=0)then
       write(message, '(a,a,a,a)' ) ch10,&
&       ' psddb8: BUG -',ch10,&
&       '  useylm must be 0 !'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
!    Read the data
     ABI_ALLOCATE(ekb0,(dimekb,dimekb))
     ekb0(:,:)=zero
     do itypat=1,ntypat
       read (nunit, '(13x,i4)' )ij
       do iproj=1,nproj
         read (nunit, '(6x,3d22.14)' )&
&         (ekb0(iproj+nproj*(ii-1),iproj+nproj*(ii-1)),ii=1,min(npsang,3))
         if(npsang>3)read (nunit, '(6x,3d22.14)' )&
&         (ekb0(iproj+nproj*(ii-1),iproj+nproj*(ii-1)),ii=4,npsang)
         do ii=1,npsang
           iekb=iproj+nproj*(ii-1)
           indlmn(1,iekb,itypat)=ii-1
           indlmn(2,iekb,itypat)=ii-1
           indlmn(3,iekb,itypat)=iproj
           indlmn(4,iekb,itypat)=ii**2-ii+1
           indlmn(5,iekb,itypat)=iekb
           indlmn(6,iekb,itypat)=1
!          For the time being, only diagonal ekb are treated in abinit v3
           ekb(iekb,itypat)=ekb0(iekb,iekb)
         end do
       end do
     end do
     ABI_DEALLOCATE(ekb0)

!    --------------------------------------------
!    -----  OTHER CASES -------------------------
!    --------------------------------------------
   else if(string==' No informat')then
     fullinit=0
   else
     write(message, '(a,a,a)' )&
&     ' psddb8 : BUG -',ch10,&
&     '  Error when reading the psp information'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  Now, the number of blocks
   read(nunit,*)
   read(nunit,*)
   read(nunit, '(24x,i4)' )nblok

!  ==================================================================================
!  Second option: read psp characteristic from file ================================
!  ==================================================================================
 else if(choice==2)then

   write(nunit, '(a)' )' '
   if (fullinit==0)then
!    This possibility is used when the DDB is initialized,
!    and the ekb s are not available from the GS input file...
     write(nunit, '(a)' )&
&     ' No information on the potentials yet '
   else

!    ==============================
!    === Norm-conserving psps =====
!    ==============================
     if (usepaw==0) then
       write(nunit, '(a)' )&
&       '  Description of the potentials (KB energies)'
       write(nunit, '(a,i6)' )&
&       '  vrsio8 (for pseudopotentials)=',vrsio8
       write(nunit, '(a,i3)' ) '  usepaw =',usepaw
       write(nunit, '(a,i3,a,i3,a,i3)' )&
&       '  dimekb =',dimekb,'       lmnmax=',lmnmax
       ABI_ALLOCATE(ekb0,(dimekb,dimekb))
       do itypat=1,ntypat
!        Compute nekb
         nekb=0
         do jlmn=1,lmnmax
           jln=indlmn(5,jlmn,itypat)
           if(jln>nekb)then
             nekb=jln
           end if
         end do
         write(nunit, '(a,i4,a,i3,a,i4)' ) &
&         '  Atom type= ',itypat,'   pspso=',pspso(itypat),'   nekb=',nekb
         write(nunit, '(a)' ) '  iln lpsang iproj  ekb(:)'
         iln0=0
         ekb0(:,:)=zero
         do ilmn=1,lmnmax
           iln =indlmn(5,ilmn,itypat)
           if (iln>iln0) then
             iln0=iln
             lpsang=indlmn(1,ilmn,itypat)
             iproj=indlmn(3,ilmn,itypat)
!            For the time being, only diagonal ekb are treated in abinit v3
             ekb0(iln,iln)=ekb(iln,itypat)
!            For non-diagonal ekb, one could use:
!            do ii=iln to nekb
!            ekb0(ii,iln)=ekb(ii+iln*(iln-1)/2,itypat)
!            end do
             write(nunit, '(3i6,3x,4es15.7)' ) iln,lpsang,iproj,&
&             (ekb0(ii,iln),ii=1,min(nekb,4))
             if(nekb>4)then
               do iekb=5,nekb,4
                 write(nunit, '(21x,4es15.7)' )&
&                 (ekb0(ii,iekb),ii=iekb,min(nekb,iekb+3))
               end do
             end if
           end if
         end do
       end do
       ABI_DEALLOCATE(ekb0)

!      ==============================
!      ============ PAW =============
!      ==============================
     else
       write(nunit, '(a)' )&
&       '  Description of the PAW dataset(s)'
       write(nunit, '(a,i6)' )&
&       '  vrsio8 (for pseudopotentials)=',vrsio8
       write(nunit, '(a,i3)' ) '  usepaw =',usepaw
       do itypat=1,ntypat
         iln0=0
         ABI_ALLOCATE(orbitals,(pawtab(itypat)%basis_size))
         do ilmn=1,pawtab(itypat)%lmn_size
           iln =indlmn(5,ilmn,itypat)
           if (iln>iln0) then
             iln0=iln;orbitals(iln)=indlmn(1,ilmn,itypat)
           end if
         end do
         write(nunit, '(a,i4,a,i3,a,i5)' ) &
&         '  Atom type=',itypat,' basis_size=',pawtab(itypat)%basis_size,&
&         '   lmn_size=',pawtab(itypat)%lmn_size
         write(nunit, '(a,50i2)' ) &
&         '    Basis functions=',orbitals(1:pawtab(itypat)%basis_size)
         write(nunit, '(a,f6.3,a,i2,a,f6.3)' ) &
&         '    r_PAW= ',pawtab(itypat)%rpaw,' shape_type= ',pawtab(itypat)%shape_type,&
&         '  r_shape= ',pawtab(itypat)%rshp
         nekb=0
         ABI_ALLOCATE(dij0,(pawtab(itypat)%lmn2_size))
         ABI_ALLOCATE(i1,(pawtab(itypat)%lmn2_size))
         ABI_ALLOCATE(i2,(pawtab(itypat)%lmn2_size))
         do jlmn=1,pawtab(itypat)%lmn_size
           ij=jlmn*(jlmn-1)/2
           do ilmn=1,jlmn
             if (abs(pawtab(itypat)%dij0(ij+ilmn))>tol16) then
               nekb=nekb+1;i1(nekb)=ilmn;i2(nekb)=jlmn
               dij0(nekb)=pawtab(itypat)%dij0(ij+ilmn)
             end if
           end do
         end do
         write(nunit,'(a,i3,a)') '    Dij0=     (only the ',nekb,' values different from zero)'
         write(nunit,'(2a)') '       i    j     Dij0        i    j     Dij0 ',&
&         '       i    j     Dij0        i    j     Dij0'
         do ii=1,nekb,4
           write(nunit,'(3x,4(1x,i4,1x,i4,1x,es12.5))') (i1(ij),i2(ij),dij0(ij),ij=ii,min(ii+3,nekb))
         end do
         ABI_DEALLOCATE(dij0)
         ABI_DEALLOCATE(i1)
         ABI_DEALLOCATE(i2)
         ABI_DEALLOCATE(orbitals)
       end do

     end if ! NCPP or PAW
   end if ! fullinit==0

!  Now, write the number of blocks
   write(nunit, '(a)' )' '
   write(nunit, '(a)' )' **** Database of total energy derivatives ****'
   write(nunit, '(a,i4)' ) ' Number of data blocks= ',nblok

 end if

end subroutine psddb8
!!***
