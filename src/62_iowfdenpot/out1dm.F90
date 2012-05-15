!{\src2tex{textfont=tt}}
!!****f* ABINIT/out1dm
!! NAME
!! out1dm
!!
!! FUNCTION
!! Output the 1 dimensional mean of potential and density
!! on the three reduced axis.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  fnameabo_app_1dm=name of the file in which the data is written, appended with _1DM
!!  natom=number of atoms in unit cell
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in unit cell.
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=real space dimensional primitive translations (bohr)
!!  typat(natom)=type integer for each atom in cell
!!  ucvol=unit cell volume in bohr**3.
!!  vtrial(nfft,nspden)=INPUT Vtrial(r).
!!  xred(3,natom)=reduced coordinates of atoms
!!  znucl(ntypat)=real(dp), nuclear number of atom type
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  data written in file fnameabo_app_1dm
!!  xred is not modified by xredxcart evenf if it is declared inout (TD).
!!
!! NOTES
!!
!! PARENTS
!!      dielmt2,ladielmt,outscfcv
!!
!! CHILDREN
!!      atmdata,leave_new,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine out1dm(fnameabo_app_1dm,natom,nfft,ngfft,nspden,ntypat,&
&  rhor,rprimd,typat,ucvol,vtrial,xred,znucl)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'out1dm'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat
 real(dp),intent(in) :: ucvol
 character(len=fnlen),intent(in) :: fnameabo_app_1dm
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3),vtrial(nfft,nspden)
 real(dp),intent(in) :: znucl(ntypat)
 real(dp),intent(inout) :: xred(3,natom)

!Local variables-------------------------------
! character(len=2), parameter :: symbol(92)=(/' H','He',        &
!&   'Li','Be',' B',' C',' N',' O',' F','Ne',   &
!&   'Na','Mg','Al','Si',' P',' S','Cl','Ar',   &
!&   ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni',&
!&        'Cu','Zn','Ga','Ge','As','Se','Br','Kr',     &
!&   'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',&
!&        'Ag','Cd','In','Sn','Sb','Te',' I','Xe',     &
!&   'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',&
!&                       'Tb','Dy','Ho','Er','Tm','Yb',&
!&             'Lu','Hf','Ta',' W','Re','Os','Ir','Pt',&
!&        'Au','Hg','Tl','Pb','Bi','Po','At','Rn',     &
!&   'Fr','Ra','Ac','Th','Pa',' U'/)
!scalars
 integer :: ia,ib,idim,ifft,islice,ispden,na,nb,ndig,nslice,nu
 real(dp) :: amu,global_den,global_pot,rcov
 character(len=2) :: symbol
 character(len=500) :: message
!arrays
 real(dp),allocatable :: lin_den(:),mean_pot(:),reduced_coord(:),xcart(:,:)
 character(len=8),allocatable :: iden(:)

! *************************************************************************

 if(nspden==4)then
   write(std_out,*)' out1dm : does not work yet for nspden=4: do nothing and return to normal execution'
   return
 end if

!Initialize the file
 write(message, '(a,a)' ) ' io1dm : about to open file ',fnameabo_app_1dm
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 open (unit=tmp_unit,file=fnameabo_app_1dm,status='unknown',form='formatted')
 rewind(tmp_unit)

 write(message, '(a,a)' ) ch10,' ABINIT package : 1DM file '
 call wrtout(tmp_unit,message,'COLL')

 write(message, '(a,a)' )ch10,' Primitive vectors of the periodic cell (bohr)'
 call wrtout(tmp_unit,message,'COLL')
 do nu=1,3
   write(message, '(1x,a,i1,a,3f10.5)' ) '  R(',nu,')=',rprimd(:,nu)
   call wrtout(tmp_unit,message,'COLL')
 end do

 write(message, '(a,a)' ) ch10,&
& ' Atom list        Reduced coordinates          Cartesian coordinates (bohr)'
 call wrtout(tmp_unit,message,'COLL')

!Set up a list of character identifiers for all atoms : iden(ia)
 ABI_ALLOCATE(iden,(natom))
 do ia=1,natom
   call atmdata(amu,rcov,symbol,znucl(typat(ia)))
   ndig=int(log10(dble(ia)+0.5_dp))+1
   if(ndig==1) write(iden(ia), '(a,a,i1,a)' )symbol,'(',ia,')   '
   if(ndig==2) write(iden(ia), '(a,a,i1,a)' )symbol,'(',ia,')  '
   if(ndig==3) write(iden(ia), '(a,a,i1,a)' )symbol,'(',ia,') '
   if(ndig==4) write(iden(ia), '(a,a,i1,a)' )symbol,'(',ia,')'
   if(ndig>4)then
     write(message, '(a,i8)' )&
&     ' out1dm : cannot handle more than 9999 atoms, while natom=',natom
     call wrtout(std_out,message,'COLL')
     close(tmp_unit)
     call leave_new('COLL')
   end if
 end do

!Compute cartesian coordinates, and print reduced and cartesian coordinates
 ABI_ALLOCATE(xcart,(3,natom))
 call xredxcart(natom,1,rprimd,xcart,xred)
 do ia=1,natom
   write(message, '(a,a,3f10.5,a,3f10.5)' ) &
&   '   ',iden(ia),xred(1:3,ia),'    ',xcart(1:3,ia)
   call wrtout(tmp_unit,message,'COLL')
 end do
 ABI_DEALLOCATE(iden)
 ABI_DEALLOCATE(xcart)

 do idim=1,3

   nslice=ngfft(idim)

!  Dummy initialisation of na, nb and ifft.
   na=0 ; nb=0; ifft=0
   select case(idim)
     case(1)
       na=ngfft(2) ; nb=ngfft(3)
     case(2)
       na=ngfft(1) ; nb=ngfft(3)
     case(3)
       na=ngfft(1) ; nb=ngfft(2)
   end select

   ABI_ALLOCATE( reduced_coord,(nslice))
   ABI_ALLOCATE(mean_pot,(nslice))
   ABI_ALLOCATE(lin_den,(nslice))

   do ispden=1,nspden

     if(ispden==1)then
       write(message, '(a,a,a)' ) ch10,'===========',&
&       '====================================================================='
       call wrtout(tmp_unit,message,'COLL')
     end if

     select case(idim)
       case(1)
         write(message, '(a)' )' Projection along the first dimension '
       case(2)
         write(message, '(a)' )' Projection along the second dimension '
       case(3)
         write(message, '(a)' )' Projection along the third dimension '
     end select
     call wrtout(tmp_unit,message,'COLL')

     if(nspden==2)then
       select case(ispden)
         case(1)
           write(message, '(a)' )' Spin up '
         case(2)
           write(message, '(a)' )' Spin down '
       end select
       call wrtout(tmp_unit,message,'COLL')
     end if

     write(message, '(2a)' ) ch10,&
&     '     Red. coord. Mean KS potential  Linear density  '
     call wrtout(tmp_unit,message,'COLL')

     write(message, '(a)' ) &
&     '                  (Hartree unit)   (electron/red. unit)'
     call wrtout(tmp_unit,message,'COLL')

     global_pot=0.0_dp
     global_den=0.0_dp
     do islice=1,nslice
       reduced_coord(islice)=(islice-1)/(dble(nslice))
       mean_pot(islice)=0.0_dp
       lin_den(islice)=0.0_dp
       do ib=1,nb
         do ia=1,na
           if(idim==1) ifft = islice + nslice*( ia    -1 + na    *(ib    -1) )
           if(idim==2) ifft = ia     + na    *( islice-1 + nslice*(ib    -1) )
           if(idim==3) ifft = ia     + na    *( ib    -1 + nb    *(islice-1) )
           mean_pot(islice)=mean_pot(islice)+vtrial(ifft,ispden)
           lin_den(islice)=lin_den(islice)+rhor(ifft,ispden)
         end do
       end do
       mean_pot(islice)=mean_pot(islice)/dble(na*nb)
       lin_den(islice)=lin_den(islice)/dble(na*nb)*ucvol
       global_pot=global_pot+mean_pot(islice)
       global_den=global_den+ lin_den(islice)
     end do
     global_pot=global_pot/dble(nslice)
     global_den=global_den/dble(nslice)

     do islice=1,ngfft(idim)
       write(message, '(i3,f10.4,es20.6,es16.6)' )&
&       islice,reduced_coord(islice),mean_pot(islice),lin_den(islice)
       call wrtout(tmp_unit,message,'COLL')
     end do

     write(message, '(a,a,es15.6,es16.6)' ) ch10,&
&     ' Cell mean       :',global_pot,global_den
     call wrtout(tmp_unit,message,'COLL')


!    End of the loop on spins
   end do

   ABI_DEALLOCATE(reduced_coord)
   ABI_DEALLOCATE(mean_pot)
   ABI_DEALLOCATE(lin_den)

!  End of the loops on the three dimensions
 end do

 write(message, '(a,a,a)' ) ch10,'===========',&
& '====================================================================='
 call wrtout(tmp_unit,message,'COLL')

 close(tmp_unit)

end subroutine out1dm
!!***
