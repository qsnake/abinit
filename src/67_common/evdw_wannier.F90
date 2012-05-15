!{\src2tex{textfont=tt}}
!!****f* ABINIT/evdw_wannier
!! NAME
!! evdw_wannier
!!
!! FUNCTION
!!  FIXME: Evaluates the van der Waals correlation energy using maximally
!!         localized Wannier funcitons (MLWF) as proposed by: 
!!         P. L. Silvestrelli in PRL 100:053002 (2008) vdw_xc=10 and 
!!         A. Ambrosetti and P. L. Silvestrelli in PRB 85:073101 (2012) vdw_xc=11.
!!
!! COPYRIGHT
!!  Copyright (C) 2010,2012 ABINIT group (CE and TR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS

!!   nsppol          = Spin polarization.
!!   nwan(nsppol)    = Total number of MLWF in the system per spin component.
!!   mwan            = max[nwan(nsppol)]     
!!   vdw_nfrag       = Number of vdW interating fragments in the unit cell.                     
!!   vdw_supercell(3)     = Distance along each rprimd components for 
!!                          which vdW interactions between MLWF will be taken into account. 
!!   vdw_typfrag(natom)   = Fragment to which each atom belongs to.
!!   vdw_xc               = vdW-WF version.  
!!   rprimd               = Real space primitive translations.    
!!   wann_centres(3,mwan,nsppol) = The centers of MLWFs  in a.u. 
!!   wann_spreads(mwan,nsppol)   = Spread of the MLWFs, in Ang**2. (from wannier90).
!!   xcart           = Coordinates of unit cell atoms in atomic units.    
!!
!! OUTPUT
!!   csix(mwan,mwan,nsppol,nsppol) = dispersion coefficient between each pair of MLWF.          
!!   corrvdw           = van der Waals correction to the energy.
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      mlwfovlp
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine evdw_wannier(csix,corrvdw,mwan,natom,nsppol,nwan,vdw_nfrag,&
& vdw_supercell,vdw_typfrag,vdw_xc,rprimd,wann_centres,wann_spreads,xcart)    

 use m_profiling

use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'evdw_wannier'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
 use interfaces_67_common, except_this_one => evdw_wannier
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: mwan,nsppol,natom,nwan(nsppol)
 integer , intent(in)  :: vdw_nfrag,vdw_supercell(3),vdw_typfrag(natom),vdw_xc
 real(dp), intent(in)  :: rprimd(3,3),wann_centres(3,mwan,nsppol),wann_spreads(mwan,nsppol)
 real(dp), intent(in)  :: xcart(3,natom)
 real(dp), intent(out) :: corrvdw
 real(dp), intent(out), allocatable :: csix(:,:,:,:) 

!Local variables-------------------------------
 integer  :: ii,inx,iny,inz,iwan,jwan,jj,ll,tnwan                      
 integer, allocatable:: ord(:,:)
 real(dp) ::fij,rij,fu,shift
 real(dp), parameter :: a = 20.d0 !Parameter related to the damping function.
 real(dp), parameter :: gama = 4.5d0/(sqrt3**3) !alpha=gama*S**3.
 real(dp), allocatable:: dcenters(:,:,:),rc(:,:),rv(:,:),wanncent(:,:,:),wannspr(:,:)
 real(dp), allocatable:: wc_rec(:,:,:),xi(:,:) 
 character(len=500) :: message                   ! to be uncommented, if needed
! *************************************************************************

!DEBUG                                           ! to be uncommented, if needed
!if(option/=1 .and. option/=2 )then
!write(message,'(a,a,a,a,a,a,i6)') ch10,&
!&  ' evdw_wannier: BUG -',ch10,&
!&  '  The argument option should be 1 or 2,',ch10,&
!&  '  however, option=',option
!call wrtout(std_out,message,'COLL')
!call leave_new('COLL')
!endif
!if(sizein<1)then
!write(message,'(a,a,a,a,a,a,i6)') ch10,&
!&  ' evdw_wannier: BUG -',ch10,&
!&  '  The argument sizein should be a positive number,',ch10,&
!&  '  however, sizein=',sizein
!call wrtout(std_out,message,'COLL')
!call leave_new('COLL')
!endif
!ENDDEBUG
 
 ABI_ALLOCATE(wanncent,(3,mwan,nsppol))
 ABI_ALLOCATE(wannspr,(mwan,nsppol))
 ABI_ALLOCATE(wc_rec,(3,mwan,nsppol))
 ABI_ALLOCATE(ord,(mwan,nsppol))

!The vdW correction is calculated in atomic units:
 wanncent(:,:,:)=wann_centres(:,:,:)/Bohr_Ang
!converting to bohr**2 and then squared 
 wannspr(:,:)=sqrt(wann_spreads(:,:)/Bohr_Ang**2)
!write(std_out,*) "spread of WF",i, "=", wann_spreads(i)
 tnwan=0
 do ii=1,nsppol
   tnwan=tnwan+nwan(ii)
 end do 
!write(std_out,*) 'Number of MLWFs:',ch10
!do ii=1,nsppol
!write(std_out,*) 'nsppol=',ii, 'nwan(nsppol)=',nwan(nsppol),ch10
!end do

 write(std_out,*) 'Original Wannier centres and spreads:',ch10
 do ii=1,nsppol
   write(std_out,*) 'nsppol=',ii,ch10 
   do iwan=1,nwan(nsppol)
     write(std_out,*) (wanncent(jj,iwan,ii),jj=1,3), wannspr(iwan,ii)
   end do
 end do

 if(vdw_nfrag>0)then
   do jj=1,nsppol
     call xredxcart(nwan(jj),-1,rprimd,wanncent(:,1:nwan(jj),jj),wc_rec(:,1:nwan(jj),jj)) 
!    got centers in reduced coor
     do iwan=1,nwan(jj)
       do ii=1,3
         if(wc_rec(ii,iwan,jj)<zero) then
           shift=REAL(CEILING(ABS(wc_rec(ii,iwan,jj))),dp)
           wc_rec(ii,iwan,jj) = wc_rec(ii,iwan,jj)+shift
         end if 
         if(wc_rec(ii,iwan,jj)>one) then                      
           shift=-REAL(INT(wc_rec(ii,iwan,jj)),dp)
           wc_rec(ii,iwan,jj) = wc_rec(ii,iwan,jj)+shift
         end if 
       end do
     end do
     call xredxcart(nwan(jj),1,rprimd,wanncent(:,1:nwan(jj),jj),wc_rec(:,1:nwan(jj),jj))
   end do

!  ====================================================================

   write(std_out,*) 'Wannier centres translated to unit cell and spr:',ch10
   do jj=1,nsppol
     write(std_out,*) 'nsppol=',jj,ch10 
     do iwan=1,nwan(jj)
       write(std_out,*) (wanncent(ii,iwan,jj),ii=1,3), wannspr(iwan,jj)
     end do
   end do
 end if !vdw_nfrag>0

 call order_wannier(mwan,natom,nwan,nsppol,ord,vdw_typfrag,wanncent,xcart)  
!Assing each MLWFs to one fragment, the same as their nearest atom.
 
 write(std_out,*) ch10,'Wannier centres and fragments',ch10
 do ll=1,abs(vdw_nfrag)
   write(std_out,*) 'MLWF centers in fragment',ll,ch10
   do jj=1,nsppol
     do iwan=1,nwan(jj)
       if (ord(iwan,jj)==ll) then
         write(std_out,*) 'X', (Bohr_Ang*wanncent(ii,iwan,jj),ii=1,3)
       end if
     end do 
   end do
 end do


!vdW-WF VERSION 1

 if(vdw_xc==10) then

   ABI_ALLOCATE(dcenters,(3,mwan,nsppol))
   ABI_ALLOCATE(rc,(mwan,nsppol))
   ABI_ALLOCATE(rv,(mwan,nsppol))
!  Calculate intermediate quantities
   do jj=1,nsppol
     do iwan=1, nwan(jj)                                              
       rc(iwan,jj)= three*(0.769d0+half*dlog(wannspr(iwan,jj)))            
!      rv(iwan,jj)= (1.475d0-half_sqrt3*dlog(wannspr(iwan,jj)))*wannspr(iwan,jj) 
       rv(iwan,jj)= (rc(iwan,jj)*wannspr(iwan,jj))/sqrt3 !r_v suggested in JPhysChemA 113:5224   
     end do                                                      
   end do
   corrvdw=0.0d0  !Initializing the vdW correction energy.
   
   ABI_ALLOCATE(csix,(mwan,mwan,nsppol,nsppol))

   do ii=1,nsppol
     do jj=1,nsppol  
       do iwan=1,nwan(ii)
         do jwan=1,nwan(jj)
           
           call getFu(wannspr(iwan,ii),wannspr(jwan,jj),rc(iwan,ii),rc(jwan,jj),fu)

           csix(iwan,jwan,ii,jj)=( ( ((wannspr(iwan,ii))**1.5d0)*&
&           (wannspr(jwan,jj)**three))/(two*(three**1.25d0) ) )*fu

         end do
       end do
     end do
   end do
   
   if (nsppol == 1) then
     csix(:,:,:,:)=sqrt2*csix(:,:,:,:)  !For non spin polarized systems
   end if 


!  DEBUG
!  write(std_out,*) ch10,'gamma=',gama,ch10
!  write(std_out,*) ch10,'C6ij coefficients',ch10,'index i(j) from mlwf in fragment 1(2):',ch10
!  do ii=1,nsppol
!  do jj=1,nsppol
!  do iwan=1,nwan(ii)
!  write(std_out,*) (csix(iwan,jwan,ii,jj),jwan=1,nwan(jj)),ch10,&
!  (csixalt(iwan,jwan,ii,jj),jwan=1,nwan(jj))
!  end do
!  end do  
!  end do 
!  END DEBUG
   
!  test   k=0   
   do ii=1,nsppol
     do iwan=1,nwan(ii)
       do inx=-abs(vdw_supercell(1)),abs(vdw_supercell(1))
         do iny=-abs(vdw_supercell(2)),abs(vdw_supercell(2))
           do inz=-abs(vdw_supercell(3)),abs(vdw_supercell(3))
             do jj=1,nsppol
               do jwan=1,nwan(jj) 
                 
                 if(inx==0.and.iny==0.and.inz==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle 
!                This avoids intrafragment vdW interactions.
                 if(vdw_supercell(1)<=0.and.inx==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
                 if(vdw_supercell(2)<=0.and.iny==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
                 if(vdw_supercell(3)<=0.and.inz==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
!                Last three conditions allow proper treatment of layered systems.
                 
                 dcenters(:,jwan,jj) = (real(inx,dp))*rprimd(:,1)+(real(iny,dp))*rprimd(:,2)+&
&                 (real(inz,dp))*rprimd(:,3)+wanncent(:,jwan,jj) 
                 rij=sqrt((dcenters(1,jwan,jj)-wanncent(1,iwan,ii))**2+&
&                 (dcenters(2,jwan,jj)-wanncent(2,iwan,ii))**2+&
&                 (dcenters(3,jwan,jj)-wanncent(3,iwan,ii))**2) 

                 fij=one/(one+exp(-a*(rij/(rv(iwan,ii)+rv(jwan,jj))-one))) !Damping function.  
                 
                 corrvdw = corrvdw - csix(iwan,jwan,ii,jj)*fij/(two*(rij**6)) !making the sum of eq(4) of 
!                JPhysChemA 113:5224-5234. Each term is divided by two because
!                we are counting twice within the unit cell, also the 
!                interactions with neighbor cells are properly acounted for in 
!                this way.
                 
!                write(std_out,*) 'i=',iwan, 'j=',jwan, 'C6ij=', csix(iwan,jwan)
!                write(std_out,*) 'inx=',inx, "iny=",iny, "inz=",inz, "Evdw=",&
!                & -(csix(iwan,jwan)*fij/(two*rij**6))*Ha_ev*ten**3 
!                write(std_out,*) 'rnl=',rnl
               end do
             end do
           end do
         end do
       end do 
     end do
   end do

   ABI_DEALLOCATE(dcenters)
   ABI_DEALLOCATE(rc)
   ABI_DEALLOCATE(rv)
   
   write(message, '(2a,i2,2a,f12.3,2a,f12.3,a)' )ch10,&
&   ' vdw_xc : ',10,ch10,&                                                              
&   ' van der Waals correction(Ha):',   corrvdw,ch10,&                            
&   ' van der Waals correction(meV):',   corrvdw*Ha_ev*ten**3,ch10                                                           
   call wrtout(std_out,message,'COLL')                                                     
   call wrtout(ab_out,message,'COLL')

 end if

!vdW-WF VERSION 2: Phys. Rev. B. 85:073101 (2012)

 if (vdw_xc==11) then
   
   ABI_ALLOCATE(dcenters,(3,mwan,nsppol))
   ABI_ALLOCATE(rv,(mwan,nsppol))
   ABI_ALLOCATE(xi,(mwan,nsppol))

!  Calculate intermediate quantities
   do jj=1,nsppol
     do iwan=1, nwan(jj)                                              
       rv(iwan,jj)= ( (1.20d0/Bohr_Ang)*wannspr(iwan,jj) )/sqrt3
       write(std_out,*) 'rv(iwan,jj)=',rv(iwan,jj),ch10
     end do                                                      
   end do

   ABI_ALLOCATE(csix,(mwan,mwan,nsppol,nsppol))
!  C6 coefficients between WF
   csix(:,:,:,:) = 0.0d0
   corrvdw = 0.0d0

   call ovlp_wann(mwan,nwan,nsppol,ord,wanncent,wannspr,xi)

!  DEBUG
   write(std_out,*)ch10,'xi(iwan,isspol)=',ch10 
   do jj=1,nsppol
     write(std_out,*) (xi(iwan,jj),iwan=1,nwan(jj))
   end do
!  END DEBUG

   do ii=1,nsppol
     do jj=1,nsppol  
       do iwan=1,nwan(ii)
         do jwan=1,nwan(jj)
           
           csix(iwan,jwan,ii,jj)=onehalf*( (wannspr(iwan,ii)*wannspr(jwan,jj))**three )*&
&           ((xi(iwan,ii)*xi(jwan,jj))*gama**onehalf)/( sqrt(xi(iwan,ii))*&
&           wannspr(iwan,ii)**onehalf + sqrt(xi(jwan,jj))*wannspr(jwan,jj)**onehalf )
           
         end do
       end do
     end do
   end do
   
   if (nsppol == 1) then
     csix(:,:,:,:)=sqrt2*csix(:,:,:,:)  !For non spin polarized systems
   end if 

!  DEBUG
   write(std_out,*) ch10,'C6ij coefficients:',ch10
   do ii=1,nsppol
     do jj=1,nsppol
       do iwan=1,nwan(ii)
         write(std_out,*) (csix(iwan,jwan,ii,jj),jwan=1,nwan(jj))
       end do
     end do  
   end do 
!  END DEBUG
   do ii=1,nsppol
     do iwan=1,nwan(ii)
       do inx=-abs(vdw_supercell(1)),abs(vdw_supercell(1))
         do iny=-abs(vdw_supercell(2)),abs(vdw_supercell(2))
           do inz=-abs(vdw_supercell(3)),abs(vdw_supercell(3))
             do jj=1,nsppol
               do jwan=1,nwan(jj) 
                 
                 if(inx==0.and.iny==0.and.inz==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle 
!                This avoids intrafragment vdW interactions.
                 if(vdw_supercell(1)<=0.and.inx==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
                 if(vdw_supercell(2)<=0.and.iny==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
                 if(vdw_supercell(3)<=0.and.inz==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
!                Last three conditions allow proper treatment of layered systems.
                 
                 dcenters(:,jwan,jj) = (real(inx,dp))*rprimd(:,1)+(real(iny,dp))*rprimd(:,2)+&
&                 (real(inz,dp))*rprimd(:,3)+wanncent(:,jwan,jj) 
                 rij=sqrt((dcenters(1,jwan,jj)-wanncent(1,iwan,ii))**2+&
&                 (dcenters(2,jwan,jj)-wanncent(2,iwan,ii))**2+&
&                 (dcenters(3,jwan,jj)-wanncent(3,iwan,ii))**2) 

                 fij=one/(one+exp(-a*(rij/(rv(iwan,ii)+rv(jwan,jj))-one))) !Damping function.  
!                DEBUG
!                write(std_out,*) 'f_i,j=',fij,ch10
!                END DEBUG               
                 corrvdw = corrvdw - csix(iwan,jwan,ii,jj)*fij/(two*(rij**6)) !making the sum of eq(4) of 
!                JPhysChemA 113:5224-5234. 
               end do
             end do
           end do
         end do
       end do 
     end do
   end do

   write(message, '(2a,i2,2a,f12.3,2a,f12.3,a)' )ch10,&
&   ' vdw_xc : ',11,ch10,&                                                              
&   ' van der Waals correction(Ha):',   corrvdw,ch10,&                            
&   ' van der Waals correction(meV):',   corrvdw*Ha_ev*ten**3,ch10                                                           
   call wrtout(std_out,message,'COLL')                                                     
   call wrtout(ab_out,message,'COLL')

   ABI_DEALLOCATE(dcenters)
   ABI_DEALLOCATE(rv)
   ABI_DEALLOCATE(xi)
 end if

 ABI_DEALLOCATE(ord)
 ABI_DEALLOCATE(wanncent)
 ABI_DEALLOCATE(wannspr)

end subroutine evdw_wannier
!!***

!!****f* ABINIT/getFu
!! NAME
!! getFu
!!
!! FUNCTION
!!  Performs double integral needed to evaluate C6 
!!  coefficients. Eq. (9) in J.Phys.Chem. 113:5224
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2012 ABINIT group (CE and TR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT

!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      evdw_wannier
!!
!! CHILDREN
!!
!! SOURCE

 subroutine getFu(sn,sl,rn,rl,fu) ! sn-->spread(n), sl-->spread(l), rn --> rc(n), rl --> rc(l) 

 use m_profiling
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getFu'
 use interfaces_32_util
!End of the abilint section

 implicit none
 real(dp),intent(in)::sn,sl,rn,rl
 real(dp),intent(out)::fu
 !local variables
 integer::nx,ny,ix,iy
 real(dp)::deltax,deltay
 real(dp)::beta,xc,yc,y,x
 real(dp),allocatable::arg1(:),res1(:),arg2(:),res2(:)

! *************************************************************************

 ny=100
 nx=100
 beta=(sn/sl)**(1.5d0)
 xc=rn
 yc=rl
 deltax=xc/(real(nx,dp)-1.d0)
 deltay=yc/(real(ny,dp)-1.d0)
 
 ABI_ALLOCATE(arg1,(ny))
 ABI_ALLOCATE(res1,(ny))
 ABI_ALLOCATE(arg2,(nx))
 ABI_ALLOCATE(res2,(nx))

 do ix=1,nx

   x=deltax*(real(ix,dp)-1.d0)
   
   do iy=1,ny
     y=deltay*(real(iy,dp)-1.d0)
     arg1(iy)=( (y**2.d0)*exp(-y) )/( (exp(-x)/beta) + exp(-y) )
   end do

   call simpson_int(ny,deltay,arg1,res1)
   arg2(ix)=(x**2.d0)*exp(-x)*res1(ny)

 end do

 call simpson_int(nx,deltax,arg2,res2)
 
 Fu = res2(nx)
 
 ABI_DEALLOCATE(arg1)
 ABI_DEALLOCATE(res1)
 ABI_DEALLOCATE(arg2)
 ABI_DEALLOCATE(res2)
end subroutine getFu
!!*** 
      
!!****f* ABINIT/order_wannier
!! NAME
!! order_wannier
!!   
!! FUNCTION     
!!  Assign each MLWF with a corresponding fragment of atoms, according 
!!  to vdw_typfrag array. Assignation is done by evaluating the distance
!!  from each MLWF center to the unit cell atoms. MLWFs belong to the 
!!  same fragment as their nearest atom.    
!! COPYRIGHT
!!  Copyright (C) 2010-2012 ABINIT group (CE and TR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT

!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      evdw_wannier
!!
!! CHILDREN
!!
!! SOURCE
subroutine order_wannier(mwan,natom,nwan,nsppol,ord,vdw_typfrag,wanncent,xcart)

 use m_profiling

   use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'order_wannier'
!End of the abilint section

   implicit none
!Arguments
   integer, intent(in)    :: mwan,natom,nwan(nsppol),nsppol,vdw_typfrag(natom)
   integer, intent(inout) :: ord(mwan,nsppol)
   real(dp),intent(in)    :: wanncent(3,mwan,nsppol),xcart(3,natom)
!Local variables
   integer :: ii,jj,ll
   real(dp):: dis,mindi

! *************************************************************************

 do ll=1,nsppol 
   do ii=1,nwan(ll)
     mindi=sqrt( dot_product(wanncent(:,ii,ll),wanncent(:,ii,ll))+dot_product(xcart(:,1),xcart(:,1))&
&     -2*(dot_product(wanncent(:,ii,ll),xcart(:,1))) ) 
     ord(ii,ll)=vdw_typfrag(1)
     do jj=2,natom
       dis=sqrt( dot_product(wanncent(:,ii,ll),wanncent(:,ii,ll))+dot_product(xcart(:,jj),xcart(:,jj))&
&       -2*(dot_product(wanncent(:,ii,ll),xcart(:,jj))) ) 
       if(dis<=mindi) then
         mindi=dis
         ord(ii,ll)=vdw_typfrag(jj)
       end if
     end do
   end do
 end do
 
 end subroutine order_wannier
!!***

!!****f* ABINIT/ovlp_wann
!! NAME
!! ovlp_wann
!!   
!! FUNCTION     
!!  Evaluate volumen reduction of MLWFs
!!  due to intrafragment overlapping 
!!  
!! COPYRIGHT
!!  Copyright (C) 2011-2012 ABINIT group (CE)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT

!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      evdw_wannier
!!
!! CHILDREN
!!
!! SOURCE
subroutine ovlp_wann(mwan,nwan,nsppol,ord,wanncent,wannspr,xi)

 use m_profiling

   use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ovlp_wann'
!End of the abilint section

   implicit none
!Arguments
   integer, intent(in)  :: mwan,nwan(nsppol),nsppol,ord(mwan,nsppol)
   real(dp),intent(in)  :: wanncent(3,mwan,nsppol),wannspr(mwan,nsppol)
   real(dp), intent(out)  :: xi(mwan,nsppol)
!Local variables
   integer :: ii,iwan,ix,iy,iz,jj,jwan,neigh,steps
   real(dp):: dis,disi,discent,veff,vfree
   integer, allocatable :: intsec(:,:,:,:)
   real(dp), allocatable :: rpoint(:)
   real(dp), parameter :: delt = 0.05d0 !Bohr, spatial mesh (1D) step 
! *************************************************************************

 ABI_ALLOCATE(intsec,(mwan,nsppol,mwan,nsppol))
 ABI_ALLOCATE(rpoint,(3))
 intsec(:,:,:,:) = 0 
 xi(:,:) = 0.0d0
!detecting WF intersecting neighbors 
 do ii=1,nsppol
   do iwan=1,nwan(ii) 
     do jj=1,nsppol
       do jwan=1,nwan(jj)
         dis = 0.0d0
         if (ord(jwan,jj)==ord(iwan,ii)) then 

           dis=sqrt(  dot_product(wanncent(:,iwan,ii),wanncent(:,iwan,ii))+&
&           dot_product(wanncent(:,jwan,jj),wanncent(:,jwan,jj))&
&           -2*( dot_product(wanncent(:,iwan,ii),wanncent(:,jwan,jj)) )  ) 
           
           if ( ii == jj ) then 
             if ( dis<=(wannspr(iwan,ii)+wannspr(jwan,jj)).and.iwan/=jwan ) then
               intsec(iwan,ii,jwan,jj) = 1
             end if
           end if
           if ( ii /= jj) then
             if ( dis<=(wannspr(iwan,ii)+wannspr(jwan,jj)) ) then
               intsec(iwan,ii,jwan,jj) = 1
             end if
           end if

         end if
       end do
     end do
   end do
 end do

!DEBUG
 write(std_out,*) 'intsec(iwan,ii,kk,jj)=',ch10
 do ii=1,nsppol
   do iwan=1,nwan(ii)
     do jj=1,nsppol
       write(std_out,*) (intsec(iwan,ii,jwan,jj),jwan=1,nwan(jj)),ch10
     end do
   end do
 end do
!END DEBUG
!Determining both free and effective volumes.
!Eqs (6) and (7) in PRB 85:073101. 
!Creation of grids around each WF centre. 
!Calculation of intersection volumes. 
 do ii = 1,nsppol 
   do iwan = 1,nwan(ii)
!    Spatial meshes and volume parameters
     steps=NINT(wannspr(iwan,ii)/delt)
     vfree = 0
     veff = 0
     rpoint(:) = 0.0d0 
     do iz=-steps,steps
       do iy=-steps,steps
         do ix=-steps,steps
           neigh = 0
           rpoint(1) = wanncent(1,iwan,ii) + ix*delt
           rpoint(2) = wanncent(2,iwan,ii) + iy*delt   
           rpoint(3) = wanncent(3,iwan,ii) + iz*delt
           discent = sqrt( dot_product(wanncent(:,iwan,ii),wanncent(:,iwan,ii))&
&           +dot_product( rpoint(:),rpoint(:) )&
&           -2*( dot_product( wanncent(:,iwan,ii),rpoint(:) ) ) ) 
           if (discent > wannspr(iwan,ii)) cycle
           if (discent <= wannspr(iwan,ii)) then

             neigh = 1
             do jj = 1,nsppol
               do jwan = 1,nwan(jj) 
                 if ( intsec(iwan,ii,jwan,jj) == 0 ) cycle 
                 if ( intsec(iwan,ii,jwan,jj) == 1 ) then
                   disi = sqrt( dot_product(rpoint(:),rpoint(:))&
&                   +dot_product( wanncent(:,jwan,jj),wanncent(:,jwan,jj) )&
&                   -2*( dot_product(rpoint(:),wanncent(:,jwan,jj)) ) )            
                   if (disi <= wannspr(jwan,jj)) then
                     neigh = neigh + 1
                   end if 
                 end if
               end do
             end do 
             if (nsppol==1) then 
               veff = veff + 1/(real(neigh,dp)**2)
             end if
             if (nsppol==2) then
               veff = veff + 1/real(neigh,dp)
             end if 
             vfree = vfree + 1/real(neigh,dp)
           end if
         end do    
       end do
     end do 
!    write(std_out,*) 'iwan=',iwan,'ii=',ii,ch10
!    write(std_out,*) 'vfree=',vfree,'neigh=',neigh,'veff=',veff,ch10
     xi(iwan,ii) = veff/vfree  
!    write(std_out,*) 'xi(iwan,ii)=',xi(iwan,ii),ch10
   end do
 end do

 ABI_DEALLOCATE(intsec)
 ABI_DEALLOCATE(rpoint)
 end subroutine ovlp_wann
!!***
