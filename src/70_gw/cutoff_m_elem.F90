!{\src2tex{textfont=tt}}
!!****f* ABINIT/cutoff_m_elem
!! NAME
!! cutoff_m_elem
!!
!! FUNCTION
!!  Main routine for the output of the optical matrix elements including a cutoff in real space.
!!
!! INPUTS
!!  ep= datatype gathering differening parameters related to the calculation of the inverse dielectric matrix
!!  ep%npwvec=dimension of igffttt
!!  energy(ep%nbnds,ep%nkibz,ep%nsppol)=KS energies
!!  ep%npwwfn=number of planewaves for wavefunctions (input variable)
!!  ep%nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(ep%nkibz,ep%nbnds,ep%nsppol)=occupation numbers, for each k point in IBZ, and each band
!!  Wf%ug(ep%npwwfn,my_minb:my_maxb,ep%nkibz,ep%nsppol)= (optional) wfs in real space, for each band treated by this processor
!!  Wf%wfr(nfftot,my_minb:my_maxb,ep%nkibz,ep%nsppol)= (optional)  wfs in G space for each band treated by this processor
!!  Wf_val%ug(ep%npwwfn,nbvw,ep%nkibz,ep%nsppol)= (optional) array containing fully and partially occupied states in G space
!!  Wf_val%wfr(nfftot,nbvw,ep%nkibz,ep%nsppol) = (optional) array containing unoccupied states in real space
!!  z0 = coordinate of the first step in the cutoff function
!!  wdth = coordinate of the second (descending) step in the cutoff
!!  function(width of the step)
!!
!! OUTPUT
!! The output will be stored in a file called output.dat
!! "direction" is a variable taken from the input files (userrc). 
!! If it is set to 1 (x axis), we consider cutoff perpendicular to the x axis; 
!! idem if is set to 2 (y axis), or 3 (z axis). If it is set to zero (default value), we consider no real space cutoff. 
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine cutoff_m_elem(ep,kmesh,gvec,Wf,energy,z0,wdth,occ,direction,gprimd)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_gwdefs

 use m_bz_mesh,  only : bz_mesh_type
 use m_wfs,      only : wfs_descriptor

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cutoff_m_elem'
 use interfaces_70_gw, except_this_one => cutoff_m_elem
!End of the abilint section

implicit none

!Arguments ------------------------------------

 type(bz_mesh_type),target,intent(in) :: kmesh
 type(wfs_descriptor),optional,intent(in) :: Wf
 type(epsilonm1_parameters),intent(in) :: ep
!scalars
 integer,intent(in) :: direction
 real(dp),intent(in) :: z0, wdth 

 complex(gwpc) :: res(3), vec_cart(3)
!arrays
 integer,intent(in) :: gvec(3,ep%npwvec) 
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: energy(ep%nbnds,ep%nkibz,ep%nsppol)
 real(dp),intent(in) :: occ(ep%nbnds,ep%nkibz,ep%nsppol)

!Local variables-------------------------------
 real(dp) :: b1(3), b2(3), b3(3)
 complex(gwpc),allocatable :: wfg1(:),wfg2(:)
 integer :: ivb,icb,ik,is,nval
 real(dp) :: diffen

!************************************************************************
!BEGIN EXECUTABLE SECTION

 b1=two_pi*gprimd(:,1)
 b2=two_pi*gprimd(:,2)
 b3=two_pi*gprimd(:,3)


! warning:the following check is valid ONLY for semiconductors!
nval=0
do ivb=1,ep%nbnds
if (abs(occ(ivb,1,1)) > 0.001) nval=nval+1
enddo

open(78,file='output_matrix_elements.dat')
ABI_ALLOCATE(wfg1,(Wf%npwwfn))
ABI_ALLOCATE(wfg2,(Wf%npwwfn))

write(78,*) 'matrixelements (NF) for a number of points'
write(78,*) 'F  NF: flag tmetal'
write(78,*) 'T  NF: flag for the long (T) or short (F) format'
if (abs(wdth-1) < 0.0001) then
 write(78,*) 'F  NF: flag for the cutoff function'
else
 write(78,*) 'T  NF: flag for the cutoff function'
endif
write(78,*) '         ', kmesh%nibz,' NF: number of groups'


!
if (direction == 1 .or.direction == 2 .or. direction == 3) then

do is=1, ep%nsppol
do ik=1,kmesh%nibz
 nval=0
 do ivb=1,ep%nbnds
 if (abs(occ(ivb,ik,1))>0.01) nval=nval+1
 enddo
 write(78,*) '1 k-points in group number  ', ik,' (NF)' 
 write(78,'(1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.5,1x, "(NF)")')  kmesh%ibz(1,ik), kmesh%ibz(2,ik), kmesh%ibz(3,ik),kmesh%wt(ik)
 write(78,*)  '   ',nval, ep%nbnds
do ivb=1,ep%nbnds
do icb=1,ep%nbnds
 if(occ(ivb,ik,is) - occ(icb,ik,is) <= 1.0e-04) CYCLE 
 diffen =  energy(icb,ik,is)*Ha_eV - energy(ivb,ik,is)*Ha_eV 
    call matrixelements(Wf%npwwfn,Wf%Wave(ivb,ik,is)%ug,Wf%Wave(icb,ik,is)%ug,gvec,kmesh%ibz(:,ik),res)
    vec_cart(:) = res(1)*b1(:) + res(2)*b2(:) + res(3)*b3(:)
    write(78,'(1x,i4,1x,i4,1x,i4,1x,e16.8,1x,e16.8,1x,e16.8)') ik,ivb,icb,diffen,energy(icb,ik,is)*Ha_eV,energy(ivb,ik,is)*Ha_eV
    write(78,'(1x,e16.8,1x,e16.8,1x,e16.8,1x,e16.8,1x,e16.8,1x,e16.8)') &
     &real(vec_cart(1)),aimag(vec_cart(1)),real(vec_cart(2)),aimag(vec_cart(2)),real(vec_cart(3)),aimag(vec_cart(3))

    call matrixelements_cutoff(Wf%npwwfn,Wf%Wave(ivb,ik,is)%ug,Wf%Wave(icb,ik,is)%ug,gvec,kmesh%ibz(:,ik),z0,wdth,direction,res)
    vec_cart(:) = res(1)*b1(:) + res(2)*b2(:) + res(3)*b3(:)
    write(78,'(1x,e16.8,1x,e16.8,1x,e16.8,1x,e16.8,1x,e16.8,1x,e16.8)') &
     &real(vec_cart(1)),aimag(vec_cart(1)),real(vec_cart(2)),aimag(vec_cart(2)),real(vec_cart(3)),aimag(vec_cart(3)) 
enddo
enddo
enddo
enddo


else  ! (abs(wdth-1) < 0.0001 .or. direction == 0) then
do is=1, ep%nsppol
do ik=1,kmesh%nibz

nval=0
do ivb=1,ep%nbnds
if (abs(occ(ivb,ik,1))>0.01) nval=nval+1
enddo

write(78,*) '1 k-points in group number  ', ik,' (NF)' 
write(78,'(1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.5,1x, "(NF)")')  kmesh%ibz(1,ik), kmesh%ibz(2,ik), kmesh%ibz(3,ik),kmesh%wt(ik)
write(78,*)  '   ',nval, ep%nbnds

do ivb=1,ep%nbnds
do icb=1,ep%nbnds
   
  if(occ(ivb,ik,is) - occ(icb,ik,is) <= 1.0e-04) CYCLE ! check in order to take into account only
                                                       ! transitions from vb to cb
diffen =  energy(icb,ik,is)*Ha_eV - energy(ivb,ik,is)*Ha_eV 
                                                       
    call matrixelements(Wf%npwwfn,Wf%Wave(ivb,ik,is)%ug,Wf%Wave(icb,ik,is)%ug,gvec,kmesh%ibz(:,ik),res)
    !   ep%npwwfn dovrebbe essere uguale a Wf%npwwfn  !?
    vec_cart(:) = res(1)*b1(:) + res(2)*b2(:) + res(3)*b3(:)
    write(78,'(1x,i4,1x,i4,1x,i4,1x,e16.8,1x,e16.8,1x,e16.8)') ik,ivb,icb,diffen,energy(icb,ik,is)*Ha_eV,energy(ivb,ik,is)*Ha_eV
    write(78,'(1x,e16.8,1x,e16.8,1x,e16.8,1x,e16.8,1x,e16.8,1x,e16.8)') &
     &real(vec_cart(1)),aimag(vec_cart(1)),real(vec_cart(2)),aimag(vec_cart(2)),real(vec_cart(3)),aimag(vec_cart(3)) 
enddo
enddo
enddo
enddo

endif

ABI_DEALLOCATE(wfg1)
ABI_DEALLOCATE(wfg2)

close(78)

end subroutine  cutoff_m_elem
!!***


!!****f* ABINIT/matrixelements
!! NAME
!! matrixelements
!!
!! FUNCTION
!! This subroutine calculates the matrix elements:
!! \hat{P}_{vc}(k) := \sum_{G} F(g_{z}) conjug(C_v(G) \hbar (k + G) C_c(G) 
!!
!! INPUTS
!!  npwwfn=number of G vectors for wavefunctions
!!  gvec(3,npwwfn)= reduced coordinates of G vectors
!!  wfg1(npwwfn),wfg1(npwwfn)= bra and ket in reciprocal space
!!  kpoint(3)= reduced coordinates of the k point of interest for the calculation of the matrix el. P^{~}_{vc}(k)
!!
!! OUTPUT
!! res =  inabla(3)= i <wfg1|\nabla|wfg2>
!!
!! PARENTS
!!      cutoff_m_elem
!!
!! CHILDREN
!!
!! SOURCE
subroutine matrixelements(npwwfn,wfg1,wfg2,gvec,kpoint,res)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'matrixelements'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
 integer,intent(in) :: npwwfn
 integer,intent(in) :: gvec(3,npwwfn)
 complex(gwpc),intent(in) :: wfg1(npwwfn),wfg2(npwwfn)
 complex(gwpc),intent(out) :: res(3)
 real(dp),intent(in) :: kpoint(3)
!Local variables-------------------------------
 integer :: ig1
 complex(gwpc) :: ct,uimag

! *************************************************************************

 uimag=dcmplx(0.d0,1.d0)
 res(:)=(0.0_gwp,0.0_gwp)
 do ig1=1,npwwfn
     ct=CONJG(wfg1(ig1))*wfg2(ig1)
   res(:)=res(:)+(gvec(:,ig1)+kpoint(:))*ct
 end do

end subroutine matrixelements
!!***

!!****f* ABINIT/matrixelements_cutoff
!! NAME
!! matrixelements_cutoff
!!
!! FUNCTION
!!  TODO: To be better described. [Ask Carlo Motta]
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!      cutoff_m_elem
!!
!! CHILDREN
!!
!! SOURCE
subroutine matrixelements_cutoff(npwwfn,wfg1,wfg2,gvec,kpoint,z0,wdth,direction,res)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'matrixelements_cutoff'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
 integer,intent(in) :: npwwfn, direction
 integer,intent(in) :: gvec(3,npwwfn)
 complex(gwpc),intent(in) :: wfg1(npwwfn),wfg2(npwwfn)
 complex(gwpc),intent(out) :: res(3)
 real(dp),intent(in) :: kpoint(3)
 real(dp), intent(in) :: z0, wdth
!Local variables-------------------------------
 integer :: ig1,ig2,a,b
 complex(gwpc) :: ct, uimag
 integer :: gz(3)
 complex(gwpc) :: Fcut

! *************************************************************************

 uimag=dcmplx(0.d0,1.d0)
 res(:)=(0.0_gwp,0.0_gwp)
  a=1+mod(3+direction-2,3)
  b=1+mod(3+direction-3,3)
 do ig1=1,npwwfn
  do ig2=1,npwwfn
   if(gvec(a,ig1)/=gvec(a,ig2) .OR. gvec(b,ig1)/=gvec(b,ig2)) CYCLE 
   gz(:)=gvec(:,ig1)-gvec(:,ig2)
   if(gz(1)==0) then
   Fcut=cmplx(wdth,0.0)
   else
     Fcut= -uimag/(cmplx(2*pi*gz(direction),0.0))*CEXP(cmplx(0.0,2*pi*gz(direction)*z0))*&
           & (CEXP(cmplx(0.0,2*pi*gz(direction)*wdth))-1.0)
   endif
   ct=CONJG(wfg1(ig1))*wfg2(ig2)
   res(:)=res(:)+cmplx((gvec(:,ig1)+kpoint(:)),0.0)*ct*Fcut
  end do
 end do

end subroutine matrixelements_cutoff
!!***

