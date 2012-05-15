!{\src2tex{textfont=tt}}
!!****f* ABINIT/nlinopt
!! NAME
!! nlinopt
!!
!! FUNCTION
!! This routine compute optical frequency dependent second harmonic generation
!! suscptibility for semiconductors
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (SSharma)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nspin = number of spins(integer)
!!  omega = crystal volume in au (real)
!!  nkpt  = total number of kpoints (integer)
!!  wkpt(nkpt) = weights of kpoints (real)
!!  nsymcrys = number of crystal symmetry operations(integer)
!!  symcrys(3,3,nsymcrys) = symmetry operations in cartisian coordinates(real)
!!  nstval = total number of valence states(integer)
!!  evalv(nstval,nspin,nkpt) = eigen value for each band in Ha(real)
!!  efermi = Fermi energy in Ha(real)
!!  pmat(nstval,nstval,nkpt,3,nspin) = momentum matrix elements in cartesian coordinates(complex)
!!                                     : changes on exit
!!  v1,v2,v3 = desired component of the dielectric function(integer) 1=x,2=y,3=z
!!  nmesh = desired number of energy mesh points(integer)
!!  de = desired step in energy(real); nmesh*de=maximum energy for plotting
!!  sc = scissors shift in Ha(real)
!!  brod = broadening in Ha(real)
!!  tol = tolerance:how close to the singularity exact exact is calculated(real)
!!  fnam=root for filenames that will contain the output  :
!!   fnam1=trim(fnam)//'-ChiTotIm.out'
!!   fnam2=trim(fnam)//'-ChiTotRe.out'
!!   fnam3=trim(fnam)//'-ChiIm.out'
!!   fnam4=trim(fnam)//'-ChiRe.out'
!!   fnam5=trim(fnam)//'-ChiAbs.out'
!!
!! OUTPUT
!!  Calculates the second harmonic generation susceptibility on a desired energy mesh and
!!  for desired direction of polarisation. The output is in files named
!!  ChiTot.out : Im\chi_{v1v2v3}(2\omega,\omega,-\omega) and Re\chi_{v1v2v3}(2\omega,\omega,-\omega)
!!  ChiIm.out  : contributions to the Im\chi_{v1v2v3}(2\omega,\omega,-\omega) from various terms
!!  ChiRe.out  : contributions to Re\chi_{v1v2v3}(2\omega,\omega,-\omega) from various terms
!!  ChiAbs.out : abs\chi_{v1v2v3}(2\omega,\omega,-\omega). The headers in these files contain
!!  information about the calculation.
!!
!! SIDE EFFECTS
!!  pmat(nstval,nstval,nkpt,3,nspin) = momentum matrix elements in cartesian coordinates(complex)
!!
!! COMMENTS
!!  Right now the routine sums over the k-points. In future linear tetrahedron method might be
!!  useful.
!!  Reference articles:
!!  1. S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl, Phys. Rev. B {\bf 67} 165332 2003
!!  2. J. L. P. Hughes and J. E. Sipe, Phys. Rev. B {\bf 53} 10 751 1996
!!  3. S. Sharma and C. Ambrosch-Draxl, Physica Scripta T 109 2004
!!  4. J. E. Sipe and Ed. Ghahramani, Phys. Rev. B {\bf 48} 11 705 1993
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine nlinopt(nspin,omega,nkpt,wkpt,nsymcrys,symcrys,nstval,evalv,efermi, &
  pmat,v1,v2,v3,nmesh,de,sc,brod,tol,fnam)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nlinopt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!no_abirules
integer, intent(in) :: nspin
real(dp), intent(in) :: omega
integer, intent(in) :: nkpt
real(dp), intent(in) :: wkpt(nkpt)
integer, intent(in) :: nsymcrys
real(dp), intent(in) :: symcrys(3,3,nsymcrys)
integer, intent(in) :: nstval
real(dp), intent(in) :: evalv(nstval,nspin,nkpt)
real(dp), intent(in) :: efermi
complex(dpc), intent(inout) :: pmat(nstval,nstval,nkpt,3,nspin)
integer, intent(in) :: v1
integer, intent(in) :: v2
integer, intent(in) :: v3
integer, intent(in) :: nmesh
real(dp), intent(in) :: de
real(dp), intent(in) :: sc
real(dp), intent(in) :: brod
real(dp), intent(in) :: tol
character(256), intent(in) :: fnam

!Local variables -------------------------
!no_abirules
! present calculation related (user specific)
integer :: iw
integer :: i1,i2,i3,i,j,k,lx,ly,lz
integer :: isp,isym,ik
integer :: ist1,ist2,istl,istn,istm
real(dp) :: ha2ev
real(dp) :: t1,t2,t3,tst
real(dp) :: ene,totre,totabs,totim
real(dp) :: e1,e2,el,en,em
real(dp) :: emin,emax
real(dp) :: const_esu,const_au,au2esu
real(dp) :: wmn,wnm,wln,wnl,wml,wlm
complex(dpc) :: idel,w,zi
complex(dpc) :: mat2w,mat1w1,mat1w2,mat2w_tra,mat1w3_tra
complex(dpc) :: b111,b121,b131,b112,b122,b132,b113,b123,b133
complex(dpc) :: b241,b242,b243,b221,b222,b223,b211,b212,b213,b231
complex(dpc) :: b311,b312,b313,b331
complex(dpc) :: b24,b21_22,b11,b12_13,b31_32
character(256) :: fnam1,fnam2,fnam3,fnam4,fnam5
! local allocatable arrays
real(dp), allocatable :: s(:,:)
real(dp), allocatable :: sym(:,:,:)
complex(dpc), allocatable :: px(:,:,:,:,:)
complex(dpc), allocatable :: py(:,:,:,:,:)
complex(dpc), allocatable :: pz(:,:,:,:,:)
complex(dpc), allocatable :: delta(:,:,:)
complex(dpc), allocatable :: inter2w(:)
complex(dpc), allocatable :: inter1w(:)
complex(dpc), allocatable :: intra2w(:)
complex(dpc), allocatable :: intra1w(:)

! *********************************************************************

!calculate the constant
 zi=(0._dp,1._dp)
 idel=zi*brod
 const_au=-2._dp/(omega*dble(nsymcrys))
 au2esu=5.8300348177d-8
 const_esu=const_au*au2esu
 ha2ev=13.60569172*2._dp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!5.8300348177d-8 : au2esu : bohr*c*10^4/4pi*2*ry2ev
!bohr: 5.2917ifc nlinopt.f907E-11
!c: 2.99792458   velocity of sound
!ry2ev: 13.60569172
!au2esu=(5.29177E-11*2.99792458*1.0E4)/(13.60569172*2)
!this const includes (e^3*hbar^3*hbar^3)/(vol*hbar^5*m_e^3)
!mass comes from converting P_mn to r_mn
!hbar^3 comes from converting all frequencies to energies in denominator
!hbar^3 comes from operator for momentum (hbar/i nabla)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!output file names
 fnam1=trim(fnam)//'-ChiTotIm.out'
 fnam2=trim(fnam)//'-ChiTotRe.out'
 fnam3=trim(fnam)//'-ChiIm.out'
 fnam4=trim(fnam)//'-ChiRe.out'
 fnam5=trim(fnam)//'-ChiAbs.out'
!fool proof:
!If there exists inversion symmetry exit with a mesg.
 tst=1.d-09
 do isym=1,nsymcrys
   t1=symcrys(1,1,isym)+1
   t2=symcrys(2,2,isym)+1
   t3=symcrys(3,3,isym)+1
!  test if diagonal elements are -1
   if (abs(t1).lt.tst.and.abs(t2).lt.tst.and.abs(t3).lt.tst) then
!    test if off-diagonal elements are zero
     if (abs(symcrys(1,2,isym)).lt.tst.and.abs(symcrys(1,3,isym)).lt.tst &
     .and.abs(symcrys(2,1,isym)).lt.tst.and.abs(symcrys(2,3,isym)).lt.tst.and.  &
     abs(symcrys(3,1,isym)).lt.tst.and.abs(symcrys(3,2,isym)).lt.tst) then
       write(std_out,*) '-----------------------------------------'
       write(std_out,*) '    the crystal has inversion symmetry   '
       write(std_out,*) '    the SHG susceptibility is zero       '
       write(std_out,*) '-----------------------------------------'
       stop
     end if
   end if
 end do
!check polarisation
 if (v1.le.0.or.v2.le.0.or.v3.le.0.or.v1.gt.3.or.v2.gt.3.or.v3.gt.3) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in nlinopt:                        '
   write(std_out,*) '    the polarisation directions incorrect    '
   write(std_out,*) '    1=x,  2=y  and 3=z                       '
   write(std_out,*) '---------------------------------------------'
   stop
 end if
!number of energy mesh points
 if (nmesh.le.0) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in nlinopt:                        '
   write(std_out,*) '    number of energy mesh points incorrect   '
   write(std_out,*) '    number has to be integer greater than 0  '
   write(std_out,*) '    nmesh*de = max energy for calculation    '
   write(std_out,*) '---------------------------------------------'
   stop
 end if
!step in energy
 if (de.le.0._dp) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in nlinopt:                        '
   write(std_out,*) '    energy step is incorrect                 '
   write(std_out,*) '    number has to real greater than 0.0      '
   write(std_out,*) '    nmesh*de = max energy for calculation    '
   write(std_out,*) '---------------------------------------------'
   stop
 end if
!broadening
 if (brod.gt.0.009) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    ATTENTION: broadening is quite high      '
   write(std_out,*) '    ideally should be less than 0.005        '
   write(std_out,*) '---------------------------------------------'
 else if (brod.gt.0.015) then
   write(std_out,*) '----------------------------------------'
   write(std_out,*) '    ATTENTION: broadening is too high   '
   write(std_out,*) '    ideally should be less than 0.005   '
   write(std_out,*) '----------------------------------------'
 end if
!tolerance
 if (tol.gt.0.006) then
   write(std_out,*) '----------------------------------------'
   write(std_out,*) '    ATTENTION: tolerance is too high    '
   write(std_out,*) '    ideally should be less than 0.004   '
   write(std_out,*) '----------------------------------------'
 end if

!fool proof ends
!
!allocate local arrays
 ABI_ALLOCATE(px,(nstval,nstval,3,3,3))
 ABI_ALLOCATE(py,(nstval,nstval,3,3,3))
 ABI_ALLOCATE(pz,(nstval,nstval,3,3,3))
 ABI_ALLOCATE(inter2w,(nmesh))
 ABI_ALLOCATE(inter1w,(nmesh))
 ABI_ALLOCATE(intra2w,(nmesh))
 ABI_ALLOCATE(intra1w,(nmesh))
 ABI_ALLOCATE(delta,(nstval,nstval,3))
 ABI_ALLOCATE(sym,(3,3,3))
 ABI_ALLOCATE(s,(3,3))
!generate the symmetrizing tensor
!TODO: check if this needs cartesian instead of reduced coordinate symops
 sym(:,:,:)=0._dp
 do isym=1,nsymcrys
   s(:,:)=symcrys(:,:,isym)
   do i=1,3
     do j=1,3
       do k=1,3
         sym(i,j,k)=sym(i,j,k)+(s(i,v1)*s(j,v2)*s(k,v3))
       end do
     end do
   end do
 end do
!initialise
 inter2w(:)=0._dp
 inter1w(:)=0._dp
 intra2w(:)=0._dp
 intra1w(:)=0._dp
 delta(:,:,:)=0._dp
 emin=0._dp
 emax=0._dp
!start loop over kpts
 do ik=1,nkpt
   write(std_out,*) ik,'of',nkpt
!  start loop over spins
   do isp=1,nspin
!    start loop over states
     do ist1=1,nstval
       e1=evalv(ist1,isp,ik)
       if (e1.lt.efermi) then   ! ist1 is a valence state
         do ist2=1,nstval
           e2=evalv(ist2,isp,ik)
           if (e2.gt.efermi) then ! ist2 is a conduction state
!            symmetrize the momentum matrix elements
             do lx=1,3
               do ly=1,3
                 do lz=1,3
                   i1=sym(lx,ly,lz)+sym(lx,lz,ly)
                   i2=sym(ly,lx,lz)+sym(ly,lz,lx)
                   i3=sym(lz,lx,ly)+sym(lz,ly,lx)
                   px(ist1,ist2,lx,ly,lz)=i1*pmat(ist1,ist2,ik,lx,isp)
                   py(ist2,ist1,lx,ly,lz)=i2*pmat(ist2,ist1,ik,lx,isp)   ! TODO: check if this should not be ly, lz on these lines
                   pz(ist2,ist1,lx,ly,lz)=i3*pmat(ist2,ist1,ik,lx,isp)
                 end do
               end do
             end do
!            end loop over states
           end if
         end do
       end if
     end do
!    calculate the energy window and \Delta_nm
     do ist1=1,nstval
       emin=min(emin,evalv(ist1,isp,ik))
       emax=max(emax,evalv(ist1,isp,ik))
       do ist2=1,nstval
         delta(ist1,ist2,1:3)=pmat(ist1,ist1,ik,1:3,isp)-pmat(ist2,ist2,ik,1:3,isp)
       end do
     end do
!    initialise the factors
!    factors are named according to the Ref. article 2.
     b111=0._dp
     b121=0._dp
     b131=0._dp
     b112=0._dp
     b122=0._dp
     b132=0._dp
     b113=0._dp
     b123=0._dp
     b133=0._dp
     b211=0._dp
     b221=0._dp
     b212=0._dp
     b222=0._dp
     b213=0._dp
     b223=0._dp
     b231=0._dp
     b241=0._dp
     b242=0._dp
     b243=0._dp
     b311=0._dp
     b312=0._dp
     b313=0._dp
     b331=0._dp
!    start the calculation
     do istn=1,nstval
       en=evalv(istn,isp,ik)
       if (en.lt.efermi) then    ! istn is a valence state
         do istm=1,nstval
           em=evalv(istm,isp,ik)
           if (em.gt.efermi) then   ! istm is a conduction state
             wmn=em+sc-en
             wnm=-wmn
!            calculate the matrix elements for two band intraband term
             mat2w_tra=0._dp
             mat1w3_tra=0._dp
             do lx=1,3
               do ly=1,3
                 do lz=1,3
                   mat2w_tra=mat2w_tra+px(istn,istm,lx,ly,lz)*pmat(istm,istn,ik,lz,isp)    &
                   *delta(istm,istn,ly)
                   mat1w3_tra=mat1w3_tra+px(istn,istm,lx,ly,lz)*pmat(istm,istn,ik,ly,isp)  &
                   *delta(istm,istn,lz)
!                  NOTE:: lx to ly m to n in pmat matrices respectively
!                  Changes are made so that this (b3) term is according to paper
!                  PRB48(Ref. 4) rather than PRB53(Ref 2) in which this term is incorrect
                 end do
               end do
             end do
             b331=mat1w3_tra/wnm
             b11=0._dp
             b12_13=0._dp
             b24=0._dp
             b31_32=0._dp
             b21_22=0._dp

             b231=8._dp*mat2w_tra/wmn
             b331=mat1w3_tra/(wnm)
!            !!!!!!!!!!!!!!!!!!!
!            istl < istn   !
!            !!!!!!!!!!!!!!!!!!!
             do istl=1,istn-1           ! istl is a valence state below istn
               el=evalv(istl,isp,ik)
               wln=el-en                ! do not add sc to the valence bands!
               wml=em+sc-el
               wnl=-wln
               wlm=-wml
!              calculate the matrix elements for three band terms
               mat2w=0._dp
               mat1w1=0._dp
               mat1w2=0._dp
               do lx=1,3
                 do ly=1,3
                   do lz=1,3

                     mat2w=mat2w+(px(istn,istm,lx,ly,lz)*pmat(istm,istl,ik,ly,isp)   &
                     *pmat(istl,istn,ik,lz,isp))

                     mat1w1=mat1w1+(py(istm,istn,lx,ly,lz)*pmat(istl,istm,ik,lz,isp) &
                     *pmat(istn,istl,ik,ly,isp))

                     mat1w2=mat1w2+(pz(istm,istn,lx,ly,lz)*pmat(istl,istm,ik,lz,isp) &
                     *pmat(istn,istl,ik,ly,isp))
                   end do
                 end do
               end do
               b111=mat2w*(1._dp/(wln+wlm))*(1._dp/wlm)
               b121=mat1w1*(1._dp/(wnm+wlm))*(1._dp/wlm)
               b131=mat1w2*(1._dp/wlm)
!              
               b221=0._dp
               b211=mat1w1/wml
               b241=-mat2w/wml
!              
               b311=mat1w2/wlm
               if (abs(wln).gt.tol) then
                 b111=b111/wln
                 b121=b121/wln
                 b131=b131/wln
                 b221=mat1w2/wln
                 b241=b241+(mat2w/wln)
                 b311=b311+(mat1w1/wln)
               else
                 b111=0._dp
                 b121=0._dp
                 b131=0._dp
                 b221=0._dp
               end if
               t1=wln-wnm
               if (abs(t1).gt.tol) then
                 b131=b131/t1
               else
                 b131=0._dp
               end if
               b11=b11-2._dp*b111
               b12_13=b12_13+b121+b131
               b21_22=-b211+b221
               b24=b24+2._dp*b241
               b31_32=b31_32+b311
!              end loop over istl
             end do
!            !!!!!!!!!!!!!!!!!!!!!!!!!!!
!            istn < istl < istm    !
!            !!!!!!!!!!!!!!!!!!!!!!!!!!!
             do istl=istn+1,istm-1
               el=evalv(istl,isp,ik)
!              calculate the matrix elements for three band terms
               mat2w=0._dp
               mat1w1=0._dp
               mat1w2=0._dp
               do lx=1,3
                 do ly=1,3
                   do lz=1,3

                     mat2w=mat2w+(px(istn,istm,lx,ly,lz)*pmat(istm,istl,ik,ly,isp)   &
                     *pmat(istl,istn,ik,lz,isp))

                     mat1w1=mat1w1+(py(istm,istn,lx,ly,lz)*pmat(istl,istm,ik,lz,isp) &
                     *pmat(istn,istl,ik,ly,isp))

                     mat1w2=mat1w2+(pz(istm,istn,lx,ly,lz)*pmat(istl,istm,ik,lz,isp) &
                     *pmat(istn,istl,ik,ly,isp))
                   end do
                 end do
               end do
               if (el.lt.efermi) then
                 wln=el-en
                 wnl=-wln
                 wml=em-el
                 wlm=-wml
               else
                 el=el+sc
                 wln=el-en
                 wnl=-wln
                 wml=em-el
                 wlm=-wml
               end if
!              
               b112=0._dp
               b122=mat1w1*(1._dp/(wnm+wlm))
               b132=mat1w2*(1._dp/(wnm+wnl))
               b242=0._dp
               b222=0._dp
               b212=0._dp
               if (abs(wnl).gt.tol) then
                 b112=mat2w/wln
                 b122=b122/wnl
                 b132=b132/wnl
                 b242=mat2w/wln
                 b222=mat1w2/wln
                 b312=mat1w1/wln
               else
                 b122=0._dp
                 b132=0._dp
               end if
               if (abs(wlm).gt.tol) then
                 b112=b112/wml
                 b122=b122/wlm
                 b132=b132/wlm
                 b242=b242-(mat2w/wml)
                 b212=mat1w1/wml
                 b312=b312+(mat1w2/wlm)
               else
                 b112=0._dp
                 b122=0._dp
                 b132=0._dp
                 b212=0._dp
               end if
               t1=wlm-wnl
               if (abs(t1).gt.tol) then
                 b112=b112/t1
               else
                 b112=0._dp
               end if
               b11=b11+2._dp*b112
               b12_13=b12_13-b122+b132
               b24=b24+2._dp*b242
               b21_22=b21_22-b212+b222
               b31_32=b31_32+b312
!              end loop over istl
             end do

!            !!!!!!!!!!!!!!!!!!!!!
!            istl > istm    !
!            !!!!!!!!!!!!!!!!!!!!!
             do istl=istm+1,nstval
               el=evalv(istl,isp,ik)+sc
               wln=el-en
               wnl=-wln
               wml=em-el
               wlm=-wml
!              calculate the matrix elements for three band terms
               mat2w=0._dp
               mat1w1=0._dp
               mat1w2=0._dp
               do lx=1,3
                 do ly=1,3
                   do lz=1,3
                     mat2w=mat2w+px(istn,istm,lx,ly,lz)*pmat(istm,istl,ik,ly,isp) &
                     *pmat(istl,istn,ik,lz,isp)

                     mat1w1=mat1w1+(py(istm,istn,lx,ly,lz)*pmat(istl,istm,ik,lz,isp) &
                     *pmat(istn,istl,ik,ly,isp))

                     mat1w2=mat1w2+(pz(istm,istn,lx,ly,lz)*pmat(istl,istm,ik,lz,isp) &
                     *pmat(istn,istl,ik,ly,isp))
                   end do
                 end do
               end do
!              
               b113=mat2w*(1._dp/(wnl+wml))*(1._dp/wnl)
               b123=mat1w1*(1._dp/wnl)
               b133=mat1w2*(1._dp/wnl)*(1._dp/(wnl+wnm))
               b243=mat2w/wln
               b223=mat1w2/wln
               b213=0._dp
               b313=-1._dp*mat1w1/wnl
               if (abs(wml).gt.tol) then
                 b113=b113/wml
                 b123=b123/wml
                 b133=b133/wml
                 b243=b243-(mat2w/wml)
                 b213=mat1w1/wml
                 b313=b313+(mat1w2/wlm)
               else
                 b113=0._dp
                 b123=0._dp
                 b133=0._dp
               end if
               t1=wnm-wml
               if (abs(t1).gt.tol) then
                 b123=b123/t1
               else
                 b123=0._dp
               end if
               b11=b11+2._dp*b113
               b12_13=b12_13+b123-b133
               b21_22=b21_22-b213+b223
               b24=b24+2._dp*b243
               b31_32=b31_32+b313
!              end loop over istl
             end do
!            
             b11=b11*zi*(1._dp/wnm)*const_esu
             b12_13=b12_13*zi*(1._dp/wnm)*const_esu
             b24=(b24+b231)*zi*(1._dp/(wnm**3))*const_esu
             b21_22=(b21_22)*zi*(1._dp/(wnm**3))*const_esu
             b31_32=(b31_32-b331)*zi*(1._dp/(wmn**3))*const_esu*0.5_dp
!            calculate over the desired energy mesh and sum over k-points
             do iw=1,nmesh
               w=(iw-1)*de+idel
               inter2w(iw)=inter2w(iw)+(wkpt(ik)*(b11/(wmn-2._dp*w)))
               inter1w(iw)=inter1w(iw)+(wkpt(ik)*(b12_13/(wmn-w)))
               intra2w(iw)=intra2w(iw)+(wkpt(ik)*(b24/(wmn-2._dp*w)))
               intra1w(iw)=intra1w(iw)+(wkpt(ik)*((b21_22+b31_32)/(wmn-w)))
             end do
!            end loop over istn and istm
           end if
         end do
       end if
     end do
!    end loop over spins
   end do
!  end loop over k-points
 end do
!
!write output in SI units and esu (esu to SI(m/v)=(value_esu)*(4xpi)/30000)
!
 open(92,file=fnam1,action='WRITE',form='FORMATTED')
 open(93,file=fnam2,action='WRITE',form='FORMATTED')
 open(94,file=fnam3,action='WRITE',form='FORMATTED')
 open(95,file=fnam4,action='WRITE',form='FORMATTED')
 open(96,file=fnam5,action='WRITE',form='FORMATTED')
!write headers
 write(92, '(a,3i3)' ) ' #calculated the component:',v1,v2,v3
 write(92, '(a,es16.6)' ) ' #tolerence:',tol
 write(92, '(a,es16.6,a)' ) ' #broadening:',brod,'Ha'
 write(92, '(a,es16.6,a)' ) ' #scissors shift:',sc,'Ha'
 write(92, '(a,es16.6,a,es16.6,a)' ) ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
 write(92, '(a)' )' # Energy      Tot-Im Chi(-2w,w,w)  Tot-Im Chi(-2w,w,w)'
 write(92, '(a)' )' # eV          *10^-7 esu        *10^-12 m/V SI units '
 write(92, '(a)' )' # '

 write(93, '(a,3i3)' ) ' #calculated the component:',v1,v2,v3
 write(93, '(a,es16.6)') ' #tolerence:',tol
 write(93, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
 write(93, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
 write(93, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
 write(93, '(a)')' # Energy      Tot-Im Chi(-2w,w,w)  Tot-Im Chi(-2w,w,w)'
 write(93, '(a)')' # eV          *10^-7 esu        *10^-12 m/V SI units '
 write(93, '(a)')' # '

 write(94, '(a,3i3)') ' #calculated the component:',v1,v2,v3
 write(94, '(a,es16.6)') ' #tolerence:',tol
 write(94, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
 write(94, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
 write(94, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
 write(94, '(a)')' # Energy(eV) Inter(2w) inter(1w) intra(2w) intra(1w)'
 write(94, '(a)')' # in esu'
 write(94, '(a)')' # '

 write(95, '(a,3i3)') ' #calculated the component:',v1,v2,v3
 write(95, '(a,es16.6)') ' #tolerence:',tol
 write(95, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
 write(95, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
 write(95, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
 write(95, '(a)')' # Energy(eV) Inter(2w) inter(1w) intra(2w) intra(1w)'
 write(95, '(a)')' # in esu'
 write(95, '(a)')' # '

 write(96, '(a,3i3)') ' #calculated the component:',v1,v2,v3
 write(96, '(a,es16.6)') ' #tolerence:',tol
 write(96, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
 write(96, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
 write(96, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
 write(96, '(a)')' # Energy(eV)  |TotChi(-2w,w,w)|   |Tot Chi(-2w,w,w)|'
 write(96, '(a)')' # eV          *10^-7 esu        *10^-12 m/V SI units '
 write(96, '(a)')' # '
!
 totim=0._dp
 totre=0._dp
 totabs=0._dp
 do iw=2,nmesh
   ene=(iw-1)*de
   ene=ene*ha2ev
!  
   totim=aimag(inter2w(iw)+inter1w(iw)+intra2w(iw)+intra1w(iw))/1.d-7
   write(92,'(f15.6,2es15.6)') ene,totim,totim*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
   totim=0._dp
!  
   totre=dble(inter2w(iw)+inter1w(iw)+intra2w(iw)+intra1w(iw))/1.d-7
   write(93,'(f15.6,2es15.6)') ene,totre,totre*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
   totre=0._dp
!  
   write(94,'(f15.6,4es15.6)') ene,aimag(inter2w(iw))/1.d-7,      &
   aimag(inter1w(iw))/1.d-7,aimag(intra2w(iw))/1.d-7, aimag(intra1w(iw))/1.d-7
!  
   write(95,'(f15.6,4es15.6)') ene,dble(inter2w(iw))/1.d-7,       &
   dble(inter1w(iw))/1.d-7,dble(intra2w(iw))/1.d-7,dble(intra1w(iw))/1.d-7
!  
   totabs=abs(inter2w(iw)+inter1w(iw)+intra2w(iw)+intra1w(iw))/1.d-7
   write(96,'(f15.6,2es15.6)') ene,totabs,totabs*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
   totabs=0._dp
 end do
!deallocate local arrays
 ABI_DEALLOCATE(px)
 ABI_DEALLOCATE(py)
 ABI_DEALLOCATE(pz)
 ABI_DEALLOCATE(sym)
 ABI_DEALLOCATE(s)
 ABI_DEALLOCATE(inter2w)
 ABI_DEALLOCATE(inter1w)
 ABI_DEALLOCATE(intra2w)
 ABI_DEALLOCATE(intra1w)
 ABI_DEALLOCATE(delta)
 close(92)
 close(93)
 close(94)
 close(95)
 close(96)
!print information
 write(std_out,*) ' '
 write(std_out,*) 'information about calculation just performed:'
 write(std_out,*) ' '
 write(std_out,*) 'calculated the component:',v1,v2,v3 ,'of second order susceptibility'
 write(std_out,*) 'tolerence:',tol
 if (tol.gt.0.008) write(std_out,*) 'ATTENTION: tolerence is too high'
 write(std_out,*) 'broadening:',brod,'Hartree'
 if (brod.gt.0.009) then
   write(std_out,*) ' '
   write(std_out,*) 'ATTENTION: broadening is quite high'
   write(std_out,*) ' '
 else if (brod.gt.0.015) then
   write(std_out,*) ' '
   write(std_out,*) 'ATTENTION: broadening is too high'
   write(std_out,*) ' '
 end if
 write(std_out,*) 'scissors shift:',sc,'Hartree'
 write(std_out,*) 'energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Hartree'
!
 return
end subroutine nlinopt
!!***
