!{\src2tex{textfont=tt}}
!!****f* ABINIT/acfd_intexact
!! NAME
!! acfd_intexact
!!
!! FUNCTION
!! In the ACFD framework, the integration over the coupling constant can be
!! performed analytically in the RPA and in the PGG approximation (for spin-
!! compensated, two-electron systems). This yields :
!!
!!   $E_c=\int_0^\infty\frac{du}{2\pi\alpha} \\
!!          Tr\{ln[1-\alpha\chi_0(iu)v]+\alpha\chi_0(iu)v\}$
!!
!! where Tr denotes the trace, $\chi_0(iu)$ is the imaginary-frequency
!! Kohn-Sham susceptibility matrix, $v$ is the Coulomb interaction matrix,
!! $\alpha=1$ in the RPA, and $\alpha=0.5$ in the PGG approximation.
!! This subroutine returns :
!!
!!  $\delta E_c=\frac{1}{\alpha}Tr\{ln[1-\alpha\chi_0(iu)v]+\alpha\chi_0(iu)v\}$
!!
!! given the Kohn-Sham susceptibility matrix $\chi_0(iu)$, and using a cut-
!! off Coulomb interaction.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, MF, XG, GMR, LSI, YMN).
!! This file is distributed under the terms of the
!! GNU General Public License,see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  freq = (imaginary) frequency of the calculation (only used for output purposes).
!!  gsq(npwdiel) = the squared norm of the planewaves.
!!  ikhxc = option for the TDDFT kernel (see defs_parameters.f).
!!  mband = maximum number of bands for each k point and spin polarization.
!!  nband(nkpt*nsppol) = number of bands for each k point and spin polarization.
!!  nkpt = number of k points.
!!  npwdiel = number of planewaves for the susceptibility matrix.
!!  nspden = number of spin-density components.
!!  nsppol = number of spin polarizations.
!!  occ(mband*nkpt*nsppol) = occupation numbers for each band at each k point.
!!  occopt = option for occupancies (needed to check input parameters).
!!  rcut_coulomb = real space cut-off radius for the Coulomb interaction in Bohr.
!!  susmat(2,npwdiel,nspden,npwdiel,nspden) = susceptibility matrix.
!!
!! OUTPUT
!!  dec = the contribution to the correlation energy.
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!! Current restrictions are:
!!  a - Spin-polarized case not tested.
!!  b - Does not work in the metallic case.
!!
!! PARENTS
!!      xcacfd
!!
!! CHILDREN
!!      k_rpa,leave_new,timab,wrtout,zgetrf,zherk,zpotrf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine acfd_intexact(dec,freq,gsq,ikhxc,mband,nband,nkpt,npwdiel,&
&                         nspden,nsppol,occ,occopt,rcut_coulomb,susmat)

 use m_profiling

 use defs_basis
 use defs_parameters

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'acfd_intexact'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_77_suscep, except_this_one => acfd_intexact
!End of the abilint section

 implicit none

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: ikhxc,mband,nkpt,npwdiel,nspden,nsppol,occopt
 real(dp),intent(in) :: freq,rcut_coulomb
 real(dp),intent(out) :: dec
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: gsq(npwdiel),occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)

!Local variables -------------------------------------------------------
 character(len = *), parameter :: fmtd = '(a,t13,5(1x,es12.5))'
 character(len = *), parameter :: fmth1 = '(12x,3(1x,i12))'
 character(len = *), parameter :: fmth2 = '(12x,1x,a12,3(1x,i12))'
 character(len = *), parameter :: fmth3 = '(12x,1x,a12,1x,a12,3(1x,i12))'
!scalars
 integer :: iband,info,ipw,ipw1,ipw2,isp
 real(dp) :: alpha,krpajj,lndetsusmatk1,lndetsusmatk2
 logical :: pggok
 character(len=500) :: message
!arrays
 integer,allocatable :: ipiv(:)
 real(dp) :: detsusmatk(2),olddet(2),susmatij(2),trsusmatk(2),tsec(2)
 real(dp),allocatable :: diag(:,:),krpa(:),susmatk(:,:,:),susmatk2(:,:,:)

!***********************************************************************

 call timab(96,1,tsec) !Use timer for dieltcel.

!Check input parameters.

 if (nspden > 2) then
   write (message,'(4a)') ch10,&
&   ' acfd_intexact: ERROR - ',ch10,&
&   '  acfd_intexact does not work yet for nspden > 2.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if ((nspden == 2).and.((occopt >= 3).and.(occopt <= 8))) then
   write (message,'(4a)') ch10,&
&   ' acfd_intexact: ERROR - ',ch10,&
&   '  acfd_intexact does not work yet in the metallic, spin-polarized case.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 select case (ikhxc)
   case (ikhxc_RPA)
     alpha = 1.0_dp
   case (ikhxc_PGG)
     alpha = 0.5_dp
     case default
     write (message,'(6a)') ch10,&
&     ' acfd_intexact: ERROR - ',ch10,&
&     '  The exact integration over the coupling constant can only be',ch10,&
&     '  performed for the RPA and spin-compensated two-electron PGG kernels.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
 end select

!PGG: Check that we are in a spin-compensated, two-electron system.

 if (ikhxc == ikhxc_PGG) then

   pggok = .false.
   if ((nspden == 1).and.(nkpt == 1)) then
     pggok = (abs(occ(1)-2._dp) < tol12)
     do iband = 2,nband(1)
       pggok = pggok.and.(abs(occ(iband)) < tol12)
     end do
   end if

   if (.not.pggok) then
     write (message,'(6a)') ch10,&
&     ' acfd_intexact: ERROR - ',ch10,&
&     '  The exact integration over the coupling constant can be performed ',ch10,&
&     '  for the PGG kernel only in spin-compensated two-electron systems.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

 end if

 write (message,'(2a)') ch10,&
& ' acfd_intexact: Performing the exact integration over the coupling constant...'
 call wrtout(std_out,message,'COLL')

!Allocate memory.

 ABI_ALLOCATE(krpa,(npwdiel))
 ABI_ALLOCATE(diag,(2,npwdiel))
 ABI_ALLOCATE(susmatk,(2,npwdiel,npwdiel))

!Compute the Hartree kernel with a cut-off in real space.

 call k_rpa(gsq,krpa,npwdiel,1,rcut_coulomb)

 krpa(:) = alpha*krpa(:)

!Calculate the product of the Kohn-Sham susceptibility matrix and of the Hartree kernel.
!Here we assume that the Kohn-Sham susceptibility matrix susmat is spin-diagonal.

 do ipw2 = 1,npwdiel

   krpajj = krpa(ipw2)

   do ipw1 = 1,npwdiel

     susmatij(:) = 0._dp

     do isp = 1,nspden
       susmatij(:) = susmatij(:)+susmat(:,ipw1,isp,ipw2,isp)
     end do

     susmatk(:,ipw1,ipw2) = -susmatij(:)*krpajj

   end do

 end do

!Calculate the trace of susmatk and replace susmatk with 1+susmatk.

 trsusmatk(:) = 0._dp

 do ipw = 1,npwdiel

   trsusmatk(:) = trsusmatk(:)-susmatk(:,ipw,ipw)

   susmatk(1,ipw,ipw) = 1._dp+susmatk(1,ipw,ipw)

 end do

!DEBUG
!write (std_out,*) ' acfd_intexact: trsusmatk = ',trsusmatk(1),'+i',trsusmatk(2)
!ENDDEBUG

!Calculate the log of the determinant of susmatk.

 if (1 == 1) then

!  Compute the LU decompostion of susmatk using ZGETRF...

!  DEBUG
!  write (std_out,*) ' acfd_intexact: call to ZGETRF to compute the LU decomposition of susmatk'
!  ENDDEBUG

   ABI_ALLOCATE(ipiv,(npwdiel))

   call ZGETRF(npwdiel,npwdiel,susmatk,npwdiel,ipiv,info)
   if (info < 0) then
     write (message,'(4a)') ch10,&
&     ' acfd_intexact: BUG - ',ch10,&
&     '  Failed to compute the LU decomposition of 1-chi_0*v.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   ABI_DEALLOCATE(ipiv)

   do ipw = 1,npwdiel
     diag(:,ipw) = susmatk(:,ipw,ipw)
   end do

 else

!  First compute susmatk*susmatk**h using ZHERK...

   ABI_ALLOCATE(susmatk2,(2,npwdiel,npwdiel))

   call ZHERK('L','N',npwdiel,npwdiel,1._dp,susmatk,npwdiel,0._dp,susmatk2,npwdiel)

!  then its Cholesky decomposition using ZPOTRF...

!  DEBUG
!  write (std_out,*) ' acfd_intexact: call to ZPOTRF to compute the Cholesky decomposition of susmatk2'
!  ENDDEBUG

   call ZPOTRF('L',npwdiel,susmatk2,npwdiel,info)
   if (info < 0) then
     write (message,'(4a)') ch10,&
&     ' acfd_intexact: BUG - ',ch10,&
&     '  Failed to compute the Cholesky decomposition of (1-chi_0*v)*(1-chi_0*v)**h.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   do ipw = 1,npwdiel
     diag(:,ipw) = susmatk2(:,ipw,ipw)
   end do

   ABI_DEALLOCATE(susmatk2)

 end if

!DEBUG
!write (std_out,*) ' acfd_intexact: diag ='
!do ipw = 1,npwdiel
!do ipw = 1,min(npwdiel,100)
!write (std_out,*) ' ',ipw,' = ',diag(1,ipw),'+i',diag(2,ipw)
!end do
!if (npwdiel > 100) write (std_out,*) ' [...]'
!ENDDEBUG

!Calculate the log of the determinant.

 detsusmatk = (/1._dp,0._dp/)

 lndetsusmatk1 = 0._dp

 do ipw = 1,npwdiel

   olddet(:) = detsusmatk(:)
   detsusmatk(1) = olddet(1)*diag(1,ipw)-olddet(2)*diag(2,ipw)
   detsusmatk(2) = olddet(1)*diag(2,ipw)+olddet(2)*diag(1,ipw)

   lndetsusmatk1 = lndetsusmatk1+log(diag(1,ipw))

 end do

 lndetsusmatk2 = log(detsusmatk(1))

!DEBUG
!write (std_out,*) ' acfd_intexact: det(susmatk)  = ',detsusmatk(1),'+i',detsusmatk(2)
!write (std_out,*) ' acfd_intexact: ln(det(susmatk)) = ',lndetsusmatk2
!write (std_out,*) '  evaluated from det(susmatk) above'
!write (std_out,*) ' acfd_intexact: ln(det(susmatk)) = ',lndetsusmatk1
!write (std_out,*) '  evaluated from sum_i ln(LU_susmatk_ii)'
!ENDDEBUG

!Calculate the contribution to the correlation energy.

 dec = (lndetsusmatk1+trsusmatk(1))/(alpha*two_pi)

 write (message,fmth2) 'frequency'
 call wrtout(std_out,message,'COLL')

 write (message,fmtd) &
& ' dEc_cut',freq,dec
 call wrtout(std_out,message,'COLL')

!Check the consistency of the result.

 if (abs(trsusmatk(2)) > tol12) then
   write (message,'(4a,es12.5,a)') ch10,&
&   ' acfd_intexact: WARNING - ',ch10,&
&   '  Im(Tr(chi_0*v)) = ',trsusmatk(2),', should be zero (up to 1.d-12).'
   call wrtout(std_out,message,'COLL')
 end if

 if (abs(detsusmatk(2)) > tol12) then
   write (message,'(4a,es12.5,a)') ch10,&
&   ' acfd_intexact: WARNING - ',ch10,&
&   '  Im(det(1-chi_0*v)) = ',detsusmatk(2),', should be zero (up to 1.d-12).'
   call wrtout(std_out,message,'COLL')
 end if

 if (abs(lndetsusmatk2-lndetsusmatk1) > tol10) then
   write (message,'(4a,es20.13,3a,es20.13,3a)') ch10,&
&   ' acfd_intexact: WARNING - ',ch10,&
&   '  ln(det(1-chi_0*v)) = ',lndetsusmatk1,' from det(1-chi_0*v).',ch10,&
&   '  ln(det(1-chi_0*v)) = ',lndetsusmatk2,' from sum_i ln(LU_ii).',ch10,&
&   '  Both evaluations differ by more than 1.d-10.'
   call wrtout(std_out,message,'COLL')
 end if

!Free memory.

 ABI_DEALLOCATE(krpa)
 ABI_DEALLOCATE(diag)
 ABI_DEALLOCATE(susmatk)

 call timab(96,2,tsec)

end subroutine acfd_intexact

!!***
