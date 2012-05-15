!{\src2tex{textfont=tt}}
!!****f* ABINIT/prteigrs
!!
!! NAME
!! prteigrs
!!
!! FUNCTION
!! Print out eigenvalues band by band and k point by k point.
!! If option=1, do it in a standard way, for self-consistent calculations.
!! If option=2, print out residuals and eigenvalues, in a format
!! adapted for nonself-consistent calculations, within the loops.
!! If option=3, print out eigenvalues, in a format
!! adapted for nonself-consistent calculations, at the end of the job.
!! If option=4, print out derivatives of eigenvalues (same format as option==3, except header that is printed)
!! If option=5, print out Fan contribution to zero-point motion correction to eigenvalues (averaged)
!!                  (same format as option==3, except header that is printed)
!! If option=6, print out DDW contribution to zero-point motion correction to eigenvalues (averaged)
!!                  (same format as option==3, except header that is printed)
!! If option=7, print out Fan+DDW contribution to zero-point motion correction to eigenvalues (averaged)
!!                  (same format as option==3, except header that is printed)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  eigen(mband*nkpt*nsppol)=eigenvalues (hartree)
!!   or, if option==4, diagonal of derivative of eigenvalues
!!   or, if option==5...7, zero-point motion correction to eigenvalues (averaged)
!!  enunit=choice parameter: 0=>output in hartree; 1=>output in eV;
!!   2=> output in both hartree and eV
!!  fermie=fermi energy (Hartree)
!!  fname_eig=filename for printing of the eigenenergies
!!  iout=unit number for formatted output file
!!  iscf=option for self-consistency
!!  kptns(3,nkpt)=k points in reduced coordinates
!!  kptopt=option for the generation of k points
!!  mband=maximum number of bands
!!  nband(nkpt)=number of bands at each k point
!!  nkpt=number of k points
!!  nnsclo_now=number of non-self-consistent loops for the current vtrial
!!    (often 1 for SCF calculation, =nstep for non-SCF calculations)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(maxval(nband(:))*nkpt*nsppol)=occupancies for each band and k point
!!  occopt=option for occupancies
!!  option= (see above)
!!  prteig=control print eigenenergies
!!  prtvol=control print volume and debugging
!!  resid(mband*nkpt*nsppol)=residuals (hartree**2)
!!  tolwfr=tolerance on band residual of wf, hartrees**2 (needed when option=2)
!!  vxcavg=average of vxc potential
!!  wtk(nkpt)=k-point weights
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      clnup1,loper3,respfn,scprqt,vtorho
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prteigrs(eigen,enunit,fermie,fname_eig,iout,iscf,kptns,kptopt,mband,nband,&
&  nkpt,nnsclo_now,nsppol,occ,occopt,option,prteig,prtvol,resid,tolwfr,vxcavg,wtk)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prteigrs'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enunit,iout,iscf,kptopt,mband,nkpt,nnsclo_now,nsppol
 integer,intent(in) :: occopt,option,prteig,prtvol
 real(dp),intent(in) :: fermie,tolwfr,vxcavg
 character(len=fnlen),intent(in) :: fname_eig
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),kptns(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),resid(mband*nkpt*nsppol)
 real(dp),intent(in) :: wtk(nkpt)

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt_max=50
 integer :: band_index,iband,ii,ikpt,isppol,nband_index,nband_k,nkpt_eff
 integer :: tmagnet,tmetal
 real(dp) :: convrt,magnet,residk,rhodn,rhoup
 character(len=39) :: kind_of_output
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,*)' prteigrs : enter, prteig,fname_eig=',prteig,',.',trim(fname_eig),'.'
!ENDDEBUG

 if (nsppol<1.or.nsppol>2) then
   write(message, '(a,a,a,a,i12,a)' ) ch10,&
&   ' prteigrs: BUG -',ch10,&
&   '  nsppol must be 1 or 2. Argument was ',nsppol,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if (prteig > 0) then
   write(message, '(a,a)' ) ' prteigrs : about to open file ',fname_eig
   call wrtout(iout,message,'COLL')
   open (unit=tmp_unit,file=fname_eig,status='unknown',form='formatted')
!  always rewind disk file and print latest eigenvalues
   rewind(tmp_unit)
 end if

 kind_of_output=              ' Eigenvalues                          '
 if(option==4) kind_of_output=' Expectation of eigenvalue derivatives'
 if(option==5) kind_of_output=' Fan corrections to eigenvalues at T=0'
 if(option==6) kind_of_output=' DDW corrections to eigenvalues at T=0'
 if(option==7) kind_of_output=' Fan+DDW corrs   to eigenvalues at T=0'

 nkpt_eff=nkpt
!DEBUG
!write(message,'(a,5i5)')' prtvol,iscf,kptopt,nkpt_eff,nkpt_max ',prtvol,iscf,kptopt,nkpt_eff,nkpt_max
!call wrtout(iout,message,'COLL')
!ENDDEBUG
 if( (prtvol==0.or.prtvol==1) .and. (iscf/=-2 .or. kptopt>0)   &
& .and. nkpt_eff>nkpt_max)nkpt_eff=nkpt_max

 if(option==1 .or. (option>=3 .and. option<=7))then

!  Print eigenvalues in hartree for enunit=0 or 2
   if (enunit==0.or.enunit==2) then

     convrt=1.0_dp
     band_index=0

     tmetal=0
     if(option==1 .and. occopt>=3 .and. occopt<=8)tmetal=1

     tmagnet=0
     if(tmetal==1 .and. nsppol==2)then
       tmagnet=1
       rhoup = 0._dp
       rhodn = 0._dp
       nband_index = 1
       do isppol=1,nsppol
         do ikpt=1,nkpt
           nband_k=nband(ikpt+(isppol-1)*nkpt)
           do iband=1,nband_k
             if(isppol==1) rhoup = rhoup + wtk(ikpt)*occ(nband_index)
             if(isppol==2) rhodn = rhodn + wtk(ikpt)*occ(nband_index)
             nband_index = nband_index + 1
           end do
         end do
       end do
       magnet = abs(rhoup - rhodn)
     end if

     if(iscf>0)then
       write(message, '(a,f10.5,a,f10.5)' ) &
&       ' Fermi (or HOMO) energy (hartree) =',fermie,'   Average Vxc (hartree)=',vxcavg
       call wrtout(iout,message,'COLL')
       if (prteig > 0) call wrtout(tmp_unit,message,'COLL')

       if(tmagnet==1)then
         write(message, '(a,es16.8,a,a,es16.8,a,es16.8)' )&
&         ' Magnetisation (Bohr magneton)=',magnet,ch10,&
&         ' Total spin up =',rhoup,'   Total spin down =',rhodn
         call wrtout(iout,message,'COLL')
         if (prteig > 0) call wrtout(tmp_unit,message,'COLL')
       end if
     end if

!    Loop over spins (suppress spin data if nsppol not 2)
     do isppol=1,nsppol

       if (nsppol==2.and.isppol==1) then
         write(message, '(a,a,i4,2x,a)' ) &
&         trim(kind_of_output),' (hartree) for nkpt=',nkpt,'k points, SPIN UP:'
       else if (nsppol==2.and.isppol==2) then
         write(message, '(a,a,i4,2x,a)' )&
&         trim(kind_of_output),' (hartree) for nkpt=',nkpt,'k points, SPIN DOWN:'
       else
         write(message, '(a,a,i4,2x,a)' )&
&         trim(kind_of_output),' (hartree) for nkpt=',nkpt,'k points:'
       end if
       call wrtout(iout,message,'COLL')
       if (prteig > 0) call wrtout(tmp_unit,message,'COLL')

       if(option>=4 .and. option<=7)then
         write(message, '(a)' )&
&         '  (in case of degenerate eigenvalues, averaged derivative)'
         call wrtout(iout,message,'COLL')
         if (prteig > 0) call wrtout(tmp_unit,message,'COLL')
       end if

       do ikpt=1,nkpt
         nband_k=nband(ikpt+(isppol-1)*nkpt)
         if(ikpt<=nkpt_eff)then
           write(message, '(a,i4,a,i3,a,f9.5,a,3f8.4,a)' ) &
&           ' kpt#',ikpt,', nband=',nband_k,', wtk=',wtk(ikpt)+tol10,', kpt=',&
&           kptns(1:3,ikpt)+tol10,' (reduced coord)'
           call wrtout(iout,message,'COLL')
           if (prteig > 0) call wrtout(tmp_unit,message,'COLL')
           do ii=0,(nband_k-1)/8
!            write(message, '(8f15.10)' ) (convrt*eigen(iband+band_index),&
             write(message, '(8f10.5)' ) (convrt*eigen(iband+band_index),&
&             iband=1+ii*8,min(nband_k,8+ii*8))
             call wrtout(iout,message,'COLL')
             if (prteig > 0) call wrtout(tmp_unit,message,'COLL')
           end do
           if(option==1 .and. occopt>=3 .and. occopt<=8)then
             write(message, '(5x,a,i4)' )  ' occupation numbers for kpt#',ikpt
             call wrtout(iout,message,'COLL')
             do ii=0,(nband_k-1)/8
               write(message, '(8f10.5)' ) (occ(iband+band_index),&
&               iband=1+ii*8,min(nband_k,8+ii*8))
               call wrtout(iout,message,'COLL')
             end do
           end if

         else
           if(ikpt==nkpt_eff+1)then
             write(message, '(a,a)' ) &
&             ' prteigrs : prtvol=0 or 1, do not print more k-points.',ch10
             call wrtout(iout,message,'COLL')
           end if
           if (prteig > 0) then
             write(message, '(a,i4,a,i3,a,f9.5,a,3f8.4,a)' ) &
&             ' kpt#',ikpt,', nband=',nband_k,', wtk=',wtk(ikpt)+tol10,', kpt=',&
&             kptns(1:3,ikpt)+tol10,' (reduced coord)'
             call wrtout(tmp_unit,message,'COLL')
             do ii=0,(nband_k-1)/8
               write(message, '(8f10.5)' ) (convrt*eigen(iband+band_index),&
&               iband=1+ii*8,min(nband_k,8+ii*8))
               call wrtout(tmp_unit,message,'COLL')
             end do
           end if
         end if
         band_index=band_index+nband_k
       end do ! do ikpt=1,nkpt
     end do ! do isppol=1,nsppol

!    End print in Hartree
   end if

!  Print in eV for enunit=1 or 2
   if (enunit==1.or.enunit==2) then
     convrt=Ha_eV
     band_index=0

     if(option==1 .and. iscf>0)then
       write(message, '(a,f10.5,a,f10.5)' ) &
&       ' Fermi (or HOMO) energy (eV) =',fermie*convrt,'   Average Vxc (eV)=',vxcavg*convrt
       call wrtout(iout,message,'COLL')
       if (prteig > 0) call wrtout(tmp_unit,message,'COLL')
     end if

!    Loop over spins (suppress spin data if nsppol not 2)
     do isppol=1,nsppol

       if (nsppol==2.and.isppol==1) then
         write(message, '(a,a,i4,2x,a)' ) &
&         trim(kind_of_output),' (   eV  ) for nkpt=',nkpt,'k points, SPIN UP:'
       else if (nsppol==2.and.isppol==2) then
         write(message, '(a,a,i4,2x,a)' ) &
&         trim(kind_of_output),' (   eV  ) for nkpt=',nkpt,'k points, SPIN DOWN:'
       else
         write(message, '(a,a,i4,2x,a)' )&
&         trim(kind_of_output),' (   eV  ) for nkpt=',nkpt,'k points:'
       end if
       call wrtout(iout,message,'COLL')
       if (prteig > 0) call wrtout(tmp_unit,message,'COLL')

       do ikpt=1,nkpt
         nband_k=nband(ikpt+(isppol-1)*nkpt)
         if(ikpt<=nkpt_eff)then
           write(message, '(a,i4,a,i3,a,f9.5,a,3f8.4,a)' ) &
&           ' kpt#',ikpt,', nband=',nband_k,', wtk=',wtk(ikpt)+tol10,', kpt=',&
&           kptns(1:3,ikpt)+tol10,' (reduced coord)'
           call wrtout(iout,message,'COLL')
           if (prteig > 0) call wrtout(tmp_unit,message,'COLL')
           do ii=0,(nband_k-1)/8
             write(message, '(8f10.5)' ) (convrt*eigen(iband+band_index),&
&             iband=1+ii*8,min(nband_k,8+ii*8) )
             call wrtout(iout,message,'COLL')
             if (prteig > 0) call wrtout(tmp_unit,message,'COLL')
           end do
         else
           if(ikpt==nkpt_eff+1)then
             write(message, '(a,a)' ) &
&             ' prteigrs : prtvol=0 or 1, do not print more k-points.',ch10
             call wrtout(iout,message,'COLL')
           end if
           if (prteig > 0) then
             write(message, '(a,i4,a,i3,a,f9.5,a,3f8.4,a)' ) &
&             ' kpt#',ikpt,', nband=',nband_k,', wtk=',wtk(ikpt)+tol10,', kpt=',&
&             kptns(1:3,ikpt)+tol10,' (reduced coord)'
             call wrtout(tmp_unit,message,'COLL')
             do ii=0,(nband_k-1)/8
               write(message, '(8f10.5)' ) (convrt*eigen(iband+band_index),&
&               iband=1+ii*8,min(nband_k,8+ii*8) )
               call wrtout(tmp_unit,message,'COLL')
             end do
           end if
         end if
         band_index=band_index+nband_k
       end do
     end do

!    End print in eV
   end if

 else if(option==2)then

   band_index=0
   do isppol=1,nsppol

     if(nsppol==2)then
       if(isppol==1)write(message, '(2a)' ) ch10,' SPIN UP channel '
       if(isppol==2)write(message, '(2a)' ) ch10,' SPIN DOWN channel '
       call wrtout(iout,message,'COLL')
       if(prteig>0)call wrtout(tmp_unit,message,'COLL')
     end if

     do ikpt=1,nkpt
       nband_k=nband(ikpt+(isppol-1)*nkpt)

       if(ikpt<=nkpt_eff)then
         write(message, '(1x,a,i5,a,f9.5,2f9.5,a)' ) &
&         'Non-SCF case, kpt',ikpt,&
&         ' (',(kptns(ii,ikpt),ii=1,3),'), residuals and eigenvalues='
         call wrtout(iout,message,'COLL')
         if (prteig > 0) then
           write(message, '(1x,a,i5,a,f9.5,2f9.5,a)' ) &
&           'Non-SCF kpt',ikpt,&
&           ' eig(',(kptns(ii,ikpt),ii=1,3),') '
           call wrtout(tmp_unit,message,'COLL')
         end if
         do ii=0,(nband_k-1)/8
           write(message, '(1p,8e10.2)' )&
&           (resid(iband+band_index),iband=1+8*ii,min(8+8*ii,nband_k))
           call wrtout(iout,message,'COLL')
         end do
         do ii=0,(nband_k-1)/6
           write(message, '(1p,6e12.4)' ) &
&           (eigen(iband+band_index),iband=1+6*ii,min(6+6*ii,nband_k))
           call wrtout(iout,message,'COLL')
           if (prteig > 0) call wrtout(tmp_unit,message,'COLL')
         end do
       else
         if(ikpt==nkpt_eff+1)then
           write(message, '(a,a)' ) &
&           ' prteigrs : prtvol=0 or 1, do not print more k-points.',ch10
           call wrtout(iout,message,'COLL')
         end if
         if (prteig > 0) then
           write(message, '(1x,a,i5,a,f9.5,2f9.5,a)' ) &
&           'Non-SCF kpt',ikpt,&
&           ' eig(',(kptns(ii,ikpt),ii=1,3),') '
           call wrtout(tmp_unit,message,'COLL')
           do ii=0,(nband_k-1)/6
             write(message, '(1p,6e12.4)' ) &
&             (eigen(iband+band_index),iband=1+6*ii,min(6+6*ii,nband_k))
             call wrtout(tmp_unit,message,'COLL')
           end do
         end if
       end if

       residk=maxval(resid(band_index+1:band_index+nband_k))
       if (residk>tolwfr) then
         write(message, '(1x,a,2i5,a,1p,e13.5)' ) &
&         ' prteigrs : nnsclo,ikpt=',nnsclo_now,ikpt,&
&         ' max resid (incl. the buffer)=',residk
         call wrtout(iout,message,'COLL')
       end if

       band_index=band_index+nband_k
     end do
   end do

 else
   write(message, '(a,a,a,a,i4,a)' ) ch10,&
&   ' prteigrs : BUG -',ch10,&
&   '  option =',option,', is not an allowed value.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if (prteig > 0) close (tmp_unit)

end subroutine prteigrs
!!***
