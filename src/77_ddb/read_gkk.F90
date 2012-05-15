!{\src2tex{textfont=tt}}
!!****f* ABINIT/read_gkk
!!
!! NAME
!! read_gkk
!!
!! FUNCTION
!! This routine reads in elphon matrix elements and completes them
!!  using the appropriate symmetries
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  elph_ds = datastructure containing elphon matrix elements
!!  Cryst<crystal_structure>=Info on the crystal unit cell.
!!  FSfullpqtofull = mapping of k+q to k
!!  n1wf = number of 1WF files to be read and analyzed
!!  nband = number of bands per kpoint
!!  phon_ds = phonon datastructure with real space interatomic force constants
!!  unitgkk = unit of GKK file for reading
!!
!! OUTPUT
!!  elph_ds = modified gkq
!!  gkk_qpt = el-ph matrix elements for irreducible qpoints and
!!    kpoints (as a function of the reduced symmetry for the qpoint)
!!  gkk_flag = flag array:
!!       -1 -> element is missing
!!        0 -> element is from symmetric qpt (Now done in complete_gkk)
!!        1 -> element is from symmetric pert
!!        2 -> element is kptsym of gkk file
!!        3 -> element was read from gkk file
!!
!! NOTES
!!
!! PARENTS
!!      get_all_gkq
!!
!! CHILDREN
!!      completeperts,gam_mult_displ,gkkfd_free,gkkfd_init,gkkfd_read_h1me
!!      hdr_clean,hdr_io,hdr_nullify,inpphon,insy3,mati3inv,nmsq_pure_gkk_sumfs
!!      normsq_gkq,phdispl_cart2red,prt_gkk_yambo,symq3,wrap2_pmhalf,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine read_gkk(elph_ds,Cryst,Bst,FSfullpqtofull,gkk_flag,n1wf,nband,phon_ds,ep_prt_yambo,unitgkk)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_elphon
 use m_errors

 use m_header,      only : hdr_clean, hdr_nullify
 use m_crystal,     only : crystal_structure
 use m_io_gkk
 use m_gamma

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_gkk'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_59_io_mpi
 use interfaces_72_response
 use interfaces_77_ddb, except_this_one => read_gkk
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1wf,nband,unitgkk,ep_prt_yambo
 type(crystal_structure),intent(in) :: Cryst
 type(bandstructure_type),intent(in) :: Bst
 type(elph_type),intent(inout) :: elph_ds
 type(phon_type),intent(inout) :: phon_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 integer,intent(out) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch,elph_ds%k_phon%nkpt,elph_ds%nsppol,elph_ds%nqptirred)

!Local variables-------------------------------
!scalars
 integer :: nsppol,nbranch,nFSband,minFSband
 integer :: fform,i1wf,ikpt_phon,iatom1,iatom2
 integer :: ib,ib1,ib2,ibb,ibranch,idir,idir1,idir2,ierr,ii,ikpt1
 integer :: ipert,ipert1,ipert2,iqptirred,iqptfull,isppol,isym1
 integer :: itim1,jkpt_phon,new
 integer :: nsym1,qtimrev,rdwr,symikpt_phon,syuse
 integer :: tdonecompl,test_flag,verify
 integer :: nqptirred_local
 real(dp) :: res,ss,timsign
 character(len=500) :: msg
 type(hdr_type) :: hdr1
!arrays
 integer :: FSirrtok(3,elph_ds%k_phon%nkpt)
 integer :: symaf1(Cryst%nsym),symq(4,2,Cryst%nsym)
 integer :: symrc1(3,3,Cryst%nsym),symrl1(3,3,Cryst%nsym)
 integer :: tmpflg(3,Cryst%natom+2,3,Cryst%natom+2)
 real(dp) :: displ_cart(2,3*Cryst%natom,3*Cryst%natom)
 real(dp) :: displ_red(2,3*Cryst%natom,3*Cryst%natom),eigval(3*Cryst%natom)
 real(dp) :: eigvec(2,3*Cryst%natom,3*Cryst%natom),kpt(3),phfrq_tmp(3*Cryst%natom),redkpt(3)
 real(dp) :: qptirred_local(3,n1wf)
 real(dp) :: tnons1(3,Cryst%nsym)
 real(dp),allocatable :: eigen1(:,:,:),gkk_qpt_tmp(:,:,:,:)
 real(dp),allocatable :: h1_mat_el(:,:,:,:,:),h1_mat_el_sq(:,:,:,:,:)
 real(dp),allocatable :: qdata(:,:,:),qdata_tmp(:,:,:,:)

!BEGIN NEW read_gkk
#if 0
 integer :: iq_ibz,jatom,nkwant,pert,pertcase,iqpt_fullbz,ikpt_phonq,ibeff,istat
 real(dp) :: sd1,sd2,e1mef,e2mef
 real(dp),allocatable :: buf_h1(:,:,:,:)
 real(dp),pointer :: kwanted(:,:)
 type(hdr_type) :: GS_Hdr
 type(gkkfd_t) :: Fd
 type(gamma_t) :: Gam

 integer,allocatable :: my_gkk_flag(:,:,:,:)
 real(dp),allocatable :: my_h1_mat_el(:,:,:,:,:),my_h1_mat_el_sq(:,:,:,:,:)
 real(dp),allocatable :: my_accum_mat(:,:,:,:)
 real(dp),allocatable :: my_accum_mat2(:,:,:,:)
 real(dp),allocatable :: gkq_sum_bands(:,:,:)
 real(dp) :: lambda_tot
 real(dp) :: lambda(elph_ds%nsppol)
 real(dp),allocatable :: tmp_mat2(:,:,:) 
 real(dp),allocatable :: zgemm_tmp_mat(:,:,:)
#endif
!END NEW read_gkk

! *************************************************************************

 ABI_UNUSED(Bst%bantot)

 nsppol    = elph_ds%nsppol
 nbranch   = elph_ds%nbranch
 nFSband   = elph_ds%nFSband
 minFSband = elph_ds%minFSband

 call hdr_nullify(hdr1)

 ABI_ALLOCATE(h1_mat_el,(2, nFSband**2, nbranch, elph_ds%k_phon%nkpt, nsppol))
 ierr = ABI_ALLOC_STAT
 if (ierr /= 0 ) then
   write (msg,'(a)')' trying to allocate array h1_mat_el '
   MSG_ERROR(msg)
 end if
 h1_mat_el= zero

 ABI_ALLOCATE(h1_mat_el_sq,(2, nFSband**2, nbranch**2,elph_ds%k_phon%nkpt, nsppol))
 ierr = ABI_ALLOC_STAT
 if (ierr /= 0 ) then
   write (msg,'(a)')' trying to allocate array h1_mat_el_sq '
   MSG_ERROR(msg)
 end if
 h1_mat_el_sq = zero

 ABI_ALLOCATE(elph_ds%qirredtofull,(elph_ds%nqptirred))

!MG array to store the e-ph quantities calculated over the input Q-grid
 ABI_ALLOCATE(qdata_tmp,(elph_ds%nqptirred,nbranch,nsppol,3))
 qdata_tmp=zero

 nqptirred_local=0 !zero number of irred q-points found
 qptirred_local(:,:)=zero

 gkk_flag = -1

 if (elph_ds%gkqwrite ==0) then
   elph_ds%gkk_qpt = zero

 else if (elph_ds%gkqwrite == 1) then
   ABI_ALLOCATE(gkk_qpt_tmp,(2,elph_ds%ngkkband*elph_ds%ngkkband,nbranch**2,nsppol))
   ierr = ABI_ALLOC_STAT
   if (ierr /= 0 ) then
     MSG_ERROR(' trying to allocate array gkk_qpt_tmp')
   end if
   gkk_qpt_tmp = zero
   do iqptirred=1,elph_ds%nqptirred*elph_ds%k_phon%nkpt
     write (elph_ds%unitgkq,REC=iqptirred) gkk_qpt_tmp
   end do
   ABI_DEALLOCATE(gkk_qpt_tmp)

 else
   write (msg,'(a,i0)')' Wrong values for gkqwrite = ',elph_ds%gkqwrite
   MSG_BUG(msg)
 end if !gkqwrite

!===========================================================
!Loop over all files we have
!read in header for perturbation
!should check that all files are complete, have same header
!(taking into account the symmetries for the qpoint),
!represent the correct qpoints ...
!MG: this task should be performed in mrggkk
!===========================================================

 ABI_ALLOCATE(eigen1,(2,nband,nband))
 do i1wf=1,n1wf

   write (msg,'(2a,i4,a,i4)')ch10,' read_gkk : reading 1WF header # ',i1wf,' /',n1wf
   call wrtout(std_out,msg,'COLL')

!  Could check for compatibility of natom, kpt grids, ecut, qpt with DDB grid...
!  MG: Also this task should be done in mrggkk

   rdwr = 5 !read without rewinding
   call hdr_io(fform,hdr1,rdwr,unitgkk)
   if (fform == 0) then
     write (msg,'(a,i0,a)')' 1WF header number ',i1wf,' was mis-read. fform == 0'
     MSG_ERROR(msg)
   end if

   write(msg,'(a,i4)')' read_gkk : have read 1WF header #',i1wf
   call wrtout(std_out,msg,'COLL')
   write (msg,'(2a,i4,a)')ch10,' read_gkk : # of kpt for this perturbation: ',hdr1%nkpt,ch10
   call wrtout(std_out,msg,'COLL')

!  Find qpoint in full grid
   new=1
   do iqptfull=1,elph_ds%nqpt_full
     kpt(:) = hdr1%qptn(:) - elph_ds%qpt_full(:,iqptfull)
     call wrap2_pmhalf(kpt(1),redkpt(1),res)
     call wrap2_pmhalf(kpt(2),redkpt(2),res)
     call wrap2_pmhalf(kpt(3),redkpt(3),res)
     ss=redkpt(1)**2+redkpt(2)**2+redkpt(3)**2
     if(ss < tol6) then
       new = 0
       exit !exit with iqptfull
     end if
   end do !iqptfull

   if (new == 1) then
!    Test should be at the end: dont care if there are additional
!    qpts in gkk file which are not on the main grid. Ignore them.
     write (msg,'(4a,3es16.6,2a)')ch10,&
&     ' read_gkk : WARNING-  ',ch10,&
&     ' qpoint = ',hdr1%qptn(:),ch10,&
&     ' not found in the input q-grid. Ignoring this point '
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,msg,'COLL')

     do isppol=1,hdr1%nsppol
       do ikpt1=1,hdr1%nkpt
         read(unitgkk) ((eigen1(:,ii,ib),ii=1,nband),ib=1,nband)
       end do
     end do
     cycle !cycle the loop on i1wf
   end if !end if (new ==1)

!  Check whether other pieces of the DDB have used this qpt already
   new=1
   do iqptirred=1,nqptirred_local
     kpt(:) = qptirred_local(:,iqptirred) - hdr1%qptn(:)
     call wrap2_pmhalf(kpt(1),redkpt(1),res)
     call wrap2_pmhalf(kpt(2),redkpt(2),res)
     call wrap2_pmhalf(kpt(3),redkpt(3),res)
     ss=redkpt(1)**2+redkpt(2)**2+redkpt(3)**2
     if(ss < tol6) then
       new=0
       exit  !MG We can use this information to avoid recalculating the dynamical matrix 
     end if !but we need to use a fixed format in GKK!
   end do !iqptirred

   if (new==1) then  !we have a new valid irreducible qpoint, add it!
     nqptirred_local = nqptirred_local+1
     if (nqptirred_local > elph_ds%nqptirred) then
       write (msg, '(a,a,a,i6,i6)') &
&       ' found too many qpoints in GKK file wrt anaddb input ', ch10, &
&       ' nqpt_anaddb nqpt_gkk = ', elph_ds%nqptirred, nqptirred_local
       MSG_ERROR(msg)
     end if
     qptirred_local(:,nqptirred_local) = hdr1%qptn(:)
     iqptirred = nqptirred_local
     tdonecompl = 0
     h1_mat_el = zero
   end if

!  now iqptirred is the index of the present qpoint in the array qptirred_local
!  and iqptfull is the index in the full qpt_full array for future reference
   elph_ds%qirredtofull(iqptirred) = iqptfull

   write (msg,'(a,i5,a,3es16.8)')&
&   ' read_gkk : full zone qpt number ',iqptfull,' is ',elph_ds%qpt_full(:,iqptfull)
   call wrtout(std_out,msg,'COLL')

!  if this perturbation has already been filled (overcomplete gkk)
!  check only 1st kpoint and spinpol, then check others
   verify = 0
   if (gkk_flag(hdr1%pertcase,hdr1%pertcase,1,1,iqptirred) /= -1) then
!    
     do isppol=1,nsppol
       do ikpt_phon=1,elph_ds%k_phon%nkpt
         if (gkk_flag(hdr1%pertcase,hdr1%pertcase,ikpt_phon,isppol,iqptirred) == -1) then
           write (std_out,*)" hdr1%pertcase,ikpt_phon,iqptirred",hdr1%pertcase,ikpt_phon,iqptirred
           MSG_ERROR('Partially filled perturbation ')
         end if
       end do ! ikpt_phon
     end do ! isppol
!    
     MSG_WARNING(' gkk perturbation is already filled')
     write(std_out,*)' hdr1%pertcase,iqptirred,iqptfull = ',hdr1%pertcase,iqptirred,iqptfull,&
&     gkk_flag(hdr1%pertcase,hdr1%pertcase,1,1,iqptirred)
     verify = 1
     write (125,*) '# matrix elements for symmetric perturbation'
!    Instead of reading eigen1 into void, verify == 1 checks them later on wrt values in memory
   end if !gkk_flag

!  Examine the symmetries of the q wavevector
!  these will be used to complete the perturbations for other atoms and idir

   call symq3(Cryst%nsym,qptirred_local(:,iqptirred),symq,Cryst%symrec,qtimrev,prtvol=0)

!  Determine dynamical matrix, phonon frequencies and displacement vector for qpoint
   call wrtout(std_out,' read_gkk : calling inpphon to calculate the dynamical matrix','COLL')

   call inpphon(displ_cart,eigval,eigvec,phfrq_tmp,phon_ds,qptirred_local(:,iqptirred))

!  Get displacement vectors for all branches in reduced coordinates
!  used in scalar product with H(1)_atom,idir  matrix elements
!  Calculate $displ_red = displ_cart \cdot gprimd$ for each phonon branch

   call phdispl_cart2red(Cryst%natom,Cryst%gprimd,displ_cart,displ_red)

!  prefactors for gk+q,n\prime;k,n matrix element
!  COMMENT : in decaft there is a weird term in the mass factor, of M-zval(species)
!  dont know why. Not needed to reproduce decaft results, though... 
!  weight is squared in evaluation of
!  gamma_{k,q,j} = 2 \pi omega_{q,j} sum_{nu,nu\prime} |g^{q,j}_{k+q,nu\prime; k,nu}|^2
!  normally cancels with the 2 \pi omega_{q,j} factor in front of the sum...

!  hdr1%pertcase = idir + (ipert-1)*3 where ipert=iatom in the interesting cases
   idir = mod (hdr1%pertcase-1,3)+1
   ipert = int(dble(hdr1%pertcase-idir)/three)+1

   write (msg,'(4a,i3,a,i3,a,i4,a)')ch10,&
&   ' read_gkk : calling insy3 to examine the symmetries of the full perturbation ',ch10,&
&   ' idir = ',idir,' ipert = ',ipert,' and Q point = ',iqptirred,ch10
   call wrtout(std_out,msg,'COLL') 

!  Examine the symmetries of the full perturbation these will be used to complete the kpoints
!  DOESNT USE TIME REVERSAL IN insy3 except for gamma

   syuse=0

   call insy3(Cryst%gprimd,idir,Cryst%indsym,ipert,Cryst%natom,Cryst%nsym,nsym1,2,Cryst%symafm,symaf1,&
&   symq,Cryst%symrec,Cryst%symrel,symrl1,syuse,Cryst%tnons,tnons1)

   do isym1=1,nsym1
     call mati3inv(symrl1(:,:,isym1),symrc1(:,:,isym1))
   end do
   FSirrtok (:,:) = 0

!  ========================================================
!  Loop over irred kpts in file, and fill the default gkk
!  ========================================================

!  MG NOTE : in the present implementation, if nsppol /=1 the code stops in rchkGSheader!
   do isppol=1,hdr1%nsppol !Loop over spins is trivial? Not tested.
     write (std_out,*) ' read_gkk : isppol = ', isppol 

     do ikpt1=1,hdr1%nkpt   !Loop over irred kpoints, WARNING  nkpt depends on qpoint and symmetry!
!      
!      this is the main read of the gkk matrix elements from the file (eigen1 arrays)
!      it has to be done exactly nsppol*nkpt times, and the kpt_phon are completed
!      where appropriate in the loop below (normally succeeding only once for each kpt)
!      
       read(unitgkk) ((eigen1(:,ii,ib),ii=1,nband),ib=1,nband)

!      Check to see if kpoint is in FS set
!      WARNING! the kpoints in the file (kptns) could be ordered arbitrarily
       do ikpt_phon=1,elph_ds%k_phon%nkpt
         kpt(:) = hdr1%kptns(:,ikpt1)-elph_ds%k_phon%kpt(:,ikpt_phon)-qptirred_local(:,iqptirred)
         call wrap2_pmhalf(kpt(1),redkpt(1),res)
         call wrap2_pmhalf(kpt(2),redkpt(2),res)
         call wrap2_pmhalf(kpt(3),redkpt(3),res)

         ss=redkpt(1)**2+redkpt(2)**2+redkpt(3)**2

!        this is not the point on the FS we are looking for
         if (ss > tol6) cycle

         if (verify == 1) then
           do ib1=1,nFSband
             do ib2=1,nFSband
               ibb = (ib1-1)*nFSband+ib2
               write (125,'(2(2E16.6,2x))') h1_mat_el(:,ibb,hdr1%pertcase,ikpt_phon,isppol),&
&               eigen1(:,minFSband-1+ib2,minFSband-1+ib1)
             end do
           end do
         end if !verify

!        if this kpoint has already been filled (overcomplete gkk)
         if (gkk_flag(hdr1%pertcase,hdr1%pertcase,ikpt_phon,isppol,iqptirred) /= -1) then
           MSG_WARNING("gkk element is already filled")
           write(std_out,*)' hdr1%pertcase,ikpt_phon,isppol,iqptirred = ',&
&           hdr1%pertcase,ikpt_phon,isppol,iqptirred,&
&           gkk_flag(hdr1%pertcase,hdr1%pertcase,ikpt_phon,isppol,iqptirred)
           exit
         end if !gkk_flag

!        save this kpoint
         do ib1=1,nFSband
           do ib2=1,nFSband
             ibb = (ib1-1)*nFSband+ib2
             
!            17.05.04 corrected ib1ib2 order in eigen1 indices...

!            real
             res=eigen1(1,minFSband-1+ib2,minFSband-1+ib1)
             h1_mat_el(1,ibb,hdr1%pertcase,ikpt_phon,isppol) = res

!            imag
             res=eigen1(2,minFSband-1+ib2,minFSband-1+ib1)
             h1_mat_el(2,ibb,hdr1%pertcase,ikpt_phon,isppol) = res
           end do !ib2
         end do !ib1
         gkk_flag(hdr1%pertcase,hdr1%pertcase,ikpt_phon,isppol,iqptirred) = 3


!        ===============================================================
!        we now have contribution to g(k+q,k; \kappa,\alpha) from one
!        kpoint,and one perturbation,
!        NB: each perturbation will contribute to all the modes later!
!        find correspondence between this kpt_phon and the others
!        provided sym conserves perturbation as well as qpoint: add to FSirrtok
!        
!        SHOULD ONLY DO THIS FOR THE SYMS NEEDED 
!        TO COMPLETE THE PERTURBATIONS!!!
!        ================================================================

         new=0
         do isym1=1,nsym1
           do itim1=0,qtimrev
             timsign=one-two*itim1
             kpt(:) = timsign*(symrc1(:,1,isym1)*elph_ds%k_phon%kpt(1,ikpt_phon)+&
&             symrc1(:,2,isym1)*elph_ds%k_phon%kpt(2,ikpt_phon)+&
&             symrc1(:,3,isym1)*elph_ds%k_phon%kpt(3,ikpt_phon))
             call wrap2_pmhalf(kpt(1),redkpt(1),res)
             call wrap2_pmhalf(kpt(2),redkpt(2),res)
             call wrap2_pmhalf(kpt(3),redkpt(3),res)
             new=1
!            FIXME: use rank scheme to avoid this loop
             do jkpt_phon=1,elph_ds%k_phon%nkpt
               ss=  (redkpt(1)-elph_ds%k_phon%kpt(1,jkpt_phon))**2&
&               +(redkpt(2)-elph_ds%k_phon%kpt(2,jkpt_phon))**2&
&               +(redkpt(3)-elph_ds%k_phon%kpt(3,jkpt_phon))**2
               if (ss < tol6) then
                 new=0
                 FSirrtok(1,jkpt_phon) = ikpt_phon
                 FSirrtok(2,jkpt_phon) = isym1
                 FSirrtok(3,jkpt_phon) = itim1
                 exit !exit jkpt_phon
               end if
             end do !jkpt_phon
             if (new == 1) then
               write (msg,'(a,3es16.6,a,i5,a,i4,a)')&
&               ' equivalent of kpt ',elph_ds%k_phon%kpt(:,ikpt_phon),' by sym ',isym1,' and itime ',itim1,' was not found'
               MSG_ERROR(msg)
             end if
           end do !itim1
         end do !isim1

!        can we safely exit here? The ikpt has been identified to ikpt_phon,
!        so no need to search the other elph_ds%k_phon%kpt
         exit
       end do !ikpt_phon
     end do !ikpt1
   end do !isppol

   if (verify == 1) cycle

   if (elph_ds%tuniformgrid == 1) then  ! check if irred kpoints found do reconstitute the FS kpts
     do ikpt_phon=1,elph_ds%k_phon%nkpt
       if (FSirrtok(1,ikpt_phon) == 0) then
         write(msg,'(a,3es16.6,2a)')&
&         ' kpt = ',elph_ds%k_phon%kpt(:,ikpt_phon),ch10,&
&         ' is not the symmetric of one of those found in the GKK file'
         MSG_ERROR(msg)
       end if
     end do !ikpt_phon
   end if

!  ===============================================================
!  We now have all irred kpoints : complete the others
!  complete gkk for symmetric ikpt_phon with sym1 which conserve
!  the full perturbation+qpoint
!  Not tested explicitly, but the results for Pb using reduced kpts look good
!  should do same RF calculation with nsym=1 and check
!  ===============================================================

   do symikpt_phon=1,elph_ds%k_phon%nkpt
!    if the element has already been filled with another sym op, cycle
     if (gkk_flag(hdr1%pertcase,hdr1%pertcase,symikpt_phon,1,iqptirred) /= -1) cycle

     ikpt_phon=FSirrtok(1,symikpt_phon)
     isym1    =FSirrtok(2,symikpt_phon)
     itim1    =FSirrtok(3,symikpt_phon)
     timsign  = one-two*itim1

!    copy kpt to symmetric kpoint
     do ibb=1,nFSband**2
       h1_mat_el(1,ibb,hdr1%pertcase,symikpt_phon,:) = h1_mat_el(1,ibb,hdr1%pertcase,ikpt_phon,:)
       h1_mat_el(2,ibb,hdr1%pertcase,symikpt_phon,:) = h1_mat_el(2,ibb,hdr1%pertcase,ikpt_phon,:)
     end do
     gkk_flag(hdr1%pertcase,hdr1%pertcase,symikpt_phon,:,iqptirred) = 2
   end do !symikpt_phon

!  do ii=1,elph_ds%k_phon%nkpt
!  do ibb=1,nFSband**2
!  write(444,*)hdr1%pertcase,ibb,ii,h1_mat_el(:,ibb,hdr1%pertcase,ii,1)
!  end do
!  end do

   if (elph_ds%tuniformgrid == 1) then
!    normally at this point we have used all the gkk for all kpoints on the FS
!    for the given irred perturbation: check
     do ikpt_phon=1,elph_ds%k_phon%nkpt
       if (gkk_flag(hdr1%pertcase,hdr1%pertcase,ikpt_phon,1,iqptirred) == -1) then
         write (msg,'(a,i3,a,3es18.6,2a,i3,a,i3,a,3es18.6,a,a,i4,a,a)')&
&         ' For full qpt ', iqptirred,') ',qptirred_local(:,iqptirred),ch10, &
&         ' the gkk element : pertcase = ',hdr1%pertcase,' ikpt = ',ikpt_phon, &
&         ' kpt = ',elph_ds%k_phon%kpt(:,ikpt_phon),ch10,&
&         ' and isppol ',1,ch10,&
&         ' was not found by symmetry operations on the irreducible kpoints given'
         MSG_ERROR(msg)
       end if
       if (FSirrtok(1,ikpt_phon) == 0) then
         write (msg,'(a)')" One of the kpoints was not equivalent to an irred one found in the gkk file "
         MSG_ERROR(msg)
       end if
     end do !ikpt_phon
   end if

   write(msg,'(a,i0)')' read_gkk : Done completing the kpoints for pertcase ',hdr1%pertcase
   call wrtout(std_out,msg,'COLL')

   tmpflg(:,:,:,:) = 0

   do idir1=1,3
     do iatom1=1,Cryst%natom
       ipert1 = (iatom1-1)*3+idir1
       do idir2=1,3
         do iatom2=1,Cryst%natom
           ipert2 = (iatom2-1)*3+idir2
           if (gkk_flag(ipert1,ipert1,1,1,iqptirred) >= 0 .and. &
&           gkk_flag(ipert2,ipert2,1,1,iqptirred) >= 0) then
             tmpflg(idir1,iatom1,idir2,iatom2) = 1
           end if
         end do
       end do
     end do
   end do


!  ===============================================
!  Full test: need all perturbations explicitly
!  ===============================================

   test_flag = 0
   if (sum(tmpflg(:,1:Cryst%natom,:,1:Cryst%natom)) == (3*Cryst%natom)**2 .and. tdonecompl == 0) test_flag = 1

   write(std_out,*)'read_gkk: tdonecompl = ', tdonecompl

!  de-activate completion of perts by symmetry for now.
!  Must be called when all irreducible perturbations are in memory!!!!
   if (test_flag == 1 .and. tdonecompl == 0) then

!    write(std_out,*) ' read_gkk : enter fxgkkphase before completeperts'
!    call fxgkkphase(elph_ds,gkk_flag,h1_mat_el,iqptirred)

     if (ep_prt_yambo==1) then
       call prt_gkk_yambo(displ_cart,displ_red,elph_ds%k_phon%kpt,h1_mat_el,iqptirred,&
&       Cryst%natom,nFSband,elph_ds%k_phon%nkpt,phfrq_tmp,hdr1%qptn)
     end if

!    ========================================================================
!    Now use more general symops to complete the other equivalent
!    perturbations: the kpoints are also shuffled by these symops
!    afterwards h1_mat_el_sq contains gamma_\tau\alpha,\tau'\alpha' in reduced coordinates
!    
!    \gamma_{\tau'\alpha',\tau\alpha} =
!    <psi_{k+q,ib2}| H(1)_{\tau'\alpha'}| psi_{k,ib1}>* \cdot  
!    <psi_{k+q,ib2}| H(1)_{\tau \alpha }| psi_{k,ib1}>
!    
!    ========================================================================
!    write(554,*)h1_mat_el

     call completeperts(Cryst,nbranch,nFSband,elph_ds%k_phon%nkpt,nsppol,&
&     gkk_flag(:,:,:,:,iqptirred),h1_mat_el,h1_mat_el_sq,qptirred_local(:,iqptirred),symq,qtimrev)

!    write(776,*)h1_mat_el_sq

     tdonecompl = 1
   end if

!  ==============================================================
!  if we have all the perturbations for this qpoint, proceed
!  with scalar product, norm squared, and add weight factors
!  
!  SHOULD HAVE A TEST SO h1_mat_el IS NOT OVERWRITTEN
!  BEFORE PREVIOUS QPOINT IS FINISHED!!!!!
!  REMARK: there is just a factor of |rprimd(:,i)| between
!  abinit and decaft matrix elements.
!  ==============================================================

   test_flag = 1
   do isppol=1,nsppol
     do ikpt_phon=1,elph_ds%k_phon%nkpt
       do ibranch=1,nbranch
         if (gkk_flag (ibranch,ibranch,ikpt_phon,isppol,iqptirred) == -1) then
           test_flag = 0
           exit
         end if
       end do
     end do
   end do

   if (test_flag /= 0) then
     call wrtout(std_out,' read_gkk : enter normsq_gkq',"COLL")

!    MG temporary array to save ph-linewidths before Fourier interpolation
     ABI_ALLOCATE(qdata,(nbranch,nsppol,3))
     qdata(:,:,:)=zero

     call normsq_gkq(displ_red,eigvec,elph_ds,FSfullpqtofull,&
&     h1_mat_el_sq,iqptirred,phfrq_tmp,qptirred_local,qdata)

     qdata_tmp(iqptirred,:,:,:)=qdata(:,:,:)
     ABI_DEALLOCATE(qdata)
   end if

   call hdr_clean(hdr1)

 end do !of i1wf 

!got all the gkk perturbations

 ABI_DEALLOCATE(eigen1)
 ABI_DEALLOCATE(h1_mat_el)
 ABI_DEALLOCATE(h1_mat_el_sq)

 if (nqptirred_local /= elph_ds%nqptirred) then
   write (msg, '(3a,i0,i0)') &
&   ' Found wrong number of qpoints in GKK file wrt anaddb input ', ch10, &
&   ' nqpt_anaddb nqpt_gkk = ', elph_ds%nqptirred, nqptirred_local
   MSG_ERROR(msg)
 end if

!normally at this point we have the gkk for all kpoints on the FS
!for all the perturbations. Otherwise a 1WF file is missing.
!NOTE: still havent checked the qpoint grid completeness
 do iqptirred=1,elph_ds%nqptirred
   do isppol=1,nsppol
     do ikpt_phon=1,elph_ds%k_phon%nkpt
       do ipert=1,nbranch
         if (gkk_flag(ipert,ipert,ikpt_phon,isppol,iqptirred) == -1) then
           write (msg,'(a,i5,1x,i5,1x,i5,1x,i5,a,a)')&
&           ' gkk element',ipert,ikpt_phon,isppol,iqptirred,' was not found by symmetry operations ',&
&           ' on the irreducible perturbations and qpoints given'
           MSG_ERROR(msg)
         end if
       end do !ipert
     end do !ikpt_phon
   end do !isppol
 end do !iqptirred

 call wrtout(std_out,'read_gkk : done completing the perturbations (and checked!)','COLL')

!MG save phonon frequencies, ph-linewidths and lambda(q,n) values before Fourier interpolation
 ABI_ALLOCATE(elph_ds%qgrid_data,(elph_ds%nqptirred,nbranch,nsppol,3))

 do iqptirred=1,elph_ds%nqptirred
   elph_ds%qgrid_data(iqptirred,:,:,:)=qdata_tmp(iqptirred,:,:,:)
 end do

 ABI_DEALLOCATE(qdata_tmp)

!BEGIN NEW_GKK
#if 0
 call hdr_nullify(GS_Hdr)

 MSG_WARNING("Entering NEW_GKK section")
 call gkkfd_init(Fd,"Dummy_file",GS_Hdr,funt=unitgkk)
!
 nkwant  =  elph_ds%k_phon%nkpt
 kwanted => elph_ds%k_phon%kpt

 ABI_ALLOCATE(my_gkk_flag,(nbranch,nbranch,nkwant,nsppol))

 ABI_ALLOCATE(my_h1_mat_el   ,(2, nFSband**2, nbranch,   nkwant, nsppol))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory my_h1_mat_el")

 ABI_ALLOCATE(my_h1_mat_el_sq,(2, nFSband**2, nbranch**2,nkwant, nsppol))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory my_h1_mat_el_sq")

!call gamma_init(Gam,Cryst,gprim,elph_ds%ep_scalprod,qptrlatt,nqibz,qibz,nqbz,qbz,nsppol,nrpt,rpt,wghatm)
!call gamma_free(Gam)

 do iq_ibz=1,Fd%nqibz
   iqptfull = elph_ds%qirredtofull(iq_ibz)

   write (msg,'(a,i5,a,3es16.8)')&
&   ' NEW_gkk : full zone qpt number ',iqptfull,' is ',elph_ds%qpt_full(:,iqptfull)
   call wrtout(std_out,msg,'COLL')

   my_h1_mat_el    = HUGE(zero)
   my_h1_mat_el_sq = HUGE(zero)
!  my_h1_mat_el    = zero
!  my_h1_mat_el_sq = zero
   my_gkk_flag     = -1
!  
!  Examine the symmetries of the q wavevector
!  these will be used to complete the perturbations for other atoms and idir
   call symq3(Cryst%nsym,Fd%qibz(:,iq_ibz),symq,Cryst%symrec,qtimrev,prtvol=0)
!  
!  Determine dynamical matrix, phonon frequencies and displacement vector for qpoint
   call inpphon(displ_cart,eigval,eigvec,phfrq_tmp,phon_ds,Fd%qibz(:,iq_ibz))
!  
!  Displacement vectors in reduced coordinates
   call phdispl_cart2red(Cryst%natom,Cryst%gprimd,displ_cart,displ_red)
!  
   do pert=1,Fd%npert_q(iq_ibz)
     pertcase = Fd%pertcase_q(pert,iq_ibz)
     ABI_ALLOCATE(buf_h1,(2,nFSband**2,nkwant,nsppol))
     istat = ABI_ALLOC_STAT
     ABI_CHECK(istat==0,"out of memory buf_h1")      

     call gkkfd_read_h1me(Fd,iq_ibz,pertcase,nFSband,minFSband,nkwant,kwanted,Cryst,buf_h1)

     do ii=1,nkwant !elph_ds%k_phon%nkpt
       do ibb=1,nFSband**2
         write(445,*)pertcase,ibb,ii,buf_h1(:,ibb,ii,1)
       end do
     end do

     my_gkk_flag(pertcase,pertcase,:,:) = 3 ! read from GKK file.

!    Save values
     do isppol=1,nsppol
       my_h1_mat_el(:,:,pertcase,:,isppol) = buf_h1(:,:,:,isppol)
     end do

     ABI_DEALLOCATE(buf_h1)
   end do
!  
!  ========================================================================
!  Now use more general symops to complete the other equivalent
!  perturbations: the kpoints are also shuffled by these symops
!  afterwards h1_mat_el_sq contains gamma_\tau\alpha,\tau'\alpha' in reduced coordinates
!  
!  gamma_{\tau'\alpha',\tau\alpha} =
!  <psi_{k+q,ib2} | H(1)_{\tau'\alpha'} | psi_{k,ib1}>*  \cdot
!  <psi_{k+q,ib2} | H(1)_{\tau \alpha } | psi_{k,ib1}>
!  ========================================================================
!  
!  write(555,*)my_h1_mat_el

   call completeperts(Cryst,nbranch,nFSband,nkwant,nsppol,&
&   my_gkk_flag,my_h1_mat_el,my_h1_mat_el_sq,Fd%qibz(:,iq_ibz),symq,qtimrev)

!  write(777,*)my_h1_mat_el_sq

   ABI_ALLOCATE(my_accum_mat ,(2,nbranch,nbranch,nsppol))
   my_accum_mat =zero
   ABI_ALLOCATE(my_accum_mat2,(2,nbranch,nbranch,nsppol))
   my_accum_mat2=zero

#if 1
   call nmsq_pure_gkk_sumFS(my_accum_mat,my_accum_mat2,displ_red,elph_ds,FSfullpqtofull,my_h1_mat_el_sq,iq_ibz)
#else
!  iqpt_fullbz = elph_ds%qirredtofull(iqptirred)
   iqpt_fullbz = elph_ds%qirredtofull(iq_ibz)
   ABI_ALLOCATE(gkq_sum_bands,(2,nbranch,nbranch))
!  
!  accum_mat and accum_mat2 are real, the imaginary part is used for debugging purpose
!  accum_mat2 is used to store the phonon-linewidhts before interpolation
!  FIXME: 
!  BE careful here, since it wont' work if kwanted /= elph_ds%k_phon%kpt
!  due to the tables and the way used to access the eigevalues.
   do isppol=1,nsppol
     do ikpt_phon=1,nkwant
!      
!      The index of k+q in the BZ.
       ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)

!      % ikptgs = irredtoGS_phon(elph_ds%k_phon%full2irr(1,ikpt_phon))
!      
!      gkq_sum_bands = 
!      \sum_{ib1,ib2} <k+q| H^{(1)}_{q,\tau_i,\alpha_i} |k> \cdot <k| H^{(1)}_{q,\tau_j,\alpha_j}|k+q>
!      
!      where ibranch = (\tau_i,\alpha_i) and  jbranch = (\tau_j,\alpha_j).
       gkq_sum_bands(:,:,:) = zero

       do ib1=1,nFSband
!        %% e1mef = Bst%eig(ikptgs,minFSband-1+ib1,isppol) - Bst%fermie 
!        %% sd1 = gaussian(e2mef,sigma)
         sd1 = elph_ds%k_phon%wtk(ib1,ikpt_phon,isppol)      !  weights for distance from the fermi surface

         do ib2=1,nFSband
!          %% e2mef = Bst%eig(ikptgs,minFSband-1+ib2,isppol) - Bst%fermie 
!          %% sd2 = gaussian(e2mef,sigma)
           sd2 = elph_ds%k_phon%wtk(ib2,ikpt_phonq,isppol)  !  weights for distance from the fermi surface
           ibeff=ib2+(ib1-1)*nFSband

           gkq_sum_bands = gkq_sum_bands + &
&           sd1*sd2*pi * reshape(my_h1_mat_el_sq(:,ibeff,:,ikpt_phon,isppol),(/2,nbranch,nbranch/))
         end do !ib2
       end do !ib1
!      
!      gamma matrix contribution in reduced coordinates (ie interpolatable form)
!      The sum over Fermi surface bands is done here, and fed into (ib1,ib2)=(1,1)
       my_h1_mat_el_sq(:,1,:,ikpt_phon,isppol) = reshape(gkq_sum_bands,(/2,nbranch**2/))

       my_accum_mat(:,:,:,isppol) = my_accum_mat(:,:,:,isppol) + gkq_sum_bands(:,:,:)
     end do ! kpt_phon
   end do ! isppol

   ABI_DEALLOCATE(gkq_sum_bands)
!  
!  MG20060603
!  scalar product wit displ_red to calculate the ph lwdth before interpolation (stored in accum_mat2)
   ABI_ALLOCATE(tmp_mat2,(2,elph_ds%nbranch,elph_ds%nbranch))
   ABI_ALLOCATE(zgemm_tmp_mat,(2,elph_ds%nbranch,elph_ds%nbranch))

   do isppol=1,nsppol
     zgemm_tmp_mat = my_accum_mat(:,:,:,isppol)
!    
     call gam_mult_displ(nbranch, displ_red, zgemm_tmp_mat, tmp_mat2)

     do ipert1=1,nbranch
       my_accum_mat2(1,ipert1,ipert1,isppol) = my_accum_mat2(1,ipert1,ipert1,isppol) + tmp_mat2(1,ipert1,ipert1)
     end do
!    
   end do

   ABI_DEALLOCATE(tmp_mat2)
   ABI_DEALLOCATE(zgemm_tmp_mat)
#endif

!  MG: values without the good prefactor
   my_accum_mat = my_accum_mat * elph_ds%occ_factor/nkwant 

!  MG: my_accum_mat2 contains the line-widhts before the Fourier interpolation
   my_accum_mat2 = my_accum_mat2 * elph_ds%occ_factor/nkwant 

!  MG20060531i
!  write e-ph quantities before Fourier interpolation
!  save e-ph values in the temporary array qdata that will be copied into elph_ds%qgrid_data

   write (msg,'(4a,3es16.6,63a)')ch10,                    &
&   ' NEW_GKK : Phonon linewidths before interpolation ',ch10,&
&   ' Q point = ',Fd%qibz(:,iq_ibz),ch10,('=',ii=1,60),ch10,  &
&   ' Mode          Frequency (Ha)  Linewidth (Ha)  Lambda '
   call wrtout(std_out,msg,'COLL')

   lambda_tot = zero
   do isppol=1,nsppol
     do ii=1,nbranch
       lambda(isppol)=zero
!      MG: the tolerance factor is somehow arbitrary
       if (abs(phfrq_tmp(ii)) > tol10) lambda(isppol)=my_accum_mat2(1,ii,ii,isppol)/ (pi*elph_ds%n0(isppol)*phfrq_tmp(ii)**2)
       lambda_tot=lambda_tot+lambda(isppol)
       write(msg,'(i8,es20.6,2es16.6)' )ii,phfrq_tmp(ii),my_accum_mat2(1,ii,ii,isppol),lambda(isppol)
       call wrtout(std_out,msg,'COLL')
!      save values
!      qdata(ii,isppol,1)=phfrq_tmp(ii)
!      qdata(ii,isppol,2)=my_accum_mat2(1,ii,ii,isppol)
!      qdata(ii,isppol,3)=lambda(isppol)
     end do !loop over branch
   end do !loop over sppol

!  normalize for number of spins
   lambda_tot = lambda_tot / nsppol

   write(msg,'(61a,44x,es16.6,62a)' )('=',ii=1,60),ch10,lambda_tot,ch10,('=',ii=1,60),ch10
   call wrtout(std_out,msg,'COLL')
!  ENDMG20060531

!  immediately calculate linewidths:
   write(std_out,*) 'summed my_accum_mat = '
   write(std_out,'(3(2E18.6,1x))') my_accum_mat(:,:,:,1)
   write(std_out,*) 'summed my_accum_mat2 = '
   write(std_out,'(3(2E18.6,1x))')  (my_accum_mat2(:,ii,ii,1),ii=1,nbranch)
   write(std_out,*) 'displ_red  = '
   write(std_out,'(3(2E18.6,1x))') displ_red

   ABI_DEALLOCATE(my_accum_mat)
   ABI_DEALLOCATE(my_accum_mat2)
 end do ! iq_ibz

 ABI_DEALLOCATE(my_gkk_flag)
 ABI_DEALLOCATE(my_h1_mat_el)
 ABI_DEALLOCATE(my_h1_mat_el_sq)

 call gkkfd_free(Fd)
#endif
!END NEW_GKK

end subroutine read_gkk
!!***
