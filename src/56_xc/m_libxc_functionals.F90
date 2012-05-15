!{\src2tex{textfont=tt}}
!!****m* ABINIT/libxc_functionals
!! NAME
!!  libxc_functionals
!!
!! FUNCTION
!!  Module containing interfaces to the LibXC library, for exchange
!!  correlation potentials and energies. The interfacing between
!!  the ABINIT and LibXC formats and datastructures happens here.
!!  Also contains basic container datatype for LibXC interfacing.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MOliveira,LHH,FL,GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module libxc_functionals

 use m_profiling

  use defs_basis
#if defined HAVE_DFT_LIBXC
  use xc_f90_types_m
  use libxc_funcs_m
  use xc_f90_lib_m
#endif

  implicit none

#if defined HAVE_DFT_LIBXC
  type libxc_functional
    private
    integer         :: family ! LDA, GGA, etc.
    integer         :: id     ! identifier
    integer         :: nspin  ! # of spin components

    type(xc_f90_pointer_t) :: conf ! the pointer used to call the library
    type(xc_f90_pointer_t) :: info ! information about the functional
  end type libxc_functional

  type(libxc_functional) :: funcs(2)

  private
  public :: libxc_functionals_init, &
&      libxc_functionals_getvxc, &
&      libxc_functionals_isgga, &
&      libxc_functionals_ismgga, &
&      libxc_functionals_nspin, &
&      libxc_functionals_end

contains
!!***

!!****f* libxc_functionals/libxc_functionals_init
!! NAME
!!  libxc_functionals_init
!!
!! FUNCTION
!!  Initialize the desired XC functional, from LibXC.
!!  * Call the LibXC initializer
!!  * Fill preliminary fields in module structures.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      wrtout,xc_f90_gga_exc_vxc,xc_f90_gga_fxc,xc_f90_gga_vxc
!!      xc_f90_lda_exc_vxc,xc_f90_lda_fxc,xc_f90_lda_kxc,xc_f90_lda_vxc
!!      xc_f90_mgga_exc_vxc,xc_f90_mgga_vxc,xc_f90_mgga_x_tb09_set_par
!!
!! SOURCE

  subroutine libxc_functionals_init(ixc,nspden)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_init'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

    implicit none

!Arguments ------------------------------------
!scalars

    integer, intent(in) :: nspden
    integer, intent(in) :: ixc

!Local variables-------------------------------
!scalars

    integer :: i, ii
    character(len=500) :: message
    type(xc_f90_pointer_t) :: str

! *************************************************************************

    funcs(1)%id = -ixc/1000
    funcs(2)%id = -ixc - funcs(1)%id*1000

    funcs(1)%nspin = nspden
    funcs(2)%nspin = nspden

    do i = 1, 2
      if (funcs(i)%id == 0) then
        funcs(i)%family = 0
        cycle
      end if

      ! Get XC functional family
      funcs(i)%family = xc_f90_family_from_id(funcs(i)%id)
      select case (funcs(i)%family)
      case (XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA)
        call xc_f90_func_init(funcs(i)%conf,funcs(i)%info,funcs(i)%id,nspden)
      case default
        write(message, '(4a,i8,2a,i8,6a)' )ch10,&
             &    ' libxc_functionals_init : ERROR -',ch10,&
             &    '  Invalid IXC = ',ixc,ch10,&
             &    '  The LibXC functional family ',funcs(i)%family,&
             &    '  is currently unsupported by ABINIT',ch10,&
             &    '  (-1 means the family is unknown to the LibXC itself)',ch10,&
             &    '  Please consult the LibXC documentation',ch10
        call wrtout(std_out,message,'COLL')
        call leave_new('COLL')
      end select

      if (funcs(i)%id == XC_LDA_C_XALPHA) then
        call xc_f90_lda_c_xalpha_set_par(funcs(i)%conf,zero)
      end if

      ! Dump functional information
      call xc_f90_info_name(funcs(i)%info,message)
      call wrtout(std_out,message,'COLL')
      ii = 0
      call xc_f90_info_refs(funcs(i)%info,ii,str,message)
      do while (ii >= 0)
        call wrtout(std_out,message,'COLL')
        call xc_f90_info_refs(funcs(i)%info,ii,str,message)
      end do
    end do

  end subroutine libxc_functionals_init
!!***

!!****f* libxc_functionals/libxc_functionals_end
!! NAME
!!  libxc_functionals_end
!!
!! FUNCTION
!!  End usage of LibXC functional. Call LibXC end function,
!!  and deallocate module contents.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      wrtout,xc_f90_gga_exc_vxc,xc_f90_gga_fxc,xc_f90_gga_vxc
!!      xc_f90_lda_exc_vxc,xc_f90_lda_fxc,xc_f90_lda_kxc,xc_f90_lda_vxc
!!      xc_f90_mgga_exc_vxc,xc_f90_mgga_vxc,xc_f90_mgga_x_tb09_set_par
!!
!! SOURCE
  subroutine libxc_functionals_end()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_end'
!End of the abilint section

    implicit none

    integer :: i

    do i = 1, 2
      if (funcs(i)%id == 0) cycle
      call xc_f90_func_end(funcs(i)%conf)
    end do

  end subroutine libxc_functionals_end
!!***

!!****f* libxc_functionals/libxc_functionals_isgga
!! NAME
!!  libxc_functionals_isgga
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  is a GGA or not
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  function libxc_functionals_isgga()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_isgga'
!End of the abilint section

    implicit none

!Arguments ------------------------------------

!Local variables-------------------------------

    logical :: libxc_functionals_isgga

! *************************************************************************

    if (any(funcs%family == XC_FAMILY_GGA)) then
      libxc_functionals_isgga = .true.
    else
      libxc_functionals_isgga = .false.
    end if

  end function libxc_functionals_isgga
!!***

!!****f* libxc_functionals/libxc_functionals_ismgga
!! NAME
!!  libxc_functionals_ismgga
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  is a Meta-GGA or not
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  function libxc_functionals_ismgga()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_ismgga'
!End of the abilint section

    implicit none

!Arguments ------------------------------------

!Local variables-------------------------------

    logical :: libxc_functionals_ismgga

! *************************************************************************

    if (any(funcs%family == XC_FAMILY_MGGA)) then
      libxc_functionals_ismgga = .true.
    else
      libxc_functionals_ismgga = .false.
    end if

  end function libxc_functionals_ismgga
!!***

!!****f* libxc_functionals/libxc_functionals_nspin
!! NAME
!!  libxc_functionals_nspin
!!
!! FUNCTION
!!  Returns the number of spin components for the XC functionals
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  function libxc_functionals_nspin()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_nspin'
!End of the abilint section

    implicit none

!Arguments ------------------------------------

!Local variables-------------------------------

    integer :: libxc_functionals_nspin

! *************************************************************************

    if (any(funcs%nspin == XC_POLARIZED)) then
      libxc_functionals_nspin = 2
    else
      libxc_functionals_nspin = 1
    end if

  end function libxc_functionals_nspin
!!***

!!****f* libxc_functionals/libxc_functionals_getvxc
!! NAME
!!  libxc_functionals_getvxc
!!
!! FUNCTION
!!  Return XC potential and energy, from input density (event gradient etc...)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!      wrtout,xc_f90_gga_exc_vxc,xc_f90_gga_fxc,xc_f90_gga_vxc
!!      xc_f90_lda_exc_vxc,xc_f90_lda_fxc,xc_f90_lda_kxc,xc_f90_lda_vxc
!!      xc_f90_mgga_exc_vxc,xc_f90_mgga_vxc,xc_f90_mgga_x_tb09_set_par
!!
!! SOURCE

  subroutine libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho,exc,vxc,grho2,vxcgr,lrho,vxclrho,tau,vxctau,dvxc,d2vxc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'libxc_functionals_getvxc'
 use interfaces_14_hidewrite
!End of the abilint section

    implicit none

!Arguments ------------------------------------

    integer, intent(in) :: ndvxc,nd2vxc,npts,nspden,order
    real(dp),intent(in)  :: rho(npts,nspden)
    real(dp),intent(out) :: vxc(npts,nspden), exc(npts)
    real(dp),intent(in),optional :: grho2(npts,2*min(nspden,2)-1)
    real(dp),intent(out),optional :: vxcgr(npts,3)
    real(dp),intent(in),optional :: lrho(npts,nspden)
    real(dp),intent(out),optional :: vxclrho(npts,nspden)
    real(dp),intent(in),optional :: tau(npts,nspden)
    real(dp),intent(out),optional :: vxctau(npts,nspden)
    real(dp),intent(out),optional :: dvxc(npts,ndvxc)
    real(dp),intent(out),optional :: d2vxc(npts,nd2vxc)

!Local variables-------------------------------

    integer  :: i, ipts
    real(dp) :: c, rhotmp(nspden), exctmp, sigma(3), vsigma(3), vxctmp(nspden)
    real(dp) :: v2rho2(3),v2rhosigma(6),v2sigma2(6),v3rho3(4)
    real(dp) :: lrhotmp(nspden), tautmp(nspden), vxclrhotmp(nspden), vxctautmp(nspden)
    real(dp), allocatable :: gnon(:)
    character(len=500) :: message

! *************************************************************************

    ! Inititalize all relevant arrays to zero
    vxc=zero
    exc=zero
    vxctmp=zero
    exctmp=zero


!LHH,FL,GMR
    v2rho2=zero
    v2rhosigma=zero
    v2sigma2=zero
    v3rho3=zero
    if (order**2 >1) dvxc=zero
    if (order**2 >4) d2vxc=zero
!LHH,FL,GMR

    if (any(funcs%family == XC_FAMILY_GGA)) vxcgr=zero
    if (any(funcs%family == XC_FAMILY_MGGA)) then
      vxcgr=zero
      vxclrho=zero
      vxctau=zero
    end if

    !The TB09 MGGA functional requires an extra quantity
    if (any(funcs%id == XC_MGGA_X_TB09)) then
      ABI_ALLOCATE(gnon,(npts))
      do ipts = 1, npts
        if (sum(rho(ipts, :)) <= 1e-7_dp) then
          gnon(ipts) = zero
        else
          if (nspden == 1) then
            gnon(ipts) = sqrt(grho2(ipts,1))/rho(ipts, 1)
          else
            gnon(ipts) = sqrt(grho2(ipts,3))/sum(rho(ipts, :))
          end if
        end if
      end do
      c = -0.012_dp + 1.023_dp*sqrt(sum(gnon)/npts)
      do i = 1, 2
        if (funcs(i)%id == XC_MGGA_X_TB09) then
          call xc_f90_mgga_x_tb09_set_par(funcs(i)%conf, c)
          write(message, '(2a,f9.6)' ) ch10,&
               &     ' In the functional TB09 c = ', c
          call wrtout(std_out,message,'COLL')
        end if
      end do
      ABI_DEALLOCATE(gnon)
    end if

    !Loop over points
    do ipts = 1, npts

      ! Convert the quantities provided by ABINIT to the ones needed by libxc
      if (nspden == 1) then
        ! ABINIT passes rho_up in the spin-unpolarized case, while the libxc
        ! expects the total density
        rhotmp(1:nspden) = two*rho(ipts,1:nspden)
      else
        rhotmp(1:nspden) = rho(ipts,1:nspden)
      end if
      if (any(funcs%family == XC_FAMILY_GGA) .or. any(funcs%family == XC_FAMILY_MGGA)) then
        sigma=zero
        if (nspden==1) then
          ! ABINIT passes |grho_up|^2 while Libxc needs |grho_tot|^2
          sigma(1) = four*grho2(ipts,1)
        else
          ! ABINIT passes |grho_up|^2, |grho_dn|^2, and |grho_tot|^2
          ! while Libxc needs |grho_up|^2, grho_up.grho_dn, and |grho_dn|^2
          sigma(1) = grho2(ipts,1)
          sigma(2) = (grho2(ipts,3) - grho2(ipts,1) - grho2(ipts,2))/two
          sigma(3) = grho2(ipts,2)
        end if
      end if
      if (any(funcs%family == XC_FAMILY_MGGA)) then
        if (nspden==1) then
          lrhotmp(1:nspden) = two*lrho(ipts,1:nspden)
          tautmp(1:nspden) = four*tau(ipts,1:nspden)
        else
          lrhotmp(1:nspden) = lrho(ipts,1:nspden)
          tautmp(1:nspden) = two*tau(ipts,1:nspden)
        end if
      end if

      !Loop over functionals
      do i = 1,2
        if (funcs(i)%id == 0) cycle

        !Get the potential (and possibly the energy)
        if (iand(xc_f90_info_flags(funcs(i)%info), XC_FLAGS_HAVE_EXC) .ne. 0) then
          select case (funcs(i)%family)
          case (XC_FAMILY_LDA)
            call xc_f90_lda_exc_vxc(funcs(i)%conf,1,rhotmp(1),exctmp,vxctmp(1))
            if (order**2 > 1) then
              call xc_f90_lda_fxc(funcs(i)%conf,1,rhotmp(1),v2rho2(1))
            endif
            if (order**2 > 4) then
              call xc_f90_lda_kxc(funcs(i)%conf,1,rhotmp(1),v3rho3(1))
            endif
          case (XC_FAMILY_GGA)
            call xc_f90_gga_exc_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),exctmp,vxctmp(1),vsigma(1))
            if (order**2 > 1) then
              call xc_f90_gga_fxc(funcs(i)%conf,1,rhotmp(1),sigma(1),v2rho2(1),v2rhosigma(1),v2sigma2(1))
            endif

          case (XC_FAMILY_MGGA)
            call xc_f90_mgga_exc_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),lrhotmp(1),&
                        tautmp(1),exctmp,vxctmp(1),vsigma(1),vxclrhotmp(1),vxctautmp(1))
          end select

        else
          exctmp=zero
          select case (funcs(i)%family)
          case (XC_FAMILY_LDA)
            call xc_f90_lda_vxc(funcs(i)%conf,1,rhotmp(1),vxctmp(1))
          case (XC_FAMILY_GGA)
            call xc_f90_gga_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),vxctmp(1),vsigma(1))


          case (XC_FAMILY_MGGA)
            call xc_f90_mgga_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),lrhotmp(1),&
                        tautmp(1),vxctmp(1),vsigma(1),vxclrhotmp(1),vxctautmp(1))
          end select
        end if

        exc(ipts) = exc(ipts) + exctmp
        vxc(ipts,1:nspden) = vxc(ipts,1:nspden) + vxctmp(1:nspden)

!LHH,FL,GMR: deal with fxc and kxc
        if (order**2>1) then
          select case (funcs(i)%family)
          case (XC_FAMILY_LDA)
            if (nspden==1) then
              if(order>=2) then
                dvxc(ipts,1)=dvxc(ipts,1)+v2rho2(1)
                if(order==3) then
                  d2vxc(ipts,1)=d2vxc(ipts,1)+v3rho3(1)
                endif
              else
                dvxc(ipts,1)=dvxc(ipts,1)+v2rho2(1)
                dvxc(ipts,2)=dvxc(ipts,2)+v2rho2(1)
              endif
            else
              dvxc(ipts,1)=dvxc(ipts,1)+v2rho2(1)
              dvxc(ipts,2)=dvxc(ipts,2)+v2rho2(2)
              dvxc(ipts,3)=dvxc(ipts,3)+v2rho2(3)
              if(order==3) then
                d2vxc(ipts,1)=d2vxc(ipts,1)+v3rho3(1)
                d2vxc(ipts,2)=d2vxc(ipts,2)+v3rho3(2)
                d2vxc(ipts,3)=d2vxc(ipts,3)+v3rho3(3)
                d2vxc(ipts,4)=d2vxc(ipts,4)+v3rho3(4)
              endif
            endif
          case (XC_FAMILY_GGA)
            if (i==1) then
              if (nspden==1) then
                dvxc(ipts,1)=v2rho2(1)*two
                dvxc(ipts,2)=dvxc(ipts,1)
                dvxc(ipts,3)=two*two*vsigma(1)
                dvxc(ipts,4)=dvxc(ipts,3)
                dvxc(ipts,5)=four*two*v2rhosigma(1)
                dvxc(ipts,6)=dvxc(ipts,5)
                dvxc(ipts,7)=two*four*four*v2sigma2(1)
                dvxc(ipts,8)=dvxc(ipts,7)
              else
                dvxc(ipts,1)=v2rho2(1)
                dvxc(ipts,2)=v2rho2(3)
                dvxc(ipts,3)=two*vsigma(1)
                dvxc(ipts,4)=two*vsigma(3)
                dvxc(ipts,5)=two*v2rhosigma(1)
                dvxc(ipts,6)=two*v2rhosigma(6)
                dvxc(ipts,7)=four*v2sigma2(1)
                dvxc(ipts,8)=four*v2sigma2(6)
              end if
            else
              if (nspden==1) then
                dvxc(ipts,9)=v2rho2(1)
                dvxc(ipts,10)=dvxc(ipts,9)
                dvxc(ipts,11)=dvxc(ipts,9)
                dvxc(ipts,12)=two*vsigma(1)
                dvxc(ipts,13)=two*v2rhosigma(1)
                dvxc(ipts,14)=dvxc(ipts,13)
                dvxc(ipts,15)=four*v2sigma2(1)
              else
                dvxc(ipts,9)=v2rho2(1)
                dvxc(ipts,10)=v2rho2(2)
                dvxc(ipts,11)=v2rho2(3)
                dvxc(ipts,12)=two*vsigma(1)
                dvxc(ipts,13)=two*v2rhosigma(1)
                dvxc(ipts,14)=two*v2rhosigma(6)
                dvxc(ipts,15)=four*v2sigma2(1)
              end if
            end if
          end select
        end if

        if (any(funcs%family == XC_FAMILY_GGA) .or. any(funcs%family == XC_FAMILY_MGGA)) then
          !Convert the quantities returned by Libxc to the ones needed by ABINIT
          if (nspden == 1) then
            vxcgr(ipts,3) = vxcgr(ipts,3) + vsigma(1)*two
          else
            vxcgr(ipts,1) = vxcgr(ipts,1) + two*vsigma(1) - vsigma(2)
            vxcgr(ipts,2) = vxcgr(ipts,2) + two*vsigma(3) - vsigma(2)
            vxcgr(ipts,3) = vxcgr(ipts,3) + vsigma(2)
          end if
        end if
        if (any(funcs%family == XC_FAMILY_MGGA)) then
          vxclrho(ipts,1:nspden) = vxclrho(ipts,1:nspden) + vxclrhotmp(1:nspden)
          vxctau(ipts,1:nspden) = vxctau(ipts,1:nspden) + two*vxctautmp(1:nspden)
        end if

      end do

    end do

  end subroutine libxc_functionals_getvxc
#endif

end module
!!***
