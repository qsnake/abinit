

AC_DEFUN([ABI_CC_HINTS],[
  dnl Init
  abi_cc_vendor_hnt="none"
  abi_cc_version_hnt="none"
  abi_sys_spec_hnt="none"

  dnl Look for hints flags
  AC_MSG_CHECKING([which cc hints to apply])

  dnl Case built from config/hints/cc_*.conf
  case "${abi_cc_vendor}" in
    intel)
      abi_cc_vendor_hnt="intel"
      abi_cc_version_hnt="default"
      abi_sys_spec_hnt="default"
      CFLAGS_HINTS='-vec-report0'
      CFLAGS_PIC='-fPIC'
      CC_LDFLAGS_HINTS='-static-libgcc -static-intel'
      ;;
    ibm)
      abi_cc_vendor_hnt="ibm"
      abi_cc_version_hnt="default"
      abi_sys_spec_hnt="default"
      CFLAGS_64BITS='-q64'
      CFLAGS_PIC='-qpic'
      ;;
    open64)
      abi_cc_vendor_hnt="open64"
      abi_cc_version_hnt="default"
      abi_sys_spec_hnt="default"
      CFLAGS_64BITS='-m64'
      CFLAGS_PIC='-fPIC'
      CFLAGS_32BITS='-m32'
      ;;
    gnu)
      abi_cc_vendor_hnt="gnu"
      abi_cc_version_hnt="default"
      case "${abi_sys_spec}" in
        irix-mips)
          abi_sys_spec_hnt="irix-mips"
          CFLAGS_64BITS='-mabi=64'
          CFLAGS_PIC='-fPIC'
          CFLAGS_32BITS='-mabi=32'
          ;;
        *)
          abi_sys_spec_hnt="default"
          CFLAGS_64BITS='-m64'
          CFLAGS_PIC='-fPIC'
          CFLAGS_32BITS='-m32'
          ;;
      esac   # [case: abi_sys_spec, indent: 2, item: True]
      ;;
    pathscale)
      abi_cc_vendor_hnt="pathscale"
      abi_cc_version_hnt="default"
      abi_sys_spec_hnt="default"
      CFLAGS_64BITS='-m64'
      CFLAGS_PIC='-fPIC'
      CFLAGS_32BITS='-m32'
      ;;
    compaq)
      abi_cc_vendor_hnt="compaq"
      abi_cc_version_hnt="default"
      abi_sys_spec_hnt="default"
      big_endian='-convert big_endian'
      ;;
  esac   # [case: abi_cc_vendor, indent: 0, item: True]

  dnl Display settings
  AC_MSG_RESULT([${abi_cc_vendor_hnt}/${abi_cc_version_hnt}/${abi_sys_spec_hnt}])

]) #ABI_CC_HINTS


AC_DEFUN([ABI_FC_HINTS],[
  dnl Init
  abi_fc_vendor_hnt="none"
  abi_fc_version_hnt="none"
  abi_sys_spec_hnt="none"

  dnl Look for hints flags
  AC_MSG_CHECKING([which fc hints to apply])

  dnl Case built from config/hints/fc_*.conf
  case "${abi_fc_vendor}" in
    absoft)
      abi_fc_vendor_hnt="absoft"
      abi_fc_version_hnt="default"
      abi_sys_spec_hnt="default"
      abi_fc_wrap='yes'
      FCFLAGS_FIXEDFORM='-ffixed'
      FCFLAGS_FREEFORM='-ffree'
      FCFLAGS_MODDIR='-p $(abinit_moddir)'
      ;;
    compaq)
      abi_fc_vendor_hnt="compaq"
      abi_fc_version_hnt="default"
      abi_sys_spec_hnt="default"
      FCFLAGS_BIGENDIAN='-convert big_endian'
      FCFLAGS_FIXEDFORM='-fixed'
      FCFLAGS_FREEFORM='-free'
      FCFLAGS_MODDIR='-module $(abinit_moddir) -I$(abinit_moddir)'
      ;;
    gnu)
      abi_fc_vendor_hnt="gnu"
      case "${abi_fc_version}" in
        4.5)
          abi_fc_version_hnt="4.5"
          abi_sys_spec_hnt="default"
          FCFLAGS_HINTS='-ffree-line-length-none -Wstrict-aliasing=2'
          FCFLAGS_64BITS='-m64'
          FCFLAGS_PIC='-fPIC'
          FCFLAGS_FREEFORM='-ffree-form'
          FCFLAGS_32BITS='-m32'
          FCFLAGS_BIGENDIAN='-fconvert=big-endian -frecord-marker=4'
          FCFLAGS_FIXEDFORM='-ffixed-form'
          FCFLAGS_MODDIR='-J$(abinit_moddir)'
          ;;
        4.1)
          abi_fc_version_hnt="4.1"
          abi_sys_spec_hnt="default"
          FCFLAGS_64BITS='-m64'
          FCFLAGS_FREEFORM='-ffree-form'
          FC_LDFLAGS_HINTS='-static-libgfortran'
          FCFLAGS_32BITS='-m32'
          FCFLAGS_FIXEDFORM='-ffixed-form'
          FCFLAGS_MODDIR='-J$(abinit_moddir) -I$(abinit_moddir)'
          ;;
        *)
          abi_fc_version_hnt="default"
          abi_sys_spec_hnt="default"
          FCFLAGS_HINTS='-ffree-line-length-none'
          FCFLAGS_64BITS='-m64'
          FCFLAGS_PIC='-fPIC'
          FCFLAGS_FREEFORM='-ffree-form'
          FCFLAGS_32BITS='-m32'
          FCFLAGS_BIGENDIAN='-fconvert=big-endian -frecord-marker=4'
          FCFLAGS_FIXEDFORM='-ffixed-form'
          FCFLAGS_MODDIR='-J$(abinit_moddir)'
          ;;
      esac   # [case: abi_fc_version, indent: 2, item: True]
      ;;
    open64)
      abi_fc_vendor_hnt="open64"
      abi_fc_version_hnt="default"
      abi_sys_spec_hnt="default"
      FCFLAGS_HINTS='-extend-source'
      FCFLAGS_64BITS='-m64 -align64'
      FCFLAGS_PIC='-fPIC'
      FCFLAGS_FREEFORM='-freeform'
      FCFLAGS_32BITS='-m32 -align32'
      FCFLAGS_BIGENDIAN='-byteswapio'
      FCFLAGS_FIXEDFORM='-fixedform'
      FCFLAGS_MODDIR='-module $(abinit_moddir)'
      ;;
    sun)
      abi_fc_vendor_hnt="sun"
      abi_fc_version_hnt="default"
      abi_sys_spec_hnt="default"
      FCFLAGS_HINTS='-e'
      FCFLAGS_64BITS='-m64'
      FCFLAGS_PIC='-KPIC'
      FCFLAGS_FREEFORM='-free'
      FCFLAGS_32BITS='-m32'
      FCFLAGS_BIGENDIAN='-xfilebyteorder=big16:%all'
      FCFLAGS_FIXEDFORM='-fixed'
      FCFLAGS_MODDIR='-moddir=$(abinit_moddir) -M$(abinit_moddir)'
      ;;
    amd)
      abi_fc_vendor_hnt="amd"
      abi_fc_version_hnt="default"
      abi_sys_spec_hnt="default"
      FCFLAGS_HINTS='-extend-source -col120'
      FCFLAGS_64BITS='-m64 -align64 -default64'
      FCFLAGS_PIC='-fPIC'
      FCFLAGS_FREEFORM='-freeform'
      FCFLAGS_32BITS='-m32 -align32'
      FCFLAGS_BIGENDIAN='-convert big_endian'
      FCFLAGS_FIXEDFORM='-fixedform'
      FCFLAGS_MODDIR='-module $(abinit_moddir) -I$(abinit_moddir)'
      ;;
    intel)
      abi_fc_vendor_hnt="intel"
      case "${abi_fc_version}" in
        10.1)
          abi_fc_version_hnt="10.1"
          case "${abi_sys_spec}" in
            *-ia64)
              abi_sys_spec_hnt="*-ia64"
              FCFLAGS_HINTS='-extend_source'
              FCFLAGS_PIC='-fPIC'
              FCFLAGS_FREEFORM='-free'
              FC_LDFLAGS_HINTS='-static-intel -static-libgcc'
              FCFLAGS_BIGENDIAN='-convert big_endian'
              FCFLAGS_FIXEDFORM='-fixed'
              FCFLAGS_MODDIR='-module $(abinit_moddir)'
              ;;
            *)
              abi_sys_spec_hnt="default"
              FCFLAGS_HINTS='-extend-source -vec-report0'
              FCFLAGS_PIC='-fPIC'
              FCFLAGS_FREEFORM='-free'
              FC_LDFLAGS_HINTS='-static-intel -static-libgcc'
              FCFLAGS_BIGENDIAN='-convert big_endian'
              FCFLAGS_FIXEDFORM='-fixed'
              FCFLAGS_MODDIR='-module $(abinit_moddir)'
              ;;
          esac   # [case: abi_sys_spec, indent: 4, item: True]
          ;;
        10.0)
          abi_fc_version_hnt="10.0"
          case "${abi_sys_spec}" in
            *-ia64)
              abi_sys_spec_hnt="*-ia64"
              FCFLAGS_HINTS='-extend_source'
              FCFLAGS_PIC='-fPIC'
              FCFLAGS_FREEFORM='-free'
              FC_LDFLAGS_HINTS='-static-intel -static-libgcc'
              FCFLAGS_BIGENDIAN='-convert big_endian'
              FCFLAGS_FIXEDFORM='-fixed'
              FCFLAGS_MODDIR='-module $(abinit_moddir)'
              ;;
            *)
              abi_sys_spec_hnt="default"
              FCFLAGS_HINTS='-extend-source -vec-report0'
              FCFLAGS_PIC='-fPIC'
              FCFLAGS_FREEFORM='-free'
              FC_LDFLAGS_HINTS='-static-intel -static-libgcc'
              FCFLAGS_BIGENDIAN='-convert big_endian'
              FCFLAGS_FIXEDFORM='-fixed'
              FCFLAGS_MODDIR='-module $(abinit_moddir)'
              ;;
          esac   # [case: abi_sys_spec, indent: 4, item: True]
          ;;
        *)
          abi_fc_version_hnt="default"
          case "${abi_sys_spec}" in
            *-ia64)
              abi_sys_spec_hnt="*-ia64"
              FCFLAGS_HINTS='-extend_source'
              FCFLAGS_PIC='-fPIC'
              FCFLAGS_FREEFORM='-free'
              FC_LDFLAGS_HINTS='-static-intel -static-libgcc'
              FCFLAGS_BIGENDIAN='-convert big_endian'
              FCFLAGS_FIXEDFORM='-fixed'
              FCFLAGS_MODDIR='-module $(abinit_moddir)'
              ;;
            *)
              abi_sys_spec_hnt="default"
              FCFLAGS_HINTS='-extend-source -vec-report0 -noaltparam -nofpscomp'
              FCFLAGS_PIC='-fPIC'
              FCFLAGS_FREEFORM='-free'
              FC_LDFLAGS_HINTS='-static-intel -static-libgcc'
              FCFLAGS_BIGENDIAN='-convert big_endian'
              FCFLAGS_FIXEDFORM='-fixed'
              FCFLAGS_MODDIR='-module $(abinit_moddir)'
              ;;
          esac   # [case: abi_sys_spec, indent: 4, item: True]
          ;;
      esac   # [case: abi_fc_version, indent: 2, item: True]
      ;;
    hitachi)
      abi_fc_vendor_hnt="hitachi"
      abi_fc_version_hnt="default"
      abi_sys_spec_hnt="default"
      FCFLAGS_HINTS='-hf95 -nosave -nohugeary'
      ;;
    g95)
      abi_fc_vendor_hnt="g95"
      abi_fc_version_hnt="default"
      abi_sys_spec_hnt="default"
      FCFLAGS_HINTS='-ffree-line-length-huge'
      FCFLAGS_64BITS='-m64'
      FCFLAGS_PIC='-fpic'
      FCFLAGS_FREEFORM='-ffree-form'
      FCFLAGS_32BITS='-m32'
      FCFLAGS_BIGENDIAN='-fconvert=big-endian -frecord-marker=4'
      FCFLAGS_FIXEDFORM='-ffixed-form'
      FCFLAGS_MODDIR='-fmod=$(abinit_moddir) -I$(abinit_moddir)'
      ;;
    fujitsu)
      abi_fc_vendor_hnt="fujitsu"
      abi_fc_version_hnt="default"
      abi_sys_spec_hnt="default"
      abi_fc_wrap='yes'
      FCFLAGS_HINTS='-Am -Ee -Ep'
      FCFLAGS_FIXEDFORM='-Fixed -X7'
      FCFLAGS_FREEFORM='-Free -X9'
      FCFLAGS_MODDIR='-M $(abinit_moddir)'
      ;;
    pathscale)
      abi_fc_vendor_hnt="pathscale"
      abi_fc_version_hnt="default"
      abi_sys_spec_hnt="default"
      FCFLAGS_HINTS='-extend-source'
      FCFLAGS_64BITS='-m64 -align64'
      FCFLAGS_PIC='-fPIC'
      FCFLAGS_FREEFORM='-freeform'
      FCFLAGS_32BITS='-m32 -align32'
      FCFLAGS_BIGENDIAN='-byteswapio'
      FCFLAGS_FIXEDFORM='-fixedform'
      FCFLAGS_MODDIR='-module $(abinit_moddir)'
      ;;
    nag)
      abi_fc_vendor_hnt="nag"
      abi_fc_version_hnt="default"
      abi_sys_spec_hnt="default"
      FCFLAGS_HINTS='-132 -english -strict95'
      FCFLAGS_64BITS='-64t'
      FCFLAGS_PIC='-PIC'
      FCFLAGS_FREEFORM='-free'
      fpp='-fpp'
      FCFLAGS_FIXEDFORM='-fixed'
      FCFLAGS_MODDIR='-mdir $(abinit_moddir) -I$(abinit_moddir)'
      ;;
    mipspro)
      abi_fc_vendor_hnt="mipspro"
      abi_fc_version_hnt="default"
      abi_sys_spec_hnt="default"
      FCFLAGS_HINTS='-extend_source -col120 -macro_expand'
      FCFLAGS_64BITS='-default64'
      FCFLAGS_PIC='-KPIC'
      FCFLAGS_FREEFORM='-freeform'
      FCFLAGS_32BITS='-32'
      abi_fc_wrap='yes'
      FCFLAGS_FIXEDFORM='-fixedform'
      ;;
    ibm)
      abi_fc_vendor_hnt="ibm"
      abi_fc_version_hnt="default"
      abi_sys_spec_hnt="default"
      FCFLAGS_HINTS='-qzerosize'
      FCFLAGS_64BITS='-q64'
      FCFLAGS_PIC='-qpic'
      FCFLAGS_FREEFORM='-qsuffix=cpp=F90:f=f90 -qfree=f90'
      FCFLAGS_32BITS='-q32'
      FCFLAGS_FIXEDFORM='-qsuffix=cpp=F:f=f -qfixed'
      FPPFLAGS_HINTS='-WF,-DHAVE_CONFIG_H'
      FCFLAGS_MODDIR='-qmoddir=$(abinit_moddir) -I$(abinit_moddir)'
      ;;
  esac   # [case: abi_fc_vendor, indent: 0, item: True]

  dnl Display settings
  AC_MSG_RESULT([${abi_fc_vendor_hnt}/${abi_fc_version_hnt}/${abi_sys_spec_hnt}])

]) #ABI_FC_HINTS


AC_DEFUN([ABI_CXX_HINTS],[
  dnl Init
  abi_cxx_vendor_hnt="none"
  abi_cxx_version_hnt="none"
  abi_sys_spec_hnt="none"

  dnl Look for hints flags
  AC_MSG_CHECKING([which cxx hints to apply])

  dnl Case built from config/hints/cxx_*.conf
  case "${abi_cxx_vendor}" in
    intel)
      abi_cxx_vendor_hnt="intel"
      abi_cxx_version_hnt="default"
      abi_sys_spec_hnt="default"
      CXXFLAGS_HINTS='-vec-report0'
      CXXFLAGS_PIC='-fPIC'
      CXX_LDFLAGS_HINTS='-static-libgcc -static-intel'
      ;;
    ibm)
      abi_cxx_vendor_hnt="ibm"
      abi_cxx_version_hnt="default"
      abi_sys_spec_hnt="default"
      CXXFLAGS_64BITS='-q64'
      CXXFLAGS_PIC='-qpic'
      ;;
    gnu)
      abi_cxx_vendor_hnt="gnu"
      abi_cxx_version_hnt="default"
      case "${abi_sys_spec}" in
        irix-mips)
          abi_sys_spec_hnt="irix-mips"
          CXXFLAGS_64BITS='-mabi=64'
          CXXFLAGS_PIC='-fPIC'
          CXXFLAGS_32BITS='-mabi=32'
          ;;
        *)
          abi_sys_spec_hnt="default"
          CXXFLAGS_64BITS='-m64'
          CXXFLAGS_PIC='-fPIC'
          CXXFLAGS_32BITS='-m32'
          ;;
      esac   # [case: abi_sys_spec, indent: 2, item: True]
      ;;
    pathscale)
      abi_cxx_vendor_hnt="pathscale"
      abi_cxx_version_hnt="default"
      abi_sys_spec_hnt="default"
      CXXFLAGS_64BITS='-m64'
      CXXFLAGS_PIC='-fPIC'
      CXXFLAGS_32BITS='-m32'
      ;;
    compaq)
      abi_cxx_vendor_hnt="compaq"
      abi_cxx_version_hnt="default"
      abi_sys_spec_hnt="default"
      big_endian='-convert big_endian'
      ;;
  esac   # [case: abi_cxx_vendor, indent: 0, item: True]

  dnl Display settings
  AC_MSG_RESULT([${abi_cxx_vendor_hnt}/${abi_cxx_version_hnt}/${abi_sys_spec_hnt}])

]) #ABI_CXX_HINTS


AC_DEFUN([ABI_AR_HINTS],[
  dnl Init
  abi_ar_vendor_hnt="none"
  abi_ar_version_hnt="none"
  abi_sys_spec_hnt="none"

  dnl Look for hints flags
  AC_MSG_CHECKING([which ar hints to apply])

  dnl Case built from config/hints/ar_*.conf
  if test "${abi_ar_vendor}" = "ibm"; then
    abi_ar_vendor_hnt="ibm"
    abi_ar_version_hnt="default"
    case "${abi_sys_spec}" in
      aix-*)
        abi_sys_spec_hnt="aix-*"
        ARFLAGS_64BITS='-X 64'
        ARFLAGS_32BITS='-X 32'
        ;;
      *)
        abi_sys_spec_hnt="default"
        
        ;;
    esac   # [case: abi_sys_spec, indent: 1, item: False]
  fi

  dnl Display settings
  AC_MSG_RESULT([${abi_ar_vendor_hnt}/${abi_ar_version_hnt}/${abi_sys_spec_hnt}])

]) #ABI_AR_HINTS


AC_DEFUN([ABI_FPP_HINTS],[
  dnl Init
  abi_fpp_vendor_hnt="none"
  abi_fpp_version_hnt="none"
  abi_sys_spec_hnt="none"

  dnl Look for hints flags
  AC_MSG_CHECKING([which fpp hints to apply])

  dnl Case built from config/hints/fpp_*.conf
  case "${abi_fpp_vendor}" in
    ibm)
      abi_fpp_vendor_hnt="ibm"
      abi_fpp_version_hnt="default"
      abi_sys_spec_hnt="default"
      FPPFLAGS_HINTS_EXT='-P -traditional-cpp'
      ;;
    *)
      abi_fpp_vendor_hnt="default"
      abi_fpp_version_hnt="default"
      abi_sys_spec_hnt="default"
      FPPFLAGS_HINTS_EXT='-P'
      ;;
  esac   # [case: abi_fpp_vendor, indent: 0, item: True]

  dnl Display settings
  AC_MSG_RESULT([${abi_fpp_vendor_hnt}/${abi_fpp_version_hnt}/${abi_sys_spec_hnt}])

]) #ABI_FPP_HINTS


AC_DEFUN([ABI_CPP_HINTS],[
  dnl Init
  abi_cpp_vendor_hnt="none"
  abi_cpp_version_hnt="none"
  abi_sys_spec_hnt="none"

  dnl Look for hints flags
  AC_MSG_CHECKING([which cpp hints to apply])

  dnl Case built from config/hints/cpp_*.conf
  case "${abi_cpp_vendor}" in
    ibm)
      abi_cpp_vendor_hnt="ibm"
      abi_cpp_version_hnt="default"
      case "${abi_sys_spec}" in
        aix-*)
          abi_sys_spec_hnt="aix-*"
          CPPFLAGS_HINTS_EXT='-P -traditional-cpp'
          ;;
        linux-*)
          abi_sys_spec_hnt="linux-*"
          CPPFLAGS_HINTS_EXT='-P -traditional-cpp'
          ;;
        *)
          abi_sys_spec_hnt="default"
          
          ;;
      esac   # [case: abi_sys_spec, indent: 2, item: True]
      ;;
    *)
      abi_cpp_vendor_hnt="default"
      abi_cpp_version_hnt="default"
      abi_sys_spec_hnt="default"
      CPPFLAGS_HINTS_EXT='-P -std=c99'
      ;;
  esac   # [case: abi_cpp_vendor, indent: 0, item: True]

  dnl Display settings
  AC_MSG_RESULT([${abi_cpp_vendor_hnt}/${abi_cpp_version_hnt}/${abi_sys_spec_hnt}])

]) #ABI_CPP_HINTS


AC_DEFUN([ABI_XPP_HINTS],[
  dnl Init
  abi_xpp_vendor_hnt="none"
  abi_xpp_version_hnt="none"
  abi_sys_spec_hnt="none"

  dnl Look for hints flags
  AC_MSG_CHECKING([which xpp hints to apply])

  dnl WARNING: no config files were found for language

  dnl Display settings
  AC_MSG_RESULT([${abi_xpp_vendor_hnt}/${abi_xpp_version_hnt}/${abi_sys_spec_hnt}])

]) #ABI_XPP_HINTS
