

AC_DEFUN([ABI_CC_OPTFLAGS],[
  dnl Init
  abi_cc_vendor_opt="none"
  abi_cc_version_opt="none"
  abi_cpu_spec_opt="none"

  dnl Look for optimizations
  AC_MSG_CHECKING([which cc optimizations to apply])

  dnl Case built from config/optim/cc_*.conf
  case "${abi_cc_vendor}" in
    intel)
      abi_cc_vendor_opt="intel"
      case "${abi_cc_version}" in
        9.0)
          abi_cc_version_opt="9.0"
          case "${abi_cpu_spec}" in
            intel_pentium4)
              abi_cpu_spec_opt="intel_pentium4"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium3)
              abi_cpu_spec_opt="intel_pentium3"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=pentiumpro -msse -xK"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=pentiumpro -msse -xK"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=pentiumpro -msse -xK"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=itanium2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=itanium2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=itanium2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=itanium"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=itanium"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=itanium"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
        9.1)
          abi_cc_version_opt="9.1"
          case "${abi_cpu_spec}" in
            intel_pentium4)
              abi_cpu_spec_opt="intel_pentium4"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium3)
              abi_cpu_spec_opt="intel_pentium3"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=pentiumpro -msse -xK"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=pentiumpro -msse -xK"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=pentiumpro -msse -xK"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=itanium2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=itanium2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=itanium2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_xeon)
              abi_cpu_spec_opt="intel_xeon"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -xT"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -xW"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -xP"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=itanium"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=itanium"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=itanium"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_coreduo)
              abi_cpu_spec_opt="intel_coreduo"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -xP"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xP"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xP"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_core2)
              abi_cpu_spec_opt="intel_core2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -msse3 -xT"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -msse3 -xT"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -msse3 -xT"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
        10.0)
          abi_cc_version_opt="10.0"
          case "${abi_cpu_spec}" in
            intel_pentium4)
              abi_cpu_spec_opt="intel_pentium4"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium3)
              abi_cpu_spec_opt="intel_pentium3"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=pentiumpro -msse -xK"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=pentiumpro -msse -xK"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=pentiumpro -msse -xK"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=itanium2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=itanium2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=itanium2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=itanium"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=itanium"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=itanium"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_coreduo)
              abi_cpu_spec_opt="intel_coreduo"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -xP"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xP"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xP"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_core2)
              abi_cpu_spec_opt="intel_core2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -msse3 -xT"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -msse3 -xT"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -msse3 -xT"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
        *)
          abi_cc_version_opt="default"
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cc_version, indent: 2, item: True]
      ;;
    ibm)
      abi_cc_vendor_opt="ibm"
      abi_cc_version_opt="default"
      case "${abi_cpu_spec}" in
        ibm_powerpc)
          abi_cpu_spec_opt="ibm_powerpc"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O4 -qarch=auto -qtune=auto -qstrict -qspill=2000"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -qarch=auto -qtune=auto -qstrict -qspill=2000"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O3 -qarch=auto -qtune=auto -qstrict -qspill=2000"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        ibm_powerpc64)
          abi_cpu_spec_opt="ibm_powerpc64"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O4 -qarch=auto -qtune=auto -qstrict -qspill=2000"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -qarch=auto -qtune=auto -qstrict -qspill=2000"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O3 -qarch=auto -qtune=auto -qstrict -qspill=2000"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    open64)
      abi_cc_vendor_opt="open64"
      abi_cc_version_opt="default"
      case "${abi_cpu_spec}" in
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3 -march=pentium4 -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3 -march=opteron -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    sun)
      abi_cc_vendor_opt="sun"
      abi_cc_version_opt="default"
      case "${abi_cpu_spec}" in
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=opteron -xarch=sse2a -m64 -xlibmil -xlibmopt -xvector=simd -xipo -xjobs=3"
              CC_LDFLAGS_OPTIM="-xvector=simd -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=opteron -xarch=sse2a -m64"
              CC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=opteron -xarch=sse2a -m64 -xvector=simd"
              CC_LDFLAGS_OPTIM="-xvector=simd"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xlibmil -xlibmopt -xvector=lib -xipo"
              CC_LDFLAGS_OPTIM="-xvector=lib -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=pentium4 -xarch=sse2"
              CC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=pentium4 -xarch=sse2"
              CC_LDFLAGS_OPTIM=""
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_pentium3)
          abi_cpu_spec_opt="intel_pentium3"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=pentium3 -xarch=sse -xlibmil -xlibmopt -xvector=lib -xipo"
              CC_LDFLAGS_OPTIM="-xvector=lib -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=pentium3 -xarch=sse"
              CC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=pentium3 -xarch=sse"
              CC_LDFLAGS_OPTIM=""
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_xeon)
          abi_cpu_spec_opt="intel_xeon"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-xO5 -fround=nearest -xarch=sse2 -xchip=native -xcache=native -xvector=lib"
              CC_LDFLAGS_OPTIM="-xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-xO2 -fround=nearest"
              CC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-xO3 -fround=nearest -xarch=sse2 -xchip=native -xcache=native -xvector=lib"
              CC_LDFLAGS_OPTIM=""
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_coreduo)
          abi_cpu_spec_opt="intel_coreduo"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xlibmil -xlibmopt -xvector=simd -xipo -xjobs=3"
              CC_LDFLAGS_OPTIM="-xvector=simd -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=pentium4 -xarch=sse2"
              CC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xvector=simd"
              CC_LDFLAGS_OPTIM="-xvector=simd"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_core2)
          abi_cpu_spec_opt="intel_core2"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xlibmil -xlibmopt -xvector=simd -xipo -xjobs=3"
              CC_LDFLAGS_OPTIM="-xvector=simd -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=pentium4 -xarch=sse2"
              CC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xvector=simd"
              CC_LDFLAGS_OPTIM="-xvector=simd"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=native -xlibmil -xlibmopt -xipo -xjobs=3"
              CC_LDFLAGS_OPTIM="-xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=native"
              CC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=native"
              CC_LDFLAGS_OPTIM=""
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    gnu)
      abi_cc_vendor_opt="gnu"
      case "${abi_cc_version}" in
        4.1)
          abi_cc_version_opt="4.1"
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        4.2)
          abi_cc_version_opt="4.2"
          case "${abi_cpu_spec}" in
            amd_athlon64)
              abi_cpu_spec_opt="amd_athlon64"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -march=athlon64"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -march=athlon64"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -march=athlon64"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            amd_opteron)
              abi_cpu_spec_opt="amd_opteron"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -march=opteron"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -march=opteron"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -march=opteron"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium4)
              abi_cpu_spec_opt="intel_pentium4"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -fforce-addr -march=pentium4 -mmmx -msse -msse2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -fforce-addr -march=pentium4 -mmmx -msse -msse2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -fforce-addr -march=pentium4 -mmmx -msse -msse2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            dec_alphaev67)
              abi_cpu_spec_opt="dec_alphaev67"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            amd_athlon)
              abi_cpu_spec_opt="amd_athlon"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -march=athlon"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -march=athlon"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -march=athlon"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            ibm_powerpc)
              abi_cpu_spec_opt="ibm_powerpc"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O4 -mpowerpc"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mpowerpc"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mpowerpc"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            dec_alphaev56)
              abi_cpu_spec_opt="dec_alphaev56"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -fforce-addr"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -fforce-addr"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -fforce-addr"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_xeon)
              abi_cpu_spec_opt="intel_xeon"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -fschedule-insns2 -march=nocona -mmmx -msse -msse2 -msse3 -mfpmath=sse"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -march=nocona -mmmx -msse -mfpmath=sse"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -fschedule-insns2 -march=nocona -mmmx -msse -msse2 -msse3 -mfpmath=sse"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -fforce-addr"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -fforce-addr"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -fforce-addr"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_coreduo)
              abi_cpu_spec_opt="intel_coreduo"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -march=prescott -mmmx -msse -msse2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium3)
              abi_cpu_spec_opt="intel_pentium3"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -fforce-addr -march=pentium3 -mmmx -msse"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -fforce-addr -march=pentium3 -mmmx -msse"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -fforce-addr -march=pentium3 -mmmx -msse"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            ibm_powerpc64)
              abi_cpu_spec_opt="ibm_powerpc64"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O4 -mpowerpc64"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -mpowerpc64"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -mpowerpc64"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_core2)
              abi_cpu_spec_opt="intel_core2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3 -march=prescott -mmmx -msse -msse2 -msse3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2 -msse3"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2 -msse3"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
        *)
          abi_cc_version_opt="default"
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3 -mtune=native -march=native"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2 -mtune=native -march=native -mfpmath=sse"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cc_version, indent: 2, item: True]
      ;;
    pathscale)
      abi_cc_vendor_opt="pathscale"
      abi_cc_version_opt="default"
      case "${abi_cpu_spec}" in
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3 -march=pentium4 -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3 -march=opteron -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    compaq)
      abi_cc_vendor_opt="compaq"
      abi_cc_version_opt="default"
      case "${abi_cpu_spec}" in
        dec_alphaev67)
          abi_cpu_spec_opt="dec_alphaev67"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3 -arch host -tune host"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -arch host -tune host"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O3 -arch host -tune host"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        dec_alphaev56)
          abi_cpu_spec_opt="dec_alphaev56"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3 -arch host -tune host"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -arch host -tune host"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O3 -arch host -tune host"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
  esac   # [case: abi_cc_vendor, indent: 0, item: True]

  dnl Display settings
  AC_MSG_RESULT([${abi_cc_vendor_opt}/${abi_cc_version_opt}/${abi_cpu_spec_opt}])

]) #ABI_CC_OPTFLAGS


AC_DEFUN([ABI_CXX_OPTFLAGS],[
  dnl Init
  abi_cxx_vendor_opt="none"
  abi_cxx_version_opt="none"
  abi_cpu_spec_opt="none"

  dnl Look for optimizations
  AC_MSG_CHECKING([which cxx optimizations to apply])

  dnl Case built from config/optim/cxx_*.conf
  case "${abi_cxx_vendor}" in
    intel)
      abi_cxx_vendor_opt="intel"
      case "${abi_cxx_version}" in
        9.0)
          abi_cxx_version_opt="9.0"
          case "${abi_cpu_spec}" in
            intel_pentium4)
              abi_cpu_spec_opt="intel_pentium4"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium3)
              abi_cpu_spec_opt="intel_pentium3"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=pentiumpro -msse -xK"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentiumpro -msse -xK"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentiumpro -msse -xK"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=itanium2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=itanium2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=itanium2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=itanium"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=itanium"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=itanium"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
        9.1)
          abi_cxx_version_opt="9.1"
          case "${abi_cpu_spec}" in
            intel_pentium4)
              abi_cpu_spec_opt="intel_pentium4"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium3)
              abi_cpu_spec_opt="intel_pentium3"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=pentiumpro -msse -xK"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentiumpro -msse -xK"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentiumpro -msse -xK"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=itanium2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=itanium2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=itanium2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_xeon)
              abi_cpu_spec_opt="intel_xeon"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -xT"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -xW"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -xP"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=itanium"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=itanium"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=itanium"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_coreduo)
              abi_cpu_spec_opt="intel_coreduo"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -xP"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xP"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xP"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_core2)
              abi_cpu_spec_opt="intel_core2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -msse3 -xT"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -msse3 -xT"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -msse3 -xT"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
        10.0)
          abi_cxx_version_opt="10.0"
          case "${abi_cpu_spec}" in
            intel_pentium4)
              abi_cpu_spec_opt="intel_pentium4"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xN"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium3)
              abi_cpu_spec_opt="intel_pentium3"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=pentiumpro -msse -xK"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentiumpro -msse -xK"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentiumpro -msse -xK"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=itanium2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=itanium2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=itanium2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=itanium"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=itanium"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=itanium"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_coreduo)
              abi_cpu_spec_opt="intel_coreduo"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -xP"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xP"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -xP"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_core2)
              abi_cpu_spec_opt="intel_core2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -mcpu=pentium4 -msse -msse2 -msse3 -xT"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -msse3 -xT"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mcpu=pentium4 -msse -msse2 -msse3 -xT"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
        *)
          abi_cxx_version_opt="default"
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cxx_version, indent: 2, item: True]
      ;;
    ibm)
      abi_cxx_vendor_opt="ibm"
      abi_cxx_version_opt="default"
      case "${abi_cpu_spec}" in
        ibm_powerpc)
          abi_cpu_spec_opt="ibm_powerpc"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O4 -qarch=auto -qtune=auto -qstrict -qspill=2000 -qessl"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2 -qarch=auto -qtune=auto -qstrict -qspill=2000 -qessl"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O3 -qarch=auto -qtune=auto -qstrict -qspill=2000 -qessl"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        ibm_powerpc64)
          abi_cpu_spec_opt="ibm_powerpc64"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O4 -qarch=auto -qtune=auto -qstrict -qspill=2000 -qessl"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2 -qarch=auto -qtune=auto -qstrict -qspill=2000 -qessl"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O3 -qarch=auto -qtune=auto -qstrict -qspill=2000 -qessl"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    sun)
      abi_cxx_vendor_opt="sun"
      abi_cxx_version_opt="default"
      case "${abi_cpu_spec}" in
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=opteron -xarch=sse2a -m64 -xlibmil -xlibmopt -xvector=simd -xipo -xjobs=3"
              CXX_LDFLAGS_OPTIM="-xvector=simd -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=opteron -xarch=sse2a -m64"
              CXX_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=opteron -xarch=sse2a -m64 -xvector=simd"
              CXX_LDFLAGS_OPTIM="-xvector=simd"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xlibmil -xlibmopt -xvector=lib -xipo"
              CXX_LDFLAGS_OPTIM="-xvector=lib -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=pentium4 -xarch=sse2"
              CXX_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=pentium4 -xarch=sse2"
              CXX_LDFLAGS_OPTIM=""
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_pentium3)
          abi_cpu_spec_opt="intel_pentium3"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=pentium3 -xarch=sse -xlibmil -xlibmopt -xvector=lib -xipo"
              CXX_LDFLAGS_OPTIM="-xvector=lib -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=pentium3 -xarch=sse"
              CXX_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=pentium3 -xarch=sse"
              CXX_LDFLAGS_OPTIM=""
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_xeon)
          abi_cpu_spec_opt="intel_xeon"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-xO5 -fround=nearest -xarch=sse2 -xchip=native -xcache=native -xvector=lib"
              CXX_LDFLAGS_OPTIM="-xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-xO2 -fround=nearest"
              CXX_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-xO3 -fround=nearest -xarch=sse2 -xchip=native -xcache=native -xvector=lib"
              CXX_LDFLAGS_OPTIM=""
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_coreduo)
          abi_cpu_spec_opt="intel_coreduo"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xlibmil -xlibmopt -xvector=simd -xipo -xjobs=3"
              CXX_LDFLAGS_OPTIM="-xvector=simd -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=pentium4 -xarch=sse2"
              CXX_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xvector=simd"
              CXX_LDFLAGS_OPTIM="-xvector=simd"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_core2)
          abi_cpu_spec_opt="intel_core2"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xlibmil -xlibmopt -xvector=simd -xipo -xjobs=3"
              CXX_LDFLAGS_OPTIM="-xvector=simd -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=pentium4 -xarch=sse2"
              CXX_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xvector=simd"
              CXX_LDFLAGS_OPTIM="-xvector=simd"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=native -xlibmil -xlibmopt -xipo -xjobs=3"
              CXX_LDFLAGS_OPTIM="-xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=native"
              CXX_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=native"
              CXX_LDFLAGS_OPTIM=""
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    gnu)
      abi_cxx_vendor_opt="gnu"
      case "${abi_cxx_version}" in
        4.1)
          abi_cxx_version_opt="4.1"
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        4.2)
          abi_cxx_version_opt="4.2"
          case "${abi_cpu_spec}" in
            amd_athlon64)
              abi_cpu_spec_opt="amd_athlon64"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -march=athlon64"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -march=athlon64"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -march=athlon64"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            amd_opteron)
              abi_cpu_spec_opt="amd_opteron"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -march=opteron"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -march=opteron"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -march=opteron"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium4)
              abi_cpu_spec_opt="intel_pentium4"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -fforce-addr -march=pentium4 -mmmx -msse -msse2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -fforce-addr -march=pentium4 -mmmx -msse -msse2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -fforce-addr -march=pentium4 -mmmx -msse -msse2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            dec_alphaev67)
              abi_cpu_spec_opt="dec_alphaev67"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            amd_athlon)
              abi_cpu_spec_opt="amd_athlon"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -march=athlon"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -march=athlon"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -march=athlon"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            ibm_powerpc)
              abi_cpu_spec_opt="ibm_powerpc"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O4 -mpowerpc"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mpowerpc"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mpowerpc"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            dec_alphaev56)
              abi_cpu_spec_opt="dec_alphaev56"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -fforce-addr"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -fforce-addr"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -fforce-addr"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_xeon)
              abi_cpu_spec_opt="intel_xeon"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -fschedule-insns2 -march=nocona -mmmx -msse -msse2 -msse3 -mfpmath=sse"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -march=nocona -mmmx -msse -mfpmath=sse"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -fschedule-insns2 -march=nocona -mmmx -msse -msse2 -msse3 -mfpmath=sse"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -fforce-addr"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -fforce-addr"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -fforce-addr"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_coreduo)
              abi_cpu_spec_opt="intel_coreduo"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -march=prescott -mmmx -msse -msse2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium3)
              abi_cpu_spec_opt="intel_pentium3"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -fforce-addr -march=pentium3 -mmmx -msse"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -fforce-addr -march=pentium3 -mmmx -msse"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -fforce-addr -march=pentium3 -mmmx -msse"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            ibm_powerpc64)
              abi_cpu_spec_opt="ibm_powerpc64"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O4 -mpowerpc64"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -mpowerpc64"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -mpowerpc64"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_core2)
              abi_cpu_spec_opt="intel_core2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3 -march=prescott -mmmx -msse -msse2 -msse3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2 -msse3"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2 -msse3"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  CXXFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  CXXFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  CXXFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
        *)
          abi_cxx_version_opt="default"
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3 -mtune=native -march=native"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O2 -mtune=native -march=native -mfpmath=sse"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cxx_version, indent: 2, item: True]
      ;;
    pathscale)
      abi_cxx_vendor_opt="pathscale"
      abi_cxx_version_opt="default"
      case "${abi_cpu_spec}" in
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3 -march=pentium4 -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3 -march=opteron -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    compaq)
      abi_cxx_vendor_opt="compaq"
      abi_cxx_version_opt="default"
      case "${abi_cpu_spec}" in
        dec_alphaev67)
          abi_cpu_spec_opt="dec_alphaev67"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3 -arch host -tune host"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2 -arch host -tune host"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O3 -arch host -tune host"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        dec_alphaev56)
          abi_cpu_spec_opt="dec_alphaev56"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3 -arch host -tune host"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2 -arch host -tune host"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O3 -arch host -tune host"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
  esac   # [case: abi_cxx_vendor, indent: 0, item: True]

  dnl Display settings
  AC_MSG_RESULT([${abi_cxx_vendor_opt}/${abi_cxx_version_opt}/${abi_cpu_spec_opt}])

]) #ABI_CXX_OPTFLAGS


AC_DEFUN([ABI_FC_OPTFLAGS],[
  dnl Init
  abi_fc_vendor_opt="none"
  abi_fc_version_opt="none"
  abi_cpu_spec_opt="none"

  dnl Look for optimizations
  AC_MSG_CHECKING([which fc optimizations to apply])

  dnl Case built from config/optim/fc_*.conf
  case "${abi_fc_vendor}" in
    compaq)
      abi_fc_vendor_opt="compaq"
      abi_fc_version_opt="default"
      case "${abi_cpu_spec}" in
        dec_alphaev67)
          abi_cpu_spec_opt="dec_alphaev67"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -arch host -tune host"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -arch host -tune host"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3 -arch host -tune host"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        dec_alphaev56)
          abi_cpu_spec_opt="dec_alphaev56"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -arch host -tune host"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -arch host -tune host"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3 -arch host -tune host"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -arch host -tune host"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -arch host -tune host"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -arch host -tune host"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    gnu)
      abi_fc_vendor_opt="gnu"
      case "${abi_fc_version}" in
        4.5)
          abi_fc_version_opt="4.5"
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -mtune=native -march=native -funroll-loops"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -mtune=native -march=native -mfpmath=sse"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        4.2)
          abi_fc_version_opt="4.2"
          case "${abi_cpu_spec}" in
            amd_athlon64)
              abi_cpu_spec_opt="amd_athlon64"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -march=athlon64"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -march=athlon64"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -march=athlon64"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            amd_opteron)
              abi_cpu_spec_opt="amd_opteron"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -march=opteron"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -march=opteron"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -march=opteron"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium4)
              abi_cpu_spec_opt="intel_pentium4"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -march=pentium4 -mmmx -msse -msse2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -march=pentium4 -mmmx -msse -msse2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -march=pentium4 -mmmx -msse -msse2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            dec_alphaev67)
              abi_cpu_spec_opt="dec_alphaev67"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            amd_athlon)
              abi_cpu_spec_opt="amd_athlon"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -march=athlon"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -march=athlon"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -march=athlon"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            ibm_powerpc)
              abi_cpu_spec_opt="ibm_powerpc"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O4 -mpowerpc"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -mpowerpc"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -mpowerpc"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            dec_alphaev56)
              abi_cpu_spec_opt="dec_alphaev56"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_xeon)
              abi_cpu_spec_opt="intel_xeon"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -fschedule-insns2 -march=nocona -mmmx -msse -msse2 -msse3 -mfpmath=sse"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -march=nocona -mmmx -msse -mfpmath=sse"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -fschedule-insns2 -march=nocona -mmmx -msse -msse2 -msse3 -mfpmath=sse"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_coreduo)
              abi_cpu_spec_opt="intel_coreduo"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -march=prescott -mmmx -msse -msse2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium3)
              abi_cpu_spec_opt="intel_pentium3"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -march=pentium3 -mmmx -msse"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -march=pentium3 -mmmx -msse"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -march=pentium3 -mmmx -msse"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            ibm_powerpc64)
              abi_cpu_spec_opt="ibm_powerpc64"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O4 -mpowerpc64"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -mpowerpc64"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -mpowerpc64"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_core2)
              abi_cpu_spec_opt="intel_core2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -march=prescott -mmmx -msse -msse2 -msse3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2 -msse3"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2 -msse3"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
        *)
          abi_fc_version_opt="default"
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -mtune=native -march=native"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -mtune=native -march=native -mfpmath=sse"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_fc_version, indent: 2, item: True]
      ;;
    open64)
      abi_fc_vendor_opt="open64"
      abi_fc_version_opt="default"
      case "${abi_cpu_spec}" in
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=pentium4 -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=opteron -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -OPT:Olimit=0 -g -ggdb"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    sun)
      abi_fc_vendor_opt="sun"
      abi_fc_version_opt="default"
      case "${abi_cpu_spec}" in
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=opteron -xarch=sse2a -m64 -xlibmil -xlibmopt -xvector=simd -xipo -xjobs=3"
              FC_LDFLAGS_OPTIM="-xvector=simd -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=opteron -xarch=sse2a -m64"
              FC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=opteron -xarch=sse2a -m64 -xvector=simd"
              FC_LDFLAGS_OPTIM="-xvector=simd"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xlibmil -xlibmopt -xvector=lib -xipo"
              FC_LDFLAGS_OPTIM="-xvector=lib -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=pentium4 -xarch=sse2"
              FC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=pentium4 -xarch=sse2"
              FC_LDFLAGS_OPTIM=""
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_pentium3)
          abi_cpu_spec_opt="intel_pentium3"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=pentium3 -xarch=sse -xlibmil -xlibmopt -xvector=lib -xipo"
              FC_LDFLAGS_OPTIM="-xvector=lib -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=pentium3 -xarch=sse"
              FC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=pentium3 -xarch=sse"
              FC_LDFLAGS_OPTIM=""
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_xeon)
          abi_cpu_spec_opt="intel_xeon"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-xO5 -fround=nearest -xarch=sse2 -xchip=native -xcache=native -xvector=lib"
              FC_LDFLAGS_OPTIM="-xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-xO2 -fround=nearest"
              FC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-xO3 -fround=nearest -xarch=sse2 -xchip=native -xcache=native -xvector=lib"
              FC_LDFLAGS_OPTIM=""
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_coreduo)
          abi_cpu_spec_opt="intel_coreduo"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xlibmil -xlibmopt -xvector=simd -xipo -xjobs=3"
              FC_LDFLAGS_OPTIM="-xvector=simd -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=pentium4 -xarch=sse2"
              FC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xvector=simd"
              FC_LDFLAGS_OPTIM="-xvector=simd"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_core2)
          abi_cpu_spec_opt="intel_core2"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xlibmil -xlibmopt -xvector=simd -xipo -xjobs=3"
              FC_LDFLAGS_OPTIM="-xvector=simd -xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=pentium4 -xarch=sse2"
              FC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=pentium4 -xarch=sse2 -xvector=simd"
              FC_LDFLAGS_OPTIM="-xvector=simd"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-xO4 -fround=nearest -xtarget=native -xlibmil -xlibmopt -xipo -xjobs=3"
              FC_LDFLAGS_OPTIM="-xipo"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-xO2 -fround=nearest -xtarget=native"
              FC_LDFLAGS_OPTIM=""
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-xO3 -fround=nearest -xtarget=native"
              FC_LDFLAGS_OPTIM=""
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    intel)
      abi_fc_vendor_opt="intel"
      case "${abi_fc_version}" in
        10.0)
          abi_fc_version_opt="10.0"
          case "${abi_cpu_spec}" in
            intel_centrino)
              abi_cpu_spec_opt="intel_centrino"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn4 -tune pn4 -xN"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xN"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xN"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            amd_opteron)
              abi_cpu_spec_opt="amd_opteron"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium4)
              abi_cpu_spec_opt="intel_pentium4"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn4 -tune pn4 -xN"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xN"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xN"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium3)
              abi_cpu_spec_opt="intel_pentium3"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn3 -tune pn3 -xK"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn3 -tune pn3 -xK"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn3 -tune pn3 -xK"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -tpp2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -tpp2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -tpp2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_xeon)
              abi_cpu_spec_opt="intel_xeon"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn4 -tune pn4 -xT"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xP"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xT"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -tpp1"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -tpp1"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O3"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_coreduo)
              abi_cpu_spec_opt="intel_coreduo"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn4 -tune pn4 -xP"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xP"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xP"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_core2)
              abi_cpu_spec_opt="intel_core2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn4 -tune pn4 -xT"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xT"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xT"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
        10.1)
          abi_fc_version_opt="10.1"
          case "${abi_cpu_spec}" in
            intel_centrino)
              abi_cpu_spec_opt="intel_centrino"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn4 -tune pn4 -xN"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xN"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xN"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            amd_opteron)
              abi_cpu_spec_opt="amd_opteron"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium4)
              abi_cpu_spec_opt="intel_pentium4"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn4 -tune pn4 -xN"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xN"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xN"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium3)
              abi_cpu_spec_opt="intel_pentium3"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn3 -tune pn3 -xK"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn3 -tune pn3 -xK"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn3 -tune pn3 -xK"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -tpp2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -tpp2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_xeon)
              abi_cpu_spec_opt="intel_xeon"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn4 -tune pn4 -ip -mcmodel=large -xT"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xT"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -ip -mcmodel=large -xT"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -tpp1"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -tpp1"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_coreduo)
              abi_cpu_spec_opt="intel_coreduo"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn4 -tune pn4 -xP"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xP"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xP"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_core2)
              abi_cpu_spec_opt="intel_core2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn4 -tune pn4 -xT"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xT"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xT"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
        9.0)
          abi_fc_version_opt="9.0"
          case "${abi_cpu_spec}" in
            amd_opteron)
              abi_cpu_spec_opt="amd_opteron"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium4)
              abi_cpu_spec_opt="intel_pentium4"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn4 -tune pn4 -tpp7 -xN"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -tpp7 -xN"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -tpp7 -xN"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium3)
              abi_cpu_spec_opt="intel_pentium3"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn3 -tune pn3 -tpp6 -xK"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn3 -tune pn3 -tpp6 -xK"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn3 -tune pn3 -tpp6 -xK"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -tpp2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -tpp2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -tpp2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -tpp1"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -tpp1"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -tpp1"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
        9.1)
          abi_fc_version_opt="9.1"
          case "${abi_cpu_spec}" in
            amd_opteron)
              abi_cpu_spec_opt="amd_opteron"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium4)
              abi_cpu_spec_opt="intel_pentium4"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn4 -tune pn4 -tpp7 -xN"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -tpp7 -xN"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -tpp7 -xN"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_pentium3)
              abi_cpu_spec_opt="intel_pentium3"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn3 -tune pn3 -tpp6 -xK"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn3 -tune pn3 -tpp6 -xK"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn3 -tune pn3 -tpp6 -xK"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -tpp2"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -tpp2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -tpp2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_xeon)
              abi_cpu_spec_opt="intel_xeon"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -arch pn4 -tune pn4 -xT"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xW"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -arch pn4 -tune pn4 -xT"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -tpp1"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -tpp1"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -tpp1"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_coreduo)
              abi_cpu_spec_opt="intel_coreduo"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -tpp7 -xP"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -tpp7 -xP"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -tpp7 -xP"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_core2)
              abi_cpu_spec_opt="intel_core2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -tpp7 -xT"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -tpp7 -xT"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -tpp7 -xT"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
        *)
          abi_fc_version_opt="default"
          case "${abi_cpu_spec}" in
            intel_itanium1)
              abi_cpu_spec_opt="intel_itanium1"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -fp-model fast=1 -fp-relaxed -ip"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -fp-model precise -fp-speculation=safe"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            intel_itanium2)
              abi_cpu_spec_opt="intel_itanium2"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -fp-model fast=1 -fp-relaxed -ip"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -fp-model precise -fp-speculation=safe"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            amd_opteron)
              abi_cpu_spec_opt="amd_opteron"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O1"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
            *)
              abi_cpu_spec_opt="default"
              case "${enable_optim}" in
                aggressive)
                  enable_optim_opt="aggressive"
                  FCFLAGS_OPTIM="-O3 -xHOST"
                  ;;
                safe)
                  enable_optim_opt="safe"
                  FCFLAGS_OPTIM="-O2 -xHost -fltconsistency -fp-model precise -fp-speculation=safe -prec-div -prec-sqrt"
                  ;;
                standard)
                  enable_optim_opt="standard"
                  FCFLAGS_OPTIM="-O2 -xHost"
                  ;;
              esac   # [case: enable_optim, indent: 6, item: True]
              ;;
          esac   # [case: abi_cpu_spec, indent: 4, item: True]
          ;;
      esac   # [case: abi_fc_version, indent: 2, item: True]
      ;;
    g95)
      abi_fc_vendor_opt="g95"
      abi_fc_version_opt="default"
      case "${abi_cpu_spec}" in
        amd_athlon64)
          abi_cpu_spec_opt="amd_athlon64"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=athlon64"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=athlon64"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=athlon64"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=opteron"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=opteron"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=opteron"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=pentium4 -mmmx -msse -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=pentium4 -mmmx -msse -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=pentium4 -mmmx -msse -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        dec_alphaev67)
          abi_cpu_spec_opt="dec_alphaev67"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_athlon)
          abi_cpu_spec_opt="amd_athlon"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=athlon"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=athlon"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=athlon"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        ibm_powerpc)
          abi_cpu_spec_opt="ibm_powerpc"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O4 -mpowerpc"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -mpowerpc"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -mpowerpc"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        dec_alphaev56)
          abi_cpu_spec_opt="dec_alphaev56"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_itanium2)
          abi_cpu_spec_opt="intel_itanium2"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_xeon)
          abi_cpu_spec_opt="intel_xeon"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=nocona -mmmx -msse -msse2 -mfpmath=sse"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=nocona -mmmx -msse -mfpmath=sse"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3 -march=nocona -mmmx -msse -msse2 -mfpmath=sse"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_itanium1)
          abi_cpu_spec_opt="intel_itanium1"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_coreduo)
          abi_cpu_spec_opt="intel_coreduo"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=prescott -mmmx -msse -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_pentium3)
          abi_cpu_spec_opt="intel_pentium3"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=pentium3 -mmmx -msse"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=pentium3 -mmmx -msse"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=pentium3 -mmmx -msse"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        ibm_powerpc64)
          abi_cpu_spec_opt="ibm_powerpc64"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O4 -mpowerpc64"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -mpowerpc64"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -mpowerpc64"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_core2)
          abi_cpu_spec_opt="intel_core2"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=prescott -mmmx -msse -msse2 -msse3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2 -msse3"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2 -msse3"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    fujitsu)
      abi_fc_vendor_opt="fujitsu"
      abi_fc_version_opt="default"
      abi_cpu_spec_opt="default"
      case "${enable_optim}" in
        aggressive)
          enable_optim_opt="aggressive"
          FCFLAGS_OPTIM="-Of -X9 -Ps -Wv,-md"
          ;;
        safe)
          enable_optim_opt="safe"
          FCFLAGS_OPTIM="-Of -X9 -Ps -Wv,-md"
          ;;
        standard)
          enable_optim_opt="standard"
          FCFLAGS_OPTIM="-Of -X9 -Ps -Wv,-md"
          ;;
      esac   # [case: enable_optim, indent: 2, item: True]
      ;;
    pathscale)
      abi_fc_vendor_opt="pathscale"
      abi_fc_version_opt="default"
      case "${abi_cpu_spec}" in
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=pentium4 -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=opteron -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    nag)
      abi_fc_vendor_opt="nag"
      abi_fc_version_opt="default"
      abi_cpu_spec_opt="default"
      case "${enable_optim}" in
        aggressive)
          enable_optim_opt="aggressive"
          FCFLAGS_OPTIM="-O4"
          ;;
        safe)
          enable_optim_opt="safe"
          FCFLAGS_OPTIM="-O2"
          ;;
        standard)
          enable_optim_opt="standard"
          FCFLAGS_OPTIM="-O3"
          ;;
      esac   # [case: enable_optim, indent: 2, item: True]
      ;;
    mipspro)
      abi_fc_vendor_opt="mipspro"
      abi_fc_version_opt="default"
      case "${abi_cpu_spec}" in
        sgi_mips4)
          abi_cpu_spec_opt="sgi_mips4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -OPT:Olimit=7168 -mips4"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -OPT:Olimit=7168 -mips4"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3 -OPT:Olimit=7168 -mips4"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        sgi_mips3)
          abi_cpu_spec_opt="sgi_mips3"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -OPT:Olimit=7168 -mips3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -OPT:Olimit=7168 -mips3"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3 -OPT:Olimit=7168 -mips3"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -OPT:Olimit=7168"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -OPT:Olimit=7168"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3 -OPT:Olimit=7168"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    ibm)
      abi_fc_vendor_opt="ibm"
      abi_fc_version_opt="default"
      case "${abi_cpu_spec}" in
        ibm_powerpc)
          abi_cpu_spec_opt="ibm_powerpc"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O4 -qmaxmem=65536 -qspill=2000 -qarch=auto -qtune=auto -qcache=auto -qstrict -qsuppress=1520-022:1520-031"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -qmaxmem=65536 -qspill=2000 -qarch=auto -qtune=auto -qcache=auto -qstrict -qsuppress=1520-022:1520-031"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3 -qmaxmem=65536 -qspill=2000 -qarch=auto -qtune=auto -qcache=auto -qstrict -qsuppress=1520-022:1520-031"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        ibm_powerpc64)
          abi_cpu_spec_opt="ibm_powerpc64"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O4 -qmaxmem=65536 -qspill=2000 -qarch=auto -qtune=auto -qcache=auto -qstrict -qsuppress=1520-022:1520-031"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -qmaxmem=65536 -qspill=2000 -qarch=auto -qtune=auto -qcache=auto -qstrict -qsuppress=1520-022:1520-031"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3 -qmaxmem=65536 -qspill=2000 -qarch=auto -qtune=auto -qcache=auto -qstrict -qsuppress=1520-022:1520-031"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
  esac   # [case: abi_fc_vendor, indent: 0, item: True]

  dnl Display settings
  AC_MSG_RESULT([${abi_fc_vendor_opt}/${abi_fc_version_opt}/${abi_cpu_spec_opt}])

]) #ABI_FC_OPTFLAGS
