#!/bin/sh

#
# Functions for the build system of ABINIT (shell-scripts)
#

# Determine how to "echo" without newline: "echo -n" or "echo ...\c"
echon()
{
  if test "${ECHON}" = ""; then
    if test "`echo -n`" = "-n"; then
      ECHON="echo"; NNL="\c"
    else
      ECHON="echo -n"; NNL=""
    fi
  fi

  ${ECHON} "${*}${NNL}"
}
