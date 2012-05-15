#!/bin/sh
#
# Copyright (C) 2011 Yann Pouillon
#
# This file is part of the Abinit GUI software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

# Stop at first error encountered
set -e

# Check that we are in the right directory
if test ! -s "./configure.ac" -o ! -s "src/Main.java"; then
  echo "[guibuild] This is not a Abinit GUI source tree - aborting now" >&2
  exit 1
fi

# Create possibly missing directories
mkdir -p gnu

# Generate M4 includes
echo "[guibuild] Generating aclocal.m4"
aclocal -I m4

# Generate configure
echo "[guibuild] Generating configure script"
autoconf

# Generate makefile inputs
# Do not use "automake --force-missing", as it overwrites the INSTALL file.
echo "[guibuild] Generating Makefile.in for each directory"
automake --add-missing --copy
