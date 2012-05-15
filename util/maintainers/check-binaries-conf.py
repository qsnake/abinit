#!/usr/bin/env python
#
# Copyright (C) 2012 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

from ConfigParser import ConfigParser,NoOptionError
from time import gmtime,strftime

import commands
import os
import re
import sys

class MyConfigParser(ConfigParser):

  def optionxform(self,option):
    return str(option)

# ---------------------------------------------------------------------------- #

#
# Auxiliary data
#

dep_levels = {
  "algo":1,
  "atompaw":3,
  "bigdft":3,
  "etsf_io":2,
  "fft":1,
  "fox":1,
  "libxc":1,
  "linalg":1,
  "math":1,
  "mpi":0,
  "gpu":1,
  "netcdf":1,
  "timer":1,
  "wannier90":3}

# ---------------------------------------------------------------------------- #

#
# Main program
#

# Check if we are in the top of the ABINIT source tree
my_name = os.path.basename(sys.argv[0])
if ( not os.path.exists("configure.ac") or
     not os.path.exists("src/98_main/abinit.F90") ):
  print "%s: You must be in the top of an ABINIT source tree." % my_name
  print "%s: Aborting now." % my_name
  sys.exit(1)

# Init
cnf_bin = MyConfigParser()
cnf_bin.read("config/specs/binaries.conf")
bin_list = cnf_bin.sections()
bin_list.sort()
dep_order = dict()
lib_order = dict()
re_num = re.compile("[0-9][0-9]_")

# Check order of dependencies and libraries
for prg in bin_list:
  if ( cnf_bin.has_option(prg,"dependencies") ):
    bin_deps = cnf_bin.get(prg,"dependencies").split()
  else:
    bin_deps = list()
  dep_old = 100
  dep_new = 100
  for dep in bin_deps:
    if ( dep in dep_levels ):
      dep_new = dep_levels[dep]
    else:
      sys.stderr.write("%s: Error: unregistered dependency '%s'\n" % \
        (my_name,dep))
      sys.exit(10)
    if ( dep_new > dep_old ):
      if ( not prg in dep_order ):
        dep_order[prg] = list()
      dep_order[prg].append(dep)
    dep_old = dep_new

  if ( cnf_bin.has_option(prg,"libraries") ):
    bin_libs = cnf_bin.get(prg,"libraries").split()
  else:
    bin_libs = list()
  lib_old = 100
  lib_new = 100
  for lib in bin_libs:
    if ( re_num.match(lib) ):
      lib_new = int(re.sub("_.*","",lib))
      if ( lib_new > lib_old ):
        if ( not prg in lib_order ):
          lib_order[prg] = list()
        lib_order[prg].append(lib)
    lib_old = lib_new

# Report any disorder
nerr = len(dep_order) + len(lib_order)

if ( nerr > 0 ):
  sys.stderr.write("%s: reporting disordered libraries\n\n" % \
    (os.path.basename(sys.argv[0])))
  sys.stderr.write("X: D=Dependency / L=Library\n\n")
  sys.stderr.write("%s  %-24s  %-24s\n" % \
    ("X","Program","Dependency/Library"))
  sys.stderr.write("%s  %s  %s\n" % ("-","-" * 24,"-" * 24))

  dep_keys = dep_order.keys()
  dep_keys.sort()
  for prg in dep_keys:
    for dep in dep_order[prg]:
      sys.stderr.write("%s  %-24s  %-24s\n" % \
      ("D",prg,dep))
  lib_keys = lib_order.keys()
  lib_keys.sort()
  for prg in lib_keys:
    for lib in lib_order[prg]:
      sys.stderr.write("%s  %-24s  %-24s\n" % \
      ("L",prg,lib))

  sys.stderr.write("\n")
  sys.exit(1)
else:
  sys.exit(0)
