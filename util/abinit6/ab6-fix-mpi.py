#!/usr/bin/env python

import os
import re
import sys

# Initial setup
my_name    = "ab6-fix-mpi"

# Check if we are in the top of the ABINIT source tree
if ( not os.path.exists("configure.ac") or
     not os.path.exists("src/main/abinit.F90") ):
 print "%s: You must be in the top of an ABINIT source tree." % my_name
 print "%s: Aborting now." % my_name
 sys.exit(1)

f90_name = re.compile("\.F90")
cpp_line = re.compile("^#")

# Fix MPI statements in source files
for root,dirs,files in os.walk("src"):
 for src in files:
  print "Checking %s" % (os.path.join(root,src))
  if ( f90_name.search(src) ):
   print "---> %s is F90" % (os.path.join(root,src))
   f90_text = file(os.path.join(root,src),"r").readlines()

   for i in range(len(f90_text)):
    if ( cpp_line.match(f90_text[i]) ):
     f90_text[i] = re.sub(" MPI"," HAVE_MPI",f90_text[i])

   file(os.path.join(root,src),"w").write("".join(f90_text))

