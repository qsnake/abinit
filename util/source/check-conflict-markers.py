#!/usr/bin/env python

import re
import os
import sys

re_markers = re.compile("^(<<<<<<< TREE|=======|>>>>>>> MERGE-SOURCE)$")
re_fbktop  = re.compile("fallbacks$")
re_fbkdir  = re.compile("(exports|sources|stamps)")
re_tmpdir  = re.compile("^tmp")
re_tmpfile = re.compile("\.(orig|rej)$")

retval = 0

for root,dirs,files in os.walk("."):

  # Ignore Makefiles
  if ( "Makefile.am" in files ):
   files.remove("Makefile.am")
  if ( "Makefile.in" in files ):
   files.remove("Makefile.in")
  if ( "Makefile" in files ):
   files.remove("Makefile")

  # Ignore Autotools subdirs
  if ( "autom4te.cache" in dirs ):
    dirs.remove("autom4te.cache")

  # Ignore Bazaar subdirs
  if ( ".bzr" in dirs ):
    dirs.remove(".bzr")

  # Ignore temporary dirs
  garb_dirs = [item for item in dirs if re_tmpdir.match(item)]
  for d in garb_dirs:
    dirs.remove(d)

  # Ignore installed fallbacks
  if ( re_fbktop.search(root) ):
    garb_dirs = [item for item in dirs if re_fbkdir.match(item)]
    for d in garb_dirs:
      dirs.remove(d)

  # Display conflict markers found
  for item in files:
    if ( not re_tmpfile.search(item) ):
      chk_data = file("%s/%s" % (root,item),"r").readlines()
      chk_stat = False
      for line in chk_data:
        if ( re_markers.match(line) ):
          chk_stat = True
          retval = 1
      if ( chk_stat ):
        sys.stderr.write("Conflict markers in %s/%s\n" % (root,item))

sys.exit(retval)
