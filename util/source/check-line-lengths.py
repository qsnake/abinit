#!/usr/bin/env python

import re
import os
import sys

# Init
re_srcfile = re.compile("\.([Ff]|[Ff]90)$")
len_limit  = 132
retval     = 0

for root,dirs,files in os.walk("src"):
  # Sort dirs
  dirs.sort()

  # Check line lengths in Fortran source files
  for item in files:
    if ( re_srcfile.search(item) ):
      lineno = 1
      for line in file("%s/%s" % (root,item),"r").readlines():
        line = re.sub("!.*","",line)
        line = re.sub("\n","",line)
        if ( len(line) > len_limit ):
          sys.stderr.write("%s/%s: line %d has more than %d characters\n" % \
            (root,item,lineno,len_limit))
          retval = 1
        lineno += 1

sys.exit(retval)
