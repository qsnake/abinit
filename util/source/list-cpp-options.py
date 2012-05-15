#!/usr/bin/env python

import os
import re

fortran = re.compile("\.([Ff]|[Ff]90)$")
cppline = re.compile("^#")
cppopts = dict()
cppkeys = ("define .*","include.*","ifdef","ifndef","elif","^if ","else","endif","defined","undef","!","&&","\|\|","\(","\)")

for root,dirs,files in os.walk('src'):
  for src in files:
    if ( fortran.search(src) ):
      code = file(os.path.join(root,src),"r").readlines()

      for line in code:
        if ( cppline.match(line) ):
          line = re.sub("^#","",line).strip()
          for kw in cppkeys:
            line = re.sub(kw,"",line)

          line = line.split()
          for item in line:
            if ( item in cppopts ):
              cppopts[item] += 1
            else:
              cppopts[item] = 1

print "Option                             Occurences"
print "--------------------------------   ----------"
names = cppopts.keys()
names.sort()
for opt in names:
  print "%-32s   %10d" % (opt,cppopts[opt])
print ""
