#!/usr/bin/env python
#
# Copyright (C) 2010-2012 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

# FIXME: detect duplicate definitions

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
# Functions
#

env_ignore = ["DEFS"]
opt_ignore = [
  "enable_config_file",
  "fcflags_opt_.*",
  "group_.*",
  "prefix",
  "with_config_file"]
val_ignore = [".*-fallback"]

def is_ignored(keyword):
  for opt in env_ignore + opt_ignore:
    if ( "*" in opt ):
      if ( re.match(opt,keyword) ):
        return True
    elif ( opt == keyword ):
        return True
  return False

def is_ignored_val(keyword):
  for val in val_ignore:
    if ( "*" in val ):
      if ( re.match(val,keyword) ):
        return True
    elif ( val == keyword ):
        return True
  return False

# ---------------------------------------------------------------------------- #

#
# Main program
#

# Check if we are in the top of the ABINIT source tree
if ( not os.path.exists("configure.ac") or
     not os.path.exists("src/98_main/abinit.F90") ):
  print "%s: You must be in the top of an ABINIT source tree." % my_name
  print "%s: Aborting now." % my_name
  sys.exit(1)

# Init
re_env = re.compile("^#[A-Z][0-9A-Z_]*=")
re_opt = re.compile("^#[a-z][0-9a-z_]*=")

# Extract environment variables from config file
cnf_env = MyConfigParser()
cnf_env.read("config/specs/environment.conf")
env_config = list()
for env in cnf_env.sections():
  if ( cnf_env.get(env,"reset") == "no" ):
    if ( not is_ignored(env) ):
      env_config.append(env)
env_config.sort()

# Extract options from config file
cnf_opt = MyConfigParser()
cnf_opt.read("config/specs/options.conf")
opt_config = list()
opt_removed = list()
for opt in cnf_opt.sections():
  tmp_sta = cnf_opt.get(opt,"status")
  if ( tmp_sta == "removed" or tmp_sta == "dropped" ):
    opt_removed.append(opt)
  elif  ( "renamed" in tmp_sta ):
    opt_removed.append(tmp_sta.split()[1])
    if ( not is_ignored(opt) ):
      opt_config.append(opt)
  elif ( tmp_sta == "hidden" ):
    opt_ignore.append(opt)
  else:
    if ( not is_ignored(opt) ):
      opt_config.append(opt)
opt_config.sort()
opt_removed.sort()

# Extract information from template
env_template = list()
opt_template = list()
for line in file("doc/config/build-config.ac","r").readlines():
  if ( re_env.match(line) ):
    tmp_env = re.sub("=.*","",line[1:-1])
    if ( not is_ignored(tmp_env) ):
      if ( not tmp_env in env_template ):
        env_template.append(tmp_env)
  if ( re_opt.match(line) ):
    tmp_opt = re.sub("=.*","",line[1:-1])
    if ( not is_ignored(tmp_opt) ):
      if ( not tmp_opt in opt_template ):
        opt_template.append(tmp_opt)
env_template.sort()
opt_template.sort()

# Check whether non-trivial option values are found in template
tpl_data = file("doc/config/build-config.ac","r").read()
opt_values = dict()
for opt in cnf_opt.sections():
  try:
    tmp_values = cnf_opt.get(opt,"values").split()
    if ( "no" in tmp_values ):
      tmp_values.remove("no")
    if ( "yes" in tmp_values ):
      tmp_values.remove("yes")
  except NoOptionError:
    tmp_values = list()

  for val in tmp_values:
    if ( (val[0] != "@") and (not is_ignored_val(val)) ):
      if ( not re.search("\\* %s" % (val),tpl_data,re.MULTILINE) ):
        if ( not opt in opt_values ):
          opt_values[opt] = list()
        opt_values[opt].append(val)
opt_valkeys = opt_values.keys()
opt_valkeys.sort()

# Compare environment and options
denv_config = [env for env in env_config if not env in env_template]
denv_template = [env for env in env_template if not env in env_config]
dopt_config = [opt for opt in opt_config if not opt in opt_template]
dopt_values = [opt for opt in opt_valkeys if opt not in dopt_config]
dopt_template = [opt for opt in opt_template if not opt in opt_config + opt_removed]
dopt_removed = [opt for opt in opt_removed if opt in opt_template]

# Report any mismatch
nerr = len(denv_config) + len(denv_template) + len(dopt_config) + \
  + len(dopt_values) + len(dopt_template) + len(dopt_removed)

if ( nerr > 0 ):
  sys.stderr.write("%s: reporting environment and option mismatches\n\n" % \
    (os.path.basename(sys.argv[0])))
  sys.stderr.write("X: D=Documentation / I=Instance / R=Removed / U=Undefined\n\n")
  sys.stderr.write("%s  %-48s  %-24s\n" % \
    ("X","Option","Outdated file"))
  sys.stderr.write("%s  %s  %s\n" % ("-","-" * 48,"-" * 24))

  for env in denv_config:
    sys.stderr.write("%s  %-48s  %-24s\n" % \
      ("I",env,"build-config.ac"))
  for env in denv_template:
    sys.stderr.write("%s  %-48s  %-24s\n" % \
      ("U",env,"environment.conf"))
  for opt in dopt_config:
    sys.stderr.write("%s  %-48s  %-24s\n" % \
      ("I",opt,"build-config.ac"))
  for opt in dopt_values:
    for val in opt_values[opt]:
      sys.stderr.write("%s  %-48s  %-24s\n" % \
      ("D","%s='%s'" % (opt,val),"build-config.ac"))
  for opt in dopt_template:
    sys.stderr.write("%s  %-48s  %-24s\n" % \
      ("U",opt,"options.conf"))
  for opt in dopt_removed:
    sys.stderr.write("%s  %-48s  %-24s\n" % \
      ("R",opt,"build-config.ac"))

  sys.stderr.write("\n")
  sys.exit(1)
else:
  sys.exit(0)
