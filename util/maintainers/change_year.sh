# Copyright (C) 1998-2012 ABINIT group (XG)
# 
# The purpose of this script is to change the copyright year
# in nearly all files in the ABINIT package. 
# First you should reapply the present script without changing it, because it has been seen that some people bring
# routines with erroneous date during a few months after the beginning of a new one ... !
#
# Then one should update the present script to put the current year (at present, valid for going from 2010 to 2011 !) !
#
# Then should be called from the top directory, with
# util/maintainers/change_year.sh */*.in */*/*.F90 */*/*.cf */*/*.html */*/*.tex */*/*.pl README */README */*/README */*/*.cnf */*/*.bat */*/*.com */*/*.src */*/*_ */*/*.mk */*/*.txt */*/*/*.pl */*/*.m4 */*/make* tests/*/*.sh */*/*.csh */*/*.py */*/*.in *.ac */*/*.am */*/*.help */*/*.c */*/*.m */*/*/*.out */*/*.h */*/*.cu */*/*.sav */*/*/*.log */*/*.conf */*/*.dep */*/*.dir */*/Makefile */*.env 
# If one works in a directory in which a compilation has already been done, then
# util/maintainers/change_year.sh */*/*.dep */*/*.dir */*/Makefile
# can be performed as well, in order to remove any problem. But these files are not archived, so this is a bit useless...
# util/maintainers/change_year.sh config/scripts/add-header-typed config/scripts/add-targets-binaries config/scripts/add-targets-libraries config/scripts/clean-build-examples config/scripts/clean-source-tree doc/manpages/abinit.1 extras/bzr_helpers/abinit-forge-branch util/developers/mkroutine.sh util/maintainers/change2.sh util/maintainers/change.sh util/misc/change_perl.sh util/misc/fixed_to_free tests/cpu/Refs/changeref 
#then
# util/maintainers/change_year.sh INSTALL config/scripts/abilint config/scripts/reshape-options-conf config/scripts/update-options-conf config/scripts/update-version-number config/wrappers/wrap-fc doc/versioning/bzr-quickref.tex util/developers/mkmodule.sh fallbacks/config/m4/conf-tricks.m4 fallbacks/config/m4/util-fixes.m4 fallbacks/config/scripts/make-macros-dirs fallbacks/config/scripts/make-macros-fallbacks fallbacks/config/scripts/make-makefiles-fallbacks doc/versioning/bzr-quickref.tex
# 
# If one works in a directory in which a compilation has already been done, then
# util/maintainers/change_year.sh */*/*.dep */*/*.dir */*/Makefile
#
# In the previous list, files without an extension are not treated (except the README files), 
# and */*/*.sh are not treated (except tests/*/*.sh), because of conflict with the present script file extension !!
# Also config/scripts/abilint cannot be treated automatically
# So, one should complement the present script with a search 
# grep 'past_year ABINIT' * */* */*/* */*/*/* */*/*/*/*
# and treat by hand the remaining files ...
#XG 100118 Remarked still many other problems with copyrights, by using the following command (replace 2011 by the present year !):
# grep -i opyright * */* */*/* */*/*/* */*/*/*/* | grep -v 2011 | grep -v '!! COPYRIGHT' | grep -v 'Oldenburg' | grep -v 'Stefan Goedecker' | grep -v 'doc/rel' | grep -v 'Remove' | grep -v 'tests/' | grep -v 'EXC group' | grep -v 'PWSCF group' | grep -v 'Makefile' | grep -v 'abinit.d' | grep -v 'fallbacks' | grep -v 'doc/features/features' | grep -v 'doc/install_notes/install' | grep -v 'COPYING'

for file in "$@"
do
 echo "working on $file"
 rm -f tmp.yr*  
 sed -e 's&Copyright (c)&Copyright (C)&' $file > tmp.yrup
 sed -e 's&(C) 1987-2011 ABINIT&(C) 1987-2012 ABINIT&' tmp.yrup > tmp.yr87
 sed -e 's&(C) 1991-2011 ABINIT&(C) 1991-2012 ABINIT&' tmp.yr87 > tmp.yr91
 sed -e 's&(C) 1992-2011 ABINIT&(C) 1992-2012 ABINIT&' tmp.yr91 > tmp.yr92
 sed -e 's&(C) 1993-2011 ABINIT&(C) 1993-2012 ABINIT&' tmp.yr92 > tmp.yr93
 sed -e 's&(C) 1996-2011 ABINIT&(C) 1996-2012 ABINIT&' tmp.yr93 > tmp.yr96
 sed -e 's&(C) 1997-2011 ABINIT&(C) 1997-2012 ABINIT&' tmp.yr96 > tmp.yr97
 sed -e 's&(C) 1998-2011 ABINIT&(C) 1998-2012 ABINIT&' tmp.yr97 > tmp.yr98
 sed -e 's&(C) 1999-2011 ABINIT&(C) 1999-2012 ABINIT&' tmp.yr98 > tmp.yr99
 sed -e 's&(C) 2000-2011 ABINIT&(C) 2000-2012 ABINIT&' tmp.yr99 > tmp.yr00
 sed -e 's&(C) 2001-2011 ABINIT&(C) 2001-2012 ABINIT&' tmp.yr00 > tmp.yr01
 sed -e 's&(C) 2002-2011 ABINIT&(C) 2002-2012 ABINIT&' tmp.yr01 > tmp.yr02
 sed -e 's&(C) 2003-2011 ABINIT&(C) 2003-2012 ABINIT&' tmp.yr02 > tmp.yr03
 sed -e 's&(C) 2004-2011 ABINIT&(C) 2004-2012 ABINIT&' tmp.yr03 > tmp.yr04
 sed -e 's&(C) 2005-2011 ABINIT&(C) 2005-2012 ABINIT&' tmp.yr04 > tmp.yr05
 sed -e 's&(C) 2006-2011 ABINIT&(C) 2006-2012 ABINIT&' tmp.yr05 > tmp.yr06
 sed -e 's&(C) 2007-2011 ABINIT&(C) 2007-2012 ABINIT&' tmp.yr06 > tmp.yr07
 sed -e 's&(C) 2008-2011 ABINIT&(C) 2008-2012 ABINIT&' tmp.yr07 > tmp.yr08
 sed -e 's&(C) 2009-2011 ABINIT&(C) 2009-2012 ABINIT&' tmp.yr08 > tmp.yr09
 sed -e 's&(C) 2010-2011 ABINIT&(C) 2010-2012 ABINIT&' tmp.yr09 > tmp.yr10
#The next lines are both needed, as some developers decide to use one, and some the other ...
 sed -e 's&(C) 2011-2011 ABINIT&(C) 2011-2012 ABINIT&' tmp.yr10 > tmp.yr11
 sed -e 's&(C) 2011 ABINIT&(C) 2011-2012 ABINIT&' tmp.yr11 > tmp.yr
 echo "changes done "
 # put the modified file at the correct place
 mv tmp.yr $file
 echo "file $file written "
done
rm -f tmp.yr*  
chmod 755 */*/*.sh */*/*.py */*/*.pl */*/*.com config/*/make* extras/*/make* 
chmod 755 config/scripts/* extras/bzr_helpers/* util/*/* tests/cpu/Refs/changeref
