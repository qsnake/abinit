#
# Top makefile for the Abinit GUI
#

ACLOCAL_AMFLAGS = -I m4

pkgdata_DATA = $(srcdir)/precompiled/config.xml

SUBDIRS = src tests

if !DO_BUILD_GUI
  pkgdata_DATA += $(srcdir)/precompiled/$(gui_progname).jar
endif

bin_SCRIPTS = abinit-gui

CLEANFILES = abinit-gui

abinit-gui:
	echo "#!$(BOURNE_SHELL)" >abinit-gui
	echo "" >>abinit-gui
	echo "JAVA='$(JAVA)'" >>abinit-gui
	echo "" >>abinit-gui
	echo "if test \"\$${JAVA}\" != \"\"; then" >>abinit-gui
	echo "  $(JAVA) -jar $(pkgdatadir)/$(gui_progname).jar" >>abinit-gui
	echo "else" >>abinit-gui
	echo "  echo 'Java was not available when Abinit was built'" >>abinit-gui
	echo "  exit 1" >>abinit-gui
	echo "fi" >>abinit-gui

