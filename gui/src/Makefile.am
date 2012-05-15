#
# Makefile for the GUI of Abinit 6
#

DIST_SUBDIRS = abinitgui com org

if DO_BUILD_GUI

pkgdata_DATA = $(gui_objdir)/$(gui_progname).jar

all_targets all: $(gui_objdir)/$(gui_progname).jar

build-stamp:
	touch build-stamp

# Compile JSCH
jsch-stamp: build-stamp
	@echo "Compilation of the JSCH library..."
	$(JAVAC) $(JAVACFLAGS) -d $(gui_objdir) $(srcdir)/com/jcraft/jzlib/*.java
	$(JAVAC) $(JAVACFLAGS) -d $(gui_objdir) $(srcdir)/com/jcraft/jsch/*.java
	$(JAVAC) $(JAVACFLAGS) -d $(gui_objdir) $(srcdir)/com/jcraft/jsch/jce/*.java
	$(JAVAC) $(JAVACFLAGS) -d $(gui_objdir) $(srcdir)/com/jcraft/jsch/jcraft/*.java
	$(JAVAC) $(JAVACFLAGS) -d $(gui_objdir) $(srcdir)/com/jcraft/jsch/jgss/*.java
	touch jsch-stamp
	@echo "JSCH library compiled."

# Compile swing-layout
swing-stamp: jsch-stamp
	@echo "Compilation of the swing-layout library..."
	$(JAVAC) $(JAVACFLAGS) -d $(gui_objdir) $(srcdir)/org/jdesktop/layout/*.java
	touch swing-stamp
	@echo "Swing-layout library compiled."

# Compile JDOM
jdom-stamp: swing-stamp
	@echo "Compilation of the JDOM library..."
	$(JAVAC) $(JAVACFLAGS) -d $(gui_objdir) -sourcepath $(srcdir) \
	  $(srcdir)/org/jdom/*.java
	$(JAVAC) $(JAVACFLAGS) -d $(gui_objdir) $(srcdir)/org/jdom/input/*.java
	touch jdom-stamp
	@echo "JDOM library compiled."

# Compile AbinitGUI
abinit-gui-stamp: jdom-stamp
	@echo "Compilation of the AbinitGUI program..."
	$(JAVAC) $(JAVACFLAGS) -d $(gui_objdir) -cp $(gui_objdir) \
	  $(srcdir)/abinitgui/*.java
	$(JAVAC) $(JAVACFLAGS) -d $(gui_objdir) -cp $(gui_objdir) \
	  $(srcdir)/Main.java
	touch abinit-gui-stamp
	@echo "AbinitGUI program compiled."

# Create AbinitGUI.jar
$(gui_objdir)/$(gui_progname).jar: abinit-gui-stamp
	@echo "JAR file creation..."
	cd $(gui_objdir) && \
	  $(JAR) $(JARFLAGS) $(gui_progname).jar Main *
	chmod +x $(gui_objdir)/$(gui_progname).jar
	@echo "JAR file created."

clean-local:
	rm -rf $(gui_objdir)/* *-stamp

endif
