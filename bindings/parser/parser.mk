#
# Makefile for the parser bindings of ABINIT
#

include ../../config.mk

VPATH = $(abinit_srcdir)/bindings/parser \
	$(abinit_srcdir)/src/42_geometry \
	$(abinit_srcdir)/src/56_recipspace

all: libabinis.a

python: ab6_invars.so

# The building of the library itself.
libabinit_tmpdir = tmp-libabinis-objects
libabinis.a: libbindings.a ../../src/libs/libab6_parser.a
	test -e "$(libabinit_tmpdir)" || \
	$(INSTALL) -d -m 755 $(libabinit_tmpdir)
	cd $(libabinit_tmpdir) && $(AR) x ../../../src/libs/libab6_parser.a
	cd $(libabinit_tmpdir) && $(AR) x ../libbindings.a
	$(AR) $(ARFLAGS) libabinis.a $(libabinit_tmpdir)/*
	$(RANLIB) libabinis.a
	rm -rf $(libabinit_tmpdir)

ab6_base.o: $(srcdir)/ab6_base.c $(srcdir)/ab6_base.h
	$(CC) $(CPPFLAGS) -I. -I$(srcdir) -I$(abinit_builddir) $(CFLAGS) -c $(abinit_srcdir)/bindings/parser/ab6_base.c

ab6_invars_c.o: $(srcdir)/ab6_invars_c.c $(srcdir)/ab6_invars_c.h $(srcdir)/ab6_invars.h
	$(CC) $(CPPFLAGS) -I. -I$(srcdir) -I$(abinit_builddir) $(CFLAGS) -c $(abinit_srcdir)/bindings/parser/ab6_invars_c.c
ab6_symmetry.o: ab6_symmetry.c ab6_symmetry.h ab6_symmetry.fortran.h
	$(CC) $(CPPFLAGS) -I. -I$(srcdir) -I$(abinit_builddir) -I$(abinit_srcdir)/src/42_geometry $(CFLAGS) -c $(abinit_srcdir)/bindings/parser/ab6_symmetry.c
ab6_kpoints.o: ab6_kpoints.c ab6_kpoints.h ab6_kpoints.fortran.h
	$(CC) $(CPPFLAGS) -I. -I$(srcdir) -I$(abinit_builddir) -I$(abinit_srcdir)/src/56_recipspace $(CFLAGS) -c $(abinit_srcdir)/bindings/parser/ab6_kpoints.c

ab6_dummy.o: $(srcdir)/ab6_dummy.f90
	$(FC) $(FCFLAGS) -c $(abinit_srcdir)/bindings/parser/ab6_dummy.f90

libbindings.a: ab6_base.o ab6_invars_c.o ab6_symmetry.o ab6_kpoints.o ab6_dummy.o
	$(AR) $(ARFLAGS) libbindings.a ab6_base.o ab6_invars_c.o ab6_symmetry.o ab6_kpoints.o ab6_dummy.o
	$(RANLIB) libbindings.a

ab6_invars_py.o: $(srcdir)/ab6_invars_py.h $(srcdir)/ab6_invars.h $(srcdir)/ab6_invars_py.c
	$(CC) $(CPPFLAGS) $(PYTHON_CPPFLAGS) -I. -I$(srcdir) -I$(abinit_builddir) $(CFLAGS) -c $(abinit_srcdir)/bindings/parser/ab6_invars_py.c


ab6_invars.so: ab6_invars_py.o
	$(FC) -shared -o ab6_invars.so ab6_invars_py.o libabinis.a

clean:
	rm -f *.o *.mod libbindings.a libabinis.a
	rm -f ab6_invars.so
	rm -f check_invars_c check_invars_f90
	rm -f check_symmetry_f90 check_symmetry_c
	rm -f *.out



# Some additional rules for make check
check_invars_c: $(srcdir)/check_invars_c.c fallbacks.o libabinis.a
	$(FC) -o check_invars_c -I. -I$(srcdir) $(abinit_srcdir)/bindings/parser/check_invars_c.c libabinis.a fallbacks.o $(lib_linalg_libs) $(lib_etsf_io_libs) $(lib_netcdf_libs) $(lib_libxc_libs)

check_invars_f90: $(srcdir)/check_invars_f90.f90 fallbacks.o libabinis.a
	$(FC) -I. -I$(srcdir) $(FCFLAGS) $(FCFLAGS_MODDIR) -o check_invars_f90 $(abinit_srcdir)/bindings/parser/check_invars_f90.f90 libabinis.a fallbacks.o $(lib_linalg_libs) $(lib_etsf_io_libs) $(lib_netcdf_libs) $(lib_libxc_libs)

fallbacks.o: $(srcdir)/fallbacks.f90
	$(FC) -I. -I$(srcdir) $(FCFLAGS) $(FCFLAGS_MODDIR) -c $(abinit_srcdir)/bindings/parser/fallbacks.f90

check_symmetry_f90: check_symmetry_f90.f90 libabinis.a fallbacks.o
	$(FC) -I. -I$(srcdir) $(FCFLAGS) $(FCFLAGS_MODDIR) -o check_symmetry_f90 $(abinit_srcdir)/bindings/parser/check_symmetry_f90.f90 -L. -labinis fallbacks.o $(lib_linalg_libs) $(lib_etsf_io_libs) $(lib_netcdf_libs) $(lib_libxc_libs)

check_symmetry_c: $(srcdir)/check_symmetry_c.c libabinis.a fallbacks.o
	$(FC) -I. -I$(srcdir) $(CFLAGS) -o check_symmetry_c $(abinit_srcdir)/bindings/parser/check_symmetry_c.c -L. -labinis fallbacks.o $(lib_linalg_libs) $(lib_etsf_io_libs) $(lib_netcdf_libs) $(lib_libxc_libs)

check: check_symmetry_f90 check_symmetry_c check_invars_f90 check_invars_c
	./check_symmetry_f90 > check_symmetry_f90.out
	diff $(abinit_srcdir)/bindings/parser/check_symmetry.ref check_symmetry_f90.out
	./check_symmetry_c  > check_symmetry_c.out
	diff $(abinit_srcdir)/bindings/parser/check_symmetry.ref check_symmetry_c.out
	./check_invars_f90 $(abinit_srcdir)/tests/v1/Input/t05.in > check_invars_f90.out
	diff $(abinit_srcdir)/bindings/parser/check_invars.ref check_invars_f90.out
	./check_invars_c $(abinit_srcdir)/tests/v1/Input/t05.in > check_invars_c.out
	diff $(abinit_srcdir)/bindings/parser/check_invars.ref check_invars_c.out
