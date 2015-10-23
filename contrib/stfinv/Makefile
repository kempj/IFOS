# this is <Makefile>
# ----------------------------------------------------------------------------
# $Id: Makefile 4968 2013-02-01 13:58:05Z lrehor $
# 
# Copyright (c) 2011 by Thomas Forbriger (BFO Schiltach) 
# 
# create library for source wavelet correction filter determination (libstfinv)
#
# ----
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version. 
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# ----
#
# See the README and the doxygen documentation for a list of libraries
# required in addition to libstfinv.
# 
# REVISIONS and CHANGES 
#    05/05/2011   V1.0   Thomas Forbriger
#
# ============================================================================
#

.PHONY: all
all: install doxydoc

.PHONY: doc
doc: doxydoc

LIBRARIES=libstfinv.a

.PHONY: install
install: $(addprefix $(LOCLIBDIR)/,$(LIBRARIES))
$(LOCLIBDIR)/%: install-include %
	mkdir -pv $(LOCLIBDIR)
	/bin/mv -fv $(word 2,$^) $(LOCLIBDIR)

# ============================================================================
#
CHECKVAR=$(if $($(1)),,$(error ERROR: missing variable $(1)))
CHECKVARS=$(foreach var,$(1),$(call CHECKVAR,$(var)))
#

$(call CHECKVARS,LOCINCLUDEDIR LOCLIBDIR LOCBINDIR)
$(call CHECKVARS,TF_BROWSER TF_WWWBASEDIR)

LIBHEADERS=$(filter-out none,$(wildcard *.h))

LIBCCSRC=$(filter-out none,$(wildcard *.cc))
LIBCSRC=$(filter-out none,$(wildcard *.c))

TESTHEADERS=$(wildcard tests/*.h)

TESTCCSRC=$(wildcard tests/*.cc)
TESTCSRC=$(wildcard tests/*.c)

TESTCCBIN=$(patsubst tests/%.cc,%,$(TESTCCSRC))
TESTCBIN=$(patsubst tests/%.c,%,$(TESTCSRC))

INCINSTALLPATH=$(LOCINCLUDEDIR)/stfinv
LIBINSTALLPATH=$(LOCLIBDIR)

# name of headers with comments stripped off
STRIPHEADER=$(addsuffix .strip,$(notdir $(LIBHEADERS)))
# name of installed (exported) header files
INSTHEADER=$(addprefix $(INCINSTALLPATH)/,$(notdir $(LIBHEADERS)))

# if defined, empty lines are kept in stripped comment headers
# to synchronize line numbers
EMPTYPRINT=1
#EMPTYPRINT=0

# define this to be cat in case you do not have remcmmnt available
ifndef TF_REMCMMNT
TF_REMCMMNT=cat
endif

#
# general part
# ------------
#

FLAGS=
FLAGS+=-fPIC $(MYFLAGS)
CXXFLAGS+=-Wall $(FLAGS)
LDFLAGS+=-L$(LOCLIBDIR)
CPPFLAGS+=-I$(LOCINCLUDEDIR) $(FLAGS)

EDITFILES=Makefile README $(wildcard *.cfg *.h *.cc *.c tests/*.cc tests/*.c) \
	    test_libstfinv/README test_libstfinv/Test_libstfinv.inp \
	    $(wildcard test_libstfinv/src/*.c) \
	    $(wildcard test_libstfinv/src/*.h) \
	    test_libstfinv/src/Makefile COPYING

flist: $(EDITFILES) 
	echo $^ | tr ' ' '\n' | sort > $@

.PHONY: edit
edit: flist; vim $<

.PHONY: clean
clean: ; 
	-find . -name \*.bak | xargs --no-run-if-empty /bin/rm -v
	-find . -name \*.o | xargs --no-run-if-empty /bin/rm -v
	-find . -name \*.d | xargs --no-run-if-empty /bin/rm -v
	-find . -name \*.hd | xargs --no-run-if-empty /bin/rm -v
	-find . -name \*.strip | xargs --no-run-if-empty /bin/rm -v
	-/bin/rm -vf flist *.a *.so *.xxx $(TESTCBIN) $(TESTCCBIN)

#======================================================================
# pattern rules
# -------------

%.d: %.cc
	$(SHELL) -ec '$(CXX) -M $(CPPFLAGS) $< \
      | sed '\''s,\($(notdir $*)\)\.o[ :]*,$(dir $@)\1.o $@ : ,g'\'' \
        > $@; \
      [ -s $@ ] || rm -f $@'

%.d: %.c
	$(SHELL) -ec '$(CC) -M $(CPPFLAGS) $< \
      | sed '\''s,\($(notdir $*)\)\.o[ :]*,$(dir $@)\1.o $@ : ,g'\'' \
        > $@; \
      [ -s $@ ] || rm -f $@'

#======================================================================
# library part
# ------------

include $(patsubst %.cc,%.d,$(LIBCCSRC) $(TESTCCSRC)) \
        $(patsubst %.c,%.d,$(LIBCSRC) $(TESTCSRC))

libstfinv.a: $(patsubst %.cc,%.o,$(LIBCCSRC)) \
             $(patsubst %.c,%.o,$(LIBCSRC))
	ar rcv $@ $^
	ranlib $@

#----------------------------------------------------------------------
# comment stripping

%.hd: %.h; echo $(notdir $<.strip): $< > $@

include $(patsubst %.h,%.hd,$(LIBHEADERS) $(TESTHEADERS))

# first awk rule assumes that the first non-comment line starts with '#'
# and that the first (copyright) comment end with pattern "^ */"
%.h.strip:
	awk 'BEGIN {hot=1;} /^ \*\// { if (hot) { hot=2; print; next;} }\
        /^#/ { hot=0; } \
        { if (hot==2) { print ""; } else if (hot) { print; } }' $< > $@
	$(TF_REMCMMNT) $< | awk 'BEGIN {hot=0;} \
        /^ *$$/ { if ((hot) && ($(EMPTYPRINT))) { print; } next; } \
        /^#/ { hot=1; } { if (hot) print; }' >> $@

#----------------------------------------------------------------------

$(INCINSTALLPATH)/%.h: %.h.strip
	mkdir -p $(INCINSTALLPATH)
	-rm -fv $@
	/bin/cp -vpd $< $@

.PHONY: install-include
install-include: $(INSTHEADER)

#======================================================================
# reinstall target
# is necessary in case of header file problems (e.g. remcmmnt not installed)
.PHONY: clean-include
clean-include: 
	/bin/rm -fv $(INSTHEADER)

.PHONY: reinstall
reinstall:
	$(MAKE) clean
	$(MAKE) clean-include
	$(MAKE) install

#======================================================================
# documentation part
# ------------------
#
# targets commonly used:
# ----------------------
#
# make doxyclean        removes all documentation
# make doxydoc          creates doxygen documentation in the DOXYWWWPATH
# make doxyview         creates doxygen documentation and launches netscape to
#                       browse in the documentation
# make doxyconf         edit the doxygen configuration file
#
# If you launch "make doxydoc" the documentation will be written to
# DOXYWWWPATH (see below). This is meant to export the documentation through
# your homepage. The doxyfull directory is just a symbolic link to this
# directory.
#

$(call CHECKVARS,TF_WWWBASEDIR TF_BROWSER)

DOXYWWWPATH=$(TF_WWWBASEDIR)/libstfinv

.PHONY: doxyclean doxyview doxydoc doxyconf

doxyclean: ;/bin/rm -rfv $(DOXYWWWPATH) doxydoc.xxx

onlinehelp.xxx: onlinehelp
	echo '/*! \page page_help Online help texts' > $@
	echo '\verbatim' >> $@
	./$< -help 2>&1 >> $@ 2>&1
	echo '\endverbatim' >> $@
	echo '*/' >> $@

DOXYSRC=$(LIBCCSRC) $(LIBCSRC) $(LIBHEADERS) README $(TESTCCSRC) $(TESTCSRC) \
  onlinehelp.xxx Makefile

PWD=$(shell env pwd)
# create doxygen intermediate configuration
doxydoc.xxx: doxydoc.cfg
	sed 's,<OUTPUTDIRECTORY>,$(DOXYWWWPATH),g;s,<STRIPFROMPATH>,$(PWD),g' \
	  $< > $@
# create commented version of doxygen configuration
doxycomm.xxx: doxydoc.cfg
	/bin/cp -vf $< $@; doxygen -u $@

$(DOXYWWWPATH)/html/index.html: doxydoc.xxx $(DOXYSRC)
	mkdir -vp $(DOXYWWWPATH)
	doxygen $<

doxydoc: $(DOXYWWWPATH)/html/index.html

doxyview: $(DOXYWWWPATH)/html/index.html
	$(TF_BROWSER) file:$< &

#======================================================================
# create package
# --------------
# is delegated to Makefile.packages
ifdef TF_MAKEPKG
.PHONY: package
package: $(TF_MAKEPKG)
	$(MAKE) -f $< \
	  PACKAGE=libstfinv \
	  PACKAGEEXPORT="trunk/src/libs/libstfinv:src" \
	  PACKAGETARGETS="src:all:" \
	  PACKAGELIBS="-"
.PHONY: srcpackage
srcpackage: $(TF_MAKEPKG)
	$(MAKE) -f $< \
	  PACKAGE=libstfinv \
	  PACKAGESRCONLY=yes \
	  PACKAGEEXPORT="trunk/src/libs/libstfinv:src" \
	  PACKAGETARGETS="src:all:" \
	  PACKAGELIBS="-"
.PHONY: fullpackage
fullpackage: $(TF_MAKEPKG)
	$(MAKE) -f $< \
	  PACKAGE=libstfinvwithlibs \
	  PACKAGEEXPORT="trunk/src/libs/libstfinv:src" \
	  PACKAGETARGETS="src:all:" \
	  PACKAGELIBS="libaff libfourier"
endif

#======================================================================
# test code
# ---------

$(TESTCCBIN): %: tests/%.o libstfinv.a
	$(CXX) $(FLAGS) -o $@ $^ -I$(LOCINCLUDEDIR) \
        -ltfxx -lfourierxx -lfftw3 -laff \
        -L$(LOCLIBDIR)

# ----- END OF Makefile ----- 
