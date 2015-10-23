# this is <Makefile>
# ----------------------------------------------------------------------------
#   ($Id: Makefile 4961 2013-02-01 12:52:09Z lrehor $)
# 
# 25/10/2000 by Thomas Forbriger (IMGF Frankfurt)
# 
# libseife Makefile
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
# This library contains code which was extracted from the program seife.f
# by Erhard Wielandt. The original version was written in 1984 at ETH Zurich.
# The code can be obtained through 
# http://www.software-for-seismometry.de/
#
# The frontpage of http://www.software-for-seismometry.de/ tells (Jan. 2013):
#
#   Welcome to Erhard Wielandt's website for free seismological software
#
#   This website offers free software for practical use, especially for the
#   design, calibration, and testing of seismographs. It also contains
#   educational software and other material used in courses and workshops on
#   seismometry. Special problems of seismometry and theoretical seismology are
#   treated in the 'misc. texts'.
#   If you are interested in software, please start from the software overview. 
#
# The code from seife.f was extracted to a library by Wolfgang Friederich in
# 1996. The current version of libseife is a an extended code collection based
# on the abovementioned code.
# 
# REVISIONS and CHANGES
#    25/10/2000   V1.0   Thomas Forbriger
#    14/12/2007   V1.1   g77 compilation is the default now
#    17/12/2007   V1.2   moved to gfortran
#    18/09/2010   V1.3   start migration to new SVN scheme
#                        discard f2c option (fallback)
#                        check variables
#    15/11/2010   V1.4   do not use fdep.sh
#    17/11/2010   V1.5   libseifemon.a is out of use
#    17/01/2011   V1.6   distinguish library creation and installation
#    10/06/2011   V1.7   preparing Makefile for export of cseife to  
#                        the DENISE code by Lisa Rehor
# 
# ============================================================================
# 

.PHONY: all
all: install doc

.PHONY: doc
doc: 

LIBRARIES=libcseife.a

.PHONY: install
install: 
	$(MAKE) $(LIBRARIES)
	 
# ============================================================================

CLIBSRC=$(wildcard cseife*.c)
HEADERS=$(shell find . -name \*.h)

.PHONY: clean edit
clean:
	-find . -name \*.d | xargs --no-run-if-empty /bin/rm -v
	-find . -name \*.o | xargs --no-run-if-empty /bin/rm -v
	-find . -name \*.bak | xargs --no-run-if-empty /bin/rm -v
	-/bin/rm -fv flist
	-/bin/rm -fv *.a

edit: flist ; vim flist
flist: $(wildcard Makefile *.c *.h) COPYING
	echo $^ | tr ' ' '\n' | sort > $@

AR=ar
AS=as
RANLIB=ranlib
CC=mpicc

FLAGS += $(MYFLAGS)
CFLAGS += -O3 -Wall $(FLAGS)
CPPFLAGS=$(addprefix -I,.) \
  $(FLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -c $(CPPFLAGS)


#======================================================================
# reinstall target
# is necessary in case of header file problems (e.g. remcmmnt not installed)
.PHONY: reinstall
reinstall:
	$(MAKE) clean
	$(MAKE) install

#======================================================================
# pattern rules
# -------------

%.d: %.c
	$(SHELL) -ec '$(CC) -M $(CPPFLAGS) $< \
      | sed '\''s,\($(notdir $*)\)\.o[ :]*,$(dir $@)\1.o $@ : ,g'\'' \
        > $@; \
      [ -s $@ ] || rm -f $@'

#======================================================================
# library part
# ------------

include $(patsubst %.c,%.d,$(CLIBSRC))

libcseife.a: $(patsubst %.c,%.o,$(CLIBSRC))

libcseife.a: 
	$(AR) rcv $@ $^
	$(RANLIB) $@

# 
# ----- END OF Makefile ----- 
