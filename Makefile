# $Id: Makefile 16563 2008-08-15 16:08:28Z bangerth $


# For the small projects Makefile, you basically need to fill in only
# four fields.
#
# The first is the name of the application. It is assumed that the
# application name is the same as the base file name of the single C++
# file from which the application is generated.


# The second field determines whether you want to run your program in
# debug or optimized mode. The latter is significantly faster, but no
# run-time checking of parameters and internal states is performed, so
# you should set this value to `on' while you develop your program,
# and to `off' when running production computations.
debug-mode = off 
#ANGIOGENESIS = 1
MECHANICS = 1
# As third field, we need to give the path to the top-level deal.II
# directory. You need to adjust this to your needs. Since this path is
# probably the most often needed one in the Makefile internals, it is
# designated by a single-character variable, since that can be
# reference using $D only, i.e. without the parentheses that are
# required for most other parameters, as e.g. in $(target).
D = /home/yangl4/deal.II_7.2.0/deal.II/


# The last field specifies the names of data and other files that
# shall be deleted when calling `make clean'. Object and backup files,
# executables and the like are removed anyway. Here, we give a list of
# files in the various output formats that deal.II supports.
clean-up-files = *gmv *gnuplot *gpl *eps *pov




#
#
# Usually, you will not need to change anything beyond this point.
#
#
# The next statement tell the `make' program where to find the
# deal.II top level directory and to include the file with the global
# settings
include $D/common/Make.global_options

#.PHONY: print_var

#print_vars:
#	echo $(INCLUDE)
# Since the whole project consists of only one file, we need not
# consider difficult dependencies. We only have to declare the
# libraries which we want to link to the object file, and there need
# to be two sets of libraries: one for the debug mode version of the
# application and one for the optimized mode. Here we have selected
# the versions for 2d. Note that the order in which the libraries are
# given here is important and that your applications won't link
# properly if they are given in another order.
#
# You may need to augment the lists of libraries when compiling your
# program for other dimensions, or when using third party libraries
libs.g   = $(lib-deal2.g) \
	   $(lib-lac.g)      \
           $(lib-base.g)
libs.o   = $(lib-deal2.o) \
	   $(lib-lac.o)      \
           $(lib-base.o)


# We now use the variable defined above which switch between debug and
# optimized mode to select the set of libraries to link with. Included
# in the list of libraries is the name of th

# Now comes the first production rule: how to link the single object
# file produced from the single C++ file into the executable. Since
# this is the first rule in the Makefile, it is the one `make' selects
# If you call it without arguments.

ifeq ($(debug-mode),on)
  woundHealing2D_V6 : woundHealing2D_V6.g.$(OBJEXT) Flist.g.$(OBJEXT) Elist.g.$(OBJEXT) ECM.g.$(OBJEXT) chemokine.g.$(OBJEXT) BMP.g.$(OBJEXT) step8_wound_mechanics.g.$(OBJEXT) $(libs.g)
	@echo ===========================Linking $@
	@$(CXX) -I $D/include/ -o $@$(EXEEXT) $^ $(LIBS) $(LDFLAGS)
else
  woundHealing2D_V6 : util.$(OBJEXT) Flist.$(OBJEXT) Elist.$(OBJEXT) ECM.$(OBJEXT) chemokine.$(OBJEXT) ADI.$(OBJEXT) BMP.$(OBJEXT) step8_wound_mechanics.$(OBJEXT) woundHealing2D_V6.$(OBJEXT) $(libs.o)
	@echo ===========================Linking $@
	g++ -std=c++11 -I $D/include/ -o $@$(EXEEXT) $^ $(LIBS) $(LDFLAGS)
endif

ifeq ($(debug-mode),on)
  woundHealing2D_V6.g.$(OBJEXT) : woundHealing2D_V6.cpp
        @echo ==============debug========= $(<F)
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
  parameters.g.$(OBJEXT) : parameters.cpp
	@echo ==============debug========= $(<F)
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
  Flist.g.$(OBJEXT) : Flist.cpp BMP.cpp chemokine.o
	@echo ==============debug========= $(<F)
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
  Elist.g.$(OBJEXT) : Elist.cpp
	@echo ==============debug========= $(<F)
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
  ECM.g.$(OBJEXT) : ECM.cpp BMP.cpp
	@echo ==============debug========= $(<F)
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
  chemokine.g.$(OBJEXT) : chemokine.cpp
	@echo ==============debug========= $(<F)
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
  BMP.g.$(OBJEXT) : BMP.cpp
	@echo ==============debug========= $(<F)
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
  step8_wound_mechanics.g.$(OBJEXT) : step8_wound_mechanics.cc  
	@echo ==============debug========= $(<F)
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@

else
  parameters.$(OBJEXT) : parameters.cpp
	@echo ==============optimized===== $(<F)
	@g++ $(CXXFLAGS.o) -c $< -o $@
  util.$(OBJEXT) : util.cpp
	@echo ==============optimized===== $(<F)
	@g++ $(CXXFLAGS.o) -c $< -o $@
  Flist.$(OBJEXT) : Flist.cpp ECM.o BMP.cpp chemokine.o
	@echo ==============optimized===== $(<F)
	@g++ $(CXXFLAGS.o) -c $< -o $@
  Elist.$(OBJEXT) : angiogenesis.cpp
	@echo ==============optimized===== $(<F)
	@g++ $(CXXFLAGS.o) -c $< -o $@ 
  ECM.$(OBJEXT) : ECM.cpp
	@echo ==============optimized===== $(<F)
	@g++ $(CXXFLAGS.o) -c $< -o $@ 
  chemokine.$(OBJEXT) : chemokine.cpp
	@echo ==============optimized===== $(<F)
	@g++ $(CXXFLAGS.o) -c $< -o $@ 
  ADI.$(OBJEXT) : ADI_Neumann.cpp
	@echo ==============optimized===== $(<F)
	@g++ $(CXXFLAGS.o) -c $< -o $@
  BMP.$(OBJEXT) : BMP.cpp
	@echo ==============optimized===== $(<F)
	@g++ $(CXXFLAGS.o) -c $< -o $@ 
  step8_wound_mechanics.$(OBJEXT) : step8_wound_mechanics.cc
	@echo ==============optimized===== $(<F)
	@g++ $(CXXFLAGS.o) -I $D/include/ -c $< -o $@
  woundHealing2D_V6.$(OBJEXT) : woundHealing2D_V6.cpp
	@echo ==============optimized===== $(<F)
	@g++ $(CXXFLAGS.o) -c $< -o $@
endif



# To make running the application somewhat independent of the actual
# program name, we usually declare a rule `run' which simply runs the
# program. You can then run it by typing `make run'. This is also
# useful if you want to call the executable with arguments which do
# not change frequently. You may then want to add them to the
# following rule:


# As a last rule to the `make' program, we define what to do when
# cleaning up a directory. This usually involves deleting object files
# and other automatically created files such as the executable itself,
# backup files, and data files. Since the latter are not usually quite
# diverse, you needed to declare them at the top of this file.
clean:
	-rm -f *.$(OBJEXT) *~ Makefile.dep $(target)$(EXEEXT) $(clean-up-files)


# Since we have not yet stated how to make an object file from a C++
# file, we should do so now. Since the many flags passed to the
# compiler are usually not of much interest, we suppress the actual
# command line using the `at' sign in the first column of the rules
# and write the string indicating what we do instead.


# The following statement tells make that the rules `run' and `clean'
# are not expected to produce files of the same name as Makefile rules
# usually do.
.PHONY: run clean


# Finally there is a rule which you normally need not care much about:
# since the executable depends on some include files from the library,
# besides the C++ application file of course, it is necessary to
# re-generate the executable when one of the files it depends on has
# changed. The following rule to created a dependency file
# `Makefile.dep', which `make' uses to determine when to regenerate
# the executable. This file is automagically remade whenever needed,
# i.e. whenever one of the cc-/h-files changed. Make detects whether
# to remake this file upon inclusion at the bottom of this file.
#
# If the creation of Makefile.dep fails, blow it away and fail
Makefile.dep: woundHealing2D_V6.cpp Makefile \
              $(shell echo $D/include/*/*/*.h)
	@echo ============================ Remaking $@
	@$D/common/scripts/make_dependencies  $(INCLUDE) -B. woundHealing2D_V6.cpp \
		> $@ \
	|| (rm -f $@ ; false)
	@if test -s $@ ; then : else rm $@ ; fi


# To make the dependencies known to `make', we finally have to include
# them:
include Makefile.dep



