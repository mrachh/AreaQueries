
# # This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create
# the file make.inc, which overrides the defaults below (which are
# for ubunutu linux/gcc system).


# compiler, and linking from C, fortran
CC=gcc
CXX=g++
FC=gfortran


# set compiler flags for c and fortran
FFLAGS= -fPIC -O3 -march=native -funroll-loops -std=legacy -w
FFLAGS_DYN= -shared -fPIC
CFLAGS= -fPIC -O3 -march=native -funroll-loops -std=gnu17 
CXXFLAGS= -std=c++11 -DSCTL_PROFILE=-1 -fPIC -O3 -march=native -funroll-loops

# set linking libraries
CLIBS = -lgfortran -lm -ldl
LIBS = -lm

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp

# flags for MATLAB MEX compilation..
MFLAGS=-compatibleArrayDims -DMWF77_UNDERSCORE1 "CFLAGS=$(CFLAGS)"
MWFLAGS=-c99complex 
MOMPFLAGS = -D_OPENMP

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../mwrap/mwrap
MEXLIBS=-lm -lstdc++ -ldl -lgfortran

AQ_INSTALL_DIR=$(PREFIX)
ifeq ($(PREFIX),)
	AQ_INSTALL_DIR = ${HOME}/lib
endif

DYLIBS = $(LIBS)

LIBNAME=libareaqueries
DYNAMICLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a
LIMPLIB = $(DYNAMICLIB)

LLINKLIB = -lareaqueries



# For your OS, override the above by placing make variables in make.inc
-include make.inc


# multi-threaded libs & flags needed
ifneq ($(OMP),OFF)
  CFLAGS += $(OMPFLAGS)
  FFLAGS += $(OMPFLAGS)
  MFLAGS += $(MOMPFLAGS)
  LIBS += $(OMPLIBS)
  DYLIBS += $(OMPLIBS)
  MEXLIBS += $(OMPLIBS)
endif


LIBDIR = lib-static

# objects to compile
#
# Common objects
COM = src/common
COMOBJS = $(COM)/cumsum.o $(COM)/sort.o $(COM)/sparse_reps.o $(COM)/prini.o 

# 3D objects
3D = src/3d
3DOBJS = $(3D)/findnear3d.o $(3D)/pts_tree3d.o $(3D)/tree_routs3d.o 

# 2D objects
2D = src/2d
2DOBJS = $(2D)/findnear2d.o $(2D)/tree_lr_2d.o 

# Test objects
TOBJS = $(COM)/hkrand.o $(COM)/dlaran.o $(COM)/prini.o


OBJS = $(COMOBJS) $(3DOBJS) $(2DOBJS) 

.PHONY: usage install lib matlab mex

default: usage

usage:
	@echo "------------------------------------------------------------------------"
	@echo "Makefile for Area Queries. Specify what to make:"
	@echo ""
	@echo "  make install        compile and install the main library"
	@echo "  make install PREFIX=(INSTALL_DIR)  "
	@echo "                      compile and install the main library at custom"
	@echo "                      location given by PREFIX"
	@echo "  make lib            compile the main library (in lib/ and lib-static/)"
	@echo "  make matlab         compile matlab interfaces"
	@echo "  make objclean       removal all object files, preserving lib & MEX"
	@echo "  make clean          also remove lib, MEX, py, and demo executables"
	@echo "  make mex            generate matlab interfaces"
	@echo "                      (for expert users only, requires mwrap)"
	@echo ""
	@echo "For faster (multicore) making, append the flag -j"
	@echo " 'make [task] OMP=OFF'        for single-threaded"
	@echo "------------------------------------------------------------------------"


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@
%.o: %.f 
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

# build the library...
lib: $(STATICLIB) $(DYNAMICLIB)
ifneq ($(OMP),OFF)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif


install: $(STATICLIB) $(DYNAMICLIB)
	echo $(AQ_INSTALL_DIR)
	mkdir -p $(AQ_INSTALL_DIR)
	cp -f lib/$(DYNAMICLIB) $(AQ_INSTALL_DIR)/
	cp -f lib-static/$(STATICLIB) $(AQ_INSTALL_DIR)/
	[ ! -f lib/$(LIMPLIB) ] || cp -f lib/$(LIMPLIB) $(AQ_INSTALL_DIR)/
	@echo "Make sure to include " $(AQ_INSTALL_DIR) " in the appropriate path variable"
	@echo "    LD_LIBRARY_PATH on Linux"
	@echo "    PATH on windows"
	@echo "    DYLD_LIBRARY_PATH on Mac OSX (not needed if default installation directory is used"
	@echo " "
	@echo "In order to link against the dynamic library, use -L"$(AQ_INSTALL_DIR) " -lfmm3d"


$(STATICLIB): $(OBJS)
	ar rcs $(STATICLIB) $(OBJS)
	mv $(STATICLIB) lib-static/
$(DYNAMICLIB): $(OBJS)
	$(FC) $(FFLAGS_DYN) $(OBJS) -o $(DYNAMICLIB) $(DYLIBS)
	mv $(DYNAMICLIB) lib/
	[ ! -f $(LIMPLIB) ] || mv $(LIMPLIB) lib/


# matlab..
MWRAPFILE = areaqueries
GATEWAY = $(MWRAPFILE)

matlab:	$(STATICLIB) matlab/$(GATEWAY).c matlab/$(GATEWAY2).c
	$(MEX) matlab/$(GATEWAY).c lib-static/$(STATICLIB) $(MFLAGS) \
	-output matlab/areaqueries $(MEXLIBS)

mex:  $(STATICLIB)
	cd matlab; $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw;\
	$(MEX) $(GATEWAY).c ../lib-static/$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE) \
	$(MEXLIBS); \

objclean:
	rm -f $(OBJS) $(COBJS) $(TOBJS)
	rm -f test/*.o 
