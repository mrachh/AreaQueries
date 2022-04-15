# user choose platform...
#HOST = osx
#HOST=linux-gfortran
HOST=linux-gfortran-openmp

ifeq ($(HOST),osx)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -c -w
FLINK = gfortran -w -o $(EXEC)
FEND = -framework accelerate
endif

ifeq ($(HOST),linux-gfortran)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math -c -w -std=legacy
FLINK = gfortran -w -o $(EXEC) 
FEND = -lblas -llapack
endif

ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math --openmp -c -w -std=legacy
FLINK = gfortran --openmp -w -o $(EXEC) 
FEND = -lblas -llapack
endif

ifeq ($(HOST),linux-ifort)
FC = ifort
FFLAGS = -O3 -c -w -xW 
FLINK = ifort -w -mkl -o $(EXEC)
WITH_SECOND = 1
endif


COM = ../src/common
3D = ../src/3d


.PHONY: all clean list

SOURCES =  test3d.f \
  $(COM)/prini_new.f \
  $(COM)/hkrand.f \
  $(COM)/sort.f \
  $(3D)/findnear3d.f \
  $(3D)/tree_lr_3d.f \
  $(COM)/dlaran.f 

ifeq ($(WITH_SECOND),1)
SOURCES += $(HELLSKITCHEN)/Common/second-r8.f
endif

OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))
EXEC = test3d

#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f
	$(FC) $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

default: run

run: $(EXEC)
	./$(EXEC)

$(EXEC): $(OBJECTS)
	$(FLINK) $(OBJECTS) $(FEND)

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



