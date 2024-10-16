### Fortran Compiler
#FC = mpif90
#FC = ftn
FC = gfortran    #To run on PC _Sahar

### Compilation flags
#FFLAGS := -O3 -ffixed-line-length-none -fdefault-real-8 
FFLAGS := -g -O0 -ffixed-line-length-none -fdefault-real-8   #Adding debugging flags _Sahar

#BIG  := #-mcmodel=medium
#DBG  := #-g #-traceback
#PROF := #-pg
OMP  := -fopenmp

### Generate script-dependent target (exe)
num := $(shell grep -o 'case_num = .*' param.f90 | sed -n 's/case_num = //p' | tr -d '"/')
TARGET := HLGD$(num)

### Source scripts
#SRC = param.f90 common.f90 zero.f90 core.f90 output.f90 main.f90      Sahar
SRC = param.f90 common.f90 output.f90 zero.f90 core.f90 main.f90

OBJ = $(SRC:.f90=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
#	$(FC) $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) -o $@ $(OBJ) $(LIB)     Sahar
	$(FC) $(FFLAGS) $(OMP) -o $@ $(OBJ)    

.PHONY: clean  
clean:
	rm -f *.o *.mod $(TARGET) *~ .depend
	mkdir data1 data1/gap0 data1/para0 data1/force0 data1/para1 data1/restart_sim data1/post_processing		#!create folders to store data1  

%.o: %.f90
#	$(FC) $(INCLUDE) -c -o $@ $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) $<   Sahar
	$(FC) $(INCLUDE) -c -o $@ $(FFLAGS) $(OMP) $<

.depend dep:
	./.makedepo $(SRC) > .depend

include .depend
