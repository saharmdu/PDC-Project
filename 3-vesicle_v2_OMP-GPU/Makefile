### Fortran Compiler
#FC = mpif90
#FC = ftn
#FC = gfortran    #To work on PC _Sahar
FC = /opt/nvidia/hpc_sdk/Linux_x86_64/24.9/compilers/bin/nvfortran

### Compilation flags
#FFLAGS := -O3 -ffixed-line-length-none -fdefault-real-8 
FFLAGS = -g -O0 -Mfree -Mextend -Mpreprocess -r8 -mp=gpu -gpu=cc89,mem:managed -Minfo=accel # Sahar

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
