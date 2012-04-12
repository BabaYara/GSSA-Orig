
# Eric's Machine
#SDIR     = .
#CPP = g++
#LDFLAGS   = -lm -framework Accelerate
#CPPFLAGS = -O2 -c -I$(SDIR)

# Keeneland with the following modules:
# Currently Loaded Modulefiles:
#   1) modules                4) gold                   7) openmpi/1.5.1-intel
#     2) torque/2.5.11          5) mkl/2011_sp1.8.273     8) PE-intel
#       3) moab/6.1.5             6) intel/2011_sp1.8.273   9) cuda/4.1
#
CPP  = icpc

# Use nvcc to compile .cu files
NVCC = nvcc
NVCCFLAGS = -arch sm_20 # For fermi's in keeneland

# Add MKL Paths
MKL_INC_PATH=/opt/intel/composer_xe_2011_sp1.8.273/mkl/include
MKL_LIB_PATH=/opt/intel/composer_xe_2011_sp1.8.273/mkl/lib/intel64

# Add MKL to the link line
LDFLAGS   = -L$(MKL_LIB_PATH) -lm -lpthread -openmp -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
#LDFLAGS   = -L$(MKL_LIB_PATH) -lm -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

# Add CUDA libraries to the link line
LDFLAGS += -lcuda -lcudart -L/sw/keeneland/cuda/4.1/linux_binary/lib64

# Include standard optimization flags and MKL headers
CPPFLAGS = -O2 -g -c -I$(MKL_INC_PATH)

# List of all the objects you need
OBJECTS  = Timer.cpp Productivity.o GH_Quadrature.o Monomials_1.o Monomials_2.o Ord_Polynomial_N.o Num_Stab_Approx.o Accuracy_Test_N.o auxfuncs.o

# Rule that tells make how to make the program from the objects
Main_GSSA_N :	Main_GSSA_N.o $(OBJECTS)
		$(CPP) -o Main_GSSA_N Main_GSSA_N.o $(OBJECTS) $(LDFLAGS) 

# Rule that tells make how to turn a .cu file into a .o
%.o: %.cu
		$(NVCC) ${NVCCFLAGS} $(CPPFLAGS) -c $<

# How does make know how to turn a .cpp into a .o?  It's built-in!
# but if you wanted to type it out it would look like:
# %.o: %.cpp
# 	$(CPP) $(CPPFLAGS) -c $<

clean :
	rm -f *.o
	rm -f core core.*

veryclean :
	rm -f *.o
	rm -f core core.*
	rm -f Main_GSSA_N
