# common
SDIR     = .

# Eric's Machine
#CPP = g++
#LDFLAGS   = -lm -framework Accelerate
#CPPFLAGS = -O2 -c -I$(SDIR)

# Keeneland with the following modules:
# Currently Loaded Modulefiles:
#   1) modules                4) gold                   7) openmpi/1.5.1-intel
#     2) torque/2.5.11          5) mkl/2011_sp1.8.273     8) PE-intel
#       3) moab/6.1.5             6) intel/2011_sp1.8.273   9) cuda/4.1
#
CPP    = icpc
MKL_INC_PATH=/opt/intel/composer_xe_2011_sp1.8.273/mkl/include
MKL_LIB_PATH=/opt/intel/composer_xe_2011_sp1.8.273/mkl/lib/intel64
LDFLAGS   = -L$(MKL_LIB_PATH) -lm -lpthread -openmp -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
#LDFLAGS   = -L$(MKL_LIB_PATH) -lm -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
CPPFLAGS = -O1 -g -c -I$(MKL_INC_PATH) -I$(SDIR)

OBJECTS  = Productivity.o GH_Quadrature.o Monomials_1.o Monomials_2.o Ord_Polynomial_N.o Num_Stab_Approx.o Accuracy_Test_N.o auxfuncs.o

Main_GSSA_N :	Main_GSSA_N.o $(OBJECTS)
		$(CPP) -o Main_GSSA_N Main_GSSA_N.o $(OBJECTS) $(LDFLAGS) 

clean :
	rm -f *.o
	rm -f core core.*

veryclean :
	rm -f *.o
	rm -f core core.*
	rm -f Main_GSSA_N
