# Use icpc to compile .cpp files
CPP  = icpc

# Use nvcc to compile .cu files
NVCC = nvcc
NVCCFLAGS = -arch sm_20 # For fermi's in keeneland

# Add MKL Paths
IMKL = /opt/intel/composer_xe_2011_sp1.8.273/mkl/include
LMKL = /opt/intel/composer_xe_2011_sp1.8.273/mkl/lib/intel64

# Add MAMGA Paths
IMAGMA = /sw/keeneland/magam/1.1/centos5.5_gnu4.4/magma_1.1.0/include
LMAGMA = /sw/keeneland/magam/1.1/centos5.5_gnu4.4/magma_1.1.0/lib

# Add MKL to the link line
LFLAGS = -L$(LMKL) -lm -lpthread -openmp -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

# Add CUDA libraries to the link line
LFLAGS += -lcuda -lcudart -L/sw/keeneland/cuda/4.1/linux_binary/lib64

# Add MAGMA libraries to the link line
LFLAGS += -L$(LMAGMA) -lmagma -lmagmablas

# Include standard optimization flags and MKL headers
CPPFLAGS = -O1 -g -c -I$(IMKL) -I$(IMAGMA)

# List of all the objects you need
OBJECTS  = Timer.o Productivity.o GH_Quadrature.o Monomials_1.o Monomials_2.o Ord_Polynomial_N.o Num_Stab_Approx.o Accuracy_Test_N.o auxfuncs.o

# Rule that tells make how to make the program from the objects
Main_GSSA_N :	Main_GSSA_N.o $(OBJECTS)
		$(CPP) -o Main_GSSA_N Main_GSSA_N.o $(OBJECTS) $(LFLAGS) 

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
