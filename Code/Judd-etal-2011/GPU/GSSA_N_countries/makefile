# Keeneland with the following modules:
# Currently Loaded Modulefiles:
#   1) modules                4) gold                   7) openmpi/1.5.1-intel
#     2) torque/2.5.11          5) mkl/2011_sp1.8.273     8) PE-intel
#       3) moab/6.1.5             6) intel/2011_sp1.8.273   9) cuda/4.1
#
CPUCC	 = icpc
GPUCC 	 = nvcc --ptxas-options=-v 
IMKL	 = /opt/intel/composer_xe_2011_sp1.8.273/mkl/include
LMKL	 = /opt/intel/composer_xe_2011_sp1.8.273/mkl/lib/intel64
ICUDA	 = /sw/keeneland/cuda/4.1/linux_binary/include
LCUDA	 = /sw/keeneland/cuda/4.1/linux_binary/lib64
LFLAGS	 = -L$(LCUDA) -lcudart -L$(LMKL) -lm -lpthread -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
CPUFLAGS = -O1 -g -c -I$(IMKL)
GPUFLAGS = -O1 -m64 -g -c -I$(IMKL) -I$(ICUDA) --gpu-architecture sm_20
HEADERS  = global.h ran.h
OBJECTS  = Timer.o Productivity.o GH_Quadrature.o Monomials_1.o Monomials_2.o Ord_Polynomial_N.o Num_Stab_Approx.o Accuracy_Test_N.o auxfuncs.o

Main_GSSA_N :	Main_GSSA_N.o $(OBJECTS) $(HEADERS)
		$(GPUCC) -m64 -o Main_GSSA_N Main_GSSA_N.o $(OBJECTS) $(LFLAGS) 

Main_GSSA_N.o :		Main_GSSA_N.cu $(HEADERS)
		$(GPUCC) $(GPUFLAGS) Main_GSSA_N.cu

Timer.o :	Timer.cpp $(HEADERS)
		$(CPUCC) $(CPUFLAGS) Timer.cpp

Productivity.o :	Productivity.cpp $(HEADERS)
		$(CPUCC) $(CPUFLAGS) Productivity.cpp

GH_Quadrature.o :	GH_Quadrature.cpp $(HEADERS)
		$(CPUCC) $(CPUFLAGS) GH_Quadrature.cpp

Monomials_1.o :	Monomials_1.cpp $(HEADERS)
		$(CPUCC) $(CPUFLAGS) Monomials_1.cpp

Monomials_2.o :	Monomials_2.cpp $(HEADERS)
		$(CPUCC) $(CPUFLAGS) Monomials_2.cpp

Ord_Polynomial_N.o :	Ord_Polynomial_N.cpp $(HEADERS)
		$(CPUCC) $(CPUFLAGS) Ord_Polynomial_N.cpp

Num_Stab_Approx.o :	Num_Stab_Approx.cpp $(HEADERS)
		$(CPUCC) $(CPUFLAGS) Num_Stab_Approx.cpp

Accuracy_Test_N.o :	Accuracy_Test_N.cpp $(HEADERS)
		$(CPUCC) $(CPUFLAGS) Accuracy_Test_N.cpp

auxfuncs.o :	auxfuncs.cpp $(HEADERS)
		$(CPUCC) $(CPUFLAGS) auxfuncs.cpp

clean :
	rm -f *.o
	rm -f core core.*

veryclean :
	rm -f *.o
	rm -f core core.*
	rm -f *.linkinfo
	rm -f Main_GSSA_N
