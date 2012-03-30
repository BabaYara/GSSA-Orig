CPUCC    = icc
SDIR     = .
HEADERS  = global.h ran.h svdcmp.h
OBJECTS  = Productivity.o GH_Quadrature.o Monomials_1.o Monomials_2.o Ord_Polynomial_N.o Num_Stab_Approx.o Accuracy_Test_N.o cholesky.o auxfuncs.o svdcmp.o qrdcmp.o simplx.o

# On my own machine (with GSL)
#IGSL    = /opt/local/include
#LGSL    = /opt/local/lib
#CPUFLAGS = -O2 -c -I$(SDIR) -I$(IGSL)
#LFLAGS   = -L$(LGSL) -lm -lgsl -lgslcblas

# On Keeneland (with MKL)
IMKL     = /opt/intel/composer_xe_2011_sp1.8.273/mkl/include
LMKL     = /opt/intel/composer_xe_2011_sp1.8.273/mkl/lib/intel64
CPUFLAGS = -g -O3 -I$(SDIR) -I$(IMKL)
LFLAGS   = -L $(MKL_LIB_PATH) -lpthread -openmp -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -Wl,--end-group


Main_GSSA_N :	Main_GSSA_N.o $(OBJECTS) $(HEADERS)
		$(CPUCC) -o Main_GSSA_N Main_GSSA_N.o $(OBJECTS) $(LFLAGS) 

Main_GSSA_N.o : 	$(SDIR)/Main_GSSA_N.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/Main_GSSA_N.cpp

Productivity.o : 	$(SDIR)/Productivity.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/Productivity.cpp

GH_Quadrature.o : 	$(SDIR)/GH_Quadrature.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/GH_Quadrature.cpp

Monomials_1.o : 	$(SDIR)/Monomials_1.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/Monomials_1.cpp

Monomials_2.o : 	$(SDIR)/Monomials_2.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/Monomials_2.cpp

Ord_Polynomial_N.o : 	$(SDIR)/Ord_Polynomial_N.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/Ord_Polynomial_N.cpp

Num_Stab_Approx.o : 	$(SDIR)/Num_Stab_Approx.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/Num_Stab_Approx.cpp

Accuracy_Test_N.o : 	$(SDIR)/Accuracy_Test_N.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/Accuracy_Test_N.cpp

cholesky.o : 	$(SDIR)/cholesky.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/cholesky.cpp

auxfuncs.o : 	$(SDIR)/auxfuncs.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/auxfuncs.cpp

svdcmp.o : 	$(SDIR)/svdcmp.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/svdcmp.cpp

qrdcmp.o : 	$(SDIR)/qrdcmp.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/qrdcmp.cpp

simplx.o : 	$(SDIR)/simplx.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/simplx.cpp

clean :
	rm -f *.o
	rm -f core core.*

veryclean :
	rm -f *.o
	rm -f core core.*
	rm -f Main_GSSA_N
