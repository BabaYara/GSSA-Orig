CPUCC    = g++
SDIR     = .
CPUFLAGS = -O2 -c -I$(SDIR)
LFLAGS   = -lm -framework Accelerate
HEADERS  = global.h ran.h svdcmp.h
OBJECTS  = Productivity.o GH_Quadrature.o Monomials_1.o Monomials_2.o Ord_Polynomial_N.o Num_Stab_Approx.o Accuracy_Test_N.o cholesky.o auxfuncs.o svdcmp.o qrdcmp.o simplx.o

Main_GSSA_N :	Main_GSSA_N.o $(OBJECTS) $(HEADERS)
		$(CPUCC) -o Main_GSSA_N Main_GSSA_N.o $(OBJECTS) $(LFLAGS) 

clean :
	rm -f *.o
	rm -f core core.*

veryclean :
	rm -f *.o
	rm -f core core.*
	rm -f Main_GSSA_N
