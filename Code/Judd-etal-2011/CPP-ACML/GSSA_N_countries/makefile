# common
SDIR     = .

# Keeneland with the following modules:
#Currently Loaded Modulefiles:
#  1) modules             3) moab/6.1.5          5) cuda/4.1            7) subversion/1.6.15   9) gcc/4.4.0
#  2) torque/2.5.11       4) gold                6) cmake/2.8.0         8) acml/4.1.0
CPP    = g++
ACML_INC_PATH=/sw/keeneland/acml/4.4.0/linux_binary/gfortran64/include
ACML_LIB_PATH=/sw/keeneland/acml/4.4.0/linux_binary/gfortran64/lib
LDFLAGS   = -L$(ACML_LIB_PATH) -lm -lacml
CPPFLAGS = -O1 -g -c -I$(ACML_INC_PATH) -I$(SDIR)
OBJECTS  = Timer.cpp Productivity.o GH_Quadrature.o Monomials_1.o Monomials_2.o Ord_Polynomial_N.o Num_Stab_Approx.o Accuracy_Test_N.o auxfuncs.o

Main_GSSA_N :	Main_GSSA_N.o $(OBJECTS)
		$(CPP) -o Main_GSSA_N Main_GSSA_N.o $(OBJECTS) $(LDFLAGS) 

clean :
	rm -f *.o
	rm -f core core.*

veryclean :
	rm -f *.o
	rm -f core core.*
	rm -f Main_GSSA_N
