# common
SDIR     = .

# Eric's Machine
CPP    = g++
LDFLAGS   = -lm -framework Accelerate
CPPFLAGS = -O1 -g -c -I$(SDIR)
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
