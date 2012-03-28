CPUCC    = g++
SDIR     = .
CPUFLAGS = -O2 -c -I$(SDIR)
LFLAGS   = -lm
HEADERS  = simplx.h
OBJECTS  = simplx.o

main :	main.o $(OBJECTS) $(HEADERS)
		$(CPUCC) -o main main.o $(OBJECTS) $(LFLAGS) 

main.o : 	$(SDIR)/main.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/main.cpp

simplx.o : 	$(SDIR)/simplx.cpp  $(HEADERS)
		$(CPUCC) $(CPUFLAGS) $(SDIR)/simplx.cpp

clean :
	rm -f *.o
	rm -f core core.*

veryclean :
	rm -f *.o
	rm -f core core.*
	rm -f main
