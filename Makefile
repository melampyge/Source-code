.PHONY: all allclean clean

CC=icc
CXX=icc

OBJS= main.o bc_func.o verlet.o io_func.o model.o divdeath.o pressure.o shutdown.o verletlist.o basics.o mtrand.o debug.o localstress.o

LD=icc

CFLAGS= -fp-model precise -openmp $(INCLUDES)
CXXFLAGS= $(CFLAGS) 
LDFLAGS= -openmp
#LDFLAGS=

all: cell_dpd

#verlet.o: verlet.cpp
#	$(CXX) $(CXXFLAGS) -c -openmp -o verlet.o verlet.cpp
#	$(CXX) $(CXXFLAGS) -c -o verlet.o verlet.cpp

#pressure.o: pressure.cpp
#	$(CXX) $(CXXFLAGS) -c -openmp -o pressure.o pressure.cpp

cell_dpd: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS) 

debug: CXX += -g -pg
debug: CC += -g -pg
debug: LDFLAGS += -pg
debug: cell_dpd
	
clean:
	rm -f *.o lock.dat

allclean:
	rm -f *.o cell_dpd lock.dat *~ *.dat *.xyz

new: | allclean all clean

newd: | allclean debug clean

