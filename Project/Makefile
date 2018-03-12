CC   = gcc
CXX  = g++

CFLAGS =-std=c++11 -lfftw3 -lm  

EXEC = HMC

SRC_CPP = Main.cpp Fourier_Transform.cpp Monte_carlo.cpp stastics.cpp functions.cpp 

OBJS_CPP = $(SRC_CPP:.cpp=.o)


all: $(EXEC)

$(EXEC): $(OBJS_CPP)
	$(CXX) -o $(EXEC) $(OBJS_CPP)  $(CFLAGS)

clean:
	-rm *.o
	-rm *~
	-rm Make.log
	

distclean: clean
	-rm $(EXEC)

.SUFFIXES: .o .cpp

.cpp.o:
	$(CXX) -c $< $(CFLAGS)
