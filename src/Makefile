CPP=g++-12
CPPFLAGS=-O3 -fopenmp  -std=gnu++17 

default: rt

all: rt

rt: rt.cpp Makefile
	$(CPP) $(CPPFLAGS) -o rt rt.cpp 
	./rt

clean:
	-rm rt