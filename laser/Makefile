CCFLAGS=-I.
LDFLAGS=-llapack -lcv
points2grid: points2grid.cpp icp.hpp arg_max.hpp detect_circle.hpp
	g++-4.4 -Wall -O3 -g0 -DNDEBUG $(CCFLAGS) $(LDFLAGS) -march=native -mtune=native --fast-math -fno-builtin -std=gnu++0x -o points2grid points2grid.cpp
run:points2grid
	time ./points2grid
