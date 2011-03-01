LDFLAGS=-L boost_1_45_0/stage/lib -l boost_system -l boost_filesystem -l lapack
BASE ?= data/L2_17aug
grid: grid.cpp grid_graph.hpp voxelgrid.hpp wurzel_tree.hpp
	g++-4.4 -Wall -O3 -g0 -DNDEBUG $(LDFLAGS) -march=native -mtune=native --fast-math -fno-builtin -std=gnu++0x -o grid grid.cpp -I boost_1_45_0 
#grid: grid.cpp grid_graph.hpp voxelgrid.hpp wurzel_tree.hpp
#    mpic++ -Wall -O3 -g0 -DNDEBUG $(LDFLAGS) -o grid grid.cpp -I boost_1_45_0 
grid_dbg: grid.cpp grid_graph.hpp voxelgrid.hpp wurzel_tree.hpp
	g++ -Wall -O0 -g3 $(LDFLAGS) -o grid_dbg grid.cpp -I boost_1_45_0 
dbg: grid_dbg
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib cgdb ./grid_dbg ../$(BASE)
runf: grid
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib time ./grid ../$(BASE) force
run: grid
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib time ./grid ../$(BASE)

test:gg
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./gg
gg:gg.cpp grid_graph.hpp
	g++ -Wall -O0 -g3 $(LDFLAGS) -o gg gg.cpp -I boost_1_45_0 
