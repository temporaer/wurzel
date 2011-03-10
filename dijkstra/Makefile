LDFLAGS=-L boost_1_45_0/stage/lib -l boost_system -l boost_filesystem -l lapack -l boost_program_options -lxerces-c
CCFLAGS=-I $(HOME)/pool/include/boost-numeric-bindings/ -Wall -O3 -g0 -DNDEBUG -march=native -mtune=native --fast-math -fno-builtin -std=gnu++0x
BASE ?= data/L2_17aug

grid: grid.o config.o config-pskel.o
	g++-4.4 $(LDFLAGS) -o grid grid.o config.o config-pskel.o
grid.o: grid.cpp grid_graph.hpp voxelgrid.hpp wurzel_tree.hpp grid_config.hpp voxel_accessors.hpp voxel_normals.hpp
	g++-4.4 $(CCFLAGS) -c -o grid.o grid.cpp -I boost_1_45_0 
config-pskel.o: config-pskel.cxx config.xsd
	g++-4.4 $(CCFLAGS) -c -o config-pskel.o config-pskel.cxx 
config.o: config.cxx config-pskel.cxx config.xsd
	g++-4.4 $(CCFLAGS) -c -o config.o config.cxx 
config-pskel.cxx: config.xsd
	rm -f config-pskel.cxx config-pskel.hxx
	xsdcxx cxx-parser config.xsd
config.cxx: config.xsd
	rm -f config.hxx config.cxx
	xsdcxx cxx-tree   config.xsd

#grid: grid.cpp grid_graph.hpp voxelgrid.hpp wurzel_tree.hpp
#    mpic++ -Wall -O3 -g0 -DNDEBUG $(LDFLAGS) -o grid grid.cpp -I boost_1_45_0 
grid_dbg: grid.cpp grid_graph.hpp voxelgrid.hpp wurzel_tree.hpp
	g++ -Wall -O0 -g3 $(LDFLAGS) -o grid_dbg grid.cpp -I boost_1_45_0 
dbg: grid_dbg
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib cgdb ./grid_dbg ../$(BASE)
runf: grid
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib time ./grid ../$(BASE) --force
run: grid
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib time ./grid ../$(BASE)

test: grid
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./grid ../data/GersteLA_128x128x410
gg:gg.cpp grid_graph.hpp
	g++ -Wall -O0 -g3 $(LDFLAGS) -o gg gg.cpp -I boost_1_45_0 

clean:
	rm -f *.o
	rm -f config-pskel.cxx config-pskel.hxx config.cxx config.hxx
