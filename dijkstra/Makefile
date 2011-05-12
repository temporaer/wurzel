LDFLAGS=-L boost_1_45_0/stage/lib -l boost_system -l boost_filesystem -l lapack -l boost_program_options -lxerces-c -lblas -L ../../../plain/levmar-2.5 -llevmar -lboost_serialization -ltbb
CCFLAGS=-I $(HOME)/pool/include/boost-numeric-bindings/ -Wall -O3 -g0 -DNDEBUG -march=native -mtune=native --fast-math -fno-builtin -std=gnu++0x -I boost_1_45_0
#CCFLAGS=-I $(HOME)/pool/include/boost-numeric-bindings/ -Wall -O1 -g3 -march=native -mtune=native --fast-math -fno-builtin -std=gnu++0x
BASE ?= data/L2_17aug
GCC ?= g++-4.5

all:grid treecompare

grid: grid.o rprop.o config.o config-pskel.o gaussian_fit.o 
	$(GCC) -o grid grid.o config.o config-pskel.o gaussian_fit.o rprop.o $(LDFLAGS) 
grid.o: grid.cpp grid_graph.hpp voxelgrid.hpp wurzel_tree.hpp grid_config.hpp voxel_accessors.hpp voxel_normals.hpp config.hxx config.xsd
	$(GCC) $(CCFLAGS) -c -o grid.o grid.cpp
rprop.o: rprop.cpp rprop.hpp
	$(GCC) $(CCFLAGS) -c -o rprop.o rprop.cpp

gaussian_fit.o: gaussian_fit.cpp gaussian_fit.hpp
	$(GCC) $(CCFLAGS) -c -o gaussian_fit.o gaussian_fit.cpp


config-pskel.o: config-pskel.cxx config.xsd
	$(GCC) $(CCFLAGS) -c -o config-pskel.o config-pskel.cxx 
config.o: config.cxx config-pskel.cxx config.xsd
	$(GCC) $(CCFLAGS) -c -o config.o config.cxx 
config-pskel.cxx: config.xsd
	rm -f config-pskel.cxx config-pskel.hxx
	xsdcxx cxx-parser config.xsd
config.cxx config.hxx: config.xsd
	rm -f config.hxx config.cxx
	xsdcxx cxx-tree   config.xsd


treecompare: treecompare.o config.o config-pskel.o
	$(GCC) -o treecompare treecompare.o config.o config-pskel.o $(LDFLAGS)
treecompare.o:treecompare.cpp treecompare_config.hpp wurzel_info.hpp wurzel_tree.hpp
	$(GCC) $(CCFLAGS) -c -o treecompare.o treecompare.cpp

#grid: grid.cpp grid_graph.hpp voxelgrid.hpp wurzel_tree.hpp
#    mpic++ -Wall -O3 -g0 -DNDEBUG $(LDFLAGS) -o grid grid.cpp -I boost_1_45_0 
grid_dbg: grid.cpp grid_graph.hpp voxelgrid.hpp wurzel_tree.hpp
	$(GCC) -Wall -O0 -g3 $(LDFLAGS) -o grid_dbg grid.cpp -I boost_1_45_0 
dbg: grid_dbg
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib cgdb ./grid_dbg ../$(BASE)
runf: grid
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib time ./grid ../$(BASE) --force
run: grid
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib time ./grid ../$(BASE)

test: grid
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./grid ../data/GersteLA_128x128x410
gg:gg.cpp grid_graph.hpp
	$(GCC) -Wall -O0 -g3 $(LDFLAGS) -o gg gg.cpp -I boost_1_45_0 

clean:
	rm -f *.o
	rm -f config-pskel.cxx config-pskel.hxx config.cxx config.hxx
	rm -f treecompare
