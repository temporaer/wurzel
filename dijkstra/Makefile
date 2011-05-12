# Release-Mode or Debug-Mode?
#OPTFLAGS=-O3 -g0 -DNDEBUG -march=native -mtune=native --fast-math -fno-builtin
OPTFLAGS=-O0 -g3

INCLUDE=-I $(HOME)/pool/include/boost-numeric-bindings/ -I boost_1_45_0 -I inc -I.
LDFLAGS=-L boost_1_45_0/stage/lib -l boost_system -l boost_filesystem -l lapack -l boost_program_options -lxerces-c -lblas -lboost_serialization -ltbb
CCFLAGS=-Wall -std=gnu++0x $(OPTFLAGS) $(INCLUDE)
BASE ?= data/L2_17aug
GCC ?= g++-4.5

all:grid treeinfo

# ##################################
#   Compile targets
# ##################################
grid: grid.o rprop.o config.o config-pskel.o 
	$(GCC) -o grid grid.o config.o config-pskel.o rprop.o $(LDFLAGS) 
grid.o: grid.cpp inc/grid_graph.hpp inc/voxelgrid.hpp inc/wurzel_tree.hpp inc/grid_config.hpp inc/voxel_accessors.hpp inc/voxel_normals.hpp config.hxx config.xsd
	$(GCC) $(CCFLAGS) -c -o grid.o grid.cpp
rprop.o: rprop.cpp inc/rprop.hpp
	$(GCC) $(CCFLAGS) -c -o rprop.o rprop.cpp

# config XML file parser
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

treeinfo: treeinfo.o config.o config-pskel.o
	$(GCC) -o treeinfo treeinfo.o config.o config-pskel.o $(LDFLAGS)
treeinfo.o:treeinfo.cpp inc/treeinfo_config.hpp inc/wurzel_info.hpp inc/wurzel_tree.hpp
	$(GCC) $(CCFLAGS) -c -o treeinfo.o treeinfo.cpp

# ##################################
#   Cleanup
# ##################################
clean:
	rm -f *.o
	rm -f config-pskel.cxx config-pskel.hxx config.cxx config.hxx
	rm -f treeinfo
	rm -f grid
