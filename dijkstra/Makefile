LDFLAGS=-L boost_1_45_0/stage/lib -l boost_system -l boost_filesystem
grid: grid.cpp grid_graph.hpp voxelgrid.hpp wurzel_tree.hpp
	g++ -Wall -O3 -g0 -DNDEBUG $(LDFLAGS) -o grid grid.cpp -I boost_1_45_0 
grid_dbg: grid.cpp grid_graph.hpp voxelgrid.hpp wurzel_tree.hpp
	g++ -Wall -O0 -g3 $(LDFLAGS) -o grid_dbg grid.cpp -I boost_1_45_0 
dbg: grid_dbg
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib cgdb ./grid_dbg
runf: grid
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib time ./grid force
run: grid
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib time ./grid
	#neato -Tpng  dijkstra-eg.dot > bla.png
	#gqview bla.png

test:gg
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./gg
gg:gg.cpp grid_graph.hpp
	g++ -Wall -O0 -g3 $(LDFLAGS) -o gg gg.cpp -I boost_1_45_0 
