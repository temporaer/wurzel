LDFLAGS=-L boost_1_45_0/stage/lib -l boost_system -l boost_filesystem
grid: grid.cpp
	g++ -Wall -O3 -g0 -DNDEBUG $(LDFLAGS) -o grid grid.cpp -I boost_1_45_0 
grid_dbg: grid.cpp
	g++ -Wall -O0 -g3 $(LDFLAGS) -o grid_dbg grid.cpp -I boost_1_45_0 
dbg: grid_dbg
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib cgdb ./grid_dbg
run: grid
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./grid
	#neato -Tpng  dijkstra-eg.dot > bla.png
	#gqview bla.png
