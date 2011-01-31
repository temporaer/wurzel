LDFLAGS=-L boost_1_45_0/stage/lib -l boost_system -l boost_filesystem
grid: grid.cpp
	g++ -O3 -g0 -DNDEBUG $(LDFLAGS) -o grid grid.cpp -I boost_1_45_0 
run: grid
	LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./grid
	#neato -Tpng  dijkstra-eg.dot > bla.png
	#gqview bla.png
