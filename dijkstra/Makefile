grid: grid.cpp
	g++ -o grid grid.cpp -I boost_1_45_0
run: grid
	./grid
	neato -Tpng  dijkstra-eg.dot > bla.png
	gqview bla.png
