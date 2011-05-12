.PHONY: dijkstra all vis treeinfo
BASE ?= data/L2_17aug
WHAT ?= mass
all: $(BASE).sato $(BASE)-wgraph.ser  treeinfo

# calculate tubularness measure
$(BASE)-upsampled.dat $(BASE).sato: $(BASE).dat
	python -O main.py $(BASE)

# calculate distance map, predecessor map, serialized graph
$(BASE)-d_map.dat $(BASE)-p_map.dat $(BASE)-wgraph.ser: $(BASE).sato $(BASE)-upsampled.dat
	make -C dijkstra grid
	cd dijkstra && LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./grid ../$(BASE) --cfg=../config.xml -s 4 -f 8

# print some info about the graph
treeinfo: $(BASE)-wgraph.ser
	make -C dijkstra treeinfo
	cd dijkstra && LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./treeinfo -b ../$(BASE) -c ../config.xml -a print

vis:
	ipython -wthread vis.py $(BASE)
imgs:
	ipython -wthread -- vis.py -t $(WHAT) $(BASE) -o

demo:
	ipython -wthread -- vis.py -t $(WHAT)  $(BASE) 
clean:
	find  -maxdepth 3 -name '*.pyc' | xargs -i rm {}
	find  -maxdepth 3 -name '*.pyo' | xargs -i rm {}
