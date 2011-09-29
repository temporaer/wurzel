.PHONY: dijkstra all vis treeinfo
DATAPATH=$(shell cat config.xml|grep datapath|perl -npe '($$_)=/>(.*)</')
BASE ?= L2_17aug
FILE  = $(DATAPATH)/$(BASE)
WHAT ?= mass
all: $(FILE)/sato.dat $(FILE)/wgraph.ser  treeinfo

# calculate tubularness measure
$(FILE)/upsampled.dat $(FILE)/sato.dat: $(FILE).dat
	mkdir -p $(FILE)
	python -O main.py $(BASE)

# calculate distance map, predecessor map, serialized graph
$(FILE)/d_map.dat $(FILE)/p_map.dat $(FILE)/wgraph.ser: $(FILE)/sato.dat $(FILE)/upsampled.dat
	make -C dijkstra grid
	echo "lupine settings!!!"
	#cd dijkstra && LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./grid $(BASE) --cfg=../config.xml -s 1 -f  5 -a 1.0 # Gerste
	#cd dijkstra && LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./grid $(BASE) --cfg=../config.xml -s 4 -f 16 -t 0.1 -a 0.4 # L2_6aug
	#cd dijkstra && LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./grid $(BASE) --cfg=../config.xml -s 4 -f  8 -t 0.2 -a 1.0 # L2_17aug
	#cd dijkstra && LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./grid $(BASE) --cfg=../config.xml -s 4 -f 8 -t 0.2 -a 1.0 # L2_17aug
	#cd dijkstra && LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./grid $(BASE) --cfg=../config.xml -s 4 -f 8 -t 0.2 -a 0.5 # L2_22aug
	cd dijkstra && LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./grid $(BASE) --cfg=../config.xml -s 4 -f 8 -t 1 -a 1.0 # Maize

# print some info about the graph
treeinfo: $(FILE)/wgraph.ser
	make -C dijkstra treeinfo
	cd dijkstra && LD_LIBRARY_PATH=boost_1_45_0/stage/lib ./treeinfo -b $(BASE) -c ../config.xml -a print

vis:
	ipython -wthread vis.py $(BASE)
imgs:
	ipython -wthread -- vis.py -t $(WHAT) $(BASE) -o

demo:
	#ipython -wthread -- vis.py -t $(WHAT)  $(BASE) 
	ipython -- vis.py -t $(WHAT)  $(BASE) 
clean:
	make -C presentation clean
	make -C dijkstra clean
	find  -maxdepth 3 -name '*.pyc' | xargs -i rm {}
	find  -maxdepth 3 -name '*.pyo' | xargs -i rm {}
