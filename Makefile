.PHONY: dijkstra all vis treeinfo
DATAPATH=$(shell cat config.xml|grep datapath|perl -npe '($$_)=/>(.*)</')
BASE ?= L2_17aug
FILE  = $(DATAPATH)/$(BASE)
WHAT ?= mass
all: sato grid treeinfo
grid: $(FILE)/d_map.dat
sato: $(FILE)/sato.dat

# calculate tubularness measure
$(FILE)/upsampled.dat $(FILE)/sato.dat: $(FILE).dat
	mkdir -p $(FILE)
	python main.py $(BASE)

# calculate distance map, predecessor map, serialized graph
$(FILE)/d_map.dat $(FILE)/p_map.dat $(FILE)/wgraph.ser: $(FILE)/sato.dat $(FILE)/upsampled.dat
	cd dijkstra && ./grid $(BASE) --cfg=../config.xml -s 2 -f 1 --no-gauss-fit # Maize

# print some info about the graph
treeinfo: $(FILE)/wgraph.ser
	make -C dijkstra treeinfo
	cd dijkstra && ./treeinfo -b $(BASE) -c ../config.xml -a print

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
