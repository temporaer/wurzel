.PHONY: dijkstraf dijkstra all sato vis
BASE ?= data/L2_17aug
WHAT ?= mass
all: sato vis
sato:
	python -O main.py $(BASE)
dijkstraf:
	make -C dijkstra runf
dijkstra:
	make -C dijkstra run
vis:
	ipython -wthread vis.py $(BASE)
imgs:
	ipython -wthread -- vis.py $(BASE) -o

demo:
	ipython -wthread -- vis.py -t $(WHAT)  $(BASE) 
clean:
	find  -maxdepth 3 -name '*.pyc' | xargs -i rm {}
	find  -maxdepth 3 -name '*.pyo' | xargs -i rm {}
