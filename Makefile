.PHONY: dijkstraf dijkstra all sato vis
BASE ?= data/L2_17aug
all: sato vis
sato:
	python -O main.py $(BASE)
dijkstraf:
	make -C dijkstra runf
dijkstra:
	make -C dijkstra run
vis:
	ipython -wthread vis.py $(BASE)
clean:
	find  -maxdepth 3 -name '*.pyc' | xargs -i rm {}
	find  -maxdepth 3 -name '*.pyo' | xargs -i rm {}
