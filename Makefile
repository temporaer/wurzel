.PHONY: dijkstraf dijkstra all sato vis
all: sato vis
sato:
	python -O main.py data/L2_17aug.dat
dijkstraf:
	make -C dijkstra runf
dijkstra:
	make -C dijkstra run
vis:
	ipython -wthread vis.py
