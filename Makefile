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
imgs:
	ipython -wthread -- vis.py $(BASE) -o

demo1:
	ipython -wthread -- vis.py -t us-vs-raw    $(BASE) 
demo2:
	ipython -wthread -- vis.py -t us-vs-ground $(BASE) 
demo3:
	ipython -wthread -- vis.py -t us-vs-ground-vs-raw $(BASE) 
clean:
	find  -maxdepth 3 -name '*.pyc' | xargs -i rm {}
	find  -maxdepth 3 -name '*.pyo' | xargs -i rm {}
