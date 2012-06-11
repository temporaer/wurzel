# Purpose

This project was created to analyze 3D MRI root data.

<img width="100" src="https://github.com/temporaer/wurzel/raw/master/images/barley-raw.jpg"/>
<img width="100" src="https://github.com/temporaer/wurzel/raw/master/images/barley-model.jpg"/>

# Guide to the files

## Vesselness measure and visualization

- `config.xml` contains meta-information about raw data files
- `main.py`    takes raw data, upsamples it and calculates vesselness measure
- `wurzel/*`   helper classes for `main.py` and `vis.py`

## Root tree generation and operations on root tree

- `dijkstra/grid.cpp`     transforms raw data and vesselness measure to root tree
- `dijkstra/treeinfo.cpp` print info about a root tree, translate to other formats

## Visualization

- `wurzel/viewer.py`  defines some default views on 3D data for convenience
- `vis.py`            visualizes results of all other components, including c++ components using `viewer.py`

# Running

## Meta-Data

For each raw data file, put a section like this in the config file:

	<datafile read-dtype="float32">
		<base-name>../data/GersteLA_96x96x410</base-name>
		<read-shape x="410" y="96" z="96"/>
		<shape x="436" y="96" z="96"/>
		<stem-plane>5</stem-plane>
		<stem-axis>0</stem-axis>
		<has-rohr>false</has-rohr>
		<scale>0.94</scale>
		<noise-cutoff>0.02</noise-cutoff>
	</datafile>

- `read-dtype` is the datatype of the elements
- `base-name` is the basename (with prepended `../`) as used lateron in BASE, which is used as a key
- `read-shape` the size of the data as stored in the raw `.dat` file
- `shape`     the size of the data after upsampling (should be isotropic then!)
- `stem-plane` index of plane in which to look for stem
- `stem-axis` number of axis the plane index refers to
- `has-rohr`  true/false, whether the dataset contains a measuring tube
- `scale`     how to get from `read-shape` resololution to millimeters (factor)
- `noise-cutoff` average noise level after normalization

## Run everything at once

- `Makefile` contains targets to run all steps separately.

  `make BASE=GersteLA_96x96x410_normal sato`

  upsamples `data/GersteLA_96x96x410_normal` and computes vesselness

  `make BASE=GersteLA_96x96x410_normal grid`

  given vesselness and upsampled raw data from previous step, find the root tree (in the graph-theoretical sense ^^)

  `make BASE=GersteLA_96x96x410_normal treeinfo`

  output some stats about the found root tree

- `Makefile` has targets for running a dataset through the whole process, e.g.

  `make BASE=GersteLA_96x96x410_normal`

  does all of the above steps.

# Parameters you can choose

- *Note* that when running everything through `make` commands as suggested above,
  you need to adjust the `dijkstra/grid` parameters in the `Makefile`.

- vesselness measure: You can mainly configure the scales at which the measure
  is calculated. Have a look at `main.py` for that purpose.

- root tree extraction: The `grid` program is configurable using command line parameters.
  Try running `grid --help` to find out more, e.g.

		Allowed options:
			--help                                produce help message
			--base                                the base name of the dataset
			-c [ --cfg ] arg (=config.xml)        the config-file containing dataset
			                                      descriptions (XML)
			--force                               force recomputation of dijkstra
			                                      algorithm
			--stem-plane arg                      plane index in which to search for stem
			                                      (read from config-file if not given)
			--stem-axis arg                       the axis of stem-plane
			                                      (read from config-file if not given)
			-v [ --max-void-dist ] arg (=5)       maximum consecutive distance (in mm) 
			                                      traveled through below-noiselevel data 
			                                      (use this to get rid of weed not 
			                                      connected to root)
			-r [ --max-radius ] arg (=1.8)        maximum root radius (in mm)
			-s [ --start-threshold ] arg (=0.1)
			                                      minimum raw value to start tracking
			-t [ --total-len-thresh ] arg (=1000000000)
			                                      maximal total length
			-d [ --dijkstra-stop-val ] arg (=1000000000)
			                                      stop dijkstra when paths longer than
			                                      this (decreases dijkstr runtime)
			-l [ --leaf-select-method ] arg (=median_raw)
			                                      method for selecting leaf candidates
			                                      [edge_detect,median_raw,subtree_weight]
			-f [ --min-flow-thresh ] arg (=0.0001)
			                                      minimal flow fraction
			--no-gauss-fit                        save some time
			--no-subpix-pos                       save some time
	
- Note to some of the parameters:

  - `force`: after Dijkstra algorithm has been executed, the result is written
    to files `d_map.dat` and `p_map.dat`. When you re-run `grid`, these files
    will be loaded instead of re-running Dijkstra. After changing scales in the
    vesselness measure or `dijkstra-stop-val`, you have to use this parameter to 
    ensure that the distance/predecessor maps are updated.
  - `dijskstra-stop-val`: 
    Most of the MRI image volume does not belong to the root. You can avoid
    determining its connectivity (which takes a lot of time) by choosing a maximum
    distance from the root you want to explore.
  - `max-void-dist`: Cut off branches where you need to cross this many
    millimeters of below-noise-level space. E.g. if you have underground weed
    that is not connected to the root, but the root is largely connected,
    you can use this parameter to remove the weed.
  - `leaf-select-method`:
    a voxel is considered to be a leaf node if it is above `start-threshold` *and*
    its distance to the root node is less than `total-len-thresh` *and*
      - `edge_detect` is the method used in the VISAPP paper for maize.
        it is true if the upstream/downstream ratio is above `min-flow-thresh`.
      - `subtree_weight`: In the subtree defined by the current node, the sum of
        the weights above threshold must be larger than `min-flow-thresh`.
      - `median_raw` In a breadth first search towards the root, and away from
        it (with same number of steps, respectively) the median of the visited
        raw-values must at least be `min-flow-thresh`.


# Notes

## Intermediate Results

All steps save intermediate outputs in the same directory as the original raw data.
Subsequent runs will make use of (at least some) cached results:

- upsampling is not repeated if avoidable
- dijkstra is only run if no d_map/p_map file exists
- treeinfo is very fast and operates only on serialized root tree

If in doubt, delete intermediate files and re-run everything.

## Parallelization

Care has been taken to parallelize as many operations as possible:

- Hessian computation of a *single* scale is parallelized with OpenMP using
  weave in Python. It will make use of all available cores.
- Optionally, you can distribute the Hessian computation of *different* scales
  over multiple computers using a running ipcluster engine and the `-p`
  parameter to `main.py`
- Dijkstra is not distributed, the parallel implementation for
  `dijkstra_shortest_path` algorithm is currently only available for adjacency
  list graphs, not the new `grid_graph`s.
- Most loops over vertices in `grid.cpp` have been parallelized using
  Intel Thread Building blocks. It automatically makes use of all available
  cores. The side-effect free action objects are implemented using
  C++0X lambda functions, which are available in `gcc` version >=4.5.

# Citing

If you use this code, please cite our publication

	Hannes Schulz, Johannes Postma, Dagmar van Dusschoten, Hanno Scharr, and Sven Behnke:
	3D Reconstruction of Plant Roots from MRI Images
	In Proceedings of International Conference on Computer Vision Theory and Applications (VISAPP), Rome, February 2012.
