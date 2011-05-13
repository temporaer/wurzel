# Purpose

This project was created to analyze 3D MRI root data.


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

- `Makefile` has targets for running a dataset through the whole process, e.g.
	
	make BASE=data/L2_6aug

  upsamples data/L2_6aug.dat, computes vesselness, finds root tree and outputs statistics about it.
  BASE should be the relative path without extention to the raw .dat file.


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
