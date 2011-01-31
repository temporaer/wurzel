#ifndef __WURZEL_TREE_HPP__
#define __WURZEL_TREE_HPP__
#include <boost/graph/adjacency_list.hpp>
#include "voxelgrid.hpp"

typedef boost::adjacency_list<
	boost::listS,
	boost::listS,
	boost::bidirectionalS,
	boost::property<boost::vertex_name_t,voxel_vertex_descriptor>
	> wurzelgraph_t;

typedef graph_traits<wurzelgraph_t> wurzelg_traits;
typedef wurzelg_traits::edge_descriptor wurzel_edge_descriptor;
typedef wurzelg_traits::edge_iterator   wurzel_edge_iterator;
typedef wurzelg_traits::in_edge_iterator   wurzel_in_edge_iterator;
typedef wurzelg_traits::vertex_descriptor wurzel_vertex_descriptor;
typedef wurzelg_traits::vertex_iterator   wurzel_vertex_iterator;


#endif
