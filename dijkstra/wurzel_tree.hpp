#ifndef __WURZEL_TREE_HPP__
#define __WURZEL_TREE_HPP__
#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "voxelgrid.hpp"

namespace ublas = boost::numeric::ublas;

struct path_orientation_t{
	typedef boost::vertex_property_tag kind;
} path_orientation;
struct marked_vertex_t{
	typedef boost::vertex_property_tag kind;
} marked_vertex;
typedef ublas::bounded_matrix<double,3,3,ublas::column_major> covmat_t;

typedef boost::adjacency_list<
	boost::listS,
	boost::listS,
	boost::bidirectionalS,
	boost::property<boost::vertex_name_t,voxel_vertex_descriptor,
	  boost::property<path_orientation_t,covmat_t,
		  boost::property<marked_vertex_t,bool
	  > > >
	> wurzelgraph_t;

typedef graph_traits<wurzelgraph_t> wurzelg_traits;
typedef wurzelg_traits::edge_descriptor wurzel_edge_descriptor;
typedef wurzelg_traits::edge_iterator   wurzel_edge_iterator;
typedef wurzelg_traits::in_edge_iterator   wurzel_in_edge_iterator;
typedef wurzelg_traits::vertex_descriptor wurzel_vertex_descriptor;
typedef wurzelg_traits::vertex_iterator   wurzel_vertex_iterator;


#endif
