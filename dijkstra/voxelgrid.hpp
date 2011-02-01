#ifndef __VOXELGRID_HPP__
#define __VOXELGRID_HPP__

#include <boost/graph/grid_graph.hpp>

#define DIMENSIONS 3
using namespace boost;


typedef unsigned int vidx_t;
typedef unsigned int eidx_t;

typedef grid_graph<DIMENSIONS,vidx_t,eidx_t> voxelgraph_t;
typedef graph_traits<voxelgraph_t> voxelg_traits;
typedef voxelg_traits::edge_descriptor voxel_edge_descriptor;
typedef voxelg_traits::edge_iterator   voxel_edge_iterator;
typedef voxelg_traits::vertex_descriptor voxel_vertex_descriptor;
typedef voxelg_traits::vertex_iterator   voxel_vertex_iterator;

struct voxel_edge_weight_map;
// ReadablePropertyGraph associated types
namespace boost {
  template<>
  struct property_map< voxelgraph_t, edge_weight_t > {
    typedef voxel_edge_weight_map type;
    typedef voxel_edge_weight_map const_type;
  };
}

/*
   Map from edges to weight values
*/
struct voxel_edge_weight_map {
  typedef double value_type;
  typedef value_type reference;
  typedef voxel_edge_descriptor key_type;
  typedef boost::readable_property_map_tag category;
  const voxelgraph_t& m_graph;
  voxel_edge_weight_map(const voxelgraph_t& g)
  :m_graph(g) { }

  // Edges have a weight equal to the average of their endpoint indexes.
  reference operator[](key_type e) const;
};

// Use these propety_map and property_traits parameterizations to refer to
// the associated property map types.
typedef boost::property_map<voxelgraph_t, boost::edge_weight_t>::const_type
        const_voxel_edge_weight_map;
typedef boost::property_traits<const_voxel_edge_weight_map>::reference
        voxel_edge_weight_map_value_type;
typedef boost::property_traits<const_voxel_edge_weight_map>::key_type
        edge_weight_map_key;

namespace boost{
	// PropertyMap valid expressions
	voxel_edge_weight_map_value_type get(const_voxel_edge_weight_map pmap, edge_weight_map_key e) {
		return pmap[e]; }
	// ReadablePropertyGraph valid expressions
	const_voxel_edge_weight_map get(boost::edge_weight_t, const voxelgraph_t&g) {
	   	return const_voxel_edge_weight_map(g); }
	voxel_edge_weight_map_value_type get(boost::edge_weight_t tag, const voxelgraph_t& g, edge_weight_map_key e) {
		return get(tag, g)[e]; }
}

#endif
