//=======================================================================
// Copyright 2009 Trustees of Indiana University.
// Authors: Michael Hansen, Andrew Lumsdaine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================

#include <iostream>
#include <fstream>
#include <map>
#include <boost/array.hpp>
#include <boost/graph/grid_graph.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#define DIMENSIONS 2
using namespace boost;

typedef grid_graph<DIMENSIONS> graph_t;
typedef graph_traits<graph_t> Traits;
typedef Traits::edge_descriptor edge_descriptor;
typedef Traits::vertex_descriptor vertex_descriptor;

struct edge_weight_map;
// ReadablePropertyGraph associated types
namespace boost {
  template<>
  struct property_map< graph_t, edge_weight_t > {
    typedef edge_weight_map type;
    typedef edge_weight_map const_type;
  };
}
/*
   Map from edges to weight values
*/
struct edge_weight_map {
  typedef double value_type;
  typedef value_type reference;
  typedef edge_descriptor key_type;
  typedef boost::readable_property_map_tag category;
  const graph_t& m_graph;
  edge_weight_map(const graph_t& g)
  :m_graph(g)
  {
  }

  // Edges have a weight equal to the average of their endpoint indexes.
  reference operator[](key_type e) const;
};

// Use these propety_map and property_traits parameterizations to refer to
// the associated property map types.
typedef boost::property_map<graph_t, boost::edge_weight_t>::const_type
        const_edge_weight_map;
typedef boost::property_traits<const_edge_weight_map>::reference
        edge_weight_map_value_type;
typedef boost::property_traits<const_edge_weight_map>::key_type
        edge_weight_map_key;

namespace boost{
	// PropertyMap valid expressions
	edge_weight_map_value_type
		get(const_edge_weight_map pmap, edge_weight_map_key e) {
			return pmap[e];
		}
	// ReadablePropertyGraph valid expressions
	const_edge_weight_map get(boost::edge_weight_t, const graph_t&g) {
		return const_edge_weight_map(g);
	}
	edge_weight_map_value_type get(boost::edge_weight_t tag,
			const graph_t& g,
			edge_weight_map_key e) {
		return get(tag, g)[e];
	}
}

  edge_weight_map::reference 
  edge_weight_map::operator[](key_type e) const {
	double d = get(vertex_index, m_graph)[target(e,m_graph)];
	//double d = get(vertex_index, m_graph)[e.first]-get(vertex_index, m_graph)[e.second];
    //return fabs(d);
    //return (get(edge_index,e.first + e.second)/2.0;
    //return 2.0;
  }

int main(int argc, char* argv[]) {

  // Define a 3x5x7 grid_graph where the second dimension doesn't wrap
  boost::array<std::size_t, 2> lengths = { { 8, 8 } };
  graph_t graph(lengths, false); // no dim is wrapped

  //shared_array_property_map<edge_weight_t,
  //                          property_map<graph_t, edge_index_t>::const_type>
  //                          weight(num_edges(graph), get(edge_index, graph));

  // Start with the first vertex in the graph
  Traits::vertex_descriptor first_vertex = vertex(0, graph);
  //print_vertex(first_vertex); // prints "(0, 0, 0)"

  const_edge_weight_map m = get(edge_weight, graph);

  shared_array_property_map<vertex_descriptor,
                            property_map<graph_t, vertex_index_t>::const_type>
                            p_map(num_vertices(graph), get(vertex_index, graph)); 
  shared_array_property_map<double,
                            property_map<graph_t, vertex_index_t>::const_type>
                            d_map(num_vertices(graph), get(vertex_index, graph)); 

  dijkstra_shortest_paths(graph, first_vertex
		  ,predecessor_map(p_map)
		  .distance_map(d_map)
		  );


  property_map<graph_t, vertex_index_t>::type name = get(vertex_index, graph);


  std::cout << "distances and parents:" << std::endl;
  Traits::vertex_iterator vi, vend;
  for (tie(vi, vend) = vertices(graph); vi != vend; ++vi) {
	  std::cout << "distance(" << name[*vi] << ") = " << d_map[*vi] << ", ";
	  std::cout << "parent(" << name[*vi] << ") = " << name[p_map[*vi]] << std::
		  endl;
  }
  std::cout << std::endl;

  std::ofstream dot_file("dijkstra-eg.dot");

  dot_file << "digraph D {\n"
	  << "  rankdir=LR\n"
	  << "  size=\"40,30\"\n"
	  << "  ratio=\"fill\"\n"
	  << "  edge[style=\"bold\"]\n" << "  node[shape=\"circle\"]\n";

  Traits::edge_iterator ei, ei_end;
  for (tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) {
	  graph_traits < graph_t >::edge_descriptor e = *ei;
	  graph_traits < graph_t >::vertex_descriptor
		  u = source(e, graph), v = target(e, graph);
	  dot_file << name[u] << " -> " << name[v]
		  << "[label=\"" << m[e] << "\"";
	  if (p_map[v] == u)
		  dot_file << ", color=\"black\"";
	  else
		  dot_file << ", color=\"grey\"";
	  dot_file << "]";
  }
  dot_file << "}";
}
