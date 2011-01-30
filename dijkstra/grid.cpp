//=======================================================================
// Copyright 2009 Trustees of Indiana University.
// Authors: Michael Hansen, Andrew Lumsdaine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <boost/array.hpp>
#include <boost/graph/grid_graph.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/multi_array.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
using namespace boost::accumulators;


#define DIMENSIONS 3
using namespace boost;

typedef boost::multi_array_ref<unsigned char, 3> array_type;
typedef array_type::index index;

array_type* g_dat;

typedef unsigned int vidx_t;
typedef unsigned int eidx_t;

typedef grid_graph<DIMENSIONS,vidx_t,eidx_t> graph_t;
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
  :m_graph(g) { }

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
	edge_weight_map_value_type get(const_edge_weight_map pmap, edge_weight_map_key e) {
		return pmap[e]; }
	// ReadablePropertyGraph valid expressions
	const_edge_weight_map get(boost::edge_weight_t, const graph_t&g) {
	   	return const_edge_weight_map(g); }
	edge_weight_map_value_type get(boost::edge_weight_t tag, const graph_t& g, edge_weight_map_key e) {
		return get(tag, g)[e]; }
}

// map from edges to weights
edge_weight_map::reference 
edge_weight_map::operator[](key_type e) const {
	//double d = get(vertex_index, m_graph)[target(e,m_graph)];
	vertex_descriptor v = target(e,m_graph);
	array_type& a = *g_dat;
	float f = exp(-(float)a[v[0]][v[1]][v[2]]/25.6f);
	return f;
}

int main(int argc, char* argv[]) {


  std::ifstream dat0("../data/L2_22aug-upsampled.dat", std::ios::in | std::ios::binary);
  unsigned char* data0 = new unsigned char[256*256*256];
  dat0.read((char*)data0,256*256*256);
  array_type A0(data0,boost::extents[256][256][256]);
  dat0.close();

  std::ifstream dat("../data/L2_22aug-preproc.dat", std::ios::in | std::ios::binary);
  unsigned char* data = new unsigned char[256*256*256];
  dat.read((char*)data,256*256*256);
  array_type A(data,boost::extents[256][256][256]);
  g_dat = &A;
  dat.close();

  unsigned char* paths = new unsigned char[256*256*256];
  std::fill(paths, paths+256*256*256, (unsigned char) 0);
  array_type B(paths,boost::extents[256][256][256]);


  // Define a 3x5x7 grid_graph where the second dimension doesn't wrap
  boost::array<vidx_t, 3> lengths = { { 256, 256, 256 } };
  boost::array<vidx_t, 3> strunk  = { { 109, 129,  24 } };
  graph_t graph(lengths, false); // no dim is wrapped

  Traits::vertex_descriptor first_vertex = vertex((vidx_t)0, graph);

  shared_array_property_map<vertex_descriptor,
                            property_map<graph_t, vertex_index_t>::const_type>
                            p_map(num_vertices(graph), get(vertex_index, graph)); 
  shared_array_property_map<double,
                            property_map<graph_t, vertex_index_t>::const_type>
                            d_map(num_vertices(graph), get(vertex_index, graph)); 

  std::cout << "Running Dijkstra..." <<std::endl;
  dijkstra_shortest_paths(graph, strunk
		  ,predecessor_map(p_map)
		  .distance_map(d_map)
		  );

  std::cout << "Determining scaling factors..." <<std::endl;
  float minp=10000000, maxp=-10000000, meanp=0;
  std::ofstream of_dist("dist.dat", std::ios::out | std::ios::binary);
  std::ofstream of_paths("paths.dat", std::ios::out | std::ios::binary);
  Traits::vertex_iterator vi, vend;
  for (tie(vi, vend) = vertices(graph); vi != vend; ++vi) {
	  //std::cout << (*vi)[0]<<","<<(*vi)[1]<<","<<(*vi)[2]<<std::endl;
	  float f = d_map[*vi];
	  if (f<minp) minp = f;
	  if (f>maxp) maxp = f;
	  meanp += f;
  }

  std::cout << "Tracing paths..." <<std::endl;
  double start_threshold       = 0.3*255;
  double total_len_perc_thresh = 0.09;

  Traits::vertices_size_type strunk_idx = boost::get(vertex_index, graph, strunk);
  accumulator_set< double, features< tag::min, tag::mean, tag::max > > pathlens;
  for (tie(vi, vend) = vertices(graph); vi != vend; ++vi) {
	  vertex_descriptor v = *vi;
	  float val = A0[v[0]][v[1]][v[2]];
	  if(val<start_threshold)
		  continue;
	  float total_dist   = d_map[*vi];
	  if((total_dist-minp)/(maxp-minp) > total_len_perc_thresh)
		  continue;
	  v = *vi;
	  unsigned int cnt = 0;
	  while(1){ 
		  if(boost::get(vertex_index,graph,v) == strunk_idx)
			  break;
		  cnt++;
		  v = p_map[v];
	  }
	  pathlens(total_dist/cnt);
  }
  std::cout << "Pathlen stats: "<< min(pathlens)<<" "<<mean(pathlens)<<" "<<max(pathlens)<<std::endl;
  for (tie(vi, vend) = vertices(graph); vi != vend; ++vi) {
	  vertex_descriptor v = *vi;
	  float val = A0[v[0]][v[1]][v[2]];
	  if(val<start_threshold)
		  continue;
	  float total_dist   = d_map[*vi];
	  if((total_dist-minp)/(maxp-minp) > total_len_perc_thresh)
		  continue;
	  v = *vi;
	  unsigned int cnt = 0;
	  while(1){ 
		  if(boost::get(vertex_index,graph,v) == strunk_idx)
			  break;
		  cnt++;
		  v = p_map[v];
	  }
	  if(total_dist/cnt > max(pathlens)*0.75)
		  continue;
	  v = *vi;
	  while(1){ 
		  if(boost::get(vertex_index,graph,v) == strunk_idx)
			  break;
		  B[v[0]][v[1]][v[2]] = 255;
		  v = p_map[v];
	  }
  }
  std::cout << "Writing results..." <<std::endl;
  for (tie(vi, vend) = vertices(graph); vi != vend; ++vi) {
	  vertex_descriptor v = *vi;
	  unsigned char d = (unsigned char)(255*(d_map[v]-minp)/maxp);
	  of_dist.write((char*)&d,sizeof(unsigned char));

	  unsigned char p = B[v[0]][v[1]][v[2]];
	  of_paths.write((char*)&p,sizeof(unsigned char));
  }
  meanp /= num_vertices(graph);
  std::cout << "shortest path found: "<<minp<<std::endl;
  std::cout << "longest path found: "<<maxp<<std::endl;
  std::cout << "average path found: "<<meanp<<std::endl;
  of_dist.close();
  of_paths.close();

}
