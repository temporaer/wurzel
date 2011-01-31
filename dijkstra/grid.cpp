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
#include <boost/tuple/tuple.hpp>
#include <boost/array.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/multi_array.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include "boost/filesystem.hpp"
#include "voxelgrid.hpp"
#include "wurzel_tree.hpp"
using namespace boost::accumulators;
namespace fs = boost::filesystem;

typedef boost::multi_array_ref<unsigned char, 3> array_type;

typedef shared_array_property_map<voxel_vertex_descriptor,
						property_map<voxelgraph_t, vertex_index_t>::const_type> predecessor_map_t;
typedef shared_array_property_map<double,
						property_map<voxelgraph_t, vertex_index_t>::const_type> distance_map_t;

array_type* g_dat;

template<class T>
struct vox2arr{
	const T& A;
	vox2arr(const T& a):A(a){}
	const typename T::element& operator[](const voxel_vertex_descriptor&v)const{
		return A[v[0]][v[1]][v[2]];
	}
};

template<class T>
vox2arr<T> make_vox2arr(const T& t){ return vox2arr<T>(t); }

// map from edges to weights
voxel_edge_weight_map::reference 
voxel_edge_weight_map::operator[](key_type e) const {
	//double d = get(vertex_index, m_graph)[target(e,m_graph)];
	voxel_vertex_descriptor v = target(e,m_graph);
	array_type& a = *g_dat;
	float f = exp(-(float)a[v[0]][v[1]][v[2]]/25.6f);
	return f;
}

static const unsigned int X=256,Y=256,Z=256,XYZ=X*Y*Z;

template<class F, class T>
void write_voxelgrid(const std::string& name, voxelgraph_t& graph, const T& map){
  std::cout << "Writing results "<<name<<"..." <<std::flush;
  std::ofstream of_dist(name.c_str(), std::ios::out | std::ios::binary);
  voxel_vertex_iterator vi, vend;
  for (tie(vi, vend) = vertices(graph); vi != vend; ++vi) {
	  voxel_vertex_descriptor v = *vi;

	  F d = (F)(map[v]);
	  of_dist.write((char*)&d,sizeof(F));
  }
  of_dist.close();
  std::cout << "done." <<std::endl;
}

void find_shortest_paths(voxelgraph_t& graph, 
		voxel_vertex_descriptor& strunk,
		predecessor_map_t&p_map, distance_map_t& d_map){

  voxel_vertex_iterator vi, vend;
  bool read_p=false, read_d=false;
  if(fs::exists("data/p_map.dat")){
	  std::cout << "Reading predecessor map from file..."<<std::endl;
	  std::ifstream ifs("data/p_map.dat");
	  for (tie(vi, vend) = vertices(graph); vi != vend; ++vi)
		  ifs.read((char*)&p_map[*vi], sizeof(voxel_vertex_descriptor));
	  read_p = true;
  }
  if(fs::exists("data/d_map.dat")){
	  std::cout << "Reading distance map from file..."<<std::endl;
	  std::ifstream ifs("data/d_map.dat");
	  for (tie(vi, vend) = vertices(graph); vi != vend; ++vi)
		  ifs.read((char*)&d_map[*vi], sizeof(double));
	  read_d = true;
  }

  if(!read_p || !read_d){
	  std::cout << "Running Dijkstra..." <<std::endl;
	  dijkstra_shortest_paths(graph, strunk
			  ,predecessor_map(p_map)
			  .distance_map(d_map)
			  );
	  write_voxelgrid<double>("data/d_map.dat", graph, d_map);
	  write_voxelgrid<voxel_vertex_descriptor>("data/p_map.dat", graph, p_map);
  }
	
}

typedef accumulator_set< double, features< tag::min, tag::mean, tag::max > > stat_t;
template<class M>
stat_t voxel_stats(voxelgraph_t& graph, const M& map){
	stat_t acc;
	voxel_vertex_iterator vi, vend;
	for (tie(vi, vend) = vertices(graph); vi != vend; ++vi)
		acc(map[*vi]);
	return acc;
}

template<class T>
void paths2adjlist(voxelgraph_t& vg, wurzelgraph_t& wg, predecessor_map_t& p_map, const T& inclusion_map){
	std::cout << "Building adjacency list tree..."<<std::flush;
	typedef std::map<voxel_vertex_descriptor,wurzel_vertex_descriptor> reverse_t;
	reverse_t reverse_lookup;
	if(1){
		voxel_vertex_iterator vi, vend;
		for (tie(vi, vend) = vertices(vg); vi != vend; ++vi){
			// add all nodes which are in the inclusion map
			if(!inclusion_map[*vi]) continue;
			wurzel_vertex_descriptor v = add_vertex(*vi,wg);    
			reverse_lookup[*vi] = v;
		}
	}
	if(1){
		wurzel_vertex_iterator vi,vii, vend, vend2;
		for (tie(vi, vend) = vertices(wg); vi != vend; ++vi){
			// add edges from the predecessor of vi to vi
			voxel_vertex_descriptor v = get(vertex_name,wg)[*vi];
			voxel_vertex_descriptor w = p_map[v];
			if(w == v) continue;
			reverse_t::iterator res = reverse_lookup.find(w);
			if(res == reverse_lookup.end()) continue;
			add_edge((*res).second,*vi,wg);
		}
	}
	std::cout <<"done."<<std::endl;
}

void erode_tree(wurzelgraph_t& wg){
	std::cout <<"Eroding tree (num_nodes: "<< num_vertices(wg)<<")..."<<std::flush;
	bool modified=true;
	wurzel_vertex_iterator vi,vii, vend, vend2, next;
	wurzel_in_edge_iterator ei,eend;
	while(modified){
		std::cout <<"."<<std::flush;
		modified = false;
		tie(vi, vend) = vertices(wg);
		for (next=vi; vi != vend; vi=next){
			next++;
			if(out_degree(*vi,wg)!=0)
				continue;
			// go up from leaf to find next fork
			tie(ei,eend) = in_edges(*vi,wg);
			if(ei==eend)
				continue;
			wurzel_vertex_descriptor pred = source(*ei,wg);
			unsigned int cnt = 0;
			static const unsigned int minlen = 10;
			while(cnt++<minlen){
				if(out_degree(pred,wg)>1)
					break;
				tie(ei,eend) = in_edges(pred,wg);
				if(ei!=eend)
					pred = source(*ei,wg);
			}
			if(cnt >= minlen)
				continue;
			remove_vertex(*vi,wg);
			modified=true;
		}
	}
	std::cout <<"done (num_nodes: "<< num_vertices(wg)<<")."<<std::endl;
}


int main(int argc, char* argv[]) {
  std::ifstream dat0("../data/L2_22aug-upsampled.dat", std::ios::in | std::ios::binary);
  unsigned char* data0 = new unsigned char[XYZ];
  dat0.read((char*)data0,X*Y*Z);
  array_type A0(data0,boost::extents[X][Y][Z]);
  dat0.close();

  std::ifstream dat("../data/L2_22aug-preproc.dat", std::ios::in | std::ios::binary);
  unsigned char* data = new unsigned char[XYZ];
  dat.read((char*)data,XYZ);
  array_type A(data,boost::extents[X][Y][Z]);
  g_dat = &A;
  dat.close();

  unsigned char* paths = new unsigned char[XYZ];
  std::fill(paths, paths+XYZ, (unsigned char) 0);
  array_type B(paths,boost::extents[X][Y][Z]);

  // Define a 3x5x7 grid_graph where the second dimension doesn't wrap
  boost::array<vidx_t, 3> lengths = { { 256, 256, 256 } };
  boost::array<vidx_t, 3> strunk  = { { 109, 129,  24 } };
  voxelgraph_t graph(lengths, false); // no dim is wrapped

  voxel_vertex_descriptor first_vertex = vertex((vidx_t)0, graph);
  voxel_vertex_iterator vi, vend;

  predecessor_map_t         p_map(num_vertices(graph), get(vertex_index, graph)); 
  distance_map_t            d_map(num_vertices(graph), get(vertex_index, graph)); 

  find_shortest_paths(graph,strunk,p_map,d_map);

  std::cout << "Determining scaling factors..." <<std::endl;
  stat_t s_allpaths = voxel_stats(graph,d_map);

  std::cout << "Tracing paths..." <<std::endl;
  double start_threshold       = 0.2*255;
  double total_len_perc_thresh = 1.00;

  voxelg_traits::vertices_size_type strunk_idx = boost::get(vertex_index, graph, strunk);
  stat_t pathlens;
  vox2arr<array_type> vox2raw(A0);
  for (tie(vi, vend) = vertices(graph); vi != vend; ++vi) {
	  voxel_vertex_descriptor v = *vi;
	  float val = vox2raw[v];
	  if(val<start_threshold)
		  continue;
	  float total_dist   = d_map[*vi];
	  v = *vi;
	  unsigned int cnt = 0;
	  while(1){ 
		  if(boost::get(vertex_index,graph,v) == strunk_idx)
			  break;
		  cnt++;
		  v = p_map[v];
	  }
	  if(cnt>0)
		  pathlens(total_dist/cnt);
  }
  std::cout << "Pathlen stats: "<< min(pathlens)<<" "<<mean(pathlens)<<" "<<max(pathlens)<<std::endl;
  for (tie(vi, vend) = vertices(graph); vi != vend; ++vi) {
	  voxel_vertex_descriptor v = *vi;
	  float val = A0[v[0]][v[1]][v[2]];
	  if(val<start_threshold)
		  continue;
	  float total_dist   = d_map[*vi];
	  v = *vi;
	  unsigned int cnt = 0;
	  while(1){ 
		  if(boost::get(vertex_index,graph,v) == strunk_idx)
			  break;
		  cnt++;
		  v = p_map[v];
	  }
	  if(total_dist/cnt > max(pathlens)*0.05)
		  continue;
	  v = *vi;
	  while(1){ 
		  if(boost::get(vertex_index,graph,v) == strunk_idx)
			  break;
		  B[v[0]][v[1]][v[2]] = 255;
		  v = p_map[v];
	  }
  }
  wurzelgraph_t wgraph;
  paths2adjlist(graph,wgraph,p_map,make_vox2arr(B));
  erode_tree(wgraph);

  // fill B again with eroded graph
  std::fill(paths, paths+XYZ, (unsigned char) 0); // == B
  wurzel_vertex_iterator wi,wend;
  for (tie(wi, wend) = vertices(wgraph); wi != wend; ++wi) {
	  voxel_vertex_descriptor w = get(vertex_name,wgraph)[*wi];
	  B[w[0]][w[1]][w[2]] = 255;
  }

  write_voxelgrid<unsigned char>("data/paths.dat",graph,make_vox2arr(B));
}
