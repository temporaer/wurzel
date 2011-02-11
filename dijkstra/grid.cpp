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
#include <vector>
#include <iterator>
#include <fstream>
#include <map>
#include <boost/tuple/tuple.hpp>
#include <boost/array.hpp>
#include <boost/range/irange.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/multi_array.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/foreach.hpp>
#include "boost/filesystem.hpp"
#include "voxelgrid.hpp"
#include "wurzel_tree.hpp"
#define SQR(X) ((X)*(X))
#define V(X) #X<<":"<<(X)<<"  "
#define foreach BOOST_FOREACH

using namespace boost::accumulators;
namespace fs = boost::filesystem;

typedef boost::multi_array_ref<float, 3> float_grid;
typedef boost::multi_array_ref<unsigned char, 3> uc_grid;

typedef shared_array_property_map<voxel_vertex_descriptor,
						property_map<voxelgraph_t, vertex_index_t>::const_type> predecessor_map_t;
typedef shared_array_property_map<double,
						property_map<voxelgraph_t, vertex_index_t>::const_type> distance_map_t;
typedef accumulator_set< double, features< tag::min, tag::mean, tag::max > > stat_t;

std::string
getfn(const std::string& base, const std::string& marker, const std::string& ext){
	std::stringstream ss;
	if(marker.length()>0){
		ss << base << "-"<<marker << "."<<ext;
	}else{
		ss << base << "."<<ext;
	}
	return ss.str();
}


std::ostream& 
operator<<(std::ostream& o, const stat_t& s) {
	o<<"min: "<<min(s)<<"  mean: "<<mean(s)<<"  max:"<<max(s);
  return o;
}
std::ostream& 
operator<<(std::ostream& o, const voxelg_traits::vertex_descriptor& vertex_to_print) {
  o << "(" << vertex_to_print[0] << ", " << vertex_to_print[1] <<
    ", " << vertex_to_print[2] << ")";
  return o;
}
template<class G, class M>
void print_neighbors(G& graph, voxelg_traits::vertex_descriptor v, const M& map){
	std::cout<< "start vertex: "<<v<<std::endl;
	std::cout << "   out-degree v: "<<out_degree(v,graph)<<std::endl;

	voxelg_traits::out_edge_iterator oei, oeend;
	unsigned int i=0;
	foreach(voxel_edge_descriptor& e, out_edges(v,graph)){
	    voxelg_traits::vertex_descriptor w = target(e, graph);
		std::cout<<"   edge"<<(++i)<<": "<<v<<"  "<<w<<" map:"<<(int)map[w]<<std::endl;
	}
}

template<class I, class J>
inline double voxdist(I a, J b){
	double s = 0;
	s += SQR(*a-*b); a++; b++;
	s += SQR(*a-*b); a++; b++;
	s += SQR(*a-*b); 
	return sqrt(s);
}


template<class T>
struct const_vox2arr{
	const T& A;
	const_vox2arr(const T& a):A(a){}
	const typename T::element& operator[](const voxel_vertex_descriptor&v)const{
		return A[v[0]][v[1]][v[2]];
	}
};

template<class T>
struct vox2arr{
	T& A;
	vox2arr(T& a):A(a){}
	typename T::element& operator[](const voxel_vertex_descriptor&v)const{
		return A[v[0]][v[1]][v[2]];
	}
};

template<class T>
vox2arr<T> make_vox2arr(T& t){ return vox2arr<T>(t); }

template<class T>
const_vox2arr<T> make_vox2arr(const T& t){ return const_vox2arr<T>(t); }

vox2arr<float_grid>* g_sato, *g_ev10, *g_ev11, *g_ev12;
boost::array<vidx_t, 3> g_strunk;

// map from edges to weights
voxel_edge_weight_map::reference 
voxel_edge_weight_map::operator[](key_type e) const {
	voxel_vertex_descriptor s = source(e,m_graph);
	voxel_vertex_descriptor t = target(e,m_graph);
	vox2arr<float_grid>& sato = *g_sato;
	vox2arr<float_grid>& ev10 = *g_ev10;
	vox2arr<float_grid>& ev11 = *g_ev11;
	vox2arr<float_grid>& ev12 = *g_ev12;

	double d = voxdist(s.begin(),t.begin());

	//double gd = ev10[s]*ev10[t] + ev11[s]*ev11[t] + ev12[s]*ev12[t];
	//double  g = (1.0 + 30.0*exp(-10.0*gd*gd));

	double v = exp( - 50.0 * (sato[t] + sato[s]) );

	//double hs = voxdist(s.begin(),g_strunk.begin());
	//double ht = voxdist(t.begin(),g_strunk.begin()); 
	//double h  = 1.0 + exp((hs-ht)/d); // target should be further away from source

	//return g * v * d * h;
	return v * d;
}

template<class F, class T>
void write_voxelgrid(const std::string& name, voxelgraph_t& graph, const T& map){
  std::cout << "Writing results "<<name<<"..." <<std::flush;
  std::ofstream of_dist(name.c_str(), std::ios::out | std::ios::binary);
  foreach (const voxel_vertex_descriptor& v,vertices(graph)) {
	  F d = (F)(map[v]);
	  of_dist.write((char*)&d,sizeof(F));
  }
  of_dist.close();
  std::cout << "done." <<std::endl;
}

void find_shortest_paths(const std::string& base, 
		voxelgraph_t& graph, 
		voxel_vertex_descriptor& strunk,
		predecessor_map_t&p_map, distance_map_t& d_map, bool force=false){

  bool read_p=false, read_d=false;
  if(fs::exists(getfn(base,"p_map","dat")) && !force){
	  std::cout << "Reading predecessor map from file..."<<std::endl;
	  std::ifstream ifs("data/p_map.dat");
	  foreach (const voxel_vertex_descriptor& v, vertices(graph))
		  ifs.read((char*)&p_map[v], sizeof(voxel_vertex_descriptor));
	  read_p = true;
  }
  if(fs::exists(getfn(base,"d_map","dat")) && !force){
	  std::cout << "Reading distance map from file..."<<std::endl;
	  std::ifstream ifs("data/d_map.dat");
	  foreach (const voxel_vertex_descriptor& v, vertices(graph))
		  ifs.read((char*)&d_map[v], sizeof(double));
	  read_d = true;
  }

  if(!read_p || !read_d){
	  std::cout << "Running Dijkstra..." <<std::endl;

	  voxel_edge_weight_map wt(graph);

	  dijkstra_shortest_paths(graph, strunk
			 ,predecessor_map(p_map)
			 .distance_map(d_map)
			 );
	  write_voxelgrid<double>(getfn(base,"d_map","dat"), graph, d_map);
	  write_voxelgrid<voxel_vertex_descriptor>(getfn(base,"p_map","dat"), graph, p_map);
  }
	
}

template<class M>
stat_t voxel_stats(voxelgraph_t& graph, const M& map){
	stat_t acc;
	voxel_vertex_iterator vi, vend;
	foreach ( const voxel_vertex_descriptor& v, vertices(graph))
		acc(map[v]);
	return acc;
}

template<class T>
void paths2adjlist(voxelgraph_t& vg, wurzelgraph_t& wg, predecessor_map_t& p_map, const T& inclusion_map){
	std::cout << "Building adjacency list tree..."<<std::flush;
	typedef std::map<voxel_vertex_descriptor,wurzel_vertex_descriptor> reverse_t;
	reverse_t reverse_lookup;
	if(1){
		foreach (const voxel_vertex_descriptor& v0, vertices(vg)){
			// add all nodes which are in the inclusion map
			if(!inclusion_map[v0]) continue;
			wurzel_vertex_descriptor v = add_vertex(v0,wg);    
			reverse_lookup[v0] = v;
		}
	}
	if(1){
		foreach(wurzel_vertex_descriptor& wv, vertices(wg)){
			// add edges from the predecessor of vi to vi
			voxel_vertex_descriptor v = get(vertex_name,wg)[wv];
			voxel_vertex_descriptor w = p_map[v];
			if(w == v) continue;
			reverse_t::iterator res = reverse_lookup.find(w);
			if(res == reverse_lookup.end()) continue;
			add_edge((*res).second,wv,wg);
		}
	}
	std::cout <<"done ("<<V(num_vertices(wg))<<")"<<std::endl;
}

template<class T, class U, class V>
void rank_op(const std::string& base, voxelgraph_t& vg, T rankmap, const U& valuemap, const V& tracesmap){
	std::ofstream ofs(getfn(base,"ranks","txt").c_str());
	ofs << "0 0 0 0"<<std::endl;
	voxelg_traits::out_edge_iterator oei,oeend;
	unsigned int cnt = 0, cnt_all=0;
	foreach (const voxel_vertex_descriptor& vd, vertices(vg)){
		if(!tracesmap[vd]) continue;
		unsigned int order = 0;
		double myval = valuemap[vd];
		foreach(const voxel_edge_descriptor& e, out_edges(vd,vg)){
		    if(valuemap[target(e,vg)] >= myval)
		        order++;
		}
		//rankmap[*vi] = 255 - 9 * order;
		rankmap[vd] = order;
		if(order==0) cnt++;

		if(order!=0){
			cnt_all++;
			continue;
		}
		ofs << vd[0]<<" "<<vd[1]<<" "<<vd[2]<<" "<<255<<std::endl;
		cnt_all++;
	}
	std::cout <<"Determined "<< cnt<<" of "<<cnt_all<<" to be maxima"<<std::endl;
}
void erode_tree(wurzelgraph_t& wg){
	std::cout <<"Eroding tree (num_nodes: "<< num_vertices(wg)<<")..."<<std::flush;
	bool modified=true;
	wurzel_vertex_iterator vi,vii, vend, vend2, next;
	wurzel_in_edge_iterator ei,eend;
	property_map<wurzelgraph_t,marked_vertex_t>::type mv_map = get(marked_vertex, wg);
	while(modified){
		std::cout <<"."<<std::flush;
		modified = false;
		foreach(wurzel_vertex_descriptor& v, vertices(wg)){
			mv_map[v] = out_degree(v,wg)==0;
		}
		tie(vi, vend) = vertices(wg);
		for (next=vi; vi != vend; vi=next){
			next++;
			if(!mv_map[*vi])
				continue;
			// go up from leaf to find next fork
			tie(ei,eend) = in_edges(*vi,wg);
			if(ei==eend)
				continue;
			voxel_vertex_descriptor     v = get(vertex_name,wg)[*vi];
			wurzel_vertex_descriptor pred = source(*ei,wg);
			unsigned int cnt = 0;
			static const unsigned int minlen = 5;
			while(cnt++<minlen){
				if(out_degree(pred,wg)>1)
					break;
				tie(ei,eend) = in_edges(pred,wg);
				if(ei!=eend)
					pred = source(*ei,wg);
			}
			if(cnt >= minlen)
				continue;
			clear_vertex(*vi,wg);
			remove_vertex(*vi,wg);
			modified=true;
		}
	}
	std::cout <<"done (num_nodes: "<< num_vertices(wg)<<")."<<std::endl;
}

template<class T>
void remove_nonmax_nodes(wurzelgraph_t& wg, const T& maxmap){
	std::cout <<"Removing nonmaximum nodes (num_nodes: "<< num_vertices(wg)<<")..."<<std::flush;
	wurzel_vertex_iterator wi,wend,next;
	wurzel_in_edge_iterator ei,eend;
	wurzelg_traits::adjacency_iterator      ai,aend;
	wurzelgraph_t::inv_adjacency_iterator  iai,iaend;
	property_map<wurzelgraph_t,path_orientation_t>::type po_map = get(path_orientation, wg);
	bool changed = true;
	while(changed){
		changed = false;
		tie(wi, wend) = vertices(wg);
		for (next=wi; wi != wend;wi=next) {
			next++;
			if(out_degree(*wi,wg) > 1)
				continue;
			if(in_degree(*wi,wg) < 1)
				continue;
			if(maxmap[get(vertex_name,wg)[*wi]]==0)
				continue;
			tie(ai,aend)   = adjacent_vertices(*wi,wg);
			tie(iai,iaend) = inv_adjacent_vertices(*wi,wg);
			clear_vertex(*wi,wg);
			remove_vertex(*wi,wg);
			if(iai!=iaend && ai!=aend)
				add_edge(*iai,*ai,wg);
			changed = true;
		}
	}
	std::cout <<"done (num_nodes: "<< num_vertices(wg)<<")."<<std::endl;
}
void merge_deg2_nodes(wurzelgraph_t& wg){
	std::cout <<"Removing deg2nodes (num_nodes: "<< num_vertices(wg)<<")..."<<std::flush;
	wurzel_vertex_iterator wi,wend,next;
	wurzel_in_edge_iterator ei,eend;
	wurzelg_traits::adjacency_iterator      ai,aend;
	wurzelgraph_t::inv_adjacency_iterator  iai,iaend;
	property_map<wurzelgraph_t,path_orientation_t>::type po_map = get(path_orientation, wg);
	bool changed = true;
	while(changed){
		changed = false;
		tie(wi, wend) = vertices(wg);
		for (next=wi; wi != wend;wi=next) {
			next++;
			if(out_degree(*wi,wg) != 1)
				continue;
			if(in_degree(*wi,wg) != 1)
				continue;
			tie(ai,aend)   = adjacent_vertices(*wi,wg);
			tie(iai,iaend) = inv_adjacent_vertices(*wi,wg);
			clear_vertex(*wi,wg);
			remove_vertex(*wi,wg);
			add_edge(*iai,*ai,wg);
			changed = true;
		}
	}
	std::cout <<"done (num_nodes: "<< num_vertices(wg)<<")."<<std::endl;
}
void merge_nodes(wurzelgraph_t& wg){
	std::cout << "Merging nodes..."<<std::flush;
	// determine local orientations
	wurzelg_traits::adjacency_iterator      ai,aend;
	wurzelgraph_t::inv_adjacency_iterator  iai,iaend;
	property_map<wurzelgraph_t,path_orientation_t>::type po_map = get(path_orientation_t(), wg);
	foreach (wurzel_vertex_descriptor& v, vertices(wg)) {
		po_map[v] = ublas::zero_matrix<double>(3,3);
	}
	typedef std::vector<wurzel_vertex_descriptor> wd_vec_t;
	wd_vec_t neighbors, neighbors2;
	foreach (const wurzel_vertex_descriptor& wv, vertices(wg)) {
		// search forward
		tie(ai,aend) = adjacent_vertices(wv,wg);
		std::copy(ai,aend,std::back_inserter(neighbors));
		foreach(const wurzel_vertex_descriptor& v2, neighbors){
			tie(ai,aend) = adjacent_vertices(v2,wg);
			std::copy(ai,aend,std::back_inserter(neighbors2));
		}
		neighbors.clear();
		// search back
		tie(iai,iaend) = inv_adjacent_vertices(wv,wg);
		std::copy(iai,iaend,std::back_inserter(neighbors));
		foreach(const wurzel_vertex_descriptor& v2, neighbors){
			tie(iai,iaend) = inv_adjacent_vertices(v2,wg);
			std::copy(iai,iaend,std::back_inserter(neighbors2));
		}

		// now add coordinates of neighbors to covariance matrix at wi
		covmat_t& cov = po_map[wv];
		foreach(const wurzel_vertex_descriptor& w, neighbors2){
			const voxel_vertex_descriptor&  v = get(vertex_name,wg)[w];
			cov(0,0) += v[0]*v[0];
			cov(0,1) += v[0]*v[1];
			cov(0,2) += v[0]*v[2];
			cov(1,1) += v[1]*v[1];
			cov(1,2) += v[1]*v[2];
			cov(2,2) += v[2]*v[2];
		}
		cov(1,0) = cov(0,1);
		cov(2,0) = cov(0,2);
		cov(2,1) = cov(1,2);
	}
	std::cout << "done."<<std::endl;
}

template<class T>
void print_wurzel_edges(const std::string& name, wurzelgraph_t& wg, T& vidx_map){
	std::ofstream ofs(name.c_str());
	foreach (const wurzel_edge_descriptor& e, edges(wg)){
		unsigned int  v = vidx_map[source(e,wg)];
		unsigned int  w = vidx_map[target(e,wg)];
		//ofs << v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<w[0]<<" "<<w[1]<<" "<<w[2]<<std::endl;
		ofs << v <<" "<< w <<std::endl;
	}
	ofs.close();
}
template<class T>
void print_wurzel_vertices(const std::string& name, wurzelgraph_t& wg, T& vidx_map){
	std::ofstream ofs(name.c_str());
	unsigned int idx=0;
	foreach (wurzel_vertex_descriptor& wd, vertices(wg)){
		voxel_vertex_descriptor  v = get(vertex_name,wg)[wd];
		vidx_map[wd] = idx++;
		unsigned int deg = out_degree(wd,wg);
		ofs << v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<deg<<" "<<(*g_ev10)[v]<<" "<<(*g_ev11)[v]<<" "<<(*g_ev12)[v]<<std::endl;
	}
	ofs.close();
}

template<class T>
boost::multi_array_ref<T, 3> 
read3darray(std::string fn, unsigned int X, unsigned int Y, unsigned int Z){
  std::cout << "Reading `"<<fn<<"', bytes="<<sizeof(T)<<std::endl;
  std::ifstream dat(fn.c_str(), std::ios::in | std::ios::binary);
  if(!dat.is_open())	{
	  std::cerr << "Could not open "<<fn<<" --> exiting."<<std::endl;
	  exit(1);
  }
  T* data = new T[X*Y*Z];
  dat.read((char*)data,X*Y*Z*sizeof(T));
  if(dat.fail()){
	  std::cerr << "Error while reading "<<fn<<" --> exiting."<<std::endl;
	  exit(1);
  }
  dat.close();
  return boost::multi_array_ref<T,3>(data,boost::extents[X][Y][Z]);
}

int main(int argc, char* argv[]) {
	static const unsigned int X=256,Y=256,Z=256,XYZ=X*Y*Z;
	bool force_recompute_dijkstra = false;
	if(argc<2) {
			std::cerr << "Usage: " << argv[0]<<" {basename} [force]"<<std::endl;
			exit(1);
		}
	const std::string base=argv[1];
	if(argc>2)
		force_recompute_dijkstra = true;
	
  // Define a 3x5x7 grid_graph where the second dimension doesn't wrap
  boost::array<vidx_t, 3> lengths = { { 256, 256, 256 } };
  boost::array<vidx_t, 3> strunk  = { { 109, 129,  24 } };
  g_strunk = strunk;
  voxelgraph_t graph(lengths, false); // no dim is wrapped

  float_grid Raw  = read3darray<float>(getfn(base,"upsampled","dat"),X,Y,Z);
  float_grid Sato = read3darray<float>(getfn(base,"","sato"),X,Y,Z); g_sato = new vox2arr<float_grid>(Sato);
  float_grid ev10 = read3darray<float>(getfn(base,"","ev10"),X,Y,Z); g_ev10 = new vox2arr<float_grid>(ev10);
  float_grid ev11 = read3darray<float>(getfn(base,"","ev11"),X,Y,Z); g_ev11 = new vox2arr<float_grid>(ev11);
  float_grid ev12 = read3darray<float>(getfn(base,"","ev12"),X,Y,Z); g_ev12 = new vox2arr<float_grid>(ev12);

  std::cout << "Sato stats in file: " << voxel_stats(graph,make_vox2arr(Sato)) <<std::endl;

  unsigned char* paths = new unsigned char[XYZ];
  std::fill(paths, paths+XYZ, (unsigned char) 0);
  uc_grid Paths(paths,boost::extents[X][Y][Z]);

  unsigned char* ranks = new unsigned char[XYZ];
  std::fill(ranks, ranks+XYZ, (unsigned char) 255);
  uc_grid Ranks(ranks,boost::extents[X][Y][Z]);

  float* flow = new float[XYZ];
  std::fill(flow, flow+XYZ, (float) 0);
  float_grid Flow(flow,boost::extents[X][Y][Z]);


  voxel_vertex_descriptor first_vertex = vertex((vidx_t)0, graph);
  voxel_vertex_iterator vi, vend;

  wurzel_vertex_iterator wi,wend;
  stat_t s_raw  = voxel_stats(graph, make_vox2arr(Raw));
  foreach(const voxel_vertex_descriptor& v, vertices(graph)) {
	  float& f = Sato[v[0]][v[1]][v[2]];
			//f = std::max(0.f, f-0.2f) * 1.f/0.7f;
			//f  = exp(-10.f * log(1.f+f)/log(2.f) );
			//f  = exp(-15.f * f );

	  float& g = Raw[v[0]][v[1]][v[2]];
	  f *= 1.0f+2.0f*g;
	  //      g  = (g-min(s_raw))/(max(s_raw)-min(s_raw));
  }
  stat_t s_sato = voxel_stats(graph, make_vox2arr(Sato));
  std::cout << "Raw:  " << s_raw <<std::endl;
  std::cout << "Sato: " << s_sato <<std::endl;

  predecessor_map_t         p_map(num_vertices(graph), get(vertex_index, graph)); 
  distance_map_t            d_map(num_vertices(graph), get(vertex_index, graph)); 

  find_shortest_paths(base,graph,strunk,p_map,d_map,force_recompute_dijkstra);
  std::cout << "Dmap: " << voxel_stats(graph,d_map) <<std::endl;

  std::cout << "Determining scaling factors..." <<std::endl;
  stat_t s_allpaths = voxel_stats(graph,d_map);

  std::cout << "Tracing paths..." <<std::endl;
  double start_threshold       = 0.1;
  double total_len_perc_thresh = 1.00;
  double avg_len_perc_thresh   = 0.10;
  double min_flow_thresh       = 0.000500;

  voxelg_traits::vertices_size_type strunk_idx = boost::get(vertex_index, graph, strunk);
  stat_t s_avg_pathlen, s_pathlen, s_cnt, s_flow;
  vox2arr<float_grid> vox2raw(Raw);
  foreach (const voxel_vertex_descriptor& v, vertices(graph)) {
	  // determine total path length statistic
	  if(vox2raw[v] < start_threshold)
		  continue;
	  s_pathlen(d_map[v]);
  }

  foreach (const voxel_vertex_descriptor& v, vertices(graph)) {
	  // determine avg path costs statistic
	  if(vox2raw[v] < start_threshold)
		  continue;
	  if(((d_map[v]-min(s_pathlen))/(max(s_pathlen)-min(s_pathlen))) > total_len_perc_thresh)
		  continue;
	  double vox_dist = 0.0;
		voxel_vertex_descriptor v2 = v;
	  while(1){ 
		  Flow[v2[0]][v2[1]][v2[2]] += vox2raw[v];
		  if(boost::get(vertex_index,graph,v2) == strunk_idx)
			  break;
		  vox_dist += voxdist(p_map[v2].begin(), v2.begin());
		  v2 = p_map[v2];
	  }
	  s_cnt(vox_dist);
	  if(vox_dist>0)
		  s_avg_pathlen(d_map[v]/vox_dist);
  }

  foreach (const voxel_vertex_descriptor& v, vertices(graph)) {
	  // determine flow statistic
	  if(vox2raw[v] < start_threshold)
		  continue; // too weak at start
	  if(((d_map[v]-min(s_pathlen))/(max(s_pathlen)-min(s_pathlen))) > total_len_perc_thresh)
		  continue; // too long
	  s_flow(Flow[v[0]][v[1]][v[2]]);
  }

  std::cout << "Total   Pathlens: "<< s_pathlen<<std::endl;
  std::cout << "Average Pathlens: "<< s_avg_pathlen<<std::endl;
  std::cout << "Hop     Pathlens: "<< s_cnt<<std::endl;

  foreach (const voxel_vertex_descriptor& v0, vertices(graph)) {
	  if(vox2raw[v0]<start_threshold) 
		  continue;                  // weak signal at start point
	  float total_dist   = d_map[v0];
	  if(((total_dist-min(s_pathlen))/(max(s_pathlen)-min(s_pathlen))) > total_len_perc_thresh)
		  continue;                  // too long
	  float flow = Flow[v0[0]][v0[1]][v0[2]];
	  if((flow-min(s_flow))/(max(s_flow)-min(s_flow)) < min_flow_thresh)
		  continue;                  // not enough mass here
	  voxel_vertex_descriptor v = v0;
	  double vox_dist = 0.0;
	  while(1){ 
		  if(boost::get(vertex_index,graph,v) == strunk_idx)
			  break;
		  vox_dist += voxdist(p_map[v].begin(), v.begin());
		  v = p_map[v];
	  }
	  if(total_dist/vox_dist > max(s_avg_pathlen)*avg_len_perc_thresh)
		  continue;
	  v = *vi;
	  while(1){ 
		  Paths[v[0]][v[1]][v[2]] = 255;
		  if(boost::get(vertex_index,graph,v) == strunk_idx)
			  break;
		  v = p_map[v];
	  }
  }
  
  // find local ranks
  rank_op(base,graph,make_vox2arr(Ranks),make_vox2arr(Sato),make_vox2arr(Paths));

  wurzelgraph_t wgraph;
  paths2adjlist(graph,wgraph,p_map,make_vox2arr(Paths));
  remove_nonmax_nodes(wgraph,make_vox2arr(Ranks));
  //erode_tree(wgraph);
  //merge_deg2_nodes(wgraph);

  std::map<wurzel_vertex_descriptor,unsigned int> idx_map;
  print_wurzel_vertices(getfn(base,"vertices","txt"),wgraph,idx_map);
  print_wurzel_edges(   getfn(base,"edges","txt"),wgraph,idx_map);
  write_voxelgrid<unsigned char>(getfn(base,"paths","dat"),graph,make_vox2arr(Paths));
  write_voxelgrid<unsigned char>(getfn(base,"ranks","dat"),graph,make_vox2arr(Ranks));
}
