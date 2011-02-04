//=======================================================================
// Copyright 2011 University of Bonn
// Author: Hannes Schulz
//=======================================================================

#include <algorithm>
#include <iostream>
#include <vector>
#include <iterator>
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

#define SQR(X) ((X)*(X))
#define V(X) #X<<":"<<(X)<<"  "

using namespace boost::accumulators;
namespace fs = boost::filesystem;

typedef boost::multi_array_ref<float, 3> float_grid;
typedef boost::multi_array_ref<unsigned char, 3> uc_grid;

typedef shared_array_property_map<voxel_vertex_descriptor,
						property_map<voxelgraph_t, vertex_index_t>::const_type> predecessor_map_t;
typedef shared_array_property_map<double,
						property_map<voxelgraph_t, vertex_index_t>::const_type> distance_map_t;
typedef accumulator_set< double, features< tag::min, tag::mean, tag::max > > stat_t;


std::ostream& 
operator<<(std::ostream& o, const stat_t& s) {
	o << "min:"<<min(s)<<"  mean:"<<mean(s)<<"  max:"<<max(s);
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
	for(tie(oei,oeend)=out_edges(v,graph);oei!=oeend;oei++){
	    voxelg_traits::vertex_descriptor w = target(*oei, graph);
		std::cout<<"   edge"<<(++i)<<": "<<v<<"  "<<w<<" map:"<<(int)map[w]<<std::endl;
	}
}

template<class I>
double voxdist(I a, I b){
	double s = 0;
	s += SQR(*a-*b); a++; b++;
	s += SQR(*a-*b); a++; b++;
	s += SQR(*a-*b); 
	return sqrt(s);
}

float_grid* g_sato;

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
	mutable T& A;
	vox2arr(T& a):A(a){}
	typename T::element& operator[](const voxel_vertex_descriptor&v)const{
		return A[v[0]][v[1]][v[2]];
	}
};

template<class T>
vox2arr<T> make_vox2arr(T& t){ return vox2arr<T>(t); }

template<class T>
const_vox2arr<T> make_vox2arr(const T& t){ return const_vox2arr<T>(t); }

// map from edges to weights    __WEIGHT__
voxel_edge_weight_map::reference 
voxel_edge_weight_map::operator[](key_type e) const {
	//double d = get(vertex_index, m_graph)[target(e,m_graph)];
	voxel_vertex_descriptor s = source(e,m_graph);
	voxel_vertex_descriptor v = target(e,m_graph);
	float_grid& a = *g_sato;
	//double f = exp(-a[v[0]][v[1]][v[2]]);
	double f = a[v[0]][v[1]][v[2]];
	return f * voxdist(s.begin(),v.begin());
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
		predecessor_map_t&p_map, distance_map_t& d_map, bool force=false){

  voxel_vertex_iterator vi, vend;
  bool read_p=false, read_d=false;
  if(fs::exists("data/p_map.dat") && !force){
	  std::cout << "Reading predecessor map from file..."<<std::endl;
	  std::ifstream ifs("data/p_map.dat");
	  for (tie(vi, vend) = vertices(graph); vi != vend; ++vi)
		  ifs.read((char*)&p_map[*vi], sizeof(voxel_vertex_descriptor));
	  read_p = true;
  }
  if(fs::exists("data/d_map.dat") && !force){
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

template<class T, class U, class V>
void rank_op(voxelgraph_t& vg, const T& rankmap, const U& valuemap, const V& tracesmap){
	voxel_vertex_iterator vi, vend;
	voxel_edge_iterator ei, eend;
	std::ofstream ofs("data/ranks.txt");
	ofs << "0 0 0 0"<<std::endl;
	voxelg_traits::out_edge_iterator oei,oeend;
	for (tie(vi, vend) = vertices(vg); vi != vend; ++vi){
		voxel_vertex_descriptor vd = *vi;
		if(!tracesmap[*vi]) continue;
		unsigned int order = 0;
		double myval = valuemap[*vi];
		tie(oei,oeend) = out_edges(*vi,vg);
		for(;oei!=oeend;++oei)
		    if(valuemap[target(*oei,vg)] <= myval)
		        order++;
		//rankmap[*vi] = 255 - 9 * order;
		rankmap[*vi] = order;

		if(order!=0)
			continue;
		voxel_vertex_descriptor v = *vi;
		ofs << v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<255<<std::endl;
	}
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
		tie(vi, vend) = vertices(wg);
		for (; vi != vend; vi++)
			mv_map[*vi] = out_degree(*vi,wg)==0;
		tie(vi, vend) = vertices(wg);
		for (next=vi; vi != vend; vi=next){
			next++;
			if(!mv_map[*vi])
				continue;
			// go up from leaf to find next fork
			tie(ei,eend) = in_edges(*vi,wg);
			voxel_vertex_descriptor v = get(vertex_name,wg)[*vi];
			if(ei==eend)
				continue;
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
	wurzel_vertex_iterator wi,wend;
	wurzel_in_edge_iterator ei,eend;
	wurzelg_traits::adjacency_iterator      ai,aend;
	wurzelgraph_t::inv_adjacency_iterator  iai,iaend;
	property_map<wurzelgraph_t,path_orientation_t>::type po_map = get(path_orientation_t(), wg);
	for (tie(wi, wend) = vertices(wg); wi != wend; ++wi) {
		po_map[*wi] = ublas::zero_matrix<double>(3,3);
	}
	typedef std::vector<wurzel_vertex_descriptor> wd_vec_t;
	wd_vec_t neighbors, neighbors2;
	for (tie(wi, wend) = vertices(wg); wi != wend; ++wi) {
		// search forward
		tie(ai,aend) = adjacent_vertices(*wi,wg);
		std::copy(ai,aend,std::back_inserter(neighbors));
		for(wd_vec_t::iterator it=neighbors.begin();it!=neighbors.end();++it){
			tie(ai,aend) = adjacent_vertices(*it,wg);
			std::copy(ai,aend,std::back_inserter(neighbors2));
		}
		neighbors.clear();
		// search back
		tie(iai,iaend) = inv_adjacent_vertices(*wi,wg);
		std::copy(iai,iaend,std::back_inserter(neighbors));
		for(wd_vec_t::iterator it=neighbors.begin();it!=neighbors.end();++it){
			tie(iai,iaend) = inv_adjacent_vertices(*it,wg);
			std::copy(iai,iaend,std::back_inserter(neighbors2));
		}

		// now add coordinates of neighbors to covariance matrix at wi
		covmat_t& cov = po_map[*wi];
		for(wd_vec_t::iterator it=neighbors2.begin();it!=neighbors2.end();++it){
			wurzel_vertex_descriptor w = *it;
			voxel_vertex_descriptor  v = get(vertex_name,wg)[w];
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
	wurzel_edge_iterator ei,eend;
	tie(ei, eend) = edges(wg);
	std::ofstream ofs(name.c_str());
	for (; ei != eend; ei++){
		unsigned int  v = vidx_map[source(*ei,wg)];
		unsigned int  w = vidx_map[target(*ei,wg)];
		//ofs << v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<w[0]<<" "<<w[1]<<" "<<w[2]<<std::endl;
		ofs << v <<" "<< w <<std::endl;
	}
	ofs.close();
}
template<class T>
void print_wurzel_vertices(const std::string& name, wurzelgraph_t& wg, T& vidx_map){
	wurzel_vertex_iterator vi,vii, vend, vend2, next;
	wurzel_in_edge_iterator ei,eend;
	tie(vi, vend) = vertices(wg);
	std::ofstream ofs(name.c_str());
	unsigned int idx=0;
	for (; vi != vend; vi++){
		voxel_vertex_descriptor  v = get(vertex_name,wg)[*vi];
		vidx_map[*vi] = idx++;
		unsigned int deg = out_degree(*vi,wg);
		ofs << v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<deg<<std::endl;
	}
	ofs.close();
}

int main(int argc, char* argv[]) {
	bool force_recompute_dijkstra = false;
  if(argc>1)
	  force_recompute_dijkstra = true;
  std::ifstream ifs_raw("../data/L2_22aug-upsampled.dat", std::ios::in | std::ios::binary);
	if(!ifs_raw.is_open()) exit(1);
  float* raw_ptr = new float[XYZ];
  ifs_raw.read((char*)raw_ptr,XYZ*sizeof(*raw_ptr));
  float_grid Raw(raw_ptr,boost::extents[X][Y][Z]);
  ifs_raw.close();

  std::ifstream ifs_sato("../data/L2_22aug.sato", std::ios::in | std::ios::binary);
	if(!ifs_sato.is_open()) exit(1);
  //std::ifstream ifs_sato("../data/L2_22aug-preproc.dat", std::ios::in | std::ios::binary);
  float* sato_ptr = new float[XYZ];
  ifs_sato.read((char*)sato_ptr,XYZ*sizeof(*sato_ptr));
  float_grid Sato(sato_ptr,boost::extents[X][Y][Z]);
  g_sato = &Sato;
  ifs_sato.close();

  unsigned char* paths = new unsigned char[XYZ];
  std::fill(paths, paths+XYZ, (unsigned char) 0);
  uc_grid Paths(paths,boost::extents[X][Y][Z]);

  unsigned char* ranks = new unsigned char[XYZ];
  std::fill(ranks, ranks+XYZ, (unsigned char) 255);
  uc_grid Ranks(ranks,boost::extents[X][Y][Z]);

  // Define a 3x5x7 grid_graph where the second dimension doesn't wrap
  boost::array<vidx_t, 3> lengths = { { 256, 256, 256 } };
  boost::array<vidx_t, 3> strunk  = { { 109, 129,  24 } };
  voxelgraph_t graph(lengths, false); // no dim is wrapped

  voxel_vertex_descriptor first_vertex = vertex((vidx_t)0, graph);
  voxel_vertex_iterator vi, vend;

  stat_t s_raw = voxel_stats(graph, make_vox2arr(Raw));
  wurzel_vertex_iterator wi,wend;
  for (tie(vi, vend) = vertices(graph); vi != vend; ++vi) {
      voxel_vertex_descriptor v = *vi;
      float& f = Sato[v[0]][v[1]][v[2]];
			f  = exp(-10.f * log(1.f+f)/log(2.f) );
			//f  = 1.f-log(1.f+f)/log(2.f);

      float& g = Raw[v[0]][v[1]][v[2]];
			g  = (g-min(s_raw))/(max(s_raw)-min(s_raw));
  }
  stat_t s_sato = voxel_stats(graph, make_vox2arr(Sato));
  std::cout << "Sato: " << s_sato <<std::endl;

  predecessor_map_t         p_map(num_vertices(graph), get(vertex_index, graph)); 
  distance_map_t            d_map(num_vertices(graph), get(vertex_index, graph)); 

  find_shortest_paths(graph,strunk,p_map,d_map,force_recompute_dijkstra);
  std::cout << "Dmap: " << voxel_stats(graph,d_map) <<std::endl;

  std::cout << "Determining scaling factors..." <<std::endl;
  stat_t s_allpaths = voxel_stats(graph,d_map);

  std::cout << "Tracing paths..." <<std::endl;
  double start_threshold       = 0.2;
  double total_len_perc_thresh = 0.50;
  double avg_len_perc_thresh   = 0.15;

  voxelg_traits::vertices_size_type strunk_idx = boost::get(vertex_index, graph, strunk);
  stat_t s_avg_pathlen, s_pathlen, s_cnt;
  vox2arr<float_grid> vox2raw(Raw);
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
	  s_cnt(cnt);
	  s_pathlen(total_dist);
	  if(cnt>0)
		  s_avg_pathlen(total_dist/cnt);
  }
  std::cout << "Total   Pathlens: "<< s_pathlen<<std::endl;
  std::cout << "Average Pathlens: "<< s_avg_pathlen<<std::endl;
  std::cout << "Hop     Pathlens: "<< s_cnt<<std::endl;
  for (tie(vi, vend) = vertices(graph); vi != vend; ++vi) {
	  voxel_vertex_descriptor v = *vi;
	  float val = Raw[v[0]][v[1]][v[2]];
	  if(val<start_threshold)
		  continue;
	  float total_dist   = d_map[*vi];
	  if(((total_dist-min(s_pathlen))/(max(s_pathlen)-min(s_pathlen))) > total_len_perc_thresh)
		  continue;
	  v = *vi;
	  unsigned int cnt = 0;
	  while(1){ 
		  if(boost::get(vertex_index,graph,v) == strunk_idx)
			  break;
		  cnt++;
		  v = p_map[v];
	  }
	  if(total_dist/cnt > max(s_avg_pathlen)*avg_len_perc_thresh)
		  continue;
	  v = *vi;
	  while(1){ 
		  Paths[v[0]][v[1]][v[2]] = 255.f * (cnt-- - min(s_cnt))/(max(s_cnt)-min(s_cnt));
		  if(boost::get(vertex_index,graph,v) == strunk_idx)
			  break;
		  v = p_map[v];
	  }
  }
  
  // find local ranks
  rank_op(graph,make_vox2arr(Ranks),make_vox2arr(Sato),make_vox2arr(Paths));

  wurzelgraph_t wgraph;
  paths2adjlist(graph,wgraph,p_map,make_vox2arr(Paths));
  remove_nonmax_nodes(wgraph,make_vox2arr(Ranks));
  //erode_tree(wgraph);
  //merge_deg2_nodes(wgraph);

  std::map<wurzel_vertex_descriptor,unsigned int> idx_map;
  print_wurzel_vertices("data/vertices.txt",wgraph,idx_map);
  print_wurzel_edges("data/edges.txt",wgraph,idx_map);
  write_voxelgrid<unsigned char>("data/paths.dat",graph,make_vox2arr(Paths));
  write_voxelgrid<unsigned char>("data/ranks.dat",graph,make_vox2arr(Ranks));
}
