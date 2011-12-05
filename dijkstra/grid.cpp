//=======================================================================
// Copyright 2011 University Bonn
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
//#include <boost/range/irange.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/multi_array.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/array.hpp>
#include "boost/filesystem.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <tbb/parallel_for_each.h>
#include <tbb/atomic.h>
#include "inc/voxelgrid.hpp"
#include "inc/wurzel_tree.hpp"
#include "inc/voxel_accessors.hpp"
#include "inc/voxel_normals.hpp"
#include "inc/wurzel_info.hpp"
#include "inc/grid_config.hpp"
//#include "gaussian_fit.hpp"
#include "inc/profiler.hpp"
#include "inc/rprop.hpp"

#define V(X) #X<<":"<<(X)<<"  "
#define foreach BOOST_FOREACH

using namespace boost::accumulators;
namespace fs = boost::filesystem;
using boost::optional;

typedef boost::multi_array_ref<float, 3> float_grid;
typedef boost::multi_array_ref<unsigned char, 3> uc_grid;

typedef shared_array_property_map<voxel_vertex_descriptor,
						property_map<voxelgraph_t, vertex_index_t>::const_type> predecessor_map_t;
typedef shared_array_property_map<double,
						property_map<voxelgraph_t, vertex_index_t>::const_type> distance_map_t;
typedef accumulator_set< double, features< tag::min, tag::mean, tag::max > > stat_t;
/*
 * Dijkstra Early Stopping
 */

class dijkstra_finish : public std::exception { };

template<class DMap>
class maxdist_visitor : public default_dijkstra_visitor {
public:
    maxdist_visitor(const DMap &d_, double m) : d(d_),maxdist(m) {}
    template <typename Vertex, typename Graph>
        void finish_vertex(const Vertex& u, const Graph& g) {
            if (d[u]> maxdist) throw dijkstra_finish();
        }
private:
    const DMap &d;
    const double maxdist;
};



/**
 * Return a file name which is built from the parts in the parameters
 * @param base    the name without extention
 * @param marker  this will be attached to base with a "-"
 * @param ext     this will be the extention
 * @return a file name which is built from the parts in the parameters
 */
std::string
getfn(const std::string& base, const std::string& marker, const std::string& ext){
	std::stringstream ss;
	if(marker.length()>0){
		ss << base << "/"<<marker << "."<<ext;
	}else{
		ss << base << "/"<<ext;
	}
	return ss.str();
}

/**
 * Print statistics to a stream for convenience
 * @param o  the stream
 * @param s  statistic
 * @return the stream
 */
std::ostream& 
operator<<(std::ostream& o, const stat_t& s) {
	o<<"min: "<<min(s)<<"  mean: "<<mean(s)<<"  max:"<<max(s);
  return o;
}
/**
 * print a vertex to a stream for convenience
 * @param o   the stream
 * @param vertex_to_print  the vertex
 * @return the stream
 */
std::ostream& 
operator<<(std::ostream& o, const voxelg_traits::vertex_descriptor& vertex_to_print) {
  o << "(" << vertex_to_print[0] << ", " << vertex_to_print[1] <<
    ", " << vertex_to_print[2] << ")";
  return o;
}

/**
 * For a given vertex, print all its neighbors to std::cout,
 * with additional integer info taken from map.
 */
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


vox2arr<float_grid>* g_sato, *g_ev10, *g_ev11, *g_ev12;
boost::array<vidx_t, 3> g_strunk;

/**
 * @return   the weight of an edge in the graph
 * @param e  the edge
 */
voxel_edge_weight_map::reference 
voxel_edge_weight_map::operator[](key_type e) const {
	voxel_vertex_descriptor s = source(e,m_graph);
	voxel_vertex_descriptor t = target(e,m_graph);
	vox2arr<float_grid>& sato = *g_sato;
	//vox2arr<float_grid>& ev10 = *g_ev10;
	//vox2arr<float_grid>& ev11 = *g_ev11;
	//vox2arr<float_grid>& ev12 = *g_ev12;

	double d = voxdist(s.begin(),t.begin());

	//double gd = ev10[s]*ev10[t] + ev11[s]*ev11[t] + ev12[s]*ev12[t];
	//double  g = (1.0 + 30.0*exp(-10.0*gd*gd));

	double v = exp( - 100.0 * (sato[t] + sato[s]) );

	//double hs = voxdist(s.begin(),g_strunk.begin());
	//double ht = voxdist(t.begin(),g_strunk.begin()); 
	//double h  = 1.0 + exp((hs-ht)/d); // target should be further away from source

	//return g * v * d * h;
	return v * d;
}

/**
 * write the data of a whole voxelgrid to a raw file
 * @param name   filename
 * @param graph  determines connectivity and order
 * @param map    the information written to the file
 */
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

/**
 * find the shortest path through a voxel graph
 * @param base   basename of the file, used to load/store cached p_map/d_map
 * @param graph  determines connectivtity
 * @param strunk the position from which we start the search
 * @param p_map  predecessor-map
 * @param d_map  distance to strunk map
 * @param force  force recomputing, even when cached d_map/p_map found
 */
void find_shortest_paths(const std::string& base, 
		voxelgraph_t& graph, 
		voxel_vertex_descriptor& strunk,
		predecessor_map_t&p_map, distance_map_t& d_map, bool force=false){

  bool read_p=false, read_d=false;
  if(fs::exists(getfn(base,"p_map","dat")) && !force){
	  std::cout << "Reading predecessor map from file..."<<std::endl;
	  std::ifstream ifs(getfn(base,"p_map","dat").c_str());
	  foreach (const voxel_vertex_descriptor& v, vertices(graph))
		  ifs.read((char*)&p_map[v], sizeof(voxel_vertex_descriptor));
	  read_p = true;
  }
  if(fs::exists(getfn(base,"d_map","dat")) && !force){
	  std::cout << "Reading distance map from file..."<<std::endl;
	  std::ifstream ifs(getfn(base,"d_map", "dat").c_str());
	  foreach (const voxel_vertex_descriptor& v, vertices(graph))
		  ifs.read((char*)&d_map[v], sizeof(double));
	  read_d = true;
  }

  if(!read_p || !read_d){
	  std::cout << "Running Dijkstra..." <<std::endl;

	  voxel_edge_weight_map wt(graph);

	  maxdist_visitor<distance_map_t> vis(d_map, 20.0);
	  try{
		  dijkstra_shortest_paths(graph, strunk
				  ,predecessor_map(p_map)
				  .distance_map(d_map)
				  .visitor(vis)
				  );
	  } catch(const dijkstra_finish&) {}
	  write_voxelgrid<double>(getfn(base,"d_map","dat"), graph, d_map);
	  write_voxelgrid<voxel_vertex_descriptor>(getfn(base,"p_map","dat"), graph, p_map);
  }
	
}

/**
 * calculate the statistics of a grid
 * @param map    where the values are stored
 * @return       statistics-object
 */
template<class Stats>
stat_t voxel_stats(const float_grid& map, Stats& acc ){
	for(int i=0;i<map.shape()[0];i++)
		for(int j=0;j<map.shape()[1];j++)
			for(int k=0;k<map.shape()[2];k++)
				acc(map[i][j][k]);
	return acc;
}
/**
 * calculate the statistics of a grid
 * @param map    where the values are stored
 * @return       statistics-object
 */
template<class M>
stat_t voxel_stats(const M& map){
	stat_t acc;
	for(unsigned int i=0;i<map.shape()[0];i++)
		for(unsigned int j=0;j<map.shape()[1];j++)
			for(unsigned int k=0;k<map.shape()[2];k++)
				acc(map[i][j][k]);
	return acc;
}
/**
 * calculate the statistics of a map
 * @param graph  determines connectivity
 * @param map    where the values are stored
 * @return       statistics-object
 */
template<class M>
stat_t voxel_stats(voxelgraph_t& graph, const M& map){
	stat_t acc;
	voxel_vertex_iterator vi, vend;
	foreach ( const voxel_vertex_descriptor& v, vertices(graph))
		acc(map[v]);
	return acc;
}

/**
 * Transform a voxelgraph (=dense grid) to a wurzelgraph (=sparse list)
 * @param vg  voxelgraph (unchanged)
 * @param wg  wurzelgraph (created)
 * @param p_map  predecessor map
 * @param inclusion_map marks nodes which are to be included in wurzelgraph
 */
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
			add_edge(res->second,wv,wg);
		}
	}
	std::cout <<"done ("<<V(num_vertices(wg))<<")"<<std::endl;
}

/**
 * rank nodes of a graph to check whether they are a local maximum
 * @param base  basename where result is stored for visualization
 * @param vg    voxelgraph
 * @param rankmap   ranks are stored here
 * @param valuemap  values are taken from here
 * @param tracesmap only nodes marked in here are considered
 */
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

/**
 * Removes edge sequences starting at leafs from a wurzelgraph 
 * which are less than minlen steps or less than max_radius_mm long long
 * @param wg       the graph
 * @param scale    factor to get from voxel to mm
 * @param max_radius_mm  minimum length of a edge sequence
 */
void erode_tree(wurzelgraph_t& wg, const double& scale, const double& max_radius_mm){
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
			static const unsigned int minlen = 9; // TODO: make a parameter!
			while(cnt++<minlen){
				if(out_degree(pred,wg)>1)
					break;
				tie(ei,eend) = in_edges(pred,wg);
				if(ei!=eend)
					pred = source(*ei,wg);
			}
			if(cnt >= minlen)
				continue;
			if(voxdist(get(vertex_name,wg)[pred].begin(),v.begin())*scale > max_radius_mm)
			        continue;
			clear_vertex(*vi,wg);
			remove_vertex(*vi,wg);
			modified=true;
		}
	}
	std::cout <<"done (num_nodes: "<< num_vertices(wg)<<")."<<std::endl;
}

/**
 * remove nodes which are not leafs and are not in the maximum map 
 * @param wg      the graph
 * @param maxmap  the map telling which nodes to include
 */
template<class T>
void remove_nonmax_nodes(wurzelgraph_t& wg, const T& maxmap){
	std::cout <<"Removing nonmaximum nodes (num_nodes: "<< num_vertices(wg)<<")..."<<std::flush;
	wurzel_vertex_iterator wi,wend,next;
	wurzel_in_edge_iterator ei,eend;
	wurzelg_traits::adjacency_iterator      ai,aend;
	wurzelgraph_t::inv_adjacency_iterator  iai,iaend;
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

/**
 * randomly remove degree-two-nodes, merging their properties into their neighbors
 * @param wg             the graph
 * @param leave_fraction fraction of nodes to be left in the graph
 */
void merge_deg2_nodes(wurzelgraph_t& wg, float leave_fraction){
	std::cout <<"Removing deg2nodes (num_nodes: "<< num_vertices(wg)<<")..."<<std::flush;
	wurzel_vertex_iterator wi,wend,next;
	wurzel_in_edge_iterator ei,eend;
	wurzelg_traits::adjacency_iterator      ai,aend;
	wurzelgraph_t::inv_adjacency_iterator  iai,iaend;

	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	property_map<wurzelgraph_t,vertex_normal_t>::type normal_map = get(vertex_normal, wg);
	property_map<wurzelgraph_t,root_stddev_t>::type stddev_map   = get(root_stddev, wg);
	property_map<wurzelgraph_t,vertex_eigenval_t>::type eigenval_map = get(vertex_eigenval, wg);
	property_map<wurzelgraph_t,vertex_param0_t>::type param0_map = get(vertex_param0, wg);

	bool             changed = true;
	static const float prob  = 0.01f;
	//int              cnt     = 0;
	unsigned int     target  = num_vertices(wg) *leave_fraction;
	while(changed && num_vertices(wg)>target){
		changed = false;
		tie(wi, wend) = vertices(wg);
		for (next=wi; wi != wend;wi=next) {
			next++;
			if(out_degree(*wi,wg) != 1)
				continue;
			if(in_degree(*wi,wg) != 1)
				continue;
			if(drand48()>prob)
				continue;
			tie(ai,aend)   = adjacent_vertices(*wi,wg);
			tie(iai,iaend) = inv_adjacent_vertices(*wi,wg);

			stddev_map[*ai] = 2.f/3.f * stddev_map[*ai] + 1.f/3.f * stddev_map[*wi];
			param0_map[*ai] = 2.f/3.f * param0_map[*ai] + 1.f/3.f * param0_map[*wi];
			pos_map[*ai]    = 2.f/3.f * pos_map[*ai]    + 1.f/3.f * pos_map[*wi];

			stddev_map[*iai] = 2.f/3.f * stddev_map[*iai] + 1.f/3.f * stddev_map[*wi];
			param0_map[*iai] = 2.f/3.f * param0_map[*iai] + 1.f/3.f * param0_map[*wi];
			pos_map[*iai]    = 2.f/3.f * pos_map[*iai]    + 1.f/3.f * pos_map[*wi];
			
			clear_vertex(*wi,wg);
			remove_vertex(*wi,wg);
			add_edge(*iai,*ai,wg);
			changed = true;
		}
	}
	std::cout <<"done (num_nodes: "<< num_vertices(wg)<<")."<<std::endl;
}
/**
 * read a 3D-map from file
 * @param fn filename
 * @param X  size of 1st dimension (innermost)
 * @param Y  size of 2nd dimension 
 * @param Z  size of 3rd dimension
 * @return a reference to the space in memory
 */
template<class T>
boost::multi_array_ref<T, 3> 
read3darray(std::string fn, unsigned int X, unsigned int Y, unsigned int Z){
  std::cout << "Reading `"<<fn<<"', bytes="<<sizeof(T)<<" size="<<X<<","<<Y<<","<<Z<<std::endl;
  std::ifstream dat(fn.c_str(), std::ios::in | std::ios::binary);
  if(!dat.is_open())	{
	  std::cerr << "Could not open `"<<fn<<"' --> exiting."<<std::endl;
	  exit(1);
  }
  T* data = new T[X*Y*Z];
  dat.read((char*)data,X*Y*Z*sizeof(T));
  if(dat.fail()){
	  std::cerr << "Error while reading `"<<fn<<"' --> exiting."<<std::endl;
	  exit(1);
  }
  dat.close();
  return boost::multi_array_ref<T,3>(data,boost::extents[X][Y][Z]);
}

/**
 * fills the vertex_position map from the vertex_name map and
 * sets the normals and eigenvalues maps to zero.
 * @param wg
 */
void
initialize_vertex_positions(wurzelgraph_t& wg){
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	property_map<wurzelgraph_t,vertex_name_t>::type     name_map = get(vertex_name, wg);

	property_map<wurzelgraph_t,vertex_normal_t>::type normal_map = get(vertex_normal, wg);
	property_map<wurzelgraph_t,vertex_eigenval_t>::type eigenval_map = get(vertex_eigenval, wg);

	foreach(wurzel_vertex_descriptor& wv, vertices(wg)){
		vec3_t&   p = pos_map[wv];
		std::copy(name_map[wv].begin(), name_map[wv].end(), p.begin());
		normal_map[wv] = ublas::scalar_matrix<double>(3,3,0);
		eigenval_map[wv] = ublas::scalar_vector<double>(3,0);
	}
}

/**
 * Determine normals and eigenvectors for all vertices. 
 * @param wg the graph with the vertices
 * @param acc an accessor which allows (interpolated) access to a 3D array
 */
template<class T>
void
determine_vertex_normals(wurzelgraph_t& wg, const T& acc){
	//std::cout << "Determining vertex normals..."<<std::flush;
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	property_map<wurzelgraph_t,vertex_normal_t>::type normal_map = get(vertex_normal, wg);
	property_map<wurzelgraph_t,vertex_eigenval_t>::type eigenval_map = get(vertex_eigenval, wg);
	tbb::parallel_for_each(vertices(wg).first,vertices(wg).second,[=](wurzel_vertex_descriptor& wv){
	//foreach(wurzel_vertex_descriptor& wv, vertices(wg)){
		get_normal(normal_map[wv],eigenval_map[wv],pos_map[wv],acc,1.0);
	});
	//std::cout << "done."<<std::endl;
}

/**
 * a visitor which keeps a median statistic and stops when fixed amount of data
 * collected
 */
template<class RawData>
struct bfs_grid_threshold_visitor:public default_bfs_visitor {
  bfs_grid_threshold_visitor(const RawData& wmap, int num, wmedstat_t* m)
  :m_stats(m),m_map(wmap), m_num(num)
  {
  }
  template < typename Vertex, typename Graph >
  void discover_vertex(Vertex u, const Graph & g) 
  {
	  (*m_stats)(m_map[u]);
	  if(count(*m_stats)==m_num){
	          throw 0.0;
	  }
  }
  wmedstat_t* m_stats;
  const RawData& m_map;
  unsigned int      m_num;
};
template<class RawData>
bfs_grid_threshold_visitor<RawData> make_bfs_grid_threshold_visitor(const RawData& wm, int i, wmedstat_t* m){
	return bfs_grid_threshold_visitor<RawData>(wm, i, m);
}

template<class Vertex>
struct implicit_bfs_grid_node{
	Vertex v;
	double cost;
	implicit_bfs_grid_node():cost(0){ }
	implicit_bfs_grid_node(const Vertex& me):v(me),cost(0){ }
	implicit_bfs_grid_node(const Vertex& me, const implicit_bfs_grid_node<Vertex>& pred){ 
		v = me;

		vox2arr<float_grid>& sato = *g_sato;
		double d = voxdist(v.begin(),pred.v.begin());
		double c = exp( - 100.0 * (sato[v] + sato[pred.v]) );
	       	cost = pred.cost + d*c; 
	}
	bool operator==(const implicit_bfs_grid_node& o){ return v==o.v; }
	bool operator<(const implicit_bfs_grid_node& o)const{ return cost>o.cost; } // priority is LESS if cost ics GREATER
};

template<bool direction,class Visitor, class GridGraph, class Vertex, class PredMap>
void implicit_bfs_on_grid(Visitor visit, GridGraph& g, Vertex start, const PredMap& pmap){
	typedef implicit_bfs_grid_node<Vertex> node;
	typedef std::priority_queue<node, std::vector<node> > implicit_bfs_pqueue;
	implicit_bfs_pqueue  queue;
	std::map<Vertex, bool> closed;
	queue.push(node(start));
	while(!queue.empty()){
		node v = queue.top();
		queue.pop();
		if(closed.find(v.v)!=closed.end())
			continue;
		if(v.v != start) // count start node once outside this function, since we start twice from here!
			visit.discover_vertex(v.v, g); 
		closed[v.v]=1;

		typename boost::graph_traits<GridGraph>::out_edge_iterator      oei,oeend;
		tie(oei,oeend)   = out_edges(v.v,g);
		for(; oei!=oeend; oei++){
			const typename boost::graph_traits<GridGraph>::vertex_descriptor adj       = target(*oei,g);
			const typename boost::graph_traits<GridGraph>::vertex_descriptor& pred_v   = pmap[v.v];
			const typename boost::graph_traits<GridGraph>::vertex_descriptor& pred_oei = pmap[adj];
			//std::cout << "         " << V(adj)<< std::endl;
			if( direction && pred_v==adj   && closed.find(adj)==closed.end()) // go towards the shoot of the plant
				queue.push(node(adj,v));
			if(!direction && pred_oei==v.v && closed.find(adj)==closed.end()) // go away from shoot of the plant
				queue.push(node(adj,v));
		}
	}
}

template<class GridGraph, class RawData, class ResultData, class PredMap>
void is_part_of_root(ResultData& res, GridGraph& g, const RawData& raw, const PredMap& pmap, const float start_thresh, const float flow_thresh){

	bool use_edge_detector = false;
	if(use_edge_detector){
		unsigned int smooth = 10;
		//unsigned int cnt=0;
		//foreach(const voxel_vertex_descriptor& v, vertices(g)){
		tbb::parallel_for_each(vertices(g).first,vertices(g).second,[&](voxel_vertex_descriptor& v){

				wmedstat_t s_median1;
				wmedstat_t s_median2;
				if(raw[v] < start_thresh)
				return;
				try{ implicit_bfs_on_grid<1>(make_bfs_grid_threshold_visitor(raw,  smooth,&s_median1), g, v, pmap); // towards shoot
				}catch(double m){}
				if(count(s_median1)!=smooth) return;// return from /lambda/
				try{ implicit_bfs_on_grid<0>(make_bfs_grid_threshold_visitor(raw,  smooth,&s_median2), g, v, pmap); // away from shoot
				}catch(double m){}
				if(count(s_median2)!=smooth) return;// return from /lambda/
				//s_median(raw[v]); // results in odd number of observations: v decides if tied!
				if(median(s_median1)/median(s_median2)>flow_thresh){
				res[v[0]][v[1]][v[2]] = 1;
				}
				//if(count(s_median)>=(unsigned int)(2*smooth) && median(s_median)>flow_thresh){
				//        res[v[0]][v[1]][v[2]] = 1;
				//}
				});
		//};
	}else{
		unsigned int smooth = 10;
		//unsigned int cnt=0;
		//foreach(const voxel_vertex_descriptor& v, vertices(g)){
		tbb::parallel_for_each(vertices(g).first,vertices(g).second,[&](voxel_vertex_descriptor& v){
				//std::cout << "," << std::flush;

				wmedstat_t s_median;
				if(raw[v] < start_thresh)
				return;
				try{ implicit_bfs_on_grid<1>(make_bfs_grid_threshold_visitor(raw,  smooth,&s_median), g, v, pmap); // towards shoot
				}catch(double m){}
				try{ implicit_bfs_on_grid<0>(make_bfs_grid_threshold_visitor(raw,  smooth,&s_median), g, v, pmap); // away from shoot
				}catch(double m){}
				if(count(s_median)<2*smooth) return;// return from /lambda/

				s_median(raw[v]); // results in odd number of observations: v decides if tied!

				if(median(s_median)>flow_thresh){
					res[v[0]][v[1]][v[2]] = 1;
				}
				});
		//};
	}
}

/**
 * a visitor which keeps a median statistic and stops when fixed amount of data
 * collected
 */
template<class WidthMap, class WeightMap>
struct bfs_scale_averager_visitor:public default_bfs_visitor {
  bfs_scale_averager_visitor(const WidthMap& wmap, const WeightMap& weights, int num, wmedstat_t* m, double dummy_scale, double weight_thresh)
  :m_stats(m),m_map(wmap), m_weights(weights), m_num(num), m_dummy_scale(dummy_scale), m_weight_thresh(weight_thresh)
  {
  }
  template < typename Vertex, typename Graph >
  void discover_vertex(Vertex u, const Graph & g) 
  {
	  property_map<wurzelgraph_t,vertex_position_t>::const_type pos_map  = get(vertex_position, g);

	  double width  = m_map[pos_map[u]];
	  double w      = m_weights[pos_map[u]];
	  if(w<m_weight_thresh)
		 width = m_dummy_scale;

	  (*m_stats)(width);
	  //(*m_stats)( width, weight=w);
	  //std::cout << V(count(*m_stats))<<V(median(*m_stats))<<std::endl;
	  if(count(*m_stats)==m_num){
	          throw median(*m_stats);
	  }
  }
  wmedstat_t* m_stats;
  const WidthMap& m_map;
  const WeightMap& m_weights;
  unsigned int      m_num;
  double m_dummy_scale, m_weight_thresh;
};
template<class WidthMap, class WeightMap>
bfs_scale_averager_visitor<WidthMap, WeightMap> make_bfs_scale_averager_visitor(const WidthMap& wm, WeightMap weights, int i, wmedstat_t* m, double dummy_scale, double weight_thresh){
	return bfs_scale_averager_visitor<WidthMap, WeightMap>(wm, weights,i, m,dummy_scale, weight_thresh);
}

/**
 * estimate the root radius for each vertex using the scale at which it was maximal
 * @param wg   the graph with the vertices
 * @param acc  an (interpolated) accessor to an array which contains the sigmas=scale where vesselness was maximal
 * @param max_radius_mm  not used
 * @param wi   contains scale so that radius can be determined in mm
 */
template<class T>
void
determine_radius_from_scales(wurzelgraph_t& wg, const T& acc, const double& max_radius_mm, const wurzel_info& wi){
	//std::cout << "Determining vertex normals..."<<std::flush;
	property_map<wurzelgraph_t,vertex_name_t>::type name_map         = get(vertex_name, wg);
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map      = get(vertex_position, wg);
	property_map<wurzelgraph_t,vertex_normal_t>::type normal_map     = get(vertex_normal, wg);
	property_map<wurzelgraph_t,vertex_eigenval_t>::type eigenval_map = get(vertex_eigenval, wg);

	foreach(wurzel_vertex_descriptor& wv, vertices(wg)){
		//double scale     = acc[name_map[wv]];
		double scale     = acc[pos_map[wv]];
		eigenval_map[wv][1] = scale*wi.scale;
		eigenval_map[wv][2] = scale*wi.scale;
	}
	//std::cout << "done."<<std::endl;
}

/**
 * estimates inertia tensor around each vertex
 * @param wg graph with the vertices
 * @param acc 3D interpolated accessor with "mass" data
 * @param max_radius_mm -- used to determine integration area
 * @param wi determines scale so that radius can be interpreted in voxels
 */
template<class T>
void
determine_inertia_tensor(wurzelgraph_t& wg, const T& acc, const double& max_radius_mm, const wurzel_info& wi){
	//std::cout << "Determining vertex normals..."<<std::flush;
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	property_map<wurzelgraph_t,vertex_normal_t>::type normal_map = get(vertex_normal, wg);
	property_map<wurzelgraph_t,vertex_eigenval_t>::type eigenval_map = get(vertex_eigenval, wg);

	const double r = max_radius_mm/wi.scale, step=1.0; // r and step in voxel

	foreach(wurzel_vertex_descriptor& wv, vertices(wg)){
		get_inertia_tensor(normal_map[wv],eigenval_map[wv],pos_map[wv],acc,r,step,wi.spross_intensity,wi.noise_cutoff);
	}
	//std::cout << "done."<<std::endl;
}

/**
 * move vertex towards the center-of-mass in its vicinity
 * @param wg   graph with vertices
 * @param acc  3D interpolated accessor to "mass" data
 */
template<class T>
void
move_vertex_in_plane(wurzelgraph_t& wg, const T& acc){
	std::cout << "moving vertices in planes..."<<std::flush;
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	property_map<wurzelgraph_t,vertex_normal_t>::type normal_map = get(vertex_normal, wg);
	property_map<wurzelgraph_t,vertex_eigenval_t>::type eigenval_map = get(vertex_eigenval, wg);
	//double sum  = 0.0;
	//int    cnt1 = 0, cnt2 = 0;
	tbb::atomic<int> cnt1,cnt2;
	double sum = 0.0;
	cnt1=0;
	cnt2=0;
	tbb::parallel_for_each(vertices(wg).first,vertices(wg).second,[&](wurzel_vertex_descriptor& wv){
	//foreach(wurzel_vertex_descriptor& wv, vertices(wg)){
		covmat_t& m = normal_map[wv];
		vec3_t&   p = pos_map[wv];

		vec3_t eigval = eigenval_map[wv] + ublas::scalar_vector<double>(3,0.0001);
		static const double thresh = 1.5;
		bool plane_defined = eigval[0] > eigval[1]*thresh
			&&           eigval[0] > eigval[2]*thresh;
		
		if(!plane_defined){
			vec3_t    avgp = ublas::scalar_vector<double>(3,0);
			double    avgw = 0.0;
			int       avgc = 0;
			foreach(const wurzel_edge_descriptor& ve, out_edges(wv,wg)){
				avgp += pos_map[target(ve,wg)];
				avgc += 1;
				avgw += acc[pos_map[target(ve,wg)]];
			}
			foreach(const wurzel_edge_descriptor& ve, in_edges(wv,wg)){
				avgp += pos_map[source(ve,wg)];
				avgc += 1;
				avgw += acc[pos_map[source(ve,wg)]];
			}
			if(avgc==1)
				return;
			avgw /= avgc;
			avgp /= avgc;
			vec3_t diff = 0.5*(avgp-p);
			p += diff;
			sum += ublas::norm_2(diff);
			++cnt1;
		}else{
			double dx1 = acc(p[0] + m(0,1), p[1] + m(1,1), p[2] + m(2,1));
			double dx2 = acc(p[0] - m(0,1), p[1] - m(1,1), p[2] - m(2,1));

			double dy1 = acc(p[0] + m(0,2), p[1] + m(1,2), p[2] + m(2,2));
			double dy2 = acc(p[0] - m(0,2), p[1] - m(1,2), p[2] - m(2,2));
			double norm = dx1+dx2+dy1+dy2 + 0.00001; // avoid div by 0
			double dx = (dx1-dx2)/norm*0.5;
			double dy = (dy1-dy2)/norm*0.5;
			p += dx * ublas::matrix_column<covmat_t>(m,1);
			p += dy * ublas::matrix_column<covmat_t>(m,2);
			sum += (SQR(dx)+SQR(dy));
			++cnt2;
		}
	});
	std::cout << "done (avg norm="<<sum/(cnt1+cnt2)<<"; "<<V(cnt1)<<V(cnt2)<<")"<<std::endl;
}

/**
 * Determine radius of root at each vertex.
 * 
 * The function fits a gaussian bell-shaped function to the data around 
 * each voxel, projected to the plane orthogonal to its "normal"
 * 
 * @param wg   graph vertices are taken from
 * @param acc  3D interpolated accessor to mass data
 * @param scale factor voxel->mm
 * @param max_radius_mm maximum radius of a root in mm
 * @param wi   additional root information
 */
template<class T>
void
wurzel_thickness(wurzelgraph_t& wg, const T& acc, const double& scale, const double max_radius_mm, const wurzel_info& wi){
	stat_t s_thickness;
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	property_map<wurzelgraph_t,vertex_normal_t>::type normal_map = get(vertex_normal, wg);
	property_map<wurzelgraph_t,root_stddev_t>::type stddev_map   = get(root_stddev, wg);
	property_map<wurzelgraph_t,vertex_eigenval_t>::type eigenval_map = get(vertex_eigenval, wg);
	property_map<wurzelgraph_t,vertex_param0_t>::type param0_map = get(vertex_param0, wg);
	const double r = max_radius_mm/scale, step=0.5; // r and step in voxel
	double sr      = r/5.0; // start radius for fitting
	std::cout << "Determining Wurzel thickness..."
		<< V(max_radius_mm)
		<< V(pow(2*r/step,2))
		<<std::flush;
	// foreach(wurzel_vertex_descriptor& wv, vertices(wg){
	tbb::parallel_for_each(vertices(wg).first,vertices(wg).second,[=](wurzel_vertex_descriptor& wv){
		covmat_t& m = normal_map[wv];
		vec3_t&   p = pos_map[wv];

		std::vector<double> values, coors;
		values.reserve(2*r/step*2*r/step*3);
		coors .reserve(2*r/step*2*r/step*3);
		ublas::vector<double> params(2);

		double centerval = acc( p[0], p[1], p[2])/wi.spross_intensity;
		params[0] = centerval;
		params[1] = 1.0/(2*sr*sr);
		//params[2] = 0;
		for(double z = -1; z <= 1; z+=1){
			for(double i = -r; i <= r; i+=step){
				for (double j = -r; j <= r; j+=step)
				{
					double val = acc(
							p[0] + z*m(0,0) + i*m(0,1) + j*m(0,2),
							p[1] + z*m(1,0) + i*m(1,1) + j*m(1,2),
							p[2] + z*m(2,0) + i*m(2,1) + j*m(2,2))/wi.spross_intensity;
					values.push_back(std::max(0.0,val-wi.noise_cutoff/4));
					//values.push_back(val<wi.noise_cutoff ? 0.0 : val);
					//values.push_back(val);
					coors.push_back(i*i + j*j);
				}
			}
		}
		/* 
		 * sadly, levmar is not thread-safe.
		 *
		 * we therefore do not use fit_gauss_curve anymore, and instead resort to
		 * an old RPROP implementation of mine, which gives the same results.
		 *
		 */
		//fit_gauss_curve(params,&values.front(),values.size(),&coors.front()); 
		//std::vector<float> sqerrs;
		RProp rprop(2);
		ublas::vector<double>& grad = rprop.getGrad();
		double last_sqerr = 0.0;
	       for(int iter=0;iter<60;iter++){
		     grad *= 0;
		     double sqerr = 0.0;
		     for(std::vector<double>::iterator vit=values.begin(),cit=coors.begin(); vit!=values.end();vit++,cit++){
			     double r = params[0]*exp(-params[1]* *cit) - *vit;
			     grad[0]  += exp(-params[1]* *cit) * r;
			     grad[1]  += -(*cit)*params[0] * exp(-params[1]* *cit) * r;
			     sqerr+= r*r;
		     }
		     grad /= values.size();
		     static const int wd = 0.01; // weight decay
		     grad[0] -= wd * (params[0]-0.08);
		     grad[1] -= wd * (params[1]-0);
		     //sqerrs.push_back(sqerr);
		     rprop.update_irprop_plus(sqerr<last_sqerr ? RProp::ARPROP_DIR_OK : RProp::ARPROP_DIR_WRONG);
		     last_sqerr = sqerr;
		     params    -= rprop.getDeltaW();
		     if(ublas::norm_2(rprop.getDeltaW()) < 0.0001)
			     break;
		     params[0]  = std::max(0.08,params[0]);  // centerscale must have minimum height to bound expscale
		     params[1]  = std::max(0.0,params[1]);  // expscale    must be at least 0
		     params[0]  = std::min(1.0,params[0]);  // centerscale must not be larger than spross_intensity
	       }

	       //std::cout << "iter: "<<0             <<"  sqerr: "<<sqerrs[0]<<std::endl
	       //         << "      "<<sqerrs.size()-1<<"  sqerr: "<<sqerrs.back()<<std::endl
	       //         << "  p0: "<<params[0]<<std::endl
	       //         << "  p1: "<<params[1]<<std::endl;
	       

		//centerscale += wi.noise_cutoff; // account for cutoff above
		double stddev = sqrt(1.0/(2*params[1])); 

		if(stddev != stddev)          stddev = 0.00001;
		if(isinf(stddev))             stddev = 0.00001;
		if(stddev > 2*max_radius_mm)  stddev = 2*max_radius_mm;
		if(stddev <        0.00001  ) stddev = 0.00001;
		stddev_map[wv] = stddev;
		param0_map[wv] = params[0];
		//s_thickness(stddev);
	});
	//foreach(wurzel_vertex_descriptor& wv, vertices(wg)){
	//        covmat_t& m = normal_map[wv];
	//        vec3_t&   p = pos_map[wv];

	//        values.clear();
	//        coors.clear();
	//        double centerval = acc( p[0], p[1], p[2])/wi.spross_intensity;
	//        params[0] = centerval;     double& centerscale = params[0];
	//        params[1] = 1.0/(2*sr*sr); double& expscale    = params[1];
	//        params[2] = 0;
	//        for(double i = -r; i <= r; i+=step){
	//                for (double j = -r; j <= r; j+=step)
	//                {
	//                        double val = acc(
	//                                        p[0]+i*m(0,1) + j*m(0,2),
	//                                         p[1]+i*m(1,1) + j*m(1,2),
	//                                         p[2]+i*m(2,1) + j*m(2,2))/wi.spross_intensity;
	//                        values.push_back(std::max(0.0,val-wi.noise_cutoff/4));
	//                        //values.push_back(val<wi.noise_cutoff ? 0.0 : val);
	//                        //values.push_back(val);
	//                        coors.push_back(i*i + j*j);
	//                }
	//        }
		
	//        fit_gauss_curve(params,&values.front(),values.size(),&coors.front());
	//        //centerscale += wi.noise_cutoff; // account for cutoff above
	//        double stddev = sqrt(1.0/(2*expscale)); 


	//        //std::cout << "-------------------------------------"<<std::endl;
	//        //for(int iter=0;iter<20;iter++){
	//        //       double g =0.0;
	//        //       for(std::vector<double>::iterator vit=values.begin(),cit=coors.begin(); vit!=values.end();vit++,cit++){
	//        //               g += -(*cit)*params[0] * exp(-params[1]* *cit) * *vit;
	//        //       }
	//        //       g /= values.size();
	//        //       params[1] -= 0.2 * g;
	//        //       //std::cout << "g: "<<g<<std::endl;
	//        //}
	//        //double stddev = sqrt(1.0/(2*params[1])) * scale; 


	//        //vec3_t& l = eigenval_map[wv];
	//        //double stddev = (sqrt(l(1))+sqrt(l(2)))/2.0;
	//        if(stddev != stddev)          stddev = 0.00001;
	//        if(isinf(stddev))             stddev = 0.00001;
	//        if(stddev > 2*max_radius_mm)  stddev = 2*max_radius_mm;
	//        if(stddev <        0.00001  ) stddev = 0.00001;
	//        //stddev *= params[0];
	//        stddev_map[wv] = stddev;
	//        param0_map[wv] = centerscale;
	//        s_thickness(stddev);
	//}
	//std::cout <<V(s_thickness)<< " done."<<std::endl;
}

/**
 * Adjust radius at each node such that the roots only get a larger or equal
 * diameter towards the root, never smaller.
 * @param wg
 */
void
smooth_thickness(wurzelgraph_t& wg){
	std::cout << "Smoothing Wurzel thickness..."<<std::flush;
	property_map<wurzelgraph_t,root_stddev_t>::type stddev_map   = get(root_stddev, wg);
	foreach(wurzel_vertex_descriptor& wv, vertices(wg)){
		if(out_degree(wv,wg)!=0)
			continue;
		wurzel_vertex_descriptor d = wv;
		double radius = stddev_map[d];
		while(in_degree(d,wg)>0){
			if(stddev_map[d]<radius)
				stddev_map[d] = radius;
			radius = std::max(stddev_map[d],radius);
			d = source(*in_edges(d,wg).first,wg);
		}
	}
}

/**
 * Determine the position of the plant stem in a given image plane
 * @param dims shape of the 3D data
 * @param spross_intensity returns the intensity at the determined position
 * @param acc  accessor to raw 3D data (not interpolated)
 * @param plane  index of plane which stem grows through
 * @param axis   index of dimension which contains plane
 * @return the position of the stem within the given plane
 */
template<class T>
boost::array<vidx_t, 3> 
locate_stem(boost::array<vidx_t,3>& dims, double& spross_intensity, const T& acc, unsigned int plane, unsigned int axis){
	std::cout <<"Locate stem..."<<std::flush;
	voxel_vertex_descriptor arg_max_val = {{0,0,0}};
	float max_val = -1E6;

	voxel_vertex_descriptor current = {{0,0,0}};
	current[axis] = plane;
	int ax1 =  axis==0            ?1:0;
	int ax2 = (axis==0 || axis==1)?2:1;

	for ( current[ax1] = 0; current[ax1] < dims[ax1]; ++current[ax1])
	{
		for ( current[ax2] = 0; current[ax2] < dims[ax2]; ++current[ax2])
		{
			double val = 0;
			const int radius = 6;
			voxel_vertex_descriptor win;
			win[axis] = plane;
			double sum = 0.0;
			for ( win[ax1]  = std::max(0,(int)current[ax1]-radius);
			      win[ax1] <= std::min(dims[ax1]-1, current[ax1]+radius);
			      ++win[ax1])
			for ( win[ax2]  = std::max(0,(int)current[ax2]-radius);
			      win[ax2] <= std::min(dims[ax2]-1, current[ax2]+radius);
			      ++win[ax2]){
				float fact = SQR(current[0]-win[0]) + SQR(current[1]-win[1]) + SQR(current[2]-win[2]);
				fact       = 1.0/(2.0*M_PI) * exp(-fact);
				      double sval = acc[win];
				      val += fact * sval;
				      sum += fact;
			      }
			val /= sum;
			if(val>max_val){
				max_val = val;
				arg_max_val = current;
			}
		}
	}
	std::cout <<"done ("<<arg_max_val[0]<<", "<<arg_max_val[1]<<", "<<arg_max_val[2]<<"), val="<<max_val<<std::endl;
	spross_intensity = max_val;
	return arg_max_val;
}


int main(int argc, char* argv[]) {
	wurzel_info info;

	// read commandline parameters
	po::variables_map vm = get_config(info,argc,argv);
	static const unsigned int X=info.X,Y=info.Y,Z=info.Z,XYZ=info.XYZ;
	bool force_recompute_dijkstra = vm.count("force");
	std::string base = (fs::path(info.directory) / vm["base"].as<std::string>()).string();

	// Define a grid_graph where dimensions don't wrap
	boost::array<vidx_t, 3> lengths = { { X, Y, Z } };
	voxelgraph_t graph(lengths, false); // no dim is wrapped

	// read upsampled raw data
	float_grid Raw  = read3darray<float>(getfn(base,"upsampled","dat"),X,Y,Z);
	boost::array<vidx_t, 3> strunk   // locate the stem
		=  locate_stem(lengths, info.spross_intensity, make_vox2arr(Raw), info.stem_plane, info.stem_axis); //{ { 109, 129,  24 } };

	//std::cout << "Barley settings!!!"<<std::endl;
	//info.spross_intensity = 425000;
	std::cout << "Lupine settings!!!"<<std::endl;
	info.spross_intensity = 0.9;

	// read additional data: Vesselness measure and (if desired) eigenvalues and scales
	float_grid Sato = read3darray<float>(getfn(base,"sato","dat"),X,Y,Z); g_sato = new vox2arr<float_grid>(Sato);
	//float_grid ev10 = read3darray<float>(getfn(base,"","ev10"),X,Y,Z); g_ev10 = new vox2arr<float_grid>(ev10);
	//float_grid ev11 = read3darray<float>(getfn(base,"","ev11"),X,Y,Z); g_ev11 = new vox2arr<float_grid>(ev11);
	//float_grid ev12 = read3darray<float>(getfn(base,"","ev12"),X,Y,Z); g_ev12 = new vox2arr<float_grid>(ev12);
	float_grid scales = read3darray<float>(getfn(base,"scales","dat"),X,Y,Z); 

	std::cout << "Sato stats in file: " << voxel_stats(graph,make_vox2arr(Sato)) <<std::endl;

	g_strunk = strunk; // ugh. ugly: use global variable.

	// allocate helper arrays: Path
	unsigned char* paths = new unsigned char[XYZ];
	std::fill(paths, paths+XYZ, (unsigned char) 0);
	uc_grid Paths(paths,boost::extents[X][Y][Z]);

	// allocate helper arrays: ranks
	unsigned char* ranks = new unsigned char[XYZ];
	std::fill(ranks, ranks+XYZ, (unsigned char) 255);
	uc_grid Ranks(ranks,boost::extents[X][Y][Z]);

	// allocate helper arrays: flow
	float* flow = new float[XYZ];
	std::fill(flow, flow+XYZ, (float) 0);
	float_grid Flow(flow,boost::extents[X][Y][Z]);


	// normalize vesselness such that it is >=0 and max is 1
	stat_t s_raw  = voxel_stats(graph, make_vox2arr(Raw));
	foreach(const voxel_vertex_descriptor& v, vertices(graph)) {
		float& f = Sato[v[0]][v[1]][v[2]];
		//f = std::max(0.f, f-0.2f) * 1.f/0.7f;
		//f  = exp(-10.f * log(1.f+f)/log(2.f) );
		//f  = exp(-15.f * f );

		float& g = Raw[v[0]][v[1]][v[2]];
		//f *= 1.0f+2.0f*g;
		//      g  = (g-min(s_raw))/(max(s_raw)-min(s_raw));
		f += std::max(g,0.f)/max(s_raw);
	}
	stat_t s_sato = voxel_stats(graph, make_vox2arr(Sato));
	std::cout << "Raw:  " << s_raw <<std::endl;
	std::cout << "Sato: " << s_sato <<std::endl;

	predecessor_map_t         p_map(num_vertices(graph), get(vertex_index, graph)); 
	distance_map_t            d_map(num_vertices(graph), get(vertex_index, graph)); 

	// Run dijkstra to find shortest paths
	{       boost::prof::profiler prof("Dijkstra");
		find_shortest_paths(base,graph,strunk,p_map,d_map,force_recompute_dijkstra);
	}
	std::cout << "Dmap: " << voxel_stats(graph,d_map) <<std::endl;

	std::cout << "Determining scaling factors..." <<std::endl;
	stat_t s_allpaths = voxel_stats(graph,d_map);

	/*******************************************************
	 *
	 * Trace paths to determine what is soil/noise and what is root
	 *
	 * *****************************************************/
	std::cout << "Tracing paths... " <<std::flush;
	double start_threshold       = vm["start-threshold"].as<double>();
	double total_len_perc_thresh = vm["total-len-frac"].as<double>();
	double avg_len_perc_thresh   = vm["avg-len-frac"].as<double>();
	double min_flow_thresh       = vm["min-flow-thresh"].as<double>() * (info.XYZ/info.read_XYZ/2);
	bool   no_gauss_fit          = vm.count("no-gauss-fit");
	bool   no_subpix_pos         = vm.count("no-subpix-pos");

	start_threshold *= info.noise_cutoff;

	voxelg_traits::vertices_size_type strunk_idx = boost::get(vertex_index, graph, strunk);
	stat_t s_avg_pathlen, s_pathlen, s_cnt, s_flow;
	vox2arr<float_grid> vox2raw(Raw);
	vox2arr<float_grid> vox2sato(Sato);
	{ boost::prof::profiler prof("Tracing");
		foreach (const voxel_vertex_descriptor& v, vertices(graph)) {
			// determine total path length statistic
			//if(vox2sato[v] < start_threshold)
			if(vox2raw[v]/info.spross_intensity < start_threshold)
				continue;
			s_pathlen(d_map[v]);
		}
		std::cout << ". " <<std::flush;

		property_map<voxelgraph_t,vertex_index_t>::type vertex_index_map = get(vertex_index, graph);
		foreach (const voxel_vertex_descriptor& v, vertices(graph)) {
			// determine avg path costs statistic
			//if(vox2sato[v] < start_threshold)
			if(vox2raw[v]/info.spross_intensity < start_threshold)
				continue;
			if(((d_map[v]-min(s_pathlen))/(max(s_pathlen)-min(s_pathlen))) > total_len_perc_thresh)
				continue;
			double vox_dist = 0.0;
			voxel_vertex_descriptor v2 = v;
			unsigned int cnt = 0;
			float flow_add = vox2raw[v]/info.spross_intensity;
			voxel_vertex_descriptor tmp;
			while(cnt++<XYZ){ 
				Flow[v2[0]][v2[1]][v2[2]] += flow_add;
				//Paths[v2[0]][v2[1]][v2[2]] = 255;
				if(vertex_index_map[v2] == strunk_idx)
					break;
				tmp = p_map[v2];
				vox_dist += voxdist(tmp.begin(), v2.begin());
				v2 = tmp;
			}
			if(cnt>=XYZ){
				std::cout << "endless loop!"<<std::endl;
				exit(1);
			}
			s_cnt(vox_dist);
			if(vox_dist>0)
				s_avg_pathlen(d_map[v]/vox_dist);
		}
		std::cout << ". " <<std::flush;
		//write_voxelgrid<unsigned char>(getfn(base,"paths1","dat"),graph,make_vox2arr(Paths));
		//std::fill(paths, paths+XYZ, (unsigned char) 0);

		foreach (const voxel_vertex_descriptor& v, vertices(graph)) {
			// determine flow statistic
			//if(vox2sato[v] < start_threshold)
			if(vox2raw[v]/info.spross_intensity < start_threshold)
				continue; // too weak at start
			if(((d_map[v]-min(s_pathlen))/(max(s_pathlen)-min(s_pathlen))) > total_len_perc_thresh)
				continue; // too long
			s_flow(Flow[v[0]][v[1]][v[2]]);
		}
	} // profiling: Tracing

	std::cout << "Total   Pathlens: "<< s_pathlen<<std::endl;
	std::cout << "Average Pathlens: "<< s_avg_pathlen<<std::endl;
	std::cout << "Hop     Pathlens: "<< s_cnt<<std::endl;

	foreach (const voxel_vertex_descriptor& v0, vertices(graph)) {
		//if(vox2sato[v0] < start_threshold)
		if(vox2raw[v0]/info.spross_intensity < start_threshold)
			continue;                  // weak signal at start point
		float total_dist   = d_map[v0];
		if(((total_dist-min(s_pathlen))/(max(s_pathlen)-min(s_pathlen))) > total_len_perc_thresh)
			continue;                  // too long
		float flow = Flow[v0[0]][v0[1]][v0[2]];
		//if((flow-min(s_flow))/(max(s_flow)-min(s_flow)) < min_flow_thresh)
		if(flow < min_flow_thresh * info.noise_cutoff)
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
		v = v0;
		while(1){ 
			Paths[v[0]][v[1]][v[2]] = 255;
			if(boost::get(vertex_index,graph,v) == strunk_idx)
				break;
			v = p_map[v];
		}
	}


	// find local ranks
	//rank_op(base,graph,make_vox2arr(Ranks),make_vox2arr(Sato),make_vox2arr(Paths));

	const double maximum_radius_mm = 1;
	wurzelgraph_t wgraph;
	// transform voxelgrid to adjacency-list graph
	{ 	boost::prof::profiler prof("voxelgraph-to-root-graph");
		paths2adjlist(graph,wgraph,p_map,make_vox2arr(Paths));
	}

	// remove short sequences at leafs
	{ 	boost::prof::profiler prof("erode tree");
		erode_tree(wgraph, info.scale, maximum_radius_mm);
	}

	// find positions of voxels which can be more precise than the 
	// voxel resolution
	{ 	boost::prof::profiler prof("subpixel-positioning");
		std::cout << "Finding subpixel vertex positions..."<<std::flush;
		initialize_vertex_positions(wgraph);
		if(!no_subpix_pos){
			for(int i=0;i<30;i++){
				determine_vertex_normals(wgraph, make_vox2arr_subpix(Sato));
				move_vertex_in_plane(wgraph, make_vox2arr_subpix(Sato));
			}
		}
		std::cout << "done."<<std::endl;
	}

	// locally fit gaussian functions to the root to determine its radius
	{ 	boost::prof::profiler prof("gauss-fitting");
		if(!no_gauss_fit)
			wurzel_thickness(wgraph, make_vox2arr_subpix(Raw), info.scale, maximum_radius_mm, info);
	}

	// substitute covariance stuff with inertia tensor, for radius estimation!
	{ 	boost::prof::profiler prof("inertia-tensor");
		//determine_inertia_tensor(wgraph, make_vox2arr_subpix(Raw), maximum_radius_mm, info);
		determine_radius_from_scales(wgraph, make_vox2arr_subpix(scales), maximum_radius_mm, info);
	}

	// reduce amount of nodes (good if the output needs to be passed on to another program)
	//merge_deg2_nodes(wgraph,0.25);

	if(1){ // serialize tree
		std::ofstream ofs_serializer(getfn(base,"wgraph","ser").c_str());
		boost::archive::text_oarchive oa(ofs_serializer);
		oa << wgraph;
	}

	//smooth_thickness(wgraph);


	//remove_nonmax_nodes(wgraph,make_vox2arr(Ranks));
	//write_voxelgrid<unsigned char>(getfn(base,"paths","dat"),graph,make_vox2arr(Paths));
	//write_voxelgrid<unsigned char>(getfn(base,"ranks","dat"),graph,make_vox2arr(Ranks));
}
