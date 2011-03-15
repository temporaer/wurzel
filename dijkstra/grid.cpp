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
#include <boost/range/irange.hpp>
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
#include "voxelgrid.hpp"
#include "wurzel_tree.hpp"
#include "voxel_accessors.hpp"
#include "voxel_normals.hpp"
#include "wurzel_info.hpp"
#include "grid_config.hpp"
#include "gaussian_fit.hpp"

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



vox2arr<float_grid>* g_sato, *g_ev10, *g_ev11, *g_ev12;
boost::array<vidx_t, 3> g_strunk;

// map from edges to weights
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
			//if(voxdist(get(vertex_name,wg)[pred].begin(),v.begin())*scale > max_radius_mm)
			//        continue;
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

template<class T>
void
determine_vertex_normals(wurzelgraph_t& wg, const T& acc){
	//std::cout << "Determining vertex normals..."<<std::flush;
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	property_map<wurzelgraph_t,vertex_normal_t>::type normal_map = get(vertex_normal, wg);
	property_map<wurzelgraph_t,vertex_eigenval_t>::type eigenval_map = get(vertex_eigenval, wg);
	foreach(wurzel_vertex_descriptor& wv, vertices(wg)){
		get_normal(normal_map[wv],eigenval_map[wv],pos_map[wv],acc);
	}
	//std::cout << "done."<<std::endl;
}

template<class T>
void
determine_inertia_tensor(wurzelgraph_t& wg, const T& acc, const double& max_radius_mm, const wurzel_info& wi){
	//std::cout << "Determining vertex normals..."<<std::flush;
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	property_map<wurzelgraph_t,vertex_normal_t>::type normal_map = get(vertex_normal, wg);
	property_map<wurzelgraph_t,vertex_eigenval_t>::type eigenval_map = get(vertex_eigenval, wg);

	const double r = max_radius_mm/wi.scale, step=0.25; // r and step in voxel

	foreach(wurzel_vertex_descriptor& wv, vertices(wg)){
		get_inertia_tensor(normal_map[wv],eigenval_map[wv],pos_map[wv],acc,r,step);
	}
	//std::cout << "done."<<std::endl;
}

template<class T>
void
move_vertex_in_plane(wurzelgraph_t& wg, const T& acc){
	std::cout << "moving vertices in planes..."<<std::flush;
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	property_map<wurzelgraph_t,vertex_normal_t>::type normal_map = get(vertex_normal, wg);
	double sum = 0.0;
	int    cnt = 0;
	foreach(wurzel_vertex_descriptor& wv, vertices(wg)){
		covmat_t& m = normal_map[wv];
		vec3_t&   p = pos_map[wv];
		
		double dx1 = acc(p[0] + m(0,1), p[1] + m(1,1), p[2] + m(2,1));
		double dx2 = acc(p[0] - m(0,1), p[1] - m(1,1), p[2] - m(2,1));
                                                               
		double dy1 = acc(p[0] + m(0,2), p[1] + m(1,2), p[2] + m(2,2));
		double dy2 = acc(p[0] - m(0,2), p[1] - m(1,2), p[2] - m(2,2));
		double norm = dx1+dx2+dy1+dy2;
		double dx = (dx1-dx2)/norm*0.5;
		double dy = (dy1-dy2)/norm*0.5;
		p += dx * ublas::matrix_column<covmat_t>(m,1);
		p += dy * ublas::matrix_column<covmat_t>(m,2);
		sum += SQR(dx)+SQR(dy);
		cnt ++;
	}
	std::cout << "done (avg norm="<<sum/cnt<<")"<<std::endl;
}

template<class T>
void
wurzel_thickness(wurzelgraph_t& wg, const T& acc, const double& scale, const double max_radius_mm, const wurzel_info& wi){
	std::cout << "Determining Wurzel thickness..."<<std::flush;
	stat_t s_thickness;
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	property_map<wurzelgraph_t,vertex_normal_t>::type normal_map = get(vertex_normal, wg);
	property_map<wurzelgraph_t,root_stddev_t>::type stddev_map   = get(root_stddev, wg);
	property_map<wurzelgraph_t,vertex_eigenval_t>::type eigenval_map = get(vertex_eigenval, wg);
	property_map<wurzelgraph_t,vertex_param0_t>::type param0_map = get(vertex_param0, wg);
	const double r = max_radius_mm/scale, step=0.25; // r and step in voxel
	double sr      = r/5.0; // start radius for fitting
	std::vector<double> values, coors;
	values.reserve(2*r/step*2*r/step);
	coors .reserve(2*r/step*2*r/step);
	double params[3]; // gaussian curve params
	foreach(wurzel_vertex_descriptor& wv, vertices(wg)){
		covmat_t& m = normal_map[wv];
		vec3_t&   p = pos_map[wv];
		values.clear();
		coors.clear();
		double centerval = acc( p[0], p[1], p[2])/wi.spross_intensity;
		params[0] = centerval;     double& centerscale = params[0];
		params[1] = 1.0/(2*sr*sr); double& expscale    = params[1];
		params[2] = 0;
		for(double i = -r; i <= r; i+=step){
			for (double j = -r; j <= r; j+=step)
			{
				double val = acc(
						p[0]+i*m(0,1) + j*m(0,2),
					 	p[1]+i*m(1,1) + j*m(1,2),
					 	p[2]+i*m(2,1) + j*m(2,2))/wi.spross_intensity;
				values.push_back(std::max(0.0,val-wi.noise_cutoff));
				//values.push_back(val<wi.noise_cutoff ? 0.0 : val);
				//values.push_back(val);
				coors.push_back(i*i + j*j);
			}
		}
		
		fit_gauss_curve(params,&values.front(),values.size(),&coors.front());
		centerscale += wi.noise_cutoff/4; // account for cutoff above
		double stddev = sqrt(1.0/(2*expscale)); 


		//std::cout << "-------------------------------------"<<std::endl;
		//for(int iter=0;iter<20;iter++){
		//       double g =0.0;
		//       for(std::vector<double>::iterator vit=values.begin(),cit=coors.begin(); vit!=values.end();vit++,cit++){
		//               g += -(*cit)*params[0] * exp(-params[1]* *cit) * *vit;
		//       }
		//       g /= values.size();
		//       params[1] -= 0.2 * g;
		//       //std::cout << "g: "<<g<<std::endl;
		//}
		//double stddev = sqrt(1.0/(2*params[1])) * scale; 


		//vec3_t& l = eigenval_map[wv];
		//double stddev = (sqrt(l(1))+sqrt(l(2)))/2.0;
		if(stddev != stddev)          stddev = 0.00001;
		if(isinf(stddev))             stddev = 0.00001;
		if(stddev > max_radius_mm)    stddev = max_radius_mm;
		if(stddev <        0.00001  ) stddev = 0.00001;
		//stddev *= params[0];
		stddev_map[wv] = stddev;
		param0_map[wv] = centerscale;
		s_thickness(stddev);
	}
	std::cout <<V(s_thickness)<< " done."<<std::endl;
}

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

	po::variables_map vm = get_config(info,argc,argv);
	static const unsigned int X=info.X,Y=info.Y,Z=info.Z,XYZ=info.XYZ;
	bool force_recompute_dijkstra = vm.count("force");
	std::string base = vm["base"].as<std::string>();
	
  // Define a 3x5x7 grid_graph where the second dimension doesn't wrap
  boost::array<vidx_t, 3> lengths = { { X, Y, Z } };
  voxelgraph_t graph(lengths, false); // no dim is wrapped

  float_grid Raw  = read3darray<float>(getfn(base,"upsampled","dat"),X,Y,Z);
  boost::array<vidx_t, 3> strunk   // locate the stem
	  =  locate_stem(lengths, info.spross_intensity, make_vox2arr(Raw), info.stem_plane, info.stem_axis); //{ { 109, 129,  24 } };
  info.spross_intensity = 425000;

  float_grid Sato = read3darray<float>(getfn(base,"","sato"),X,Y,Z); g_sato = new vox2arr<float_grid>(Sato);
  //float_grid ev10 = read3darray<float>(getfn(base,"","ev10"),X,Y,Z); g_ev10 = new vox2arr<float_grid>(ev10);
  //float_grid ev11 = read3darray<float>(getfn(base,"","ev11"),X,Y,Z); g_ev11 = new vox2arr<float_grid>(ev11);
  //float_grid ev12 = read3darray<float>(getfn(base,"","ev12"),X,Y,Z); g_ev12 = new vox2arr<float_grid>(ev12);

  std::cout << "Sato stats in file: " << voxel_stats(graph,make_vox2arr(Sato)) <<std::endl;

  g_strunk = strunk;

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
	  //f *= 1.0f+2.0f*g;
	  //      g  = (g-min(s_raw))/(max(s_raw)-min(s_raw));
		f += g/max(s_raw);
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

  std::cout << "Tracing paths... " <<std::flush;
  double start_threshold       = vm["start-threshold"].as<double>();
  double total_len_perc_thresh = vm["total-len-frac"].as<double>();
  double avg_len_perc_thresh   = vm["avg-len-frac"].as<double>();
  double min_flow_thresh       = vm["min-flow-thresh"].as<double>();

  start_threshold *= info.noise_cutoff;

  voxelg_traits::vertices_size_type strunk_idx = boost::get(vertex_index, graph, strunk);
  stat_t s_avg_pathlen, s_pathlen, s_cnt, s_flow;
  vox2arr<float_grid> vox2raw(Raw);
  vox2arr<float_grid> vox2sato(Sato);
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
	  float flow_add = vox2raw[v];
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

  const double maximum_radius_mm = 4;
  wurzelgraph_t wgraph;
  paths2adjlist(graph,wgraph,p_map,make_vox2arr(Paths));
  erode_tree(wgraph, info.scale, maximum_radius_mm);
  initialize_vertex_positions(wgraph);

  std::cout << "Finding subpixel vertex positions..."<<std::flush;
  for(int i=0;i<30;i++){
          determine_vertex_normals(wgraph, make_vox2arr_subpix(Sato));
          move_vertex_in_plane(wgraph, make_vox2arr_subpix(Sato));
  }
  std::cout << "done."<<std::endl;
  wurzel_thickness(wgraph, make_vox2arr_subpix(Raw), info.scale, maximum_radius_mm, info);

  // substitute covariance stuff with inertia tensor, for radius estimation!
  determine_inertia_tensor(wgraph, make_vox2arr_subpix(Raw), maximum_radius_mm, info);

  if(1){
	  // serialize tree
	  std::ofstream ofs_serializer(getfn(base,"wgraph","ser").c_str());
	  boost::archive::text_oarchive oa(ofs_serializer);
	  oa << wgraph;
  }

  //smooth_thickness(wgraph);


  //remove_nonmax_nodes(wgraph,make_vox2arr(Ranks));

  //merge_deg2_nodes(wgraph);

  //write_voxelgrid<unsigned char>(getfn(base,"paths","dat"),graph,make_vox2arr(Paths));
  //write_voxelgrid<unsigned char>(getfn(base,"ranks","dat"),graph,make_vox2arr(Ranks));
}
