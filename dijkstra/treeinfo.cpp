#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/unordered_map.hpp>

#include <boost/multi_array.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/covariance.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/reverse_graph.hpp>

#include <boost/filesystem.hpp>

#include <boost/numeric/ublas/io.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "inc/voxel_accessors.hpp"
#include "inc/wurzel_tree.hpp"
#include "inc/treeinfo_config.hpp"

#define V(X) #X<<":"<<(X)<<"  "

namespace fs = boost::filesystem;
using namespace boost::accumulators;
typedef accumulator_set< double, features< tag::min, tag::mean, tag::max, tag::variance, tag::count > > stat_t;
typedef accumulator_set< double, features< tag::median, tag::count, tag::mean > > medstat_t;

std::map<void*,double> g_diameter;
typedef boost::multi_array_ref<float, 3> float_grid;


// (a2b3 − a3b2, a3b1 − a1b3, a1b2 − a2b1)
vec3_t cross_product(
        const vec3_t &a, const vec3_t &b)
{
        vec3_t c(3);
        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];
        return c;
}

/**
 * construct a file name from the parts in the parameters
 * @param base   basename of fiel w/o ext
 * @param marker is attached to base with a dash in between
 * @param ext    is the extention of the generated filename
 * @return new complete filename
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
 * maps edges to weights
 * this is only here because otherwise we can't link
 */
voxel_edge_weight_map::reference 
voxel_edge_weight_map::operator[](key_type e) const {
	voxel_vertex_descriptor s = source(e,m_graph);
	voxel_vertex_descriptor t = target(e,m_graph);
	return 0.0;
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
  std::cout << "# Reading `"<<fn<<"', bytes="<<sizeof(T)<<" size="<<X<<","<<Y<<","<<Z<<std::endl;
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

// for ground truth graphs
void determine_mass_from_vol(std::string basename,wurzelgraph_t& g, const wurzel_info& info){
	auto   edge_mass_map = get(edge_mass, g);
	auto      radius_map = get(vertex_radius, g);
	auto        pos_map  = get(vertex_position, g);

	foreach(const wurzel_edge_descriptor& e, edges(g)){
		edge_mass_map[e] = 0.0;
		const vec3_t& pes = pos_map[source(e,g)]; // source of e
		const vec3_t& pet = pos_map[target(e,g)]; // target of e
		vec3_t e_dir      = pet-pes;               // normalized direction of e
		double h = ublas::norm_2(e_dir);
		double R = radius_map[source(e,g)];
		double r = radius_map[target(e,g)];
		double V = h*M_PI/3.0 * (R*R + R*r + r*r);
		edge_mass_map[e] = V;
	}
}
void determine_mass_from_raw(std::string basename,wurzelgraph_t& g, const wurzel_info& info){
	// NOTE: this assumes that the root has VOXEL coordinates, NOT mm!
	static const unsigned int X=info.X,Y=info.Y,Z=info.Z;
	//boost::array<vidx_t, 3> lengths = { { X, Y, Z } };

	fs::path datadir = info.directory;
	float_grid Raw  = read3darray<float>(getfn((datadir / basename).string(),"upsampled","dat"),X,Y,Z);
	property_map<wurzelgraph_t,vertex_name_t>::type name_map  = get(vertex_name, g);
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, g);
	property_map<wurzelgraph_t,vertex_radius_t>::type radius_map = get(vertex_radius, g);
	property_map<wurzelgraph_t,marked_vertex_t>::type   mark_map = get(marked_vertex, g);
	property_map<wurzelgraph_t,edge_mass_t>::type   edge_mass_map = get(edge_mass, g);
	auto raw_acc = make_vox2arr(Raw);
	static const double sample_step_len = 0.2; // mm
	foreach(const wurzel_edge_descriptor& e, edges(g)){
		edge_mass_map[e] = 0.0;
		const vec3_t& pes = pos_map[source(e,g)]; // source of e
		const vec3_t& pet = pos_map[target(e,g)]; // target of e
		vec3_t e_dir      = pet-pes;               // normalized direction of e
		e_dir            /= ublas::norm_2(e_dir);

		std::vector<std::pair<vec3_t,vec3_t> > sub_edges;
		{
			double len_e     = ublas::norm_2(pes-pet);
			double first_dist = std::min(len_e, sample_step_len); // if shorter than sample_step_len, use whole edge
			double dist       = first_dist;
			vec3_t prev_pv   = pes;                    // the `last' vertex we compared to
			for(; dist<len_e; dist+=sample_step_len){
				const vec3_t pv = pes + dist * e_dir;
				sub_edges .push_back(std::make_pair(prev_pv,pv));
				prev_pv = pv;
			}
			const vec3_t pv = pes + std::max(dist-sample_step_len,0.0) * e_dir;
			sub_edges.push_back(std::make_pair(pv,pet));
		}
		// calculate two vectors orthogonal to e
		vec3_t tmp = e_dir;
		tmp[0] += 1;
		vec3_t v1 = cross_product(tmp,e_dir);
		vec3_t v2 = cross_product(v1,e_dir);
		v1 /= ublas::norm_2(v1);
		v2 /= ublas::norm_2(v2);

		double edgemass_sum = 0.;
		for(unsigned int subedge=0;subedge<sub_edges.size();subedge++){
			const vec3_t prev_pv = sub_edges[subedge].first;
			const vec3_t p      = sub_edges[subedge].second;
			// p now is the center of a slice we want to sample
			double dist   = ublas::norm_2(p-pes)/ublas::norm_2(pet-pes);
			double radius =    dist  * radius_map[source(e,g)]  +
					(1-dist) * radius_map[target(e,g)];

			
			for(double dx=-sqrt(2.0)*radius; dx<=sqrt(2.0)*radius;dx+=0.5){
				for(double dy=-sqrt(2.0)*radius; dy<=sqrt(2.0)*radius;dy+=0.5){
					if(sqrt(dx*dx + dy*dy)>1.1*radius)
						continue;
					vec3_t pr = p + dx*v1+dy*v2;
					// for rounding: add 0.5 to all components
					pr[0] += 0.5;
					pr[1] += 0.5;
					pr[2] += 0.5;
                    pr[0]  = std::max(0.0, std::min(X-1.0, pr[0]));
                    pr[1]  = std::max(0.0, std::min(Y-1.0, pr[1]));
                    pr[2]  = std::max(0.0, std::min(Z-1.0, pr[2]));
					float& r = Raw[(unsigned int)pr[0]][(unsigned int)pr[1]][(unsigned int)pr[2]];
					edgemass_sum += r;
					r = 0;
				}
			}
		}
		edge_mass_map[e] = edgemass_sum;
	}
}

/**
 * load a serialized tree from a file
 */
void loadtree(std::string basename, wurzelgraph_t& g){
	std::ifstream ifs((basename+"/wgraph.ser").c_str());
	boost::archive::text_iarchive ia(ifs);
	ia >> g;
}
/**
 * load a ground truth tree from a file
 */
void load_gt_tree(std::string basename, wurzelgraph_t& g, const wurzel_info& wi){
	property_map<wurzelgraph_t,vertex_position_t>::type  pos_map = get(vertex_position, g);
	property_map<wurzelgraph_t,vertex_radius_t>::type radius_map = get(vertex_radius, g);

	fs::path p(basename);

	std::string fn = (p / std::string("rootsegmentdata.dat")).string();
	std::cout << "# reading data from "<<fn<<std::endl;
	std::ifstream ifs(fn.c_str());
	char str[512];
	ifs.getline(str,512); // header
	typedef std::map<wurzel_vertex_descriptor,int> next_t;
	typedef std::map<int,wurzel_vertex_descriptor> index_t;
	next_t next_map;
	index_t index_map;
	int idx = -1;
	bool is_stub=false;
	bool is_leaf=false;

	vec3_t offset;
	double read_xdim   = wi.read_XYZ / 832 / 64;
	double offset_fact = (256.0 - read_xdim)/(256-64); // 1 if read_xdim==64, 0 if read_xdim==256
	if(read_xdim != 64)
		offset_fact = 0.0;
	offset[0] = offset_fact * 0.60 / wi.scale; // in voxels
	offset[1] = offset_fact * 0.00 / wi.scale; // in voxels
	offset[2] = offset_fact * 0.60 / wi.scale; // in voxels

	std::cout << "# offset: "<<offset<<std::endl;

	while(!ifs.eof()){
		idx++;
		float x,y,z, diameter,length, nx,ny,nz, len2;
		int prev, next;
		ifs  >> x >> y >> z;
		ifs  >> prev >> next;
		ifs  >> nx >> ny >> nz;
		ifs  >> diameter >> length >> len2;
		//std::cout << V(x) << V(y)<<V(z)<<V(next)<<V(diameter)<<V(length)<<std::endl;
		is_leaf = next==-1;
		if(is_leaf && is_stub){
			is_stub = false;
			continue;
		}
		
		if(length<0.01 && !is_leaf){
			is_stub = true;
			continue; // less than a millimeter long roots are stubs that should be ignored
		}


		wurzel_vertex_descriptor v = add_vertex(g);    
		pos_map[v][0] = z*10./wi.scale + 17 + 5.5*10./wi.scale; // cm --> mm, swap x/z, -->voxel
		pos_map[v][1] = y*10./wi.scale + 6  + 40.5*10./wi.scale; // cm --> mm, -->voxel
		pos_map[v][2] = x*10./wi.scale + 17 + 5.5*10./wi.scale; // cm --> mm, swap x/z, -->voxel
		pos_map[v] += offset;
		radius_map[v] = diameter / 2. * 10. / wi.scale; // cm --> mm, -->voxel
		next_map[v] = next;
		index_map[idx] = v;
	}
	for (index_t::iterator it = index_map.begin(); it != index_map.end(); ++it)
	{
		//int current_idx = it->first;
		wurzel_vertex_descriptor current_v = it->second;
		if(!current_v)
			continue;
		int next_idx = next_map[current_v];
		if(next_idx < 0) 
			continue;
		wurzel_vertex_descriptor next_v = index_map[next_idx];
		if(next_v)
			//add_edge(next_v,current_v,g);
			add_edge(current_v,next_v,g);
	}
	std::cout << "# done "<<fn<<std::endl;

	auto inside_bb_check 	= [&](const wurzel_vertex_descriptor& v)->bool{
		double x = pos_map[v][0];
		double y = pos_map[v][1];
		double z = pos_map[v][2];
		if(x<0   +17) return false;
		if(x>=221+17) return false;
		if(y<0   + 6) return false;
		if(y>=821+ 6) return false;
		if(z<0   +17) return false;
		if(z>=221+17) return false;
		return true;
	};

	bool changed=true;
	while(changed){
		changed = false;
		foreach(wurzel_edge_descriptor e, edges(g)){
			bool is_s_inside = inside_bb_check(source(e,g));
			bool is_t_inside = inside_bb_check(target(e,g));
			if(!is_s_inside || !is_t_inside){
				if(!is_s_inside){
					clear_vertex(source(e,g),g);
					remove_vertex(source(e,g),g);
				}
				if(!is_t_inside){
					clear_vertex(target(e,g),g);
					remove_vertex(target(e,g),g);
				}
				changed = true;
				break;
			}
		}
	}
}

/**
 * key-value pair for XML
 */
template<class T>
struct kv_helper{
  std::string k;
  T v;
};

/**
 * key and value constructor for XML
 *
 * this is abbreviated, imagine it is called make_kv_helper(...)
 * to match with the standard naming convention :)
 */
template<class T>
kv_helper<T> kv(std::string s, T o){ kv_helper<T> x;x.k=s; x.v=o; return x;}

/**
 * print a key-value pair to a stream
 */
template<class T>
std::ostream& operator<<(std::ostream& o, kv_helper<T> kvh){
	o << " "<<kvh.k<<"=\""<<kvh.v<<"\" ";
	return o;
}

/**
 * print a vertex to an XML file, includin properties
 */
void to_xml_node(std::ofstream& o, unsigned int& idx, const wurzel_vertex_descriptor& v, const wurzelgraph_t& g){
	property_map<wurzelgraph_t,vertex_index_t>::const_type  index_map = get(vertex_index, g);
	property_map<wurzelgraph_t,vertex_position_t>::const_type  pos_map = get(vertex_position, g);
	property_map<wurzelgraph_t,root_stddev_t>::const_type   stddev_map = get(root_stddev, g);
	property_map<wurzelgraph_t,edge_mass_t>::const_type       mass_map = get(edge_mass, g);
	property_map<wurzelgraph_t,vertex_param0_t>::const_type param0_map = get(vertex_param0, g);

	o << "<Node "
	  << kv("id", idx)
	  << kv("bo", 1)
	  << kv("rad", g_diameter[v]/200.0) // meters, i suppose...
	  << kv("x", pos_map[v][2]/100.0-0.5)
	  << kv("y", pos_map[v][1]/100.0-0.5)
	  << kv("z", pos_map[v][0]/100.0-0.5)
	  << ">"<<std::endl;
	foreach(const wurzel_vertex_descriptor& w, adjacent_vertices(v,g)){
		to_xml_node(o,++idx,w,g);
	}
	o << "</Node>"<<std::endl;
}

/**
 * print a wurzel graph (=a single tree) to an XML-Forest file
 */
void to_xml_forest(std::string basename, const wurzelgraph_t& g){
	// find root
	wurzel_vertex_descriptor root;
	bool found=false;
	foreach(const wurzel_vertex_descriptor &v, vertices(g)){
		if(in_degree(v,g)==0){
			root = v;
			found = true;
			break;
		}
	}
	assert(found);

	std::ofstream ofs((basename+".xml").c_str());
	ofs <<"<?xml version=\"1.0\" ?>"<<std::endl;
	ofs << std::setprecision(10)
	    <<"<Forest "
            <<kv("id",848)
            <<kv("bo",95) 
	    <<kv("cm",1023)
	    <<kv("cn",8)
	    <<kv("d",1.000000)
	    <<kv("dlx",132.0)
	    <<kv("dly",100.0)
	    <<kv("dlz",100.000000)
	    <<kv("sl",1.000000)
	    <<kv("zFac",-1.000000)
	    <<kv("iso",40.0)
	    <<kv("dataset",basename)
	    <<kv("transferFunc","data/lup_red.xml")
	    <<">"
	    <<std::endl;
	ofs << "<Tree id=\"1\">"<<std::endl;
	unsigned int idx=0;
	to_xml_node(ofs,idx,root,g);
	ofs << "</Tree>"<<std::endl;
	ofs<<"/Forest>"<<std::endl;
}

/**
 * determine extent of bounding box for g
 */
void graph_stats(wurzelgraph_t& g){
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, g);

	vec3_t xmax = ublas::scalar_vector<double>(3,-1E6);
	vec3_t xmin = ublas::scalar_vector<double>(3,1E6);
	vec3_t ymax = ublas::scalar_vector<double>(3,-1E6);
	vec3_t ymin = ublas::scalar_vector<double>(3,1E6);
	vec3_t zmax = ublas::scalar_vector<double>(3,-1E6);
	vec3_t zmin = ublas::scalar_vector<double>(3,1E6);
	foreach(const wurzel_vertex_descriptor& v1, vertices(g)){
		const vec3_t& pv1 = pos_map[v1];
		if(xmax[0]<pv1[0])xmax = pv1;
		if(xmin[0]>pv1[0])xmin = pv1;
		if(ymax[1]<pv1[1])ymax = pv1;
		if(ymin[1]>pv1[1])ymin = pv1;
		if(zmax[2]<pv1[2])zmax = pv1;
		if(zmin[2]>pv1[2])zmin = pv1;
	}
	std::cout << "xmin= "<< xmin<<std::endl;
	std::cout << "xmax= "<< xmax<<std::endl;
	std::cout << "ymin= "<< ymin<<std::endl;
	std::cout << "ymax= "<< ymax<<std::endl;
	std::cout << "zmin= "<< zmin<<std::endl;
	std::cout << "zmax= "<< zmax<<std::endl << std::endl;
}
/**
 * calulate mass statistics
 */
double total_mass(wurzelgraph_t& wg, const wurzel_info& wi){
	double sum = 0.0;
	double minpos = /*wi.stem_plane*/14 * wi.X/410.0; // TODO: make this work for other datasets!!
	property_map<wurzelgraph_t,vertex_position_t>::type  pos_map = get(vertex_position, wg);
	property_map<wurzelgraph_t,root_stddev_t>::type   stddev_map = get(root_stddev, wg);
	property_map<wurzelgraph_t,edge_mass_t>::type       mass_map = get(edge_mass, wg);
	property_map<wurzelgraph_t,vertex_param0_t>::type param0_map = get(vertex_param0, wg);
	double sum_horiz = 0.0, sum_vert = 0.0;
	foreach(const wurzel_edge_descriptor& e, edges(wg)){
		const wurzel_vertex_descriptor& s = source(e,wg);
		const wurzel_vertex_descriptor& t = target(e,wg);
		if(pos_map[s][wi.stem_axis] < minpos)
			continue;                                      
		if(pos_map[t][wi.stem_axis] < minpos)
			continue;
		const double length  = voxdist(pos_map[s].begin(),pos_map[t].begin());          // pos_map is in mm

		sum += mass_map[e] * length;

		const vec3_t& spos = pos_map[s];
		const vec3_t& tpos = pos_map[t];
		const vec3_t  diff = tpos-spos;
		bool is_vert = diff[0]*diff[0] > (diff[1]*diff[1]+diff[2]*diff[2]);
		if(is_vert) sum_vert  += mass_map[e];
		else        sum_horiz += mass_map[e];
	}
	double weight_scale     = 1.0000; // milli-gram
	//weight_scale /= wi.X*wi.Y*wi.Z  / 192.0 / 192.0 / 410.0; // change 0.27 so that it fits the resolution
	std::cout << "mass_total= "<<sum * weight_scale <<std::endl;
	std::cout << "mass_total_horiz= "<<sum_horiz * weight_scale <<std::endl;
	std::cout << "mass_total_mass= "<<sum_vert  * weight_scale <<std::endl;
	std::cout << "mass_V_H_ratio= "<<sum_vert/sum_horiz<<std::endl;
	return sum*weight_scale;
}

/**
 * a segment is a sequence of connected edges without crossings
 */
struct wurzel_segment{
	const wurzelgraph_t& wg;
	wurzel_segment(const wurzelgraph_t& g):wg(g){}
	std::list<wurzel_edge_descriptor> edges;
	void mark_non_ending_vertices(wurzelgraph_t& g, bool average_radius_in_segment){
		property_map<wurzelgraph_t,vertex_position_t>::type pos_map   = get(vertex_position, g);
		property_map<wurzelgraph_t,marked_vertex_t>::type   mark_map  = get(marked_vertex, g);
		property_map<wurzelgraph_t,vertex_radius_t>::type radius_map  = get(vertex_radius, g);
		property_map<wurzelgraph_t,edge_mass_t>::type mass_map = get(edge_mass, g);
		std::vector<wurzel_vertex_descriptor> ends;
		foreach(const wurzel_edge_descriptor &e, edges){
			const wurzel_vertex_descriptor& s = source(e,g);
			const wurzel_vertex_descriptor& t = target(e,g);

			mark_map[s] = false;
			mark_map[t] = false;

			if(out_degree(s,g)+in_degree(s,g)!=2) { ends.push_back(s); }
			if(out_degree(t,g)+in_degree(t,g)!=2) { ends.push_back(t); }
		}
		if(ends.size() != 2){
			std::cout << "#  Weeeeeeeird segment!"<<std::endl;
			return;
		}
		foreach(const wurzel_edge_descriptor &e, edges){
			const wurzel_vertex_descriptor& s = source(e,g);
			const wurzel_vertex_descriptor& t = target(e,g);

			foreach(const wurzel_vertex_descriptor& v, ends){
				double ds = voxdist(pos_map[s].begin(),pos_map[v].begin()); // assumed to be in mm here!!
				double dt = voxdist(pos_map[t].begin(),pos_map[v].begin()); // assumed to be in mm here!!
				if(ds>2) mark_map[s]=true;
				if(dt>2) mark_map[t]=true;
			}
		}
		
		medstat_t stats;
		foreach(const wurzel_edge_descriptor &e, edges){
			const wurzel_vertex_descriptor& s = source(e,g);
			const wurzel_vertex_descriptor& t = target(e,g);
			if(mark_map[s]) stats(radius_map[s]);
			if(mark_map[t]) stats(radius_map[t]);
		}
		foreach(const wurzel_edge_descriptor &e, edges){
			const wurzel_vertex_descriptor& s = source(e,g);
			const wurzel_vertex_descriptor& t = target(e,g);
			if(count(stats)>0){
				bool set_s=false;
				bool set_t=false;
				foreach(const wurzel_vertex_descriptor& v, ends){
					if(v==s){
						radius_map[s] = std::max((double)radius_map[s], (double)mean(stats));
						set_s = true;
					}
					if(v==t){
						radius_map[t] = std::max((double)radius_map[t], (double)mean(stats));
						set_t = true;
					}
				}
				if(!set_s)
					radius_map[s] = mean(stats);
				if(!set_t)
					radius_map[t] = mean(stats);
			}
		}
	}
	void average_mass_in_segment(wurzelgraph_t& g){
		double edge_mass_sum = 0.0;
		double total_len     = 0.0;
		property_map<wurzelgraph_t,vertex_position_t>::type pos_map   = get(vertex_position, g);
		property_map<wurzelgraph_t,edge_mass_t>::type mass_map = get(edge_mass, g);
		foreach(const wurzel_edge_descriptor &e, edges){
			const wurzel_vertex_descriptor& s = source(e,g);
			const wurzel_vertex_descriptor& t = target(e,g);

			edge_mass_sum += mass_map[e];
			total_len     += voxdist(pos_map[s].begin(),pos_map[t].begin());
		}
		foreach(const wurzel_edge_descriptor &e, edges){
			//const wurzel_vertex_descriptor& s = source(e,g);
			//const wurzel_vertex_descriptor& t = target(e,g);
			//double len = voxdist(pos_map[s].begin(),pos_map[t].begin());
			//mass_map[e] = edge_mass_sum / total_len * len;
			mass_map[e] = edge_mass_sum / total_len; // converts mass to units per length(!)
		}
	}
	void add(const wurzel_edge_descriptor& e){
		edges.push_back(e);
	}
	bool has(const wurzel_edge_descriptor& e){
		return edges.end()==std::find(edges.begin(),edges.end(),e);
	}
	bool has(const wurzel_vertex_descriptor& v){
		foreach(const wurzel_edge_descriptor &e, edges){
			if(source(e,wg)==v) return true;
			if(target(e,wg)==v) return true;
		}
		return false;
	}
	bool is_vert(){
		property_map<wurzelgraph_t,vertex_position_t>::const_type pos_map  = get(vertex_position, wg);
		const wurzel_vertex_descriptor& s = source(edges.front(),wg);
		const wurzel_vertex_descriptor& t = target(edges.back(),wg);
		const vec3_t& spos = pos_map[s];
		const vec3_t& tpos = pos_map[t];
		const vec3_t& diff = tpos-spos;
		return diff[0]*diff[0] > (diff[1]*diff[1]+diff[2]*diff[2]);
	}
	double length(){
		double sum = 0.0;
		property_map<wurzelgraph_t,vertex_position_t>::const_type pos_map  = get(vertex_position, wg);
		foreach(const wurzel_edge_descriptor &e, edges){
			const wurzel_vertex_descriptor& s = source(e,wg);
			const wurzel_vertex_descriptor& t = target(e,wg);
			double length = voxdist(pos_map[s].begin(),pos_map[t].begin()); // assumed to be in mm here!!
			sum += length;
		}
		return sum;
	}
	bool has_leaf(){
		foreach(const wurzel_edge_descriptor &e, edges){
			const wurzel_vertex_descriptor& s = source(e,wg);
			const wurzel_vertex_descriptor& t = target(e,wg);
			if(out_degree(s,wg)==0) return true;
			if(out_degree(t,wg)==0) return true;
		}
		return false;
	}
	int pieces(){
		return edges.size();
	}
};

/**
 * recursively find all segments in a wurzelgraph
 */
void build_seg(const wurzel_edge_descriptor& e, wurzel_segment seg, std::list<wurzel_segment>& segments, const wurzelgraph_t& wg){
	const wurzel_vertex_descriptor& t = target(e,wg);
	if(out_degree(t,wg)==0){// ends in leaf
		seg.add(e);
		segments.push_back(seg);
		return;
	}
	if(out_degree(t,wg)>1){ // ends in furcation
		seg.add(e);
		segments.push_back(seg);
		// continue building new segments from nodes after t
		foreach(const wurzel_edge_descriptor& e2, out_edges(t,wg)){
			build_seg(e2,wurzel_segment(wg),segments,wg);
		}
		return;
	}
	// now out_degree==1
	seg.add(e);
	build_seg(*out_edges(t,wg).first, seg,segments,wg);
}

void determine_nodes_far_away_from_crossings(wurzelgraph_t& wg, bool average_radius_in_segment, bool average_mass_in_segment){
	std::list<wurzel_segment> segments;
	std::map<wurzel_edge_descriptor,bool> foundmap;
	wurzel_vertex_descriptor root;
	bool found = false;
	foreach(const wurzel_vertex_descriptor& v, vertices(wg)){
		if(in_degree(v,wg)==0){
			root = v;
			found = true;
			break;
		}
	}
	if(!found){
		std::cout << "Could not find root node...."<<std::endl;
		exit(1);
	}
	foreach(const wurzel_edge_descriptor& e2, out_edges(root,wg)){
		build_seg(e2,wurzel_segment(wg), segments, wg);
	}
	foreach(wurzel_segment& seg,segments){
		seg.mark_non_ending_vertices(wg,average_radius_in_segment);
		if(average_mass_in_segment)
			seg.average_mass_in_segment(wg);
	}
	
}

/**
 * determine average segment length
 */
double avg_len_btw_splits(wurzelgraph_t& wg){
	std::list<wurzel_segment> segments;
	std::map<wurzel_edge_descriptor,bool> foundmap;
	//property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	wurzel_vertex_descriptor root;
	bool found = false;
	foreach(const wurzel_vertex_descriptor& v, vertices(wg)){
		if(in_degree(v,wg)==0){
			root = v;
			found = true;
			break;
		}
	}
	if(!found){
		std::cout << "Could not find root node...."<<std::endl;
		exit(1);
	}
	foreach(const wurzel_edge_descriptor& e2, out_edges(root,wg)){
		build_seg(e2,wurzel_segment(wg), segments, wg);
	}
	
	
	std::cout << "num_segments= "<<segments.size()<<std::endl;
	stat_t s_lens, s_lens_hor, s_lens_ver;
	stat_t s_nums;
	foreach(wurzel_segment& seg,segments){
		if(seg.has_leaf())
			continue;
		s_lens(seg.length());
		s_nums(seg.pieces());
		if(seg.is_vert())
			s_lens_ver(seg.length());
		else
			s_lens_hor(seg.length());
	}
	std::cout <<"seg_avg_len_mm= "<< mean(s_lens)<<std::endl;
	std::cout <<"seg_avg_len_std_mm= "<<sqrt(variance(s_lens))<<std::endl;
	std::cout <<"seg_avg_len_vertical_mm=   "<< mean(s_lens_ver)<<std::endl;
	std::cout <<"seg_avg_len_vertical_std_mm = "<<sqrt(variance(s_lens_ver))<<std::endl;
	std::cout <<"seg_avg_len_horizontal_mm= "<< mean(s_lens_hor)<<std::endl;
	std::cout <<"seg_avg_len_horizontal_std_mm="<<sqrt(variance(s_lens_hor))<<std::endl;
	std::cout <<"seg_avg_pieces= "<< mean(s_nums)<<std::endl;
	std::cout <<"seg_avg_std_pieces="<<sqrt(variance(s_nums))<<std::endl;
	return mean(s_lens);
}

/**
 * determine total length
 */
double total_length(wurzelgraph_t& wg){
	double sum = 0.0, sum_horiz=0.0, sum_vert=0.0;
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	foreach(const wurzel_edge_descriptor& e, edges(wg)){
		const wurzel_vertex_descriptor& s = source(e,wg);
		const wurzel_vertex_descriptor& t = target(e,wg);
		double length = voxdist(pos_map[s].begin(),pos_map[t].begin()); // assumed to be in mm here!!
		sum += length;

		const vec3_t& spos = pos_map[s];
		const vec3_t& tpos = pos_map[t];
		const vec3_t& diff = tpos-spos;
		bool is_vert = diff[0]*diff[0] > (diff[1]*diff[1]+diff[2]*diff[2]);
		if(is_vert) sum_vert  += length;
		else        sum_horiz += length;
	}
	std::cout << "len_total= "<<sum <<std::endl;
	std::cout << "len_total_horiz= "<<sum_horiz <<std::endl;
	std::cout << "len_total_vertical= "<<sum_vert  <<std::endl;
	std::cout << "len_V_H_ratio= "<<sum_vert/sum_horiz<<std::endl;

	return sum;
}
/**
 * determine distance between two graphs
 * @param g1 1st graph
 * @param g2 2nd graph
 * @param tolerance  distance up to which nodes are considered "the same"
 * @param verbose show progress
 */
void distance(const std::string& ident, std::ofstream& ofs, std::ofstream& scat, wurzelgraph_t& g1, wurzelgraph_t& g2, const double& tolerance, bool verbose){
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map1  = get(vertex_position, g1);
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map2  = get(vertex_position, g2);
	property_map<wurzelgraph_t,vertex_radius_t>::type radius_map1 = get(vertex_radius, g1);
	property_map<wurzelgraph_t,marked_vertex_t>::type   mark_map1 = get(marked_vertex, g1);
	property_map<wurzelgraph_t,marked_vertex_t>::type   mark_map2 = get(marked_vertex, g2);
	property_map<wurzelgraph_t,vertex_radius_t>::type radius_map2 = get(vertex_radius, g2);
	auto   edge_mass_map1 = get(edge_mass, g1);
	auto   edge_mass_map2 = get(edge_mass, g2);
	int no_counterpart = 0;
	int v_cnt = 0;
	stat_t s_distances, s_masses;
	accumulator_set<double, stats<tag::variance, tag::covariance<double,tag::covariate1> > > s_radius1;
	accumulator_set<double, stats<tag::variance, tag::covariance<double,tag::covariate1> > > s_radius2;
	double sum_length_no_counterpart=0.0;
	double sum_length_matched       =0.0;
	double sum_mass_no_counterpart=0.0;
	double sum_mass_matched       =0.0;
	vec3_t offsets = ublas::scalar_vector<double>(3,0);

	static const double sample_step_len = 0.2; // mm
	foreach(const wurzel_edge_descriptor& e1, edges(g1)){
		const vec3_t& pe1s = pos_map1[source(e1,g1)]; // source of e1
		const vec3_t& pe1t = pos_map1[target(e1,g1)]; // target of e1
		vec3_t e1_dir      = pe1t-pe1s;               // normalized direction of e1
		e1_dir            /= ublas::norm_2(e1_dir);
		bool use_radius1 = mark_map1[source(e1,g1)] && mark_map1[target(e1,g1)];

		std::vector<std::pair<vec3_t,vec3_t> > sub_edges;
		{
			double len_e1     = ublas::norm_2(pe1s-pe1t);
			double first_dist = std::min(len_e1, sample_step_len); // if shorter than sample_step_len, use whole edge
			double dist       = first_dist;
			vec3_t prev_pv1   = pe1s;                    // the `last' vertex we compared to
			for(; dist<len_e1; dist+=sample_step_len){
				const vec3_t pv1 = pe1s + dist * e1_dir;
				sub_edges .push_back(std::make_pair(prev_pv1,pv1));
				prev_pv1 = pv1;
			}
			const vec3_t pv1 = pe1s + std::max(dist-sample_step_len,0.0) * e1_dir;
			sub_edges.push_back(std::make_pair(pv1,pe1t));
			//std::cout << "splitting a line of len="<<len_e1<< " to piecnum="<<sub_edges.size()<< V(first_dist)<<std::endl;
		}
		for(unsigned int subedge=0;subedge<sub_edges.size();subedge++){
			const vec3_t prev_pv1 = sub_edges[subedge].first;
			const vec3_t pv1      = sub_edges[subedge].second;
			double min_d          = 1E6;
			double radius2        = 1E6;
			bool use_radius2      = false;
			double mass_radius2          = 1E6;
			wurzel_edge_descriptor best_e2;
			foreach(const wurzel_edge_descriptor& e2, edges(g2)){
				const vec3_t& pe2s = pos_map2[source(e2,g2)];
				const vec3_t& pe2t = pos_map2[target(e2,g2)];
				double dp          = ublas::inner_prod(pv1-pe2s, pe2t-pe2s);
				double  t          = dp / ublas::inner_prod(pe2t-pe2s,pe2t-pe2s);
				double d, r;
				if(t<0){
					d = ublas::norm_2(pv1-pe2s);
					if(min_d>d){
						r = radius_map2[source(e2,g2)];
					}
				}
				else if(t>1){
					d = ublas::norm_2(pv1-pe2t);
					if(min_d>d){
						r = radius_map2[target(e2,g2)];
					}
				}
				else{
					vec3_t ul = (pe2t-pe2s)/ublas::norm_2(pe2t-pe2s);
					vec3_t w  = pv1-pe2s;
					d = ublas::norm_2( w - ul * ublas::inner_prod(w,ul));
					if(min_d>d){
						double dist2 = ublas::norm_2(pv1-pe1s)/ublas::norm_2(pe1t-pe1s);
						r =        dist2  * radius_map2[source(e2,g2)]  +
							(1-dist2) * radius_map2[target(e2,g2)];
					}
				}
				if(min_d>d){
					min_d   = d;
					radius2 = r;
					use_radius2 = mark_map2[source(e2,g2)] && mark_map2[target(e2,g2)];
					best_e2 = e2;

					double m = edge_mass_map2[e2];
					mass_radius2   = sqrt(m / M_PI); // m has been divided by len already!
				}
			}
			double dist = ublas::norm_2(pv1-pe1s)/ublas::norm_2(pe1t-pe1s);
			double radius1 =    dist  * radius_map1[source(e1,g1)]  +
					(1-dist) * radius_map1[target(e1,g1)];
			double mass_radius1 = sqrt(edge_mass_map1[e1]/M_PI);
			if(min_d>tolerance){
				ofs << pv1[0]<<"\t"<<pv1[1]<<"\t"<<pv1[2]<<"\t"<<prev_pv1[0]<<"\t"<<prev_pv1[1]<<"\t"<<prev_pv1[2]<<"\t"<<radius1<<"\t"<<std::endl;
				no_counterpart ++;
				double len = ublas::norm_2(pv1-prev_pv1);
				sum_length_no_counterpart += len;
				sum_mass_no_counterpart += len * edge_mass_map1[e1];
			}else{
				if(ublas::norm_2(pv1-pos_map2[source(best_e2,g2)]) < ublas::norm_2(pv1-pos_map2[target(best_e2,g2)]))
					offsets += pv1 - pos_map2[source(best_e2,g2)];
				else
					offsets += pv1 - pos_map2[target(best_e2,g2)];
				double length_matched = ublas::norm_2(pv1-prev_pv1);
				sum_length_matched   += length_matched;
				sum_mass_matched     += length_matched * edge_mass_map1[e1];
				s_distances(min_d);
				if(use_radius1 && use_radius2){
					scat << radius1 << "\t"<<radius2 <<"\t" << length_matched << "\t"<< length_matched*M_PI*radius1*radius1 << "\t"<<length_matched*M_PI*radius2*radius2 << "\t" << mass_radius1 << "\t"<<mass_radius2  << std::endl;
					s_radius1(radius1, boost::accumulators::covariate1=radius2);
					s_radius2(radius2, boost::accumulators::covariate1=radius1);
				}
			}
			v_cnt ++;
		}

	}
	offsets /= count(s_distances);
	std::cout << ident<<"_tolerance = "<<tolerance<<std::endl;
	std::cout << ident<<"_vertices_wo_counterpart = "<<no_counterpart<<std::endl;
	std::cout << ident<<"_mass_wo_counterpart = "<<sum_mass_no_counterpart<<std::endl;
	std::cout << ident<<"_mass_matched = "<<sum_mass_matched<<std::endl;
	std::cout << ident<<"_length_wo_counterpart = "<<sum_length_no_counterpart<<std::endl;
	std::cout << ident<<"_length_matched = "<<sum_length_matched<<std::endl;
	std::cout << ident<<"_distance_mean = "<<mean(s_distances)<<std::endl;
	std::cout << ident<<"_distance_std = "<<sqrt(variance(s_distances))<<std::endl;
	std::cout << ident<<"_offsets = "<<offsets[0]<<", "<<offsets[1]<<", "<<offsets[2]<<std::endl;

	
	std::cout << ident<<"_radius_corr1 = "<<covariance(s_radius1)/(sqrt(variance(s_radius1))*sqrt(variance(s_radius2)))<<std::endl;
	std::cout << ident<<"_radius_corr2 = "<<covariance(s_radius2)/(sqrt(variance(s_radius1))*sqrt(variance(s_radius2)))<<std::endl;
	std::cout << ident<<"_cnt          = "<<v_cnt<<std::endl;
}

/**
 * scale all lengths to mm
 */
void scale_to_mm(wurzelgraph_t& g, const wurzel_info& wi){
	std::cout << "# scaling: "<< wi.scale <<std::endl;
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, g);
	property_map<wurzelgraph_t,vertex_radius_t>::type radius_map = get(vertex_radius, g);
	property_map<wurzelgraph_t,root_stddev_t>::type stddev_map = get(root_stddev, g);
	auto   edge_mass_map = get(edge_mass, g);
	foreach(wurzel_vertex_descriptor& v, vertices(g)){
		pos_map[v]    *= wi.scale;
		radius_map[v] *= wi.scale;
		stddev_map[v] *= wi.scale;
	}
	foreach(wurzel_edge_descriptor e, edges(g)){
		edge_mass_map[e] *= wi.scale * wi.scale * wi.scale;
	}
}

/**
 * print the total length
 */
void action_length(std::vector<wurzel_info>& wis, std::vector<std::string>& bases, const po::variables_map& vm){
	std::string base = bases[0];

	wurzelgraph_t g;
	loadtree(base, g);

	scale_to_mm(g,wis[0]);

	std::cout << "total length("<<base<<"): "<< total_length(g)<<std::endl;

}
/**
 * print the total mass
 */
void action_mass(std::vector<wurzel_info>& wis, std::vector<std::string>& bases, const po::variables_map& vm){
	std::string base = bases[0];

	wurzelgraph_t g;
	loadtree(base, g);

	scale_to_mm(g,wis[0]);

	std::cout << "total mass ("<<base<<"): "<< total_mass(g, wis[0])<<std::endl;

}
/**
 * print the distance btw two graphs (non-symmetric!!!)
 */
void action_distance(std::vector<wurzel_info>& wis, std::vector<std::string>& bases, const po::variables_map& vm){
	fs::path datadir = wis[0].directory;
	std::string base1 = bases[0];
	std::string base2 = bases[1];
	double tolerance = vm["tolerance"].as<double>();

	wurzelgraph_t g1, g2;
	if(base1.find("snr")==std::string::npos) load_gt_tree((datadir/base1).string(),g1,wis[0]);
	else                                        loadtree((datadir/base1).string(), g1);
	if(base2.find("snr")==std::string::npos) load_gt_tree((datadir/base2).string(),g2,wis[1]);
	else                                        loadtree((datadir/base2).string(), g2);

	determine_mass_from_vol(bases[0],g1,wis[0]); // do this before scaling to mm!
	determine_mass_from_raw(bases[1],g2,wis[1]); // do this before scaling to mm!
	scale_to_mm(g1,wis[0]);
	scale_to_mm(g2,wis[1]);

	int n1 = num_vertices(g1);
	int n2 = num_vertices(g2);
	std::cout<<"# Verticex count: "<<n1 << " vs "<< n2<<std::endl;
	std::cout<<"# Tolerance     : "<<tolerance << std::endl;

	property_map<wurzelgraph_t,marked_vertex_t>::type   mark_map  = get(marked_vertex, g1);
	foreach (wurzel_vertex_descriptor& wd, vertices(g1))
		mark_map[wd] = true; // use all data in ground truth for radius estimation

	// estimate: avg in segments and convert to units/length
	determine_nodes_far_away_from_crossings(g2,true,true); // use only things far away from crossings for radius estimation, average radius in segments
	// ground truth: convert to units/length
	{
		auto pos_map  = get(vertex_position, g1);
		auto mass_map = get(edge_mass, g1);
		foreach(const wurzel_edge_descriptor &e, edges(g1)){
			const wurzel_vertex_descriptor& s = source(e,g1);
			const wurzel_vertex_descriptor& t = target(e,g1);
			double len = ublas::norm_2(pos_map[s]-pos_map[t]);
			mass_map[e] /= len; // converts mass to units per length(!)
		}
	}

	{
		// assume g1 is ground truth!
		std::ofstream  ofs( (datadir/base2/"missing.txt").string().c_str() );
		std::ofstream  scat( (datadir/base2/"missing_scatter.txt").string().c_str() );
		distance("missing",ofs,scat, g1,g2,tolerance, vm["verbose"].as<bool>());
	}
	{
		// assume g1 is ground truth!
		std::ofstream  ofs( (datadir/base2/"toomuch.txt").string().c_str() );
		std::ofstream  scat( (datadir/base2/"toomuch_scatter.txt").string().c_str() );
		distance("toomuch",ofs,scat, g2,g1,tolerance, vm["verbose"].as<bool>());
	}
}

/**
 * print all edges to a file for viewing
 */
template<class T>
void print_wurzel_edges(const std::string& name, wurzelgraph_t& wg, T& vidx_map, const wurzel_info& wi){
	std::ofstream ofs(name.c_str());
	property_map<wurzelgraph_t,root_stddev_t>::type stddev_map = get(root_stddev, wg);
	property_map<wurzelgraph_t,edge_mass_t>::type mass_map = get(edge_mass, wg);
	property_map<wurzelgraph_t,vertex_radius_t>::type radius_map = get(vertex_radius, wg);
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	foreach (const wurzel_edge_descriptor& e, edges(wg)){
		const wurzel_vertex_descriptor& s = source(e,wg);
		const wurzel_vertex_descriptor& t = target(e,wg);
		double length = voxdist(pos_map[s].begin(),pos_map[t].begin()); // assumed to be in mm here!!

		const vec3_t &spos = pos_map[s];
		const vec3_t &tpos = pos_map[t];
		const vec3_t  diff = tpos-spos;
		double        ang  = atan2(fabs(diff[0]),sqrt(diff[1] *diff[1]+diff[2]*diff[2]));

		unsigned int  v = vidx_map[s];
		unsigned int  w = vidx_map[t];
		//ofs << v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<w[0]<<" "<<w[1]<<" "<<w[2]<<std::endl;
		//double thickness = stddev_map[source(e,wg)]+stddev_map[target(e,wg)];
		double diameter = 0.5 * (2*radius_map[s]+2*radius_map[t]);
		double mass      = mass_map[e];
		ofs << v <<" "<< w <<" "<<diameter<<" "<<mass<< " "<<ang <<" "<<0.5*(spos[0]+tpos[0])<<std::endl;
	}
	ofs.close();
}

/**
 * a visitor which keeps a median statistic and stops when fixed amount of data
 * collected
 */
template<class WidthMap>
struct bfs_median_visitor:public default_bfs_visitor {
  bfs_median_visitor(WidthMap wmap, int num, medstat_t* m):m_stats(m),m_map(wmap), m_num(num){ }
  template < typename Vertex, typename Graph >
  void discover_vertex(Vertex u, const Graph & g) 
  {
	  (*m_stats)( m_map[u]);
	  if(count(*m_stats)==m_num){
	          throw median(*m_stats);
	  }
  }
  medstat_t* m_stats;
  WidthMap m_map;
  unsigned int      m_num;
};
template<class WidthMap>
bfs_median_visitor<WidthMap> make_bfs_median_visitor(WidthMap wm, int i, medstat_t* m){
	return bfs_median_visitor<WidthMap>(wm, i, m);
}

/**
 * print vertices to a file and add additional information for viewing
 *
 * uses the median filter above to smooth the radius
 */
template<class T>
void print_wurzel_vertices(const std::string& name, wurzelgraph_t& wg, T& vidx_map){
	std::ofstream ofs(name.c_str());
	unsigned int idx=0;
	stat_t s_diameter, s_mass;
	property_map<wurzelgraph_t,vertex_index_t>::type index_map = get(vertex_index, wg);
	property_map<wurzelgraph_t,root_stddev_t>::type stddev_map = get(root_stddev, wg);
	property_map<wurzelgraph_t,edge_mass_t>::type mass_map = get(edge_mass, wg);
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	property_map<wurzelgraph_t,vertex_eigenval_t>::type ev_map = get(vertex_eigenval, wg);
	property_map<wurzelgraph_t,vertex_radius_t>::type radius_map = get(vertex_radius, wg);
	foreach (wurzel_vertex_descriptor& wd, vertices(wg)){
		voxel_vertex_descriptor  v = get(vertex_name,wg)[wd];
		const vec3_t&            p = pos_map[wd];
		vidx_map[wd] = idx++;
		double mass=0.0;
		double mass_radius = 0.0;
		if(in_degree(wd,wg)>0){
			wurzel_edge_descriptor e = *in_edges(wd,wg).first;

			vec3_t s = pos_map[source(e,wg)];
			vec3_t t = pos_map[target(e,wg)];

			mass_radius = sqrt(mass_map[e] / M_PI); // has been divided by length already

			mass        = mass_map[e];
		}else if(out_degree(wd,wg)>0){
			wurzel_edge_descriptor e = *out_edges(wd,wg).first;

			vec3_t s = pos_map[source(e,wg)];
			vec3_t t = pos_map[target(e,wg)];

			mass_radius = sqrt(mass_map[e] / M_PI); // has been divided by length already

			mass        = mass_map[e];
		}
		s_mass(mass);

		//double d1 = pow(ev_map[wd][1], 2.00);
		//double d2 = pow(ev_map[wd][2], 2.00);

		double d1 = stddev_map[wd];
		double d2 = stddev_map[wd];
		medstat_t s_median;
		std::map<wurzel_vertex_descriptor,default_color_type> vertex2color;
		boost::associative_property_map< std::map<wurzel_vertex_descriptor, default_color_type> >
			    cmap(vertex2color);
		int smooth = 13;
		try{ breadth_first_visit(wg, wd,
				       	visitor(make_bfs_median_visitor(stddev_map,   smooth, &s_median)).
				       	color_map(cmap));
		}catch(double m){}
		try{ breadth_first_visit(make_reverse_graph(wg), wd,
				       	visitor(make_bfs_median_visitor(stddev_map,  2*smooth, &s_median)).
				       	color_map(cmap));
		}catch(double m){}
		//d1 = d2 = 2.0 * median(s_median);

		d1 = d2 = radius_map[wd] * 2;

		//double d1 = pow(ev_map[wd][1], 0.25);
		//double d2 = pow(ev_map[wd][2], 0.25);
		//std::cout << "d1: "<<d1<< " d2: "<< d2<<std::endl;
		double diameter = 0.5 * (d1+d2);

		ofs << p[0]<<" "<<p[1]<<" "<<p[2]<<" "<<mass_radius<<" "<<diameter<<" "<<0<<" "<<0<<std::endl;
		g_diameter[wd]=diameter;
		s_diameter(diameter);
	}
	ofs.close();
	std::cout << "diameter_avg_mm= "<<mean(s_diameter)<<std::endl;
	std::cout << "diameter_max_mm= "<<max(s_diameter)<<std::endl;
	std::cout << "diameter_min_mm= "<<min(s_diameter)<<std::endl;
	std::cout << "mass_avg_mg= "<<mean(s_mass)<<std::endl;
	std::cout << "mass_max_mg= "<<max(s_mass)<<std::endl;
	std::cout << "mass_min_mg= "<<min(s_mass)<<std::endl;
}

/**
 * print some interesting information about a graph
 */
void action_print(std::vector<wurzel_info>& wis, std::vector<std::string>& shortbases, const po::variables_map& vm){
  std::map<wurzel_vertex_descriptor,unsigned int> idx_map;

  std::string base = (fs::path(wis[0].directory) / shortbases[0]).string();

  wurzelgraph_t g;
  bool is_ground_truth = base.find("snr")==std::string::npos && base.find("Gerste")==std::string::npos;
  if(is_ground_truth)
	  load_gt_tree(base,g,wis[0]);
  else
	  loadtree(base, g);

  if (is_ground_truth)
	  determine_mass_from_vol(shortbases[0],g,wis[0]);
  else
	  determine_mass_from_raw(shortbases[0],g,wis[0]);

  scale_to_mm(g,wis[0]);
  if(!is_ground_truth){
	  determine_nodes_far_away_from_crossings(g,true,true); // average /radii/ in segments
  }
  else{
	// ground truth: convert to units/length
	  auto pos_map  = get(vertex_position, g);
	  auto mass_map = get(edge_mass, g);
	  foreach(const wurzel_edge_descriptor &e, edges(g)){
		  const wurzel_vertex_descriptor& s = source(e,g);
		  const wurzel_vertex_descriptor& t = target(e,g);
		  double len = ublas::norm_2(pos_map[s]-pos_map[t]);
		  mass_map[e] /= len; // converts mass to units per length(!)
	  }
  }
  total_mass(g, wis[0]);
  total_length(g);
  avg_len_btw_splits(g);
  print_wurzel_vertices(getfn(base,"vertices","txt"),g,idx_map);
  print_wurzel_edges(   getfn(base,"edges","txt"),g,idx_map,wis[0]);
  to_xml_forest(base,g);
}

int main(int argc, char* argv[]){
	std::vector<wurzel_info> wis;

	// read commandline params
	po::variables_map vm = get_config(wis,argc,argv);
	std::vector<std::string> actions = vm["action"].as<std::vector<std::string> >();
	std::vector<std::string> bases        = vm["base"].as<std::vector<std::string> >();
	std::vector<std::string> shortbases   = vm["base"].as<std::vector<std::string> >();

	for(unsigned int i=0;i<wis.size(); i++){
		bases[i] = (fs::path(wis[i].directory) / bases[i]).string();
	}

	//graph_stats(g1);
	//graph_stats(g2);
	
	// select action
	foreach(const std::string& a, actions){
		if(a=="distance"){
			action_distance(wis,shortbases,vm);
		}else if(a == "mass"){
			action_mass(wis,bases,vm);
		}else if(a == "length"){
			action_length(wis,bases,vm);
		}else if(a == "print"){
			action_print(wis,shortbases,vm);
		}
	}
	
}
