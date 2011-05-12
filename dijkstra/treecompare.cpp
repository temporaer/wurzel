#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/unordered_map.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/reverse_graph.hpp>

#include <boost/numeric/ublas/io.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "voxel_accessors.hpp"
#include "wurzel_tree.hpp"
#include "treecompare_config.hpp"

using namespace boost::accumulators;
typedef accumulator_set< double, features< tag::min, tag::mean, tag::max, tag::variance, tag::count > > stat_t;
typedef accumulator_set< double, features< tag::median, tag::count > > medstat_t;

std::map<void*,double> g_radius;


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


voxel_edge_weight_map::reference 
voxel_edge_weight_map::operator[](key_type e) const {
	voxel_vertex_descriptor s = source(e,m_graph);
	voxel_vertex_descriptor t = target(e,m_graph);
	return 0.0;
}


void loadtree(std::string basename, wurzelgraph_t& g){
	std::ifstream ifs((basename+"-wgraph.ser").c_str());
	boost::archive::text_iarchive ia(ifs);
	ia >> g;
}

template<class T>
struct kv_helper{
  std::string k;
  T v;
};
template<class T>
kv_helper<T> kv(std::string s, T o){ kv_helper<T> x;x.k=s; x.v=o; return x;}
template<class T>
std::ostream& operator<<(std::ostream& o, kv_helper<T> kvh){
	o << " "<<kvh.k<<"=\""<<kvh.v<<"\" ";
	return o;
}

void to_xml_node(std::ofstream& o, unsigned int& idx, const wurzel_vertex_descriptor& v, const wurzelgraph_t& g){
	property_map<wurzelgraph_t,vertex_index_t>::const_type  index_map = get(vertex_index, g);
	property_map<wurzelgraph_t,vertex_position_t>::const_type  pos_map = get(vertex_position, g);
	property_map<wurzelgraph_t,root_stddev_t>::const_type   stddev_map = get(root_stddev, g);
	property_map<wurzelgraph_t,edge_mass_t>::const_type       mass_map = get(edge_mass, g);
	property_map<wurzelgraph_t,vertex_param0_t>::const_type param0_map = get(vertex_param0, g);

	o << "<Node "
	  << kv("id", idx)
	  << kv("bo", 1)
	  << kv("rad", g_radius[v]/200.0) // meters, i suppose...
	  << kv("x", pos_map[v][2]/100.0-0.5)
	  << kv("y", pos_map[v][1]/100.0-0.5)
	  << kv("z", pos_map[v][0]/100.0-0.5)
	  << ">"<<std::endl;
	foreach(const wurzel_vertex_descriptor& w, adjacent_vertices(v,g)){
		to_xml_node(o,++idx,w,g);
	}
	o << "</Node>"<<std::endl;
}

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
	std::cout << "xmin: "<< xmin<<std::endl;
	std::cout << "xmax: "<< xmax<<std::endl;
	std::cout << "ymin: "<< ymin<<std::endl;
	std::cout << "ymax: "<< ymax<<std::endl;
	std::cout << "zmin: "<< zmin<<std::endl;
	std::cout << "zmax: "<< zmax<<std::endl << std::endl;
}
double total_mass(wurzelgraph_t& wg, const wurzel_info& wi){
	double sum = 0.0;
	double minpos = wi.scale * /*wi.stem_plane*/14 * wi.X/410.0; // TODO: make this work for other datasets!!
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

		static const double k=1.2;
		const double param1s = 1.0/(2.0*stddev_map[s]*stddev_map[s]*wi.scale*wi.scale);                              // determine param in exp of gauss again
		const double param1t = 1.0/(2.0*stddev_map[t]*stddev_map[t]*wi.scale*wi.scale);                              // determine param in exp of gauss again

		static const double vox_to_mm3 = wi.scale*wi.scale*wi.scale;
		const double param0s = param0_map[s]/vox_to_mm3;
		const double param0t = param0_map[t]/vox_to_mm3;

		//double masss   = -M_PI *param0s *(exp(-param1s *stddevs*stddevs)-1)/param1s;
		//double masst   = -M_PI *param0t *(exp(-param1t *stddevt*stddevt)-1)/param1t;
		const double masss = M_PI*(2*param0s*(1-exp(-k*k/2.0)) + k*k*wi.noise_cutoff/4)/(2*param1s);
		const double masst = M_PI*(2*param0t*(1-exp(-k*k/2.0)) + k*k*wi.noise_cutoff/4)/(2*param1t);

		//mass_map[e]    = masss*length;
		mass_map[e] = length / 3.0 * (masss +masst +sqrt(masss *masst));
		sum += mass_map[e];

		const vec3_t& spos = pos_map[s];
		const vec3_t& tpos = pos_map[t];
		const vec3_t  diff = tpos-spos;
		bool is_vert = diff[0]*diff[0] > (diff[1]*diff[1]+diff[2]*diff[2]);
		if(is_vert) sum_vert  += mass_map[e];
		else        sum_horiz += mass_map[e];
	}
	double weight_scale     = 0.2712; // milli-gram
	weight_scale /= wi.X*wi.Y*wi.Z  / 192.0 / 192.0 / 410.0; // change 0.27 so that it fits the resolution
	std::cout << "Total mass: "<<sum * weight_scale <<std::endl;
	std::cout << "Horiz mass: "<<sum_horiz * weight_scale <<std::endl;
	std::cout << "Vertc mass: "<<sum_vert  * weight_scale <<std::endl;
	std::cout << "V/H   mass: "<<sum_vert/sum_horiz<<std::endl;
	return sum*weight_scale;
}

struct wurzel_segment{
	const wurzelgraph_t& wg;
	wurzel_segment(const wurzelgraph_t& g):wg(g){}
	std::list<wurzel_edge_descriptor> edges;
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
	
	
	std::cout << "Found #segments: "<<segments.size()<<std::endl;
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
	std::cout <<"Avg seg len: "<< mean(s_lens)<<" +/- "<<sqrt(variance(s_lens))<< " mm"<<std::endl;
	std::cout <<"        ver: "<< mean(s_lens_ver)<<" +/- "<<sqrt(variance(s_lens_ver))<< " mm"<<std::endl;
	std::cout <<"        hor: "<< mean(s_lens_hor)<<" +/- "<<sqrt(variance(s_lens_hor))<< " mm"<<std::endl;
	std::cout <<"Avg seg num: "<< mean(s_nums)<<" +/- "<<sqrt(variance(s_nums))<< " pieces"<<std::endl;
	return mean(s_lens);
}
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
	std::cout << "Total len: "<<sum <<std::endl;
	std::cout << "Horiz len: "<<sum_horiz <<std::endl;
	std::cout << "Vertc len: "<<sum_vert  <<std::endl;
	std::cout << "V/H   len: "<<sum_vert/sum_horiz<<std::endl;

	return sum;
}
void distance(wurzelgraph_t& g1, wurzelgraph_t& g2, const double& tolerance, bool verbose){
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map1  = get(vertex_position, g1);
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map2  = get(vertex_position, g2);
	int no_counterpart = 0;
	int v_cnt = 0;
	stat_t s_distances;

	foreach(const wurzel_vertex_descriptor& v1, vertices(g1)){
		const vec3_t& pv1 = pos_map1[v1];
		double min_d = 1E6;
		foreach(const wurzel_edge_descriptor& e2, edges(g2)){
			const vec3_t& pe2s = pos_map2[source(e2,g2)];
			const vec3_t& pe2t = pos_map2[target(e2,g2)];
			double dp = ublas::inner_prod(pe2s-pv1, pe2t-pe2s);
			double  t = dp / pow(ublas::norm_2(pe2t-pe2s),2);
			//if(t<-tolerance) continue;
			if(t<0){
				min_d = std::min(min_d, ublas::norm_2(pv1-pe2s));
				continue;
			}
			//if(t-1>tolerance) continue;
			else if(t>1){
				min_d = std::min(min_d, ublas::norm_2(pv1-pe2t));
				continue;
			}
			else{
				double d = ublas::norm_2(cross_product(pv1-pe2s,pv1-pe2t))/ublas::norm_2(pe2t-pe2s);
				//if(d<tolerance)
				min_d = std::min(min_d, d);
			}
		}
		if(min_d>tolerance)
			no_counterpart ++;
		s_distances(min_d);
		v_cnt ++;
		if(v_cnt % 1000==0 && verbose)
			std::cout <<"\rFraction not found: "<<((float)no_counterpart/v_cnt)<<" count: "<<v_cnt<< " mind: "<< min_d<<std::flush;
	}
	std::cout <<std::endl<< "Vertices w/o counterpart: "<<no_counterpart<<" at tolerance "<<tolerance<<std::endl;
	std::cout <<std::endl<< "Average distance: "<<mean(s_distances)<< " var:"<<variance(s_distances)<<std::endl;
	
}

void scale_to_mm(wurzelgraph_t& g, const wurzel_info& wi){
	std::cout << "scaling: "<< wi.scale <<std::endl;
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, g);
	foreach(wurzel_vertex_descriptor& v, vertices(g)){
		vec3_t& p = pos_map[v];
		p *= wi.scale;
	}
}

void action_length(std::vector<wurzel_info>& wis, std::vector<std::string>& bases, const po::variables_map& vm){
	std::string base = bases[0];

	wurzelgraph_t g;
	loadtree(base, g);

	scale_to_mm(g,wis[0]);

	std::cout << "total length("<<base<<"): "<< total_length(g)<<std::endl;

}
void action_mass(std::vector<wurzel_info>& wis, std::vector<std::string>& bases, const po::variables_map& vm){
	std::string base = bases[0];

	wurzelgraph_t g;
	loadtree(base, g);

	scale_to_mm(g,wis[0]);

	std::cout << "total mass ("<<base<<"): "<< total_mass(g, wis[0])<<std::endl;

}
void action_distance(std::vector<wurzel_info>& wis, std::vector<std::string>& bases, const po::variables_map& vm){
	std::string base1 = bases[0];
	std::string base2 = bases[1];
	double tolerance = vm["tolerance"].as<double>();

	wurzelgraph_t g1, g2;
	loadtree(base1, g1);
	loadtree(base2, g2);

	scale_to_mm(g1,wis[0]);
	scale_to_mm(g2,wis[1]);

	int n1 = num_vertices(g1);
	int n2 = num_vertices(g2);
	std::cout<<"Verticex count: "<<n1 << " vs "<< n2<<std::endl;
	distance(g1,g2,tolerance, vm["verbose"].as<bool>());
}

template<class T>
void print_wurzel_edges(const std::string& name, wurzelgraph_t& wg, T& vidx_map, const wurzel_info& wi){
	std::ofstream ofs(name.c_str());
	property_map<wurzelgraph_t,root_stddev_t>::type stddev_map = get(root_stddev, wg);
	property_map<wurzelgraph_t,edge_mass_t>::type mass_map = get(edge_mass, wg);
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
		double thickness = mass_map[e]/length;
		double mass      = mass_map[e]/length;
		ofs << v <<" "<< w <<" "<<thickness<<" "<<mass<< " "<<ang <<" "<<0.5*(spos[0]+tpos[0])<<std::endl;
	}
	ofs.close();
}

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

template<class T>
void print_wurzel_vertices(const std::string& name, wurzelgraph_t& wg, T& vidx_map){
	std::ofstream ofs(name.c_str());
	unsigned int idx=0;
	property_map<wurzelgraph_t,vertex_index_t>::type index_map = get(vertex_index, wg);
	property_map<wurzelgraph_t,root_stddev_t>::type stddev_map = get(root_stddev, wg);
	property_map<wurzelgraph_t,edge_mass_t>::type mass_map = get(edge_mass, wg);
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	property_map<wurzelgraph_t,vertex_eigenval_t>::type ev_map = get(vertex_eigenval, wg);
	foreach (wurzel_vertex_descriptor& wd, vertices(wg)){
		voxel_vertex_descriptor  v = get(vertex_name,wg)[wd];
		const vec3_t&            p = get(vertex_position,wg)[wd];
		vidx_map[wd] = idx++;
		//unsigned int deg = out_degree(wd,wg);
		double mass = 0, cnt = 0;
		if(in_degree(wd,wg)>0){
			wurzel_edge_descriptor e = *in_edges(wd,wg).first;
			double l = voxdist(pos_map[source(e,wg)].begin(),pos_map[target(e,wg)].begin());
			mass += mass_map[e]/l;
			cnt ++;
		}
		if(out_degree(wd,wg)>0){
			wurzel_edge_descriptor e = *out_edges(wd,wg).first;
			double l = voxdist(pos_map[source(e,wg)].begin(),pos_map[target(e,wg)].begin());
			mass += mass_map[*out_edges(wd,wg).first]/l;
			cnt ++;
		}
		if(cnt>0) mass /= cnt;
		mass = std::max(mass, 0.05); // just for visualization!
		mass = std::min(mass, 6.00); // just for visualization!

		//double d1 = pow(ev_map[wd][1], 2.00);
		//double d2 = pow(ev_map[wd][2], 2.00);

		double d1 = stddev_map[wd];
		double d2 = stddev_map[wd];
		medstat_t s_median;
		std::map<wurzel_vertex_descriptor,default_color_type> vertex2color;
		boost::associative_property_map< std::map<wurzel_vertex_descriptor, default_color_type> >
			    cmap(vertex2color);
		int smooth = 16;
		try{ breadth_first_visit(wg, wd,
				       	visitor(make_bfs_median_visitor(stddev_map,   smooth, &s_median)).
				       	color_map(cmap));
		}catch(double m){}
		try{ breadth_first_visit(make_reverse_graph(wg), wd,
				       	visitor(make_bfs_median_visitor(stddev_map,  2*smooth, &s_median)).
				       	color_map(cmap));
		}catch(double m){}
		d1 = d2 = 2.0 * median(s_median);

		//double d1 = pow(ev_map[wd][1], 0.25);
		//double d2 = pow(ev_map[wd][2], 0.25);
		//std::cout << "d1: "<<d1<< " d2: "<< d2<<std::endl;
		double radius = 0.5 * (d1+d2)    * 0.8;
		ofs << p[0]<<" "<<p[1]<<" "<<p[2]<<" "<<mass<<" "<<radius<<" "<<0<<" "<<0<<std::endl;
		g_radius[wd]=radius;
	}
	ofs.close();
}

void action_print(std::vector<wurzel_info>& wis, std::vector<std::string>& bases, const po::variables_map& vm){
  std::map<wurzel_vertex_descriptor,unsigned int> idx_map;

  wurzelgraph_t g;
  loadtree(bases[0], g);

  scale_to_mm(g,wis[0]);
  total_mass(g, wis[0]);
  total_length(g);
  avg_len_btw_splits(g);
  print_wurzel_vertices(getfn(bases[0],"vertices","txt"),g,idx_map);
  print_wurzel_edges(   getfn(bases[0],"edges","txt"),g,idx_map,wis[0]);
  to_xml_forest(bases[0],g);
}

int main(int argc, char* argv[]){
	std::vector<wurzel_info> wis;
	po::variables_map vm = get_config(wis,argc,argv);
	std::vector<std::string> actions = vm["action"].as<std::vector<std::string> >();
	std::vector<std::string> bases   = vm["base"].as<std::vector<std::string> >();

	//graph_stats(g1);
	//graph_stats(g2);
	
	foreach(const std::string& a, actions){
		if(a=="distance"){
			action_distance(wis,bases,vm);
		}else if(a == "mass"){
			action_mass(wis,bases,vm);
		}else if(a == "length"){
			action_length(wis,bases,vm);
		}else if(a == "print"){
			action_print(wis,bases,vm);
		}
	}
	
}
