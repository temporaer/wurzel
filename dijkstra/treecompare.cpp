#include <iostream>
#include <fstream>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include <boost/numeric/ublas/io.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "voxel_accessors.hpp"
#include "wurzel_tree.hpp"
#include "treecompare_config.hpp"

using namespace boost::accumulators;
typedef accumulator_set< double, features< tag::min, tag::mean, tag::max, tag::variance > > stat_t;


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
	double minpos = wi.scale * /*wi.stem_plane*/15 * wi.Z/410; // TODO: make this work for other datasets!!
	property_map<wurzelgraph_t,vertex_position_t>::type  pos_map = get(vertex_position, wg);
	property_map<wurzelgraph_t,root_stddev_t>::type   stddev_map = get(root_stddev, wg);
	property_map<wurzelgraph_t,edge_mass_t>::type       mass_map = get(edge_mass, wg);
	property_map<wurzelgraph_t,vertex_param0_t>::type param0_map = get(vertex_param0, wg);
	foreach(const wurzel_edge_descriptor& e, edges(wg)){
		const wurzel_vertex_descriptor& s = source(e,wg);
		const wurzel_vertex_descriptor& t = target(e,wg);
		if(pos_map[s][wi.stem_axis] < minpos)
			continue;                                      
		if(pos_map[t][wi.stem_axis] < minpos)
			continue;
		double length  = voxdist(pos_map[s].begin(),pos_map[t].begin());                   // assume scaled to mm in advance

		double radiuss = 1.00 *stddev_map[s];                                              // already in mm
		double radiust = 1.00 *stddev_map[t];                                              // already in mm

		double param1s = 1.0/2.0/stddev_map[s]/stddev_map[s];                              // determine param in exp of gauss again
		double param1t = 1.0/2.0/stddev_map[t]/stddev_map[t];                              // determine param in exp of gauss again

		double masss   = -M_PI *param0_map[s] *(exp(-param1s *radiuss*radiuss)-1)/param1s;
		double masst   = -M_PI *param0_map[s] *(exp(-param1t *radiust*radiust)-1)/param1t;

		mass_map[e]    = length *masss;
		//mass_map[e] = length / 3.0 * (masss +masst +sqrt(masss *masst));
		sum += mass_map[e];
	}
	const double weight_scale     = 0.2712; // milli-gram
	return sum*weight_scale;
}
double total_length(wurzelgraph_t& wg){
	double sum = 0.0;
	property_map<wurzelgraph_t,vertex_position_t>::type pos_map  = get(vertex_position, wg);
	foreach(const wurzel_edge_descriptor& e, edges(wg)){
		const wurzel_vertex_descriptor& s = source(e,wg);
		const wurzel_vertex_descriptor& t = target(e,wg);
		sum += voxdist(pos_map[s].begin(),pos_map[t].begin()); // assumed to be in mm here!!
	}
	return sum;
}
void distance(wurzelgraph_t& g1, wurzelgraph_t& g2, const double& tolerance){
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
	std::cout <<"Large against small: "<<std::endl;
	if(n1>n2)
		distance(g1,g2,tolerance);
	else
		distance(g2,g1,tolerance);
	std::cout <<"Small against large: "<<std::endl;
	if(n1>n2)
		distance(g2,g1,tolerance);
	else
		distance(g1,g2,tolerance);
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
		}
	}
	
}