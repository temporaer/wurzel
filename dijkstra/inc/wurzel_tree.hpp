//=======================================================================
// Copyright 2011 University Bonn
// Author: Hannes Schulz
//=======================================================================
#ifndef __WURZEL_TREE_HPP__
#define __WURZEL_TREE_HPP__
#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include "voxelgrid.hpp"

namespace ublas = boost::numeric::ublas;

struct vertex_normal_t{
	typedef boost::vertex_property_tag kind;
} vertex_normal;

struct vertex_position_t{
	typedef boost::vertex_property_tag kind;
} vertex_position;

struct vertex_eigenval_t{
	typedef boost::vertex_property_tag kind;
} vertex_eigenval;

struct vertex_radius_t{
	typedef boost::vertex_property_tag kind;
} vertex_radius;

struct marked_vertex_t{
	typedef boost::vertex_property_tag kind;
} marked_vertex;

struct root_stddev_t{
	typedef boost::vertex_property_tag kind;
} root_stddev;

struct vertex_param0_t{
	typedef boost::vertex_property_tag kind;
} vertex_param0;

struct edge_mass_t{
	typedef boost::edge_property_tag kind;
} edge_mass;

typedef ublas::bounded_matrix<double,3,3,ublas::column_major> covmat_t;
typedef ublas::bounded_vector<double,3>                       vec3_t;

typedef boost::adjacency_list<
	boost::listS,
	boost::listS,
	boost::bidirectionalS,
	boost::property<boost::vertex_name_t,voxel_vertex_descriptor,
	  boost::property<vertex_position_t,vec3_t,
		  boost::property<vertex_normal_t,covmat_t,
		        boost::property<vertex_eigenval_t,vec3_t,
		  	     boost::property<root_stddev_t,double,
		  	          boost::property<vertex_param0_t,double,
			               boost::property<marked_vertex_t,bool,
					       boost::property<vertex_radius_t,double
	  > > > > > > > >,
	boost::property<edge_mass_t,double>
	> wurzelgraph_t;

typedef graph_traits<wurzelgraph_t> wurzelg_traits;
typedef wurzelg_traits::edge_descriptor wurzel_edge_descriptor;
typedef wurzelg_traits::edge_iterator   wurzel_edge_iterator;
typedef wurzelg_traits::in_edge_iterator   wurzel_in_edge_iterator;
typedef wurzelg_traits::vertex_descriptor wurzel_vertex_descriptor;
typedef wurzelg_traits::vertex_iterator   wurzel_vertex_iterator;


#endif
