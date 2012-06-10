//=======================================================================
// Copyright 2011 University Bonn
// Author: Hannes Schulz
//=======================================================================
#ifndef __GRID_CONFIG_HPP__
#define __GRID_CONFIG_HPP__

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <inc/wurzel_tree.hpp>
#include "wurzel_info.hpp"
#include "../config.hxx"
#define foreach BOOST_FOREACH

namespace po = boost::program_options;

/// method for selecting leaf candidates
enum POR_METHOD{ 
    /// a leaf candidate has 
    /// - raw value > start_thresh
    /// - total path len < total_len_thresh
    /// - median(downstream) / median(upstream) > flow_thresh
    /// This is the method described in the VISAPP Paper
    POR_EDGE_DETECT, 

    /// a leaf candidate has 
    /// - raw value > start_thresh
    /// - total path len < total_len_thresh
    /// - median(downstream,downstream) > flow_thresh
    POR_MEDIAN_RAW, 

    /// a leaf candidate has 
    /// - raw value > start_thresh
    /// - total path len < total_len_thresh
    /// - the subtree has at least weight flow_thresh
    /// This is the method originally used for Barley
    POR_SUBTREE_WEIGHT
};

po::variables_map
get_config(wurzel_info& wi, int argc, char* argv[]){
	po::variables_map vm;
	try{
		po::options_description desc("Allowed options");
		po::positional_options_description p;
		p.add("base", 1);
		desc.add_options()
			("help",  "produce help message")
			("base",  "the base name of the dataset")
			("cfg,c", po::value<std::string>()->default_value("config.xml"),"the config-file containing dataset descriptions (XML)")
			("force", "force recomputation of dijkstra algorithm")
			("stem-plane", po::value<int>(),"plane index in which to search for stem")
			("stem-axis",  po::value<int>(),"the axis of stem-plane")
			("max-radius,r", po::value<double>()->default_value(1.8),    "maximum root radius (in mm)")
			("start-threshold,s", po::value<double>()->default_value(0.1),    "minimum raw value to start tracking")
			("total-len-thresh,t", po::value<double>()->default_value(1E9),   "maximal total length")
			("dijkstra-stop-val,d", po::value<double>()->default_value(1E9),   "stop dijkstra when paths longer than this (decreases dijkstr runtime)")
			("avg-len-frac,a", po::value<double>()->default_value(0.20),      "maximal average length fraction")
            ("leaf-select-method,l", po::value<std::string>()->default_value("median_raw"), "method for selecting leaf candidates [edge_detect,median_raw,subtree_weight]")
			("min-flow-thresh,f", po::value<double>()->default_value(0.0001), "minimal flow fraction")
			("no-gauss-fit",  "save some time")
			("no-subpix-pos", "save some time")
			;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);
		if(vm.count("help")){
			std::cout << desc << std::endl;
			exit(0);
		}
	}catch(std::exception& e){
		std::cerr << "error: "<<e.what() << std::endl;
		exit(1);
	}catch(...){
		std::cerr << "unknown error"<<std::endl;
		exit(1);
	}
	std::string base = vm["base"].as<std::string>();
	std::string config_file = vm["cfg"].as<std::string>();
	std::cout << "Reading config from file `"<<config_file<<"'."<<std::endl;

	std::ifstream config_ifs(config_file.c_str());
	std::auto_ptr<config_t> cfg;
	try{
		cfg = config(config_ifs);
	}catch(xml_schema::exception& e){
		std::cerr << e << std::endl;
		exit(1);
	}catch (const xml_schema::properties::argument&) {
		std::cerr << "invalid property argument (empty namespace or location)" << std::endl;
		exit(1);
	}
	catch (const xsd::cxx::xml::invalid_utf16_string&) {
		std::cerr << "invalid UTF-16 text in DOM model" << std::endl;
		exit(1);
	}
	catch (const xsd::cxx::xml::invalid_utf8_string&) {
		std::cerr << "invalid UTF-8 text in object model" << std::endl;
		exit(1);
	}

	int found = 0;
	foreach(datafile_t& df,cfg->datafiles().datafile()){
		if(df.base_name()==base){
			wi.X = df.shape().x();
			wi.Y = df.shape().y();
			wi.Z = df.shape().z();
			wi.read_XYZ = df.read_shape().x() 
				*     df.read_shape().y()
				*     df.read_shape().z() ;
			wi.stem_plane = df.stem_plane();
			wi.stem_axis  = df.stem_axis();
			wi.scale      = df.scale();
			wi.noise_cutoff = df.noise_cutoff();
			wi.directory  = cfg->datapath();
			found++;
		}
	}
	if(!found){
		std::cerr << "could not find in config-file: `"<<base<<"'."<<std::endl;
		exit(1);
	}
	if(found>1){
		std::cerr << "could more than once in config-file: `"<<base<<"'."<<std::endl;
		exit(1);
	}

	if(vm.count("stem-plane"))
		wi.stem_plane = vm["stem-plane"].as<int>();
	if(vm.count("stem-axis"))
		wi.stem_axis  = vm["stem-axis"].as<int>();

	wi.XYZ=wi.X*wi.Y*wi.Z;

	return vm;
}

#endif /* __GRID_CONFIG_HPP__ */
