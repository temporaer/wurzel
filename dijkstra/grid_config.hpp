#ifndef __GRID_CONFIG_HPP__
#define __GRID_CONFIG_HPP__

#include <iostream>
#include <boost/program_options.hpp>
#include "wurzel_info.hpp"

namespace po = boost::program_options;

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
			("force", "force recomputation of dijkstra algorithm")
			("stem-plane", po::value<int>(),"plane index in which to search for stem")
			("stem-axis",  po::value<int>(),"the axis of stem-plane")
			("start-threshold,s", po::value<double>()->default_value(0.1),    "minimum raw value to start tracking")
			("total-len-frac,t", po::value<double>()->default_value(1.0),     "maximal total length fraction")
			("avg-len-frac,a", po::value<double>()->default_value(0.20),      "maximal average length fraction")
			("min-flow-thresh,f", po::value<double>()->default_value(0.0001), "minimal flow fraction")
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
	if(base.find("Gerste")>=0){
		std::cout <<"Detected instance of Gerste"<<std::endl;
		wi.X = 872;
		wi.Y = 192;
		wi.Z = 192;
		if(vm.count("stem-plane")==0)
			wi.stem_plane = 5;
		if(vm.count("stem-plane")==0)
			wi.stem_axis  = 0;
	}
	else if(base.find("L2")>=0){
		std::cout <<"Detected instance of Lupine"<<std::endl;
		wi.X = 256;
		wi.Y = 256;
		wi.Z = 256;
		if(vm.count("stem-plane")==0)
			wi.stem_plane = 24;
		if(vm.count("stem-plane")==0)
			wi.stem_axis  = 2;
	}
	if(vm.count("stem-plane"))
		wi.stem_plane = vm["stem-plane"].as<int>();
	if(vm.count("stem-axis"))
		wi.stem_axis = vm["stem-axis"].as<int>();

	wi.XYZ=wi.X*wi.Y*wi.Z;

	return vm;
}

#endif /* __GRID_CONFIG_HPP__ */
