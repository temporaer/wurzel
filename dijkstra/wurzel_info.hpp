#ifndef __WURZEL_INFO_HPP__
#define __WURZEL_INFO_HPP__
#include <boost/array.hpp>

struct wurzel_info
{
	unsigned int X,Y,Z;   ///< dimensions of data cube
	unsigned int XYZ;     ///< == X*Y*Z
	unsigned int stem_plane; ///< plane in which to search for stem
	unsigned int stem_axis;  ///< axis of stem_plane
	double       scale;   ///< multiply with this to get milli meter
	double       spross_intensity;   ///< normalize data by this
	double       noise_cutoff; ///< signal below this is "noise"
	boost::array<vidx_t, 3> strunk;
};

#endif /* __WURZEL_INFO_HPP__ */
