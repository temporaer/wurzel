#include <fstream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/foreach.hpp>
#include "icp.hpp"
#define foreach BOOST_FOREACH

namespace ublas=boost::numeric::ublas;
namespace acc=boost::accumulators;
typedef ublas::bounded_vector<float,3> point_type;
typedef util::ICP<3,point_type>        icp_type;
typedef acc::accumulator_set< double, acc::features< acc::tag::min, acc::tag::mean, acc::tag::max > > stat_t;

typedef ublas::bounded_matrix<double,3,3,ublas::column_major> mat_t;

//#include "detect_circle.hpp"

struct tensor_visitor{
	tensor_visitor(const point_type& _c){
		m = ublas::scalar_matrix<double>(3,3,0.0);
		sum = 0;
		c = _c;
	}
	mat_t m;
	point_type c;
	double sum;
	void operator()(const point_type& q){
		point_type p = q-c;
		m(0,0) += p[0] * p[0];
		m(1,1) += p[1] * p[1];
		m(2,2) += p[2] * p[2];
		m(0,1) += p[0] * p[1];
		m(1,2) += p[1] * p[2];
		sum += 1;
	}
	void finish(){
		m /= sum;
	}
	float plateness(){
		ublas::bounded_vector<double,3> lambda;
		//m(0,0) += 0.01;
		//m(1,1) += 0.01;
		//m(2,2) += 0.01;
		int info = boost::numeric::bindings::lapack::syev( 'N', 'U', m, lambda, boost::numeric::bindings::lapack::minimal_workspace() );
		if(info!=0){ 
			return 0.f; // could be point, line or plane... 
		}
		std::sort(lambda.begin(), lambda.end());
		return exp(-lambda[2]/lambda[1]) * exp(-fabs(lambda[0]));
	}
};


// inside-cylinder-test
// taken from
// http://www.flipcode.com/archives/Fast_Point-In-Cylinder_Test.shtml
float CylTest_CapsFirst( const point_type & pt1, const point_type & pt2, float lengthsq, float radius_sq, const point_type & testpt )
{
	float dx, dy, dz;	// vector d  from line segment point 1 to point 2
	float pdx, pdy, pdz;	// vector pd from point 1 to test point
	float dot, dsq;

	dx = pt2[0] - pt1[0];	// translate so pt1 is origin.  Make vector from
	dy = pt2[1] - pt1[1];     // pt1 to pt2.  Need for this is easily eliminated
	dz = pt2[2] - pt1[2];

	pdx = testpt[0] - pt1[0];		// vector from pt1 to test point.
	pdy = testpt[1] - pt1[1];
	pdz = testpt[2] - pt1[2];

	// Dot the d and pd vectors to see if point lies behind the 
	// cylinder cap at pt1[0], pt1[1], pt1[2]

	dot = pdx * dx + pdy * dy + pdz * dz;

	// If dot is less than zero the point is behind the pt1 cap.
	// If greater than the cylinder axis line segment length squared
	// then the point is outside the other end cap at pt2.

	if( dot < 0.0f || dot > lengthsq )
	{
		return( -1.0f );
	}
	else 
	{
		// Point lies within the parallel caps, so find
		// distance squared from point to line, using the fact that sin^2 + cos^2 = 1
		// the dot = cos() * |d||pd|, and cross*cross = sin^2 * |d|^2 * |pd|^2
		// Carefull: '*' means mult for scalars and dotproduct for vectors
		// In short, where dist is pt distance to cyl axis: 
		// dist = sin( pd to d ) * |pd|
		// distsq = dsq = (1 - cos^2( pd to d)) * |pd|^2
		// dsq = ( 1 - (pd * d)^2 / (|pd|^2 * |d|^2) ) * |pd|^2
		// dsq = pd * pd - dot * dot / lengthsq
		//  where lengthsq is d*d or |d|^2 that is passed into this function 

		// distance squared to the cylinder axis:

		dsq = (pdx*pdx + pdy*pdy + pdz*pdz) - dot*dot/lengthsq;

		if( dsq > radius_sq )
		{
			return( -1.0f );
		}
		else
		{
			return( dsq );		// return distance squared to axis
		}
	}
}

std::ostream&
operator<<(std::ostream& o, const stat_t& s) {
        o<<"min: "<<acc::min(s)<<"  mean: "<<acc::mean(s)<<"  max:"<<acc::max(s);
  return o;
}

void
read_points(const std::string& fn, std::vector<point_type>& rec, unsigned int ss=1){
	std::ifstream is(fn.c_str());
	float a,b,c;
	unsigned int cnt = 0;
	while(is.good()){
		is >> a>>b>>c;
		if((cnt++%ss)!=0) continue;
		point_type pt;
		pt[0] = a;
		pt[1] = b;
		pt[2] = c;
		rec.push_back(pt);
	}
}

point_type g_orig_pt1, g_orig_pt2;
double g_orig_radius, g_orig_len;
void
get_tube(const icp_type::TQuat& q, const icp_type::TVec& trans, const icp_type::TVec& scale, std::vector<point_type>& rec){
	const float res = 200;
	rec.reserve(res*res);
	//icp_type::VectorRotator rot;
	float heightfact = 1.0f;
	g_orig_radius = scale[0]*0.35f;
	g_orig_len    = scale[2]*heightfact;
       	for (float i = 0; i < 1.f; i+=1.f/(res))
	{
		for(float ang=0.f;ang<=2*M_PI;ang+=2*M_PI/res){
			point_type pt;
			pt[0] = sin(ang) * g_orig_radius;
			pt[1] = cos(ang) * g_orig_radius;
			pt[2] = i        * g_orig_len;
			//pt = rot(pt, q) + trans;
			pt = pt + trans;
			rec.push_back(pt);
		}
	}

	g_orig_pt1[0] = 0;
	g_orig_pt1[1] = 0;
	g_orig_pt1[2] = 0;
	g_orig_pt1   += trans;

	g_orig_pt2[0] = 0;
	g_orig_pt2[1] = 0;
	g_orig_pt2[2] = scale[2]* heightfact;
	g_orig_pt2   += trans;

	rec.push_back(g_orig_pt1);
	rec.push_back(g_orig_pt2);
}

int main(){
	std::vector<point_type> rec;
	std::vector<point_type> tube;
	rec.reserve(7000000);
	read_points("../data/reispflanze_wurzeln.txt", rec, 8);
	stat_t dim0,dim1,dim2;
	foreach(const point_type& p, rec){
		dim0(p[0]);
		dim1(p[1]);
		dim2(p[2]);
	}

	float mindim2 = acc::min(dim2)+0.1*(acc::max(dim2)-acc::min(dim2));
	std::cout << "stats0: "<<dim0<<std::endl;
	std::cout << "stats1: "<<dim1<<std::endl;
	std::cout << "stats2: "<<dim2<<std::endl;
	icp_type::TVec trans;
	trans[0] = acc::mean(dim0);
	trans[1] = acc::mean(dim1);
	trans[2] = acc::min(dim2);
	icp_type::TVec scale;
	scale[0] = acc::max(dim0)-acc::min(dim0);
	scale[1] = acc::max(dim1)-acc::min(dim1);
	scale[0] = scale[1] = std::min(scale[0], scale[1]);
	scale[2] = 0.5f*acc::max(dim2)-acc::min(dim2);

	std::vector<point_type> rec2;
	rec2.reserve(rec.size());

	foreach(const point_type& p, rec){ // permute 
		if(p[2]<mindim2)
			continue;
		rec2.push_back(p);
		//do{
		//        rec2.back()[2] = rec[drand48()*rec.size()][2];
		//}while(rec2.back()[2]<mindim2);
	}

	std::cout << "Registering..."<<std::endl;
	icp_type icp(1000);
	icp.setVerbose(true);
	icp.registerModel(rec2.begin(), rec2.end());
	icp_type::TVec modelCentroid = icp.getModelCentroid();

	std::cout << "Get Tube..."<<std::endl;
	get_tube(icp_type::TQuat(1,0,0,0),trans,scale,tube);

	std::cout << "Matching..."<<std::endl;
	icp.match(tube.begin(), tube.end());
	std::cout << "Trans: "<< icp.getTrans()<<std::endl;
	std::cout << "Rot:   "<< icp.getRot()<<std::endl;
	std::cout << "Scale: "<< icp.getScale()<<std::endl;

	// translate rec to centered coordinates
	foreach(point_type& p, rec){
		p[0]-=modelCentroid[0];
		p[1]-=modelCentroid[1];
		p[2]-=modelCentroid[2];
	}

	// transform the cylinder according to ICP outcome
	point_type tube_centroid; int n;
	boost::tie(tube_centroid, n) = icp.getCentroid(tube.begin(),tube.end());
	point_type pt1  = icp_type::VectorRotator()(g_orig_pt1-tube_centroid,icp.getRot())*icp.getScale() + icp.getTrans();
	point_type pt2  = icp_type::VectorRotator()(g_orig_pt2-tube_centroid,icp.getRot())*icp.getScale() + icp.getTrans();
	pt1 = tube[tube.size()-2];
	pt2 = tube[tube.size()-1];

	double       r2 = g_orig_radius * icp.getScale(); r2 *= r2;
	double       l2 = g_orig_len    * icp.getScale(); l2 *= l2;

	std::cout << "Pt1: "<<pt1<<std::endl;
	std::cout << "Pt2: "<<pt2<<std::endl;
	std::cout << "R2 : "<<r2<<std::endl;
	std::cout << "L2 : "<<l2<<std::endl;

	// output
	unsigned int cnt=0;
	std::ofstream osrec("rec.txt");
	std::ofstream osrec_del("rec_del.txt");
	std::ofstream osrec_out("rec_out.txt");
	stat_t s_plateness;
	unsigned int progress=0;
	unsigned int total = rec.size();
	foreach(const point_type& p, rec){
		progress++;
		if(progress%1000==0){
			std::cout << "\r";
			int i=0;
			for(i=0;i<progress; i+=total/100){
				std::cout<<"#";
			}
			for(;i<total; i+=total/100){
				std::cout<<"-";
			}
			std::cout<<std::flush;
		}

		float d2 = CylTest_CapsFirst(pt1,pt2,l2,r2,p);
		if(d2<0) {
			osrec_out << p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
			continue;
		}
		if((r2-d2)/r2 < 0.1){
		        tensor_visitor tv = icp.modelTree().visit_within_range(p, 8, tensor_visitor(p));
		        tv.finish();
		        float pn = tv.plateness();
		        s_plateness(pn);
		        if(pn > 0.10){
		                osrec_del << p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
		                continue;
		        }
		}
		osrec << p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
		cnt++;
	}
	std::cout<<"cnt: "<< cnt<<std::endl;
	std::cout<<"pn : "<< s_plateness<<std::endl;
	std::ofstream ostube("tube.txt");
	foreach(const point_type& p, tube){
		ostube << p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
	}
}
