#ifndef __VOXEL_NORMALS_HPP__
#define __VOXEL_NORMALS_HPP__


#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "wurzel_tree.hpp"
#include "voxel_accessors.hpp"

#define SQR(X) ((X)*(X))
template<class T>
void
get_normal(covmat_t&M,	vec3_t& l, const vec3_t& v, const T& acc){
	namespace ublas = boost::numeric::ublas;
	namespace lapack = boost::numeric::bindings::lapack;
	using namespace std;
	using ublas::prod;
	using ublas::trans;
	//typedef ublas::bounded_matrix<double,3,3,ublas::column_major> mat;
	ublas::bounded_vector<double,3> lambda;
	ublas::bounded_vector<double,34*3> work; // the 34*width is from syev.hpp
	covmat_t A;
	double x=v[0],y=v[1],z=v[2];
	double sum = 0;
	const int r = 3;
	for(int i=-r;i<=r;i++){
		for(int j=-r;j<=r;j++){
			for(int k=-r;k<=r;k++){
				double d = acc(x+i,y+j,z+k);

				// covariance matrix
				A(0,0) += d * SQR(i);
				A(1,1) += d * SQR(j);
				A(2,2) += d * SQR(k);

				A(0,1) += d * i * j;
				A(0,2) += d * i * k;
				A(1,2) += d * j * k;

				sum += d;
			}
		}
	}
	A /= sum;

	int info = lapack::syev('V', 'U', A, lambda, lapack::optimal_workspace());
	if(info!=0)
		return;

	unsigned int idx = distance(lambda.begin(),max_element(lambda.begin(),lambda.end()));
	if(idx!=0){
		std::swap(A(0,0),A(0,idx));
		std::swap(A(1,0),A(1,idx));
		std::swap(A(2,0),A(2,idx));
		std::swap(lambda(0),lambda(idx));
	}
	M = A;
	l = lambda;
}

template<class T>
void
get_inertia_tensor(covmat_t&M,	vec3_t& l, const vec3_t& v, const T& acc, const double& r, const double& step){
	namespace ublas = boost::numeric::ublas;
	namespace lapack = boost::numeric::bindings::lapack;
	using namespace std;
	using ublas::prod;
	using ublas::trans;
	//typedef ublas::bounded_matrix<double,3,3,ublas::column_major> mat;
	ublas::bounded_vector<double,3> lambda;
	ublas::bounded_vector<double,34*3> work; // the 34*width is from syev.hpp
	covmat_t A;
	double x=v[0],y=v[1],z=v[2];
	double sum = 0;
	for(double i=-r;i<=r;i+=step){
		for(double j=-r;j<=r;j+=step){
			for(double k=-r;k<=r;k+=step){
				double d = acc(x+i,y+j,z+k);
				// inertia tensor
				A(0,0) += d * (SQR(j)+SQR(k));
				A(1,1) += d * (SQR(i)+SQR(k));
				A(2,2) += d * (SQR(i)+SQR(j));

				A(0,1) -= d * i * j;
				A(0,2) -= d * i * k;
				A(1,2) -= d * j * k;
				sum += d;
			}
		}
	}
	//A /= sum;

	int info = lapack::syev('V', 'U', A, lambda, lapack::optimal_workspace());
	if(info!=0)
		return;

	unsigned int idx = distance(lambda.begin(),max_element(lambda.begin(),lambda.end()));
	if(idx!=0){
		std::swap(A(0,0),A(0,idx));
		std::swap(A(1,0),A(1,idx));
		std::swap(A(2,0),A(2,idx));
		std::swap(lambda(0),lambda(idx));
	}
	M = A;
	l = lambda;
}

#endif /* __VOXEL_NORMALS_HPP__ */
