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
get_normal(covmat_t&M,	const vec3_t& v, const T& acc){
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
	for(int i=-1;i<=1;i++){
		for(int j=-1;j<=1;j++){
			for(int k=-1;k<=1;k++){
				double d = acc(x+i,y+j,z+k);
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

	int info = lapack::syev('V', 'U', A, lambda, lapack::optimal_workspace());
	if(info!=0)
		return;

	unsigned int idx = distance(lambda.begin(),max_element(lambda.begin(),lambda.end()));
	if(idx!=0){
		std::swap(A(0,0),A(0,idx));
		std::swap(A(1,0),A(1,idx));
		std::swap(A(2,0),A(2,idx));
	}
	M = A;
}

#endif /* __VOXEL_NORMALS_HPP__ */
