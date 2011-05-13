#ifndef __VOXEL_ACCESSORS_HPP__
#define __VOXEL_ACCESSORS_HPP__
//=======================================================================
// Copyright 2011 University Bonn
// Author: Hannes Schulz
//=======================================================================

#include "voxelgrid.hpp"
#include "wurzel_tree.hpp"

#define USE_MANUAL_ACCESSOR 0
#if USE_MANUAL_ACCESSOR
#  define ACC(A,v) A.data()[v[0] * (A.shape()[1]*A.shape()[2]) + v[1]*A.shape()[2] + v[2]]
#else
#  define ACC(A,v) A[v[0]][v[1]][v[2]]
#endif
#define SQR(X) ((X)*(X))

template<class T>
struct const_vox2arr{
	const T& A;
	const_vox2arr(const T& a):A(a){}
	const typename T::element& operator[](const voxel_vertex_descriptor&v)const{
		return ACC(A,v);
	}
};

template<class T>
struct vox2arr{
	T& A;
	vox2arr(T& a):A(a){}
	typename T::element& operator[](const voxel_vertex_descriptor&v)const{
		return ACC(A,v);
	}
};

template<class T>
vox2arr<T> make_vox2arr(T& t){ return vox2arr<T>(t); }

template<class T>
const_vox2arr<T> make_vox2arr(const T& t){ return const_vox2arr<T>(t); }

template<class T>
typename T::element 
get_subpixel_value(const T& A, const float& ffx, const float& ffy, const float& ffz){
			const typename T::size_type* shape = A.shape();
		  // taken from CIMg.h by David Tschumperle
      const float fx = ffx<0?0:(ffx>shape[0]-1?shape[0]-1:ffx),
					        fy = ffy<0?0:(ffy>shape[1]-1?shape[1]-1:ffy),
								 	fz = ffz<0?0:(ffz>shape[2]-1?shape[2]-1:ffz);
      const unsigned int x = (unsigned int)fx, y = (unsigned int)fy, z = (unsigned int)fz;
      const float dx = fx-x, dy = fy-y, dz = fz-z;
      const unsigned int nx = dx>0?x+1:x, ny = dy>0?y+1:y, nz = dz>0?z+1:z;
      const typename T::element
        &Iccc = A[x][y][z],  &Incc = A[nx][y][z],  &Icnc = A[x][ny][z],  &Innc = A[nx][ny][z],
        &Iccn = A[x][y][nz], &Incn = A[nx][y][nz], &Icnn = A[x][ny][nz], &Innn = A[nx][ny][nz];
      return Iccc + dx*(Incc-Iccc) + dy*(Icnc-Iccc) + dz*(Iccn-Iccc) +
        dx*dy*(Iccc+Innc-Icnc-Incc) + dx*dz*(Iccc+Incn-Iccn-Incc) + dy*dz*(Iccc+Icnn-Iccn-Icnc) +
        dx*dy*dz*(Iccn+Innn+Icnc+Incc-Icnn-Incn-Iccc-Innc);
}

template<class T>
struct const_vox2arr_subpix{
	const T& A;
	const_vox2arr_subpix(const T& a):A(a){}
	const typename T::element operator[](const vec3_t&v)const{
		const float ffx = v[0], ffy = v[1], ffz = v[2];
		return get_subpixel_value(A, ffx, ffy, ffz);
	}
	const typename T::element operator()(const float& ffx, const float& ffy, const float& ffz)const{
		return get_subpixel_value(A, ffx, ffy, ffz);
	}
};

template<class T>
struct vox2arr_subpix{
	T& A;
	vox2arr_subpix(T& a):A(a){}
	typename T::element operator[](const vec3_t&v)const{
		const float ffx = v[0], ffy = v[1], ffz = v[2];
		return get_subpixel_value(A, ffx, ffy, ffz);
	}
	typename T::element operator()(const float& ffx, const float& ffy, const float& ffz)const{
		return get_subpixel_value(A, ffx, ffy, ffz);
	}
};

template<class T>
vox2arr_subpix<T> make_vox2arr_subpix(T& t){ return vox2arr_subpix<T>(t); }

template<class T>
const_vox2arr_subpix<T> make_vox2arr_subpix(const T& t){ return const_vox2arr_subpix<T>(t); }

template<class I, class J>
inline double voxdist(I a, J b){
	double s = 0;
	s += SQR(*a-*b); a++; b++;
	s += SQR(*a-*b); a++; b++;
	s += SQR(*a-*b); 
	return sqrt(s);
}

#endif /* __VOXEL_ACCESSORS_HPP__ */
