# vim:ts=4:sw=4:sts=4:et:ai
import numpy as np
import scipy.ndimage as nd
from scipy.weave import inline, converters

def grad(img, ax):
    return nd.convolve1d(img, weights=[-.5,.5], axis=ax)

def hessian(img, sigma=-1):
    assert img.ndim == 3
    Dx  = grad(img, 0)
    Dy  = grad(img, 1)
    Dz  = grad(img, 2)

    Dxx = grad(Dx, 0)
    Dyy = grad(Dy, 1)
    Dzz = grad(Dz, 2)

    Dxy = grad(Dx, 1)
    Dxz = grad(Dx, 2)

    Dyz = grad(Dy, 2)
    D = {"Dxx": Dxx, "Dyy": Dyy, "Dzz": Dzz, "Dxy": Dxy, "Dxz": Dxz, "Dyz": Dyz}
    if sigma>0:
        for k,v in D.items():
            D[k] = nd.gaussian_filter(v, sigma)
    return D

def eig3x3(hessian):
    nx,ny,nz = hessian["Dxx"].shape
    lambda1  = np.empty(hessian["Dxx"].shape)
    lambda2  = np.empty(hessian["Dxx"].shape)
    lambda3  = np.empty(hessian["Dxx"].shape)
    code = """
      #line 37 "structure_tensor.py"
      namespace ublas = boost::numeric::ublas;
      namespace lapack = boost::numeric::bindings::lapack;
      using namespace std;
      using ublas::prod;
      using ublas::trans;
      typedef ublas::bounded_matrix<double,3,3,ublas::column_major> mat;
      ublas::bounded_vector<double,3> lambda;

      //typedef ublas::matrix<double,ublas::column_major> mat;
      //ublas::vector<double> lambda(3);
      mat A;
      int i,j,k;
      #pragma omp parallel for private(i,j,k,A,lambda)
      for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            for(k=0;k<nz;k++){
              A(0,0) = Dxx(i,j,k);
              A(1,1) = Dyy(i,j,k);
              A(2,2) = Dzz(i,j,k);
              A(0,1) = Dxy(i,j,k);
              A(0,2) = Dxz(i,j,k);
              A(1,2) = Dyz(i,j,k);

              lapack::syev('N','U',A,lambda,lapack::minimal_workspace());
              sort(lambda.begin(),lambda.end());

              lambda1(i,j,k) = lambda(0);
              lambda2(i,j,k) = lambda(1);
              lambda3(i,j,k) = lambda(2);
            }
          }
        }
        """
    V = vars()
    variables = "nx ny nz lambda1 lambda2 lambda3".split()
    variables.extend(hessian.keys())
    map(lambda x:V.__setitem__(x[0],x[1]), hessian.items())
    inline(code, variables,
                 compiler="gcc",
                 extra_compile_args =['-O3 -fopenmp'],
                 extra_link_args=['-lgomp'],
                 headers=[
                          '<boost/mpl/or.hpp>',
                          '<boost/numeric/ublas/fwd.hpp>',
                          '<boost/numeric/bindings/lapack/syev.hpp>',
                          '<boost/numeric/bindings/traits/ublas_matrix.hpp> ',
                          '<boost/numeric/bindings/traits/ublas_vector.hpp> ',
                          '<boost/numeric/ublas/matrix.hpp>',
                          '<boost/numeric/ublas/banded.hpp>',
                          '<boost/numeric/ublas/vector.hpp>',
                          '<boost/numeric/ublas/io.hpp>',
                          '<iostream>',
                          '<cmath>'],
                 type_converters=converters.blitz,
                 libraries=["lapack"])
    return lambda1, lambda2, lambda3

def get_ev_of_hessian(D,sigma=-1):
    print "Calculating Hessian"
    res = hessian(D, sigma)
    print "Calculating Symev"
    l1,l2,l3 = eig3x3(res)
    print "done"
    return l1,l2,l3
