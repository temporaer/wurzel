# vim:ts=4:sw=4:sts=4:et:ai
import os
import numpy as np
import scipy.ndimage as nd
from scipy.weave import inline, converters

def grad(img, ax):
    return nd.convolve1d(img, weights=[-.5,.0,.5], axis=ax)

def hessian(img, sigma,gamma=1):
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
    for k in D:
        v = D[k]
        v[:] *= sigma**(2*gamma)  # we derived twice (!)
    return D

def eig3x3(hessian):
    nx,ny,nz = hessian["Dxx"].shape
    lambda1  = np.empty(hessian["Dxx"].shape,dtype="float32")
    lambda2  = np.empty(hessian["Dxx"].shape,dtype="float32")
    lambda3  = np.empty(hessian["Dxx"].shape,dtype="float32")

    ev10     = np.empty(hessian["Dxx"].shape,dtype="float32")
    ev11     = np.empty(hessian["Dxx"].shape,dtype="float32")
    ev12     = np.empty(hessian["Dxx"].shape,dtype="float32")
    code = """
      #line 40 "structure_tensor.py"
      namespace ublas = boost::numeric::ublas;
      namespace lapack = boost::numeric::bindings::lapack;
      using namespace std;
      using ublas::prod;
      using ublas::trans;
      typedef ublas::bounded_matrix<double,3,3,ublas::column_major> mat;
      ublas::bounded_vector<double,3> lambda;
      ublas::bounded_vector<double,34*3> work; // the 34*width is from syev.hpp

      mat A;
      int i,j,k,idx;
      #pragma omp parallel for private(i,j,k,A,idx,lambda,work)
      for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            for(k=0;k<nz;k++){
              A(0,0) = Dxx(i,j,k);
              A(1,1) = Dyy(i,j,k);
              A(2,2) = Dzz(i,j,k);
              A(0,1) = Dxy(i,j,k);
              A(0,2) = Dxz(i,j,k);
              A(1,2) = Dyz(i,j,k);

              //lapack::syev('V','U',A,lambda,lapack::optimal_workspace());  // V/N compute/donotcompute eigenvectors 
                lapack::syev('N','U',A,lambda,lapack::workspace(work));
              idx = distance(lambda.begin(),min_element(lambda.begin(),lambda.end(), cmp_abs()));
              ev10(i,j,k)    = A(0,idx);
              ev11(i,j,k)    = A(1,idx);
              ev12(i,j,k)    = A(2,idx);

              sort(lambda.begin(),lambda.end(),cmp_abs());

              lambda1(i,j,k) = lambda(0);  // eigenvalues in ascending order
              lambda2(i,j,k) = lambda(1);
              lambda3(i,j,k) = lambda(2);
            }
          }
        }
        """
    sc = """
    struct cmp_abs{
      bool operator()(double a, double b){
          return fabs(a)<fabs(b);
      }
    };
    """
    V = vars()
    variables = "nx ny nz lambda1 lambda2 lambda3 ev10 ev11 ev12".split()
    variables.extend(hessian.keys())
    map(lambda x:V.__setitem__(x[0],x[1]), hessian.items())
    inline(code, variables,
                 verbose=2,
                 compiler="gcc",
                 support_code=sc,
                 extra_compile_args =['-O3 -fopenmp'],
                 #extra_compile_args =['-O3'],
                 extra_link_args=['-lgomp'],
                 include_dirs=["%s/pool/include/boost-numeric-bindings"%os.environ["HOME"]],
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
                          '<cmath>',
                          '<algorithm>',
                          '<vector>',
                          '<cmath>'],
                 type_converters=converters.blitz,
                 libraries=["lapack"])
    assert np.isnan(ev10).sum()==0
    assert np.isnan(ev11).sum()==0
    assert np.isnan(ev12).sum()==0
    return {"lambda1":lambda1,"lambda2": lambda2, "lambda3":lambda3, "ev10": ev10, "ev11": ev11, "ev12": ev12}

def get_ev_of_hessian(D,sigma,gamma=1):
    print "Calculating Hessian"
    res = hessian(D, sigma, gamma)
    print "Calculating Symev"
    eig = eig3x3(res)
    print "done"
    return eig
