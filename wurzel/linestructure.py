import numpy as np

# Copyright 2011 University of Bonn
# Author: Hannes Schulz

def cnan(x):
    """ check for not-a-number in parameter x """
    if np.isnan(x).sum()>0:
        import pdb
        pdb.set_trace()

def get_curve_3D(eig, alpha=0.25,g23=0.5,g12=0.5): # renumerated according to sato et al: l3 is smallest
    """
    return the vesselness measure based on eigenvalues of normalized hessians in eig
    """
    #return sato(eig,alpha,g23, g12)
    return frangi(eig)

def frangi_memory_hungry(eig):
    print "Finding 3D curvilinear structures using Frangi et al..."
    l1     = eig["lambda1"]  # in fangi et al, we have |l1| <= |l2| <= |l3|
    l2     = eig["lambda2"]
    l3     = eig["lambda3"]
    alpha  = 0.5
    beta   = 0.5
    RB     = np.abs(l1) / np.sqrt(np.abs(l2 * l3))
    RA     = np.abs(l2) / np.abs(l3)
    S      = np.sqrt(l1**2+l2**2+l3**2)
    c      = np.max(S)/3.   # half of maximum frobenius norm of hessian says fangi et al.
    V      = (1-np.exp(RA**2)/(-2. * alpha**2)) * np.exp(RB**2/(-2. * beta**2)) * (1-np.exp(S**2/(-2. * c**2)))
    P      = (l2 >= 0.) + (l3 >= 0.)
    V[P]   = 0
    V[V<0] = 0
    print "done."
    return V

def frangi(eig):
    print "Finding 3D curvilinear structures using Frangi et al..."
    l1     = eig["lambda1"]  # in fangi et al, we have |l1| <= |l2| <= |l3|
    l2     = eig["lambda2"]
    l3     = eig["lambda3"]
    alpha  = 0.5
    beta   = 0.5
    RB     = np.abs(l1)
    RB    /= np.sqrt(np.abs(l2*l3))

    RA     = np.abs(l2)
    RA    /= np.abs(l3)

    S      = l1.copy()
    S    **=     2
    S     += l2**2
    S     += l3**2
    np.sqrt(S,S)
    c      = np.max(S)/3.   # half of maximum frobenius norm of hessian says fangi et al.
    V      = 1-np.exp(RA**2)
    V     /= -2. * alpha**2
    V     *= np.exp(RB**2/(-2. * beta**2))
    V     *= 1-np.exp(S**2/(-2. * c**2))
    #V      = (1-np.exp(RA**2)/(-2. * alpha**2)) * np.exp(RB**2/(-2. * beta**2)) * (1-np.exp(S**2/(-2. * c**2)))
    P      = (l2 >= 0.) + (l3 >= 0.)
    V[P]   = 0
    V[V<0] = 0
    print "done."
    return V

def sato(eig, alpha=0.25,g23=0.5,g12=0.5): # renumerated according to sato et al: l3 is smallest
    print "Finding 3D curvilinear structures..."
    l3      = eig["lambda1"]  # in sato et al, l3 < l2 < l1; but lapack gives them in ascending order(!)
    l2      = eig["lambda2"]
    l1      = eig["lambda3"]
    #gc.collect()
    #alpha  = 0.25      # 0.0: cannot deal well with curved lines. Sato: 0.25
    #g23    = 0.5       # Sato: robust btw. 0.5,1.0. sharpness of
    #                        #       selectivity for cross-section isotropy
    #g12    = 0.5       # Sato: robust btw. 0.5,1.0. gives w12
    #                        #       asymmetrical characteristic in neg/pos regions of l1
    #                        #       l1<0 and |l1| large: blob, not line!

    l2a     = np.abs(l2)

    factA   = np.abs(l3) * (l2/l3)**g23
    factB   = (1+l1/l2a)**g12
    factC   = (1-alpha*l1/l2a)**g12

    A       = (l3<l2) * (l2<l1) * (l1 <= 0)
    B       = (l3<l2) * (l2<0) * (0<l1) * (l1<l2a/alpha)
    l123    = np.zeros(A.shape, dtype                                                                    = "float32")
    l123[A] = (factA*factB)[A]
    l123[B] = (factA*factC)[B]
    cnan(l123)

    print "Done."

    return l123

