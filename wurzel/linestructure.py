import gc
import numpy as np

def cnan(x):
    if np.isnan(x).sum()>0:
        import pdb
        pdb.set_trace()

def get_curve_3D(l3,l2,l1, alpha=0.25,g23=0.5,g12=0.5): # renumerated according to sato et al: l3 is smallest
    print "Finding 3D curvilinear structures..."
    gc.collect()
    #alpha       = 0.25      # 0.0: cannot deal well with curved lines. Sato: 0.25
    #g23         = 0.5       # Sato: robust btw. 0.5,1.0. sharpness of
    #                        #       selectivity for cross-section isotropy
    #g12         = 0.5       # Sato: robust btw. 0.5,1.0. gives w12
    #                        #       asymmetrical characteristic in neg/pos regions of l1
    #                        #       l1<0 and |l1| large: blob, not line!

    l2a         = np.abs(l2)

    factA       = np.abs(l3) * (l2/l3)**g23
    factB       = (1+l1/l2a)**g12
    factC       = (1-alpha*l1/l2a)**g12

    A           = (l3<l2) * (l2<l1) * (l1 <= 0)
    B           = (l3<l2) * (l2<0) * (0<l1) * (l1<l2a/alpha)
    l123        = np.zeros(A.shape)
    l123[A]     = (factA*factB)[A]
    l123[B]     = (factA*factC)[B]
    cnan(l123)

    print "Done."

    return l123

