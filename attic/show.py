# vim:ts=4:sw=4:sts=4:et
import matplotlib.pyplot as plt
import cPickle
#import sys
import numpy as np
import gc
#import matplotlib.pyplot as plt
#from scipy.io.numpyio import fwrite, fread
from scipy.ndimage import gaussian_filter
from scipy.signal import resample

from enthought.mayavi import mlab
#from enthought.mayavi.modules.scalar_cut_plane import ScalarCutPlane

class DatasetProcessor(object):
    def __init__(self, datafile):
        if datafile.find(".dat")>0:
            with open(datafile) as fd:
                self.D = np.fromfile(file=fd, dtype=np.uint8).reshape((256,256,120) ).astype("float32")/255.0
            self.D = self.D[50:200,50:200,10:80]
            self.median_filter()
            self.upsample()
            self.save(datafile.replace(".dat",".pickle"))
        else:
            with open(datafile) as f:
                self.D = cPickle.load(f)
    def save(self, datafile):
        with open(datafile,"wb") as f:
            cPickle.dump(self.D, f, cPickle.HIGHEST_PROTOCOL)
    def upsample(self):
        print "Resampling 3rd axis..."
        print "mm: 100 x 100 x 131"
        print "Dims:", self.D.shape
        from scipy.ndimage.interpolation import zoom
        self.D = zoom(self.D, [1,1,256.0/120.0 * 100.0/100.0])
        #self.D = resample(self.D, 120 * (256.0/120.0 * 131.0/100.0), axis=2, window=("gaussian", 40))
        #self.D = resample(self.D, 120 * (256.0/120.0 * 100.0/100.0), axis=2, window=10)
        print "Dims:", self.D.shape
        print "done."

    def median_filter(self):
        D  = self.D
        x = np.median(np.median(D,axis=1),axis=1)
        for i in xrange(len(x)):
            D[i,:,:] -= x[i]
        self.D = D

    def show_found(self, cm="Spectral"):
        print "Show structure"
        D = self.S
        src = mlab.pipeline.scalar_field(D)
        #if D.shape[2]<D.shape[1]:
        #    src.image_data.spacing = (1.0, 1.0, 256.0/120.0*131.0/100.0)
        #mlab.pipeline.volume(src, vmin=D.min()+0.1*D.ptp(),vmax=D.max()-0.1*D.ptp())

        #mlab.contour3d(D)
        #mlab.pipeline.volume(src, vmin=D.min()+0.2*D.ptp(),vmax=D.max()-0.2*D.ptp())
        mlab.pipeline.iso_surface(src, contours=[D.min()+0.1*D.ptp(), ], opacity=0.3, colormap=cm)
        #mlab.pipeline.iso_surface(src, contours=[D.max()-0.2*D.ptp(), ], opacity=0.2)
        #mlab.pipeline.image_plane_widget(src,
        #                        plane_orientation='z_axes',
        #                        slice_index=80,
        #                    )
        #mlab.pipeline.image_plane_widget(src,
        #                        plane_orientation='x_axes',
        #                        slice_index=10,
        #                    )
        #mlab.outline()
        print "done"
    def show_imgs(self):
        D = self.D
        #print "zooming in"
        #from scipy.ndimage.morphology import morphological_gradient
        #from scipy.ndimage.interpolation import zoom
        #D = zoom(self.D, [1,1,256.0/120.0], order=1)
        #D = morphological_gradient(D, size=(3,3,3))
        #D = np.log(D+1)
        #D = gaussian_filter(D, 2.0)
        print "show data, shape=",D.shape, D.flags
        src = mlab.pipeline.scalar_field(D)
        #if D.shape[2]<D.shape[1]:
        #    src.image_data.spacing = (1.0, 1.0, 256.0/120.0*131.0/100.0)
        #mlab.pipeline.volume(src, vmin=D.min()+0.15*D.ptp(),vmax=D.max()-0.4*D.ptp())

        #mlab.contour3d(D)
        #mlab.pipeline.volume(src, vmin=D.min()+0.2*D.ptp(),vmax=D.max()-0.2*D.ptp())
        mlab.pipeline.iso_surface(src, contours=[D.min()+0.20*D.ptp(), ], opacity=0.5, colormap="bone")
        #mlab.pipeline.iso_surface(src, contours=[D.max()-0.2*D.ptp(), ], opacity=0.2)
        #mlab.pipeline.image_plane_widget(src,
        #                        plane_orientation='z_axes',
        #                        slice_index=80,
        #                    )
        #mlab.pipeline.image_plane_widget(src,
        #                        plane_orientation='x_axes',
        #                        slice_index=10,
        #                    )
        #mlab.outline()
        print "done"
def get_ev(D, scale=3):
    from structure_tensor import hessian, eig3x3
    print "Calculating Hessian"
    res = hessian(D, scale)
    print "Calculating Symev"
    l1,l2,l3 = eig3x3(res)
    print "done"
    return l1,l2,l3

def save_lambdas(L):
    print "Saving Lambdas"
    import cPickle
    with open("lambdas.pickle", "wb") as f:
        cPickle.dump(L, f, cPickle.HIGHEST_PROTOCOL)
def load_lambdas():
    print "Loading Lambdas"
    import cPickle
    with open("lambdas.pickle") as f:
        return cPickle.load(f)

def cnan(x):
    if np.isnan(x).sum()>0:
        import pdb
        pdb.set_trace()
def select_struct(l3,l2,l1, alpha=0.25,g23=0.5,g12=0.5): # renumerated according to sato et al: l3 is smallest
    print "Selecting Structure Elements..."
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

def run():
    dp = DatasetProcessor("L2_22aug.dat")
    #dp = DatasetProcessor("L2_22aug.pickle")
    dp.D -= dp.D.min()

    mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    mlab.clf()
    #dp.show_imgs()

    sigmas  = [1.5, 3.0]
    sigmas2 = [0.0, 0.0]
    cms     = ["bone", "Spectral"]
    for sigma, sigma2, cm in zip(sigmas,sigmas2,cms):
        gc.collect()
        if True:
            D = gaussian_filter(dp.D, sigma)
            l1,l2,l3 = get_ev(D, sigma2)
            save_lambdas([l1,l2,l3])
        else:
            l1,l2,l3 = load_lambdas()
        dp.S  = select_struct(l1,l2,l3,0.25,0.5,0.5)
        gc.collect()
        dp.S -= dp.S.min()
        dp.show_found(cm)
    plt.show()

if __name__ == "__main__":
    run()
