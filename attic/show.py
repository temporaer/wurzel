# vim:ts=4:sw=4:sts=4:et
#import matplotlib.pyplot as plt
#import cPickle
#import sys
import numpy as np
#import matplotlib.pyplot as plt
#from scipy.io.numpyio import fwrite, fread
from scipy.ndimage import gaussian_filter

from enthought.mayavi import mlab
#from enthought.mayavi.modules.scalar_cut_plane import ScalarCutPlane

class DatasetProcessor(object):
    def __init__(self, videofile):
        import conv
        self.D = conv.mhd2npy(videofile)
        #with open(videofile) as fd:
        #    self.D = np.fromfile(file=fd, dtype=np.uint8).reshape((256,256,120) ).astype("float32")/255.0

    def show_spectral(self):
        from scipy.ndimage.interpolation import zoom
        from scipy.ndimage import gaussian_filter
        self.D = zoom(self.D, 0.5)
        self.D = gaussian_filter(self.D, [2,2,1])
        from scikits.learn.feature_extraction import image
        from scikits.learn.cluster import spectral_clustering
        mask = self.D>(self.D.min()+0.0*self.D.ptp())
        graph = image.img_to_graph(self.D, mask=mask)
        graph.data = np.exp(-graph.data/graph.data.std())
        labels = spectral_clustering(graph, k=3)
        label_im = -np.ones(mask.shape)
        label_im[mask] = labels

        D = label_im
        src = mlab.pipeline.scalar_field(D)
        src.image_data.spacing = (1.0, 1.0, 2.0)
        mlab.pipeline.volume(src, vmin=D.min()+0.1*D.ptp(),vmax=D.max()-0.1*D.ptp())

        #mlab.pipeline.iso_surface(src, contours=[0.5])
        #mlab.pipeline.iso_surface(src, contours=[1.5])
        #mlab.pipeline.iso_surface(src, contours=[2.5])

    def show_meanshift(self):
        print "Constructing Matrix..."
        #from scikits.learn.cluster import MeanShift, estimate_bandwidth
        from scikits.learn.cluster import MeanShift
        from scipy.ndimage.interpolation import zoom
        from scipy.ndimage import gaussian_filter
        self.D = zoom(self.D, 0.2)
        self.D = gaussian_filter(self.D, [2,2,1])
        D = self.D
        X,Y,Z = np.mgrid[:self.D.shape[0]:1,:self.D.shape[1]:1,:self.D.shape[2]:1]
        D = np.vstack([x.flatten() for x in [D,X,Y,Z]]).T
        #print "Estimating bandwidth"
        #bandwidth = estimate_bandwidth(D, quantile=0.3)
        bandwidth = 1
        print "Running MeanShift"
        ms = MeanShift(bandwidth=bandwidth)
        ms.fit(D)
        print "Done."
        labels = ms.labels_
        #cluster_centers = ms.cluster_centers_
        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)
        print "number of estimated clusters : %d" % n_clusters_

        D = labels.reshape(D.shape).copy("F")
        src = mlab.pipeline.scalar_field(D)
        src.image_data.spacing = (1.0, 1.0, 2.0)
        mlab.pipeline.volume(src, vmin=D.min()+0.1*D.ptp(),vmax=D.max()-0.1*D.ptp())


    def show_found_overlay(self):
        print "read found"
        import conv
        D = conv.mhd2npy("/home/local/l_schulz/tmp/build-itk/InsightApplications-3.20.0/build/Overlay.mhd").copy("C")
        print "show found, shape=", D.shape, D.flags
        src = mlab.pipeline.scalar_field(D)

        from enthought.tvtk.util.ctf import ColorTransferFunction
        ctf = ColorTransferFunction()
        ctf.add_rgb_point(0, 1, 1, 1)
        ctf.add_rgb_point(1, 0, 0, 0)

        #vol            = mlab.pipeline.volume(src, vmin=0.1,vmax=0.9)
        #vol._volume_property.set_color(ctf)
        #vol._ctf       = ctf
        #vol.update_ctf = True
        mlab.pipeline.iso_surface(src, contours=[0.2,], opacity=0.5, colormap="Spectral")
        print "done"
    def show_found(self):
        print "read found"
        D = np.zeros((256,256,256))
        with open("ImageSpace-EigenSpaceMap.txt", "r") as f:
            import re
            for s in f:
                r = r'\[(\d+), (\d+), (\d+)\].*\[(.*), (.*), (.*)\]'
                m = re.search(r,s)
                I = [int(i)   for i in m.groups()[:3]]
                #V = [float(f) for f in m.groups()[3:]]
                D[I[0],I[1],I[2]] = 1
        print "show found"
        src = mlab.pipeline.scalar_field(D)
        #src.image_data.spacing = (1.0, 1.0, 2.0)

        from enthought.tvtk.util.ctf import ColorTransferFunction
        ctf = ColorTransferFunction()
        ctf.add_rgb_point(0, 1, 1, 1)
        ctf.add_rgb_point(1, 0, 0, 0)

        #vol            = mlab.pipeline.volume(src, vmin=0.1,vmax=0.9)
        #vol._volume_property.set_color(ctf)
        #vol._ctf       = ctf
        #vol.update_ctf = True
        mlab.pipeline.iso_surface(src, contours=[0.2,], opacity=0.3)
        print "done"

    def show_imgs(self):
        #print "zooming in"
        #from scipy.ndimage.morphology import morphological_gradient
        #from scipy.ndimage.interpolation import zoom
        #D = zoom(self.D, [1,1,256.0/120.0], order=1)
        #D = morphological_gradient(D, size=(3,3,3))
        #D = self.D
        #D = np.log(D+1)
        #D = gaussian_filter(D, 2.0)
        print "show data, shape=",D.shape, D.flags
        src = mlab.pipeline.scalar_field(D)
        #src.image_data.spacing = (1.0, 1.0, 256.0/120.0)
        #mlab.pipeline.volume(src, vmin=D.min()+0.1*D.ptp(),vmax=D.max()-0.1*D.ptp())

        #mlab.contour3d(D)
        #mlab.pipeline.volume(src, vmin=D.min()+0.2*D.ptp(),vmax=D.max()-0.2*D.ptp())
        mlab.pipeline.iso_surface(src, contours=[D.min()+0.15*D.ptp(), ], opacity=0.5, colormap="bone")
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

def select_struct(l3,l2,l1): # renumerated according to sato et al: l3 is smallest
    A           = -((l3<l2) * (l2<0))
    lmin23      = -np.maximum(l2,l3)
    lmin23[A]   = 0

    lgmean23    = np.sqrt(l2*l3)
    lgmean23[A] = 0

    g23         = 0.5
    l23         = np.abs(l3) * (l2/l3)**g23
    l23[A]      = 0

    alpha       = 0.25
    g12         = 0.5
    l2a         = np.abs(l2)
    l123        = np.zeros(l23.shape)
    B           = l1 <= 0
    l123[B]     = ((1+l1/l2a)**g12)[B]
    B           = (l2a/alpha>l1) * (l1 > 0)
    l123[B]     = ((1-alpha * l1/l2a)**g12)[B]
    l123       *= l23

    return l123

if __name__ == "__main__":
    dp = DatasetProcessor("/home/local/l_schulz/tmp/build-itk/InsightApplications-3.20.0/build/bla.mhd")
    D  = dp.D
    x = np.median(np.median(D,axis=1),axis=1) 
    for i in xrange(len(x)):
        D[i,:,:] -= x[i]
    if True:
        D = gaussian_filter(D, 2.1)
        l1,l2,l3 = get_ev(D)
        save_lambdas([l1,l2,l3])
    else:
        l1,l2,l3 = load_lambdas()
    D -= D.min()
    S = select_struct(l1,l2,l3)
    dp.D = S
    #plt.matshow(D.sum(axis=0))
    #plt.matshow(D.sum(axis=1))
    #plt.matshow(D.sum(axis=2))
    #plt.plot(x,"+-")
    #plt.show()
    mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    mlab.clf()
    dp.show_imgs()
    #dp.show_found()
    #dp.show_found_overlay()
    #dp.show_spectral()
    #dp.show_meanshift()
