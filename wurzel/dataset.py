# vim:ts=4:sw=4:sts=4:et:ai
import os
import numpy as np
import cPickle
from scipy.ndimage import gaussian_filter
from config import config

class WurzelType:
    lupine, gerste, reis = range(3)
class WurzelInfo:
    def __init__(self, fn):
        cfg = config("dijkstra/config.xml")

        basename,ext = os.path.splitext(fn)
        self.read_shape = cfg.read_shape(basename)
        self.shape      = cfg.shape(basename)
        self.read_dtype = cfg.read_dtype(basename)
        self.has_rohr   = cfg.has_rohr(basename)

        print "WurzelInfo: ", self.read_shape, self.shape, self.read_dtype, self.has_rohr

        upsampled = fn.find("-upsampled")>=0
        sato      = fn.find(".sato")>=0
        if upsampled or sato:
            self.read_shape = self.shape
            self.read_dtype = "float32"

class dataset(object):
    def __init__(self,datafile,crop=False,usepickled=True,upsample=None,medianfilt=True,dz=120, remove_rohr=False):
        if not isinstance(datafile,str):
            self.D = datafile
            return
        self.info = WurzelInfo(datafile)
        info = self.info
        if not info.has_rohr: remove_rohr = False
        if not info.has_rohr: medianfilt = False

        picklename = datafile.replace(".dat",".pickle")
        if usepickled and os.path.exists(picklename):
            self.load(picklename)
            assert all([x==y for x,y in zip(self.D.shape, info.shape )])
        else:
            with open(datafile) as fd:
                self.D = np.fromfile(file=fd, dtype=info.read_dtype).reshape(info.read_shape).astype("float32")
            if info.read_dtype in [np.uint8, "uint8"]:
                self.D /= 255.0
            if medianfilt:  self.median_filter()
            if remove_rohr: self.get_rid_of_roehrchen()
            self.upsample(upsample)
            if medianfilt or remove_rohr or upsample:
                self.save(picklename)

    def median_filter(self):
        print "Median-Filtering..."
        D  = self.D
        x = np.median(np.median(D,axis=1),axis=1)
        for i in xrange(len(x)):
            D[i,:,:] -= x[i]
        self.D = D
        print "done."
    def upsample(self, method):
        from scipy.signal import resample
        from scipy.ndimage.interpolation import zoom
        #print "mm: 100 x 100 x 131"
        #print "Dims:", self.D.shape
        fact = np.array(self.info.shape).astype("float32") / np.array(self.info.read_shape).astype("float32")
        if method == "zoom":
            print "Resampling..."
            self.D = zoom(self.D, fact).astype("float32")
        elif method == "resample":
            print "Resampling..."
            a = self.info.resample_ax
            s = self.info.shape[a]
            self.D = resample(self.D, s, axis=a, window=10).astype("float32")
        elif method == None:
            pass
        else:
            raise NotImplementedError("Unknown upsampling method: %s" % method)
        #print "Dims:", self.D.shape
        print "done."
    def get_smoothed(self, sigma):
        return dataset(gaussian_filter(self.D,sigma))

    def load(self, picklename):
        with open(picklename) as f:
            self.D = cPickle.load(f)
    def save(self, datafile):
        with open(datafile,"wb") as f:
            cPickle.dump(self.D, f, cPickle.HIGHEST_PROTOCOL)
        upsname = datafile.replace(".pickle","-upsampled.dat")
        print "Saving to upsampled: ", self.D.shape
        self.D.tofile(upsname)
    def get_rid_of_roehrchen(self):
        D = self.D
        X = D.sum(axis=2)
        Y = gaussian_filter(X, 5)
        p = np.argmax(Y)
        py = p % X.shape[0]
        px = p / X.shape[0]

        S  = 0

        for i in xrange(5,20):
            ix, iy = np.mgrid[0:X.shape[0],0:X.shape[1]]
            iy -= py
            ix -= px
            circ = np.sqrt(iy**2+ix**2) < i
            if i == 5:
                S = Y[circ].sum()
                continue
            S2 = Y[circ].sum()
            #print "i:", i, " S2: ", S2, " S:", S, " ratio: ", (S2-S)/S
            if (S2-S)/S < 0.1:
                for z in xrange(self.D.shape[2]):
                    self.D[:,:,z][circ] = self.D.mean()
                break
            S = S2
