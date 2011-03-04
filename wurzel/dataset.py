# vim:ts=4:sw=4:sts=4:et:ai
import os
import numpy as np
import cPickle
from scipy.ndimage import gaussian_filter

class WurzelType:
    lupine, gerste, reis = range(3)
class WurzelInfo:
    def __init__(self, fn):
        self.wurzel_type = self.gettype(fn)
        wt = self.wurzel_type

        if wt==WurzelType.lupine:
            self.read_shape = (256,256,120)
            self.shape      = (256,256,256)
            self.read_dtype = "uint8"
            self.resample_ax= 2
        elif wt==WurzelType.gerste:
            self.read_shape = (410,192,192)
            self.shape      = (872,192,192)
            self.read_dtype = "float32"
            self.resample_ax= 0

        upsampled = fn.find("-upsampled")>=0
        if upsampled:
            self.read_shape = self.shape
            self.read_dtype = "float32"
    def gettype(self,fn):
        if fn.find("reis")>=0:
            return WurzelType.reis
        if fn.find("Gerste")>=0:
            return WurzelType.gerste
        if fn.find("L2")>=0:
            return WurzelType.lupine

class dataset(object):
    def __init__(self,datafile,crop=False,usepickled=True,upsample="zoom",medianfilt=True,dz=120, remove_rohr=False):
        if not isinstance(datafile,str):
            self.D = datafile
            return
        self.info = WurzelInfo(datafile)
        info = self.info
        if info.wurzel_type != WurzelType.lupine: remove_rohr = False
        if info.wurzel_type != WurzelType.lupine: medianfilt = False

        picklename = datafile.replace(".dat",".pickle")
        if usepickled and os.path.exists(picklename):
            self.load(picklename)
        else:
            with open(datafile) as fd:
                self.D = np.fromfile(file=fd, dtype=info.read_dtype).reshape(info.read_shape).astype("float32")
            if info.read_dtype in [np.uint8, "uint8"]:
                self.D /= 255.0
            if medianfilt:  self.median_filter()
            if remove_rohr: self.get_rid_of_roehrchen()
            self.upsample(upsample)
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
