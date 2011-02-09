# vim:ts=4:sw=4:sts=4:et:ai
import os
import numpy as np
import cPickle
from scipy.ndimage import gaussian_filter

class dataset(object):
    def __init__(self,datafile,crop=False,usepickled=True,upsample=None,dtype=np.uint8,medianfilt=True,dz=120, remove_rohr=False):
        if not isinstance(datafile,str):
            self.D = datafile
            return

        picklename = datafile.replace(".dat",".pickle")
        if usepickled and os.path.exists(picklename):
            self.load(picklename)
        else:
            with open(datafile) as fd:
                self.D = np.fromfile(file=fd, dtype=dtype).reshape((256,256,dz) ).astype("float32")
            if dtype in [np.uint8, "uint8"]:
                self.D /= 255.0
            if medianfilt:
                self.median_filter()
            if remove_rohr:
                self.get_rid_of_roehrchen()
            self.upsample(upsample)
            self.save(picklename)
        if crop:
            self.D = self.D[50:200,50:200,10:80]
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
        if method == "zoom":
            print "Resampling 3rd axis..."
            self.D = zoom(self.D, [1,1,256.0/120.0 * 100.0/100.0]).astype("float32")
        elif method == "resample":
            print "Resampling 3rd axis..."
            self.D = resample(self.D, 120 * (256.0/120.0 * 100.0/100.0), axis=2, window=10).astype("float32")
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
        self.D.tofile(upsname)
    def get_rid_of_roehrchen(self):
        D = self.D
        X = D.sum(axis=2)
        Y = gaussian_filter(X, 5)
        p = np.argmax(Y)
        py = p % X.shape[0]
        px = p / X.shape[0]

        mv = Y[px,py]
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
