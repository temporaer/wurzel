# vim:ts=4:sw=4:sts=4:et:ai
import os
import numpy as np
import cPickle
from scipy.ndimage import gaussian_filter

class dataset(object):
    def __init__(self,datafile,crop=False,usepickled=True,upsample=None,dtype=np.uint8,medianfilt=True,dz=120):
        if not isinstance(datafile,str):
            self.D = datafile
            return

        picklename = datafile.replace(".dat",".pickle")
        if usepickled and os.path.exists(picklename):
            self.load(picklename)
        else:
            with open(datafile) as fd:
                self.D = np.fromfile(file=fd, dtype=dtype).reshape((256,256,dz) ).astype("float32")/255.0
            if medianfilt:
                self.median_filter()
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
        print "Resampling 3rd axis..."
        #print "mm: 100 x 100 x 131"
        #print "Dims:", self.D.shape
        if method == "zoom":
            self.D = zoom(self.D, [1,1,256.0/120.0 * 100.0/100.0]).astype("float32")
        elif method == "resample":
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
