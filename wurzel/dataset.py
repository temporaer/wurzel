# vim:ts=4:sw=4:sts=4:et:ai

# Copyright 2011 University of Bonn
# Author: Hannes Schulz

import os
import numpy as np
import cPickle
from scipy.ndimage import gaussian_filter
from config import config

class WurzelInfo:
    """ determines and provides access to commonly used features of root files.
       The data is determined based on the basename of the file as a key
       to the information stored in config.xml.
    """
    def __init__(self, fn, config_file="config.xml"):
        # load config file
        cfg = config(config_file)
        self.datapath = cfg.datapath.text

        # remove extensions
        basename,ext = os.path.splitext(fn)
        # remove common suffixes
        basename = basename.replace("-upsampled","")
        basename = basename.replace("-vertices","")
        try:
            self.read_shape = cfg.read_shape(basename)
            self.shape      = cfg.shape(basename)
            self.read_dtype = cfg.read_dtype(basename)
            self.has_rohr   = cfg.has_rohr(basename)
            self.scale      = cfg.scale(basename)
        except:
            print "Could not find dataset `%s' in config!"%basename
            import sys
            sys.exit(1)

        print "WurzelInfo: ", self.read_shape, self.shape, self.read_dtype, self.has_rohr

        upsampled = fn.find("-upsampled")>=0
        sato      = fn.find(".sato")>=0
        if upsampled or sato:
            self.read_shape = self.shape
            self.read_dtype = "float32"

class dataset(object):
    def __init__(self,datafile,crop=False,usepickled=True,upsample=None,medianfilt=True,remove_rohr=False):
        """
        Loads and stores the actual data and info about it.
        @param crop         shape to crop this to (if you want to try algorithm a lot, make dataset smaller to reduce time)
        @param usepickled   when true, do not load original data, instead take previously loaded+preprocessed pickled data
        @param upsample     whether to upsample the data according to information in WurzelInfo
        @param medianfilt   remove "sheet" structures occurring in Lupine data
        @param remove_rohr  remove reference tube in Lupine data
        """
        if not isinstance(datafile,str):
            self.D = datafile
            return
        self.info = WurzelInfo(datafile)
        info = self.info
        if not info.has_rohr: remove_rohr = False
        if not info.has_rohr: medianfilt = False

        picklename = os.path.join(info.datapath, datafile.replace(".dat",".pickle"))
        if usepickled and os.path.exists(picklename):
            self.load(picklename)
            if not all([x==y for x,y in zip(self.D.shape, info.shape )]):
                print "After loading pickle, dimensions do not match: ", self.D.shape, info.shape
                import sys
                sys.exit(1)
        else:
            with open(os.path.join(info.datapath, datafile)) as fd:
                self.D = np.fromfile(file=fd, dtype=info.read_dtype).reshape(info.read_shape).astype("float32")
            if info.read_dtype in [np.uint8, "uint8"]:
                self.D /= 255.0
            if medianfilt:  self.median_filter()
            if remove_rohr: self.get_rid_of_roehrchen()
            #assert self.D.min()>= 0
            self.D[self.D<0]=0
            self.upsample(upsample)
            if not medianfilt:
                cnt = (self.D<0).sum()
                print "fraction below zero: ", cnt/np.prod(self.D.shape)
                self.D[self.D<0]=0 # this is an upsampling-artefact (hopefully)
            if not all([x==y for x,y in zip(self.D.shape, info.shape )]):
                print "After resampling, dimensions do not match: ", self.D.shape, info.shape
                import sys
                sys.exit(1)

            if medianfilt or remove_rohr or upsample:
                self.save(picklename)

    def median_filter(self):
        """ filter sheet artifacts in Lupine data """
        print "Median-Filtering..."
        D  = self.D
        x = np.median(np.median(D,axis=1),axis=1)
        for i in xrange(len(x)):
            D[i,:,:] -= x[i]
        self.D = D
        print "done."
    def upsample(self, method):
        """ change resolution to enable processing with isotropic filters """
        from scipy.signal import resample
        from scipy.ndimage.interpolation import zoom
        #print "mm: 100 x 100 x 131"
        #print "Dims:", self.D.shape
        fact = np.array(self.info.shape).astype("float32") / np.array(self.info.read_shape).astype("float32")+0.00001 # hrmpf!!
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
        """ return a gauss-filtered version of the data """
        return dataset(gaussian_filter(self.D,sigma))

    def load(self, picklename):
        """ load pickled data """
        with open(picklename) as f:
            self.D = cPickle.load(f)
    def save(self, datafile):
        """ save a pickle, and an -upsampled.dat raw data file """
        with open(datafile,"wb") as f:
            cPickle.dump(self.D, f, cPickle.HIGHEST_PROTOCOL)
        upsname = datafile.replace(".pickle","-upsampled.dat")
        print "Saving to upsampled: ", self.D.shape
        self.D.tofile(upsname)
    def get_rid_of_roehrchen(self):
        """ remove measureing tube present in Lupine data """
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
