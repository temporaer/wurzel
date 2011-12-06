# vim:ts=4:sw=4:sts=4:et:ai

# Copyright 2011 University of Bonn
# Author: Hannes Schulz


import gc, os, sys
from wurzel import linestructure
from wurzel import img3dops
from wurzel.dataset import dataset
import numpy as np


#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#root = MPI.ROOT


def vesselness(d, sigma):
    """ calculate vesselness of dataset d, after smoothing with gaussian of width sigma """
    print "Vesselness with sigma", sigma
    d           = d.get_smoothed(sigma)
    gc.collect()
    eig         = img3dops.get_ev_of_hessian(d.D,sigma,gamma   = 0.8)
    S           = linestructure.get_curve_3D(eig,0.25,0.5,0.5)
    res         = {}
    res["sato"] = S.astype("float32")
    res["ev10"] = eig["ev10"]
    res["ev11"] = eig["ev11"]
    res["ev12"] = eig["ev12"]
    return res

def loaddata(fn,slave=True):
    """" load a dataset from filename fn, use pickled data if slave=True. 
        for each dataset, this has to be called at least once with slave=False 
        to create the pickled version.
    """
    print "Loading dataset `%s'"%fn
    D = dataset("%s.dat"%fn,crop=False,upsample="zoom",remove_rohr=True,usepickled=slave)
    return D

def parallel_run(args):
        basename, sigma = args
        import main
        D = main.loaddata(basename)
        res = main.vesselness(D,sigma)
        sato = res['sato']
        return sato

if __name__ == "__main__":
    parallelize = False
    try:
        basename = sys.argv[1]
        if len(sys.argv)>2:
            parallelize   = sys.argv[2]=="-p"
    except:
        print "usage: ", sys.argv[0], "[basename] [-p]"
        print "basename     is the base name w/o extension of the raw MRI data"
        print "-p           if this is given, parallelize over multiple computers using a running ipcluster"
        sys.exit()

    # when parallelizing over multiple computers, assure here that data is upsampled /exactly/ once!
    D = loaddata(basename,slave=[True,False][parallelize])
    #D = loaddata(basename,slave=[True,True][parallelize])

    # ##########################################
    #  Determine set of scales
    # ##########################################
    # Barley
    num_scales = 20
    minimum_radius_mm = 0.08 / 2
    maximum_radius_mm = 2.5  / 2
    sigma0 = minimum_radius_mm
    s      = (maximum_radius_mm/minimum_radius_mm)**(1./num_scales)
    sigmas = []
    sigmas.extend([sigma0 * s**i for i in xrange(0,num_scales)])
    print "Sigmas (mm) : ", sigmas
    sigmas = [x / D.info.scale for x in sigmas] # convert from mm to vox
    print "Sigmas (vox): ", sigmas

    # ##########################################
    #  Prepare remote systems (load data, ...)
    # ##########################################
    if parallelize:
        from IPython.parallel import Client
        rc = Client(profile="wurzel_cluster")
        #rc = Client(profile_dir="/home/VI/staff/schulz/.ipython/profile_default")
        rc[:].execute("import os; os.chdir(\"%s\")" % os.getcwd())
        lview = rc.load_balanced_view()
        lview.block = True
        



    # ##########################################
    #  Execute cmdl remotely or locally
    # ##########################################
    print "Executing..."
    if parallelize:
        print "in parallel..."
        sato = lview.map(parallel_run, [(basename,s) for s in sigmas] )
        print "done."
    else:
        print "Sequential, and finding maxima"
        D = loaddata(basename)
        sato0 = None
        for s in sigmas:
            res  = vesselness(D,s)
            sato = res['sato']
            if sato0==None:
                sato0  = sato.astype("float32")
                scales = np.ones(sato.shape)*s
                continue
            arg         = sato0 < sato
            sato0[arg]  = sato[arg]
            scales[arg] = s
            

    sato0 -= sato0.min();
    sato0 /= sato0.max()
    basename = os.path.join(D.info.datapath, basename)
    print "Saving to ", os.path.join(basename,"sato.dat"), "data range:", sato0.min(),sato0.max()
    sato0.tofile(os.path.join(basename,"sato.dat"))
    scales.tofile(os.path.join(basename,"scales.dat"))

    print "Done"
