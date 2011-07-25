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
    eig         = img3dops.get_ev_of_hessian(d.D,sigma,gamma   = 1)
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

    # ##########################################
    #  Determine set of scales
    # ##########################################
    # Barley
    s = 1.20
    sigma0 = 1.0
    # Lupine
    s = 1.30
    sigma0 = 1.0
    sigmas = []
    sigmas.extend([sigma0 * s**i for i in xrange(0,10)])
    print "Sigmas: ", sigmas
    sigmas = [x * 0.52 / D.info.scale for x in sigmas]  # 0.52 is the scale of the dataset which the params above were adjusted for

    # ##########################################
    #  Prepare remote systems (load data, ...)
    # ##########################################
    if parallelize:
        from IPython.kernel import client
        mec = client.MultiEngineClient()
        mec.activate()
        mec.flush()
        print "running on ids ", mec.get_ids()

        print "Importing on remote systems"
        mec.block = True
        print mec.execute("import os")
        print mec.execute("os.system('rm -f /tmp/*.so')")
        print mec.execute("os.chdir(\"%s\")"%os.path.realpath("."))
        print mec.execute("from main import *")
        print mec.execute("from wurzel.dataset import dataset")
        print mec.execute("D = loaddata('%s')"%basename)
        #print mec.push(dict(D=D))

        print "Scattering sigmas"
        mec.scatter('sigmas',   sigmas)

    # ##########################################
    #  Determine what is to be done on each host
    # ##########################################
    cmdl = ["res = [vesselness(D,s) for s in sigmas]",
            "sato = [x['sato'] for x in res]",
            #"ev10 = [x['ev10'] for x in res]",
            #"ev11 = [x['ev11'] for x in res]",
            #"ev12 = [x['ev12'] for x in res]",
            ]

    # ##########################################
    #  Execute cmdl remotely or locally
    # ##########################################
    print "Executing..."
    if parallelize:
        for c in cmdl:
            print mec.execute(c)
        print "Gathering..."
        mec.block=True
        sato = mec.gather("sato")
        #ev10 = mec.gather("ev10")
        #ev11 = mec.gather("ev11")
        #ev12 = mec.gather("ev12")
    else:
        D = loaddata(basename)
        for c in cmdl:
            print "Executing ", c
            c = compile(c, '<string>', 'exec')
            eval(c,globals(), locals())

    # ##########################################
    #  Find max and arg-max of sato data
    # ##########################################
    print "Finding maxima"
    scales = np.zeros(sato[0].shape).astype("float32")
    scales[:] = sigmas[0]
    #import pdb; pdb.set_trace()
    for s in xrange(1,len(sigmas)):
        arg = sato[0] < sato[s]
        sato[0][arg] = sato[s][arg]
        scales[arg]  = sigmas[s]
        #ev10[0][arg] = ev10[s][arg]
        #ev11[0][arg] = ev11[s][arg]
        #ev12[0][arg] = ev12[s][arg]
    x = sato[0]
    x -= x.min();
    x /= x.max()
    basename = os.path.join(D.info.datapath, basename)
    print "Saving to ", basename+".sato", "data range:", x.min(),x.max()
    x.tofile(basename+".sato")
    scales.tofile(basename+".scales")
    #ev10[0].tofile(basename+".ev10")
    #ev11[0].tofile(basename+".ev11")
    #ev12[0].tofile(basename+".ev12")

    print "Done"
