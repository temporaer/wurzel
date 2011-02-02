# vim:ts=4:sw=4:sts=4:et:ai
import gc, os, sys
from wurzel import linestructure
from wurzel import img3dops
from wurzel.dataset import dataset
import numpy as np


#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#root = MPI.ROOT


def c3d(d, sigma):
    print "c3d with sigma", sigma
    d = d.get_smoothed(sigma)
    gc.collect()
    l1,l2,l3 = img3dops.get_ev_of_hessian(d.D)
    S = linestructure.get_curve_3D(l1,l2,l3,0.25,0.5,0.5)
    S *= sigma**2
    #print "Saving S"
    #np.save("data/S-%02.01d.npy"%sigma, S)
    #comm.Send(S, dest=0, tag=77)
    return S.astype("float32")

def loaddata(fn,slave=True):
    print "Loading dataset"
    D = dataset("%s.dat"%fn,crop=False,upsample="zoom",usepickled=slave)
    return D

if __name__ == "__main__":
    from IPython.kernel import client
    try:
        filename = sys.argv[1]
    except:
        print "usage: ", sys.argv[0], "[filename]"
        sys.exit()
    basename = filename[:-4] # cut off ".dat"
    #D = loaddata(basename,slave=False)

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
    s = 1.1
    sigma0 = 1.1
    sigmas = []
    sigmas.extend([sigma0 * s**i for i in xrange(1,4)])
    mec.scatter('sigmas',   sigmas)
    print mec.execute("print sigmas")
    print "Executing..."
    print mec.execute("res = [c3d(D,s) for s in sigmas]")

    print "Gathering..."
    mec.block=True
    res = mec.gather("res")
    #import pdb; pdb.set_trace()
    x = reduce(np.maximum, res).astype("float32")
    x -= x.min(); x /= x.max()
    np.save("res.npy", x)
    x.tofile(basename+".sato")
    #os.system("scp res.npy l_schulz@gitta:checkout/git/wurzel/data/neu.npy")
    #os.execl("/usr/bin/ssh", "ssh", "-X", "l_schulz@gitta","cd checkout/git/wurzel; ipython -wthread vis.py")
    print "Done"
