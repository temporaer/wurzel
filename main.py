# vim:ts=4:sw=4:sts=4:et:ai
import gc, os
from wurzel import linestructure
from wurzel import img3dops
from wurzel.dataset import dataset
import numpy as np

from IPython.kernel import client

#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#root = MPI.ROOT


def c3d(d, sigma):
    d = d.get_smoothed(sigma)
    gc.collect()
    l1,l2,l3 = img3dops.get_ev_of_hessian(d.D)
    S = linestructure.get_curve_3D(l1,l2,l3).astype("float32")
    S *= sigma
    #comm.Send(S, dest=0, tag=77)
    return S

def loaddata():
    print "Loading dataset"
    D = dataset("data/L2_22aug.dat",crop=False,upsample="zoom",usepickled=True)
    return D

if __name__ == "__main__":

    mec = client.MultiEngineClient()
    mec.activate()
    print "running on ids ", mec.get_ids()

    print "Importing on remote systems"
    mec.block = True
    print mec.execute("import os")
    print mec.execute("os.chdir(\"%s\")"%os.path.realpath("."))
    print mec.execute("from main import *")
    print mec.execute("D = loaddata()")
    #print mec.execute("""D = dataset("data/L2_22aug.dat",crop=False,upsample="zoom",usepickled=True)""")
    print mec.execute("print D")

    print "Scattering sigmas"
    mec.scatter('sigmas',   [0.5, 1.5, 2.5, 3.5, 4.5])
    print "Executing..."
    mec.execute("res = [c3d(D,s) for s in sigmas]")
    print "Gathering..."
    #print comm.Recv(D,1,tag=77)
    #comm.Gather()
    res = mec.gather("res")
    x = reduce(np.maximum, res).astype("float32")
    np.save("res.npy", x)
    print "Done"
