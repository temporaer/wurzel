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
    eig = img3dops.get_ev_of_hessian(d.D)
    S = linestructure.get_curve_3D(eig,0.25,0.5,0.5)
    #S *= sigma**2 # as in Sato et al
    #print "Saving S"
    #np.save("data/S-%02.01d.npy"%sigma, S)
    #comm.Send(S, dest=0, tag=77)
    res = {}
    res["sato"] = S.astype("float32")
    res["ev10"] = eig["ev10"]
    res["ev11"] = eig["ev11"]
    res["ev12"] = eig["ev12"]
    return res

def loaddata(fn,slave=True):
    print "Loading dataset `%s'"%fn
    D = dataset("%s.dat"%fn,crop=False,upsample="zoom",remove_rohr=True,usepickled=slave)
    return D

if __name__ == "__main__":
    try:
        filename = sys.argv[1]
    except:
        print "usage: ", sys.argv[0], "[filename]"
        sys.exit()
    basename = filename
    D = loaddata(basename,slave=True)

    s = 1.20
    sigma0 = 1.0
    sigmas = []
    sigmas.extend([sigma0 * s**i for i in xrange(0,6)])
    print "Sigmas: ", sigmas

    use_ip = True
    if use_ip:
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
    cmdl = ["res = [c3d(D,s) for s in sigmas]",
            "sato = [x['sato'] for x in res]",
            "ev10 = [x['ev10'] for x in res]",
            "ev11 = [x['ev11'] for x in res]",
            "ev12 = [x['ev12'] for x in res]",]

    print "Executing..."
    if use_ip:
        for c in cmdl:
            print mec.execute(c)
        print "Gathering..."
        mec.block=True
        sato = mec.gather("sato")
        ev10 = mec.gather("ev10")
        ev11 = mec.gather("ev11")
        ev12 = mec.gather("ev12")
    else:
        D = loaddata(basename)
        for c in cmdl:
            print "Executing ", c
            c = compile(c, '<string>', 'exec')
            eval(c,globals(), locals())

    print "Finding maxima"
    #import pdb; pdb.set_trace()
    for s in xrange(1,len(sigmas)):
        arg = sato[0] < sato[s]
        sato[0][arg] = sato[s][arg]
        ev10[0][arg] = ev10[s][arg]
        ev11[0][arg] = ev11[s][arg]
        ev12[0][arg] = ev12[s][arg]
    x = sato[0]
    x -= x.min();
    x /= x.max()
    print "Saving to ", basename+".sato", "range:", x.min(),x.max()
    print x.dtype
    print x.flags
    np.save("res.npy", x)
    x.tofile(basename+".sato")
    ev10[0].tofile(basename+".ev10")
    ev11[0].tofile(basename+".ev11")
    ev12[0].tofile(basename+".ev12")
    
    #os.system("scp res.npy l_schulz@gitta:checkout/git/wurzel/data/neu.npy")
    #os.execl("/usr/bin/ssh", "ssh", "-X", "l_schulz@gitta","cd checkout/git/wurzel; ipython -wthread vis.py")
    print "Done"
