# vim:ts=4:sw=4:sts=4:et:ai
import numpy as np
from enthought.mayavi import mlab
from enthought.tvtk.util.ctf import ColorTransferFunction

def show_vectorfield(S, U,V,W):
    print "Show Vectors"
    mins = S.min()
    ptps = S.ptp()
    mlab.quiver3d(U,V,W, scalars=S,scale_mode="scalar",vmin=mins+0.2*ptps,vmax=mins+0.8*ptps)
    print "done."

def show_volume(D, cm="Spectral", minfact=0.1, maxfact=0.9,visible=True, normalize=True):
    print "Show volume"
    mind = D.min()
    D -= D.min()

    src = mlab.pipeline.scalar_field(D)
    ptpd = D.ptp()
    if normalize:
        R = (mind+minfact*ptpd, mind+maxfact*ptpd)
    else:
        R = (minfact-mind, maxfact-mind)


    v = mlab.pipeline.volume(src, vmin=R[0],vmax=R[1])

    if not (cm == "Spectral"):
        ctf = ColorTransferFunction()
        ctf.range = R
        ctf.add_rgb_point(mind, 1,1,1)
        ctf.add_rgb_point(R[0], 1,1,1)
        ctf.add_rgb_point(R[1], 0,0,0)
        ctf.add_rgb_point(mind+ptpd, 0,0,0)
        v._volume_property.set_color(ctf)
        v._ctf = ctf
        v.update_ctf = True

    from enthought.tvtk.util.ctf import PiecewiseFunction
    otf = PiecewiseFunction()
    otf.add_point(mind, 0)
    otf.add_point(R[0], 0)
    #otf.add_point(R[0]+0.1*ptpd, 0.1)
    #otf.add_point(R[0]+0.2*ptpd, 0.3)
    #otf.add_point(R[0]+0.3*ptpd, 0.5)
    otf.add_point(R[0]+0.2*ptpd, 0.7)
    otf.add_point(R[1], 1.0)
    otf.add_point(mind+ptpd, 1.0)
    v._otf = otf
    v._volume_property.set_scalar_opacity(otf)
    v.update_ctf = True

    
    print "done"

def show_iso(D,fact=0.2, cm="bone",opacity=0.5,visible=True,normalize=True):
    print "Show Iso"
    Dmin = D.min()
    D -= Dmin
    Dptp = D.ptp()
    src = mlab.pipeline.scalar_field(D)
    #mlab.contour3d(D)
    #mlab.pipeline.volume(src, vmin=D.min()+0.2*D.ptp(),vmax=D.max()-0.2*D.ptp())
    if fact.__class__ == float:
        if normalize: c = fact * Dptp
        else:         c = fact - Dmin
        print "Thresholding at ", c
        mlab.pipeline.iso_surface(src, contours=[c, ], opacity=opacity, colormap=cm)
    else:
        for i in fact:
            if normalize: c = i * Dptp
            else:         c = i - Dmin
            print "Thresholding at ", c
            mlab.pipeline.iso_surface(src, contours=[c,], opacity=opacity, colormap=cm)
    #mlab.pipeline.iso_surface(src, contours=[D.max()-0.2*D.ptp(), ], opacity=0.2)
    #from enthought.mayavi.modules.streamline import Streamline
    #mlab.add_module(s)
    print "done"

def show_points(fn,fne=None,cm="Blues",mode="2dtriangle",color=None):
    L = []
    S = []
    D = []
    print "Show Point3D"
    with open(fn) as f:
        for line in f.readlines():
            line = line.split()
            if len(line)<7:
                continue
            L.append([float(x) for x in line[:3]])
            S.append(float(line[3]))
            D.append([float(x) for x in line[4:]])
    S = np.array(S)
    L = np.vstack(L)
    D = np.vstack(D)
    #mlab.quiver3d(L[:,0],L[:,1],L[:,2], D[:,0],D[:,1],D[:,2], scale_factor=3.)

    if (L<0).sum() > 0:
        """ this is in cm """
        L[:,0] = (-L[:,0])/100.*256 + 120
        L[:,1] = (L[:,1])/100.*256  + 110
        L[:,2] = (L[:,2])/131.*256  + -0
        #L[:,0] = (-L[:,0])/100.*256 + 110
        #L[:,1] = (L[:,1])/100.*256  + 120
        #L[:,2] = (L[:,2])/131.*256  + -5
        gt = True
    else:
        L[:] += 0.5
        gt = False

    #mlab.points3d(L[:,0].squeeze(),L[:,1].squeeze(),L[:,2].squeeze(),S,scale_factor=1.5,colormap="bone",scale_mode="none",mode="2dtriangle")

    #pts = mlab.points3d(L[:,0],L[:,1],L[:,2],S,scale_factor=0.8,colormap=cm,scale_mode="none",resolution="20",mode=mode)
    #pts = mlab.points3d(L[:,0],L[:,1],L[:,2],S,scale_factor=0.1,colormap=cm,scale_mode="none",resolution="20",mode=mode)
    pts = mlab.points3d(L[:,0],L[:,1],L[:,2],S,scale_factor=10,colormap=cm,resolution="20",mode=mode,scale_mode='scalar')
    print "Done."
    if None==fne: return
    print "Show Edges3D"
    L = []
    with open(fne) as f:
        for line in f.readlines():
            vec = line.split()
            line = [int(x) for x in vec[:-1]]
            L.append(line)

        L = np.vstack(L)
        pts.mlab_source.dataset.lines = L
        tube = mlab.pipeline.tube(pts,tube_sides=8)
        tube.filter.radius_factor = 100.
        tube.filter.vary_radius = 'vary_radius_by_scalar'
        if color==None:
            color = (1,0,0) if gt else (0,1,0)
        #color = (0,1,0) if gt else (1,0,0)
        mlab.pipeline.surface(tube,color=color)
    print "Done."
