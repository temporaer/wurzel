# vim:ts=4:sw=4:sts=4:et:ai
import numpy as np
from enthought.mayavi import mlab
from enthought.tvtk.util.ctf import ColorTransferFunction

def show_volume2(D, cm="Spectral", minfact=0.1, maxfact=0.9,visible=True):
    print "Show volume"
    D -= D.min()

    src = mlab.pipeline.scalar_field(D)
    mind = D.min()
    ptpd = D.ptp()
    R = (mind+minfact*ptpd, mind+maxfact*ptpd)


    v = mlab.pipeline.volume(src, vmin=mind+minfact*ptpd,vmax=mind+maxfact*ptpd)

    ptpR = R[1]-R[0]
    if not (cm == "Spectral") or True:
        ctf = ColorTransferFunction()
        ctf.range = R
        ctf.add_rgb_point(mind, 0,0,0)
        ctf.add_rgb_point(R[0], 0,0,0)
        ctf.add_rgb_point(R[0]+0.25*ptpR, 0.500,0.500,0.500)
        ctf.add_rgb_point(R[0]+0.50*ptpR, 0.750,0.750,0.750)
        ctf.add_rgb_point(R[0]+0.75*ptpR, 0.875,0.875,0.875)
        ctf.add_rgb_point(R[1]          , 0.0,0.0,0.0)
        ctf.add_rgb_point(mind+ptpd, 1,1,1)
        v._volume_property.set_color(ctf)
        v._ctf = ctf
        v.update_ctf = True

    from enthought.tvtk.util.ctf import PiecewiseFunction
    otf = PiecewiseFunction()
    otf.add_point(mind, 1.0)
    otf.add_point(R[0], 1.0)
    #otf.add_point(R[0]+0.25*ptpR, 0.5)
    #otf.add_point(R[0]+0.50*ptpR, 0.25)
    otf.add_point(R[0]+0.90*ptpR, 0.800)
    otf.add_point(R[0]+0.95*ptpR, 0.500)
    otf.add_point(R[1], 0.0)
    otf.add_point(mind+ptpd, 0.0)
    v._otf = otf
    v._volume_property.set_scalar_opacity(otf)
    v.update_ctf = True

    
    print "done"
def show_volume(D, cm="Spectral", minfact=0.1, maxfact=0.9,visible=True):
    print "Show volume"
    D -= D.min()

    src = mlab.pipeline.scalar_field(D)
    mind = D.min()
    ptpd = D.ptp()
    R = (mind+minfact*ptpd, mind+maxfact*ptpd)


    v = mlab.pipeline.volume(src, vmin=mind+minfact*ptpd,vmax=mind+maxfact*ptpd)

    if not (cm == "Spectral") or True:
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
    otf.add_point(R[0]+0.4*ptpd, 0.7)
    otf.add_point(R[1], 1.0)
    otf.add_point(mind+ptpd, 1.0)
    v._otf = otf
    v._volume_property.set_scalar_opacity(otf)
    v.update_ctf = True

    
    print "done"

def show_iso(D,fact=0.2, cm="bone",opacity=0.5,visible=True):
    print "Show Iso"
    D -= D.min()
    src = mlab.pipeline.scalar_field(D)
    #mlab.contour3d(D)
    #mlab.pipeline.volume(src, vmin=D.min()+0.2*D.ptp(),vmax=D.max()-0.2*D.ptp())
    if fact.__class__ == float:
        mlab.pipeline.iso_surface(src, contours=[D.min()+fact*D.ptp(), ], opacity=opacity, colormap=cm)
    else:
        for i in fact:
            mlab.pipeline.iso_surface(src, contours=[D.min()+i*D.ptp(), ], opacity=opacity, colormap=cm)
    #mlab.pipeline.iso_surface(src, contours=[D.max()-0.2*D.ptp(), ], opacity=0.2)
    #from enthought.mayavi.modules.streamline import Streamline
    #mlab.add_module(s)
    print "done"

def show_points(fn,fne=None,cm="Blues",mode="2dtriangle"):
    L = []
    S = []
    print "Show Point3D"
    with open(fn) as f:
        for line in f.readlines():
            line = [int(x) for x in line.split()]
            if len(line)!=4:
                continue
            L.append(line[:-1])
            S.append(line[-1])
    L = np.vstack(L)
    pts = mlab.points3d(L[:,0],L[:,1],L[:,2],S,scale_factor=0.8,colormap=cm,scale_mode="none",resolution="20",mode=mode)
    print "Done."
    if None==fne: return
    print "Show Edges3D"
    L = []
    with open(fne) as f:
        for line in f.readlines():
            line = [int(x) for x in line.split()]
            L.append(line)
        L = np.vstack(L)
        pts.mlab_source.dataset.lines = L
        tube = mlab.pipeline.tube(pts,tube_radius=0.1)
        #tube.filter.vary_radius = 'vary_radius_by_scalar'
        mlab.pipeline.surface(tube,color=(0.8,0.8,0.8))
    print "Done."
