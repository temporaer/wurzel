# vim:ts=4:sw=4:sts=4:et:ai
from enthought.mayavi import mlab
from enthought.tvtk.util.ctf import ColorTransferFunction

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
        ctf.add_rgb_point(R[0], 1,1,1)
        ctf.add_rgb_point(R[1], 0,0,0)
        v._volume_property.set_color(ctf)
        v._ctf = ctf
        v.update_ctf = True

    from enthought.tvtk.util.ctf import PiecewiseFunction
    otf = PiecewiseFunction()
    otf.add_point(mind, 0)
    otf.add_point(R[0], 0)
    otf.add_point(R[0]+0.1*ptpd, 0.1)
    otf.add_point(R[0]+0.2*ptpd, 0.3)
    otf.add_point(R[0]+0.3*ptpd, 0.5)
    otf.add_point(R[0]+0.4*ptpd, 0.7)
    otf.add_point(R[1], 1)
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
    mlab.pipeline.iso_surface(src, contours=[D.min()+fact*D.ptp(), ], opacity=opacity, colormap=cm)
    #mlab.pipeline.iso_surface(src, contours=[D.max()-0.2*D.ptp(), ], opacity=0.2)
    print "done"
