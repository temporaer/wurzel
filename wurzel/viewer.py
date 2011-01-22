# vim:ts=4:sw=4:sts=4:et:ai
from enthought.mayavi import mlab

def show_volume(D, cm="Spectral", minfact=0.1, maxfact=0.9):
    print "Show volume"
    src = mlab.pipeline.scalar_field(D)
    mind = D.min()
    ptpd = D.ptp()
    mlab.pipeline.volume(src, vmin=mind+minfact*ptpd,vmax=mind+maxfact*ptpd)
    print "done"

def show_iso(D,fact=0.2):
    print "Show Iso"
    src = mlab.pipeline.scalar_field(D)
    #mlab.contour3d(D)
    #mlab.pipeline.volume(src, vmin=D.min()+0.2*D.ptp(),vmax=D.max()-0.2*D.ptp())
    mlab.pipeline.iso_surface(src, contours=[D.min()+fact*D.ptp(), ], opacity=0.5, colormap="bone")
    #mlab.pipeline.iso_surface(src, contours=[D.max()-0.2*D.ptp(), ], opacity=0.2)
    print "done"
