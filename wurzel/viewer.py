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

def xyz2pol(L):
    cy = np.min(L[:,1]) + L[:,1].ptp()/2
    cz = np.min(L[:,2]) + L[:,2].ptp()/2
    A = np.arctan2(L[:,2]-cz,L[:,1]-cy)
    M = np.vstack((L[:,0],A,L[:,0])).T
    M[:,1] *= L[:,0].ptp()/2/np.pi/2
    M[:,2] = 0
    return M

def xyz2pol_grid(L):
    cy = L.shape[1]/2
    cz = L.shape[2]/2
    G  = np.mgrid[0:L.shape[1],0:L.shape[2]]
    G0 = np.mgrid[0:L.shape[1],0:L.shape[2]]
    G[0] -= cy
    G[1] -= cz

    ares = 500

    A = np.arctan2(G[1],G[0])
    M = np.zeros((L.shape[0],ares))
    for i in xrange(L.shape[1]):
        for j in xrange(L.shape[2]):
            M[:,int(A[i,j]/np.pi/2.0*ares)] = L[:,i,j]





def show_points(fn,fne=None,cm="Blues",mode="2dtriangle",color=None,swap=True,scaled=True,dscale=1,what=3):
    L = []
    S = []
    D = []
    if scaled:
        from dataset import WurzelInfo
        info = WurzelInfo(fn)
        print "Scale: ", info.scale
    wireframe = False
    if what == "wireframe":
        wireframe = True
        what = 3
    print "Show Point3D `%s'" % fn
    with open(fn) as f:
        for line in f.readlines():
            line = line.split()
            #if len(line)<7:
            #    continue
            L.append([float(x) for x in line[:3]]) # 0..2: coordinate axes
            S.append(float(line[what]))           # 3: mass, 4: diameter
            D.append([float(x) for x in line[4:]]) # D: not used
    S = np.array(S)
    L = np.vstack(L)
    D = np.vstack(D)
    if wireframe:
        S[:] = 1
    print "NumPoints: ", L.shape[0]
    if scaled:
        L /= info.scale  # L is in mm, now it matches raw data again
    S *= dscale
    S[S<0.1]=0.8
    #mlab.quiver3d(L[:,0],L[:,1],L[:,2], D[:,0],D[:,1],D[:,2], scale_factor=3.)

    if (L<0).sum() > 0 and False:
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
    #L = xyz2pol(L)

    #mlab.points3d(L[:,0].squeeze(),L[:,1].squeeze(),L[:,2].squeeze(),S,scale_factor=1.5,colormap="bone",scale_mode="none",mode="2dtriangle")

    #pts = mlab.points3d(L[:,0],L[:,1],L[:,2],S,scale_factor=0.8,colormap=cm,scale_mode="none",resolution="20",mode=mode)
    #pts = mlab.points3d(L[:,0],L[:,1],L[:,2],S,scale_factor=0.1,colormap=cm,scale_mode="none",resolution="20",mode=mode)
    pts = mlab.points3d(L[:,0],L[:,1],L[:,2],S,scale_factor=.01,colormap=cm,resolution="20",mode=mode,scale_mode='scalar')
    print "Done."
    if None==fne: return
    print "Show Edges3D"
    E = []
    thresh = 100
    with open(fne) as f:
        for line in f.readlines():
            vec = line.split()
            line = [int(x) for x in vec[:2]]
            #if np.linalg.norm(L[line[0],:] - L[line[1],:]) > thresh:
            #    continue
            E.append(line)

        E = np.vstack(E)
        pts.mlab_source.dataset.lines = E
        tube = mlab.pipeline.tube(pts,tube_sides=7,tube_radius=1, name="root tubes")
        #tube.filter.radius_factor = 1.
        tube.filter.vary_radius = 'vary_radius_by_absolute_scalar'
        #tube.filter.vary_radius = 'vary_radius_by_scalar'
        tube.filter.capping = True
        #tube.filter.use_default_normal = True

        if color==None:
            color = (1,0,0) if gt else (0,1,0)
        #color = (0,1,0) if gt else (1,0,0)
        mlab.pipeline.surface(tube,color=color)
    print "Done."

def show_laserpoints(L,cm="Blues",mode="2dtriangle",color=(0,0,0), ss=8):
    print "Show LaserPoints"
    print L.shape
    mlab.points3d(L[::ss,0].squeeze(),L[::ss,1].squeeze(),L[::ss,2].squeeze(),scale_factor=0.1,colormap=cm,scale_mode="none",mode=mode,resolution=4,color=color)
    print "Done."
