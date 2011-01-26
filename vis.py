# vim:ts=4:sw=4:sts=4:et:ai
import sys
import numpy as np
import wurzel.viewer as viewer
from wurzel.dataset import dataset
#from enthought import mayavi 
from enthought.mayavi import mlab
#from scipy.ndimage.morphology import grey_closing
#from scipy.ndimage.measurements import label

#@mlab.animate(delay=500, ui=False)
def mkimg(fig,name):
    scene = fig.scene
    #print "Rendering..."
    #scene.camera.position = [340.40906071071839, 40.846184594957535, 187.06626070910627]
    #scene.camera.focal_point = [185.219910913566, 123.09733843648607, 108.41372549551366]
    #scene.camera.view_angle = 30.0
    #scene.camera.view_up = [-0.42032768067096582, 0.071883920708805144, 0.9045205043586888]
    #scene.camera.clipping_range = [0.75260136121495613, 752.60136121495611]
    #scene.camera.compute_view_plane_normal()
    #scene.render()
    #print "Rendering..."
    #mlab.savefig(name+"-1.png",magnification=2)
    #scene.camera.position = [83.953707876827011, 311.92359815432201, 125.50936268716015]
    #scene.camera.focal_point = [89.856495599231664, 153.7986247717472, 229.42062081848633]
    #scene.camera.view_angle = 30.0
    #scene.camera.view_up = [-0.8681739651971645, 0.24949546570724399, 0.42898249235296071]
    #scene.camera.clipping_range = [0.59227624649847499, 592.27624649847496]
    #scene.camera.compute_view_plane_normal()
    #scene.render()
    #mlab.savefig(name+"-2.png",magnification=2)
    scene.camera.position = [-282.14857380964736, 128.5, 128.5]
    scene.camera.focal_point = [128.5, 128.5, 128.5]
    scene.camera.view_angle = 30.0
    scene.camera.view_up = [0.0, 0.0, 1.0]
    scene.camera.clipping_range = [152.81708807155087, 736.19580241679205]
    scene.camera.compute_view_plane_normal()
    scene.render()
    mlab.savefig(name+"-3.png", magnification=3)
    scene.camera.position = [128.5, 591.45266609081978, 128.5]
    scene.camera.focal_point = [128.5, 128.5, 128.5]
    scene.camera.view_angle = 30.0
    scene.camera.view_up = [-0.014778103364283094, 0.0, 0.99989079786792456]
    scene.camera.clipping_range = [204.59813942991155, 789.28445608218215]
    scene.camera.compute_view_plane_normal()
    scene.render()
    mlab.savefig(name+"-4.png", magnification=3)

def dump_uint8(D, vmin, vmax):
    dmin = D.min()
    dptp = D.ptp()
    vmin = dmin + vmin * dptp
    vmax = dmin + vmax * dptp

    D[D<vmin] = vmin
    D[D>vmax] = vmax
    D        -= D.min()
    D        *= 255.0 / D.max()
    D.astype("uint8").tofile("dump.dat")

if __name__ == "__main__":
  offscreen = False
  dump      = False
  for s in sys.argv[1:]:
    if s == "-o": offscreen = True
    if s == "-d": dump     = True
  if offscreen:
    mlab.options.offscreen = True
  fig = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
  fig.scene.anti_aliasing_frames=0
  mlab.clf()

  Draw = dataset("data/L2_22aug.dat", upsample="zoom", crop=False,usepickled=True).D
  #Draw = dataset("dump.dat", upsample=None, crop=False,usepickled=False).D
  Dsato = np.load("data/neu.npy").astype("float32")
  if dump:
    dump_uint8(Dsato, 0.02, 0.2)
    sys.exit()
  #dmin = Dsato.min()
  #dptp = Dsato.ptp()
  #minv = dmin+0.15*dptp
  #Dsato[Dsato<minv]=0
  #L, nf = label(Dsato)
  #print "Number of connected components: ", nf

  viewer.show_iso(Draw, 0.20 , "bone", 1.0)   # need 0.2 to get rid of noise
  if offscreen:
    mkimg(fig, "raw")
    mlab.clf()

  #viewer.show_volume(Draw, "bone", 0.13,0.3)   # 0.13, 0.3 works well
  #if offscreen:
  #  mkimg(fig, "raw")
  #  mlab.clf()

  #from scipy.ndimage.interpolation import zoom
  #Dsato = zoom(Dsato,2)
  #viewer.show_iso(Dsato, 0.05 , "Spectral", 1.0)  
  #if offscreen:
  # mkimg(fig, "sato-iso")
  # mlab.clf()

  #viewer.show_volume(Dsato, "Spectral", 0.02, 0.2)  # 0.03, 0.2 works well
  #if offscreen:
  #  mkimg(fig, "sato")
  #  mlab.clf()


