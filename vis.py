# vim:ts=4:sw=4:sts=4:et:ai

# Copyright 2011 University of Bonn
# Author: Hannes Schulz

import sys
import numpy as np
import wurzel.viewer as viewer
from wurzel.dataset import dataset, WurzelInfo
#from enthought import mayavi 
from enthought.mayavi import mlab
#from scipy.ndimage.morphology import grey_closing
#from scipy.ndimage.measurements import label

#@mlab.animate(delay=500, ui=False)
def mkimg(fig,name):
    scene = fig.scene
    #print "Rendering..."
    #scene.camera.position = [463.9984860311863, 1809.2988117125938, 109.73671134438729]
    #scene.camera.focal_point = [421.51009012013674, 95.419905483722687, 94.895551919937134]
    #scene.camera.view_angle = 30.0
    #scene.camera.view_up = [0.99969286540503455, -0.02478022775248526, -0.00033936824643516789]
    #scene.camera.clipping_range = [1482.2672553487132, 2010.9001920572837]
    #scene.camera.compute_view_plane_normal()
    #scene.render()
    #mlab.savefig(name+"-plain.png",magnification=2)
    scene.camera.position = [421.63833162281662, 0.41217041015625, -1602.0273769039982]
    scene.camera.focal_point = [421.63833162281662, 0.41217041015625, 0.0]
    scene.camera.view_angle = 30.0
    scene.camera.view_up = [0.99991970691306742, -0.012672005637842045, 0.0]
    scene.camera.clipping_range = [1574.0671031349582, 1641.0877875575584]
    scene.camera.compute_view_plane_normal()
    scene.render()
    mlab.savefig(name+"-polar.png",magnification=2)


if __name__ == "__main__":
  offscreen = False
  dump      = False
  basename  = None
  token     = ""
  for idx, s in enumerate(sys.argv[1:]):
    if   s == "-o": offscreen = True
    elif s == "-d": dump      = True
    elif s == "-t": token     = sys.argv[idx+2]
    else: basename = s
  if not basename:
      print "Usage: ",sys.argv[0],"[-o] [-d] [-t token] {basename}"
      sys.exit(1)

  if offscreen:
    mlab.options.offscreen = True
  fig = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(2048,1792))
  #fig = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(400,300))
  fig.scene.anti_aliasing_frames=1
  mlab.clf()

  #ev10  = dataset("data/L2_22aug.ev10", dz=256, dtype="float32", upsample=None, crop=False,usepickled=False,medianfilt=False).D
  #ev11  = dataset("data/L2_22aug.ev11", dz=256, dtype="float32", upsample=None, crop=False,usepickled=False,medianfilt=False).D
  #ev12  = dataset("data/L2_22aug.ev12", dz=256, dtype="float32", upsample=None, crop=False,usepickled=False,medianfilt=False).D

  #Draw  = dataset(basename+".dat", upsample="zoom", crop=False,usepickled=True, remove_rohr=True).D
  #Ddist1 = dataset(basename+"-paths1.dat", dz=256,upsample=None, crop=False,usepickled=False,medianfilt=False).D.swapaxes(0,2)
  #Ddist = dataset(basename+"-paths.dat", dz=256,upsample=None, crop=False,usepickled=False,medianfilt=False).D.swapaxes(0,2)
  #Ddist = dataset(basename+"-d_map.dat", dz=256, dtype="float64", upsample=None, crop=False,usepickled=False,medianfilt=False).D.swapaxes(0,2)
  Dsato = dataset(basename+".sato", upsample=None, crop=False,usepickled=False,medianfilt=False).D
  Draw  = dataset(basename+"-upsampled.dat", upsample=None, crop=False,usepickled=False,medianfilt=False).D
  #Draw  = np.fromfile("/home/local/cuv/build/nlmeanresult.dat", dtype="float32").reshape(410,192,192)
  #dmin = Dsato.min()
  #dptp = Dsato.ptp()
  #minv = dmin+0.15*dptp
  #Dsato[Dsato<minv]=0
  #L, nf = label(Dsato)
  #print "Number of connected components: ", nf

  #viewer.show_points(basename+"-ranks.txt",cm="Spectral", mode="sphere")
  #viewer.show_points(basename+"-vertices.txt", basename+"-edges.txt")
  #viewer.show_points( basename+"-vertices.txt", basename+"-edges.txt")
  #viewer.show_iso(255.-Ddist, 1./26., "jet", 0.7)  
  #viewer.show_iso(255-Ddist, [1/26.,2/26.], "RdGy", 0.2)
  #viewer.show_iso(255-Ddist, 250., "Spectral", 0.7)
  #viewer.show_volume(Dsato, "Spectral", 0.2, 0.40)
  #if offscreen:
  #  mkimg(fig, "path-ours")
  #  #mlab.clf()

  #viewer.show_points( "data/GersteLA_192x192x410_normal-vertices.txt", "data/GersteLA_192x192x410_normal-edges.txt")
  #viewer.show_points( "data/GersteLA_128x128x410-vertices.txt", "data/GersteLA_128x128x410-edges.txt", color=(0,0,1))
  if token == "raw":
      viewer.show_iso(Draw, 0.015 , "bone", 0.15)
      if offscreen:
        mkimg(fig, "raw")
        mlab.clf()
  if token == "rawvol":
      viewer.show_volume(Draw, "bone", 0.01, 0.1)
      if offscreen:
        mkimg(fig, "rawvol")
        mlab.clf()

  if token == "mass":
      viewer.show_points( basename+"-vertices.txt", basename+"-edges.txt", dscale=1,what=3)
      #viewer.show_iso(Draw, 0.015 , "bone", 0.15)   # need 0.2 to get rid of noise
      if offscreen:
        mkimg(fig, "mass")
        mlab.clf()

  if token == "scales":
      info = WurzelInfo(basename+".scales")
      Dscales = np.fromfile(basename+".scales", dtype="float32").reshape(info.shape)
      #viewer.show_iso(Dscales, 0.3, "bone", 0.15)
      arg0=Dsato<0.01
      Dscales[arg0] = 0
      viewer.show_volume(Dscales, "bone", 0.01, 0.8)  # 0.03, 0.2 works well
      if offscreen:
        mkimg(fig, "scales")
        mlab.clf()

  if token == "diameter":
      viewer.show_points( basename+"-vertices.txt", basename+"-edges.txt", dscale=1,what=4)
      viewer.show_iso(Draw, 0.150 , "bone", 0.40)   # lupine
      if offscreen:
        mkimg(fig, "diameter")
        mlab.clf()

  if token == "wireframe":
      viewer.show_points( basename+"-vertices.txt", basename+"-edges.txt", dscale=0.5,what="wireframe")
      #viewer.show_iso(Dsato, 0.002 , "bone", 0.1)
      #viewer.show_iso(Draw, 0.015 , "bone", 0.15)   # barley
      viewer.show_iso(Draw, 0.150 , "bone", 0.40)   # lupine
      #viewer.show_volume(Draw, "bone", 0.01, 0.1)  # 0.03, 0.2 works well
      if offscreen:
        mkimg(fig, "wireframe")
        mlab.clf()

  if token == "satoiso":
      viewer.show_iso(Dsato, 0.002 , "Spectral", 0.2)
      if offscreen:
        mkimg(fig, "satoiso")
        mlab.clf()

  if token == "satovol":
      #viewer.show_volume(Dsato, "bone", 0.00001, 0.01) # barley
      viewer.show_volume(Dsato, "bone", 0.01, 0.04) # lupine
      if offscreen:
        mkimg(fig, "satovol")
        mlab.clf()

  #DReis = np.load("data/reispflanze_wurzeln-laser.npy")
  #viewer.show_laserpoints(DReis);

  #DReis = np.loadtxt("laser/rec_del.txt")
  #viewer.show_laserpoints(DReis,"bone",ss=1,color=(1,0,0));
  #DReis = np.loadtxt("laser/rec_out.txt")
  #viewer.show_laserpoints(DReis,"bone",ss=1,color=(0,0,1));
  #DReis = np.loadtxt("laser/rec.txt")
  #viewer.show_laserpoints(DReis,"bone",ss=1,color=(0,1,0));
  #DReis = np.loadtxt("laser/tube.txt")
  #viewer.show_laserpoints(DReis,"Spectral",color=(0,0,1));

  #viewer.show_iso(Draw, 0.19 , "bone", 0.1)   # need 0.2 to get rid of noise
  #if offscreen:
  #  mkimg(fig, "us-vs-ground")
  #  mlab.clf()

  #viewer.show_iso(Ddist1, 0.05 , "RdGy", 1.0)
  #viewer.show_iso(Draw, 0.19 , "bone", 0.17)   # need 0.2 to get rid of noise
  #if offscreen:
  #  mkimg(fig, "rendered-paths1")
  #  mlab.clf()

  #viewer.show_iso(Ddist, 0.05 , "RdGy", 1.0)
  #viewer.show_iso(Draw, 0.19 , "bone", 0.17)   # need 0.2 to get rid of noise
  #if offscreen:
  #  mkimg(fig, "rendered-paths")
  #  mlab.clf()

  #viewer.show_volume(Draw, "bone", 0.15,1.0)   # 0.13, 0.3 works well
  #if offscreen:
  #  mkimg(fig, "raw")
  #  mlab.clf()

  #viewer.show_iso(Dsato, 0.002 , "Spectral", 0.2)  
  #if offscreen:
  #   mkimg(fig, "sato-iso")
  #   mlab.clf()

  #viewer.show_volume(Dsato, "bone", 0.005, 0.2)  # 0.03, 0.2 works well
  #if offscreen:
  #  mkimg(fig, "sato-bone")
  #  mlab.clf()


