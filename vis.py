# vim:ts=4:sw=4:sts=4:et:ai

# Copyright 2011 University of Bonn
# Author: Hannes Schulz

import sys, os
import numpy as np
import wurzel.viewer as viewer
from wurzel.dataset import dataset, WurzelInfo
#from enthought import mayavi 
from enthought.mayavi import mlab
#from scipy.ndimage.morphology import grey_closing
#from scipy.ndimage.measurements import label
#import IPython

#@mlab.animate(delay=500, ui=False)
def mkimg(fig,name):
    scene = fig
    scene.scene.camera.position = [-296.99079988325883, 407.90377529750242, -1481.1370972660488]
    scene.scene.camera.focal_point = [128.50921821594238, 421.38778924942017, 127.63644504547119]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [-0.0014367787462688198, 0.99996695718867323, -0.0080012622541963223]
    scene.scene.camera.clipping_range = [1391.5929806947984, 2011.2498060590401]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.render()
    mlab.savefig(os.path.join("plots", basename + "--" + name+"-1.png"),magnification=2)

    scene.scene.camera.position = [111.58714792682832, 509.05217371314114, -40.52978953059533]
    scene.scene.camera.focal_point = [10.281192710643207, 445.46754443016596, 124.36301423212157]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [0.084694202542404673, 0.91110955015930073, 0.40336866470292765]
    scene.scene.camera.clipping_range = [0.73702606876179377, 737.02606876179379]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.render()
    mlab.savefig(os.path.join("plots", basename + "--" + name+"-0.png"),magnification=2)

    scene.scene.camera.position = [23.11240194582102, 682.93747207629042, -317.94910456889863]
    scene.scene.camera.focal_point = [139.60023362585179, 686.62894911692354, 122.47988165121504]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [-0.0014367787462688198, 0.99996695718867323, -0.0080012622541963223]
    scene.scene.camera.clipping_range = [195.12015622774646, 784.56302137817386]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.render()
    mlab.savefig(os.path.join("plots", basename + "--" + name+"-2.png"),magnification=2)



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

  info = WurzelInfo(basename+".dat")
  datapath = info.datapath

  if offscreen:
    mlab.options.offscreen = True
  #fig = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(2048,1792))
  fig = mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(500,500))
  fig.scene.anti_aliasing_frames=1
  mlab.clf()
  token = token.split("-")
  if "raw" in token:
      Draw  = np.load(os.path.join(datapath, basename, "upsampled.pickle"))
      viewer.show_iso(Draw, 4 , "bone", 0.15)
      if offscreen:
        mkimg(fig, "raw")
        mlab.clf()

  if "raw_before" in token:
      Draw  = np.fromfile(os.path.join(datapath, basename + ".dat"),dtype="float32").reshape(256,832,256)
      viewer.show_iso(Draw, 4 , "bone", 0.15)
      if offscreen:
        mkimg(fig, "raw")
        mlab.clf()


  if "rawvol" in token:
      import re
      snr   = re.search(r'snr(\d+)', basename).group(1)
      Draw  = np.load(os.path.join(datapath, basename, "upsampled.pickle"))
      #viewer.show_volume(Draw, "bone", 0.004, 0.1)
      viewer.show_volume(Draw, "bone", 0.004*500./float(snr), 0.1)
      if offscreen:
        mkimg(fig, "rawvol")
        mlab.clf()

  if "mass" in token:
      #viewer.show_points( basename+"/vertices.txt", basename+"/edges.txt", dscale=0.5,what="wireframe")
      #viewer.show_points( basename+"/vertices.txt", basename+"/edges.txt", dscale=1,what=3)
      viewer.show_points( basename+"/vertices.txt", basename+"/edges.txt", dscale=4,what=3, color=(.1,.1,.8),opacity=0.3)
      #viewer.show_iso(Draw, 0.015 , "bone", 0.15)   # need 0.2 to get rid of noise
      if offscreen:
        mkimg(fig, "mass")
        mlab.clf()
  if "mass_gt" in token:
      gt_name = basename.split("roots")[0]
      viewer.show_points( gt_name+"/vertices.txt", gt_name+"/edges.txt", dscale=4,what=3,color=(.8,.3,.3), opacity=0.3)
      if offscreen:
        mkimg(fig, "mass_gt")
        mlab.clf()
  if "mass_cmp" in token:
      viewer.show_points( basename+"/vertices.txt", basename+"/edges.txt", dscale=4,what=3, color=(.1,.1,.8),opacity=0.3)
      gt_name = basename.split("roots")[0]
      viewer.show_points( gt_name+"/vertices.txt", gt_name+"/edges.txt", dscale=4,what=3,color=(.8,.3,.3), opacity=0.3)
      if offscreen:
        mkimg(fig, "mass_cmp")
        mlab.clf()

  if "scales" in token:
      info = WurzelInfo(basename+".scales")
      Dscales = np.fromfile(basename+".scales", dtype="float32").reshape(info.shape)
      if offscreen:
        mkimg(fig, "scales")
        mlab.clf()

  if "diameter" in token:
      viewer.show_points( basename+"/vertices.txt", basename+"/edges.txt", dscale=1,what=4)
      if offscreen:
        mkimg(fig, "diameter")
        mlab.clf()

  if "diameter_cmp" in token:
      viewer.show_points( basename+"/vertices.txt", basename+"/edges.txt", dscale=1,what=4, color=(.1,.1,.8),opacity=0.3)
      gt_name = basename.split("roots")[0]
      viewer.show_points( gt_name+"/vertices.txt", gt_name+"/edges.txt", dscale=1,what=4,color=(.8,.3,.3), opacity=0.3)
      if offscreen:
        mkimg(fig, "diameter_cmp")
        mlab.clf()

  if "wireframe_cmp" in token:
      viewer.show_points( basename+"/vertices.txt", basename+"/edges.txt", dscale=0.5,what="wireframe")
      L = basename.split("roots")
      viewer.show_points( L[0]+"/vertices.txt", L[0]+"/edges.txt", dscale=0.5,what="wireframe",color=(1,0,0))
      if offscreen:
        mkimg(fig, "wireframe_cmp")
        mlab.clf()

  if "wireframe_gt" in token:
      gt_name = basename.split("roots")[0]
      viewer.show_points( gt_name+"/vertices.txt", gt_name+"/edges.txt", dscale=0.5,what="wireframe", color=(.8,.3,.3))
      if offscreen:
          mkimg(fig, "wireframe_gt")
          mlab.clf()

  if "wireframe" in token:
      viewer.show_points( basename+"/vertices.txt", basename+"/edges.txt", dscale=0.5,what="wireframe", color=(.1,.1,.8))
      #viewer.show_iso(Dsato, 0.002 , "bone", 0.1)
      #Draw  = np.load(os.path.join(datapath, basename, "upsampled.pickle"))
      #viewer.show_iso(Draw, 10, "bone", 0.40)
      #viewer.show_volume(Draw, "bone", 0.01, 0.1)  # 0.03, 0.2 works well
      if offscreen:
        mkimg(fig, "wireframe")
        mlab.clf()

  if "wireframe_diff" in token:
      gt_name = basename.split("roots")[0]
      viewer.show_points( gt_name+"/vertices.txt", gt_name+"/edges.txt", dscale=0.5, what="wireframe", opacity=0.2,color=(0.8,0.3,0.3))
      viewer.show_points( basename+"/vertices.txt", basename+"/edges.txt", dscale=0.5*.85, what="wireframe", opacity=0.2,color=(0.1,0.1,0.8))
      viewer.show_points_onefile( basename+"/toomuch.txt", dscale=0.5, what="wireframe", color=(.1,.1,1),dumb=True)
      viewer.show_points_onefile( basename+"/missing.txt", dscale=0.5, what="wireframe", color=(1,.3,.3),dumb=True, opacity=0.3)
      if offscreen:
        mkimg(fig, "wireframe_diff")
        mlab.clf()
      

  if "satoiso" in token:
      Dsato = dataset(basename+"/sato.dat", upsample=None, crop=False,usepickled=False,medianfilt=False).D
      viewer.show_iso(Dsato, 4, "Spectral", 0.2)
      if offscreen:
        mkimg(fig, "satoiso")
        mlab.clf()

  if "satovol" in token:
      Dsato  = dataset(basename+"/sato.dat")
      #viewer.show_volume(Dsato, "bone", 0.00001, 0.01) # barley
      viewer.show_volume(Dsato, "bone", 0.01, 0.04) # lupine
      if offscreen:
        mkimg(fig, "satovol")
        mlab.clf()

  if "dmap" in token:
      info = WurzelInfo(basename+".dat")
      dmap  = np.fromfile(os.path.join(datapath, basename, "d_map.dat"),dtype="float64").reshape(info.shape)
      #import pdb; pdb.set_trace()
      viewer.show_iso(dmap, .003 , "bone", 0.15, absolutethresh=True)
      if offscreen:
        mkimg(fig, "dmap")
        mlab.clf()

#import pdb; pdb.set_trace()
if not offscreen:
    import IPython; IPython.frontend.terminal.embed.embed()
else:
    mlab.show()

