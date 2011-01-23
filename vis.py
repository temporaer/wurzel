# vim:ts=4:sw=4:sts=4:et:ai
import numpy as np
import wurzel.viewer as viewer
from wurzel.dataset import dataset
from enthought.mayavi import mlab
from scipy.ndimage.morphology import grey_closing
from scipy.ndimage.measurements import label

if __name__ == "__main__":
  mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
  mlab.clf()
  #D0 = 0.5**2 * np.load("data/S-0.npy").astype("float32")
  #D1 = 1.5**2 * np.load("data/S-1.npy").astype("float32")
  #D2 = 2.5**2 * np.load("data/S-2.npy").astype("float32")
  #D3 = 3.5**2 * np.load("data/S-4.npy").astype("float32")
  #D = reduce(np.maximum,[D0 , D1 , D2 , D3])

  Do = dataset("data/L2_22aug.dat", upsample="zoom", crop=False,usepickled=True).D
  Dr = np.load("data/res.npy").astype("float32")
  #dmin = Dr.min()
  #dptp = Dr.ptp()
  #minv = dmin+0.15*dptp
  #Dr[Dr<minv]=0
  #L, nf = label(Dr)
  #print "Number of connected components: ", nf
  viewer.show_iso(Do, 0.20 , "bone", 1.0)
  viewer.show_iso(Dr, 0.07 , "Spectral", 0.2)

