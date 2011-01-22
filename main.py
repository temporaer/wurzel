# vim:ts=4:sw=4:sts=4:et:ai
from wurzel import linestructure
from wurzel import img3dops
from wurzel.dataset import dataset
import wurzel.viewer as viewer

if __name__ == "__main__":
    d = dataset("data/L2_22aug.dat",crop=True)
    l1,l2,l3 = img3dops.get_ev_of_hessian(d.D)
    S = linestructure.get_curve_3D(l1,l2,l3)
    viewer.show_iso(S)
