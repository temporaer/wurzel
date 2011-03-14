import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
from wurzel.dataset import WurzelInfo
from mpl_toolkits.axes_grid1 import make_axes_locatable


info = WurzelInfo(sys.argv[1])
D = np.loadtxt(sys.argv[1] + "-edges.txt")

mass = D[:,3]
ang  = D[:,4]
Z    = D[:,5]
print "Z: ", min(Z), max(Z)
print "a: ", min(ang), max(ang)

# ignore stuff above 14th plane
b = Z > 20 * info.shape[0] / 410.0 * info.scale

Z    = Z[b]
mass = mass[b]
ang  = ang[b]

resx  = 20
resy  = 40
img  = np.ones((resy,resx))
imgn = np.ones((resy,resx))

min_mass = min(mass)
max_mass = max(mass)
min_ang  = min(ang)
max_ang  = max(ang)
min_Z    = min(Z)
max_Z    = max(Z)

assert min_mass > 0
assert np.abs(min_ang) < 0.01
assert np.abs(max_ang-1.57) < 0.1

#mass = ((mass-np.min(mass))/np.ptp(mass) * (resx-1)).astype("int")
ang  = ((ang -np.min(ang ))/np.ptp(ang ) * (resx-1)).astype("int")
Z    = ((Z   -np.min(Z   ))/np.ptp(Z   ) * (resy-1)).astype("int")

fig = plt.figure(figsize=(4,4))
fig.subplots_adjust(hspace=0.1,wspace=0.1,left=0.07,bottom=0.07,right=1,top=0.85)
if True:
	for i in xrange(Z.shape[0]):
		x = ang[i]
		y = Z[i]
		z = mass[i]
		img[y,x] += z

	ax  = fig.add_subplot(121)
	#cax = ax.matshow(img/imgn,cmap="binary",norm=LogNorm(vmin=min_mass,vmax=max_mass))
	cax = ax.matshow(img/imgn,cmap="binary")
	ax.set_xticks([x for x in np.arange(0,resx+.01,resx/2)])
	ax.set_yticks([y for y in np.arange(0,resy+.01,resy/10)])
	ax.set_xticklabels(["%1.2f"%x for x in np.arange(min_ang,max_ang+.01,(max_ang-min_ang)/2)])
	ax.set_yticklabels(["%3.0f"%y for y in np.arange(min_Z,max_Z+.01,(max_Z-min_Z)/10)])
	ax.set_xlabel('angle to vertical (rad)')
	ax.set_ylabel('depth (mm)')
	ax.set_title("Mass distribution")
	divider   = make_axes_locatable(ax)
        axMarg = divider.append_axes("bottom", 1.2, pad=0.1, sharex=ax)

	axMarg.plot(np.arange(min_ang,max_ang,(max_ang-min_ang)/img.shape[1]),img.sum(axis=0))


if True:
	for i in xrange(Z.shape[0]):
		x = ang[i]
		y = Z[i]
		z = mass[i]
		img[y,x] += z
		imgn[y,x]+= 1

	ax  = fig.add_subplot(122)
	#cax = ax.matshow(img/imgn,cmap="binary",norm=LogNorm(vmin=min_mass,vmax=max_mass))
	cax = ax.matshow(img/imgn,cmap="binary")
	ax.set_xticks([x for x in np.arange(0,resx+.01,resx/3)])
	#ax.set_yticks([y for y in np.arange(0,resy+.01,resy/10)])
	ax.set_xticklabels(["%1.2f"%x for x in np.arange(min_ang,max_ang+.01,(max_ang-min_ang)/3)])
	#ax.set_yticklabels(["%3.0f"%y for y in np.arange(min_Z,max_Z+.01,(max_Z-min_Z)/10)])
	ax.set_xlabel('angle to vertical (rad)')
	plt.setp(ax.get_yticklabels(), visible=False)
	#ax.set_ylabel('depth (mm)')
	ax.set_title("Expected Mass")

	#ax  = fig.add_subplot(224)
	#ax.plot(np.arange(min_ang,max_ang,(max_ang-min_ang)/img.shape[1]),img.sum(axis=0))
	#ax.set_aspect(2)
	#ax.apply_aspect()

plt.savefig("massdistrib.pdf", bbox_inches='tight')
plt.show()
