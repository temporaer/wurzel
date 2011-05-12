import numpy as np
import matplotlib.pyplot as plt
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

fig = plt.figure(figsize=(6,8))
fig.subplots_adjust(hspace=0.1,wspace=0.1,left=0.10,bottom=0.07,right=1,top=0.85)
xticklab = ["horizontal", "vertical"]
if True:
	for i in xrange(Z.shape[0]):
		x = ang[i]
		y = Z[i]
		z = mass[i]
		img[y,x] += z

	ax  = fig.add_subplot(121)
	#cax = ax.matshow(img/imgn,cmap="binary",norm=LogNorm(vmin=min_mass,vmax=max_mass))
	cax = ax.matshow(img/imgn,cmap="binary")
	#ax.set_xticks([x for x in np.arange(0,resx+.01,resx/3)])
	ax.set_yticks([y for y in np.arange(0,resy+.01,resy/10)])
	#xticklab = ["%1.2f"%x for x in np.arange(min_ang,max_ang+.01,(max_ang-min_ang)/3)]
	ax.set_yticklabels(["%3.0f"%y for y in np.arange(min_Z,max_Z+.01,(max_Z-min_Z)/10)])
	ax.set_ylabel('depth (mm)')
	ax.set_title("Mass distribution")
	ax.set_ylim(resy-1.5,-0.5)
	plt.setp(ax.get_xticklabels(), visible=False)

	divider = make_axes_locatable(ax)
        axMarg  = divider.append_axes("bottom", 1.2, pad=0.1, sharex=ax)
	plt.setp(axMarg.get_yticklabels(), visible=False)
	axMarg.plot(img.sum(axis=0))
	axMarg.set_xticks([3.5, resx-5.5])
	axMarg.set_xticklabels(xticklab)
	axMarg.set_aspect("auto","box-forced")
	axMarg.set_xticklabels(xticklab)
	axMarg.set_xlabel('angle to vertical (rad)')
	axMarg.set_xlim(-0.5,resx-1.5)
	#axMarg.apply_aspect()
	#axMarg.set_ylim(


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
	ax.set_ylim(resy-1.5,-0.5)
	#ax.set_yticks([y for y in np.arange(0,resy+.01,resy/10)])
	ax.set_xticklabels(["%1.2f"%x for x in np.arange(min_ang,max_ang+.01,(max_ang-min_ang)/3)])
	#ax.set_yticklabels(["%3.0f"%y for y in np.arange(min_Z,max_Z+.01,(max_Z-min_Z)/10)])
	plt.setp(ax.get_yticklabels(), visible=False)
	#ax.set_ylabel('depth (mm)')
	ax.set_title("Expected Mass")
	plt.setp(ax.get_xticklabels(), visible=False)

	divider   = make_axes_locatable(ax)
        axMarg = divider.append_axes("bottom", 1.2, pad=0.1, sharex=ax)
	plt.setp(axMarg.get_yticklabels(), visible=False)
	axMarg.plot((img/imgn).sum(axis=0))
	axMarg.set_xticklabels(xticklab)
	axMarg.set_xticks([3.5, resx-5.5])
	axMarg.set_xlabel('angle to vertical (rad)')
	axMarg.set_xlim(-0.5,resx-1.5)

plt.savefig("massdistrib.pdf", bbox_inches='tight')
plt.show()
