import numpy as np

def getu(i):
	base  = "data/GersteLA_192x192x410_normal-upsampled.dat"
	shape = (872,192,192)
	thresh = 0.04
	if i==0: return base,shape,thresh

	base  = "data/GersteLA_64x64x410-upsampled.dat"
	shape = (410,90,90)
	thresh = 0.04  * 64 / 192
	if i==1: return base,shape,thresh

	base  = "data/GersteLA_96x96x410.dat"
	shape = (436,96,96)
	thresh = 0.04  * 96 / 192
	if i==2: return base,shape,thresh

	base  = "data/GersteLA_128x128x410-upsampled.dat"
	shape = (582,128,128)
	thresh = 0.04  * 128 / 192
	if i==3: return base,shape,thresh

	base  = "data/GersteLA_256x256x410-upsampled.dat"
	shape = (1163,256,256)
	thresh = 0.04  * 256 / 192
	if i==4: return base,shape,thresh
def get(i):
	base  = "data/GersteLA_192x192x410_normal.dat"
	shape = (410,192,192)
	thresh = 0.04
	if i==0: return base,shape,thresh

	base  = "data/GersteLA_64x64x410.dat"
	shape = (410,64,64)
	thresh = 0.04  * 64 / 192
	if i==1: return base,shape,thresh

	base  = "data/GersteLA_96x96x410.dat"
	shape = (410,96,96)
	thresh = 0.04  * 96 / 192
	if i==2: return base,shape,thresh

	base  = "data/GersteLA_128x128x410.dat"
	shape = (410,128,128)
	thresh = 0.04  * 128 / 192
	if i==3: return base,shape,thresh

	base  = "data/GersteLA_256x256x410.dat"
	shape = (410,256,256)
	thresh = 0.04  * 256 / 192
	if i==4: return base,shape,thresh




for i in xrange(5):
	base,shape,thresh = getu(i)

	x=np.fromfile(base,dtype="float32").reshape(shape)

	s = x[5,:,:]
	b = s>(s.max()/1.2)
	t = np.median(s[b])
	t = 450000

	print "------- ", base, "--------"
	print "Stem: ", t

	y = x[15:,:,:] / t
	b = y>thresh
	print y[b].sum() *0.27 / np.prod(np.array(shape))*np.prod(np.array([410,192,192]))


