import matplotlib.pyplot as plt
import numpy as np
import numpy.random as rnd
from scipy.stats.distributions import norm, gamma

def rejection_sample(proposal_sampler,proposal_pdf,target,k):
    while True:
		g = proposal_sampler()
		u = rnd.uniform()
		if u < target(g)/(k*proposal_pdf(g)):
			return g


u,v    = 5, 1
sigma  = 1.0
lm, ls = 3, 1

def get_max():
		L = []
		x0 = -(np.sqrt(u**2+2*u-4*sigma**2+1)-u+1)/2
		x1 =  (np.sqrt(u**2+2*u-4*sigma**2+1)+u-1)/2
		if u > x0:
			L.append(pl_no_gamma(x0))
		if u > x1:
			L.append(pl_no_gamma(x1))
		L.append(pl_no_gamma(u))
		L.append(pl_no_gamma(0))
		return np.max(L)

l = np.arange(100)*0.1

def pl_no_gamma(l):
	if u > l:
		return   1/(l+1.0)*norm.pdf(np.linalg.norm([u-l,v]),0,sigma)
	return     1/(l+1.0)*norm.pdf(v,0,sigma)

def pl(l):
	if u > l:
		return   1/(l+1.0)*norm.pdf(np.linalg.norm([u-l,v]),0,sigma)*gamma.pdf(l,lm,ls)
	return     1/(l+1.0)*norm.pdf(v,0,sigma)*gamma.pdf(l,lm,ls)

L = []
def rs():
	cnt = 0
	while True:
		g = rnd.gamma(lm, ls)
		u = rnd.uniform()
		if u < pl(g)/(1*gamma.pdf(g,lm,ls)):
			L.append(cnt)
			return g
		cnt += 1

#samples = [rs() for x in xrange(5000)]
#print "average number of rejects: ", np.mean(L)
#n, bins, patches = plt.hist(samples, 100, normed=1, facecolor='green', alpha=0.75)
#plt.close()
#y = [40*pl(b) for b in bins]
#y = [gamma.pdf(b,lm,ls) for b in bins]
m = get_max()
print "m=",m
plt.plot(l,gamma.pdf(l,lm,ls),'b--',linewidth=1)
plt.plot(l,[pl(x)/m for x in l],'r--',linewidth=1)
#plt.plot(l,[pl_no_gamma(x) for x in l],'r--',linewidth=1)
plt.plot(l,[m for x in l],'g--',linewidth=1)
#plt.plot(bins,y,'r--',linewidth=1)
plt.show()


