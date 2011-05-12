# vim:ts=4:sw=4:sts=4:et:ai
import numpy as np
import numpy.random as rnd
import pydot
from progressbar import ProgressBar, Percentage, Bar, ETA, \
     RotatingMarker
from scipy.stats.distributions import norm, vonmises, gamma, expon
import matplotlib.pyplot as plt

def rejection_sample(proposal_sampler,proposal_pdf,target,k):
    while True:
		g = proposal_sampler()
		u = rnd.uniform()
		if u < target(g)/(k*proposal_pdf(g)):
			return g

def rotmat(a):
    ca = np.cos(a)
    sa = np.sin(a)
    return np.array([[ca,-sa],
                     [sa,ca]])
def sample_pa(x,d,sigma):
    first=(x[1]-d[1]).mean()
    second=(x[0]-d[0]).mean()
    print "x, d, sigma ", x, d, sigma
    cnt = 0
    L = []
    z = first / np.sqrt(first**2+second**2) 
    L.append(-np.arccos(-z))
    L.append(-np.arccos( z))
    L.append( np.arccos(-z))
    L.append( np.arccos( z))
    Lmax = max([norm.pdf(np.cos(a)*first + np.sin(a)*second,0,sigma) for a in L])
    while True:
        u=rnd.uniform(0,Lmax)
        a=rnd.uniform(-np.pi/2.,np.pi/2.)
        normal = norm.pdf(np.cos(a)*first + np.sin(a)*second,0,sigma) # right?
        if u < normal:
            break
        cnt += 1
    print cnt
    return a

def post_l(l,x,d,alpha,sigma):
    #normal = [norm.pdf(np.cos(alpha)*(x[1]-d[1]) + np.sin(alpha)*(x[0]-d[0])) for d in data]\
    mean=(np.cos(alpha)*(d[1]) + np.sin(alpha)*(d[0])).mean()
    std=np.sqrt(np.cos(alpha)*(x[1]) + np.sin(alpha)*(x[0]))
    normal = norm.pdf(l,mean,std) # sigma?
    #normal = norm.pdf(np.cos(alpha)*(x[1]-d[:,1]) + np.sin(alpha)*(x[0]-d[:,0])) 
    return np.prod(normal)

class edge2d_base_distrib(object):
    def __init__(self,length_a0,length_b0,angle_kappa):
        self.length_a0 = length_a0
        self.length_b0 = length_b0
        self.angle_kappa = angle_kappa
    def __call__(self, parent, psi=None, nu=None):
        """ sample from prior """
        newname = [i for i in parent.name]
        newname.append(len(parent.subnodes))

        parlen   = parent.length
        dist     = rnd.uniform(0,parlen)
        paralpha = parent.angle
        beta     = np.pi / 4.
        pos      = parent.pos + dist * np.array([np.cos(paralpha), np.sin(paralpha)])
        p2       = pos + parlen * np.array([np.cos(paralpha+beta), np.sin(paralpha+beta)])
        angle    = sample_pa(pos,p2,parent.width) # TODO ugly width
        k        = 4.
        length   = parlen + parlen * ((1.-rnd.uniform())/(k-1.))**(k-1.)
        n = edge2Dnode(parent.depth+1, newname,pos, length,angle,parent,nu,psi)

        return n
    def adjust(self, node):
        """ sample from posterior """
        D  = np.vstack([d.pos for d in node.data ])
        old=np.array([node.get_likelihood_angle(d) for d in node.data ])
        old=np.log(old).mean()
        #import ipdb as pdb; pdb.set_trace()
        D -= node.parent.pos
        D  = D.T

        old_angle = node.angle
        # determine posterior of angle given prior and data
        x_ang = np.arctan2(D[1,:],D[0,:])
        R1 = self.angle_kappa * np.cos(node.parent.angle) + np.sum(np.cos(x_ang))
        R2 = self.angle_kappa * np.sin(node.parent.angle) + np.sum(np.sin(x_ang))
        #R1 =  np.sum(np.cos(x_ang))
        #R2 =  np.sum(np.sin(x_ang))
        mu = np.arctan2(R2,R1)
        Rn = R1 / np.cos(mu)
        node.angle = vonmises.rvs(Rn,loc=mu)
        #node.angle = mu

        X = D.copy()

        # rotate around newly drawn angle
        D  = np.dot(rotmat(-node.angle), D)
        assert D.shape[1] == len(node.data)


        # determine posterior of length given prior and data
        aN = self.length_a0 + D.shape[1]/2.0
        bN = self.length_b0 + D.shape[1]/2.0 * np.var(D[0,:])
        #aN =  D.shape[1]/2.0
        #bN =  D.shape[1]/2.0 * np.var(D[0,:])
        #import ipdb; ipdb.set_trace()
        node.length = gamma.rvs(aN,loc=1.0/bN)
        #node.length = np.var(D[0,:])

        old_pos = node.pos.copy()

        node.update_position()

        def f():
            plt.plot(X[0,:],X[1,:],".b")
            #plt.plot(D[0,:],D[1,:],".r")
            plt.plot(node.pos[0]-node.parent.pos[0],node.pos[1]-node.parent.pos[1],"*r")
            plt.plot(old_pos[0]-node.parent.pos[0],old_pos[1]-node.parent.pos[1],"*b")
            #plt.plot(node.parent.pos[0],node.parent.pos[1],"*r")
            plt.show()

        new=np.array([node.get_likelihood_angle(d) for d in node.data ])
        new=np.log(new).mean()
        print('old likelihood: %f     new likelihood: %f'%(old,new))
        #f()

class datum(object):
    def __init__(self, id):
        self.id   = id
        self.pos  = None
        self.node = None

class hyperparam(object):
    def __init__(self,alpha0,Lambda,gamma):
        self.alpha0 = alpha0
        self.Lambda = Lambda
        self.gamma  = gamma

class adams_tree(object):
    def __init__(self,hp):
        self.hp = hp
    def alpha(self,depth):
        return self.hp.Lambda**depth * self.hp.alpha0

class tree_sampler(adams_tree):
    def __init__(self,hp,base_distrib):
        adams_tree.__init__(self,hp)
        self.base_distrib = base_distrib
    def get_psi_nu(self, depth):
        psi = rnd.beta(1,self.hp.gamma)
        nu  = rnd.beta(1,self.alpha(depth))
        return psi, nu
    def generate(self,data, root):
        for d in data:
            d.node = root
        for i in xrange(1):
            self.__sample_assignments(data,root)
            self.__sample_sticks(root)
            self.__sample_nodes(root)
    def clear_data(self, n):
        del n.data[:]
        for x in n.subnodes:
            self.clear_data(x)
    def __sample_assignments(self,data,root):
        self.clear_data(root)
        widgets = ['Allocating: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ', ETA()]
        #pbar = ProgressBar(widgets=widgets, maxval=len(data)).start()
        for idx, d in enumerate(data):
            #pbar.update(idx)
            pslice = rnd.uniform(0,d.node.get_likelihood(d))
            umin, umax = 0.0, 1.0
            while True:
                u = rnd.uniform(umin,umax)
                n = self.__find_node(u,root)
                p = n.get_likelihood(d)
                if p >= pslice:
                    d.node = n
                    n.add_data(d)
                    break
                if n < d.node: umin = u
                else:          umax = u
        #pbar.finish()
        root.update_lengths()
    def __sample_sticks(self,node):
        node.nu = rnd.beta(node.lendata+1,node.lennodes+self.alpha(node.depth))
        for idx, n in enumerate(node.subnodes):
            n.psi = rnd.beta(node.lennodes+1, self.hp.gamma + np.sum([x.lennodes for x in node.subnodes[idx+1:]]))
        for n in node.subnodes:
            self.__sample_sticks(n)
    def __sample_nodes(self,node):
        if node.lendata > 0:
            self.base_distrib.adjust(node)

        for n in node.subnodes:
            self.__sample_nodes(n)
    def __find_node(self,u,n):
        if u < n.nu: return n
        u = (u-n.nu)/(1-n.nu)
        while u > 1.0-np.prod([1.0-s.psi for s in n.subnodes]) or len(n.subnodes)==0:
            # draw a new psi-stick
            psi,nu = self.get_psi_nu(n.depth+1)
            n.subnodes.append(self.base_distrib(n,psi=psi,nu=nu))
        e = [0]
        for x in n.subnodes:
            e.append( (1-(1-e[-1])*(1-x.psi)) )
        for idx,ei in enumerate(e):
            if ei > u:
                u = (u-e[idx-1])/(ei-e[idx-1])
                m = n.subnodes[idx-1]
                break
        return self.__find_node(u,m)



class tree_process(adams_tree):
    def __init__(self,gen, hp):
        adams_tree.__init__(self,hp)
        self.gen  = gen
    def sift(self, n, d):
        b = rnd.uniform() < 1.0/(self.alpha(n.depth)+1)
        if b:
            n.add_data(d)
            n.sample(d)
            return
        N = [x.lennodes + x.lendata for x in n.subnodes]
        N.append(self.hp.gamma)
        N = np.array(N).astype("float")
        N /= sum(N)
        #print "  "*n.depth, "prob: ",  N, " x=",1.0/(self.alpha(n)+1)
        s = rnd.multinomial(1,N,1)[0]
        res = np.where(s)[0]
        if res == len(n.subnodes):
            n.subnodes.append(self.gen(n))
            self.sift(n.subnodes[-1],d)
        else:
            self.sift(n.subnodes[res],d)
        n.lennodes += 1

class node(object):
	def __init__(self, depth,name,nu=None,psi=None):
        self.data     = []
        self.subnodes = []
        self.depth    =  depth
        self.name     =  name
        self.nu       =  nu
        self.psi      =  psi
        self.lendata  =  0
        self.lennodes =  0
    def __lt__(self,other):
        return self.name < other.name
    def add_data(self,d):
        self.data.append(d)
        self.lendata += 1
    def update_lengths(self):
        self.lendata = len(self.data)
        for s in self.subnodes:
            s.update_lengths()
        self.lennodes = sum([x.lendata+x.lennodes for x in self.subnodes])

class drawablenode(node):
	def __init__(self, depth, name,nu=None,psi=None):
        node.__init__(self,depth,name,nu,psi)
    def namestr(self):
        return "N" + ".".join(map(str,self.name))
    def __str__(self):
        desc = "[%d, %d] "%(self.lendata, self.lennodes)
        s = "  " * self.depth + desc + self.namestr() + "; " +",".join([str(d.id) for d in self.data]) + "\n"
        for n in self.subnodes:
            s += str(n)
        return s
    #def dot(self, g, root):
    #    for d in self.data:
    #        n = pydot.Node(str(d.id))
    #        n.set_style('filled')
    #        n.set_fillcolor('red')
    #        g.add_node(n)
    #        g.add_edge(pydot.Edge(root,n))
    #    for n in self.subnodes:
    #        n2 = pydot.Node(n.namestr(), label="XXX")
    #        g.add_node(n2)
    #        g.add_edge(pydot.Edge(root,n2))
    #        n.dot(g,n2)

class edge2Dnode(drawablenode):
    def __init__(self, depth, name, pos, length,angle,parent=None,nu=None,psi=None):
        drawablenode.__init__(self,depth, name,nu,psi)
        self.width  = 0.1
        self.length = length
        self.angle  = angle
        self.pos    = pos
        self.parent = parent
        self.vonmisesscale = 100

    def get_likelihood(self, d):
        """" sample a data point for this edge """
        pos = d.pos - self.parent.pos
        pos = np.dot(rotmat(-self.angle), pos)
        lik = halfnorm.pdf(pos[0],scale=self.length) * \
              vonmises.pdf(np.arctan2(pos[1],pos[0]),self.vonmisesscale,loc=self.angle)
        #assert lik!=0.0
        return lik

    def get_likelihood_angle(self, d):
        """" sample a data point for this edge """
        pos = d.pos - self.parent.pos
        lik = vonmises.pdf(np.arctan2(pos[1],pos[0]),self.vonmisesscale,loc=self.angle)
        return lik
    def update_position(self):
        self.pos = self.parent.pos + 3.0/4.0 * np.dot(rotmat(self.angle), [self.length,0])

    def dot(self, g, parent=None):
        pos = self.pos
        color = "red" if len(self.data)>0 else "black"
        #me  = pydot.Node(self.namestr(), label="%d"%len(self.data), pos="%2.2f, %2.2f!"%(pos[0],pos[1]),color=color)
        me  = pydot.Node(self.namestr(), label=self.namestr(), pos="%2.2f, %2.2f!"%(pos[0],pos[1]),color=color)
        g.add_node(me)
        if parent:
            g.add_edge(pydot.Edge(parent,me))
        for d in self.data:
            dn = pydot.Node(str(d.id), label="X", pos="%2.2f, %2.2f!"%(d.pos[0],d.pos[1]),color="green")
            g.add_node(dn)
            g.add_edge(pydot.Edge(me,dn,color="gray"))

        for n in self.subnodes:
            if n.lennodes>0 or True:
                n.dot(g,me)

    def sample(self, d):
        """" sample a data point for this edge """
        dist     = rnd.uniform(0,self.length)
        w        = rnd.normal(0,self.width)
        d.pos    = np.dot(rotmat(self.angle), [dist, w]) + self.pos
        d.ownpos = self.pos

def plot_tree(name,N):
    G = pydot.Dot('Tree', graph_type="digraph")
    N.dot(G)
    G.write_png(name,prog='neato')

def gentree(name,kappa=8,samples=1000,alpha0=5,Lambda=.5,gamma=.25):
    base_distrib = edge2d_base_distrib(3,1,kappa)

    hp = hyperparam(alpha0,Lambda,gamma)
    tp = tree_process(base_distrib,hp)

    #def __init__(self, depth, name, pos, length,angle,parent=None,nu=None,psi=None):
    N = edge2Dnode(0,[],np.zeros(2),10,0)

    L = []
    for i in xrange(samples):
        L.append(datum(i))
        tp.sift(N,L[-1])

    import matplotlib.pyplot as plt
    np.save("G.txt",np.array([a.pos for a in L]))
    x = [a.pos[0] for a in L]
    y = [a.pos[1] for a in L]
    plt.plot(x,y, ".b")
    x = [a.ownpos[0] for a in L]
    y = [a.ownpos[1] for a in L]
    plt.plot(x,y, "*c")
    plt.axis('equal')

    #plot_tree("G.png", N)


    plt.show()

    return N, L, hp, base_distrib

def test_inference():
    N,data,hp,base_distrib,virtual_parent = gentree("G.png")
    ts = tree_sampler(hp,base_distrib)
    psi,nu = ts.get_psi_nu(0)
    N2 = base_distrib(virtual_parent,nu=nu)
    virtual_parent.subnodes = []
    virtual_parent.subnodes.append(N2)
    ts.generate(data,N2)
    plot_tree("H.png", virtual_parent)


if __name__ == "__main__":
    gentree("G.png")
    #test_inference()

