# vim:ts=4:sw=4:sts=4:et:ai
import numpy as np
import numpy.random as rnd
import pydot
import matplotlib.pyplot as plt
from scipy.stats.distributions import norm

def rotmat(a):
    ca = np.cos(a)
    sa = np.sin(a)
    return np.array([[ca,-sa],
                     [sa,ca]])

class edge2dgenerator(object):
    def __init__(self, kappa, minlen):
        self.kappa = kappa
        self.minlen = minlen
    def __call__(self, parent):
        newname = [i for i in parent.name]
        newname.append(parent.lennodes)

        angle  = rnd.vonmises(parent.angle, self.kappa)
        length = self.minlen+rnd.uniform()

        n = edge2Dnode(parent.depth+1, newname,angle,length)
        n.pos = n.position(parent.pos)
        return n

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

class tree_sampler(object):
    def __init__(self,hp,gen):
        self.hp = hp
        self.gen = gen
    def generate(self,data, root):
        for d in data:
            d.node = root
            root.add_data(d)
        self.__sample_assignments(data,root)
    def __sample_assignments(self,data,root):
        for d in data:
            pslice = rnd.uniform(0,d.node.get_likelihood(d))
            umin, umax = 0, 1
            while True:
                u = rnd.uniform(umin,umax)
                n = self.__find_node(u,root)
                p = n.get_likelihood(d)
                if p > pslice:
                    d.node = n
                    n.add_data(d)
                    break
                if n < d.node: umin = u
                else:          umax = u
        root.update_lengths()
    def __sample_sticks(self,root):
        pass
    def __sample_nodes(self,root):
        pass
    def __find_node(self,u,n):
        if u < n.nu: return n
        u = (u-n.nu)/(1-n.nu)
        while u<1-np.prod(lambda x:1-x, n.psi):
            # draw a new psi-stick
            n.psi.append(rnd.beta(1,self.hp.gamma))
            n.subnodes.append(self.gen(n))
        e = [0]
        for x in n.psi:
            e = (1-(1-e[-1])*(1-x))
        for idx,ei in enumerate(e):
            if ei > u:
                u = (u-ei)/(e[idx+1]-ei)
                n = n.subnodes[idx]
                break
        self.__find_node(u,n)



class tree_process(object):
    def __init__(self,gen, hp):
        self.gen  = gen
        self.hp   = hp
    def alpha(self,n):
        return self.hp.Lambda**n.depth * self.hp.alpha0
    def sift(self, n, d):
        b = rnd.uniform() < 1.0/(self.alpha(n)+1)
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
	def __init__(self, depth,name):
        self.data     = []
        self.subnodes = []
        self.depth    =  depth
        self.name     =  name
        self.nu       = 0
        self.psi      = []
        self.lendata  =  0
        self.lennodes =  0
    def __lt__(self,other):
        return self.name < other.name
    def update_nu_psi(self):
        pass
    def add_data(self,d):
        self.data.append(d)
        self.lendata += 1
    def update_lengths(self):
        self.lendata = len(self.data)
        for s in self.subnodes:
            s.update_lengths()
        self.lennodes = sum([x.lendata+x.lennodes for x in self.subnodes])

class drawablenode(node):
	def __init__(self, depth, name):
        node.__init__(self,depth,name)
    def namestr(self):
        return "N" + ".".join(map(str,self.name))
    def __str__(self):
        desc = "[%d, %d] "%(self.lendata, self.lennodes)
        s = "  " * self.depth + desc + self.namestr() + "; " +",".join([str(d.id) for d in self.data]) + "\n"
        for n in self.subnodes:
            s += str(n)
        return s
    def dot(self, g, root):
        for d in self.data:
            n = pydot.Node(str(d.id))
            n.set_style('filled')
            n.set_fillcolor('red')
            g.add_node(n)
            g.add_edge(pydot.Edge(root,n))
        for n in self.subnodes:
            n2 = pydot.Node(n.namestr(), label="XXX")
            g.add_node(n2)
            g.add_edge(pydot.Edge(root,n2))
            n.dot(g,n2)

class edge2Dnode(drawablenode):
    def __init__(self, depth, name, angle, length):
        drawablenode.__init__(self,depth, name)
        self.angle  = angle
        self.length = length
        self.width  = 0.01
        self.pos    = np.zeros(2)

    def get_likelihood(self, d):
        pos = d.pos.copy()
        pos -= self.pos
        pos = np.dot(rotmat(-self.angle),pos)
        if pos[0]<0 or pos[0]>self.length: return 0
        return norm.pdf(pos[1],0,self.width)/self.length

    def position(self,ppos):
        pos = ppos + self.length * np.array((np.cos(self.angle), np.sin(self.angle)))
        return pos
    def dot(self, g, root, ppos):
        for n in self.subnodes:
            pos = self.position(ppos)
            pos2 = pos
            color = "red" if len(n.data)>0 else "black"
            n2 = pydot.Node(n.namestr(), label="%d"%len(n.data), pos="%2.2f, %2.2f!"%(pos2[0],pos2[1]),color=color)
            g.add_node(n2)
            g.add_edge(pydot.Edge(root,n2))
            n.dot(g,n2,pos)
    def sample(self, d):
        distfromline = rnd.normal(0,self.width)
        lenonline    = rnd.uniform(0,self.length)
        d.pos = self.pos + np.dot(rotmat(np.pi + self.angle),np.array((lenonline,distfromline)))

def gentree(name,kappa=3,minlen=1,samples=1000,alpha0=5,Lambda=.5,gamma=.25):
    gen = edge2dgenerator(kappa,minlen)

    hp = hyperparam(alpha0,Lambda,gamma)
    tp = tree_process(gen,hp)

    N = edge2Dnode(0,[],0,1)
    L = []
    for i in xrange(samples):
        L.append(datum(i))
        tp.sift(N,L[-1])
    print N

    #np.save("G.txt",np.array([a.pos for a in L]))
    x = [a.pos[0] for a in L]
    y = [a.pos[1] for a in L]
    plt.plot(x,y, ".")

    #G = pydot.Dot('Tree', graph_type="digraph")
    #root = pydot.Node("root")
    #G.add_node(root)
    #N.dot(G,root,np.zeros(2))
    #G.write_png(name,prog='neato')

    plt.show()

if __name__ == "__main__":
    gentree("G.png")

