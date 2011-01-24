# vim:ts=4:sw=4:sts=4:et:ai
import numpy as np
import numpy.random as rnd
import pydot

alpha0 = 5
Lambda = 1
gamma  = 0.2

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
        return n

class datum(object):
    def __init__(self, id):
        self.id = id

class tree_process(object):
    def __init__(self,gen, alpha0, Lambda, gamma):
        self.gen    = gen
        self.alpha0 = alpha0
        self.Lambda = Lambda
        self.gamma  = gamma
    def alpha(self,n):
        return Lambda**n.depth * alpha0
    def sift(self, n, d):
        b = rnd.uniform() < 1.0/(self.alpha(n)+1)
        if b:
            n.data.append(d)
            n.lendata += 1
            return
        N = [x.lennodes + x.lendata for x in n.subnodes]
        N.append(gamma)
        N = np.array(N).astype("float")
        N /= sum(N)
        s = rnd.multinomial(1,N,1)[0]
        res = np.where(s)[0]
        if res == len(n.subnodes):
            n.subnodes.append(self.gen(n))
            self.sift(n.subnodes[-1],d)
        else:
            self.sift(n.subnodes[res],d)
        n.lennodes += 1

class node(object):
	def __init__(self, depth, name):
        self.data     = []
        self.subnodes = []
        self.name     = name
        self.lendata  =  0
        self.lennodes =  0
        self.depth    =  depth
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

class edge2Dnode(node):
    def __init__(self, depth, name, angle, length):
        node.__init__(self,depth, name)
        self.angle  = angle
        self.length = length
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


def gentree(name,kappa=0.5,minlen=1,samples=30,alpha0=10,Lambda=.5,gamma=.2):
    gen = edge2dgenerator(8,1)

    tp = tree_process(gen,alpha0,Lambda,gamma)

    N = edge2Dnode(0,[],0,1)
    for i in xrange(samples):
        tp.sift(N,datum(i))
    print N

    G = pydot.Dot('Tree', graph_type="digraph")
    root = pydot.Node("root")
    G.add_node(root)
    N.dot(G,root,np.zeros(2))
    G.write_png(name,prog='neato')

if __name__ == "__main__":
    gentree("G.png")

