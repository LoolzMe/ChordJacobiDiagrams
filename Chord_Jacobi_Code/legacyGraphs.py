   
#OLD CODE IS GIVEN BELOW. SOME OF IT MAY STILL BE USEFUL



from sympy import flatten
from sympy import Matrix
from sympy import sympify
#import random

import time

from LinAlgebra import terms
from LinAlgebra import simplify
from LinAlgebra import mult
from LinAlgebra import toVect
from LinAlgebra import toLC


from math import factorial as fact

from itertools import chain

from itertools import product
"""
generates all elements in the Cartersian product
"""


from itertools import permutations as perm
from itertools import combinations as comb
"""
comb(L,r) generates combinations of length r of elements from L
(no repetitions).
"""
from itertools import combinations_with_replacement as combR
"""
for example, combR((1,2,3),2) generates (1,1), (1,2), (1,3), 
(2,2), (2,3), (3,3).
"""



from chords import graph2seq
from chords import Ch_diag
from chords import can_form_lc
from chords import gen_4T
from chords import circ_diag


"""
AUXILIARY FUNCTIONS
"""

#tested
def remove_rep(L):
    """
    L is a list or tuple (possibly with repetitions). remove_rep returns the 
    corresponding list of elements without repetitions. The order of 
    elements remains the same. 
    """
    out = []
    for el in L:
        if el not in out:
            out.append(el)
    return out

"""
Here is how one can test whether a list is sorted or not:...
 
L == sorted(L)

This command does NOT change the list L.
"""

#tested
def act_edge(t,ed):
    """
    t is a permutation (given as a tuple) and ed is an edge.
    the function returns the result of the action of this permutation on the edge ed.
    """
    if t[ed[0]-1]<t[ed[1]-1]:
        return (t[ed[0]-1], t[ed[1]-1])
    else: 
        return (t[ed[1]-1], t[ed[0]-1])


#tested
def act(t,G):
    """
    t is a permutation in S_n and G is a graph with n vertices.
    The function returns the result of the action of t on vertices of G
    The output IS SORTED.
    """
    out = [act_edge(t,e) for e in G]; out.sort()
    return out


#tested
def sgn(t):
    s = 1; n=len(t) # t is an element of S_n 
    for c in comb(range(1,n+1),2):
        i,j = c
        if t[i-1]>t[j-1]:
            s=s*(-1)
    return s


#
def to_vertices(G):
    """
    A uni-trivalent graph G is given in the edge format;
    the function returns the vertex presentation of G with some 
    choice of cyclic orders in the neighborhoods of vertices; 
    for example, to_vertices([(3, 4),(3, 4),(1,3),(2,4)]) returns 
    ([(1, 2, 3), (1, 2, 4)], (3, 4));
    we assume that the set of vertices is {1,2,...,n};
    note each neighborhood of a vertex is given by a "sorted" tuple.
    """
    n_v = max(flatten(G))
    vert = []
    for a in range(n_v):
        vert.append([])
    #entries of vert will be neighborhoods of vertices
    for k in range(len(G)):
        i,j=G[k]
        vert[i-1].append(k+1)
        vert[j-1].append(k+1)
    hr=[]; V=[]
    for v in vert:
        if len(v)==1:
            hr.append(v[0])
        else:
            V.append(tuple(v))
    hr.sort(); V.sort()
    return (V,tuple(hr))


#tested
def val(i,G):
    """
    returns the valency of the vertex i in the graph G (G is in the edge format)
    """
    vert= flatten(G)
    return vert.count(i)

#tested
def val_prof(G):
    """
    G is a graph in the edge format;
    the function returns the valency profile of G;
    for instance, every uni-trivalent graph has the valency 
    profile (3,3,...,3,1,...,1);
    we assume that the set of vertices is {1,2,...,n}
    """
    n = max(flatten(G)) #number of vertices
    vL=sorted([val(i,G) for i in range(1,n+1)])
    return tuple(vL[::-1])

#tested
def to_Vstandard(G):
    """
    G is a graph in the edge format;
    the function returns the isomorphic graph in the 
    valency standard format; for the valency standard format, 
    vertices with higher valency have smaller indices
    """
    vert=flatten(G); n = max(vert) #n is the number of vertices
    def fun(i):
        """
        fun(i) returns (-1)*valency of i in G 
        """
        return -vert.count(i)
    L = list(range(1, n+1))
    L.sort(key=fun)
    t=tuple(L.index(j)+1 for j in range(1,n+1)) #we need the inverse of the permutation
    return act(t, G)


#tested
def attach_hr(G,e):
    """
    G is a graph in the edge format and e in an edge of G;
    the function returns the new graph which is obtained from 
    G by attaching a "hair" to the edge e
    """
    if e not in G:
        print('Sorry, the graph does not have edge ',e)
        return None
    (i,j)=e
    n = max(flatten(G)) #the number of vertices of G
    out=G[:]; out.remove(e)
    out+=[(i,n+1),(j,n+1),(n+1,n+2)]
    return sorted(out)
    
#tested
def body(r):
    """
    returns the 'body' of the caterpillar with r rungs;
    the output is in the edge format;
    if r==1, the function returns the theta graph
    """
    if r==1:
        return [(1,2) for i in range(3)]
    out=[(1,2),(2*r-1,2*r)]
    out+=[(2*i-1,(2*i)) for i in range(1,r+1)]
    out+=[(2*i-1,2*i+1) for i in range(1,r)]
    out+=[(2*i,2*i+2) for i in range(1,r)]
    return sorted(out)


"""
jac_diag is the class for uni-trivalent graphs (i.e. Jacobi diagrams).  
CONVENTIONS: our graph can have multiple edges, but they cannot have loops.
A trivalent vertex can be incident to at most one hair. 

The first n1 vertices 1,2,...,n1 are trivalent they are incident to hairs. 
The last n2 trivalent vertices n1+1,... n1+n2 are trivalent and 
none of them is incident to a hair. Vertices with labels 
n1+n2+1,... 2*n1+n2 are univalent and hairs are 
(1,n1+n2+1), ..., (1, 2*n1+n2).   

G (in the input of the jac_diag) is the list of edges which 
connect two trivalent vertices.  

So G is the "shaved" version of Jacobi diagram: 
it has n1 + n2 vertices, the first n1 vertices 
are bivalent and the last n2 vertices are trivalent.

2 n1 + 3 n2 = 2e 
n1 + n2 = n


"""
#sgn_aut will be tested, even should be tested...
class Jac_diag:
    def __init__(self, G):
        self.shaved = G
        self.e = len(G) #the number of edges which are not incident to univalent vertices
        self.n = max(flatten(G)) #the number of trivalent vertices of the Jacobi diagram
        #self.vert3 = tuple(range(1,self.n+1))#this is the tuple of vertices of the shaved graph
        self.n1 = 3*self.n - 2*self.e# the number of trivalent vertices which are incident to hairs
        self.n2 = 2*(self.e - self.n)# the number of trivalent vertices which are not incident to hairs.
        
    def __eq__(self, other):
        G1=self.shaved[:]; G1.sort()
        G2=other.shaved[:]; G2.sort() #the order of edges should not play a role
        return G1==G2
    
    def full(self):
        """
        returns the graph with all its hairs
        """
        n1=self.n1; n2 =self.n2
        out = self.shaved + [(i,n1+n2+i) for i in range(1,n1+1)]
        return sorted(out)
    
    def is_aut_sh(self,t):
        """
        t is a permutation in S_{n1}\times S_{n2}
        the function returns True if t is an automorphism of the shaved graph.
        Otherwise, it returns False.
        """
        G=self.shaved
        for ed in G:
            ee = act_edge(t,ed)
            if G.count(ed)!=G.count(ee):
                return False
        return True
    
    def aut(self):
        """
        generates all automorphisms of the shaved graph.
        NOTE that if our graph has multiple edges then the automorphism 
        group is the product what is generated by this function and 
        S_{2}\times S_{2} \times ...
        """
        n1 = self.n1; n = self.n
        tup1 = tuple(range(1,n1+1))
        tup2 = tuple(range(n1+1,n+1))
        for s in perm(tup1):
            for t in perm(tup2):
                if self.is_aut_sh(s+t):
                    u = tuple(n + s[i-1] for i in range(1,n1+1))# the action on univalent vertices
                    yield s+t+u
    
    def sgn_aut(self,t):
        """
        returns the sign of the automorphism t of the corresponding Jacobi diagram.
        It is assumed that t is an automorphism of the graph with hairs.
        """
        G = self.full(); s=1
        for ed in G:
            if t[ed[0]-1] > t[ed[1]-1]:
                s = s*(-1)
        return s * sgn(t)
    
    def even(self):
        """
        returns True if the Jacobi diagram is even. 
        Otherwise, False.
        """
        for t in self.aut():
            if self.sgn_aut(t)==-1:
                return False
        return True
            
                    
    def act(self,t):
        """
        t is a permutation in S_{n1} \times S_{n2}
        act(self,t) returns the result of the action of this permutation 
        on the shaved graph.
        """
        return [act_edge(t,e) for e in self.shaved]
    
    def can_form(self):
        """
        returns the canonical form of the shaved graph (we use S_{n_1} \times S_{n_2})
        """
        out = self.shaved[:]; out.sort()
        tup1 = tuple(range(1,self.n1+1))
        tup2 = tuple(range(self.n1+1,self.n+1))
        for s in perm(tup1):
            for t in perm(tup2):
                new = self.act(s+t)
                new.sort()
                if new < out:
                    out = new[:]
        return out
    
    def isom(self,other):
        """
        returns True if the corresponding graphs are isomorphic. 
        Otherwise, it returns False. 
        """
        if (self.n, self.e)!=(other.n, other.e):
            return False
        return self.can_form()==other.can_form()
    
    def nb_vert(self,v):
        """
        Returns the list of vertices adjacent to v (no repetitions).
        We assume that our graph does not have loops. 
        VD: Do we need this function?
        """
        adj = []
        for e in self.full():
            if v in e:
                vv = [j for j in e if j!=v]
                if vv[0] not in adj:
                    adj.append(vv[0])
        return adj
            
        
        
"""
functions for producing examples: 
"""
#tested
def wh(n):
    """
    returns the cycle graph with n vertices (i.e. the shaved version of the wheel)
    """
    return [(i,i+1) for i in range(1,n)] + [(1,n)]

#tested
def whJD(n):
    """
    returns the wheel with n hairs (as the instance of Jac_diag)
    """
    return Jac_diag(wh(n))

#tested
def ladd(r):
    """
    returns the shaved ladder with r rungs, r should be > 1. 
    """
    if r == 2:
        return wh(4)
    else:
        out = [(1,2),(3,4),(1,5),(2,6),(3,2*r-1),(4,2*r)]
        out = out + [(2*j+5,2*j+6) for j in range (r-2)]# adding r-2 rungs
        out = out + [(2*j+5,2*j+7) for j in range (r-3)]# additional left vertical edges 
        out = out + [(2*j+6,2*j+8) for j in range (r-3)]# additional right vertical edges 
    return out

#tested
def laddJD(r):
    """
    returns the ladder with r rungs (r >1) as the instance of Jac_diag
    """
    return Jac_diag(ladd(r))


"""
Functions for generating trivalent and bi-trivalent graphs.
A graph is represented by the (sorted) list of its edges
Every edge is an ordered tuple (i,j) with i<j

For example, the tetrahedron (K_4) 
is [(1,2),(1,3),(1,4), (2,3), (2,4),(3,4)] 
We may have loops and multiple edges. 

Bivalent vertices (if any) have smaller labels than those 
for trivalent ones. 
"""

"""
The command max(flatten(G)) returns the number of vertices of G. 
We assume that G does not have vertices of valency = 0.
set(flatten(G)) gives us the set of vertices of G
"""

"""
Note that two bi-trivalent graphs are isomorphic iff 
they have the same canonical form.
"""            
   
#tested
def nb_vert(G,v):
    
    """
    Returns the set of vertices adjacent to v.
    """
    adj=set()
    for e in G:
        if v in e:
            if e[0]==v:
                adj=adj.union({e[1]})
            else:
                adj=adj.union({e[0]})
    return adj 


def nb_vert_l(G,v):
    """
    Returns the list of vertices adjacent to v
    YW: Loops are counted in this one
    """
    adj = [] #MW: I changed adj from a set to a list, since it was unable to output list(adj)
    vv = [e for e in G if v in e]
    for e in vv:
        adj.append(e[1])
        adj.append(e[0])
    return remove_rep(adj) 



#VD: Here is another function which tests connectivity
def conn(G):
    """
    Returns True if G is a connected graph
    and False otherwise. It is assumed that G does not 
    have vertices of valency 0.
    """
    orb={1}; add_prev={1}
    while True:
        add_next = set()
        for i in add_prev:
            for j in nb_vert(G,i):
                if j not in orb:
                    add_next=add_next.union({j})
        if add_next==set():
            break
        orb = orb.union(add_next); add_prev=add_next
        
    return len(orb)==max(flatten(G))



"""
There are exactly two isomorphism classes of trivalent 
graphs with 2 vertices: Theta and dmbbll.
Both graphs are connected and dmbbll has 2 loops.  
"""

#tested
def add_vert(G,j):
    """
    G is a bi-trivalent graph and 
    1 <= j <= len(G) is one of its edges.  
    The function returns the bi-trivalent graph by inserting 
    an extra vertex in the "middle of" edge j
   
    YW: will insert loop instead if j = 0
    """
    if j == 0: #the zero flag means you're going to add a loop.
        out = [(i+1,k+1) for (i,k) in G]
        out = [(1,1)]+out
    else: 
        out = [(i+1, k+1) for (i,k) in G[:j-1]+ G[j:]]
        (i,k) = G[j-1]; out = [(1,i+1),(1,k+1)]+out
    out.sort()
    return out


#tested
def add_edge(G):
    """
    G is a bi-trivalent graph with the valency profile (2,n-2)
    The output is obtained from G by adding extra edge (1,2)
    """
    if val_prof(G)[0]!=2:
        return('Your graph has a wrong valency profile')
    out = [(1,2)]+G; out.sort()
    return out


    
"""
Suppose that I have a list of all isomorphism classes of 
trivalent graphs with n-2 vertices.
"""  

#tested for trivalent graphs with 4 vertices and with 6 vertices.
def gener_tri(G):
    """
    Returns a list of 
    only isomorphism classes of graphs.
    """
    for j in range(1, len(G)+1):
        out = add_vert(G,j)
        for k in range(0, len(out)+1): #JW: there appears to have been an indentation error in previous verison
            out1 = add_vert(out,k)
            out1 = add_edge(out1)
            out1.sort()
            yield can_form(out1)


         
"""
Various examples of bi-trivalent graphs. 
I checked that all the lists are sorted. 
"""
tetra = [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]

Th = [(1,2),(1,2),(1,2)] #this is the Theta-graph(shaving is not necessary)

dmbbll=[(1,1),(1,2),(2,2)]
# two loops connected by an edge.


diamond = [(1,3),(1,4),(2,3),(2,4),(3,4)]
# BTW, the corresponding Jacobi diagram is even. 

house = [(1,4),(1,5),(2,3),(2,4),(3,5),(4,5)]
#the house does not coincide with its canonical form.

dmbbll1 = [(1,3),(1,3),(2,4),(2,4),(3,4)]

util = [e for e in product((1,2,3),(4,5,6))] # gives us the utility graph
#the utility graph is not in the canonical form.

tritri= [(1,2),(1,5),(2,5), (3,4),(3,6),(4,6),(5,6)]
# two triangles attached by one edge

L2 = [Th,dmbbll] # all isom classes of trivalent graphs with 2 vertices



"""
We need to deal with (linearly) closed Jacobi diagrams.
A closed Jacobi diagram is a uni-trivalent graph with a TOTAL 
order on the set of its univalent vertices. 
Equivalently, we can fix the "orientation" on the closed 
Jacobi diagram by chosing a cyclic order in the neighborhood of 
each trivalent vertex. 

VD: I suggest to represent closed Jacobi diagrams as 
(sorted) lists of edges of the corresponding uni-trivalent graph. 
The trivalent vertices are labeled by 1,..., n and 
univalent vertices are labeled by n+1, ..., n+n1

I suggest to introduce a class for closed Jacobi diagrams.


For example, the Mercedes Benz sign is represented by 
[(1,2),(1,3),(1,4)]

The closure of the wheel with 2 hairs has two terms:

wh2 = [(1,2),(1,2),(1,3),(2,4)] 

len(J) is the number of edges 

max(flatten(J)) is the total number of vertices. 

n + n1 = num_v
3*n + n_1 = 2 num_e 
"""


#tested
def num_triv(J):
    """
    returns the number of trivalent vertices of the closed Jacobi diagram J
    """
    num_v= max(flatten(J))
    return (2*len(J)-num_v)//2 

#tested
def num_univ(J):
    """
    returns the number of univalent vertices of the closed Jacobi diagram J
    """
    num_v= max(flatten(J))
    return num_v - num_triv(J)

#used for testing
def num_univ_t(J):
    """
    returns the number of univalent vertices of the closed Jacobi diagram J
    """
    num_v= max(flatten(J))
    return (3*num_v-2*len(J))//2 


                    


"""
The stuff below is for testing:
"""

#mb = [(1,2),(1,3),(1,4)]

#wh2 = [(1,2),(1,2),(1,3),(2,4)] 


#tested
def all_gra(n,e):
    """
    Generates all possible graphs with n vertices and e edges.
    Multiple edges and loops are allowed.
    Be careful! It is an "insane" generator...
    """
    Edges = list(combR(range(1,n+1),2)) 
    # Edges is the list of all possible edges. It is automatically sorted.
    for c in combR(Edges,e):
        yield list(c)


#tested
def trivalent(G):
    """
    Returns true if a given graph is trivalent, false otherwise.
    It does NOT work for the graph with the empty set of edges. 
    """
    vert= flatten(G) # the list of vertices with meaningful repetitions
    n = max(vert)
    if 2*len(G) != 3*n:
        return False
    for i in range(1,n+1):
        if vert.count(i)!=3:
            return False
    return True
            

def conn_list(L):
    """
    L is a list (or tuple) of graphs. The function returns the list 
    [...,True, False,...] whose ith entry is True in L[i] is connected 
    and False otherwise.
    """
    t0 = time.time()
    out = [conn(G) for G in L]
    print(time.time()-t0)
    return out


def conn_list_M(L):
    """
    L is a list (or tuple) of graphs. The function returns the list 
    [...,True, False,...] whose ith entry is True in L[i] is connected 
    and False otherwise.
    """
    t0 = time.time()
    out = [connected(G) for G in L]
    print(time.time()-t0)
    return out

"""
VD: I tested function conn and connected using the above functions 
conn_list and conn_list_M on relatively large lists of graphs. 
The results are the same, however, conn_list is (approx.) 2 times faster.  
"""


