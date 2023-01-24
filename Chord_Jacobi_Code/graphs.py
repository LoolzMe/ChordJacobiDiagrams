#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The contributors to this package were (reversed alphabetic order): Mathhew Wynne, 
Yeahuay Wu, Vasily Dolgushev and Elif Altinay-Ozaslan.
"""

"""
Graphs CAN have mulitple edges but they cannot have loops. 


We use two formats to represent graphs: edge format and vertex format. 

In the EDGE FORMAT, a graph G is presented by the SORTED list of its edges. 
Each edge is represent as the tuple of vertices (i,j) and we assume that 
i < j. We think of (i,j) as the edge from i to j.  
For example, the unitrivalent graph with three vertices of val = 3
and three hairs is represented as
G = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6)]. 
Note that the list is sorted. 

When we work with unitrivalent graphs (in the edge format), it is 
good to assume that 1,2,...,n1 are trivalent and n1+1,n1+2, ..., n1+n2
are univalent.

In the VERTEX FORMAT, a graph G is presented by the sorted list 
of "neighborhoods" of vertices, i.e. list of tuples; for instance
(i1, i2,...,iv) represents the vertex incident to 
edges i1, i2,... ,iv.
For example, the above graph 
G = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6)] is represented 
in the vertex format by the list 
[(1,2,3), (1,4,5), (2,4,6), (3,), (5,),(6,)]. 
Another example is the Mercedes Benz graph [(1,2), (1,3), (1,4)]. 
Its vertex representation is [(1,), (1, 2, 3), (2,), (3,)]. 

When we deal with unitrivalent graphs, it makes sense to omit tuples 
corresponding to univalent vertices. For example, the Mercedes Benz graph
is represented by [(1,2,3)].
"""

#from queue import * # do we really need queue? 


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
ABOUT IMPORTED FUNCTION

mult(k,v) returns the result of multiplication a linear combination v by k;

For graph2seq(G), G is a graph (in edge format) which represents a chord diagram;
The function returns the sequence representation of this chord diagram. 

Instanced of the class Ch_diag are (circular) chord diagrams.

"""

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
END OF AUXILIARY FUNCTIONS
"""

#tested
def is_aut(t, G):
    """
    a tuple t represents a permutation of {1,2,...,n}
    and G is a graph whose set of vertices is {1,2,...,n};
    the command returns True if t is an automorphism of G;
    otherwise the command returns False;
    we assume that G has no double edges
    """
    for e in G:
        ee=act_edge(t,e)
        if ee not in G:
            return False
    return True
        


"""
A closed Jacobi diagram is represent by the tuple (V,hr) where V 
is the sorted list of (neighborhoods of) trivalent vertices and hr is the tuple 
of hairs (i.e. edges incident to univalent vertices). 
Now every trivalent vertex is a tuple
of length 3 (i,j,k) where i,j,k are the edges incident to this vertex.
The set of edges is the set of e distinct integers.

It is possible that some hairs are not incident to trivalent vertices. 
In fact, if the Jacobi diagram does not have trivalent vertices then 
it is a chord diagram. In this case, every hair is an edge which connects 
two univalent vertices.  

We say that a tuple (V,hr) is groomed if the set of edges is {1,2..., e}
where e is the number of edges. 
"""


#tested
def groom(J):
    """
    J = (V,hr) represents a closed Jacobi diagram and (V,hr) is possibly 
    ungroomed. The function returns the groomed representation of this 
    closed Jacobi diagram  
    """
    V,hr = J[0], J[1]
    L = remove_rep(flatten(V)+list(hr)); L.sort()
    def groom_it(t):
        return tuple(L.index(a)+1 for a in t)
    newV = [groom_it(v) for v in V]
    return (newV,groom_it(hr))


"""
Instances of Jac_cl represent (linearly) closed Jacobi diagrams. 
V is the sorted list of trivalent vertices represented like this 
(i,j,k) where i,j,k are edges incident to vertex (i,j,k). 
The cyclic order of the subset of edges incident to the vertex 
is determined by the total order of the tuple. 
CONVENTION: i is the smallest index in {i,j,k} in the 
sense of the usual order. 

We assume that the set of edges is {1,2,...,e}

hr is the tuple which represents the total order on the set of hairs 
(i.e. edges incident to univalent vertices)

var = (V,hr)

The current version of __eq__ is the naive one. 
It is possible that two different instances of the 
class Jac_cl represent isomorphic closed Jacobi diagrams. 

What should we do with the orientation? 
It tried to introduce Or as a default variable and it does not 
work in __init__


"""
class Jac_cl:
    def __init__(self, var):
        self.var = var # the instance variable 
        self.V= var[0] # sorted list of triv vertices
        self.hr = var[1] # the total order on the set of hairs
        self.vert = self.V + [(i,) for i in self.hr] # the list of all vertices
        #note that the output of self.vert is often NOT sorted 
        self.univ = len(var[1]) # number of univalent vertices
        self.triv = len(var[0]) # number of trivalent vertices
        self.num_v = len(var[1])+ len(var[0]) # total number of vertices
        self.e = (3*len(var[0]) + len(var[1]))//2 #number of edges
        self.inc_triv = set(flatten(var[0])) # the set of edges incident to trivalent vertices
    #tested
    def __eq__(self, other):
        return self.var==other.var
    #tested
    def good_hr(self):
        """
        returns the tuple of hairs which are indicent to trivalent 
        vertices
        """
        return tuple(i for i in self.hr if i in self.inc_triv)
    #tested
    def not_chdiag(self):
        """
        returns True if self is NOT a chord diagram, i.e. it has 
        at least one trivalent vertex.
        """
        return len(self.V) > 0
    #tested
    def edge(self,i):
        """
        i belongs to the set of edges of self;
        the method returns the corresponding edge as the 
        tuple of numbers (j_1,j_2); j_1 and j_2 are incident to edge i
        """
        n = self.num_v
        return tuple(j+1 for j in range(n) if i in self.vert[j])
        
    #tested
    def to_edges(self):
        """
        returns the representation of 'self' in the edge format;
        note that the information about the orientation is lost
        """
        e = self.e
        out = [self.edge(i) for i in range(1,e+1)]
        return sorted(out)
    #tested
    def triv_inc_hr(self,h):
        """ 
        returns the trivalent vertex incident to a given good hair;
        if h is not a good hair then the method returns nothing
        """
        for v in self.V:
            if h in v:
                return v
    
    #tested
    def stu(self):
        """
        'applies' STU relation to the left most good hair
        of self. The output is a linear combination of closed 
        Jacobi diagrams
        """
        i = self.good_hr()[0] # the left most good hair
        v = self.triv_inc_hr(i) # the unique trivalent vertex incident to hair i
        if i==v[0]:
            j,k=v[1],v[2]
        else:
            if i==v[2]:
                j,k=v[0], v[1]
            else:
                j,k = v[2], v[0]
        newV = self.V[:] # [:] is important. We do not want to change the instance self.
        newV.remove(v) # since self.V is sorted, newV is also sorted
        ind = self.hr.index(i)
        hr1 = self.hr[:ind] + (k,j) + self.hr[ind+1:]
        hr2 = self.hr[:ind] + (j,k) + self.hr[ind+1:]
        term1 = (1, groom((newV,hr1)))
        term2 = (-1,groom((newV,hr2)))
        return [term1, term2]


"""
Examples of instances of Jac_cl:...
"""

T = Jac_cl(([(1,2,3),(1,4,5),(2,5,6)],(3,4,6)))

MB = Jac_cl(([(1,2,3)],(1,2,3)))

TT = Jac_cl(([(1, 2, 3),(1, 11, 7), (2, 9, 5), (3, 4, 12),(4, 5, 6),(7, 8, 9),(10, 11, 12)],(6,8,13,13,10)))

# two version of the theta graph; note that the theta graph
# has no univalent vertices
Theta = Jac_cl(([(1,2,3),(1,3,2)],()))
Theta1 = Jac_cl(([(1,2,3),(1,2,3)],()))

#tested
def stu_lc(v):
    """
    v is a linear combination of linearly closed Jacobi diagrams;
    each graph in this linear combination has the same number of
    trivalent vertices and it is given in the vertex format; 
    the function returns the result of 
    applying stu to this linear combination
    """
    def stu_term(tt):
        """
        tt is a tuple (coefficent, graph);
        the function returns the coefficient*stu(graph)
        """
        return mult(tt[0],Jac_cl(tt[1]).stu())
    out =[]
    for t in v:
        out+=stu_term(t)
    return out

#
def to_lc_CDs(Vh):
    """
    Vh is a closed Jacobi diagram given in the vertex format; 
    the function returns the corresponding linear combination 
    of chord diagram (in the sequence format)
    """
    num3=len(Vh[0])#number of trivalent vertices
    if num3==0:
        return [(1,list(Vh[1]))]
    v = Jac_cl(Vh).stu() #stu is applied once
    for i in range(num3-1):
        v=stu_lc(v)
    return [(t[0],list(t[1][1])) for t in v]
    

#
def cl_Jac_to_CD(Vh):
    """
    Vh is a closed Jacobi diagram given in the vertex format;
    the function averages over all total orders of hairs and
    returns the corresponding linear combination 
    of chord diagram (in the canonical sequence format)
    """
    V, hr = Vh; num3 = len(V) #num3 is the number of 3-valent vertices
    lc =[]
    for p in perm(hr):
        lc+=Jac_cl((V,p)).stu()
    #we averaged over the total orders on hairs and applied the STU once
    for i in range(num3-1):
        lc=stu_lc(lc)
    return can_form_lc([(t[0],list(t[1][1])) for t in lc])

 


    """ 
    A class that represent an open Jacobi diagrams.
    Open Jacobi diagrams are quite similar to the closed ones,
    although it can have graph structures outside the main circle.
    
    There is orientation 
    """

# Transfered to surface_dynamics package for sage
# class OpenJDiag:
#     def __init__(self, G):
#         self.G = G
#         self.n_hairs = countHairs(G)


#     # def ApplyAS(self, n):
#     #     pass

 

# def countHairs(G):
#     fl_G = flatten(G)
#     dict = { }
#     for i in fl_G:
#         if i in fl_G:
#             dict[i] += 1
#         else:
#             dict[i] = 1

#     n_hairs = 0
#     for i in dict.items:
#         if i > 1:
#             n_hairs += 1


#     return n_hairs