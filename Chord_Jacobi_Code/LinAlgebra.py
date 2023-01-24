"""
Written by Vasily Dolgushev, Yeahuay (Joie) Wu and Matthew Wynne. 
Some functions are taken from a package written by Vasily Dolgushev 
and Geoffrey Schneider. 
"""


import itertools
import time
import math

from sympy import Matrix, zeros
from sympy import sympify
from math import factorial as fact

from itertools import combinations_with_replacement as cr

from itertools import combinations as comb
# comb(seq, r) is the generator for 'seq choose r'


"""
binomC(n,k) returns `n choose k'
It is just for testing.
"""
#tested
def binomC(n,k):
    return fact(n)//(fact(k)*fact(n-k))


"""
FUNCTIONS FOR WORKING WITH LINEAR COMBINATIONS:

A linear combination is typically represented like as a list 
[(c1,O1), (c2,O2), ...], where c1,c2,... are integers or SymPy rationals.
O1, O2 are objects of the same type. It is possible that 
Oi = Oj for i != j. Such linear combination is NOT simplified. 
We think of objects O1, O2,... as elements of a basis.  
"""

#tested
def terms(vec):
    """
    terms returns the list of "terms" of a linear combination vec.
    We need "terms" for "simplify".
    Note that there are no duplicates in terms(vec).
    The output of terms is SORTED. So we have to be careful.
    If elements of the basis are instances of a class then
    the computer should know how to compare them.
    """
    out = []
    for el in vec:
        if el[0]!=0 and el[1] not in out:
            out.append(el[1])
    return sorted(out)
    

#tested
def simplify(vec):
    """"
    simplifies a linear combination vec
    """
    out = [] 
    for ent in terms(vec):
        coeff = tuple(el[0] for el in vec if el[1] == ent)
        c = sum(coeff) # sum of all coefficients in front of ent
        if c != 0:
            out.append((c,ent))
    return out


#tested
def mult(k,x):
    """
    multiplies the linear combination x by the number k.
    Note that if x is not simplified then the output may be 
    non-simplified. 
    """
    if k ==0:
        return []
    else:
        return [(k*el[0], el[1]) for el in x]

#tested
def lsum(L):
    """
    L is a list or tuple of linear combinations.
    lsum returns the (simplified) sum of the linear combinations in the list L
    """
    newl = L[0]
    for k in range(1,len(L)):
        newl = newl + L[k]
    return simplify(newl)


#tested
def linExt(bfunc,x):
    """
    bfunc a function taking basis vectors to arbitrary linear combinations 
    of basis vectors. linExt allows us to extend it to the case when the argument 
    is a linear combination. Note that the output is simplified!
    """
    if x ==[]:
        return []
    l = [mult(k[0],bfunc(k[1])) for k in x]
    return lsum(l)


#tested
def bilinExt(bfunc, x, y):
    """
    bfunc a function taking basis vectors to arbitrary linear combinations 
    of basis vectors. "bilinExt" gives us the bilinear externsion of bfunc.
    Here x and y are linear combinations.
    The outout of bilinExt(bfunc,x,y) is of course a linear combination.
    Note: The simplification IS USED!
    Note that B(0,?) = B(?,0) = 0 for any bilinear function B.
    """
    if x == [] or y == []:
        return []
    l = []
    for k in x:
        for j in y:
            l.append(mult(k[0]*j[0],bfunc(k[1],j[1])))
    return lsum(l)


#tested    
def toVect(x, B):
    """
    toVect converts a simplified (!) linear combination with elements
    in the basis B into the corresponding coordinate vector.
    Note: x must be simplified and all terms of x must be multiples of 
    vectors from B.
    """
    def coeff(e):
        for t in x:
            if t[1] == e:
                return t[0]
        return 0
    return [coeff(e) for e in B]
    
        


#tested
def toLC(v,B):
    """
    toLC converts a coordinate vector v (with respect to the basis B) into 
    the corresponding linear combination.
    Note: v and B must have the same length!
    """
    if len(v)!=len(B):
        print('The coordinate vector has a wrong length!')
        return None
    return [(v[i],B[i]) for i in range(len(v)) if v[i]!=0]


#tested
def are_mult(v,w):
    """
    v,w are simplified linear combinations;
    the function returns True if v and w are linear 
    dependent; note that every coefficient in v and w is 
    non-zero
    """
    if v==[] and w==[]:
        return True
    n =len(v)
    if len(w)!=n:
        return False
    for i in range(n):
        if v[i][1]!=w[i][1]:
            return False
    c_v=v[0][0]; c_w=w[0][0] #the first coefficients
    vv = tuple(c_w*a[0] for a in v)
    ww = tuple(c_v*a[0] for a in w)
    return vv==ww
    
    

    """ Code that was successfully stolen from surface_dynamics
        Self-stealing isn't a thing, r8?
    """

class LinearSpace(object):
    def __init__(self, bases=[]) -> None:
        self._space = bases
        self._bases = bases
        self._matrix = zeros(1, len(bases))

    def bases(self):
        return self._bases

    def span(self, objects) -> None:
        self._space = objects
        self._bases = objects

class LinearObject(object):
    def __init__(self, space) -> None:
        self._bases = space.bases()
        self._scalars = [0] * len(self._bases)




"""
The following linear combinations are coordinate vectors were used for testing.
 

v = [(-3,(5,)), (2,(1,)), (0,(2,)), (-2,(1,))]

vv = [(7, (3,)), (-5,(1,)), (3,(8,)), (0,(3,))]

x = [0, 0, -4, 0, 2, 1,0, 0,0,0]

B = [(i,) for i in range(1,11)] # an example of a basis

"""
