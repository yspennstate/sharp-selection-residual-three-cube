"""Enclosing-ball centers and pairwise separated tetrahedra of sign vectors."""
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
from check_grid import SIGNS, dot, sub, cross, det, positive, vertices

TETS=[t for t in combinations(range(8),4)
      if det(*(sub(SIGNS[i],SIGNS[t[0]]) for i in t[1:]))]
VOLUME=[abs(det(*(sub(SIGNS[i],SIGNS[t[0]]) for i in t[1:])))//8 for t in TETS]
AXES=[v for v in product(range(-2,3),repeat=3) if v!=(0,0,0) and positive(v)==v]
INTERVALS=[[(min(dot(s,SIGNS[i]) for i in t),max(dot(s,SIGNS[i]) for i in t))
            for s in AXES] for t in TETS]

def separated(i,j):
    return any(b<=c or d<=a for (a,b),(c,d) in zip(INTERVALS[i],INTERVALS[j]))

COMPAT=[sum(1<<j for j in range(i+1,len(TETS)) if separated(i,j)) for i in range(len(TETS))]

def fills(mask,need=6):
    if not need:
        return True
    while mask:
        bit=mask & -mask
        i=bit.bit_length()-1
        mask^=bit
        if VOLUME[i]<=need and fills(mask & COMPAT[i],need-VOLUME[i]):
            return True
    return False

def centers(points):
    """Candidate rational ball centers as integer numerator and positive-or-negative denominator."""
    p=points[0]
    disp=[sub(q,p) for q in points[1:]]
    if not disp:
        return p,1
    if len(disp)==1:
        return tuple(x+y for x,y in zip(p,points[1])),2
    if len(disp)==2:
        u,v=disp
        a,b,c=dot(u,u),dot(u,v),dot(v,v)
        den=2*(a*c-b*b)
        if not den:
            return None
        num=tuple(den*x+c*(a-b)*y+a*(c-b)*z for x,y,z in zip(p,u,v))
        return num,den
    cof=[cross(disp[(i+1)%3],disp[(i+2)%3]) for i in range(3)]
    den=2*dot(disp[0],cof[0])
    if not den:
        return None
    norms=[dot(q,q) for q in disp]
    return tuple(den*x+dot(col,norms) for x,col in zip(p,zip(*cof))),den

def admissible_tetrahedra(cols,limit):
    pts=vertices(cols)
    @cache
    def contained(ids):
        center=centers([pts[i] for i in ids])
        if center is None:
            return frozenset()
        num,den=center
        return frozenset(i for i,p in enumerate(pts) if
                         sum((den*x-y)**2 for x,y in zip(p,num))*limit.denominator
                         <= den*den*limit.numerator)
    good=0
    for i,t in enumerate(TETS):
        if any(all(j in contained(s) for j in t)
               for count in (2,3,4) for s in combinations(t,count)):
            good|=1<<i
    return good
