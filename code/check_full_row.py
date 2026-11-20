"""Subdivide the four-parameter families with a full spatial row, using rational certificates."""
from fractions import Fraction as F
from functools import cache
from itertools import permutations
from math import gcd
from collections import Counter
from check_grid import dot
from check_cover import admissible_tetrahedra,fills

def verify(depth,alternate=False,bound=F(27,4)):
    unit=1 << (depth+1)
    count=Counter()
    @cache
    def covered(mask):
        return fills(mask)
    @cache
    def evaluate(second,x):
        a,b,c,d=x
        rows=[(unit,unit,unit),(unit,a,b),(unit,c,d) if second==0 else (c,unit,d)]
        g=gcd(unit,*x)
        cols=[tuple(v//g for v in col) for col in zip(*rows)]
        limit=bound*(unit//g)**2
        return (sum(dot(v,v) for v in cols)<=limit,admissible_tetrahedra(cols,limit))
    def check(second,verts,left):
        data=[evaluate(second,x) for x in verts]
        if all(z[0] for z in data):
            count['independent']+=1
            return True
        common=data[0][1]
        for _,mask in data[1:]:
            common &=mask
        if covered(common):
            count['cover']+=1
            return True
        if left==0:
            count['unproved']+=1
            return False
        edges=[(sum((a-b)**2 for a,b in zip(x,y)),i,j)
               for i,x in enumerate(verts) for j,y in enumerate(verts[:i])]
        longest=max(e[0] for e in edges)
        choices=[(i,j) for length,i,j in edges if length==longest]
        i,j=choices[0] if alternate else choices[-1]
        mid=tuple((a+b)//2 for a,b in zip(verts[i],verts[j]))
        assert all(2*m==a+b for m,a,b in zip(mid,verts[i],verts[j]))
        children=[]
        for k in (i,j):
            child=list(verts)
            child[k]=mid
            children.append(check(second,child,left-1))
        return all(children)
    for second in (0,1):
        for perm in permutations(range(4)):
            point=[unit,-unit,unit,-unit] if alternate else [-unit]*4
            verts=[tuple(point)]
            for i in perm:
                point[i]=-point[i]
                verts.append(tuple(point))
            check(second,verts,depth)
    count.setdefault('unproved',0)
    print('bound',bound,'depth',depth,'alternate',alternate,'counts',dict(count),'matrices',evaluate.cache_info().currsize)
    return count['unproved']==0

if __name__=='__main__':
    assert verify(24)
    assert verify(28,alternate=True)
    assert verify(28,bound=F(20,3))
