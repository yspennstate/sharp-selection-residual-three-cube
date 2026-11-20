"""Exact helpers on the sign vectors of the cube: products, determinants and the vertices of a column triple."""
from itertools import product

SIGNS = tuple(product((-1, 1), repeat=3))

def dot(a, b):
    return sum(x*y for x, y in zip(a, b))

def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def det(a, b, c):
    return dot(a, cross(b,c))

def positive(v):
    for x in v:
        if x:
            return v if x>0 else tuple(-y for y in v)
    return v

def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))

def vertices(cols):
    return [tuple(dot(row,s) for row in zip(*cols)) for s in SIGNS]
