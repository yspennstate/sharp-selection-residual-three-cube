"""Exact three-dimensional grid check. Integer coordinates are m times sign-cube units."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from collections import Counter
import argparse

SIGNS = tuple(product((-1, 1), repeat=3))
PAIRS = tuple(combinations(range(3), 2))

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

def orbit(cols):
    # Row sign changes and row permutations act on the three columns.
    for p in permutations(range(3)):
        for s in SIGNS[:4]:  # simultaneous reversal is removed by positive()
            yield tuple(sorted(positive(tuple(v[p[i]]*s[i] for i in range(3))) for v in cols))

def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))

def vertices(cols):
    return [tuple(dot(row,s) for row in zip(*cols)) for s in SIGNS]

def certificate(cols, threshold, counts):
    pts = vertices(cols)
    norms = [dot(p,p) for p in pts]
    # If all generators are orthogonal, the center is farthest, by independent signs.
    if all(dot(cols[i], cols[j]) == 0 for i,j in PAIRS):
        return F(norms[0]), (0,0,0), 1
    best = threshold
    witness = None
    dual = [cross(cols[(i+1)%3], cols[(i+2)%3]) for i in range(3)]
    volume = abs(dot(cols[0],dual[0]))
    for ids in combinations(range(8),4):
        if tuple(sorted(7-i for i in ids)) < ids:
            continue
        counts['tetrahedra tested'] += 1
        a,*others = ids
        rows = [sub(pts[i],pts[a]) for i in others]
        rhs = [norms[i]-norms[a] for i in others]
        cof = [cross(rows[(i+1)%3],rows[(i+2)%3]) for i in range(3)]
        den = 2*dot(rows[0],cof[0])
        if not den:
            continue
        num = tuple(dot(z,rhs) for z in zip(*cof))
        if any(abs(dot(w,num)) > volume*abs(den) for w in dual):
            continue
        # Offset-free comparison of squared distances to the sphere center.
        powers = [den*den*q-2*den*dot(num,p) for q,p in zip(norms,pts)]
        if min(powers) < powers[a]:
            continue
        counts['empty spheres in box'] += 1
        radius = F(powers[a]+dot(num,num),den*den)
        if radius > best:
            best = radius
            witness = num,den
    if witness is None:
        return threshold, None, None
    return best,*witness

def solve(m):
    counts = Counter()
    dirs = sorted(v for v in product(range(-m,m+1),repeat=3) if v!=(0,0,0) and positive(v)==v)
    # The Hamming layer and the diagonal attain 17/3 using at most two varying groups.
    best = F(17*m*m,3)
    witness = 'diagonal pair and Hamming layer'
    if m%2==0:
        best = F(27*m*m,4)
        witness = 'perpendicular triple with diagonal -1/2 and off-diagonal 1'
    seen = set()
    lengths = {v:dot(v,v) for v in dirs}
    for cols in combinations(dirs,3):
        counts['direction triples'] += 1
        if sum(lengths[v] for v in cols) <= best or cols in seen:
            continue
        if not det(*cols):
            continue
        seen.update(orbit(cols))
        counts['full rank orbits above running variance cutoff'] += 1
        radius,num,den = certificate(cols,best,counts)
        if radius>best:
            best = radius
            witness = (cols,tuple(F(z,den) for z in num))
    print(f'm={m}, directions={len(dirs)}, unscaled squared radius={best/(m*m)}, three-group unit-cube squared residual={best/(36*m*m)}')
    print(f'  witness: {witness}')
    print('  '+', '.join(f'{label}={n}' for label,n in sorted(counts.items())))
    return best/(m*m)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('m',nargs='*',type=int,default=[1,2])
    args=parser.parse_args()
    for m in args.m:
        if m<1:
            parser.error('positive grid denominator required')
        solve(m)
