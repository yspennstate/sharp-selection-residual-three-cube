# The Sharp Selection Residual for Set-Valued Maps on the Three-Cube

Yitzchak Shmalo

For upper semicontinuous maps $F_1,\dots,F_n$ from $[0,1]^3$ to its nonempty compact subsets there are a point $x$ and
points $y_i\in F_i(x)$ with $\|x-\frac1n\sum_i y_i\|\le 3\sqrt3/(4n)$, and for every $n\ge3$ the constant is attained by
maps whose values have at most two points. The Shapley-Folkman-Starr theorem gives the constant $3/2$ in place of $3\sqrt3/4$.

`main.tex` and `main.pdf` are the paper.

## The full-row certificate

Section 5 of the paper proves the case of a full row by an exact rational certificate. The program uses only Python's
standard library:

    cd code
    python check_full_row.py

`check_full_row.py` imports `check_grid.py` and `check_cover.py`; each of these two modules also has a main program of
its own, for an earlier grid check that the paper does not use. It prints three lines, recorded in
`code/full_row_checks.txt`:

    bound 27/4 depth 24 alternate False counts {'cover': 1776, 'independent': 20, 'unproved': 0} matrices 414
    bound 27/4 depth 28 alternate True counts {'cover': 1632, 'independent': 22, 'unproved': 0} matrices 411
    bound 20/3 depth 28 alternate False counts {'cover': 2834, 'unproved': 0} matrices 668

The first line is the run the paper uses. The second is an independent run, with the starting simplices reflected and
the other choice among tied longest edges; it also leaves no simplex uncertified. The third is the run with $27/4$
replaced by $20/3$ that Remark 5.2 uses.

## License

MIT for the code; the paper is CC BY 4.0.
