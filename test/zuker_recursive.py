"""This is a file only used for testing purposes.
Implements Zuker recursive algorithm"""

from typing import Callable

class RecursiveSolver:
    def __init__(self, seq: str, 
        m: int,
        a: int,
        b: int,
        c: int,
        eH: Callable[[str, int, int], int],
        eS: Callable[[str, int, int], int],
        eL: Callable[[str, int, int, int, int], int]
    ):
        self.seq = seq
        self.m = m # minimum seperation for a pair
        self.a = a
        self.b = b
        self.c = c
        self.eH = eH
        self.eS = eS
        self.eL = eL
    

    def match(self, i: int, j: int):
        return (self.seq[i], self.seq[j]) in [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')]

    def W(self, i: int, j: int):
        # base
        if i == j or i-1 == j:
            return 0
        # case 1: i is not paired
        r1 = self.W(i+1, j)
        # case 2: i is paired with k
        r2 = float('inf')
        for k in range(i, j+1):
            r2 = min(self.V(i, k) + self.W(k+1, j), r2)
        return min(r1, r2)
    
    def V(self, i: int, j: int):
        # base
        if j-i < self.m or not self.match(i,j):
            return float('inf')
        # Case 1: Hairpin loop
        r1 = self.eH(self.seq, i, j)
        # Case 2: Stacking loop
        r2 = self.V(i+1, j-1) + self.eS(self.seq, i,j)
        # Case 3: Internal loop
        # TODO: bottle neck, maybe use heuristic to limit interior loop size
        r3 = float('inf')
        for i2 in range(i+1, j):
            for j2 in range(i2+1, j):
                if i2 == i+1 and j2 == j-1:
                    # stacking loop has already been considered in Case 2
                    continue
                r3 = min(self.V(i2, j2)+ self.eL(self.seq, i,j,i2,j2), r3)  
        # Case 4: Multiloop
        r4 = self.a + self.WM2(i+1, j-1)  
        return min(r1, r2, r3, r4) 
    
    def WM2(self, i: int, j: int):
        # base
        if i == j or i-1 == j:
            return float('inf')
        # Case 1: i is not paired
        r1 = self.WM2(i+1, j) + self.c
        # Case 2: i is paired with k
        # Note: i cannot pair with j
        r2 = float('inf')
        for k in range(i+1, j):
            energy = self.V(i,k) + self.b + self.WM(k+1,j)
            r2 = min(energy, r2)
        return min(r1, r2)
    
    def WM(self, i: int, j: int):
        # base
        if i == j or i-1 == j:
            return float('inf')
        # Case 1: i is not paired
        r1 = self.WM(i+1, j) + self.c
        # Case 2: i is paired with k and there are more loops in [k+1, j]
        r2 = float('inf')
        for k in range(i+1, j+1):
            energy = self.V(i, k) + self.b + self.WM(k+1,j)
            r2 = min(energy, r2)
        # Case 3: j is paired with k and there are no more loops in [k+1, j]
        # Note: i can pair with j
        r3 = float('inf')
        for k in range(i+1, j+1):
            energy = self.V(i, k) + self.b + (j-k)*self.c
            r3 = min(energy, r3)
        return min(r1, r2, r3)

    def solve(self):
        N = len(self.seq)
        return self.W(0, N-1)




