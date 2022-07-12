from typing import Callable

class Solver:
    def __init__(self, seq: str, 
        m: int,
        a: int,
        b: int,
        c: int,
        eH: Callable[[int, int], int],
        eS: Callable[[int, int], int],
        eL: Callable[[int, int], int]
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
        if (not self.match(i,j)) or j-i <= self.m:
            return float('inf')
        # Case 1: Hairpin loop
        r1 = self.eH(i, j)
        # Case 2: Stacking loop
        r2 = self.V(i+1, j-1) + self.eS(i,j)
        # Case 3: Internal loop
        # TODO: bottle neck, maybe use heuristic to limit interior loop size
        r3 = float('inf')
        for i2 in range(i+1, j):
            for j2 in range(i2+1, j):
                if i2 == i+1 and j2 == j-1:
                    # stacking loop has already been considered in Case 2
                    continue
                r3 = min(self.V(i2, j2)+ self.eL(i,j,i2,j2), r3)      
        # Case 4: Multiloop
        r4 = self.a + self.WM2(i+1, j-1)  
        return min(r1, r2, r3, r4) 
    
    def WM2(self, i: int, j: int):
        # base
        if i == j or i-1 == j:
            return float('inf')
        # Case 1: j is not paired
        r1 = self.WM2(i+1, j) + self.c
        # Case 2: j is paired with k
        # Note: j cannot pair with i
        r2 = float('inf')
        for k in range(i+1, j):
            energy = self.V(i,k) + self.b + self.WM(k+1,j)
            r2 = min(energy, r2)
        return min(r1, r2)
    
    def WM(self, i: int, j: int):
        # base
        if i == j or i-1 == j:
            return float('inf')
        # Case 1: j is not paired
        r1 = self.WM(i+1, j) + self.c
        # TODO: we can combine case 2 and case 3
        # Case 2: j is paired with k and there are more loops outside of [k,j]
        r2 = float('inf')
        # Note: j can pair with i
        for k in range(i+1, j+1):
            energy = self.V(i, k) + self.b + self.WM(k+1,j)
            r2 = min(energy, r2)
        # Case 3: j is paired with k and there are no more loops
        r3 = float('inf')
        for k in range(i+1, j+1):
            energy = self.V(i, k) + self.b + (j-k)*self.c
            r3 = min(energy, r3)
        return min(r1, r2, r3)

    def solve(self):
        N = len(self.seq)
        return self.V(0, N-1)

m = 1
a, b, c = -10, -1, -1
eH = lambda i, j: i-j
eS = lambda i, j: i-j
eL = lambda i, j, i2, j2: i-j + i2-j2

rna = "AUCUCUAGGCU"
solver = Solver(rna, m, a, b, c, eH, eS, eL)


def debug(solver: Solver):
    print()
    print(' ' * 5, end='')
    for i in range(len(rna)):
        print(f'{rna[i]:>5}', end='')
    print()
    print(' ' * 5, end='')
    for i in range(len(rna)):
        print(f'{i:>5}', end='')
    print()
    for i in range(len(rna)):
        print(f'{i:>5}', end="")
        for j in range(len(rna)):
            if not j - i >= -1:
                print(' '*5, end="")
                continue
            result = solver.W(i, j)
            print(f'{result:>5}', end="")
        print()

# debug(solver)



