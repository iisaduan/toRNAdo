from itertools import product

from zuker_backtrack import *
import random

import zuker

# try all parameters
def test_params():
    count = 0
    for vals in product(range(-3, 1, 1), repeat=7):
        m, a, b, c = vals[:4]
        eH = lambda i, j: vals[4]
        eS = lambda i, j: vals[5]
        eL = lambda i, j, i2, j2: vals[6]
        for _ in range(3):
            rna = generate_rna(10)
            dpsolver = Solver(rna, m, a, b, c, eH, eS, eL)
            dpsolver.fill_table()
            bfsolver = zuker.Solver2(rna, m, a, b, c, eH, eS, eL)
            for i in range(len(rna)):
                for j in range(len(rna)):
                    if not j - i >= -1: continue
                    dpans = dpsolver.W[i][j].val
                    bfans = bfsolver.W(i, j)
                    if dpans == bfans:
                        count += 1
                        if count % 1000 == 0:
                            print('.', end='', flush=True)
                    else:
                        debug(dpsolver, bfsolver)
                        raise AssertionError("values not equal")             
    print()
    print(f"All tests {count} passed")

m = 0
a, b, c = -1, -1, -1
eH = lambda i, j: i-j
eS = lambda i, j: i-j
eL = lambda i, j, i2, j2: 0

rna = "ACGCGU"

new_solver = Solver(rna, m, a, b, c, eH, eS, eL)
new_solver.fill_table()
print(new_solver.solve())

def generate_rna(length: int) -> str:
    s = []
    for _ in range(length):
        c = random.choice(['A', 'U', 'C', 'G'])
        s.append(c)
    return ''.join(s)



def debug(solver1: Solver, solver2: zuker.Solver):
    print()
    print(' ' * 10, end='')
    for i in range(len(rna)):
        print(f'{rna[i]:>10}', end='')
    print()
    print(' ' * 10, end='')
    for i in range(len(rna)):
        print(f'{i:>10}', end='')
    print()
    for i in range(len(rna)):
        print(f'{i:>10}', end="")
        for j in range(len(rna)):
            if not j - i >= -1:
                print(' '*10, end="")
                continue
            result1 = solver1.W[i][j].val
            result2 = solver2.W(i,j) if solver2 else ""
            s = str(result1) + '|' + str(result2)
            print(f'{s:>10}', end="")
        print()

if __name__ == '__main__':
    # solver2 = zuker.Solver2(rna, m, a, b, c, eH, eS, eL)
    # debug(new_solver, solver2)
    test_params()