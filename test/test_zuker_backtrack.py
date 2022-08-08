"""Test file for Zuker DP algorithm and Zuker max distance and distance vector algorithms"""

from functools import reduce
from itertools import product
from algs.zuker_backtrack import *
from algs.zuker_distance import DistanceSolver
import test.zuker_recursive as zuker_recursive
import random

def test_default_energy_functions():
    """Test the default parameters and energy functions"""
    for N in range(20,50):
        for _ in range(1):
            rna = generate_rna(N)
            # run Zuker algorithm to get optimal solutions
            solver = Solver(seq=rna)
            solver.fill_table()
            min_energy = solver.solve()
            print("Min energy is: ", min_energy)
            assert min_energy <= 0

            # a random fold to compare against
            afold = gen_random_fold(rna)

            # run distance solver to get max distance and vector
            distance_solver = DistanceSolver(rna, afold, solver.W, solver.V,
                                    solver.WM, solver.WM2)
            distance_solver.fill_distance_table()
            distance_solver.fill_vector_table()
            max_distance_dp = distance_solver.solve()[0]   
            vector_dp = distance_solver.solve()[1]
            print("Max distance is: ", max_distance_dp)
            print("The distance vector is: ", vector_dp)
            num_folds_dp = reduce((lambda x, y: x + y), vector_dp)
            print("The number of optimal folds is: ", num_folds_dp)
        print(N)
    print()
    print("All tests passed!")


def test_max_distance():
    """
    Test Zuker max distance algorithm
    by comparing with brute force algorithm
    """
    print (". = 1000 tests")
    count = 0
    for vals in product(range(-3, 2, 1), repeat=7):
        m, a, b, c = vals[:4]
        if m <= 0: continue
        eH = lambda s, i, j: vals[4]
        eS = lambda s, i, j: vals[5]
        eL = lambda s, i, j, i2, j2: vals[6]
        # try different values for bounding internal loop size. None means no bound
        L = list(range(0, 10))
        L.append(None)
        for bound in L:
            for N in range(1,10):
                for _ in range(2):
                    rna = generate_rna(N)
                    solver = Solver(rna, m, a, b, c, eH, eS, eL, bound)
                    solver.fill_table()
                    optimal_folds = backtrack(solver)

                    # a random fold to compare against
                    afold = gen_random_fold(rna)

                    max_distance = len(distancevector(afold, optimal_folds)) - 1

                    # result from solver
                    distance_solver = DistanceSolver(rna, afold, solver.W, solver.V,
                                            solver.WM, solver.WM2)
                    distance_solver.fill_distance_table()
                    max_distance_dp = distance_solver.solve()[0]  
                    distance_solver.fill_vector_table()
                    if max_distance != max_distance_dp:
                        max_distance_solution = distance_solver.get_one_max_dist_solution()
                        print("Given folding is: ", afold)
                        print("Ouputed farthest solution: ", max_distance_solution)
                        print(f"Max distance does not match for m={m} a={a} b={b} c={c} eH={vals[4]} eS={vals[5]} eL={vals[6]}, rna={rna} afold={afold} | expected={max_distance} got={max_distance_dp}")
                        debug2(rna, distance_solver.maxdist_W)
                        debug2(rna, distance_solver.maxdist_V)
                        debug2(rna, distance_solver.maxdist_WM2)
                        debug2(rna, distance_solver.maxdist_WM)
                    assert max_distance == max_distance_dp, \
                    f"Max distance does not match for m={m} a={a} b={b} c={c} eH={vals[4]} eS={vals[5]} eL={vals[6]}, rna={rna} afold={afold} | expected={max_distance} got={max_distance_dp}"
                    count += 1
                    if count % 1000 == 0:
                        print('.', end='', flush=True)
    print()
    print("All tests passed!")

def test_vector():
    """
    Test Zuker distance vector algorithm
    by comparing with brute force algorithm.
    """
    print (". = 100 tests")
    count = 0
    for vals in product(range(-3, 2, 1), repeat=7):
        m, a, b, c = vals[:4]
        if m <= 0: continue
        eH = lambda s, i, j: vals[4]
        eS = lambda s, i, j: vals[5]
        eL = lambda s, i, j, i2, j2: vals[6]
        # try different values for bounding internal loop size. None means no bound
        L = list(range(0, 10))
        L.append(None)
        for bound in L:
            for N in range(10, 30):
                for _ in range(2):
                    rna = generate_rna(N)
                    solver = Solver(rna, m, a, b, c, eH, eS, eL, bound)
                    solver.fill_table()
                    optimal_folds = backtrack(solver)

                    # a random fold to compare against
                    afold = gen_random_fold(rna)

                    vector = distancevector(afold, optimal_folds)

                    distance_solver = DistanceSolver(rna, afold, solver.W, solver.V,
                                                solver.WM, solver.WM2)
                    distance_solver.fill_vector_table()
                    vector_dp = distance_solver.solve()[1]

                    assert vector == vector_dp, \
                    f"Vector does not match for m={m} a={a} b={b} c={c} eH={vals[4]} eS={vals[5]} eL={vals[6]}, rna={rna} afold={afold} | expected={vector} got={vector_dp}"
                    num_folds_b = len(optimal_folds)
                    num_folds_dp = reduce((lambda x, y: x + y), vector_dp)
                    assert num_folds_b == num_folds_dp

                    count += 1
                    if count % 100 == 0:
                        print('.', end='', flush=True)

def gen_random_fold_h(rna: str, i: int, j: int, fold: list[int]):
    if i == j: return

    potential_indices = []
    for k in range(i+1, j):
        if match_of(rna[i]) == rna[k]:
            potential_indices.append(k)
    
    potential_indices.append(None)
    match_result = random.choice(potential_indices)

    if match_result is None:
        return

    k = match_result
    
    fold[i] = k
    fold[k] = i

    gen_random_fold_h(rna, i+1, k, fold)
    gen_random_fold_h(rna, k+1, j, fold)

def gen_random_fold(rna: str):
    """
    Given an rna, generate a random fold.
    """
    fold = [i for i in range(len(rna))]
    gen_random_fold_h(rna, 0, len(rna), fold)  
    return fold

matches = {'A': 'U', 'U': 'A', 'C':'G', 'G':'C'}
def match_of(c):
    """Return the complementary base"""
    return matches[c]

def backtrack(solver: Solver) -> list[int]:
    """
    Given a solver where .fill_table() has been called
    generate all foldings from that solver.

    Note: if you want to know the details of specific foldings,
    look into the optimal_folds2 list. It's the list of foldings with
    more details.
    """
    i, j = 0, len(solver.seq)-1
    initial_fold = [None for _ in range(len(solver.seq))]
    initial_fold2 = set()
    optimal_folds, optimal_folds2 = [], []
    for a, b in backtrackh(solver, 'W', i, j, initial_fold, initial_fold2):
        optimal_folds.append(a)
        optimal_folds2.append(b)

    
    # make sure that they are unique
    optimal_folds_u = set( tuple(fold) for fold in optimal_folds )
    assert len(optimal_folds) == len(optimal_folds_u), "Generated folds not unique"

    # convert None to self match
    for fold in optimal_folds:
        for i in range(len(fold)):
            if fold[i] is None:
                fold[i] = i

    return optimal_folds

def backtrackh(solver: Solver, tablename: str, i: int, j: int, fold: list, fold2: set):
    # size = j-i+1
    # if size <= 1:
    #     yield fold.copy()
    #     return
    
    T = get_table(solver, tablename)

    if len(T[i][j].eList) == 0:
        yield fold.copy(), fold2.copy()

    for eventname, match, eventlist in T[i][j].eList:
        if match is not None:
            folds, folde = match
            fold[folds] = folde
            fold[folde] = folds
            fold2.add((tablename, eventname, folds, folde))

        if len(eventlist) == 0:
            yield fold.copy(), fold2.copy()
        
        if len(eventlist) == 1:
            ntablename, ni, nj = eventlist[0]
            for f in backtrackh(solver, ntablename, ni, nj, fold, fold2):
                yield f
        
        if len(eventlist) == 2:
            ntablename, ni, nj = eventlist[0]
            for f1 in backtrackh(solver, ntablename, ni, nj, fold, fold2):
                ntablename2, ni2, nj2 = eventlist[1]
                for f2 in backtrackh(solver, ntablename2, ni2, nj2, fold, fold2):
                    yield f2

        if match is not None:
            fold[folds] = None
            fold[folde] = None
            fold2.remove((tablename, eventname, folds, folde))


def get_table(solver, tablename):
    """Get the corresponding DP table in the solver of the given name"""
    if tablename == 'W':
        return solver.W
    if tablename == 'V':
        return solver.V
    if tablename == 'WM2':
        return solver.WM2
    if tablename == 'WM':
        return solver.WM

def test_backtrack():
    """Test backtrack generates unique solutions"""
    print (". = 1000 tests")
    count = 0
    for vals in product(range(-3, 2, 1), repeat=7):
        m, a, b, c = vals[:4]
        if m <= 0: continue
        eH = lambda s, i, j: vals[4]
        eS = lambda s, i, j: vals[5]
        eL = lambda s, i, j, i2, j2: vals[6]
        # try different values for bounding internal loop size. None means no bound
        L = list(range(0, 10))
        L.append(None)
        for bound in L:
            for _ in range(2):
                rna = generate_rna(9)
                solver = Solver(seq=rna, m=m, a=a, b=b, c=c, eH=eH, eS=eS, eL=eL, internal_loop_size=bound)
                solver.fill_table()
                backtrack(solver)
                count += 1
                if count % 1000 == 0:
                    print('.', end='', flush=True)
    print()
    print("All tests passed!")

def test_params():
    """
    Test Zuker DP algorithm for different parameters
    by comparing DP and recursive algorithm
    """
    print (". = 1000 tests")
    count = 0
    for vals in product(range(-3, 2, 1), repeat=7):
        m, a, b, c = vals[:4]
        if m <= 0: continue
        eH = lambda s, i, j: vals[4]
        eS = lambda s, i, j: vals[5]
        eL = lambda s, i, j, i2, j2: vals[6]
        for _ in range(3):
            rna = generate_rna(5)
            dpsolver = Solver(seq=rna, m=m, a=a, b=b, c=c, eH=eH, eS=eS, eL=eL)
            dpsolver.fill_table()
            bfsolver = zuker_recursive.RecursiveSolver(rna, m, a, b, c, eH, eS, eL)
            for i in range(len(rna)):
                for j in range(len(rna)):
                    if not j - i >= -1: continue
                    dpans = dpsolver.W[i][j].val
                    bfans = bfsolver.W(i, j)
                    if dpans != bfans:
                        debug(rna, dpsolver, bfsolver)
                        raise AssertionError("values not equal")   
            count += 1
            if count % 1000 == 0:
                print('.', end='', flush=True)
    print()
    print(f"All tests {count} passed")

def generate_rna(length: int) -> str:
    """Generates an RNA of given length"""
    s = []
    for _ in range(length):
        c = random.choice(['A', 'U', 'C', 'G'])
        s.append(c)
    return ''.join(s)

def debug2(rna: str, table, table2=None):
    """Compare one table against another"""
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
            result1 = table[i][j].val
            result2 = table2[i][j] if table2 else ""
            s = str(result1) + '|' + str(result2)
            print(f'{s:>10}', end="")
        print()

def debug(rna: str, solver1: Solver, solver2: zuker_recursive.RecursiveSolver):
    """
    Compare the Zuker DP algorithm against Zuker recursive algorithm by
    comparing their W entries.
    """
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

def distancevector(fold: list, allfolds: list):
    """Compute the distance vector throught a brute force algorithm"""
    vector = [0] * len(fold)*2

    for other in allfolds:
        dist = foldcompare(fold, other)
        vector[dist] += 1

    # trim the end so that it looks clean
    while vector and vector[-1] == 0:
        vector.pop()
    return vector

def foldcompare(f1, f2):
    """
    Compute the distance between two folds.
    ([0, 1], [0, 1]) -> 0
    ([1, 0], [1, 0]) -> 0
    ([1, 0, 2], [2, 1, 0]) -> 2
    ([1, 0, 2], [0, 2, 1]) -> 2
    """
    count = 0
    # fold in f1 but not in f2
    for i in range(len(f1)):
        if f1[i] > i and f1[i] != f2[i]:
            count += 1
    # fold in f2 but not in f1
    for i in range(len(f2)):
        if f2[i] > i and f1[i] != f2[i]:
            count += 1
    return count

if __name__ == '__main__':
    # test_backtrack()
    # test_params()
    # test_max_distance()
    # test_vector()
    # test_default_energy_functions()
    pass
