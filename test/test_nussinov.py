from algs.nussinov import *

import random

def test_max_distance_dp():
    """
    Generate rnas of different sizes,
    and test the DP algorithm against the
    brute force algorithm.
        - Test the max difference
        - Test the vector
    """
    test_range = range(5, 25)
    item_per_range = 1000
    for length in test_range:
        print(f"Testing length = {length}...")
        for _ in range(item_per_range):
            rna = generate_rna(length)

            _, opt_solutions = nussinov_dp(rna)
            folding = construct_one_opt_solution(rna, opt_solutions)

            maxdiff, _ = max_distance_dp(folding, opt_solutions)

            optimal_folds = backtrack(rna, opt_solutions)
            v = distancevector(folding, optimal_folds)
            maxdiff_bf = len(v)-1
            assert maxdiff == maxdiff_bf

            # Test histogram/vector calculations
            vdp = max_distance_vec(folding, opt_solutions).v
            assert vdp == v

def test_suboptimal_fold():
    """
    Generate rnas of different sizes,
    and test the DP algorithm against the
    brute force algorithm.
        - Test the max difference
        - Test the vector
    """
    test_range = range(5, 6)
    item_per_range = 10
    for length in test_range:
        print(f"Testing length = {length}...")
        for _ in range(item_per_range):
            rna = generate_rna(length)

            dist, opt_solutions = nussinov_dp(rna)
            folding = gen_random_fold(rna)

            maxdiff, _ = max_distance_dp(folding, opt_solutions)

            optimal_folds = backtrack(rna, opt_solutions)
            v = distancevector(folding, optimal_folds)
            maxdiff_bf = len(v)-1
            assert maxdiff == maxdiff_bf

            # Test histogram/vector calculations
            vdp = max_distance_vec(folding, opt_solutions).v
            assert vdp == v

            print(f'rna: {rna} | fold: {folding}')
            print(f'optdist: {dist}')
            for o in optimal_folds:
                print('-', o)
            print(f'vec: {vdp} ({maxdiff})')

matches = {'A': 'U', 'U': 'A', 'C':'G', 'G':'C'}
def match_of(c):
    return matches[c]

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

def test_nussinov():
    for _ in range(10):
        rna = generate_rna(6)
        opt, solgraph = nussinov_dp(rna)
        allopts = backtrack(rna, solgraph)
        print(f'rna: {rna} | {opt}')
        for o in allopts:
            print(o)

def generate_rna(length: int) -> str:
    s = []
    for _ in range(length):
        c = random.choice(['A', 'U', 'C', 'G'])
        s.append(c)
    return ''.join(s)

def backtrack(rna: str, opt_solutions: list) -> list:
    """
    Given the opt_solutions table from nussinov algorithm,
    return a list of all optimal RNA foldings.

    Example:
    Original String: CGCG
    opt_solutions Dict:
    {
        (0, 1): [('M', (0, 1))], # <unreachable>
        (1, 2): [('M', (1, 2))],
        (2, 3): [('M', (2, 3))],
        (0, 2): [('L', 0), ('M', (0, 1))], # <unreachable>
        (1, 3): [('L', 1), ('M', (1, 2))], # <unreachable>
        (0, 3): [('M', (0, 1)), ('M', (0, 3))]
    }

    Result:
    [
        [1, 0, 3, 2],
        [3, 2, 1, 0]
    ]
    """
    start = (0, len(rna)-1)
    initial_fold = [None for _ in range(len(rna))]
    optimal_folds = backtrackh(opt_solutions, start, initial_fold)

    # make sure that they are unique
    optimal_folds2 = set( tuple(fold) for fold in optimal_folds )
    assert len(optimal_folds) == len(optimal_folds2)

    # convert None to self match
    for fold in optimal_folds:
        for i in range(len(fold)):
            if fold[i] is None:
                fold[i] = i

    return optimal_folds

def foldcompare(f1, f2):
    """
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

assert foldcompare([0, 1], [0, 1]) == 0
assert foldcompare([1, 0], [1, 0]) == 0
assert foldcompare([1, 0, 2], [2, 1, 0]) == 2
assert foldcompare([1, 0, 2], [0, 2, 1]) == 2

def backtrackh(opt_solutions: list, curr: tuple, folding: list):
    currsize = curr[1] - curr[0] + 1
    if currsize <= 1:
        return [folding.copy()]

    foldings = []
    events = opt_solutions[curr[0]][curr[1]]
    for ename, e in events:
        if ename == 'L':
            newcurr = curr[0]+1, curr[1]
            for result in backtrackh(opt_solutions, newcurr, folding):
                foldings.append(result.copy())
        elif ename == 'M':
            foldstart, foldend = e # foldstart == curr[0]+1
            newcurr1 = curr[0]+1, foldend-1
            newcurr2 = foldend+1, curr[1]

            folding[foldstart] = foldend
            folding[foldend] = foldstart

            for partial_folding in backtrackh(opt_solutions, newcurr1, folding):
                for result in backtrackh(opt_solutions, newcurr2, partial_folding):
                    foldings.append(result.copy())

            folding[foldstart] = None
            folding[foldend] = None

        else:
            raise NotImplementedError()
    return foldings

def distancevector(fold: list, allfolds: list):
    vector = [0] * len(fold)*2

    for other in allfolds:
        dist = foldcompare(fold, other)
        vector[dist] += 1

    # trim the end so that it looks clean
    while vector and vector[-1] == 0:
        vector.pop()
    return vector


if __name__ == '__main__':
    test_suboptimal_fold()
    test_max_distance_dp()
    print("All tests passed!")