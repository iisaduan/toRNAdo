matchings_pairs = [('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')]
def is_base_pair(b1, b2):
    return (b1, b2) in matchings_pairs

def nussinov_dp(rna: str) -> tuple[int, list[list[int]]]:
    """
    Given rna sequence, compute the maximum number of paired bases.

    Return an int result and 2D table opt_solutions.
    result = the maximum number of paired bases (int) for rna
    opt_solutions[i][j] = a list of choices for i that produce optimal 
                       parings for rna[i:j+1]. A choice is a tuple
                       ('L', i) for a loss or ('M', (i, k)) for a pairing.
    
    nussinov_dp("CGCG") -> (2, opt_solutions)
    opt_solutions:
        [0][1]: [('M', (0, 1))]
        [0][2]: [('L', 0), ('M', (0, 1))]
        [0][3]: [('M', (0, 1)), ('M', (0, 3))]
        [1][2]: [('M', (1, 2))]
        [1][3]: [('L', 1), ('M', (1, 2))]
        [2][3]: [('M', (2, 3))]
    }
    """
    N = len(rna)
    if len(rna) == 0: 
        return 0
    
    # dp[i][j] = the maximum number of paired bases (int) for rna rna[i:j+1]
    # Note that dp[i+1][i] denotes the empty string base case, so we have to store
    # the cell dp[N][N-1]. This means we have to store (N+1)*N matrix to store that cell.
    dp = [[None for j in range(N)] for i in range(N+1)]
    opt_solutions = [[None for j in range(N)] for i in range(N+1)]

    # base cases dp[i+1][i] = 0 (size 0), dp[i][i] = 0 (size 1)
    for i in range(0, N):
        dp[i][i] = 0
        opt_solutions[i][i] = []
        dp[i+1][i] = 0
        opt_solutions[i+1][i] = []
    
    for start in reversed(range(0, N)):
        for end in range(start+1, N):
            # initialize maximum number of folds and optimal choices to the lose it case
            max_folds = dp[start+1][end]
            opt_choices_start = [("L", start)]
            # now consider the different options of pairing the nucleobase at start 
            # with another nucleobase at certain index in the string
            for altpair in range(start+1, end+1):
                if is_base_pair(rna[start], rna[altpair]):
                    num_folds = 1 + dp[start+1][altpair-1] + dp[altpair+1][end]
                    if num_folds > max_folds:
                        max_folds = num_folds
                        opt_choices_start = [("M", (start, altpair))]
                    elif num_folds == max_folds:
                        opt_choices_start.append(("M", (start, altpair)))
                    else:
                        continue
            dp[start][end] = max_folds
            opt_solutions[start][end] = opt_choices_start
    return dp[0][N-1], opt_solutions

def construct_one_opt_solution_helper(start, end, opt_solutions, solution):
    """return a list where the ith index contains the index that i is paired with.
    If i is unpaired, we say i is paired with itself."""
    if end <= start:
        return
    else:
        choice = opt_solutions[start][end][0]
        if choice[0] == "L":
            construct_one_opt_solution_helper(start+1, end, opt_solutions, solution)
        elif choice[0] == "M":
            (index, pair_index) = choice[1]
            if index!= start:
                raise ValueError()
            solution[start] = pair_index
            solution[pair_index] = start
            construct_one_opt_solution_helper(start+1, pair_index-1, opt_solutions, solution)
            construct_one_opt_solution_helper(pair_index+1, end, opt_solutions, solution)

def construct_one_opt_solution(rna, opt_solutions=None):
    N = len(rna)
    constructed_solution = [i for i in range(N)]
    if opt_solutions is None:
        _, opt_solutions = nussinov_dp(rna)
    construct_one_opt_solution_helper(0, N-1, opt_solutions, constructed_solution)
    return constructed_solution

def convert_DBN_to_folding(dbn: str) -> list:
    dbn_length = len(dbn)
    folding = [i for i in range(dbn_length)]
    count = 0
    index_array = [-1 for i in range(dbn_length)]
    for index in range(dbn_length):
        if dbn[index] == "(":
            count += 1
            index_array[count] = index
        elif dbn[index] == ")":
            pair_index = index_array[count]
            assert pair_index != -1
            folding[index] = pair_index
            folding[pair_index] = index
            count -= 1
        else:
            pass
    return folding

def min_similarities_dp(folding, opt_solutions):
    N = len(folding)
    # dp[i][j] = the minimum number of paired bases that an optimal folding 
    #            shares with the given folding in the range [i:j+1]
    # Note that dp[i+1][i] denotes the empty string base case, so we have to store
    # the cell dp[N][N-1]. This means we have to store (N+1)*N matrix to store that cell.
    dp = [[None for j in range(N)] for i in range(N+1)]
    breadcrumbs = [[None for j in range(N)] for i in range(N+1)]

    # base cases dp[i+1][i] = 0 (size 0), dp[i][i] = 0 (size 1)
    for i in range(0, N):
        dp[i][i] = 0
        breadcrumbs[i][i] = []
        dp[i+1][i] = 0
        breadcrumbs[i+1][i] = []

    for start in reversed(range(0, N)):
        for end in range(start+1, N):
            pair_index = folding[start]
            min_similarities = N
            min_similarities_choices = []
            # TODO: in the DP solution, we need to check that this dictionary entry is not empty
            opt_choices = opt_solutions[start][end]

            current_similarities = 0
            current_choice = []
            for choice in opt_choices:
                if choice[0] == "L":
                    current_similarities = dp[start+1][end]
                    current_choice = ("L", start)
                elif choice[0] == "M":
                    # get the basepair matching in the optimal choice
                    _, opt_pair_index = choice[1]  
                    current_similarities = dp[start+1][opt_pair_index-1] \
                                            + dp[opt_pair_index+1][end] \
                                            + int(pair_index == opt_pair_index)
                                        
                    current_choice = ("M", (start, opt_pair_index))

                if current_similarities < min_similarities:
                    min_similarities = current_similarities
                    min_similarities_choices = [current_choice]
                elif current_similarities == min_similarities:
                    min_similarities_choices.append(current_choice)

            dp[start][end] = min_similarities
            breadcrumbs[start][end] = min_similarities_choices
    
    return dp[0][N-1], breadcrumbs

from itertools import zip_longest
class V:
    """
    Similarity Vector, for use with min_similarity algorithm.
    """

    def __init__(self, v=None):
        """
        Create an empty vector, 
        or create a vector from the given list v.
        """
        if v is None:
            self.v = []
        else:
            self.v = v
    
    def __rshift__(self, n):
        """
        Shift the vector right by n.
        V([1]) >> 2 == V([0,0,1])
        """
        v = [0] * n
        v.extend(self.v)
        return V(v)
    
    def __add__(self, other: 'V'):
        """
        Add two vectors
        V([1]) + V([1,3]) == V([2,3])
        """
        v = []
        for x, y in zip_longest(self.v, other.v, fillvalue=0):
            v.append(x+y)
        return V(v)
    
    def __mul__(self, scalar: int):
        """
        Multiply by scalar
        V([1]) * 2 == V([2])
        """
        return V([ scalar * x for x in self.v])
    
    def __xor__(self, other: 'V'):
        """
        Vector coalesce
        V([1,2]) ^ V([2,3]) == V([2,7,6])
        """
        vec = V()
        for shift, coeff in enumerate(other.v):
            vec += (self * coeff) >> shift
        return vec
    
    def __eq__(self, other: 'V'):
        return self.v == other.v
    
    def to_diff_vec(self, maxfold: int):
        """
        Convert from similarities vector (self) to
        difference vector (list) using the formula
        diff = 2*maxfold - 2*sim.
        """
        v = self.v
        diffv = [0] * len(v) * 2
        for i, coeff in enumerate(v):
            newi = 2*maxfold - 2*i
            diffv[newi] = coeff
        
        while diffv and diffv[-1] == 0:
            diffv.pop()

        return diffv


assert V() * 2  == V()
assert V([1]) * 2 == V([2])
assert V([1]) >> 2 == V([0, 0, 1])
assert V([1]) ^ V([1]) == V([1])
assert V([1]) ^ V([2]) == V([2])
assert V([2]) ^ V([1]) == V([2])
assert V([1,2]) ^ V([2,3]) == V([2,7,6])

def min_similarities_vec(folding, opt_solutions):
    N = len(folding)

    dp = [[None for j in range(N)] for i in range(N+1)]
    for i in range(0, N):
        dp[i][i] = V([1])
        dp[i+1][i] = V([1])

    for start in reversed(range(0, N)):
        for end in range(start+1, N):
            pair_index = folding[start]
            opt_choices = opt_solutions[start][end]

            vec = V()
            for choice in opt_choices:
                if choice[0] == "L":
                    vec += dp[start+1][end]
                elif choice[0] == "M":
                    _, opt_pair_index = choice[1]  
                    vec += (dp[start+1][opt_pair_index-1] ^ dp[opt_pair_index+1][end]) \
                            >> int(pair_index == opt_pair_index)

            dp[start][end] = vec
    
    return dp[0][N-1]



