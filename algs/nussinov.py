from algs.distance_vector import V
from algs.utils import is_base_pair

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
    if N == 0: 
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

def construct_one_opt_solution(rna, opt_solutions):
    """
    opt_solutions: breadcrumbs for either the Nussinov DP algorithm or the Max distance algorithm
    Returns a folding from the set of optimal solutions
    """
    N = len(rna)
    constructed_solution = [i for i in range(N)]
    construct_one_opt_solution_helper(0, N-1, opt_solutions, constructed_solution)
    return constructed_solution

def count_matched_basepairs(folding) -> list:
    """
    Return a 2D list dp.
    dp[i][j] = the number of basepairs in folding that are entirely contained in rna[i,j+1]
    """
    N = len(folding)
    dp = [[0 for j in range(N)] for i in range(N+1)]
    for start in reversed(range(0, N)):
        for end in range(0, N):
            pair_index = folding[start]
            is_start_pair_contained = start < pair_index <= end
            dp[start][end] = dp[start+1][end] + int(is_start_pair_contained)
    return dp        

def max_distance_dp(folding, opt_solutions):
    N = len(folding)
    matched_pairs_counts = count_matched_basepairs(folding)
    # dp[i][j] = the maximum distance that any optimal folding of rna[i:j+1] has from
    #            the given folding if we only consider basepairs contained in rna[i:j+1].
    #            The distance metric used here is symmetric set difference of the basepairs.
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
            pair_in_range = int(start < pair_index <= end)
            max_distance = 0
            max_distance_choices = []
            opt_choices = opt_solutions[start][end]

            current_distance = 0
            current_choice = []
            for choice in opt_choices:
                if choice[0] == "L":
                    current_distance = dp[start+1][end] + (1 if pair_in_range else 0)
                    current_choice = ("L", start)
                elif choice[0] == "M":
                    # get the basepair matching in the optimal choice
                    _, opt_pair_index = choice[1]

                    current_distance = dp[start+1][opt_pair_index-1] \
                                        + dp[opt_pair_index+1][end] \
                                        + matched_pairs_counts[start+1][end] \
                                        - matched_pairs_counts[start+1][opt_pair_index-1] \
                                        - matched_pairs_counts[opt_pair_index+1][end]

                    if opt_pair_index != pair_index:
                        if pair_in_range:
                            current_distance += 2
                        else:
                            current_distance += 1

                    current_choice = ("M", (start, opt_pair_index))
                if current_distance > max_distance:
                    max_distance = current_distance
                    max_distance_choices = [current_choice]
                elif current_distance == max_distance:
                    max_distance_choices.append(current_choice)

            dp[start][end] = max_distance
            breadcrumbs[start][end] = max_distance_choices
    return dp[0][N-1], breadcrumbs


def max_distance_vec(folding, opt_solutions):
    N = len(folding)
    matched_pairs_counts = count_matched_basepairs(folding)
    dp = [[None for j in range(N)] for i in range(N+1)]

    for i in range(0, N):
        dp[i][i] = V([1])
        dp[i+1][i] = V([1])

    for start in reversed(range(0, N)):
        for end in range(start+1, N):
            pair_index = folding[start]
            pair_in_range = int(start < pair_index <= end)
            vec = V()
            opt_choices = opt_solutions[start][end]
            for choice in opt_choices:
                if choice[0] == "L":
                    vec += dp[start+1][end] >> (1 if pair_in_range else 0)
                elif choice[0] == "M":
                    _, opt_pair_index = choice[1]    
                    vec += (dp[start+1][opt_pair_index-1]
                            ^ dp[opt_pair_index+1][end]) \
                            >> (2 * int(pair_index != opt_pair_index)
                                - (1 if not pair_in_range else 0)
                                + matched_pairs_counts[start+1][end] \
                                - matched_pairs_counts[start+1][opt_pair_index-1] \
                                - matched_pairs_counts[opt_pair_index+1][end])

            dp[start][end] = vec
    return dp[0][N-1]