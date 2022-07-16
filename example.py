from functools import reduce
from zuker_backtrack import *
from zuker_distance import DistanceSolver
from test_zuker_backtrack import gen_random_fold
from nussinov import *

def run_Zuker_functions(rna, folding=None):
    print("-----Running Zuker algorithm-----")
    # run Zuker algorithm to get optimal folds
    solver = Solver(rna)
    solver.fill_table()
    min_energy = solver.solve()
    print("The minimum energy is: ", min_energy)

    # choose a folding to compare against
    afold = folding if folding is not None else gen_random_fold(rna)

    # run distance solver to get max distance and vector
    distance_solver = DistanceSolver(rna, afold, solver.W, solver.V,
                                solver.WM, solver.WM2)
    distance_solver.fill_distance_table()
    distance_solver.fill_vector_table()
    max_distance = distance_solver.solve()[0]   
    distance_vector = distance_solver.solve()[1]
    print("Max distance is: ", max_distance)
    print("The distance vector is: ", distance_vector)
    num_folds_dp = reduce((lambda x, y: x + y), distance_vector)
    print("The number of optimal folds is: ", num_folds_dp)
    print()


def run_Nussinov_functions(rna, folding=None):
    print("-----Running Nussinov algorithm-----")
    # run Nussinov algorithm to get optimal folds
    opt_val, opt_solutions = nussinov_dp(rna)
    print("The maximum number of folds is: ", opt_val)

    # choose a folding to compare against
    afold = folding if folding is not None else gen_random_fold(rna)

    max_distance, _ = max_distance_dp(afold, opt_solutions)
    distance_vector = max_distance_vec(afold, opt_solutions).v
    print("Max distance is: ", max_distance)
    print("The distance vector is: ", distance_vector)
    num_folds_dp = reduce((lambda x, y: x + y), distance_vector)
    print("The number of optimal folds is: ", num_folds_dp)
    print()


if __name__ == '__main__':
    rna = "CUUCCCAGGUAACAAACCAACCAACUUUCGAUCUCUUGUAGAUCUGUUCUCUAAACGAACUUUAAAAUCUGUGUGGCUGUCACUCGGCUGCAUGCUUAGUGCACUCACGCAGUAUAAUUA"
    run_Nussinov_functions(rna)
    run_Zuker_functions(rna)
