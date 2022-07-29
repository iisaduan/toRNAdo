from functools import reduce
from algs.utils import convert_folding_to_DBN, display_histogram
from algs.zuker_backtrack import *
from algs.zuker_distance import DistanceSolver
from algs.nussinov import *

PLOT_FILENAME = "plot"

def run_Zuker_functions(rna: str, folding: list = None, internal_loop_size: int = None,\
    display: bool = False, plot: bool = False, plot_file: str = PLOT_FILENAME):

    print("-----Running Zuker's algorithm-----")
    # run Zuker algorithm to get optimal folds
    solver = Solver(seq = rna, internal_loop_size = internal_loop_size)
    solver.fill_table()
    min_energy = solver.solve()
    print("The minimum energy is: ", min_energy)

    # choose a folding to compare against
    afold = None
    folding_given = False
    if folding:
        # use the given folding
        afold = folding
        folding_given = True
    else:
        # select an arbitrary folding
        afold = solver.get_one_solution()

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
    if display:
        # display an optimal fold that is the most distant from the given fold
        max_distance_fold = distance_solver.get_one_max_dist_solution()
        print("One of the farthest foldings is: ")
        print(convert_folding_to_DBN(max_distance_fold))
    if plot:
        # plot the histogram for the distance vector
        display_histogram(distance_vector, 'Z', plot_file, folding_given)
    print()


def run_Nussinov_functions(rna: str, folding: list = None, display: bool = False,\
    plot: bool = False, plot_file: str = PLOT_FILENAME):
    print("-----Running Nussinov's algorithm-----")
    # run Nussinov algorithm to get optimal folds
    opt_val, opt_solutions = nussinov_dp(rna)
    print("The maximum number of matched pairs is: ", opt_val)

    # choose a folding to compare against
    afold = None
    folding_given = False
    if folding:
        # use the given folding
        afold = folding
        folding_given = True
    else:
        # select an arbitrary folding
        afold = construct_one_opt_solution(rna, opt_solutions)

    max_distance, max_distance_solutions = max_distance_dp(afold, opt_solutions)
    distance_vector = max_distance_vec(afold, opt_solutions).v
    print("Max distance is: ", max_distance)
    print("The distance vector is: ", distance_vector)
    num_folds_dp = reduce((lambda x, y: x + y), distance_vector)
    print("The number of optimal folds is: ", num_folds_dp)
    if display:
        # display an optimal fold that is the most distant from the given fold
        max_distance_fold = construct_one_opt_solution(rna, max_distance_solutions)
        print("One of the farthest foldings is: ")
        print(convert_folding_to_DBN(max_distance_fold))
    if plot:
        # plot the histogram for the distance vector
        display_histogram(distance_vector, 'N', plot_file, folding_given)
    print()


if __name__ == '__main__':
    rna = "CUUCCCAGGUAACAAACCAACCAACUUUCGAUCUCUUGUAGAUCUGUUCUCUAAACGCUUCCCAGGUAACAAACCAACCAACUUUCGAUCUCUUGUAGAUCUGUUCUCUAAACG"
    run_Nussinov_functions(rna, display=True, plot=True, plot_file="plot_Nussinov")
    run_Zuker_functions(rna, display=True, plot=True, plot_file="plot_Zuker", internal_loop_size=None)
    run_Zuker_functions(rna, display=True, plot=True, plot_file="plot_Zuker_30", internal_loop_size=30)
