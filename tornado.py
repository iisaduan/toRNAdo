"The command line tool for distance profiling"

import argparse
from core import run_Nussinov_functions, run_Zuker_functions
from algs.utils import RNA_BASES, convert_DBN_to_folding

# N stands for Nussinov algorithm
# Z stands for Zuker algorithm
ALG_OPTIONS = ['N', 'Z']


def check_rna(rna: str):
    """Check if the given RNA file contains valid characters.

    Args:
        rna (str): a string

    Raises:
        ValueError: if a character in the string is not an RNA base, then this
            error is raised.
    """
    for c in rna:
        if c not in RNA_BASES:
            raise ValueError("Input RNA string is not valid")

def get_distance_statistics(rna_file: str, dbn_file: str, alg: str, internal_loop_size: int, \
    display: bool, plot: bool, plot_file: str):
    """Provide the distance statistics according to the command line options that
    the user selects.

    Args:
        rna_file (str): a filepath that stores the RNA sequence
        dbn_file (str): a filepath that stores a folding for the RNA sequence
            in dot-bracket notation
        alg (str): the algorithm the user chooses for RNA folding
        internal_loop_size (int): the maximum size for an internal loop in the
            Zuker algorithm.
        display (bool): whether to display a folding of max distance from the
            base folding.
        plot (bool): whether to plot the histogram for the distance vector
        plot_file (str): the filepath for storing the histogram plot.

    Raises:
        ValueError: raised if rna_file is not provided
        ValueError: raised if RNA sequence is empty
        ValueError: raised if the provided folding is not a valid folding for
            the provided RNA sequence
        ValueError: raised if the choice of algorithm is invalid
    """
    try:
        # Check the validity of inputs
        if not rna_file:
            raise ValueError("Need to specify the RNA file name")
        else:
            with open(rna_file, 'r') as f:
                rna = f.read().replace(' ', '').replace('\n', '').upper()
            if len(rna) == 0:
                raise ValueError("Input RNA string cannot be empty")
            check_rna(rna)
        if not dbn_file:
            folding = None
        else:
            with open(dbn_file, 'r') as f:
                dbn = f.read().replace(' ', '').replace('\n', '')
            if len(dbn) != len(rna):
                raise ValueError("The dbn string length does not match the RNA string length")
            folding = convert_DBN_to_folding(rna, dbn)
        if alg not in ALG_OPTIONS:
            raise ValueError("Choice of algorithm is not valid")           
    except ValueError as e:
        print(e)
        return
    else:
        if alg == 'N':
            run_Nussinov_functions(rna, folding, display, plot, plot_file)
        if alg == "Z":
            run_Zuker_functions(rna, folding, internal_loop_size, display, plot, plot_file)
        return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("alg", type=str, choices=['N', 'Z'], \
        help="select an algorithm to get optimal foldings: N for Nussinov, Z for Zuker")
    parser.add_argument("rna_file", type=str, help="filename for RNA string")
    parser.add_argument("-s","--internal_loop_size", type=int, help="bound the maximum size for an internal loop in Zuker's algorithm for speedup")
    parser.add_argument("-dbn", "--dbn_file", type=str, help="filename for dbn string that represents a valid folding")
    parser.add_argument("-d", "--display_max_distance", action="store_true", help="display a folding of max distance from the given folding")
    parser.add_argument("-p", "--plot", action="store_true", help="plot the histogram for the distance vector")
    parser.add_argument("-f", "--plot_file", type=str, help="filename for storing the histogram plot")
    args = parser.parse_args()
    get_distance_statistics(args.rna_file, args.dbn_file, args.alg, args.internal_loop_size, args.display_max_distance, args.plot, args.plot_file)
