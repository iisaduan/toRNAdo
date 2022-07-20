import argparse
from core import run_Nussinov_functions, run_Zuker_functions
from algs.utils import RNA_BASES, convert_DBN_to_folding

# N stands for Nussinov algorithm
# Z stands for Zuker algorithm
ALG_OPTIONS = ['N', 'Z']


def check_rna(rna):
    for c in rna:
        if c not in RNA_BASES:
            raise ValueError("Input RNA string is not valid")

def get_distance_statistics(rna_file: str, dbn_file: str, alg: str, d: bool, h: bool, plot_file: str):
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
            run_Nussinov_functions(rna, folding, d, h, plot_file)
        if alg == "Z":
            run_Zuker_functions(rna, folding, d, h, plot_file)
        return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("rna_file", type=str, help="filename for RNA string")
    parser.add_argument("alg", type=str, choices=['N', 'Z'], \
        help="select an algorithm to get optimal foldings: N for Nussinov, Z for Zuker")
    parser.add_argument("--dbn_file", type=str, help="filename for dbn string that represents a valid folding")
    parser.add_argument("--d", action="store_true", help="display a folding of max distance from the given folding")
    parser.add_argument("--h", action="store_true", help="create the histogram for the distance vector")
    parser.add_argument("--f", type=str, help="filename for storing the histogram plot")
    args = parser.parse_args()
    get_distance_statistics(args.rna_file, args.dbn_file, args.alg, args.d, args.h, args.f)
