import argparse
from core import run_Nussinov_functions, run_Zuker_functions
from utils import RNA_BASES, convert_DBN_to_folding

# N stands for Nussinov algorithm
# Z stands for Zuker algorithm
ALG_OPTIONS = ['N', 'Z']


def check_rna(rna):
    for c in rna:
        if c not in RNA_BASES:
            raise ValueError("Input RNA string is not valid")

def get_distance_statistics(rna_file: str, dbn_file: str, alg: str):
    try:
        if not rna_file:
            raise ValueError("Need to specify the RNA file name")
        if not dbn_file:
            raise ValueError("Need to specify the dbn file name")
        with open(rna_file, 'r') as f:
            rna = f.read().replace(' ', '').replace('\n', '').upper()
        with open(dbn_file, 'r') as f:
            dbn = f.read().replace(' ', '').replace('\n', '')

        if len(rna) == 0:
            raise ValueError("Input RNA string cannot be empty")
        if len(rna) != len(dbn):
            raise ValueError("The dbn string length does not match the RNA string length")
        check_rna(rna)
        folding = convert_DBN_to_folding(rna, dbn)
        if alg not in ALG_OPTIONS:
            raise ValueError("Choice of algorithm is not valid")           
    except ValueError as e:
        print(e)
        return
    else:
        if alg == 'N':
            run_Nussinov_functions(rna, folding)
        if alg == "Z":
            run_Zuker_functions(rna, folding)
        return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("rna_file", type=str, help="filename for RNA string")
    parser.add_argument("dbn_file", type=str, help="filename for dbn string that represents a valid folding")
    parser.add_argument("alg", type=str, choices=['N', 'Z'], \
        help="select an algorithm to get optimal foldings: N for Nussinov, Z for Zuker")
    args = parser.parse_args()
    get_distance_statistics(args.rna_file, args.dbn_file, args.alg)
