"""Provides common functionalities to both Nussinov and Zuker algorithms"""

import matplotlib.pyplot as plt

RNA_BASES = ['A', 'U', 'C', 'G']
MATCHING_PAIRS = [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')]

def is_base_pair(b1: str, b2: str) -> bool:
    """Determine whether two bases are complementary

    Args:
        b1 (str): a base in RNA sequence
        b2 (str): another base in RNA sequence

    Raises:
        ValueError: the base b1 or b2 is not an RNA base

    Returns:
        bool: True if b1 and b2 are complementary bases. False otherwise.
    """

    if b1 not in RNA_BASES or b2 not in RNA_BASES:
        raise ValueError("Input is not an RNA string")
    return (b1, b2) in MATCHING_PAIRS

def convert_DBN_to_folding(rna: str, dbn: str) -> list:
    """Convert a dot-bracket string to a list that represents a folding

    Args:
        rna (str): an RNA sequence
        dbn (str): a folding represented by the dot-bracket notation

    Raises:
        ValueError: dot-bracket notation is invalid

    Returns:
        list: the same folding represented by a list where folding[i] stores the
            index that i is matched to. folding[i] = i if i is unmatched.
    """
    dbn_length = len(dbn)
    folding = [i for i in range(dbn_length)]
    count = 0
    index_array = [-1 for i in range(dbn_length)]
    for index in range(dbn_length):
        if dbn[index] == "(":
            count += 1
            index_array[count] = index
        elif dbn[index] == ")":
            if count == 0:
                raise ValueError("The dbn representation is not valid. Extra ) detected")
            pair_index = index_array[count]
            if not is_base_pair(rna[index], rna[pair_index]):
                raise ValueError("The dbn representation is not valid. RNA bases that do not form a basepair are paired")
            assert pair_index != -1
            folding[index] = pair_index
            folding[pair_index] = index
            count -= 1
        else:
            pass
    if count != 0:
        raise ValueError("The dbn representation is not valid. Extra ( detected")
    return folding

def convert_folding_to_DBN(folding: list) -> str:
    """Convert a list that represents a folding to a dot-bracket string that
    represents the same folding.

    Args:
        folding (list): list that represents a folding

    Returns:
        str: the same folding in dot-bracket notation
    """
    dbn = ""
    for i,j in enumerate(folding):
        if i < j:
            dbn += "("
        elif i > j:
            dbn += ")"
        else:
            dbn += "."
    return dbn

def display_histogram(vec: list, alg: str, plot_file: str, folding_given: bool):
    """Save a histogram plot for the distance vector

    Args:
        vec (list): the distance vector to be displayed
        alg (str): 'N' if we are comparing against all Nussinov-optimal
            solutions and 'Z' if we are comparing against all Zuker-optimal
            solutions
        plot_file (str): the filepath for saving the plot
        folding_given (bool): whether the folding is provided by the user or
            randomly selected
    """
    num_bins = len(vec)
    plt.bar(list(range(num_bins)), list(vec))
    # Force y-axis to use scientific notation
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    if alg == 'N': 
        alg_name = "Nussinov"
    elif alg == 'Z':
        alg_name = "Zuker"
    
    plt.ylabel(f"Number of {alg_name}-Optimal Solutions")
    plt.xlabel("Distance")
    if folding_given:
        plt.title(f"Distances from the Given Folding", y=1.08, fontsize=18)
    else:
        plt.title(f"Distances from a Selected Random Folding", y=1.08, fontsize=18)
    plt.savefig(plot_file, bbox_inches='tight')
    plt.clf()


if __name__ == '__main__':
    assert convert_DBN_to_folding("ACG", "...") == [0,1,2]
    assert convert_DBN_to_folding("AU", "()") == [1,0]
    assert convert_DBN_to_folding("ACAUU", "(.())") == [4,1,3,2,0]
    assert convert_folding_to_DBN([]) == ""
    assert convert_folding_to_DBN([0,1,2]) == "..."
    assert convert_folding_to_DBN([1,0]) == "()"
    assert convert_folding_to_DBN([4,1,3,2,0]) == "(.())"