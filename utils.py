
RNA_BASES = ['A', 'U', 'C', 'G']
MATCHING_PAIRS = [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')]

def is_base_pair(b1, b2):
    if b1 not in RNA_BASES or b2 not in RNA_BASES:
        raise ValueError("Input is not an RNA string")
    return (b1, b2) in MATCHING_PAIRS

def convert_DBN_to_folding(rna: str, dbn: str) -> list:
    """Convert a dot-bracket string to a folding"""
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
    """Convert a folding to a dot-bracket string"""
    dbn = ""
    for i,j in enumerate(folding):
        if i < j:
            dbn += "("
        elif i > j:
            dbn += ")"
        else:
            dbn += "."
    return dbn

if __name__ == '__main__':
    assert convert_DBN_to_folding("ACG", "...") == [0,1,2]
    assert convert_DBN_to_folding("AU", "()") == [1,0]
    assert convert_DBN_to_folding("ACAUG", "(.())") == [4,1,3,2,0]
    assert convert_folding_to_DBN([]) == ""
    assert convert_folding_to_DBN([0,1,2]) == "..."
    assert convert_folding_to_DBN([1,0]) == "()"
    assert convert_folding_to_DBN([4,1,3,2,0]) == "(.())"