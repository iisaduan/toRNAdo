
def convert_DBN_to_folding(dbn: str) -> list:
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
            pair_index = index_array[count]
            assert pair_index != -1
            folding[index] = pair_index
            folding[pair_index] = index
            count -= 1
        else:
            pass
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

assert convert_DBN_to_folding("") == []
assert convert_DBN_to_folding("...") == [0,1,2]
assert convert_DBN_to_folding("()") == [1,0]
assert convert_DBN_to_folding("(.())") == [4,1,3,2,0]
assert convert_folding_to_DBN([]) == ""
assert convert_folding_to_DBN([0,1,2]) == "..."
assert convert_folding_to_DBN([1,0]) == "()"
assert convert_folding_to_DBN([4,1,3,2,0]) == "(.())"