from zuker_backtrack import EntryTable

ThingToRecurseOn = list[tuple[str, int, int]]
E = tuple[str, tuple[int, int], list[ThingToRecurseOn]]

class DistanceEntry:
    def __init__(self,
        val: int,
        eList: list[E] = None
    ):
        if eList is None:
            eList = []
        self.val = val
        self.eList = eList

    def update(self, val: float, e: E = None):
        # e: tuple['Event', tuple[i,j]/None, list[ThingToRecurseOn]]
        # ThingToRecurseOn: tuple['TableName', i, j]
        if val > self.val:
            self.val = val
            self.eList = [e]
        elif val == self.val and self.val != float('-inf'):
            self.eList.append(e)
DistanceEntryTable = list[list[DistanceEntry]]

class DistanceSolver:
    def __init__(self,
    seq: str,
    folding: list,
    opt_W: EntryTable,
    opt_V: EntryTable,
    opt_WM: EntryTable,
    opt_WM2: EntryTable
    ):
        self.seq = seq
        self.folding = folding
        # pass in DP tables from running Zuker's algorithm
        self.opt_W = opt_W
        self.opt_V = opt_V
        self.opt_WM = opt_WM
        self.opt_WM2 = opt_WM2
        # precompute counts of base pairs that are contained in certain range
        # in the provided folding
        self.contained_counts = self.count_contained_basepairs()
        # create DP tables for calculating max distance
        N = len(self.seq)
        self.W =   [[DistanceEntry(None) for j in range(N)] for i in range(N+1)]
        self.V =   [[DistanceEntry(None) for j in range(N)] for i in range(N+1)]
        self.WM =  [[DistanceEntry(None) for j in range(N)] for i in range(N+1)]
        self.WM2 = [[DistanceEntry(None) for j in range(N)] for i in range(N+1)]
    
    def fill_table(self):
        N = len(self.seq)
        for i in range(0, N):
            self.compute("W", i, i)
            self.compute("V", i, i)
            self.compute("WM", i, i)
            self.compute("WM2", i, i)
            
            self.compute("W", i+1, i)
            self.compute("V", i+1, i)
            self.compute("WM", i+1, i)
            self.compute("WM2", i+1, i)

        for i in reversed(range(0, N)):
            for j in range(i+1, N):
                self.compute("V", i, j)
                self.compute("W", i, j) # relies on V[i][j]
                self.compute("WM", i, j) # relies on V[i][j]
                self.compute("WM2", i, j)
    
    def get_table(self, tablename):
        if tablename == 'W':
            return self.W
        elif tablename == 'V':
            return self.V
        elif tablename == 'WM2':
            return self.WM2
        elif tablename == 'WM':
            return self.WM
    
    def get_opt_table(self, tablename):
        if tablename == 'W':
            return self.opt_W
        elif tablename == 'V':
            return self.opt_V
        elif tablename == 'WM2':
            return self.opt_WM2
        elif tablename == 'WM':
            return self.opt_WM
    
    def count_contained_basepairs(self) -> list:
        """
        Return a 2D list dp.
        dp[i][j] = the number of basepairs in self.folding that are entirely contained in rna[i,j+1]
        """
        N = len(self.seq)
        dp = [[0 for j in range(N)] for i in range(N+1)]
        for start in reversed(range(0, N)):
            for end in range(0, N):
                pair_index = self.folding[start]
                is_start_pair_contained = start < pair_index <= end
                dp[start][end] = dp[start+1][end] + int(is_start_pair_contained)
        return dp
    
    def compute(self, tablename, i, j):
        table = self.get_table(tablename)
        opt_table = self.get_opt_table(tablename)
        # base
        if i == j or i-1 == j:
            if tablename == "W":
                table[i][j] = DistanceEntry(0)
            else:
                table[i][j] = DistanceEntry(float('-inf'))
            return
        entry = DistanceEntry(float('-inf')) 
        opt_choices = opt_table[i][j].eList
        for choice in opt_choices:
            distance: int = 0
            _, match, recurseL = choice
            # total number of basepairs in folding[i:j+1]
            distance = self.contained_counts[i][j]
            for newtableName, new_i, new_j in recurseL:
                distance += self.get_table(newtableName)[new_i][new_j].val
                # subtract the basepairs that are considered in recursive calls
                distance -= self.contained_counts[new_i][new_j]
            if match:
                # (i,j) is the match in the current optimal choice
                # check the pairing of i in the given folding
                pair_index = self.folding[i]
                # i is unpaired in folding
                if i == pair_index:
                    distance += 1
                # i is paired to the same index
                elif pair_index == j:
                    distance -= 1
                # i is paired to a different index
                elif i < pair_index < j:
                    distance += 1
                # i is paired to an index not in the range
                else:
                    distance += 1
            entry.update(distance, choice)
        table[i][j] = entry
    
    def solve(self):
        N = len(self.seq)
        return self.W[0][N-1].val
    
    def get_one_solution_h(self, i, j, table, solution):
        if j <= i:
            return
        # choose the first solution
        _, match, recurseL = table[i][j].eList[0]
        if match:
            i,j = match
            solution[i] = j
            solution[j] = i
        for tablename, i, j in recurseL:
            self.get_one_solution_h(i, j, self.get_table(tablename), solution)
    
    # TODO: allow choosing a randomly selected solution/median solution etc.
    def get_one_solution(self):
        """get an arbitrary solution"""
        N = len(self.seq)
        solution = [i for i in range(N)]
        self.get_one_solution_h(0, N-1, self.W, solution)
        return solution



