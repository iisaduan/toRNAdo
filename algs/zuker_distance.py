from algs.zuker_backtrack import EntryTable
from algs.distance_vector import V

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

    def update(self, val: int, e: E = None):
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
        self.maxdist_W =   [[DistanceEntry(None) for j in range(N)] for i in range(N+1)]
        self.maxdist_V =   [[DistanceEntry(None) for j in range(N)] for i in range(N+1)]
        self.maxdist_WM =  [[DistanceEntry(None) for j in range(N)] for i in range(N+1)]
        self.maxdist_WM2 = [[DistanceEntry(None) for j in range(N)] for i in range(N+1)]
        # create DP tables for calculating distance vectors
        self.vec_W =   [[None for j in range(N)] for i in range(N+1)]
        self.vec_V =   [[None for j in range(N)] for i in range(N+1)]
        self.vec_WM =  [[None for j in range(N)] for i in range(N+1)]
        self.vec_WM2 = [[None for j in range(N)] for i in range(N+1)]
    
    def fill_distance_table(self):
        N = len(self.seq)
        for i in range(0, N):
            for tablename in ["W", "V", "WM", "WM2"]:
                self.compute_max_distance(tablename, i, i)
                self.compute_max_distance(tablename, i+1, i)

        for i in reversed(range(0, N)):
            for j in range(i+1, N):
                self.compute_max_distance("V", i, j)
                self.compute_max_distance("W", i, j) # relies on V[i][j]
                self.compute_max_distance("WM", i, j) # relies on V[i][j]
                self.compute_max_distance("WM2", i, j)
    
    def fill_vector_table(self):
        N = len(self.seq)
        for i in range(0, N):
            for tablename in ["W", "V", "WM", "WM2"]:
                self.compute_vec(tablename, i, i)
                self.compute_vec(tablename, i+1, i)

        for i in reversed(range(0, N)):
            for j in range(i+1, N):
                self.compute_vec("V", i, j)
                self.compute_vec("W", i, j) # relies on V[i][j]
                self.compute_vec("WM", i, j) # relies on V[i][j]
                self.compute_vec("WM2", i, j)
    
    def get_opt_table(self, tablename):
        if tablename == 'W':
            return self.opt_W
        elif tablename == 'V':
            return self.opt_V
        elif tablename == 'WM2':
            return self.opt_WM2
        elif tablename == 'WM':
            return self.opt_WM

    def get_max_dist_table(self, tablename):
        if tablename == 'W':
            return self.maxdist_W
        elif tablename == 'V':
            return self.maxdist_V
        elif tablename == 'WM2':
            return self.maxdist_WM2
        elif tablename == 'WM':
            return self.maxdist_WM
    
    def get_vec_table(self, tablename):
        if tablename == 'W':
            return self.vec_W
        elif tablename == 'V':
            return self.vec_V
        elif tablename == 'WM2':
            return self.vec_WM2
        elif tablename == 'WM':
            return self.vec_WM
    
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
    
    def compute_max_distance(self, tablename, i, j):
        table = self.get_max_dist_table(tablename)
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
                distance += self.get_max_dist_table(newtableName)[new_i][new_j].val
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
        max_distance, vec = None, None
        if self.maxdist_W[0][N-1]:
            max_distance = self.maxdist_W[0][N-1].val
        if self.vec_W[0][N-1]:
            vec = self.vec_W[0][N-1].v
        return max_distance, vec
    
    def get_one_max_dist_solution_h(self, i, j, table, solution):
        if j <= i:
            return
        # choose the first solution
        _, match, recurseL = table[i][j].eList[0]
        if match:
            i,j = match
            solution[i] = j
            solution[j] = i
        for tablename, i, j in recurseL:
            self.get_one_max_dist_solution_h(i, j, \
                self.get_max_dist_table(tablename), solution)
    
    def get_one_max_dist_solution(self):
        """get an arbitrary solution"""
        N = len(self.seq)
        solution = [i for i in range(N)]
        self.get_one_max_dist_solution_h(0, N-1, self.maxdist_W, solution)
        return solution
    
    def compute_vec(self, tablename, i, j):
        table = self.get_vec_table(tablename)
        opt_table = self.get_opt_table(tablename)
        # base
        if i == j or i-1 == j:
            if tablename == "W":
                # one solution of distance 0
                # and no solution of other distance
                table[i][j] = V([1])
            else:
                # no solution of any distance 
                # since there is no valid solution
                table[i][j] = V([0])
            return
        opt_choices = opt_table[i][j].eList
        vec = V([0])
        for choice in opt_choices:
            new_vec = V([1])
            _, match, recurseL = choice
            # total number of basepairs in folding[i:j+1]
            shift = self.contained_counts[i][j]
            for newtableName, new_i, new_j in recurseL:
                new_vec = new_vec ^ self.get_vec_table(newtableName)[new_i][new_j]
                # subtract the basepairs that are considered in recursive calls
                shift -= self.contained_counts[new_i][new_j]
            if match:
                # (i,j) is the match in the current optimal choice
                # check the pairing of i in the given folding
                pair_index = self.folding[i]
                # i is unpaired in folding
                if i == pair_index:
                    shift += 1
                # i is paired to the same index
                elif pair_index == j:
                    shift -= 1
                # i is paired to a different index
                elif i < pair_index < j:
                    shift += 1
                # i is paired to an index not in the range
                else:
                    shift += 1
            new_vec = (new_vec >> shift)
            vec += new_vec
        table[i][j] = vec





