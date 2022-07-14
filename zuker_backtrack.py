from typing import Callable

ThingToRecurseOn = list[tuple[str, int, int]]
E = tuple[str, tuple[int, int], list[ThingToRecurseOn]]

class Entry:
    def __init__(self,
        val: float,
        eList: list[E] = None
    ):
        if eList is None:
            eList = []
        self.val = val
        self.eList = eList

    def update(self, val: float, e: E = None):
        # e: tuple['Event', tuple[i,j]/None, list[ThingToRecurseOn]]
        # ThingToRecurseOn: tuple['TableName', i, j]
        if val < self.val:
            self.val = val
            self.eList = [e]
        elif val == self.val and self.val != float('inf'):
            self.eList.append(e)


class Solver:
    def __init__(self, seq: str, 
        m: int,
        a: int,
        b: int,
        c: int,
        eH: Callable[[int, int], int],
        eS: Callable[[int, int], int],
        eL: Callable[[int, int], int]
    ):
        self.seq = seq
        self.m = m # minimum seperation for a pair
        self.a = a
        self.b = b
        self.c = c
        self.eH = eH
        self.eS = eS
        self.eL = eL

        # create table
        N = len(self.seq)
        self.W =   [[Entry(None) for j in range(N)] for i in range(N+1)]
        self.V =   [[Entry(None) for j in range(N)] for i in range(N+1)]
        self.WM =  [[Entry(None) for j in range(N)] for i in range(N+1)]
        self.WM2 = [[Entry(None) for j in range(N)] for i in range(N+1)]
    
    def fill_table(self):
        N = len(self.seq)
        for i in range(0, N):
            self.compute_W(i, i)
            self.compute_V(i, i)
            self.compute_WM(i, i)
            self.compute_WM2(i, i)
            
            self.compute_W(i+1, i)
            self.compute_V(i+1, i)
            self.compute_WM(i+1, i)
            self.compute_WM2(i+1, i)

        for i in reversed(range(0, N)):
            for j in range(i+1, N):
                self.compute_V(i, j)
                self.compute_W(i, j) # relies on V[i][j]
                self.compute_WM(i, j) # relies on V[i][j]
                self.compute_WM2(i, j)
                

    def match(self, i: int, j: int):
        return (self.seq[i], self.seq[j]) in [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')]

    def compute_W(self, i: int, j: int):
        # base
        if i == j or i-1 == j:
            self.W[i][j] = Entry(0)
            return
        entry = Entry(float('inf'))
        # case 1: i is not paired
        r1 = self.W[i+1][j].val
        breadcrumb = ('L', None, [('W', i+1, j)])
        entry.update(r1, breadcrumb)

        # case 2: i is paired with k
        for k in range(i, j+1):
            r2 = self.V[i][k].val + self.W[k+1][j].val
            breadcrumb = ('M', (i,k), [('V', i, k), ('W', k+1, j)])
            entry.update(r2, breadcrumb)
        self.W[i][j] = entry
    
    def compute_V(self, i: int, j: int):
        # base
        if j-i <= self.m or not self.match(i,j):
            self.V[i][j] = Entry(float('inf'))
            return
        entry = Entry(float('inf'))
        # Case 1: Hairpin loop
        r1 = self.eH(i, j)
        breadcrumb = ('H', None, [])
        entry.update(r1, breadcrumb)
        # Case 2: Stacking loop
        r2 = self.V[i+1][j-1].val + self.eS(i,j)
        breadcrumb = ('S', (i+1,j-1), [('V', i+1, j-1)])
        entry.update(r2, breadcrumb)
        # Case 3: Internal loop
        # TODO: bottle neck, maybe use heuristic to limit interior loop size
        for i2 in range(i+1, j):
            for j2 in range(i2+1, j):
                if i2 == i+1 and j2 == j-1:
                    # stacking loop has already been considered in Case 2
                    continue
                r3 = self.V[i2][j2].val+ self.eL(i,j,i2,j2)
                breadcrumb = ('I', (i2,j2), [('V', i2, j2)])
                entry.update(r3, breadcrumb)     
        # Case 4: Multiloop
        r4 = self.a + self.WM2[i+1][j-1].val
        breadcrumb = ('ML', None, [('WM2', i+1, j-1)])
        entry.update(r4, breadcrumb)

        self.V[i][j] = entry
    
    def compute_WM2(self, i: int, j: int):
        # base
        if i == j or i-1 == j:
            self.WM2[i][j] = Entry(float('inf'))
            return
        entry = Entry(float('inf'))
        # Case 1: i is not paired
        r1 = self.WM2[i+1][j].val + self.c
        breadcrumb = ('L', None, [('WM2', i+1, j)])
        entry.update(r1, breadcrumb)
        # Case 2: i is paired with k
        # Note: i cannot pair with j
        for k in range(i+1, j):
            r2 = self.V[i][k].val + self.b + self.WM[k+1][j].val
            breadcrumb = ('M', (i,k), [('V', i, k), ("WM", k+1, j)])
            entry.update(r2, breadcrumb)
        self.WM2[i][j] = entry
    
    def compute_WM(self, i: int, j: int):
        # base
        if i == j or i-1 == j:
            self.WM[i][j] = Entry(float('inf'))
            return
        entry = Entry(float('inf'))
        # Case 1: i is not paired
        r1 = self.WM[i+1][j].val + self.c
        breadcrumb = ('L', None, [('WM', i+1, j)])
        entry.update(r1, breadcrumb)
        # Case 2: i is paired with k and there are more loops in [k+1,j]
        for k in range(i+1, j+1):
            r2 = self.V[i][k].val + self.b + self.WM[k+1][j].val
            breadcrumb = ('M-WM', (i,k), [('V', i, k), ("WM", k+1, j)])
            entry.update(r2, breadcrumb)
        # Case 3: i is paired with k and there are no more loops in [k+1,j]
        # Note: i can pair with j
        for k in range(i+1, j+1):
            r3 = self.V[i][k].val + self.b + (j-k)*self.c
            breadcrumb = ('M-None', (i,k), [('V', i, k)])
            entry.update(r3, breadcrumb)
        self.WM[i][j] = entry

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
        for tName, i, j in recurseL:
            if tName == "W":
                self.get_one_solution_h(i, j, self.W, solution)
            elif tName == "V":
                self.get_one_solution_h(i, j, self.V, solution)
            elif tName == "WM2":
                self.get_one_solution_h(i, j, self.WM2, solution)
            elif tName == "WM":
                self.get_one_solution_h(i, j, self.WM, solution)
    
    # TODO: allow choosing a randomly selected solution/median solution etc.
    def get_one_solution(self):
        """get an arbitrary solution"""
        N = len(self.seq)
        solution = [i for i in range(N)]
        self.get_one_solution_h(0, N-1, self.W, solution)
        return solution

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

class DistanceCalc:
    def __init__(self,
    seq: str,
    folding: list,
    opt_W: DistanceEntryTable,
    opt_V: DistanceEntryTable,
    opt_WM: DistanceEntryTable,
    opt_WM2: DistanceEntryTable
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
            self.compute_W(i, i)
            self.compute_V(i, i)
            self.compute_WM(i, i)
            self.compute_WM2(i, i)
            
            self.compute_W(i+1, i)
            self.compute_V(i+1, i)
            self.compute_WM(i+1, i)
            self.compute_WM2(i+1, i)

        for i in reversed(range(0, N)):
            for j in range(i+1, N):
                self.compute_V(i, j)
                self.compute_W(i, j) # relies on V[i][j]
                self.compute_WM(i, j) # relies on V[i][j]
                self.compute_WM2(i, j)
    
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
    
    def compute_W(self, i: int, j: int):
        # base
        if i == j or i-1 == j:
            # the only optimal solution has no basepair
            self.W[i][j] = DistanceEntry(0)
            return
        entry = DistanceEntry(float('-inf'))
        pair_index = self.folding[i]
        pair_in_range = int(i < pair_index <= j)  
        opt_choices = self.opt_W[i][j].eList

        for choice in opt_choices:
            r: int = 0
            breadcrumb: tuple = None
            # case 1: i is not paired
            if choice[0] == 'L':
                r = self.W[i+1][j].val + (1 if pair_in_range else 0)
                breadcrumb = ('L', None, [('W', i+1, j)])
            # case 2: i is paired with k
            elif choice[0] == "M":
                # get the basepair matching in the optimal choice
                _, k = choice[1]  
                # count the number of basepairs that are not accounted for by recursive calls
                # specially, since we are recursing on ('V', i, k), ('W', k+1, j),
                # we need to count the # of basepairs in folding that are contained in [i+1, j]
                # but not in [i+1, k-1] or [k+1, j]
                count = self.contained_counts[i+1][j] \
                        - self.contained_counts[i+1][k-1] \
                        - self.contained_counts[k+1][j]
                # compare the basepair (i, k) with the pairing of i in folding 
                r = self.V[i][k].val + self.W[k+1][j].val \
                    + 2 * int(pair_index != k) \
                    - (1 if not pair_in_range else 0) \
                    + count
                breadcrumb = ('M', (i,k), [('V', i, k), ('W', k+1, j)])
            entry.update(r, breadcrumb)
        self.W[i][j] = entry
    
    def compute_V(self, i: int, j: int):
        # base
        if i == j or i-1 == j:
            # no optimal solution exists
            self.V[i][j] = DistanceEntry(float('-inf'))
            return
        entry = DistanceEntry(float('-inf'))  
        opt_choices = self.opt_V[i][j].eList
        # TODO: maybe no need to break by cases,
        # just look at whether there is a matching in the breadcrumb
        # and the recursive calls
        for choice in opt_choices:
            r: int = 0
            breadcrumb: tuple = None
            # Case 1: Hairpin loop
            if choice[0] == 'H':
                # distance is the number of basepairs contained in folding in [i+1,j-1]
                r = self.contained_counts[i+1][j-1]
                breadcrumb = ('H', None, [])
            # Case 2: Stacking loop
            elif choice[0] == "S":
                # compare the basepair (i+1, j-1) with the pairing of i+1 in folding
                pair_index = self.folding[i+1]
                pair_in_range = int(i+1 < pair_index <= j-1)  
                count = self.contained_counts[i+2][j-1] \
                        - self.contained_counts[i+2][j-2]
                r = self.V[i+1][j-1].val \
                    + 2 * int(pair_index != j-1) \
                    - (1 if not pair_in_range else 0) \
                    + count
                breadcrumb = ('S', (i+1,j-1), [('V', i+1, j-1)])
            # Case 3: Internal loop
            elif choice[0] == "I":
                # get the optimal interior match
                (i2, j2) = choice[1]
                # compare the basepair (i2,j2) with the pairing of i2 in folding
                # add the number of pairs in folding in [i+1,j-1] that are not
                # entirely contained in [i2+1, j2-1] to the difference
                pair_index = self.folding[i2]
                pair_in_range = int(i2 < pair_index <= j2)  
                # count the number of differences that are missed in recursive cases
                count = self.contained_counts[i+1,j-1] - self.contained_counts[i2+1, j2-1]
                r = self.V[i2][j2].val \
                    + 2 * int(pair_index != j2) \
                    - (1 if not pair_in_range else 0) \
                    + count
                breadcrumb = ('I', (i2,j2), [('V', i2, j2)])
            # Case 4: Multiloop
            elif choice[0] == "ML":
                r = self.WM2[i+1][j-1].val
                breadcrumb = ('ML', None, [('WM2', i+1, j-1)])
            entry.update(r, breadcrumb)
        self.V[i][j] = entry  
   
    def compute_WM2(self, i: int, j: int):
        # base
        if i == j or i-1 == j:
            self.WM2[i][j] = DistanceEntry(float('-inf'))
            return
        entry = DistanceEntry(float('-inf'))
        pair_index = self.folding[i]
        pair_in_range = int(i < pair_index <= j)  
        opt_choices = self.opt_WM2[i][j].eList
        for choice in opt_choices:
            r: int = 0
            breadcrumb: tuple = None
            # Case 1: i is not paired
            if choice[0] == 'L':
                r = self.WM2[i+1][j].val + (1 if pair_in_range else 0)
                breadcrumb = ('L', None, [('WM2', i+1, j)])
            # Case 2: i is paired with k
            elif choice[0] == "M":
                # get the basepair matching in the optimal choice
                _, k = choice[1]
                # count the number of differences that are missed in recursive cases
                count = self.contained_counts[i+1][j] \
                        - self.contained_counts[i+1][k-1] \
                        - self.contained_counts[k+1][j]
                r = self.V[i][k].val + self.WM[k+1][j].val \
                    + 2 * int(pair_index != k) \
                    - (1 if not pair_in_range else 0) \
                    + count
                breadcrumb = ('M', (i,k), [('V', i, k), ("WM", k+1, j)])
            entry.update(r, breadcrumb)
        self.WM2[i][j] = entry
    
    def compute_WM(self, i: int, j: int):
        # base
        if i == j or i-1 == j:
            self.WM[i][j] = DistanceEntry(float('-inf'))
            return
        entry = DistanceEntry(float('-inf'))
        pair_index = self.folding[i]
        pair_in_range = int(i < pair_index <= j)  
        opt_choices = self.opt_WM[i][j].eList
        for choice in opt_choices:
            r: int = 0
            breadcrumb: tuple = None
            # Case 1: i is not paired
            if choice[0] == 'L':
                r = self.WM[i+1][j].val + (1 if pair_in_range else 0)
                breadcrumb = ('L', None, [('WM', i+1, j)])
            # Case 2: i is paired with k and there are more loops in [k+1,j]
            elif choice[0] == "M-WM":
                # get the basepair matching in the optimal choice
                _, k = choice[1]  
                # count the number of differences that are missed in recursive cases
                count = self.contained_counts[i+1][j] \
                        - self.contained_counts[i+1][k-1] \
                        - self.contained_counts[k+1][j]
                r = self.V[i][k].val + self.WM[k+1][j].val \
                    + 2 * int(pair_index != k) \
                    - (1 if not pair_in_range else 0) \
                    + count
                breadcrumb = ('M-WM', (i,k), [('V', i, k), ("WM", k+1, j)])
            # Case 3: i is paired with k and there are no more loops in [k+1,j]
            elif choice[0] == "M-None":
                # get the basepair matching in the optimal choice
                _, k = choice[1]
                # count the number of differences that are missed in other cases
                count = self.contained_counts[i+1][j] - self.contained_counts[i+1][k-1]
                r = self.V[i][k].val \
                    + 2 * int(pair_index != k) \
                    - (1 if not pair_in_range else 0) \
                    + count
                breadcrumb = ('M-None', (i,k), [('V', i, k)])
            entry.update(r, breadcrumb)
        self.WM[i][j] = entry

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
        for tName, i, j in recurseL:
            if tName == "W":
                self.get_one_solution_h(i, j, self.W, solution)
            elif tName == "V":
                self.get_one_solution_h(i, j, self.V, solution)
            elif tName == "WM2":
                self.get_one_solution_h(i, j, self.WM2, solution)
            elif tName == "WM":
                self.get_one_solution_h(i, j, self.WM, solution)
    
    # TODO: allow choosing a randomly selected solution/median solution etc.
    def get_one_solution(self):
        """get an arbitrary solution"""
        N = len(self.seq)
        solution = [i for i in range(N)]
        self.get_one_solution_h(0, N-1, self.W, solution)
        return solution


def debug(table, table2):
    print()
    print(' ' * 10, end='')
    for i in range(len(rna)):
        print(f'{rna[i]:>10}', end='')
    print()
    print(' ' * 10, end='')
    for i in range(len(rna)):
        print(f'{i:>10}', end='')
    print()
    for i in range(len(rna)):
        print(f'{i:>10}', end="")
        for j in range(len(rna)):
            if not j - i >= -1:
                print(' '*10, end="")
                continue
            result1 = table[i][j].val
            result2 = table2[i][j] if table2 else ""
            s = str(result1) + '|' + str(result2)
            print(f'{s:>10}', end="")
        print()

if __name__ == '__main__':
    m = 0
    a, b, c = -1, -1, -1
    eH = lambda i, j: i-j
    eS = lambda i, j: i-j
    eL = lambda i, j, i2, j2: 0

    rna = "ACGU"

    new_solver = Solver(rna, m, a, b, c, eH, eS, eL)
    new_solver.fill_table()
    print("minimum energy is: ", new_solver.solve())


    folding = new_solver.get_one_solution()
    print("given folding is: ", folding)
    distance_calc = DistanceCalc(rna, folding,
    new_solver.W, new_solver.V, new_solver.WM, new_solver.WM2)
    distance_calc.fill_table()
    print("max distance is: ", distance_calc.solve())
    folding_2 = distance_calc.get_one_solution()
    print("farthest optimal folding is: ", folding_2)

    # debug(distance_calc.W, None)
    # debug(distance_calc.V, None)
    # debug(distance_calc.WM, None)
    # debug(distance_calc.WM2, None)
