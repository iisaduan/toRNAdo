"""Implements the Zuker DP algorithm"""

from decimal import Decimal, getcontext
from typing import Callable
from algs.seqfold.fold import _hairpin
from algs.seqfold.fold import _stack
from algs.seqfold.fold import _internal_loop
from algs.seqfold.rna import RNA_ENERGIES
from algs.utils import is_base_pair

RecursiveCall = tuple[str, int, int]
Traceback = list[RecursiveCall]
# Example Choice: ('S', (i, j), [('V', i+1, j-1)])
# interpreted as: Stacking loop, matches i and j, recurse on V[i+1][j-1]
Choice = tuple[str, tuple[int, int], Traceback]

# set precision of decimal.Decimal class
getcontext().prec = 2

class Entry:
    """An entry in the Zuker algorithm DP tables
    """
    def __init__(self,
        val: Decimal,
        eList: list[Choice] = None
    ):
        """Initialize an entry.

        Args:
            val (Decimal): optimal value for this entry
            eList (list[Choice], optional): a list of optimal choices that
                gives the same optimal value for this entry. Defaults to None.
        """
        if eList is None:
            eList = []
        self.val = val
        self.eList = eList

    def update(self, val: Decimal, e: Choice = None):
        """Update an entry with a new value and the corresponding choice.

        Args:
            val (Decimal): new value for this entry
            e (Choice, optional): associated breadcrumb. Defaults to None.
        """
        if val < self.val:
            self.val = val
            self.eList = [e]
        elif val == self.val and self.val != Decimal('Infinity'):
            self.eList.append(e)

EntryTable = list[list[Entry]]

class Solver:
    """Implements the Zuker DP algorithm
    """
    def __init__(self, seq: str,
        m: int = None,
        a: float = None,
        b: float = None,
        c: float = None,
        eH: Callable[[str, int, int], float] = None,
        eS: Callable[[str, int, int], float] = None,
        eL: Callable[[str, int, int, int, int], float] = None,
        internal_loop_size: int = None,
    ):
        """Initialize the Solver with parameters and functions used
        in the algorithm and energy calculations for different loops.

        Args:
            seq (str): RNA sequence
            m (int, optional): minimum size of a closed loop. For example,
                a loop (i,j) must satisfy j-i >= m. Defaults to None.
            a (float, optional): parameter for multiloop energy calculation. Defaults to None.
            b (float, optional): parameter for multiloop energy calculation. Defaults to None.
            c (float, optional): parameter for multiloop energy calculation. Defaults to None.
            eH (Callable[[str, int, int], float], optional): energy function for Hairpin loops. Defaults to None.
            eS (Callable[[str, int, int], float], optional): energy function for Stacking loops. Defaults to None.
            eL (Callable[[str, int, int, int, int], float], optional): energy function for Internal loops. Defaults to None.
            internal_loop_size (int, optional): maximum size of an internal loop.
                For example, we require an internal loop (i,j,i',j') to satisfy
                j'-i'+1 <= internal_loop_size. Defaults to None.
        """
        self.seq = seq
        self.internal_loop_size = internal_loop_size
        # if any of the following parameters is None, we import the value from
        # the seqfold module copied from https://github.com/Lattice-Automation/seqfold.
        self.m = m if m is not None else 4
        self.a = a if a is not None else Decimal(RNA_ENERGIES.MULTIBRANCH[0])
        self.b = b if b is not None else Decimal(RNA_ENERGIES.MULTIBRANCH[1])
        self.c = c if c is not None else Decimal(RNA_ENERGIES.MULTIBRANCH[2])
        self.eH = eH if eH is not None else \
            lambda s, i, j: Decimal(_hairpin(s, i, j, float('37.0'), RNA_ENERGIES))
        self.eS = eS if eS is not None else \
            lambda s, i, j: Decimal(_stack(self.seq, i, i+1, j, j-1, float('37.0'), RNA_ENERGIES))
        self.eL = eL if eL is not None else \
            lambda s, i, j, i2, j2: Decimal(_internal_loop(self.seq, i, i2, j, j2, float('37.0'), RNA_ENERGIES))

        # create tables
        N = len(self.seq)
        self.W =   [[Entry(None) for j in range(N)] for i in range(N+1)]
        self.V =   [[Entry(None) for j in range(N)] for i in range(N+1)]
        self.WM =  [[Entry(None) for j in range(N)] for i in range(N+1)]
        self.WM2 = [[Entry(None) for j in range(N)] for i in range(N+1)]      
    
    def fill_table(self):
        """Fill the tables in the order of V, W, WM, WM2
        """
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
                

    def match(self, i: int, j: int) -> bool:
        """Determine if the bases at index i and j are complementary

        Args:
            i (int): index into self.seq
            j (int): another index into self.seq

        Returns:
            bool: True if the bases at the indices are complementary and False otherwise.
        """
        return is_base_pair(self.seq[i], self.seq[j])

    def compute_W(self, i: int, j: int):
        """Compute W[i][j].

        Args:
            i (int): row index
            j (int): column index
        """
        # base
        if i == j or i-1 == j:
            self.W[i][j] = Entry(Decimal('0'))
            return
        entry = Entry(Decimal('Infinity'))
        # case 1: i is not paired
        r1 = self.W[i+1][j].val
        breadcrumb = ('L', None, [('W', i+1, j)])
        entry.update(r1, breadcrumb)

        # case 2: i is paired with k
        for k in range(i, j+1):
            r2 = self.V[i][k].val + self.W[k+1][j].val
            breadcrumb = ('M', None, [('V', i, k), ('W', k+1, j)])
            entry.update(r2, breadcrumb)
        self.W[i][j] = entry
    
    def compute_V(self, i: int, j: int):
        """Compute V[i][j].

        Args:
            i (int): row index
            j (int): column index
        """
        # base
        if j-i < self.m or not self.match(i,j):
            self.V[i][j] = Entry(Decimal('Infinity'))
            return
        entry = Entry(Decimal('Infinity'))
        # Case 1: Hairpin loop
        r1 = self.eH(self.seq, i, j)
        breadcrumb = ('H', (i, j), [])
        entry.update(r1, breadcrumb)
        # Case 2: Stacking loop
        if self.match(i+1, j-1):
            r2 = self.V[i+1][j-1].val + self.eS(self.seq, i, j)
            breadcrumb = ('S', (i, j), [('V', i+1, j-1)])
            entry.update(r2, breadcrumb)
        # Case 3: Internal loop
        # bottle neck, if internal loop size is specified,
        # we use it to bound the loop size for speedup
        for i2 in range(i+1, j):
            # determine the range for j2 of the internal loop (i, j, i_2, j_2)
            # if internal loop size is unbounded, we only require j_2 < j
            end = j
            if self.internal_loop_size:
                end = i2 + self.internal_loop_size
            for j2 in range(i2+1, end):
                if j2 > j-1:
                    continue
                if i2 == i+1 and j2 == j-1:
                    # stacking loop has already been considered in Case 2
                    continue
                if self.match(i2, j2):
                    r3 = self.V[i2][j2].val + self.eL(self.seq, i, j, i2, j2)
                    breadcrumb = ('I', (i, j), [('V', i2, j2)])
                    entry.update(r3, breadcrumb)     
        # Case 4: Multiloop
        r4 = self.a + self.WM2[i+1][j-1].val
        breadcrumb = ('ML', (i, j), [('WM2', i+1, j-1)])
        entry.update(r4, breadcrumb)

        self.V[i][j] = entry
    
    def compute_WM2(self, i: int, j: int):
        """Compute WM2[i][j].

        Args:
            i (int): row index
            j (int): column index
        """
        # base
        if i == j or i-1 == j:
            self.WM2[i][j] = Entry(Decimal('Infinity'))
            return
        entry = Entry(Decimal('Infinity'))
        # Case 1: i is not paired
        r1 = self.WM2[i+1][j].val + self.c
        breadcrumb = ('L', None, [('WM2', i+1, j)])
        entry.update(r1, breadcrumb)
        # Case 2: i is paired with k
        # Note: i cannot pair with j
        for k in range(i+1, j):
            r2 = self.V[i][k].val + self.b + self.WM[k+1][j].val
            breadcrumb = ('M', None, [('V', i, k), ("WM", k+1, j)])
            entry.update(r2, breadcrumb)
        self.WM2[i][j] = entry
    
    def compute_WM(self, i: int, j: int):
        """Compute WM[i][j].

        Args:
            i (int): row index
            j (int): column index
        """
        # base
        if i == j or i-1 == j:
            self.WM[i][j] = Entry(Decimal('Infinity'))
            return
        entry = Entry(Decimal('Infinity'))
        # Case 1: i is not paired
        r1 = self.WM[i+1][j].val + self.c
        breadcrumb = ('L', None, [('WM', i+1, j)])
        entry.update(r1, breadcrumb)
        # Case 2: i is paired with k and there are more loops in [k+1,j]
        for k in range(i+1, j+1):
            r2 = self.V[i][k].val + self.b + self.WM[k+1][j].val
            breadcrumb = ('M-WM', None, [('V', i, k), ("WM", k+1, j)])
            entry.update(r2, breadcrumb)
        # Case 3: i is paired with k and there are no more loops in [k+1,j]
        # Note: i can pair with j
        for k in range(i+1, j+1):
            r3 = self.V[i][k].val + self.b + (j-k)*self.c
            breadcrumb = ('M-None', None, [('V', i, k)])
            entry.update(r3, breadcrumb)
        self.WM[i][j] = entry

    def solve(self) -> Decimal:
        """Return the minimum energy of a folding of self.seq.

        Returns:
            Decimal: the minimum energy of a folding of self.seq
        """
        N = len(self.seq)
        return self.W[0][N-1].val

    def get_table(self, tablename: str) -> EntryTable:
        """Get the table object corresponding to the tablename.

        Args:
            tablename (str): the name of the table

        Returns:
            EntryTable: the table object corresponding to the name
        """
        if tablename == 'W':
            return self.W
        elif tablename == 'V':
            return self.V
        elif tablename == 'WM2':
            return self.WM2
        elif tablename == 'WM':
            return self.WM
    
    def get_one_solution_h(self, i: int, j: int, table: EntryTable, solution: list):
        """Helper function for getting one optimal solution for the Zuker algorithm.

        Args:
            i (int): start index
            j (int): end index
            table (EntryTable): the current table we are recursing into
            solution (list): stores the resulting optimal solution
        """
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
    
    def get_one_solution(self) -> list:
        """Get an arbitrary optimal solution for the Zuker algorithm.

        Returns:
            list: the resulting optimal solution/folding
        """
        N = len(self.seq)
        solution = [i for i in range(N)]
        self.get_one_solution_h(0, N-1, self.W, solution)
        return solution



if __name__ == '__main__':
    # Example of how to use this solver
    m = 0
    a, b, c = -1, -1, -1
    eH = lambda rna, i, j: i-j
    eS = lambda rna, i, j: i-j
    eL = lambda rna, i, j, i2, j2: 0

    rna = "ACGU"
    new_solver = Solver(rna, m, a, b, c, eH, eS, eL)
    new_solver.fill_table()
    print("minimum energy is: ", new_solver.solve())
    folding = new_solver.get_one_solution()
    print("given folding is: ", folding)
