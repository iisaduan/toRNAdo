from itertools import zip_longest
class V:
    """
    Distance Vector
    """

    def __init__(self, v=None):
        """
        Create an empty vector, 
        or create a vector from the given list v.
        """
        if v is None:
            self.v = []
        else:
            self.v = v
    
    def __rshift__(self, n):
        """
        Shift the vector right by n.
        V([1]) >> 2 == V([0,0,1])
        """
        v = [0] * n
        v.extend(self.v)
        return V(v)
    
    def __add__(self, other: 'V'):
        """
        Add two vectors
        V([1]) + V([1,3]) == V([2,3])
        """
        v = []
        for x, y in zip_longest(self.v, other.v, fillvalue=0):
            v.append(x+y)
        return V(v)
    
    def __mul__(self, scalar: int):
        """
        Multiply by scalar
        V([1]) * 2 == V([2])
        """
        return V([ scalar * x for x in self.v])
    
    def __xor__(self, other: 'V'):
        """
        Vector coalesce
        V([1,2]) ^ V([2,3]) == V([2,7,6])
        """
        vec = V()
        for shift, coeff in enumerate(other.v):
            vec += (self * coeff) >> shift
        return vec
    
    def __eq__(self, other: 'V'):
        return self.v == other.v

assert V() * 2  == V()
assert V([1]) * 2 == V([2])
assert V([1]) >> 2 == V([0, 0, 1])
assert V([1]) ^ V([1]) == V([1])
assert V([1]) ^ V([2]) == V([2])
assert V([2]) ^ V([1]) == V([2])
assert V([1,2]) ^ V([2,3]) == V([2,7,6])