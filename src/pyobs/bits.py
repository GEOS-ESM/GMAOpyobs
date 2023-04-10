"""
  Bit manipulation.
"""
class BITS(object):
    """
    This class allow you to access the individuals bits of a number
    with a protocol similar to list indexing and slicing.
    """
    def __init__(self,value=0):
        self._d = value

    def __getitem__(self, index):
        if isinstance(index, slice):
            start = index.start
            end   = index.stop
            mask = 2**(end - start) -1
            return (self._d >> start) & mask
        elif isinstance(index, int):
            return (self._d >> index) & 1 

    def __setitem__(self,index,value):
        if isinstance(index, slice):
            start = index.start
            end   = index.end
            mask = 2**(end - start) -1
            value = (value & mask) << start
            mask = mask << start
            self._d = (self._d & ~mask) | value
            return (self._d >> start) & mask
        elif isinstance(index, int):
            value    = (value&1)<<index
            mask     = (1)<<index
            self._d  = (self._d & ~mask) | value



    def __int__(self):
        return self._d

        
