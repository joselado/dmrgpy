import sympy 

class MultiOperator(sympy.Symbol):
    """Multioperator class"""
    def __init__(self,name):
        super().__init__(name) # initialize



#
#def obj2MO(a,name="multioperator"):
#    """
#    Convert an input in a multioperator
#    """
#    if isinstance(a, collections.Iterable): # if it is a tuple
#        mo = 0 # initialize
#        mo = MultiOperator(name=name) # create object
#        for ia in a:
#            mo.add_operator(ia[0],ia[1])
#        return mo
#    elif type(a)==MultiOperator:
#        a.name = name
#        return a
#    else: raise # unidentified input
#


