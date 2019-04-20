
import collections

# this class allows to define operators of the form
# A_0@A_1@....
# for mpscpp.x


class MultiOperator():
    """
    Object to deal with multioperators in mpscpp.x
    """
    def __init__(self,name="multioperator"): # do nothing
        self.op = [] # empty list
        self.name = name
    def add_operator(self,name,i):
        """Store operator"""
        self.op.append([name,str(i)])
    def get_dict(self):
        """Return the dictionary to be used in tasks.in"""
        d = dict()
        d[self.name+"_n"] = len(self.op)
        for i in range(len(self.op)): # loop
            o = self.op[i]
            name0 = self.name+"_operator_"+str(i)+"_name"
            name1 = self.name+"_operator_"+str(i)+"_site"
            d[name0] = o[0]
            d[name1] = o[1]
        return d



def obj2MO(a,name="multioperator"):
    """
    Convert an input in a multioperator
    """
    if isinstance(a, collections.Iterable): # if it is a tuple
        mo = MultiOperator(name=name) # create object
        for ia in a:
            mo.add_operator(ia[0],ia[1])
        return mo
    elif type(a)==MultiOperator:
        a.name = name
        return a
    else: raise # unidentified input




