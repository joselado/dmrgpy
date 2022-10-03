# library for immutable operators,
# they can act over a wavefunction, but they do not have
# algebra

from ..mps import MPS

class StaticOperator():
    def __init__(self,MO,MBO):
        """Init, takes as input a multioperator and the MBO"""
        self.MBO = MBO # store the many-body object
        self.SO = generate_SO(MO,MBO) # generate the static operator
    def __mul__(self,v):
        if type(v)==MPS: # input is an MPS
            return pure_applyoperator_dmrg(self.MBO,self.SO,v)
        else: raise



def pure_applyoperator_dmrg(self,A,wf):
    """Apply a pure operator to a many body wavefunction"""
    self.execute(lambda: wf.write()) # write WF
    task = {"pureapplyoperator":"true",
            "pureapplyoperator_wf0":wf.name,
            "pureapplyoperator_operator":"pureapplyoperator_operator.mpo",
            "pureapplyoperator_wf1":"pureapplyoperator_wf1.mps",
            }
    self.execute(lambda: open(self.path+"/pureapplyoperator_operator.mpo","wb").write(A))
    self.task = task
    self.execute(lambda : self.run()) # run calculation
    return MPS(self,name="pureapplyoperator_wf1.mps").copy() # copy


def generate_SO(A,MBO):
    """Generate a static many-body object"""
    task = {"gen_pureoperator":"true",
            "gen_pureoperator_operator_in":"gen_pureoperator_operator.in",
            "gen_pureoperator_operator_out":"gen_pureoperator_operator.mpo",
            }
    MBO.execute(lambda: A.write(name="pureapplyoperator_operator.in"))
    MBO.task = task
    MBO.execute(lambda : MBO.run()) # run calculation
    return open(MBO.path+"/gen_pureoperator_operator.mpo","rb").read() 


