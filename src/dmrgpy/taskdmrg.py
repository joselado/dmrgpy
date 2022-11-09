from __future__ import print_function
import numpy as np


def setup_task(self,mode="GS",task=dict()):
  select_task(self,mode=mode)
  for key in task: self.task[key] = task[key] # additional arguments
  write_tasks(self) # write the tasks



def select_task(self,mode="GS"):
  """Setup the sweep parameters"""
  task = dict() # dictionary
  if mode=="GS": # default mode
    task["GS"] = "true"
  elif mode=="excited": # default mode
    task["excited"] = "true"
  elif mode=="correlator": # default mode
    task["GS"] = "true"
    task["correlator"] = "true"
  elif mode=="entropy": # default mode
    task["GS"] = "true"
    task["entropy"] = "true"
  elif mode=="dos": # default mode
    task["dos"] = "true"
  elif mode=="spismj": # default mode
    task["spismj"] = "true"
  elif mode=="overlap": # default mode
    task["overlap"] = "true"
  elif mode=="dynamical_correlator": # default mode
    task["dynamical_correlator"] = "true"
  elif mode=="time_evolution": # default mode
    task["time_evolution"] = "true"
#    task["orthogonal_kpm"] = "true"
  else: raise
  self.task = task # initialize




def write_tasks(self):
  fo = open("tasks.in","w")
  fo.write("tasks\n{\n")
  #
  if self.use_ampo_hamiltonian: 
      fo.write(" use_ampo_hamiltonian = true\n")
  if self.gs_from_file and self.wf0 is not None: 
      self.execute(self.wf0.write) # write WF
      fo.write(" gs_from_file = true\n")
      fo.write(" starting_file_gs = "+self.wf0.name+"\n") # starting WF
      fo.write(" skip_dmrg_gs = "+obj2str(self.skip_dmrg_gs)+"\n") # starting WF
  else: fo.write(" gs_from_file = false\n")
  for key in self.task:
    fo.write(key+" = "+obj2str(self.task[key])+"\n")
  if self.sites_from_file: 
      fo.write(" sites_from_file = true\n")
#("GS = true\ngap = false\ncorrelator = false\n}\n")
  # parameters of dmrg algorithm
  fo.write(" maxm = "+str(self.maxm)+"\n") # maximum bond dimension
  mpomaxm = max([self.maxm,self.mpomaxm]) # maximum
  fo.write(" mpomaxm = "+str(mpomaxm)+"\n") # maximum bond dimension
  fo.write(" noise = "+str(self.noise)+"\n") # maximum bond dimension
  fo.write(" cutoff = "+str(self.cutoff)+"\n") # maximum discarded weight
  fo.write(" nsweeps = "+str(self.nsweeps)+"\n") # maximum discarded weight
  ### this is a special addition to allow for generic interactions ###
  fo.write("}\n")
  fo.close()



def obj2str(a):
    """
    Convert into a string
    """
    if a is None: return "None"
    elif type(a)==int: return str(a)
    elif type(a)==bool:
        if a: return "true"
        else: return "false"
    elif type(a)==str: return a
    else: return str(a)

