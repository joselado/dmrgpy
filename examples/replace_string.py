import os
ds = os.walk(os.getcwd())



def modify(ls):
    """Modify an input file to add the path"""
    return ls.replace("time.clock","time.perf_counter")

ds = [d[0] for d in ds] # loop

for d in ds:
  os.chdir(d) # go to that directory
  if os.path.isfile("main.py"):
      ls = open("main.py").read() # read all the lines
      print(d)
      open("main.py","w").write(modify(ls)) # write file
#  os.system("rm -f *.OUT")

