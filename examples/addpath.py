import os
ds = os.walk(os.getcwd())



def modify(ls):
    """Modify an input file to add the path"""
    ls = ls.split("\n") # split
    lo = "# Add the root path of the dmrgpy library\n"
    lo += "import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')\n"
    lo += "\n"
    for l in ls: # loop over lines
      if "import os" in l or "import sys" in l or "sys.path.append" in l:
          continue
      if "__future__" in l: continue
      if "# Add the root path" in l: continue
      else: 
          if l!="": lo += l + "\n" # add line
    return lo

ds = [d[0] for d in ds] # loop

for d in ds:
  os.chdir(d) # go to that directory
  if os.path.isfile("main.py"):
      ls = open("main.py").read() # read all the lines
      print(d)
      open("main.py","w").write(modify(ls)) # write file
#  os.system("rm -f *.OUT")

