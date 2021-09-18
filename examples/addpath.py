import os
ds = os.walk(os.getcwd())
d0 = os.getcwd()



def modify(ls,depth=1):
    """Modify an input file to add the path"""
    ls = ls.split("\n") # split
    lo = "# Add the root path of the dmrgpy library\n"
    dd = ""
    for i in range(depth): dd += "../"
    lo += "import os ; import sys ; sys.path.append(os.getcwd()+'/../"+dd+"src')\n"
#    lo += "\n"
    for l in ls: # loop over lines
      if "import os" in l or "import sys" in l or "sys.path.append" in l:
          continue
      if "__future__" in l: continue
      if "# Add the root path" in l: continue
      else: 
          lo += l + "\n" # add line
    return lo

ds = [d[0] for d in ds] # loop

def get_depth(d):
  d1 = d.replace(d0,"") # remove the parent directory
  o = 0
  for i in d1:
      if i=="/": o+=1
  return o

for d in ds:
  os.chdir(d) # go to that directory
  print(d)
  if os.path.isfile("main.py"):
      ls = open("main.py").read() # read all the lines
      print(d)
      open("main.py","w").write(modify(ls,depth=get_depth(d))) # write file
#  os.system("rm -f *.OUT")

