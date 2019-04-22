#!/usr/bin/python

from __future__ import print_function
import os



# different files for Linux and Mac
import platform


def addbashrc():
    """
    Add program to the .bashrc
    """
    pwd = path = os.path.dirname(os.path.realpath(__file__))+"/../" 
    if platform.system()=="Linux":
      bashrc = os.environ["HOME"]+"/.bashrc" # path to .bashrc
      print("Detected Linux system")
      f = open(bashrc).read() # read bashrc
    else:
      # Michael's fix for Mac
      f = ''
      for profileString in ["/.bash_profile", "/.bash_login","/.profile"]:
        bashrc = os.environ["HOME"]+profileString # path to .bashrc
        try:
          print(profileString)
          f = open(bashrc).read() # read bashrc
          break
        except FileNotFoundError:
          pass
      print("Detected Mac system")
    print("Adding DMRGROOT to your ",bashrc)
    
    route = "\n###############################\n"
    route += "# Added by spinflare\n"
    route += "###############################\n"
    route += "  export PATH=\""+pwd+"/bin:$PATH\"\n"
    route += "###############################\n"
    
    open(bashrc,"w").write(f+route) # write in the bashrc
    
    print("Added \n"+route+"\n route to ",bashrc)





if __name__=="__main__":
    addbashrc() # add to the .bashrc


