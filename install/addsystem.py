#!/usr/bin/python

from __future__ import print_function
import os




# different files for Linux and Mac
import platform


def addbashrc():
    """
    Add program to the .bashrc
    """
    shell = os.environ["SHELL"]
    if "bash" in shell: shellrc = ".bashrc"
    elif "zsh" in shell: shellrc = ".zshrc"
    else:
        print("Unrecognise shell",name)
        exit()
    pwd = os.path.dirname(os.path.realpath(__file__))+"/../" 
    if platform.system()=="Linux":
      shellrc = os.environ["HOME"]+"/"+shellrc # path to .bashrc
      print("Detected Linux system")
      f = open(shellrc).read() # read bashrc
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
    print("Adding DMRGROOT to your ",shellrc)
    
    route = "\n\n###############################\n"
    route += "  export DMRGROOT=\""+pwd+"/src\"\n"
    route += "###############################\n"
    
    open(shellrc,"w").write(f+route) # write in the bashrc
    
    print("Added \n"+route+"\n route to ",shellrc)





if __name__=="__main__":
    addbashrc() # add to the .bashrc


