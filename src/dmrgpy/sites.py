# initialize the sites for the C++ calculation
from .writemps import write_sites
import os

def initialize(self,**kwargs):
    self.path = os.getcwd()+"/.mpsfolder/" # folder of the calculations
    self.clean() # clean calculation
    self.inipath = os.getcwd() # original folder
    os.system("mkdir -p "+self.path) # create folder for the calculations
    self.sites_from_file = False
    self.task = {"write_sites":"true"}
    self.execute(lambda: write_sites(self)) # write the different sites
    self.run() # run the calculation
    from .mode import get_mode
    if not get_mode(self)=="ED":
        self.bin_sites = open(self.path+"/sites.sites","rb").read()
    self.sites_from_file = True




