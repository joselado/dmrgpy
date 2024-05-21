# this is a special library to launch parallel calculations using slurm

import dill as pickle
import os
from . import filesystem as fs
import signal
import subprocess

pickle.settings['recurse'] = True

def pcall(fin,xs,batch_size=1,**kwargs):
    """Wrapper to allow for a batch size"""
    if batch_size==1: 
        return pcall_rerun_failed(fin,xs,**kwargs)
    else: 
        nx = len(xs) # number of xs
        xsn = [] # empty list
        o = []
        for i in range(len(xs)):
            o.append(xs[i]) # store
            if i%batch_size==0: # reached the limit
                xsn.append(o) # store
                o = [] # reset
        def fnew(y): return [fin(x) for x in y] # call this batch
        outs = pcall_single(fnew,xsn,**kwargs) # call the inputs
        out = []
        for o in outs: out += o # add
        return out



def pcall_single(fin,xs,time=10,error=None):
    """Run a parallel calculation with slurm"""
    # create the parallel folder
    pfolder = os.getcwd()+"/.parallel"
    os.system("rm -rf "+pfolder) # remove directory
    fs.mkdir(pfolder) # create directory
    # now go with the rest
    n = len(xs) # number of calculations
#    f = lambda x: fin(x)
    f = fin
    from .manybodychain import dmrgpath
    realdmrgpath = dmrgpath+"/../"
    main = "import sys ; sys.path.append('"+realdmrgpath+"')\n"
    main += "import dill as pickle\nimport os\n"
    main += "try: ii = int(os.environ['SLURM_ARRAY_TASK_ID'])\n"
    main += "except: ii = 0\n"
    main += "f = pickle.load(open('function.obj','rb'))\n"
    main += "v = pickle.load(open('array.obj','rb'))\n"
    main += "folder = 'folder_'+str(ii)\n"
    main += "os.system('mkdir '+folder)\n"
    main += "os.chdir(folder)\n"
    main += "pwd = os.getcwd()\n"
#    main += "out = f(v[ii])\n"
    main += "try: out = f(v[ii])\n"
    main += "except: out = None\n"
    main += "os.chdir(pwd)\n"
    main += "print(out)\n"
    main += "pickle.dump(out,open('out.obj','wb'))\n"
    main += "os.system('touch DONE')\n"

    pickle.dump(f,open(pfolder+"/function.obj","wb")) # write function
    pickle.dump(xs,open(pfolder+"/array.obj","wb")) # write object
    open(pfolder+"/run.py","w").write(main) # write script
    hours = str(int(time)) # hours
    mins = int((time-int(time))*60)
    mins = str(max([mins,1])) # at least 1 minute
    runsh = "#!/bin/bash\n#SBATCH -n 1\n#SBATCH -t "+str(int(time))+":"+str(mins)+":00\n"
    runsh += "#SBATCH --mem-per-cpu=5000\n"
    runsh += "#SBATCH --array=0-"+str(n-1)+"\n"
    runsh += "srun python run.py\n"
    open(pfolder+"/run.sh","w").write(runsh) # parallel file
    pwd = os.getcwd() # current directory 
    os.chdir(pfolder) # go to the folder
#    os.system("sbatch run.sh >> run.out") # run calculation
    env = get_env() # get the cleaned environment
    out,err = subprocess.Popen(["sbatch","run.sh"],stdout=subprocess.PIPE,env=env).communicate()
    job = job_number(out) # job number
    jobkill(job) # kill the job if exiting
    os.chdir(pwd) # back to main
    import time
    from os import path
    time.sleep(0.5) # wait half a second
    while True:
        time.sleep(1.) # wait one second
        finished = True
        for i in range(n): # check all the folders
            if not path.exists(pfolder+"/folder_"+str(i)+"/DONE"):
                finished = False
        if path.exists(pfolder+"/STOP"): # stop by force .parallel
            finished = True
        if finished: break # exit the loop
    kill_process(job) # kill all the processes just in case
    # get all the data
    ys = []
    for i in range(n): # loop over outputs
        folder = pfolder+"/folder_"+str(i)+"/"
        try: y = pickle.load(open(folder+'out.obj','rb')) # get the output
        except: y = None # in case there are errors
        if y is None: y = error # use this as backup variable
        ys.append(y) # store in the output
    time.sleep(1.) # wait one second
    return ys # return all the data 




def rerun_failed(fun,ys,xs,**kwargs):
    """Rerun those calculations that failed"""
    nc = len(ys) # number of calculations
    dd = dict() # dictionary
    for i in range(nc): # loop over calculations
        if ys[i] is None: # this one failed
            dd[i] = xs[i] # store this input to run
            print("Job ",i,xs[i],"failed, rerunning")
    # run again all the calculations that failed
    xsnew = [dd[key] for key in dd] # inputs that failed
    if len(xsnew)==0: # all calculations ok
        return ys # return the outputs
    else: # some have failed
        ysnew = pcall_single(fun,xsnew,**kwargs) # call the function
        # store in the array
        ii = 0 # counter
        for key in dd: # loop over keys
            ys[key] = ysnew[ii] # store
            ii += 1
        return ys # return outputs


def pcall_rerun_failed(fun,xs,call_tries=10,**kwargs):
    """Call for all the inputs, and rerun those that failed"""
    ys = pcall_single(fun,xs,**kwargs) # call the function
    for ic in range(call_tries): # call several times
        ys = rerun_failed(fun,ys,xs,**kwargs) # rerun those that failed
    return ys # return all




def get_env():
    """Get a cleaned up environment for slurm"""
    env = os.environ # dictionary
    envout = {} # output dictionary
    for key in env: # loop over environmental variables
        if "SLURM_" not in key and "SBATCH_" not in key:
           envout[key] = env[key] # store
    return envout # return this dictionary


def job_number(out):
    """Get the job number"""
    out = str(out)
    print(out)
    out = out.split("Submitted batch job")[1]
    print(out)
    out = out.split("\\n")[0]
    return int(out) # return the job


def jobkill(n):
    """Kill the job when the program is killed"""
    def killf(*args):
      subprocess.Popen(["scancel",str(n)],stdout=subprocess.PIPE).communicate()
      print("Job killed")
      exit()
    signal.signal(signal.SIGINT, killf)
    signal.signal(signal.SIGTERM, killf)



def kill_process(n):
    """Command to kill the prcesses"""
    subprocess.Popen(["scancel",str(n)],stdout=subprocess.PIPE).communicate()
    print("Job killed")



