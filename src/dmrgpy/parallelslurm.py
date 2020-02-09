# this is a special library to launch parallel calculations using slurm

import dill as pickle
import os


def pcall(fin,xs):
    """Run a parallel calculation with slurm"""
    n = len(xs) # number of calculations
    f = lambda x: fin(x)
    main = "import dill as pickle\nimport os\n"
    main += "import sys ; sys.path.append(os.environ['DMRGROOT'])\n"
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
    pfolder = ".parallel"
    os.system("rm -rf "+pfolder) # create directory
    os.system("mkdir "+pfolder) # create directory

    pickle.dump(f,open(pfolder+"/function.obj","wb")) # write function
    pickle.dump(xs,open(pfolder+"/array.obj","wb")) # write object
    open(pfolder+"/run.py","w").write(main) # write script
    runsh = "#!/bin/bash\n#SBATCH -n 1\n#SBATCH -t 80:00:00\n"
    runsh += "#SBATCH --mem-per-cpu=5000\n"
    runsh += "#SBATCH --array=0-"+str(n-1)+"\n"
    runsh += "srun python run.py\n"
    open(pfolder+"/run.sh","w").write(runsh) # parallel file
    pwd = os.getcwd() # current directory 
    os.chdir(pfolder) # go to the folder
    os.system("sbatch run.sh >> run.out") # run calculation
    os.chdir(pwd) # back to main
    import time
    from os import path
    time.sleep(0.5) # wait half a second
    while True:
        finished = True
        for i in range(n):
            if not path.exists(pfolder+"/folder_"+str(i)+"/DONE"):
                finished = False
        if finished: break
    # get all the data
    ys = []
    for i in range(n):
        folder = pfolder+"/folder_"+str(i)+"/"
        y = pickle.load(open(folder+'out.obj','rb'))
        ys.append(y)
    return ys






