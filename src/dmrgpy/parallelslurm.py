# this is a special library to launch parallel calculations using slurm

import dill as pickle
import os


def pcall(fin,xs):
    """Run a parallel calculation with slurm"""
    n = len(xs) # number of calculations
    f = lambda x: fin(x)
    pwd = os.getcwd() # current directory 
    main = "import dill as pickle\nimport os\n"
    main += "try: ii = int(os.environ['SLURM_ARRAY_TASK_ID'])\n"
    main += "except: ii = 0\n"
    main += "f = pickle.load(open('function.obj','rb'))\n"
    main += "v = pickle.load(open('array.obj','rb'))\n"
    main += "folder = 'folder_'+str(ii)\n"
    main += "os.system('mkdir '+folder)\n"
    main += "os.chdir(folder)\n"
    main += "out = f(v[ii])\n"
#    main += "try: out = f(v[ii])\n"
#    main += "except: out = None\n"
    main += "print(out)\n"
    main += "pickle.dump(out,open('out.obj','wb'))\n"
    pfolder = ".parallel"
    os.system("rm -rf "+pfolder) # create directory
    os.system("mkdir "+pfolder) # create directory

    pickle.dump(f,open(pfolder+"/function.obj","wb")) # write function
    pickle.dump(xs,open(pfolder+"/array.obj","wb")) # write object
    open(pfolder+"/run.py","w").write(main) # write script
    runsh = "#!/bin/bash\n#SBATCH -n 1\n#SBATCH -t 20:00:00"
    runsh += "#SBATCH --mem-per-cpu=5000\n"
    runsh += "#SBATCH --array=0-"+str(n-1)+"\n"
    runsh += "srun python main.py\n"
    open(pfolder+"/run.sh","w").write(runsh) # parallel file
    os.chdir(pfolder) # go to the folder
#    os.system("sbatch run.sh") # run calculation





