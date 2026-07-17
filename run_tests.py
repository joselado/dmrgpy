#!/usr/bin/python

# Run the full pytest suite under tests/, then clean up the working
# directories (.mpsfolder, .pychainfolder, .dmrgfolder, ...) that DMRG/ED
# runs leave behind, the same way clean.py does after running examples.

import os
import subprocess
import sys

pwd = os.path.dirname(os.path.abspath(__file__))

ret = subprocess.call([sys.executable, "-m", "pytest", "tests"] + sys.argv[1:],
                       cwd=pwd)

os.system(f"python \"{os.path.join(pwd, 'clean.py')}\"")

sys.exit(ret)
