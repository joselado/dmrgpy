# Single source of truth for which C++ DMRG backend `install.py` compiles
# by default. Kept in sync by hand with src/dmrgpy/cppext.py's own
# DEFAULT_ITENSOR_VERSION (the two can't share one constant: this script
# runs before src/ is on the Python path) -- change both together.
DEFAULT_ITENSOR_VERSION = 3
