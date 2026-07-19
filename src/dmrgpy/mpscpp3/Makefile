LIBRARY_DIR=ITensor/

#################################################################
#################################################################
#################################################################
#################################################################

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

#Targets -----------------

build: pybind

clean:
	rm -fr .debug_objs *.o _dmrgcpp*.so

#################################################################
# pybind11 extension module (see bindings.cc). Builds _dmrgcpp<ext-suffix>.so
# against the vendored ITensor v3 static library.
#
# Same rationale as mpscpp2/Makefile (see the long comment there): built
# with a *different* compiler than the one that builds ITensor itself when
# Python's own environment provides one (e.g. conda's own gxx_linux-64),
# so the extension's libstdc++/BLAS matches the one already loaded by
# numpy/scipy in the same process. The only difference from mpscpp2's copy
# of this Makefile is the language standard: ITensor v3 requires C++17 with
# the concepts TS extension (-fconcepts), not v2's C++14 -- see
# ITensor/options.mk.sample's CCCOM line.
#################################################################

PYTHON=python3
PYBIND_INCLUDES=$(shell $(PYTHON) -m pybind11 --includes)
# Ask the interpreter itself for the suffix (not $(PYTHON)-config):
# conda ships python3-config but no python-config, so the -config
# form breaks whenever PYTHON is a bare 'python'.
EXT_SUFFIX=$(shell $(PYTHON) -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")
EXT_MODULE=_dmrgcpp$(EXT_SUFFIX)

CONDA_GXX=$(shell d=$$(dirname $$(command -v $(PYTHON)) 2>/dev/null); ls $$d/*-conda*-g++ 2>/dev/null | head -1)
ifeq ($(CONDA_GXX),)
PYBIND_CCCOM=$(CCCOM)
PYBIND_LIBFLAGS=$(LIBFLAGS)
else
PYBIND_CCCOM=$(CONDA_GXX) -m64 -std=c++17 -fconcepts -fPIC
# BLAS/LAPACK flags come from options.mk (BLAS_LAPACK_LIBFLAGS), so this
# branch stays in sync with whatever backend ITensor itself was built
# against (plain lapack, OpenBLAS, MKL, ...) instead of hardcoding
# -lblas -llapack, which silently ignored that configuration.
PYBIND_LIBFLAGS=-L$(ITENSOR_LIBDIR) -litensor $(BLAS_LAPACK_LIBFLAGS)
endif

pybind: $(EXT_MODULE)

$(EXT_MODULE): bindings.cc $(ITENSOR_LIBS)
	$(PYBIND_CCCOM) -shared $(CCFLAGS) $(PYBIND_INCLUDES) bindings.cc -o $(EXT_MODULE) $(PYBIND_LIBFLAGS)

