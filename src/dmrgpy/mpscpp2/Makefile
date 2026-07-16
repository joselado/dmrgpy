LIBRARY_DIR=ITensor/

APP=mpscpp

HEADERS=

CCFILES=$(APP).cc

#################################################################
#################################################################
#################################################################
#################################################################

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

#Mappings --------------
OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

#Rules ------------------

%.o: %.cc $(HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: $(APP)
debug: $(APP)-g

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o $(APP) $(LIBFLAGS)

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS)

clean:
	rm -fr .debug_objs *.o $(APP) $(APP)-g _dmrgcpp*.so

mkdebugdir:
	mkdir -p .debug_objs

#################################################################
# pybind11 extension module (see bindings.cc). Builds _dmrgcpp<ext-suffix>.so
# against the same vendored ITensor static library used by $(APP).
#
# Built with a *different* compiler than $(APP) when Python's own
# environment provides one (PYBIND_CCCOM below), specifically to avoid a
# real, reproducible segfault: $(CCCOM) is whatever system compiler builds
# ITensor/$(APP), which never has to coexist with anything else in a
# process. This extension does -- it's loaded into the same Python process
# as numpy/scipy, which (at least under conda) load their own bundled
# libstdc++/BLAS. Compiling this extension with the *system* compiler needs
# a libstdc++ newer than conda's; a couple of things were tried and
# rejected before landing here:
#   - static-linking libstdc++ into the extension: crashes once numpy/scipy
#     load a second, different libstdc++ into the same process (duplicate
#     globals/exception tables/RTTI).
#   - preloading the newer system libstdc++ at Python import time
#     (ctypes.CDLL(...,RTLD_GLOBAL)): only works if nothing else loaded a
#     conflicting libstdc++ first, which numpy alone already does in
#     practice (see git history for cppext.py) -- not reliable enough.
#   - renaming/isolating BLAS+LAPACK+libstdc++ into private-SONAME copies
#     via patchelf: works for BLAS/LAPACK, but doing the same to libstdc++
#     hits the same duplicate-runtime crash as static linking.
# The fix that actually works: build with a compiler whose own libstdc++
# *is* the one Python's environment already provides (e.g. conda's own
# gxx_linux-64 package), so there's only ever one libstdc++/BLAS image in
# the process, regardless of import order. Confirmed empirically to survive
# the worst case (scipy.linalg/scipy.sparse imported before dmrgpy).
# Falls back to $(CCCOM) if no such compiler is found -- same-compiler
# builds (e.g. a plain venv without a bundled libstdc++/BLAS) don't have
# this problem in the first place.
#################################################################

PYTHON=python3
PYBIND_INCLUDES=$(shell $(PYTHON) -m pybind11 --includes)
EXT_SUFFIX=$(shell $(PYTHON)-config --extension-suffix)
EXT_MODULE=_dmrgcpp$(EXT_SUFFIX)

CONDA_GXX=$(shell d=$$(dirname $$(command -v $(PYTHON)) 2>/dev/null); ls $$d/*-conda*-g++ 2>/dev/null | head -1)
ifeq ($(CONDA_GXX),)
PYBIND_CCCOM=$(CCCOM)
PYBIND_LIBFLAGS=$(LIBFLAGS)
else
PYBIND_CCCOM=$(CONDA_GXX) -m64 -std=c++14 -fPIC
PYBIND_LIBFLAGS=-L$(ITENSOR_LIBDIR) -litensor -lpthread -L/usr/lib/x86_64-linux-gnu -lblas -llapack
endif

pybind: $(EXT_MODULE)

$(EXT_MODULE): bindings.cc $(ITENSOR_LIBS)
	$(PYBIND_CCCOM) -shared $(CCFLAGS) $(PYBIND_INCLUDES) bindings.cc -o $(EXT_MODULE) $(PYBIND_LIBFLAGS)

