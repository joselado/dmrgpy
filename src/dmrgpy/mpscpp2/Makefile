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
# Phase 1 scaffold: pybind11 extension module (see bindings.cc). Builds
# _dmrgcpp<ext-suffix>.so against the same vendored ITensor static library
# used by $(APP).
#
# Deliberately dynamically linked against libstdc++ (NOT -static-libstdc++):
# an earlier version of this target used -static-libstdc++/-static-libgcc to
# work around conda Python bundling an older libstdc++ than this system
# compiler needs. That "fix" was itself the cause of a hard segfault (not
# just a missing-symbol error): once any *other* C++ extension in the same
# process (e.g. scipy.linalg/scipy.sparse, which dmrgpy itself imports)
# loads its own dynamic libstdc++, having a second, statically-embedded copy
# of libstdc++'s globals/exception tables/RTTI alongside it crashes. Plain
# dynamic linking avoids that; the original missing-symbol problem is
# instead handled at import time by cppext.py, which preloads the correct
# libstdc++ globally before anything else gets a chance to load a different
# one (see dmrgpy/__init__.py, which imports cppext first for exactly this
# reason).
#################################################################

PYTHON=python3
PYBIND_INCLUDES=$(shell $(PYTHON) -m pybind11 --includes)
EXT_SUFFIX=$(shell $(PYTHON)-config --extension-suffix)
EXT_MODULE=_dmrgcpp$(EXT_SUFFIX)

pybind: $(EXT_MODULE)

$(EXT_MODULE): bindings.cc $(ITENSOR_LIBS)
	$(CCCOM) -shared $(CCFLAGS) $(PYBIND_INCLUDES) bindings.cc -o $(EXT_MODULE) $(LIBFLAGS)

