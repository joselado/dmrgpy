# Import cppext first, before any submodule below has a chance to import
# scipy: cppext preloads the system libstdc++ the optional in-process C++
# extension (mpscpp2/_dmrgcpp) needs, and that preload has to win the race
# against scipy's own libstdc++ or using the extension later can segfault
# (see cppext.py). Harmless if the extension is never built/used.
from . import cppext
