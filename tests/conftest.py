import os
import sys

# Make "import dmrgpy" work when running pytest from a source checkout,
# without requiring the package to be installed / PYTHONPATH-configured.
_SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
