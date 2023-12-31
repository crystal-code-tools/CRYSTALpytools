from CRYSTALpytools import adsorb
from CRYSTALpytools import base
from CRYSTALpytools import calculate
from CRYSTALpytools import convert
from CRYSTALpytools import execute
from CRYSTALpytools import geometry
from CRYSTALpytools import crystal_io
from CRYSTALpytools import plot
from CRYSTALpytools import thermodynamics
from CRYSTALpytools import units
from CRYSTALpytools import utils

import sys

if sys.version_info >= (3, 8):
    from importlib import metadata
else:
    from importlib_metadata import metadata

__author__ = "CRYSTALpytools Development Team"
__email__ = "crystalcodetools@gmail.com"
__version__ = metadata.version('CRYSTALpytools')
