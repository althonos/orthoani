import warnings

import pyorthoani
from pyorthoani import *

__version__ = pyorthoani.__version__
__author__ = pyorthoani.__author__
__license__ = pyorthoani.__license__


warnings.warn(
    "orthoani was renamed to pyorthoani starting from v0.6.1, "
    "please update your code and dependencies",
    DeprecationWarning,
    stacklevel=2
)
