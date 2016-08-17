
#
# This is the astrobject library
#


#import collections
#import spectr
#from photometry import *
#from instruments.instrument import *

__version__ = "0.3.5"

from baseobject   import *


from photometry   import *
from spectroscopy import *

# don't invert collections and collection
from collections import *
from collection  import *
from instruments.instrument import *



