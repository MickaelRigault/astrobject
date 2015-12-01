#
# Init the astro objects
#
from baseobject   import *
from spectroscopy import *
from transient    import *
from photometry   import *
from multiband    import *
import instruments

__all__ = ['baseobject','spectroscopy','photometry',
           'multiband',"transient","instruments"]
