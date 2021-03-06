# Create a new Instrument

You want to create a new instrument called `newinst`.



## File name and mandatory functions

**The filename** of the new instrument must be called `newinst.py` and should be located with the other instruments inside `astrobject/instruments/`

There are 3 **mandarory functions** inside the `newinst.py` file:

- `newinst(*args, **kwargs)`: It returns the object corresponding the the newinst (see following section)
- `is_newinst_file(filename)`: Tests if the given file belongs to this instrument

- `which_band_is_file(filename)`: Returns the band associated to the given filename.

## Build the Instrument Class

you must import the mother object `Instrument`:
``` python
from .baseinstrument import Instrument
```

Create the class:
```
class NewInst( Instrument ):
	instrument_name=newinst
```


The other methods depends how the input data are build. For standard format, most of the work consist of telling the software where to find the information inside the header.
The imporante one are:

- observation date:
```
@property
def mjd(self):
	""" The observation data in Modified Julian Date"""
	return TELL_HERE_HOW_TO_GET_THE_DATE_IN_MJD
```
- bandname
```
@property
def bandname(self):
	""" The name of the band """
	return TELL_HERE_HOW_TO_GET_THE_BANDNAME_IN_A_WAY_SNCOSMO_KNOWS
```

- zeropoint
```
@property
def mab0(self):
	""" The zeropoint"""
	return TELL_HERE_HOW_TO_GET_THE_ZEROPOINT_HEADER_KEY_OR_FIXED_VALUE?
```

- (gain, could be important)
```
@property
def _gain(self):
	""" The Gain of the instrument """
	return TELL_HERE_HOW_TO_GET_GAIN_IF_ANY__OTHERWISE_RETURN None
```

### The Variance

The Variance is a key parameter and the way it is stored is really instrument dependent.
If there is a simple way to get the variance from the loaded file (stored as `self.fits`) than simple do:
```
@property
def var(self):
	""" The variance associated with the data"""
	return TELL_HERE_HOW_TO_GET_THE_VARIANCE.
	# For instance if the variance is the second entry of the fitsfile:
	return self.fits[1].data
	# If the error is stored as the second entry
	return self.fits[1].data**2
```

Things can be more complicated, and sometime the variance is stored in a second file. Then I suggests defining a new method call `set_variance(self, filename)` such that:

```
from astropy.io import fits as pf
def set_variance(self, filename):
	""" Set the variance associated with the data """
	self._properties["var"] = pf.getdata(filename)
```
This will do the job since by default `self.var` returns  `self._properties['var']` (see `photometry.py`)

### Data mask

Some instrument have masked pixels. Do the same way as for the variance, but this time with the key:
```
@property
def datamask(self):
	...
```
and the associated default access in `self._side_properties["datamask"]`


### Data given in counts and not counts per second

Some instruments provide data in counts and not counts\_per\_second. I think we should be consistent in the soft and already work in counts\_per\_second. To do that simply do the following:

```
@property
def rawdata(self):
	""" """
	return super(NEWINST, self).rawdata / self.exposuretime
```
(this requires that you have defined exposuretime. Not that here, exposuretime could be an 2D array)

Also, be careful that the zeropoint does follow and is the zeropoint for counts\_per\_second.

## What if the data are not in the first fits entry?

If the data entry is, for instance in the DATAINDEX entry of the fits file, then overwrite the `__build__()` method as follow:
```
def __build__(self,**kargs):
        """  """
        super(NewInst,self).__build__(**kwargs)
        self._build_properties = kwargs_update(self._build_properties,
                                               **dict(
                    data_index = DATAINDEX,
                    ))
```

# Final Step

Once you are happy, simply add your new instrument to the `astrobject/instruments/instrument.py`. To do so:
add the import:
```
import newinst
```
and add it to the list of `KNOWN_INSTRUMENTS`

# How to know more?

Check other files. The hst one is pretty easy while panstarrs has a lot a tricky things like external variance map and data given in counts.

The inheritance channel is as follow:

Image (from photometry.py) -> Instrument (from instruments/baseinstruments.py) -> NEWINST


# Example WISE:
```
__all__  = ["wise", "WISE_INFO"]

WISE_INFO = {"bands":["w1","w2","w3","w4"],
"w1":{}}

# -------------------- #
# - Instrument Info  - #
# -------------------- #

import numpy as np
from astropy.io import fits as pf
from .baseinstrument import Instrument
from ..utils.tools import kwargs_update


DATAINDEX = 0


def wise(*args,**kwargs):
    return WISE(*args,**kwargs)

def is_wise_file(filename):
    """This test if the input file is a WISE one"""
    return pf.getheader(filename).get("ORIGIN") == "SDSS"

def which_band_is_file(filename):
    """This resuts the band of the given file if it is a
    sdss one"""
    if not is_wise_file(filename):
        return None
    return pf.getheader(filename).get("FILTER")

def which_obs_mjd(filename):
    """ read the sdss-filename and return the
    modified julian date """
    if not is_wise_file(filename):
        return None
    return get_mjd(pf.getheader(filename))

#
#   def _get_default_background_(self,*args,**kwargs):
            return np.zeros(np.shape(self.rawdata))

#  data, rawdata, background
#  data = rawdata - background


class WISE( Instrument ):
    """ """
    instrument_name = "WISE"	
    @property
    def var(self):
        """ - Not mandatory this is the variance image - """
	    return YOUCOULDFINDAWAYDATAINDEX
	# You could add any other thing this way.
	# E.G. "sky"   :
	# At the beginning of the class add PROPERTIES = ["sky"]
    # @property
	# def sky(self):
	#   """ This is the sky 2D image"""
	#   return self._properties['sky']
    # def set_sky(self, 2Diamge):
	#       """ """
	#       self._properties['sky'] = 2Diamge


    # -- Derived values
	@property
	def exposuretime(self):
	     return FINDAWAY

	@property
    def mab0(self):
          return FINDAWAY

    @property
    def bandname(self):
	     if self._properties['bandname'] is None:
                  self._properties['bandname'] = FINDAWAY
	     return self._properties['bandname']

    @property
    def mjd(self):
        """This is the Modify Julien Date at the start of the Exposure"""
        return self.fits[0].header["EXPSTART"]


    @property
    def _gain(self):
        return 1.

    @property
    def _dataunits_to_electron(self):
        """The gain converts ADU->electron. The Data shouls be in ADU/s"""
        return self._gain * self.exposuretime
    
    @property
    def mjd(self):
        raise NotImplementedError("'mjd' (Modified Julian Date) must be implemented")
	
```
