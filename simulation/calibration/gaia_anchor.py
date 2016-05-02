#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" Based on provided simumation this module enables to retreive the input values """

"""
Simulated values input for the fitting techniques:

* pixel coordinates: x, y
* observed magnitudes: g, r, g_err, r_err (i, i_err)
* ref magnitudes (gaia): gref, rref, gref_err, rref_err (iref,iref_err)
* observational condition: airmass
* Bandpasses: ZTF, Gaia
"""



import numpy as np
import matplotlib.pyplot as mpl
