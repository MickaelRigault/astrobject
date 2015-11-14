#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""This module gather the basic tools used in a lot of methods"""

import numpy as np

__all__ = ["kwargs_update",
           "load_pkl","dump_pkl"]


def kwargs_update(default,**kwargs):
    """
    """
    k = default.copy()
    for key,val in kwargs.iteritems():
        k[key] = val
        
    return k

def load_pkl(filename):
    """
    """
    import cPickle as pkl
    try:
        pkl_file = open(filename,'rb')
    except:
        raise IOError("The given file does not exist %s"%filename)
    
    return pkl.load(pkl_file)


def dump_pkl(data,filename,**kwargs):
    """
    """
    from cPickle import dump
    if len(filename.split("."))>1 and filename.split(".")[-1]=="pkl":
        outfile =  open(filename,"wb")
    else:
        outfile =  open(filename+".pkl","wb")
    
    dump(data, outfile,**kwargs)
    outfile.close()


def shape_ajustment(X,Y,model_X,k=4,verbose=False):
    """ DOC TO BE DONE
    return  Y /w the same binning as model_X  in their commun wavelength
    """
    def  commun_wavelength(x1,x2):
        """
        """
        flagx1 = (x1>=x2.min()) & (x1<=x2.max())
        flagx2 = (x2>=x1.min()) & (x2<=x1.max())
        return flagx1,flagx2
        
    from scipy.interpolate import UnivariateSpline
    flagX,flagmodel = commun_wavelength(X,model_X)
    Yrebin = UnivariateSpline(X[flagX], Y[flagX],k=k,s=0)(model_X)
    
    if len(Yrebin)==len(model_X):
        return Yrebin
    else:
        if verbose:
            print 'WARNING [shape_adjustment] non-mached shape ... I am fixing that'
        
        YrebinOK = np.empty((len(Yrebin)+1),)
        YrebinOK[1:] = Yrebin
        YrebinOK[0]  = Yrebin[0]
        
        return YrebinOK
