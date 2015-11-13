#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""This module gather the basic tools used in a lot of methods"""


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
