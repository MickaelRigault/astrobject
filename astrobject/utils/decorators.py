#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""This module gather the customed decorators of astrobject"""

__all__ = ["_autogen_docstring_inheritance"]


#####################################
# == Stoolen from Matplotlib     == #
#####################################
def _autogen_docstring_inheritance(base):
    """Autogenerated wrappers will get their docstring from a base function
    with an addendum."""
    return base

def make_method(obj):
    """Decorator to make the function a method of *obj*.

    In the current context::
      @make_method(Axes)
      def toto(ax, ...):
          ...
    makes *toto* a method of `Axes`, so that one can directly use::
      ax.toto()
    COPYRIGHT: from Yannick Copin
    """

    def decorate(f):
        setattr(obj, f.__name__, f)
        return f

    return decorate


def speed_test(func):
    def wrapper(*args, **kwargs):
        t1 = time.time()
        for x in xrange(5000):
            results = func(*args, **kwargs)
        t2 = time.time()
        print('%s took %0.3f ms' % (func.func_name, (t2-t1)*1000.0))
        return results
    return wrapper
